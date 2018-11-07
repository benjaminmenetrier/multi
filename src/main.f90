program main

use algo
use bmatrix
use fft
use hmatrix
use interp
use netcdf
use rand

implicit none

! Parameters
integer,parameter :: no = 3
integer,parameter :: ni = 5
integer,parameter :: n = 128
real(8),parameter :: sigmao = 0.1
logical :: full_res = .false.
logical :: specific_interp = .true.
logical :: no_update = .false.
character(len=1024) :: lmp = 'spectral'

! Local variables
integer :: nobs,io,jo,ii
integer,allocatable :: fac(:),nn(:),nprec(:)
real(8),allocatable :: yo(:),d(:),hxg(:)
real(8),allocatable :: xb(:),xg(:),xa(:),dxa(:),dxbbar(:,:),dxabar(:,:),dxabar_interp(:)
real(8),allocatable :: sigmab(:,:),spvar(:,:)
real(8),allocatable :: jb_planczos(:,:),jo_planczos(:,:),jb_planczosif(:,:),jo_planczosif(:,:)
real(8),allocatable :: xa_mse_planczos(:,:),xa_mse_planczosif(:,:)
real(8),allocatable :: dx_planczos(:,:,:),dx_planczosif(:,:,:)
real(8),allocatable :: ritzval(:,:),dxbarritzvec(:,:,:),dxritzvec(:,:,:)
real(8),allocatable :: vb(:),dvb(:,:),dva(:,:),dva_interp(:),vritzvec(:,:,:)
real(8),allocatable :: precvec(:,:,:,:)

! Set seed
call set_seed(.false.)

! FFT test
call fft_test(n)

! Interpolation test
call interp_test(n)

! Number of observations
nobs = n/2**no
write(*,'(a,i4)') 'Number of observations:                 ',nobs
if (n/=nobs*2**no) then
   write(*,'(a)') 'Error: n should be related to the number of outer iterations'
   stop
end if

! H matrix test
call h_test(n,nobs)

! Allocation
allocate(fac(no))
allocate(nn(no))
allocate(nprec(no))
allocate(yo(nobs))
allocate(d(nobs))
allocate(hxg(nobs))
allocate(sigmab(n,no))
allocate(spvar(n,no))
allocate(xb(n))
allocate(xg(n))
allocate(xa(n))
allocate(dxa(n))
allocate(dxbbar(n,no))
allocate(dxabar(n,no))
allocate(dxabar_interp(n))
allocate(dx_planczos(n,ni,no))
allocate(dx_planczosif(n,ni,no))
allocate(jb_planczos(0:ni,no))
allocate(jo_planczos(0:ni,no))
allocate(jb_planczosif(0:ni,no))
allocate(jo_planczosif(0:ni,no))
allocate(xa_mse_planczos(0:ni,no))
allocate(xa_mse_planczosif(0:ni,no))
allocate(ritzval(ni,no))
allocate(dxbarritzvec(n,ni,no))
allocate(dxritzvec(n,ni,no))
allocate(vb(n))
allocate(dvb(n,no))
allocate(dva(n,no))
allocate(dva_interp(n))
allocate(vritzvec(n,ni,no))
allocate(precvec(n,ni,no,2))

! Setup observations
call rand_normal(nobs,yo)
yo = yo*sigmao

! Full resolution background error parameters
call gp_variance(n,sigmab(:,no))
call sp_variance(n,spvar(:,no))

! B matrix test
call b_test(n,sigmab(:,no),spvar(:,no),n/nobs)

! Background state
call rand_normal(n,vb)
call apply_u(n,sigmab(:,no),spvar(:,no),vb,xb)

! Multi-incremental PLanczos
write(*,'(a)') 'Multi-incremental PLanczos'
do io=1,no
   ! Resolution
   if (full_res) then
      fac(io) = 1
   else
      fac(io) = 2**(no-io)
   end if
   nn(io) = n/fac(io)
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,', resolution: ',fac(io)

   ! B parameters
   call gp_variance(nn(io),sigmab(1:nn(io),io))
   call sp_variance(nn(io),spvar(1:nn(io),io))

   ! Initialization
   dvb(:,io) = 0.0
   dva(:,io) = 0.0

   ! Guess
   if ((io==1).or.no_update) then
      xg(1:nn(io)) = xb(1:nn(io))
   else
      xa(1:nn(io-1)) = xg(1:nn(io-1))+dxa(1:nn(io-1))
      call interp_gp(nn(io-1),xa(1:nn(io-1)),nn(io),.true.,xg(1:nn(io)))
   end if

   ! Innovation
   call apply_h(nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   ! Background increment
   if (.not.no_update) then
      do jo=1,io-1
         call interp_sp(nn(jo),dva(1:nn(jo),jo),nn(io),.true.,dva_interp(1:nn(io)))
         dvb(1:nn(io),io) = dvb(1:nn(io),io)-dva_interp(1:nn(io))
      end do
   end if

   if (io==1) then
      ! No preconditioning vectors
      nprec(io) = 0
   else
      ! Preconditioning vectors
      select case (trim(lmp))
      case ('none')
         ! No preconditioning vectors
         nprec(io) = 0
      case ('spectral')
         do jo=1,io-1
            nprec(jo) = 0
            do ii=1,ni
               if (ritzval(ii,jo)>1.0+1.0e-3) then
                  nprec(jo) = nprec(jo)+1
                  call interp_sp(nn(jo),(1.0-1.0/ritzval(ii,jo))*vritzvec(1:nn(jo),ii,jo),nn(io),.false.,precvec(1:nn(io),nprec(jo),jo,1))
                  call interp_sp(nn(jo),vritzvec(1:nn(jo),ii,jo),nn(io),.false.,precvec(1:nn(io),nprec(jo),jo,2))
               end if
            end do
         end do
      case default
         write(*,'(a)') 'Error: wrong LMP'
         stop
      end select
   end if

   ! Minimization
   call planczos(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),sigmao,dvb(1:nn(io),io),nobs,d,ni,io,nprec(1:io),precvec(1:nn(io),:,1:io,:),dva(1:nn(io),io),dx_planczos(1:nn(io),:,io),jb_planczos(:,io),jo_planczos(:,io),ritzval(:,io),vritzvec(1:nn(io),:,io))

   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',jb_planczos(ii,io)+jo_planczos(ii,io),' = ',jb_planczos(ii,io),' + ',jo_planczos(ii,io)
   end do

   ! New analysis increment
   call apply_u(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),dva(1:nn(io),io),dxa(1:nn(io)))

   ! MSE
   xa_mse_planczos(0,io) = sum(xg(1:nn(io))**2)/real(nn(io))
   do ii=1,ni
      xa_mse_planczos(ii,io) = sum((xg(1:nn(io))+dx_planczos(1:nn(io),ii,io))**2)/real(nn(io))
   end do
end do
write(*,*) 

! Multi-incremental PLanczosIF
write(*,'(a)') 'Multi-incremental PLanczosIF'
do io=1,no
   ! Resolution
   if (full_res) then
      fac(io) = 1
   else
      fac(io) = 2**(no-io)
   end if
   nn(io) = n/fac(io)
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)

   ! B parameters
   call gp_variance(nn(io),sigmab(1:nn(io),io))
   call sp_variance(nn(io),spvar(1:nn(io),io))

   ! Initialization
   dxbbar(:,io) = 0.0
   dxabar(:,io) = 0.0

   ! Guess
   if ((io==1).or.no_update) then
      xg(1:nn(io)) = xb(1:nn(io))
   else
      xa(1:nn(io-1)) = xg(1:nn(io-1))+dxa(1:nn(io-1))
      call interp_gp(nn(io-1),xa(1:nn(io-1)),nn(io),.true.,xg(1:nn(io)))
   end if

   ! Innovation
   call apply_h(nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   ! Background increment
   if (.not.no_update) then
      do jo=1,io-1
         if (specific_interp) then
            call interp_gp_b_inv(nn(jo),sigmab(1:nn(jo),jo),spvar(1:nn(jo),jo),dxabar(1:nn(jo),jo),nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),.true.,dxabar_interp(1:nn(io)))
         else
            call interp_gp(nn(jo),dxabar(1:nn(jo),jo),nn(io),.true.,dxabar_interp(1:nn(io)))
         end if
         dxbbar(1:nn(io),io) = dxbbar(1:nn(io),io)-dxabar_interp(1:nn(io))
      end do
   end if

   if (io==1) then
      ! No preconditioning vectors
      nprec(io) = 0
   else
      ! Preconditioning vectors
      select case (trim(lmp))
      case ('none')
         ! No preconditioning vectors
         nprec(io) = 0
      case ('spectral')
         do jo=1,io-1
            nprec(jo) = 0
            do ii=1,ni
               if (ritzval(ii,jo)>1.0+1.0e-3) then
                  nprec(jo) = nprec(jo)+1
                  if (specific_interp) then
                     call interp_gp_b(nn(jo),sigmab(1:nn(jo),jo),spvar(1:nn(jo),jo),(1.0-1.0/ritzval(ii,jo))*dxritzvec(1:nn(jo),ii,jo),nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),.false.,precvec(1:nn(io),nprec(jo),jo,1))
                     call interp_gp_b_inv(nn(jo),sigmab(1:nn(jo),jo),spvar(1:nn(jo),jo),dxbarritzvec(1:nn(jo),ii,jo),nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),.false.,precvec(1:nn(io),nprec(jo),jo,2))
                  else
                     call interp_gp(nn(jo),(1.0-1.0/ritzval(ii,jo))*dxritzvec(1:nn(jo),ii,jo),nn(io),.false.,precvec(1:nn(io),nprec(jo),jo,1))
                     call interp_gp(nn(jo),dxbarritzvec(1:nn(jo),ii,jo),nn(io),.false.,precvec(1:nn(io),nprec(jo),jo,2))
                  end if
               end if
            end do
         end do
      case default
         write(*,'(a)') 'Error: wrong LMP'
         stop
      end select
   end if

   ! Minimization
   call planczosif(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),sigmao,dxbbar(1:nn(io),io),nobs,d,ni,io,nprec(1:io),precvec(1:nn(io),:,1:io,:),dxabar(1:nn(io),io),dx_planczosif(1:nn(io),:,io),jb_planczosif(:,io),jo_planczosif(:,io),ritzval(:,io),dxbarritzvec(1:nn(io),:,io),dxritzvec(1:nn(io),:,io))

   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',jb_planczosif(ii,io)+jo_planczosif(ii,io),' = ',jb_planczosif(ii,io),' + ',jo_planczosif(ii,io)
   end do

   ! New analysis increment
   call apply_b(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),dxabar(1:nn(io),io),dxa(1:nn(io)))

   ! MSE
   xa_mse_planczosif(0,io) = sum(xg(1:nn(io))**2)/real(nn(io))
   do ii=1,ni
      xa_mse_planczosif(ii,io) = sum((xg(1:nn(io))+dx_planczosif(1:nn(io),ii,io))**2)/real(nn(io))
   end do
end do
write(*,*) 

! PLanczos-PLanczosIF comparison
write(*,'(a)') 'PLanczos-PLanczosIF comparison:'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' J=Jb+Jo:',jb_planczos(ii,io)+jo_planczos(ii,io)-(jb_planczosif(ii,io)+jo_planczosif(ii,io)),' = ',jb_planczos(ii,io)-jb_planczosif(ii,io),' + ',jo_planczos(ii,io)-jo_planczosif(ii,io)
   end do
end do
write(*,*) 

! PLanczos-PLanczosIF score
write(*,'(a)') 'PLanczos-PLanczosIF score:'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' MSE planczos / planczosif:',xa_mse_planczos(ii,io),' / ',xa_mse_planczosif(ii,io)
   end do
end do
write(*,*) 

end program main
