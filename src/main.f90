program main

use algo
use bmatrix
use fft
use hmatrix
use interp
use type_lmp
use netcdf
use rand

implicit none

! Parameters
integer,parameter :: no = 3
integer,parameter :: ni = 5
integer,parameter :: n = 128
real(8),parameter :: sigmao = 0.1
logical :: full_res = .true.
logical :: lanczos_from_planczosif = .true.
logical :: planczosif_from_lanczos = .false.
logical :: no_update = .false.
character(len=1024) :: lmp = 'spectral'
real(8),parameter :: ritzvalmin = 1.000001

! Local variables
integer :: nobs,io,jo,ip,ii
integer :: fac(no),nn(no)
real(8),allocatable :: yo(:),d(:),hxg(:)
real(8),allocatable :: xb(:),xg(:),xa(:),dxa(:),dxbbar(:,:),dxabar(:,:),dxabar_interp(:)
real(8),allocatable :: sigmab(:,:),spvar(:,:)
real(8),allocatable :: jb_lanczos(:,:),jo_lanczos(:,:),jb_planczosif(:,:),jo_planczosif(:,:)
real(8),allocatable :: xa_mse_lanczos(:,:),xa_mse_planczosif(:,:)
real(8),allocatable :: dx_lanczos(:,:,:),dx_planczosif(:,:,:)
real(8),allocatable :: ritzval(:,:),ritzvec(:,:,:)
real(8),allocatable :: vb(:),dvb(:,:),dva(:,:),dva_interp(:)
real(8),allocatable :: xtmp1(:),xtmp2(:)
type(lmp_type) :: lmp_lanczos(no,no),lmp_planczosif(no,no)
type(lmp_type) :: lmp_lanczos_tmp(no,no),lmp_planczosif_tmp(no,no)

! Check parameters
if (lanczos_from_planczosif.and.planczosif_from_lanczos) then
   write(*,'(a)') 'lanczos_from_planczosif and planczosif_from_lanczos should not be true together'
end if

! Set seed
call set_seed(.false.)

! FFT test
call fft_test(n)

! Number of observations
nobs = n/2**(no-1)
write(*,'(a,i4)') 'Number of observations:                 ',nobs
if (n/=nobs*2**(no-1)) then
   write(*,'(a)') 'Error: n should be related to the number of outer iterations'
   stop
end if

! H matrix test
call h_test(n,nobs)

! Set resolutions
do io=1,no
   if (full_res) then
      fac(io) = 1
   else
      fac(io) = 2**(no-io)
   end if
   nn(io) = n/fac(io)
end do

! Allocation
do io=1,no
   do jo=1,io
      call lmp_alloc(lmp_lanczos(jo,io),nn(io),ni)
      if (planczosif_from_lanczos) call lmp_alloc(lmp_lanczos_tmp(jo,io),nn(io),ni) 
      call lmp_alloc(lmp_planczosif(jo,io),nn(io),ni)
      if (lanczos_from_planczosif) call lmp_alloc(lmp_planczosif_tmp(jo,io),nn(io),ni) 
   end do
end do
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
allocate(dx_lanczos(n,ni,no))
allocate(dx_planczosif(n,ni,no))
allocate(jb_lanczos(0:ni,no))
allocate(jo_lanczos(0:ni,no))
allocate(jb_planczosif(0:ni,no))
allocate(jo_planczosif(0:ni,no))
allocate(xa_mse_lanczos(0:ni,no))
allocate(xa_mse_planczosif(0:ni,no))
allocate(ritzval(ni,no))
allocate(ritzvec(n,ni,no))
allocate(vb(n))
allocate(dvb(n,no))
allocate(dva(n,no))
allocate(dva_interp(n))
allocate(xtmp1(n))
allocate(xtmp2(n))

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

! Multi-incremental planczosif
write(*,'(a)') 'Multi-incremental planczosif'
do io=1,no
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
      call interp_gp(nn(io-1),xa(1:nn(io-1)),nn(io),xg(1:nn(io)))
   end if

   ! Innovation
   call apply_h(nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   ! Background increment
   if (.not.no_update) then
      do jo=1,io-1
         call interp_planczosif_incr(nn(jo),sigmab(1:nn(jo),jo),spvar(1:nn(jo),jo),dxabar(1:nn(jo),jo),nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),planczosif_from_lanczos,dxabar_interp(1:nn(io)))
         dxbbar(1:nn(io),io) = dxbbar(1:nn(io),io)-dxabar_interp(1:nn(io))
      end do
   end if

   ! Preconditioning vectors
   select case (trim(lmp))
   case ('none')
      do jo=1,io
         ! No preconditioning vectors
         lmp_planczosif(jo,io)%np = 0
      end do
   case ('spectral')
      if (io==1) then
         ! No preconditioning vectors
         lmp_planczosif(io,io)%np = 0
      else
         ! Select preconditioning vectors
         ip = 0
         do ii=1,ni
            if (ritzval(ii,io-1)>ritzvalmin) then
               ! Update
               ip = ip+1

               ! Vectors
               lmp_planczosif(io,io)%ritzvec_trunc(1:nn(io-1),ip) = ritzvec(1:nn(io-1),ii,io-1)

               ! Value
               lmp_planczosif(io,io)%val(ip) = ritzval(ii,io-1)
            end if
         end do
         lmp_planczosif(io,io)%np = ip
      end if

      ! Interpolate preconditioning vectors
      if (planczosif_from_lanczos) then
         call interp_planczosif_lmp(ni,io,nn(1:io),sigmab(1:nn(io),1:io),spvar(1:nn(io),1:io),lmp_planczosif(1:io,1:io),lmp_lanczos_tmp(1:io,1:io))
      else
         call interp_planczosif_lmp(ni,io,nn(1:io),sigmab(1:nn(io),1:io),spvar(1:nn(io),1:io),lmp_planczosif(1:io,1:io))
      end if
   case default
      write(*,'(a)') 'Error: wrong LMP'
      stop
   end select

   ! Minimization
   call planczosif(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),sigmao,dxbbar(1:nn(io),io),nobs,d,ni,io,lmp_planczosif(1:io,io),dxabar(1:nn(io),io),dx_planczosif(1:nn(io),:,io),jb_planczosif(:,io),jo_planczosif(:,io),ritzval(:,io),ritzvec(1:nn(io),:,io))

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

! Multi-incremental lanczos
write(*,'(a)') 'Multi-incremental lanczos'
do io=1,no
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
      call interp_gp(nn(io-1),xa(1:nn(io-1)),nn(io),xg(1:nn(io)))
   end if

   ! Innovation
   call apply_h(nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   ! Background increment
   if (.not.no_update) then
      do jo=1,io-1
         call interp_lanczos_incr(nn(jo),sigmab(1:nn(jo),jo),spvar(1:nn(jo),jo),dva(1:nn(jo),jo),nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),lanczos_from_planczosif,dva_interp(1:nn(io)))
         dvb(1:nn(io),io) = dvb(1:nn(io),io)-dva_interp(1:nn(io))
      end do
   end if

   ! Preconditioning vectors
   select case (trim(lmp))
   case ('none')
      ! No preconditioning vectors
      do jo=1,io
         lmp_lanczos(jo,io)%np = 0
      end do
   case ('spectral')
      if (io==1) then
         ! No preconditioning vectors
         lmp_lanczos(io,io)%np = 0
      else
         ! Select preconditioning vectors
         ip = 0
         do ii=1,ni
            if (ritzval(ii,io-1)>ritzvalmin) then
               ! Update
               ip = ip+1

               ! Vectors
               lmp_lanczos(io,io)%ritzvec_trunc(1:nn(io-1),ip) = ritzvec(1:nn(io-1),ii,io-1)

               ! Value
               lmp_lanczos(io,io)%val(ip) = ritzval(ii,io-1)
            end if
         end do
         lmp_lanczos(io,io)%np = ip
      end if

      ! Interpolate preconditioning vectors
      if (lanczos_from_planczosif) then
         call interp_lanczos_lmp(ni,io,nn(1:io),sigmab(1:nn(io),1:io),spvar(1:nn(io),1:io),lmp_lanczos(1:io,1:io),lmp_planczosif(1:io,1:io))
      else
         call interp_lanczos_lmp(ni,io,nn(1:io),sigmab(1:nn(io),1:io),spvar(1:nn(io),1:io),lmp_lanczos(1:io,1:io))
      end if
   case default
      write(*,'(a)') 'Error: wrong LMP'
      stop
   end select

   ! Minimization
   call lanczos(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),sigmao,dvb(1:nn(io),io),nobs,d,ni,io,lmp_lanczos(1:io,io),dva(1:nn(io),io),dx_lanczos(1:nn(io),:,io),jb_lanczos(:,io),jo_lanczos(:,io),ritzval(:,io),ritzvec(1:nn(io),:,io))

   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',jb_lanczos(ii,io)+jo_lanczos(ii,io),' = ',jb_lanczos(ii,io),' + ',jo_lanczos(ii,io)
   end do

   ! New analysis increment
   call apply_u(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),dva(1:nn(io),io),dxa(1:nn(io)))

   ! MSE
   xa_mse_lanczos(0,io) = sum(xg(1:nn(io))**2)/real(nn(io))
   do ii=1,ni
      xa_mse_lanczos(ii,io) = sum((xg(1:nn(io))+dx_lanczos(1:nn(io),ii,io))**2)/real(nn(io))
   end do
end do
write(*,*) 



! lanczos-planczosif comparison
write(*,'(a)') 'lanczos-planczosif comparison:'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' J=Jb+Jo:',jb_lanczos(ii,io)+jo_lanczos(ii,io)-(jb_planczosif(ii,io)+jo_planczosif(ii,io)),' = ',jb_lanczos(ii,io)-jb_planczosif(ii,io),' + ',jo_lanczos(ii,io)-jo_planczosif(ii,io)
   end do
end do
write(*,*) 

! lanczos-planczosif score
write(*,'(a)') 'lanczos-planczosif score:'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' MSE lanczos / planczosif:',xa_mse_lanczos(ii,io),' / ',xa_mse_planczosif(ii,io)
   end do
end do
write(*,*) 


end program main
