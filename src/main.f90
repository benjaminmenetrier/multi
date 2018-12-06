program main

use algo
use fft
use interp
use type_bmatrix
use type_hmatrix
use type_lmp
use type_rmatrix
use netcdf
use rand

implicit none

! Parameters
integer,parameter :: no = 3                            ! Number of outer iterations
integer,parameter :: ni = 5                            ! Number of inner iterations
integer,parameter :: n = 128                           ! Full resolution
logical,parameter :: full_res = .true.                 ! All outer iterations at full resolution if true
logical,parameter :: lanczos_from_planczosif = .false. ! Lanczos method using PLanczosIF interpolation
logical,parameter :: planczosif_from_lanczos = .false. ! PLanczosIF method using Lanczos interpolation
logical,parameter :: lmp = .true.                      ! Use LMP
real(8),parameter :: ritzvalmin = 1.000001             ! Minimum Ritz value for preconditioning vectors
logical,parameter :: new_seed = .true.                 ! New random seed
integer,parameter :: ntest = 1                         ! Number of tests for MSE

! Local variables
integer :: nobs,io,jo,ip,ii,itest
integer :: fac(no),nn(no)
real(8) :: xb(n),xg(n),dxb(n),dxbbar(n,no),dxabar(n,no),dxabar_interp(n)
real(8) :: vb(n),dvb(n,no),dva(n,no),dva_interp(n)
real(8) :: dx_lanczos(n,ni,no),dx_planczosif(n,ni,no)
real(8) :: jb_lanczos(0:ni,no),jo_lanczos(0:ni,no),jb_planczosif(0:ni,no),jo_planczosif(0:ni,no)
real(8) :: xa_mse_lanczos(0:ni,no,ntest),xa_mse_planczosif(0:ni,no,ntest)
real(8) :: ritzval(ni,no),ritzvec(n,ni,no)
real(8),allocatable :: nu(:),yo(:),d(:),hxg(:)
type(bmatrix_type) :: bmatrix_full,bmatrix(no)
type(hmatrix_type) :: hmatrix_full,hmatrix(no)
type(lmp_type) :: lmp_lanczos(no),lmp_planczosif(no)
type(lmp_type) :: lmp_lanczos_tmp(no),lmp_planczosif_tmp(no)
type(rmatrix_type) :: rmatrix

! Check parameters
if (lanczos_from_planczosif.and.planczosif_from_lanczos) then
   write(*,'(a)') 'lanczos_from_planczosif and planczosif_from_lanczos should not be true together'
end if

! Set seed
call set_seed(new_seed)

! FFT test
call fft_test(n)

! Number of observations
nobs = n/2**(no-1)
write(*,'(a,i4)') 'Number of observations:                 ',nobs
if (n/=nobs*2**(no-1)) then
   write(*,'(a)') 'Error: n should be related to the number of outer iterations'
   stop
end if

! Setup full resolution H matrix
call hmatrix_setup(hmatrix_full,n,nobs)

! Test full resolution H matrix test
call hmatrix_test(hmatrix_full,n,nobs)

! Setup R matrix
call rmatrix_setup(rmatrix,nobs)

! Set resolutions
do io=1,no
   if (full_res) then
      fac(io) = 1
   else
      fac(io) = 2**(no-io)
   end if
   nn(io) = n/fac(io)
end do

! Setup full resolution B matrix
call bmatrix_setup(bmatrix_full,n)

! Test full resolution B matrix test
call bmatrix_test(bmatrix_full,n,n/nobs)

! Allocation
do io=1,no
   call lmp_alloc(lmp_lanczos(io),nn(io),ni,io)
   if (planczosif_from_lanczos) call lmp_alloc(lmp_lanczos_tmp(io),nn(io),ni,io)
   call lmp_alloc(lmp_planczosif(io),nn(io),ni,io)
   if (lanczos_from_planczosif) call lmp_alloc(lmp_planczosif_tmp(io),nn(io),ni,io)
end do
allocate(nu(nobs))
allocate(yo(nobs))
allocate(d(nobs))
allocate(hxg(nobs))

do itest=1,ntest
   ! Set seed
   call set_seed(new_seed)
   
   ! Setup observations
   call rand_normal(nobs,nu)
   call rmatrix_apply_sqrt(rmatrix,nobs,nu,yo)
   
   ! Generate background state
   call rand_normal(n,vb)
   call bmatrix_apply_sqrt(bmatrix_full,n,vb,xb)
   
   ! Multi-incremental Lanczos
   if (ntest==1) write(*,'(a)') 'Multi-incremental Lanczos'
   do io=1,no
      if (ntest==1) write(*,'(a,i2,a,i2)') '   Outer iteration ',io,', resolution: ',fac(io)
   
      ! Setup H matrix
      call hmatrix_setup(hmatrix(io),nn(io),nobs)
   
      ! Setup B matrix
      call bmatrix_setup(bmatrix(io),nn(io))
   
      ! Compute background increment
      dvb(1:nn(io),io) = 0.0
      do jo=1,io-1
         call interp_lanczos_incr(nn(jo),bmatrix(jo),dva(1:nn(jo),jo),nn(io),bmatrix(io),lanczos_from_planczosif,dva_interp(1:nn(io)))
         dvb(1:nn(io),io) = dvb(1:nn(io),io)-dva_interp(1:nn(io))
      end do
   
      ! Compute guess
      call bmatrix_apply_sqrt(bmatrix(io),nn(io),dvb(1:nn(io),io),dxb(1:nn(io)))
      xg(1:nn(io)) = xb(1:nn(io))-dxb(1:nn(io))
   
      ! Compute innovation
      call hmatrix_apply(hmatrix(io),nn(io),xg(1:nn(io)),nobs,hxg)
      d = yo-hxg
   
      ! Compute preconditioning vectors
      if (lmp) then
         if (io==1) then
            ! No preconditioning vectors
            lmp_lanczos(1)%outer(1)%np = 0
         else
            ! Select preconditioning vectors
            ip = 0
            do ii=1,ni
               if (ritzval(ii,io-1)>ritzvalmin) then
                  ! Update
                  ip = ip+1
   
                  ! Vectors
                  lmp_lanczos(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip) = ritzvec(1:nn(io-1),ii,io-1)
   
                  ! Value
                  lmp_lanczos(io)%outer(io)%val(ip) = ritzval(ii,io-1)
               end if
            end do
            lmp_lanczos(io)%outer(io)%np = ip
         end if
   
         ! Interpolate preconditioning vectors
         if (lanczos_from_planczosif) then
            call interp_lanczos_lmp(ni,io,nn(1:io),bmatrix(1:io),lmp_lanczos(1:io),lmp_planczosif_tmp(1:io))
         else
            call interp_lanczos_lmp(ni,io,nn(1:io),bmatrix(1:io),lmp_lanczos(1:io))
         end if
      else
         ! No preconditioning vectors
         do jo=1,io
            lmp_lanczos(io)%outer(jo)%np = 0
         end do
      end if
   
      ! Minimization
      call lanczos(nn(io),bmatrix(io),hmatrix(io),rmatrix,dvb(1:nn(io),io),nobs,d,ni,lmp_lanczos(io),dva(1:nn(io),io),dx_lanczos(1:nn(io),:,io),jb_lanczos(:,io),jo_lanczos(:,io),ritzval(:,io),ritzvec(1:nn(io),:,io))
   
      ! Result
      do ii=0,ni
         if (ntest==1) write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',jb_lanczos(ii,io)+jo_lanczos(ii,io),' = ',jb_lanczos(ii,io),' + ',jo_lanczos(ii,io)
      end do
   
      ! MSE
      xa_mse_lanczos(0,io,itest) = sum(xg(1:nn(io))**2)/real(nn(io))
      do ii=1,ni
         xa_mse_lanczos(ii,io,itest) = sum((xg(1:nn(io))+dx_lanczos(1:nn(io),ii,io))**2)/real(nn(io))
      end do
   end do
   if (ntest==1) write(*,*) 
   
   ! Multi-incremental PLanczosIF
   if (ntest==1) write(*,'(a)') 'Multi-incremental PLanczosIF'
   do io=1,no
      if (ntest==1) write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   
      ! Setup H matrix
      call hmatrix_setup(hmatrix(io),nn(io),nobs)
   
      ! Setup B matrix
      call bmatrix_setup(bmatrix(io),nn(io))
   
      ! Compute background increment
      dxbbar(1:nn(io),io) = 0.0
      do jo=1,io-1
         call interp_planczosif_incr(nn(jo),bmatrix(jo),dxabar(1:nn(jo),jo),nn(io),bmatrix(io),planczosif_from_lanczos,dxabar_interp(1:nn(io)))
         dxbbar(1:nn(io),io) = dxbbar(1:nn(io),io)-dxabar_interp(1:nn(io))
      end do
   
      ! Compute guess
      call bmatrix_apply(bmatrix(io),nn(io),dxbbar(1:nn(io),io),dxb(1:nn(io)))
      xg(1:nn(io)) = xb(1:nn(io))-dxb(1:nn(io))
   
      ! Compute innovation
      call hmatrix_apply(hmatrix(io),nn(io),xg(1:nn(io)),nobs,hxg)
      d = yo-hxg
   
      ! Preconditioning vectors
      if (lmp) then
         if (io==1) then
            ! No preconditioning vectors
            lmp_planczosif(io)%outer(io)%np = 0
         else
            ! Select preconditioning vectors
            ip = 0
            do ii=1,ni
               if (ritzval(ii,io-1)>ritzvalmin) then
                  ! Update
                  ip = ip+1
   
                  ! Vectors
                  lmp_planczosif(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip) = ritzvec(1:nn(io-1),ii,io-1)
   
                  ! Value
                  lmp_planczosif(io)%outer(io)%val(ip) = ritzval(ii,io-1)
               end if
            end do
            lmp_planczosif(io)%outer(io)%np = ip
         end if
   
         ! Interpolate preconditioning vectors
         if (planczosif_from_lanczos) then
            call interp_planczosif_lmp(ni,io,nn(1:io),bmatrix(1:io),lmp_planczosif(1:io),lmp_lanczos_tmp(1:io))
         else
            call interp_planczosif_lmp(ni,io,nn(1:io),bmatrix(1:io),lmp_planczosif(1:io))
         end if
      else
         ! No preconditioning vectors
         do jo=1,io
            lmp_planczosif(io)%outer(jo)%np = 0
         end do
      end if
   
      ! Minimization
      call planczosif(nn(io),bmatrix(io),hmatrix(io),rmatrix,dxbbar(1:nn(io),io),nobs,d,ni,lmp_planczosif(io),dxabar(1:nn(io),io),dx_planczosif(1:nn(io),:,io),jb_planczosif(:,io),jo_planczosif(:,io),ritzval(:,io),ritzvec(1:nn(io),:,io))
   
      ! Result
      do ii=0,ni
         if (ntest==1) write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',jb_planczosif(ii,io)+jo_planczosif(ii,io),' = ',jb_planczosif(ii,io),' + ',jo_planczosif(ii,io)
      end do
   
      ! MSE
      xa_mse_planczosif(0,io,itest) = sum(xg(1:nn(io))**2)/real(nn(io))
      do ii=1,ni
         xa_mse_planczosif(ii,io,itest) = sum((xg(1:nn(io))+dx_planczosif(1:nn(io),ii,io))**2)/real(nn(io))
      end do
   end do
   if (ntest==1) write(*,*) 
   
   if (ntest==1) then
      ! Lanczos-PLanczosIF comparison
      write(*,'(a)') 'Lanczos-PLanczosIF comparison:'
      do io=1,no
         write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
         do ii=0,ni
            write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' J=Jb+Jo:',jb_lanczos(ii,io)+jo_lanczos(ii,io)-(jb_planczosif(ii,io)+jo_planczosif(ii,io)),' = ',jb_lanczos(ii,io)-jb_planczosif(ii,io),' + ',jo_lanczos(ii,io)-jo_planczosif(ii,io)
         end do
      end do
   end if
end do

if (ntest>1) then
   ! Lanczos-PLanczosIF score
   write(*,'(a)') 'Lanczos-PLanczosIF score:'
   do io=1,no
      write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
      do ii=0,ni
         write(*,'(a,i3,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' MSE lanczos / planczosif:',sum(xa_mse_lanczos(ii,io,:))/real(ntest,8),' / ',sum(xa_mse_planczosif(ii,io,:))/real(ntest,8)
      end do
   end do
end if

end program main
