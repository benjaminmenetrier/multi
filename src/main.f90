program main

use fft
use interp
use netcdf
use rand
use type_algo
use type_bmatrix
use type_hmatrix
use type_lmp
use type_rmatrix


implicit none

! Parameters
integer             :: n                             ! Full resolution
integer             :: no                            ! Number of outer iterations
integer             :: ni                            ! Number of inner iterations
character(len=1024) :: lmp_mode                      ! LMP mode ('none', 'spectral', 'ritz')
logical             :: new_seed                      ! New random seed
logical             :: full_res                      ! All outer iterations at full resolution if true
! Other parameters used in modules
real(8)             :: sigma_obs                     ! Observation Error standard error
real(8)             :: Lb                            ! Correlation length-scale
real(8)             :: sigmabvar                     ! Grid-point standard deviation variations amplitude

! Local variables
integer                        :: nobs,io,jo,ii,iii
integer,allocatable            :: fac(:),nn(:)
real(8),allocatable            :: xb(:),xg(:),dxb(:),dxbbar(:,:),dxabar(:,:),dxabar_interp(:)
real(8),allocatable            :: vb(:),dvb(:,:),dva(:,:),dva_interp(:)
real(8),allocatable            :: nu(:),yo(:),d(:),hxg(:)
type(bmatrix_type),allocatable :: bmatrix(:)
type(hmatrix_type),allocatable :: hmatrix(:)
type(bmatrix_type)             :: bmatrix_full
type(hmatrix_type)             :: hmatrix_full
type(rmatrix_type)             :: rmatrix
type(algo_type),allocatable    :: algo_lanczos(:),algo_planczosif(:)
type(lmp_type),allocatable     :: lmp_lanczos(:),lmp_planczosif(:)

! /!\ Create later a module that read the inputs and call it here. /!\

! Read the parameters from the standard input
read(*,*) n, no, ni, lmp_mode, sigma_obs, sigmabvar, Lb, full_res, new_seed

! ! Read the parameters from the file 'parameters.dat'
! open(111,file='parameters.dat', action ='read')
! read(111,*) n
! read(111,*) no
! read(111,*) ni
! read(111,*) lmp_mode
! read(111,*) full_res
! read(111,*) new_seed
! ! Other parameters used in modules
! read(111,*) sigma_obs
! read(111,*) sigmabvar
! read(111,*) Lb
! close(111)

! Allocations of tables
allocate(fac(no),nn(no))
allocate(xb(n),xg(n),dxb(n),dxbbar(n,no),dxabar(n,no),dxabar_interp(n))
allocate(vb(n),dvb(n,no),dva(n,no),dva_interp(n))
allocate(algo_lanczos(no),algo_planczosif(no))
allocate(bmatrix(no))
allocate(hmatrix(no))
allocate(lmp_lanczos(no),lmp_planczosif(no))

! Set seed
call set_seed(new_seed)

! FFT test
call fft_test(n)


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
call rmatrix_setup(rmatrix,nobs,sigma_obs)

! Set resolutions
! do io=1,no
!    if (full_res) then
!       fac(io) = 1
!    else
!       fac(io) = 2**(no-io)
!    end if
!    nn(io) = n/fac(io)
! end do


do io=1,no
   if (full_res) then
      fac(io) = 1
   else
      fac(io) = 1.2**(no-io)
   end if
   nn(io) = n/fac(io)
end do


! Setup full resolution B matrix
call bmatrix_setup(bmatrix_full,n,n,sigmabvar,Lb)

! Test full resolution B matrix test
call bmatrix_test(bmatrix_full,n,n/nobs)

do io=1,no
   ! Setup other resolutions B matrix
   call bmatrix_setup(bmatrix(io),nn(io),n,sigmabvar,Lb)
end do

! Test equivalence condition for interpolators
call interp_test(no,nn,bmatrix)

! Allocation
do io=1,no
   call algo_alloc(algo_lanczos(io),nn(io),ni)
   call lmp_alloc(lmp_lanczos(io),nn(io),ni,io,lmp_mode,'control')
   call algo_alloc(algo_planczosif(io),nn(io),ni)
   call lmp_alloc(lmp_planczosif(io),nn(io),ni,io,lmp_mode,'model')
end do
allocate(nu(nobs))
allocate(yo(nobs))
allocate(d(nobs))
allocate(hxg(nobs))

! Setup observations
call rand_normal(nobs,nu)
call rmatrix_apply_sqrt(rmatrix,nobs,nu,yo)

! Generate background state
call rand_normal(n,vb)
call bmatrix_apply_sqrt(bmatrix_full,n,vb,xb)

! Multi-incremental Lanczos in control space

! Result file:
open(42,file='results/lanczos_control_space.dat')
write(42,'(a)') '# Outer iteration , resolution , Inner iteration , J=Jb+Jo , Jb , Jo'

write(*,'(a)') 'Multi-incremental Lanczos in control space'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,', resolution: ',fac(io)

   ! Setup H matrix
   call hmatrix_setup(hmatrix(io),nn(io),nobs)

   ! Compute background increment
   dvb(1:nn(io),io) = 0.0
   
   do jo=1,io-1
      call interp_incr_control(nn(jo),dva(1:nn(jo),jo),nn(io),dva_interp(1:nn(io)))
      dvb(1:nn(io),io) = dvb(1:nn(io),io)-dva_interp(1:nn(io))
   end do

   ! Compute guess
   call bmatrix_apply_sqrt(bmatrix(io),nn(io),dvb(1:nn(io),io),dxb(1:nn(io)))
   xg(1:nn(io)) = xb(1:nn(io))-dxb(1:nn(io))

   ! Compute innovation
   call hmatrix_apply(hmatrix(io),nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   if (io>1) then
      ! Copy preconditioning vectors
      select case (trim(lmp_mode))
      case ('spectral','ritz')
         ! Copy eigenpairs
         lmp_lanczos(io)%outer(io)%eigenval = algo_lanczos(io-1)%eigenval
         lmp_lanczos(io)%outer(io)%eigenvec = algo_lanczos(io-1)%eigenvec

         ! Copy Lanczos vectors
         lmp_lanczos(io)%outer(io)%lancvec_trunc(1:nn(io-1),:) = algo_lanczos(io-1)%lancvec(1:nn(io-1),:)

         ! Compute omega
         lmp_lanczos(io)%outer(io)%omega = algo_lanczos(io-1)%eigenvec(ni,:)*algo_lanczos(io-1)%lastbeta/algo_lanczos(io-1)%eigenval
      end select

      ! Interpolate preconditioning vectors
      call interp_lmp_control(ni,io,nn(1:io),lmp_lanczos(1:io))
   end if

   ! Minimization
   call algo_apply_lanczos(algo_lanczos(io),nn(io),bmatrix(io),hmatrix(io),rmatrix,dvb(1:nn(io),io),nobs,d,ni,lmp_lanczos(io),dva(1:nn(io),io))

   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',algo_lanczos(io)%jb(ii)+algo_lanczos(io)%jo(ii),' = ',algo_lanczos(io)%jb(ii),' + ',algo_lanczos(io)%jo(ii)
      ! Write the results in a file:
      write(42,'(i2,a,i2,a,i3,a,e15.8,a,e15.8,a,e15.8)') io,' ',fac(io),' ',ii,' ',algo_lanczos(io)%jb(ii)+algo_lanczos(io)%jo(ii),' ',algo_lanczos(io)%jb(ii),' ',algo_lanczos(io)%jo(ii)
   end do
end do
write(*,'(a)') '' 
close(42)

! Multi-incremental PLanczosIF in model space

! Result file:
open(43,file='results/PlanczosIF_model_space.dat')
write(43,'(a)') '# Outer iteration , resolution , Inner iteration , J=Jb+Jo , Jb , Jo'

write(*,'(a)') 'Multi-incremental PLanczosIF in model space'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)

   ! Setup H matrix
   call hmatrix_setup(hmatrix(io),nn(io),nobs)

   ! Setup B matrix
   call bmatrix_setup(bmatrix(io),nn(io),n,sigmabvar,Lb)

   ! Compute background increment
   dxbbar(1:nn(io),io) = 0.0
   do jo=1,io-1
      call interp_incr_model(nn(jo),dxabar(1:nn(jo),jo),nn(io),dxabar_interp(1:nn(io)))
      dxbbar(1:nn(io),io) = dxbbar(1:nn(io),io)-dxabar_interp(1:nn(io))
   end do

   ! Compute guess
   call bmatrix_apply(bmatrix(io),nn(io),dxbbar(1:nn(io),io),dxb(1:nn(io)))
   xg(1:nn(io)) = xb(1:nn(io))-dxb(1:nn(io))

   ! Compute innovation
   call hmatrix_apply(hmatrix(io),nn(io),xg(1:nn(io)),nobs,hxg)
   d = yo-hxg

   if (io>1) then
      ! Copy preconditioning vectors
      select case (trim(lmp_mode))
      case ('spectral','ritz')
         ! Copy eigenpairs
         lmp_planczosif(io)%outer(io)%eigenval = algo_planczosif(io-1)%eigenval
         lmp_planczosif(io)%outer(io)%eigenvec = algo_planczosif(io-1)%eigenvec

         ! Copy Lanczos vectors
         lmp_planczosif(io)%outer(io)%lancvec_trunc(1:nn(io-1),:) = algo_planczosif(io-1)%lancvec(1:nn(io-1),:)

         ! Compute omega
         lmp_planczosif(io)%outer(io)%omega = algo_planczosif(io-1)%eigenvec(ni,:)*algo_planczosif(io-1)%lastbeta/algo_planczosif(io-1)%eigenval
      end select

      ! Interpolate preconditioning vectors
      call interp_lmp_model(ni,io,nn(1:io),bmatrix(1:io),lmp_planczosif(1:io))
   end if

   ! Minimization
   call algo_apply_planczosif(algo_planczosif(io),nn(io),bmatrix(io),hmatrix(io),rmatrix,dxbbar(1:nn(io),io),nobs,d,ni,lmp_planczosif(io),dxabar(1:nn(io),io))

   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,', J=Jb+Jo: ',algo_planczosif(io)%jb(ii)+algo_planczosif(io)%jo(ii),' = ',algo_planczosif(io)%jb(ii),' + ',algo_planczosif(io)%jo(ii)
      ! Write the results in a file:
      write(43,'(i2,a,i2,a,i3,a,e15.8,a,e15.8,a,e15.8)') io,' ',fac(io),' ',ii,' ',algo_planczosif(io)%jb(ii)+algo_planczosif(io)%jo(ii),' ',algo_planczosif(io)%jb(ii),' ',algo_planczosif(io)%jo(ii)
   end do
end do
close(43)
write(*,'(a)') '' 

! Lanczos-PLanczosIF comparison

! Result file:
open(44,file='results/lanczos_control_vs_PlanczosIF_model.dat')
write(44,'(a)') '# Outer iteration , resolution , Inner iteration , delta_J , delta_Jb , delta_Jo'

write(*,'(a)') 'Lanczos-PLanczosIF comparison:'
do io=1,no
   write(*,'(a,i2,a,i2)') '   Outer iteration ',io,' resolution: ',fac(io)
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Inner iteration ',ii,' J=Jb+Jo:',algo_lanczos(io)%jb(ii)+algo_lanczos(io)%jo(ii)-(algo_planczosif(io)%jb(ii)+algo_planczosif(io)%jo(ii)),' = ',algo_lanczos(io)%jb(ii)-algo_planczosif(io)%jb(ii),' + ',algo_lanczos(io)%jo(ii)-algo_planczosif(io)%jo(ii)
      ! Write the results in a file:
      write(44,'(i2,a,i2,a,i3,a,e15.8,a,e15.8,a,e15.8)') io,' ',fac(io),' ',ii,' ',algo_lanczos(io)%jb(ii)+algo_lanczos(io)%jo(ii)-(algo_planczosif(io)%jb(ii)+algo_planczosif(io)%jo(ii)),' ',algo_lanczos(io)%jb(ii)-algo_planczosif(io)%jb(ii),' ',algo_lanczos(io)%jo(ii)-algo_planczosif(io)%jo(ii)
   end do
end do
close(44)

end program main
