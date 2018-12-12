!----------------------------------------------------------------------
! Module: type_bmatrix
! Purpose: B matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_bmatrix

use fft

implicit none

type bmatrix_type
   real(8),allocatable :: sigmab(:)
   real(8),allocatable :: spvar(:)
end type bmatrix_type

real(8),parameter :: pi = acos(-1.0)
real(8),parameter :: Lb = 0.5e-1
real(8),parameter :: spvarmin = 1.0e-4
real(8),parameter :: tol = 1.0e-5
integer,parameter :: niter = 50

contains

!----------------------------------------------------------------------
! Subroutine: bmatrix_setup
! Purpose: setup B matrix
!----------------------------------------------------------------------
subroutine bmatrix_setup(bmatrix,nn)

implicit none

! Passed variables
type(bmatrix_type),intent(inout) :: bmatrix
integer,intent(in) :: nn

! Local variables
integer :: i
real(8) :: x(nn)
real(8) :: v(nn)

! Release memory
if (allocated(bmatrix%sigmab)) deallocate(bmatrix%sigmab)
if (allocated(bmatrix%spvar)) deallocate(bmatrix%spvar)

! Allocation
allocate(bmatrix%sigmab(nn))
allocate(bmatrix%spvar(nn))

! Compute grid-point variance
do i=1,nn
   bmatrix%sigmab(i) = 1.0+0.5*sin(2*pi*real(i-1,8)/real(nn,8))
end do

! Compute spectral variance
bmatrix%spvar(1) = 1.0
do i=2,nn/2
   bmatrix%spvar(2*(i-1)) = 2.0*max(exp(-2.0*(pi*real(i-1,8)*Lb)**2),spvarmin)
   bmatrix%spvar(2*(i-1)+1) = bmatrix%spvar(2*(i-1))
end do
bmatrix%spvar(nn) = max(exp(-2.0*(pi*real(nn/2,8)*Lb)**2),spvarmin)

! Compute normalization
x = 0.0
x(1) = 1.0
call gp2sp(nn,x,v)
v = v*bmatrix%spvar
call gp2sp_ad(nn,v,x)
bmatrix%spvar = bmatrix%spvar/x(1)

end subroutine bmatrix_setup

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt
! Purpose: apply square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt(bmatrix,nn,v,x)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
real(8),intent(in) :: v(nn)
real(8),intent(out) :: x(nn)

! Local variables
real(8) :: vtmp(nn)

! Apply spectral standard-deviation
vtmp = v*sqrt(bmatrix%spvar)

! Adjoint FFT
call gp2sp_ad(nn,vtmp,x)

! Apply grid-point standard-deviation
x = x*bmatrix%sigmab

end subroutine bmatrix_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt_ad
! Purpose: apply adjoint of the square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt_ad(bmatrix,nn,x,v)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
real(8),intent(in) :: x(nn)
real(8),intent(out) :: v(nn)

! Local variables
real(8) :: xtmp(nn)

! Apply grid-point standard-deviation
xtmp = x*bmatrix%sigmab

! FFT
call gp2sp(nn,xtmp,v)

! Apply spectral standard-deviation
v = v*sqrt(bmatrix%spvar)

end subroutine bmatrix_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply
! Purpose: bmatrix_apply B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply(bmatrix,nn,x,bx)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
real(8),intent(in) :: x(nn)
real(8),intent(out) :: bx(nn)

! Local variables
real(8) :: v(nn)

! Apply adjoint of the square-root of the B matrix
call bmatrix_apply_sqrt_ad(bmatrix,nn,x,v)

! Apply square-root of the B matrix
call bmatrix_apply_sqrt(bmatrix,nn,v,bx)

end subroutine bmatrix_apply

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_inv_precond
! Purpose: apply B matrix inverse preconditioner
!----------------------------------------------------------------------
subroutine bmatrix_apply_inv_precond(bmatrix,nn,bx,x)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
real(8),intent(in) :: bx(nn)
real(8),intent(out) :: x(nn)

! Local variables
real(8) :: v(nn)

! Apply grid-point standard-deviation inverse
x = bx/bmatrix%sigmab

! Adjoint inverse FFT
call sp2gp_ad(nn,x,v)

! Apply spectral standard-deviation inverse
v = v/bmatrix%spvar

! Inverse FFT
call sp2gp(nn,v,x)

! Apply grid-point standard-deviation inverse
x = x/bmatrix%sigmab

end subroutine bmatrix_apply_inv_precond

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_inv
! Purpose: inverse B matrix with a PLanczos algorithm
!----------------------------------------------------------------------
subroutine bmatrix_apply_inv(bmatrix,nn,b,guess,x)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
real(8),intent(in) :: b(nn)
real(8),intent(in) :: guess(nn)
real(8),intent(out) :: x(nn)

! Local variables
integer :: ii,ierr,nitermax
real(8) :: cost(0:niter),y(niter),rhs(niter),subdiag(niter-1),eigenval(niter),eigenvec(niter,niter),work(max(1,2*niter-2))
real(8) :: alpha(0:niter),beta(0:niter+1)
real(8) :: q(nn,niter),r0(nn),t(nn,0:niter),u(nn,0:niter),v(nn,0:niter+1),w(nn,0:niter),z(nn,0:niter+1),bx(nn)
real(8) :: rmse_rat

! Initialization
v(:,0) = 0.0
u(:,0) = guess
call bmatrix_apply(bmatrix,nn,u(:,0),r0)
r0 = b-r0
if (any(abs(r0)>0.0)) then
   call bmatrix_apply_inv_precond(bmatrix,nn,r0,t(:,0))
   beta(0) = sqrt(sum(t(:,0)*r0))
   v(:,1) = r0/beta(0)
   z(:,1) = t(:,0)/beta(0)
   beta(1) = 0.0
   rhs = 0.0
   rhs(1) = beta(0)
   cost(0) = -0.5*sum(u(:,0)*r0)
   nitermax = niter

   do ii=1,niter
      ! Update
      call bmatrix_apply(bmatrix,nn,z(:,ii),q(:,ii))
      q(:,ii) = q(:,ii)-beta(ii)*v(:,ii-1)
      alpha(ii) = sum(q(:,ii)*z(:,ii))
      w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
      call bmatrix_apply_inv_precond(bmatrix,nn,w(:,ii),t(:,ii))
      beta(ii+1) = sqrt(sum(t(:,ii)*w(:,ii)))
      v(:,ii+1) = w(:,ii)/beta(ii+1)
      z(:,ii+1) = t(:,ii)/beta(ii+1)

      ! Compute eigenpairs
      if (ii==1) then
         eigenval(1) = alpha(1)
         eigenvec(1,1) = 1.0
      else
         eigenval(1:ii) = alpha(1:ii)
         subdiag(1:ii-1) = beta(2:ii)
         call dsteqr('I',ii,eigenval(1:ii),subdiag(1:ii-1),eigenvec(1:ii,1:ii),ii,work(1:2*ii-2),ierr)
         if (ierr/=0) then
            write(*,'(a)') 'Error in dsteqr'
            stop
         end if
         eigenval(1:ii) = max(eigenval(1:ii),1.0)
      end if

      ! Update increment
      y(1:ii) = matmul(eigenvec(1:ii,1:ii),matmul(transpose(eigenvec(1:ii,1:ii)),rhs(1:ii))/eigenval(1:ii))
      u(:,ii) = u(:,0)+matmul(z(:,1:ii),y(1:ii))

      ! Cost function
      cost(ii) = -0.5*sum(u(:,ii)*r0)
      if (cost(ii)>cost(ii-1)) then
         nitermax = ii-1
         exit
      end if
   end do

   ! Copy final iterate
   x = u(:,nitermax)
else
   ! Copy guess
   x = guess
end if

! Check inversion accuracy
call bmatrix_apply(bmatrix,nn,x,bx)
rmse_rat = sqrt(sum((b-bx)**2)/sum(b**2))
if (rmse_rat>tol) then
   write(*,*) '      B inversion RMSE ratio:',rmse_rat
   stop
end if

end subroutine bmatrix_apply_inv

!----------------------------------------------------------------------
! Subroutine: bmatrix_test
! Purpose: test B matrix
!----------------------------------------------------------------------
subroutine bmatrix_test(bmatrix,nn,dobs)

implicit none

! Passed variables
type(bmatrix_type),intent(in) :: bmatrix
integer,intent(in) :: nn
integer,intent(in) :: dobs

! Local variables
integer :: i
real(8) :: gp1(nn),gp2(nn)
real(8) :: sp1(nn),sp2(nn)
real(8) :: sumgp,sumsp

! Initialization
call random_number(gp1)
call random_number(gp2)
call gp2sp(nn,gp2,sp2)

! Direct U + adjoint U test
call bmatrix_apply_sqrt_ad(bmatrix,nn,gp1,sp1)
call bmatrix_apply_sqrt(bmatrix,nn,sp2,gp2)
sumgp = sum(gp1*gp2)
sumsp = sum(sp1*sp2)
write(*,'(a,e15.8)') 'Direct U + adjoint U test:               ',sumgp-sumsp

! Print other parameters
write(*,'(a,e15.8)') 'Correlation conditioning number:         ',maxval(bmatrix%spvar)/minval(bmatrix%spvar)
write(*,'(a,e15.8)') 'Correlation at obs separation:           ',gp2(dobs)
write(*,'(a)',advance='no') 'Correlation shape:                       '
gp1 = 0.0
gp1(1) = 1.0
call bmatrix_apply(bmatrix,nn,gp1,gp2)
gp2 = gp2/(bmatrix%sigmab(1)*bmatrix%sigmab)
do i=1,nn/2
   if (abs(gp2(i))<5.0e-3) exit
   write(*,'(f5.2)',advance='no') gp2(i)
end do
write(*,*)

end subroutine bmatrix_test

end module type_bmatrix
