!----------------------------------------------------------------------
! Module: bmatrix
! Purpose: B matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module bmatrix

use fft

implicit none

real(8),parameter :: pi = acos(-1.0)
real(8),parameter :: Lb = 0.5e-1
real(8),parameter :: spvarmin = 1.0e-4
logical,parameter :: direct_inverse = .false.

contains

!----------------------------------------------------------------------
! Subroutine: gp_variance
! Purpose: compute grid-point variance
!----------------------------------------------------------------------
subroutine gp_variance(nn,sigmab)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(out) :: sigmab(nn)

! Local variables
integer :: i

! Compute grid-point variance
do i=1,nn
   sigmab(i) = 1.0!+0.5*sin(2*pi*real(i-1,8)/real(nn,8))
end do

end subroutine gp_variance

!----------------------------------------------------------------------
! Subroutine: sp_variance
! Purpose: compute spectral variance
!----------------------------------------------------------------------
subroutine sp_variance(nn,spvar)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(out) :: spvar(nn)

! Local variables
integer :: i
real(8) :: x(nn)
real(8) :: v(nn)

! Compute spectral variance
spvar(1) = 1.0
do i=2,nn/2
   spvar(2*(i-1)) = 2.0*max(exp(-2.0*(pi*real(i-1,8)*Lb)**2),spvarmin)
   spvar(2*(i-1)+1) = spvar(2*(i-1))
end do
spvar(nn) = max(exp(-2.0*(pi*real(nn/2,8)*Lb)**2),spvarmin)

! Compute normalization
x = 0.0
x(1) = 1.0
call gp2sp(nn,x,v)
v = v*spvar
call gp2sp_ad(nn,v,x)
spvar = spvar/x(1)

end subroutine sp_variance

!----------------------------------------------------------------------
! Subroutine: apply_u
! Purpose: apply square-root of the B matrix
!----------------------------------------------------------------------
subroutine apply_u(nn,sigmab,spvar,v,x)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: v(nn)
real(8),intent(out) :: x(nn)

! Local variables
real(8) :: vtmp(nn)

! Apply spectral standard-deviation
vtmp = v*sqrt(spvar)

! Adjoint FFT
call gp2sp_ad(nn,vtmp,x)

! Apply grid-point standard-deviation
x = x*sigmab

end subroutine apply_u

!----------------------------------------------------------------------
! Subroutine: apply_u_ad
! Purpose: apply adjoint of the square-root of the B matrix
!----------------------------------------------------------------------
subroutine apply_u_ad(nn,sigmab,spvar,x,v)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: v(nn)

! Local variables
real(8) :: xtmp(nn)

! Apply grid-point standard-deviation
xtmp = x*sigmab

! FFT
call gp2sp(nn,xtmp,v)

! Apply spectral standard-deviation
v = v*sqrt(spvar)

end subroutine apply_u_ad

!----------------------------------------------------------------------
! Subroutine: apply_b
! Purpose: apply B matrix
!----------------------------------------------------------------------
subroutine apply_b(nn,sigmab,spvar,x,bx)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: bx(nn)

! Local variables
real(8) :: v(nn)

! Apply adjoint of the square-root of the B matrix
call apply_u_ad(nn,sigmab,spvar,x,v)

! Apply square-root of the B matrix
call apply_u(nn,sigmab,spvar,v,bx)

end subroutine apply_b

!----------------------------------------------------------------------
! Subroutine: apply_b_inv_precond
! Purpose: apply B matrix inverse preconditioner
!----------------------------------------------------------------------
subroutine apply_b_inv_precond(nn,sigmab,spvar,bx,x)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: bx(nn)
real(8),intent(out) :: x(nn)

! Local variables
real(8) :: v(nn)

! Apply grid-point standard-deviation inverse
x = bx/sigmab

! Adjoint inverse FFT
call sp2gp_ad(nn,x,v)

! Apply spectral standard-deviation inverse
v = v/spvar

! Inverse FFT
call sp2gp(nn,v,x)

! Apply grid-point standard-deviation inverse
x = x/sigmab

end subroutine apply_b_inv_precond

!----------------------------------------------------------------------
! Subroutine: inverse_b
! Purpose: inverse B matrix
!----------------------------------------------------------------------
subroutine inverse_b(nn,sigmab,spvar,b,guess,x)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: b(nn)
real(8),intent(in) :: guess(nn)
real(8),intent(out) :: x(nn)

! Parameter
integer,parameter :: ni = 50

! Local variables
integer :: ii,ji,ierr,nimax
real(8) :: cost(0:ni),y(ni),rhs(ni),subdiag(ni-1),eigenval(ni),eigenvec(ni,ni),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1)
real(8) :: q(nn,ni),r0(nn),t(nn,0:ni),u(nn,0:ni),v(nn,0:ni+1),w(nn,0:ni),z(nn,0:ni+1),bx(nn)
real(8) :: rmse_rat

! Initialization
v(:,0) = 0.0
u(:,0) = guess
call apply_b(nn,sigmab,spvar,u(:,0),r0)
r0 = b-r0
call apply_b_inv_precond(nn,sigmab,spvar,r0,t(:,0))
beta(0) = sqrt(sum(t(:,0)*r0))
v(:,1) = r0/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)
cost(0) = -0.5*sum(u(:,0)*r0)
nimax = ni

do ii=1,ni
   ! Update
   call apply_b(nn,sigmab,spvar,z(:,ii),q(:,ii))
   q(:,ii) = q(:,ii)-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call apply_b_inv_precond(nn,sigmab,spvar,w(:,ii),t(:,ii))
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
   write(*,*) 'cost',cost(ii)
   if (cost(ii)>cost(ii-1)) then
      nimax = ii-1
      exit
   end if
end do

! Copy final iterate
x = u(:,nimax)

! Check RMSE ratio
call apply_b(nn,sigmab,spvar,x,bx)
rmse_rat = sqrt(sum((b-bx)**2)/sum((b)**2))
if (rmse_rat>1.0e-6) then
   write(*,*) '      B inversion RMSE ratio:',rmse_rat
   stop
end if

end subroutine inverse_b

!----------------------------------------------------------------------
! Subroutine: b_test
! Purpose: test B matrix
!----------------------------------------------------------------------
subroutine b_test(nn,sigmab,spvar,dobs)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
integer,intent(in) :: dobs

! Local variables
integer :: i
real(8) :: gp1(nn),gp2(nn)
real(8) :: sp1(nn),sp2(nn)
real(8) :: sumgp,sumsp
real(8) :: sigmab_one(nn)

! Initialization
call random_number(gp1)
call random_number(gp2)
call gp2sp(nn,gp2,sp2)

! Direct U + adjoint U test
call apply_u_ad(nn,sigmab,spvar,gp1,sp1)
call apply_u(nn,sigmab,spvar,sp2,gp2)
sumgp = sum(gp1*gp2)
sumsp = sum(sp1*sp2)
write(*,'(a,e15.8)') 'Direct U + adjoint U test:               ',sumgp-sumsp

! B normalization test
gp1 = 0.0
gp1(1) = 1.0
sigmab_one = 1.0
call apply_b(nn,sigmab_one,spvar,gp1,gp2)
write(*,'(a,e15.8)') 'B normalization test:                    ',gp2(1)
write(*,'(a,e15.8)') 'Correlation conditioning number:         ',maxval(spvar)/minval(spvar)
write(*,'(a,e15.8)') 'Correlation at obs separation:           ',gp2(dobs)
write(*,'(a)',advance='no') 'Correlation shape:                       '
do i=1,nn/2
   if (abs(gp2(i))<5.0e-3) exit
   write(*,'(f5.2)',advance='no') gp2(i)
end do
write(*,*)

end subroutine b_test

end module bmatrix
