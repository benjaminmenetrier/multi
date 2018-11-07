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
real(8),parameter :: Lb = 1.0e-1

private
public :: sp_variance,gp_variance,apply_u,apply_u_ad,apply_b,b_test

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
   sigmab(i) = 1.0+0.5*sin(2*pi*real(i-1,8)/real(nn,8))
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
   spvar(2*(i-1)) = 2.0*exp(-2.0*(pi*real(i-1,8)*Lb)**2)
   spvar(2*(i-1)+1) = spvar(2*(i-1))
end do
spvar(nn) = exp(-2.0*(pi*real(nn/2,8)*Lb)**2)

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
write(*,'(a,e15.8)') 'Correlation at obs separation:           ',gp2(1+dobs)
write(*,*) 

end subroutine b_test

end module bmatrix
