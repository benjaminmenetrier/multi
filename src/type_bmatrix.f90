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
real(8),parameter :: Lb = 5.0e-2
real(8),parameter :: sigmabvar = 0.0

contains

!----------------------------------------------------------------------
! Subroutine: bmatrix_setup
! Purpose: setup B matrix
!----------------------------------------------------------------------
subroutine bmatrix_setup(bmatrix,nn,nnmax)

implicit none

! Passed variables
type(bmatrix_type),intent(inout) :: bmatrix
integer,intent(in) :: nn
integer,intent(in) :: nnmax

! Local variables
integer :: i
real(8) :: spvar(nnmax),x(nnmax),v(nnmax)

! Check nnmax
if (nn>nnmax) then
   write(*,'(a)') 'ERROR: nn > nnmax in bmatrix_setup'
   stop
end if

! Release memory
if (allocated(bmatrix%sigmab)) deallocate(bmatrix%sigmab)
if (allocated(bmatrix%spvar)) deallocate(bmatrix%spvar)

! Allocation
allocate(bmatrix%sigmab(nn))
allocate(bmatrix%spvar(nn))

! Compute grid-point variance
do i=1,nn
   bmatrix%sigmab(i) = 1.0+sigmabvar*sin(2*pi*real(i-1,8)/real(nn,8))
end do

! Compute spectral variance
spvar(1) = 1.0
do i=2,nnmax/2
   spvar(2*(i-1)) = 2.0*exp(-2.0*(pi*real(i-1,8)*Lb)**2)
   spvar(2*(i-1)+1) = spvar(2*(i-1))
end do
spvar(nnmax) = exp(-2.0*(pi*real(nnmax/2,8)*Lb)**2)

! Compute normalization
x = 0.0
x(1) = 1.0
call gp2sp(nnmax,x,v)
v = v*spvar
call gp2sp_ad(nnmax,v,x)
spvar = spvar/x(1)

! Copy spectral variance
bmatrix%spvar(1:nn) = spvar(1:nn)

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
write(*,'(a)')

end subroutine bmatrix_test

end module type_bmatrix
