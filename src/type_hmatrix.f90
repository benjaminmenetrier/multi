!----------------------------------------------------------------------
! Module: type_hmatrix
! Purpose: H matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_hmatrix

implicit none

type hmatrix_type
   integer,allocatable :: ind(:)
end type hmatrix_type

contains

!----------------------------------------------------------------------
! Subroutine: hmatrix_setup
! Purpose: setup H matrix
!----------------------------------------------------------------------
subroutine hmatrix_setup(hmatrix,nn,nobs)

implicit none

! Passed variables
type(hmatrix_type),intent(inout) :: hmatrix
integer,intent(in) :: nn
integer,intent(in) :: nobs

! Local variables
integer :: iobs

! Release memory
if (allocated(hmatrix%ind)) deallocate(hmatrix%ind)

! Allocation
allocate(hmatrix%ind(nobs))

! Set indices
do iobs=1,nobs
   hmatrix%ind(iobs) = (iobs-1)*nn/nobs+1
end do

end subroutine hmatrix_setup

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply
! Purpose: apply H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply(hmatrix,nn,x,nobs,y)

implicit none

! Passed variables
type(hmatrix_type),intent(in) :: hmatrix
integer,intent(in) :: nn
real(8),intent(in) :: x(nn)
integer,intent(in) :: nobs
real(8),intent(out) :: y(nobs)

! Local variables
integer :: iobs

! Apply observation operator
do iobs=1,nobs
   y(iobs) = x(hmatrix%ind(iobs))
end do

end subroutine hmatrix_apply

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply_ad
! Purpose: apply adjoint of the H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply_ad(hmatrix,nobs,y,nn,x)

implicit none

! Passed variables
type(hmatrix_type),intent(in) :: hmatrix
integer,intent(in) :: nobs
real(8),intent(in) :: y(nobs)
integer,intent(in) :: nn
real(8),intent(out) :: x(nn)

! Local variables
integer :: iobs

! Apply observation operator
x = 0.0
do iobs=1,nobs
   x(hmatrix%ind(iobs)) = y(iobs)
end do

end subroutine hmatrix_apply_ad

!----------------------------------------------------------------------
! Subroutine: hmatrix_test
! Purpose: test H matrix
!----------------------------------------------------------------------
subroutine hmatrix_test(hmatrix,nn,nobs)

implicit none

! Passed variables
type(hmatrix_type),intent(in) :: hmatrix
integer,intent(in) :: nn
integer,intent(in) :: nobs

! Local variables
real(8) :: x1(nn),x2(nn)
real(8) :: y1(nobs),y2(nobs)
real(8) :: sumx,sumy

! Initialization
call random_number(x1)
call random_number(y2)

! Direct H + adjoint H test
call hmatrix_apply(hmatrix,nn,x1,nobs,y1)
call hmatrix_apply_ad(hmatrix,nobs,y2,nn,x2)
sumx = sum(x1*x2)
sumy = sum(y1*y2)
write(*,'(a,e15.8)') 'Direct H + adjoint H test:               ',sumx-sumy

end subroutine hmatrix_test

end module type_hmatrix
