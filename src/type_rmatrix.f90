!----------------------------------------------------------------------
! Module: type_rmatrix
! Purpose: R matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_rmatrix

implicit none

type rmatrix_type
   real(8),allocatable :: sigmao(:)
end type rmatrix_type

real(8),parameter :: sigmao = 0.1 ! Observation error standard-error

contains

!----------------------------------------------------------------------
! Subroutine: rmatrix_setup
! Purpose: setup R matrix
!----------------------------------------------------------------------
subroutine rmatrix_setup(rmatrix,nobs)

implicit none

! Passed variables
type(rmatrix_type),intent(inout) :: rmatrix
integer,intent(in) :: nobs

! Local variables
integer :: iobs

! Release memory
if (allocated(rmatrix%sigmao)) deallocate(rmatrix%sigmao)

! Allocation
allocate(rmatrix%sigmao(nobs))

! Set observation error standard deviation
do iobs=1,nobs
   rmatrix%sigmao(iobs) = sigmao
end do

end subroutine rmatrix_setup

!----------------------------------------------------------------------
! Subroutine: rmatrix_apply_sqrt
! Purpose: apply R matrix square-root
!----------------------------------------------------------------------
subroutine rmatrix_apply_sqrt(rmatrix,nobs,yin,yout)

implicit none

! Passed variables
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: nobs
real(8),intent(in) :: yin(nobs)
real(8),intent(out) :: yout(nobs)

! Multiply by standard-deviation
yout = yin*rmatrix%sigmao

end subroutine rmatrix_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: rmatrix_apply_inv
! Purpose: apply R matrix inverse
!----------------------------------------------------------------------
subroutine rmatrix_apply_inv(rmatrix,nobs,yin,yout)

implicit none

! Passed variables
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: nobs
real(8),intent(in) :: yin(nobs)
real(8),intent(out) :: yout(nobs)

! Divide by variance
yout = yin/rmatrix%sigmao**2

end subroutine rmatrix_apply_inv

end module type_rmatrix
