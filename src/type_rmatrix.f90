!----------------------------------------------------------------------
! Module: type_rmatrix
! Purpose: R matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_rmatrix

use tools_rand

implicit none

type rmatrix_type
   integer :: nobs
   real(8),allocatable :: sigmao(:)
contains
   procedure :: setup => rmatrix_setup
   procedure :: apply_sqrt => rmatrix_apply_sqrt
   procedure :: apply_inv => rmatrix_apply_inv
   procedure :: randomize => rmatrix_randomize
end type rmatrix_type

contains

!----------------------------------------------------------------------
! Subroutine: rmatrix_setup
! Purpose: setup R matrix
!----------------------------------------------------------------------
subroutine rmatrix_setup(rmatrix,nobs,sigma_obs)

implicit none

! Passed variables
class(rmatrix_type),intent(inout) :: rmatrix
integer,intent(in) :: nobs
real(8),intent(in) :: sigma_obs

! Local variables
integer :: iobs

! Copy dimension
rmatrix%nobs = nobs

! Release memory
if (allocated(rmatrix%sigmao)) deallocate(rmatrix%sigmao)

! Allocation
allocate(rmatrix%sigmao(nobs))

! Set observation error standard deviation
do iobs=1,nobs
   rmatrix%sigmao(iobs) = sigma_obs
end do

end subroutine rmatrix_setup

!----------------------------------------------------------------------
! Subroutine: rmatrix_apply_sqrt
! Purpose: apply R matrix square-root
!----------------------------------------------------------------------
subroutine rmatrix_apply_sqrt(rmatrix,yin,yout)

implicit none

! Passed variables
class(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: yin(rmatrix%nobs)
real(8),intent(out) :: yout(rmatrix%nobs)

! Multiply by standard-deviation
yout = yin*rmatrix%sigmao

end subroutine rmatrix_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: rmatrix_apply_inv
! Purpose: apply R matrix inverse
!----------------------------------------------------------------------
subroutine rmatrix_apply_inv(rmatrix,yin,yout)

implicit none

! Passed variables
class(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: yin(rmatrix%nobs)
real(8),intent(out) :: yout(rmatrix%nobs)

! Divide by variance
yout = yin/rmatrix%sigmao**2

end subroutine rmatrix_apply_inv

!----------------------------------------------------------------------
! Subroutine: rmatrix_randomize
! Purpose: randomize the R matrix
!----------------------------------------------------------------------
subroutine rmatrix_randomize(rmatrix,yout)

implicit none

! Passed variables
class(rmatrix_type),intent(in) :: rmatrix
real(8),intent(out) :: yout(rmatrix%nobs)

! Local variable
real(8) :: nu(rmatrix%nobs)

! Gaussian random vector
call rand_normal(rmatrix%nobs,nu)

! Apply R matrix square-root
call rmatrix%apply_sqrt(nu,yout)

end subroutine rmatrix_randomize

end module type_rmatrix
