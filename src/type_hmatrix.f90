!----------------------------------------------------------------------
! Module: type_hmatrix
! Purpose: H matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_hmatrix

use type_geom

implicit none

type hmatrix_type
   integer             :: nobs
   real(8),allocatable :: x_obs(:)
   real(8),allocatable :: y_obs(:)
   real(8),allocatable :: yo(:)
contains
   procedure :: setup => hmatrix_setup
   procedure :: apply => hmatrix_apply
   procedure :: apply_ad => hmatrix_apply_ad
   procedure :: test => hmatrix_test
end type hmatrix_type

contains

!----------------------------------------------------------------------
! Subroutine: hmatrix_setup
! Purpose: setup observations locations
!----------------------------------------------------------------------
subroutine hmatrix_setup(hmatrix,nobs)

implicit none

! Passed variables
class(hmatrix_type),intent(inout) :: hmatrix
integer,intent(in) :: nobs

! Copy dimension
hmatrix%nobs = nobs

! Allocation
allocate(hmatrix%x_obs(hmatrix%nobs))
allocate(hmatrix%y_obs(hmatrix%nobs))
allocate(hmatrix%yo(hmatrix%nobs))

! Grid coordinates
call random_number(hmatrix%x_obs)
call random_number(hmatrix%y_obs)

end subroutine hmatrix_setup

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply
! Purpose: apply H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply(hmatrix,geom,x,y)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in) :: geom
real(8),intent(in) :: x(geom%nh)
real(8),intent(out) :: y(hmatrix%nobs)

! Local variables
integer :: iobs

! Apply observation operator
do iobs=1,hmatrix%nobs
   call geom%interp_gp(hmatrix%x_obs(iobs),hmatrix%y_obs(iobs),x,y(iobs))
end do

end subroutine hmatrix_apply

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply_ad
! Purpose: apply adjoint of the H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply_ad(hmatrix,geom,y,x)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in) :: geom
real(8),intent(in) :: y(hmatrix%nobs)
real(8),intent(out) :: x(geom%nh)

! Local variables
integer :: iobs

! Apply observation operator adjoint
x = 0.0
do iobs=1,hmatrix%nobs
   call geom%interp_gp_ad(hmatrix%x_obs(iobs),hmatrix%y_obs(iobs),y(iobs),x)
end do

end subroutine hmatrix_apply_ad

!----------------------------------------------------------------------
! Subroutine: hmatrix_test
! Purpose: setup H matrix
!----------------------------------------------------------------------
subroutine hmatrix_test(hmatrix,geom)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in) :: geom

! Local variables
real(8) :: x1(geom%nh),x2(geom%nh)
real(8) :: y1(hmatrix%nobs),y2(hmatrix%nobs)
real(8) :: dp1,dp2

! Initialization
call random_number(x1)
call random_number(y2)

! Direct H + adjoint H test
call hmatrix%apply(geom,x1,y1)
call hmatrix%apply_ad(geom,y2,x2)
dp1 = sum(x1*x2)
dp2 = sum(y1*y2)
write(*,'(a,e15.8)') '         Direct H + adjoint H test: ',2.0*abs(dp1-dp2)/abs(dp1+dp2)

end subroutine hmatrix_test

end module type_hmatrix
