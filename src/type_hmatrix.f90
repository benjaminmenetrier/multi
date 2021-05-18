!----------------------------------------------------------------------
! Module: type_hmatrix
! Purpose: H matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_hmatrix

use netcdf
use tools_netcdf
use type_geom

implicit none

type hmatrix_type
   integer             :: nobs
   real(8),allocatable :: x_obs(:)
   real(8),allocatable :: y_obs(:)
   real(8),allocatable :: yo(:)
   real(8),allocatable :: Hnl_coeff
contains
   procedure :: setup => hmatrix_setup
   procedure :: write => hmatrix_write
   procedure :: measure => hmatrix_measure
   procedure :: measure_nl => hmatrix_measure_nl
   procedure :: apply => hmatrix_apply
   procedure :: apply_ad => hmatrix_apply_ad
   procedure :: apply_nl => hmatrix_apply_nl
   procedure :: test => hmatrix_test
end type hmatrix_type

contains

!----------------------------------------------------------------------
! Subroutine: hmatrix_setup
! Purpose: setup H matrix
!----------------------------------------------------------------------
subroutine hmatrix_setup(hmatrix,nobs,Hnl_coeff,x_obs,y_obs)

implicit none

! Passed variables
class(hmatrix_type),intent(inout) :: hmatrix
integer,intent(in)                :: nobs
real(8),intent(in)                :: Hnl_coeff
real(8),intent(in)                :: x_obs(nobs)
real(8),intent(in)                :: y_obs(nobs)

! Copy dimension
hmatrix%nobs = nobs
hmatrix%Hnl_coeff = Hnl_coeff

! Allocation
allocate(hmatrix%x_obs(hmatrix%nobs))
allocate(hmatrix%y_obs(hmatrix%nobs))
allocate(hmatrix%yo(hmatrix%nobs))

! Old method:
! Grid coordinates
!call random_number(hmatrix%x_obs)
!call random_number(hmatrix%y_obs)

hmatrix%x_obs = x_obs
hmatrix%y_obs = y_obs

end subroutine hmatrix_setup

!----------------------------------------------------------------------
! Subroutine: hmatrix_write
! Purpose: write H matrix
!----------------------------------------------------------------------
subroutine hmatrix_write(hmatrix,ncid,x_obs_id,y_obs_id,nobs_id)

implicit none

! Passed variables
class(hmatrix_type),intent(inout) :: hmatrix
integer,intent(in) :: ncid
integer,intent(inout) :: nobs_id,x_obs_id,y_obs_id

! Create dimensions
call ncerr('hmatrix_write',nf90_def_dim(ncid,'nobs',hmatrix%nobs,nobs_id))

! Create variables
call ncerr('hmatrix_write',nf90_def_var(ncid,'x_obs',nf90_double,(/nobs_id/),x_obs_id))
call ncerr('hmatrix_write',nf90_def_var(ncid,'y_obs',nf90_double,(/nobs_id/),y_obs_id))

! Write variables
call ncerr('hmatrix_write',nf90_put_var(ncid,x_obs_id,hmatrix%x_obs))
call ncerr('hmatrix_write',nf90_put_var(ncid,y_obs_id,hmatrix%y_obs))

end subroutine hmatrix_write

!----------------------------------------------------------------------
! Subroutine: hmatrix_measure_nl
! Purpose: non linear measure of the state field.
!----------------------------------------------------------------------
subroutine hmatrix_measure_nl(hmatrix,geom,x,x_nl)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in)     :: geom
real(8),intent(out)            :: x_nl(geom%nh)
real(8),intent(in)             :: x(geom%nh)

! Local variables
integer :: inh

do inh=1,geom%nh
   x_nl(inh) = (1-hmatrix%Hnl_coeff)*x(inh) + hmatrix%Hnl_coeff*x(inh)*x(inh)*x(inh)
end do

end subroutine hmatrix_measure_nl

!----------------------------------------------------------------------
! Subroutine: hmatrix_measure
! Purpose: measure of the state field linearized around the guess
!----------------------------------------------------------------------
subroutine hmatrix_measure(hmatrix,geom,xg,x,x_tl)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in)     :: geom
real(8),intent(out)            :: x_tl(geom%nh)
real(8),intent(in)             :: xg(geom%nh), x(geom%nh)

! Local variables
integer :: inh

do inh=1,geom%nh
   x_tl(inh) = (1-hmatrix%Hnl_coeff)*x(inh) + hmatrix%Hnl_coeff*3*xg(inh)*xg(inh)*x(inh)
end do

end subroutine hmatrix_measure

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply
! Purpose: apply H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply(hmatrix,geom,xg,x,y)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in)     :: geom
real(8),intent(in)             :: x(geom%nh), xg(geom%nh)
real(8),intent(out)            :: y(hmatrix%nobs)

! Local variables
integer :: iobs, inh
real(8) :: x_tl(geom%nh)

! Observation operator linearized around xg
call hmatrix%measure(geom,xg,x,x_tl)

! Apply observation operator
y = 0.0
do iobs=1,hmatrix%nobs
   call geom%interp_bilinear(hmatrix%x_obs(iobs),hmatrix%y_obs(iobs),x_tl,y(iobs))
end do
 
end subroutine hmatrix_apply

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply_ad
! Purpose: apply adjoint of the H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply_ad(hmatrix,geom,xg,y,x_tl)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in) :: geom
real(8),intent(in) :: y(hmatrix%nobs),xg(geom%nh)
real(8),intent(out) :: x_tl(geom%nh)

! Local variables
integer :: iobs, inh
real(8) :: yg(hmatrix%nobs),y_tl(hmatrix%nobs), x(geom%nh)

! Apply observation operator adjoint
x = 0.0
do iobs=1,hmatrix%nobs
   call geom%interp_bilinear_ad(hmatrix%x_obs(iobs),hmatrix%y_obs(iobs),y(iobs),x)
end do

! Observation operator linearized around the guess
call hmatrix%measure(geom,xg,x,x_tl)

end subroutine hmatrix_apply_ad

!----------------------------------------------------------------------
! Subroutine: hmatrix_apply_nl
! Purpose: apply non linear H matrix
!----------------------------------------------------------------------
subroutine hmatrix_apply_nl(hmatrix,geom,x,y)

implicit none

! Passed variables
class(hmatrix_type),intent(in) :: hmatrix
type(geom_type),intent(in)     :: geom
real(8),intent(in)             :: x(geom%nh)
real(8),intent(out)            :: y(hmatrix%nobs)

! Local variables
integer :: iobs, inh
real(8) :: xnl(geom%nh)

! Non linear observation operator (cubic):
call hmatrix%measure_nl(geom,x,xnl)

! Apply observation operator
y = 0.0
do iobs=1,hmatrix%nobs
   call geom%interp_bilinear(hmatrix%x_obs(iobs),hmatrix%y_obs(iobs),xnl,y(iobs))
end do

end subroutine hmatrix_apply_nl

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
real(8) :: xg(geom%nh)
real(8) :: y1(hmatrix%nobs),y2(hmatrix%nobs)
real(8) :: dp1,dp2

! Initialization
call random_number(x1)
call random_number(y2)
call random_number(xg)

! Direct H + adjoint H test
call hmatrix%apply(geom,xg,x1,y1)
call hmatrix%apply_ad(geom,xg,y2,x2)
dp1 = sum(x1*x2)
dp2 = sum(y1*y2)
write(*,'(a,e15.8)') '         Direct H + adjoint H test: ',2.0*abs(dp1-dp2)/abs(dp1+dp2)

end subroutine hmatrix_test

end module type_hmatrix
