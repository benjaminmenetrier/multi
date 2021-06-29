!----------------------------------------------------------------------
! Module: type_hoperator
! Purpose: H matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_hoperator

use tools_kinds
use tools_netcdf
use type_geom

implicit none

type hoperator_type
   integer :: nobs
   real(kind_real),allocatable :: x_obs(:)
   real(kind_real),allocatable :: y_obs(:)
   real(kind_real),allocatable :: yo(:)
   real(kind_real),allocatable :: Hnl_coeff
contains
   procedure :: setup => hoperator_setup
   procedure :: apply => hoperator_apply
   procedure :: apply_tl => hoperator_apply_tl
   procedure :: apply_ad => hoperator_apply_ad
   procedure :: test => hoperator_test
end type hoperator_type

contains

!----------------------------------------------------------------------
! Subroutine: hoperator_setup
! Purpose: setup H matrix
!----------------------------------------------------------------------
subroutine hoperator_setup(hoperator,nobs,Hnl_coeff)

implicit none

! Passed variables
class(hoperator_type),intent(inout) :: hoperator
integer,intent(in) :: nobs
real(kind_real),intent(in) :: Hnl_coeff

! Copy dimension
hoperator%nobs = nobs
hoperator%Hnl_coeff = Hnl_coeff

! Allocation
allocate(hoperator%x_obs(hoperator%nobs))
allocate(hoperator%y_obs(hoperator%nobs))
allocate(hoperator%yo(hoperator%nobs))

! Grid coordinates
call random_number(hoperator%x_obs)
call random_number(hoperator%y_obs)

end subroutine hoperator_setup

!----------------------------------------------------------------------
! Subroutine: hoperator_apply
! Purpose: apply H operator
!----------------------------------------------------------------------
subroutine hoperator_apply(hoperator,geom,x,y)

implicit none

! Passed variables
class(hoperator_type),intent(in) :: hoperator
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: y(hoperator%nobs)

! Local variables
integer :: iobs,inh
real(kind_real) :: xtmp(geom%nh)

! Apply nonlinear variable transform
do inh=1,geom%nh
   xtmp(inh) = (one-hoperator%Hnl_coeff)*x(inh)+hoperator%Hnl_coeff*x(inh)**3
!   xtmp(inh) = x(inh)+hoperator%Hnl_coeff*sin(x(inh))
end do

! Apply bilinear interpolation
do iobs=1,hoperator%nobs
   call geom%interp_bilinear(hoperator%x_obs(iobs),hoperator%y_obs(iobs),xtmp,y(iobs))
end do

end subroutine hoperator_apply

!----------------------------------------------------------------------
! Subroutine: hoperator_apply_tl
! Purpose: apply H operator, tangent linear
!----------------------------------------------------------------------
subroutine hoperator_apply_tl(hoperator,geom,xg,x,y)

implicit none

! Passed variables
class(hoperator_type),intent(in) :: hoperator
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: xg(geom%nh)
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: y(hoperator%nobs)

! Local variables
integer :: iobs,inh
real(kind_real) :: xtmp(geom%nh)

! Apply linearized variable transform
do inh=1,geom%nh
   xtmp(inh) = ((one-hoperator%Hnl_coeff)+3.0_kind_real*hoperator%Hnl_coeff*xg(inh)**2)*x(inh)
!   xtmp(inh) = (one+hoperator%Hnl_coeff*cos(xg(inh)))*x(inh)
end do

! Apply bilinear interpolation
do iobs=1,hoperator%nobs
   call geom%interp_bilinear(hoperator%x_obs(iobs),hoperator%y_obs(iobs),xtmp,y(iobs))
end do

end subroutine hoperator_apply_tl

!----------------------------------------------------------------------
! Subroutine: hoperator_apply_ad
! Purpose: apply H operator, adjoint
!----------------------------------------------------------------------
subroutine hoperator_apply_ad(hoperator,geom,xg,y,x)

implicit none

! Passed variables
class(hoperator_type),intent(in) :: hoperator
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: xg(geom%nh)
real(kind_real),intent(in) :: y(hoperator%nobs)
real(kind_real),intent(out) :: x(geom%nh)

! Local variables
integer :: iobs,inh
real(kind_real) :: xtmp(geom%nh)

! Apply bilinear interpolation, adjoint
xtmp = zero
do iobs=1,hoperator%nobs
   call geom%interp_bilinear_ad(hoperator%x_obs(iobs),hoperator%y_obs(iobs),y(iobs),xtmp)
end do

! Apply linearized variable transform, adjoint
do inh=1,geom%nh
   x(inh) = ((one-hoperator%Hnl_coeff)+3.0_kind_real*hoperator%Hnl_coeff*xg(inh)**2)*xtmp(inh)
!   x(inh) = (one+hoperator%Hnl_coeff*cos(xg(inh)))*xtmp(inh)
end do

end subroutine hoperator_apply_ad

!----------------------------------------------------------------------
! Subroutine: hoperator_test
! Purpose: setup H matrix
!----------------------------------------------------------------------
subroutine hoperator_test(hoperator,geom)

implicit none

! Passed variables
class(hoperator_type),intent(in) :: hoperator
type(geom_type),intent(in) :: geom

! Local variables
real(kind_real) :: x1(geom%nh),x2(geom%nh)
real(kind_real) :: xg(geom%nh)
real(kind_real) :: y1(hoperator%nobs),y2(hoperator%nobs)
real(kind_real) :: dp1,dp2

! Initialization
call random_number(x1)
call random_number(y2)
call random_number(xg)

! Direct H + adjoint H test
call hoperator%apply_tl(geom,xg,x1,y1)
call hoperator%apply_ad(geom,xg,y2,x2)
dp1 = sum(x1*x2)
dp2 = sum(y1*y2)
write(*,'(a,e15.8)') '         Direct H + adjoint H test: ',two*abs(dp1-dp2)/abs(dp1+dp2)

end subroutine hoperator_test

end module type_hoperator
