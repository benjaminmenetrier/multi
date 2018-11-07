!----------------------------------------------------------------------
! Module: hmatrix
! Purpose: H matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module hmatrix

implicit none

private
public :: apply_h,apply_h_ad,h_test

contains

!----------------------------------------------------------------------
! Subroutine: apply_h
! Purpose: apply H matrix
!----------------------------------------------------------------------
subroutine apply_h(nn,x,nobs,y)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: x(nn)
integer,intent(in) :: nobs
real(8),intent(out) :: y(nobs)

! Local variables
integer :: iobs,i

! Apply observation operator
do iobs=1,nobs
   i = (iobs-1)*nn/nobs+1
   y(iobs) = x(i)
end do

end subroutine apply_h

!----------------------------------------------------------------------
! Subroutine: apply_h_ad
! Purpose: apply adjoint of the H matrix
!----------------------------------------------------------------------
subroutine apply_h_ad(nobs,y,nn,x)

implicit none

! Passed variables
integer,intent(in) :: nobs
real(8),intent(in) :: y(nobs)
integer,intent(in) :: nn
real(8),intent(out) :: x(nn)

! Local variables
integer :: iobs,i

! Apply observation operator
x = 0.0
do iobs=1,nobs
   i = (iobs-1)*nn/nobs+1
   x(i) = y(iobs)
end do

end subroutine apply_h_ad

!----------------------------------------------------------------------
! Subroutine: h_test
! Purpose: test H matrix
!----------------------------------------------------------------------
subroutine h_test(nn,nobs)

implicit none

! Passed variables
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
call apply_h(nn,x1,nobs,y1)
call apply_h_ad(nobs,y2,nn,x2)
sumx = sum(x1*x2)
sumy = sum(y1*y2)
write(*,'(a,e15.8)') 'Direct H + adjoint H test:               ',sumx-sumy

end subroutine h_test

end module hmatrix
