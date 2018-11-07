!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use fft

implicit none

private
public :: interp_sp,interp_gp,interp_gp_b_inv,interp_gp_b,interp_test

contains

!----------------------------------------------------------------------
! Subroutine: interp_sp
! Purpose: spectral interpolation with zero padding
!----------------------------------------------------------------------
subroutine interp_sp(nntrunc,sptrunc,nn,lnorm,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
logical,intent(in) :: lnorm
real(8),intent(out) :: sp(nn)

! Local variables
real(8) :: norm

if (nntrunc==nn) then
   sp = sptrunc
elseif (nntrunc<nn) then
   ! Initialize
   sp = 0.0

   ! Copy 
   sp(1:nntrunc) = sptrunc

   if (lnorm) then
      ! Normalize
      norm = sqrt(real(nn,8)/real(nntrunc,8))
      sp = sp*norm
   end if

   ! Halve Nyquist frequency
   sp(nntrunc) = 0.5*sp(nntrunc)
end if

end subroutine interp_sp

!----------------------------------------------------------------------
! Subroutine: interp_gp
! Purpose: grid-point interpolation with zero padding in spectral space
!----------------------------------------------------------------------
subroutine interp_gp(nntrunc,gptrunc,nn,lnorm,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
logical,intent(in) :: lnorm
real(8),intent(out) :: gp(nn)

! Local variables
real(8) :: sptrunc(nntrunc),sp(nn)

! FFT
call gp2sp(nntrunc,gptrunc,sptrunc)

! Spectral interpolation
call interp_sp(nntrunc,sptrunc,nn,lnorm,sp)

! Inverse FFT
call sp2gp(nn,sp,gp)

end subroutine interp_gp

!----------------------------------------------------------------------
! Subroutine: interp_gp_b_inv
! Purpose: grid-point interpolation with zero padding in spectral space, specific for B^{-1} space
!----------------------------------------------------------------------
subroutine interp_gp_b_inv(nntrunc,sigmabtrunc,spvartrunc,gptrunc,nn,sigmab,spvar,lnorm,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: lnorm
real(8),intent(out) :: gp(nn)

! Local variables
integer :: i
real(8) :: gptrunctmp(nntrunc)
real(8) :: sptrunc(nntrunc),sp(nn)

! Apply grid-point standard-deviation
gptrunctmp = gptrunc*sigmabtrunc

! FFT
call gp2sp(nntrunc,gptrunctmp,sptrunc)

! Apply spectral standard-deviation
sptrunc = sptrunc*sqrt(spvartrunc)

! Spectral interpolation
call interp_sp(nntrunc,sptrunc,nn,lnorm,sp)

! Apply spectral standard-deviation inverse
do i=1,nn
   if (spvar(i)>tiny(1.0)) then
      sp(i) = sp(i)/sqrt(spvar(i))
   else
      sp(i) = 0.0
   end if
end do

! Inverse FFT
call sp2gp(nn,sp,gp)

! Apply grid-point standard-deviation inverse
gp = gp/sigmab

end subroutine interp_gp_b_inv

!----------------------------------------------------------------------
! Subroutine: interp_gp_b
! Purpose: grid-point interpolation with zero padding in spectral space, specific for B space
!----------------------------------------------------------------------
subroutine interp_gp_b(nntrunc,sigmabtrunc,spvartrunc,gptrunc,nn,sigmab,spvar,lnorm,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: lnorm
real(8),intent(out) :: gp(nn)

! Local variables
integer :: i
real(8) :: gptrunctmp(nntrunc)
real(8) :: sptrunc(nntrunc),sp(nn)

! Apply grid-point standard-deviation inverse
gptrunctmp = gptrunc/sigmabtrunc

! Inverse adjoint FFT
call sp2gp_ad(nntrunc,gptrunctmp,sptrunc)

! Apply spectral standard-deviation inverse
do i=1,nntrunc
   if (spvartrunc(i)>tiny(1.0)) then
      sptrunc(i) = sptrunc(i)/sqrt(spvartrunc(i))
   else
      sptrunc(i) = 0.0
   end if
end do

! Spectral interpolation
call interp_sp(nntrunc,sptrunc,nn,lnorm,sp)

! Apply spectral standard-deviation
sp = sp*sqrt(spvar)

! Adjoint FFT
call gp2sp_ad(nn,sp,gp)

! Apply grid-point standard-deviation
gp = gp*sigmab

end subroutine interp_gp_b

!----------------------------------------------------------------------
! Subroutine: interp_test
! Purpose: test interpolation
!----------------------------------------------------------------------
subroutine interp_test(nn)

implicit none

! Passed variables
integer,intent(in) :: nn

! Local variables
real(8) :: gp(nn),gptrunc(nn/2)
real(8) :: sp(nn),sptrunc(nn/2)

! Interpolation test
call random_number(gptrunc)
call gp2sp(nn/2,gptrunc,sptrunc)
call interp_sp(nn/2,sptrunc,nn,.true.,sp)
call sp2gp(nn,sp,gp)
call sp2gp(nn/2,sptrunc,gptrunc)
write(*,'(a,e15.8)') 'Interpolation test:                      ',maxval(abs(gp(1:nn-1:2)-gptrunc))

end subroutine interp_test

end module interp
