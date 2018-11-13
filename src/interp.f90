!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use bmatrix
use fft

implicit none

contains

!----------------------------------------------------------------------
! Subroutine: interp_sp
! Purpose: spectral interpolation (zero padding)
!----------------------------------------------------------------------
subroutine interp_sp(nntrunc,sptrunc,nn,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
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

   ! Normalize
   norm = sqrt(real(nn,8)/real(nntrunc,8))
   sp = sp*norm

   ! Halve Nyquist frequency
   sp(nntrunc) = 0.5*sp(nntrunc)
end if

end subroutine interp_sp

!----------------------------------------------------------------------
! Subroutine: interp_gp
! Purpose: grid-point interpolation (linear)
!----------------------------------------------------------------------
subroutine interp_gp(nntrunc,gptrunc,nn,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: gp(nn)

! Local variables
integer :: rat,itrunc,j,i
real(8) :: r

if (nntrunc==nn) then
   gp = gptrunc
elseif (nntrunc<nn) then
   ! Initialize
   rat = nn/nntrunc

   ! Interpolate
   do itrunc=1,nntrunc
      do j=0,rat-1
         i = (itrunc-1)*rat+j+1
         r = real(j,8)/real(rat,8)
         if (itrunc<nntrunc) then
            gp(i) = r*gptrunc(itrunc+1)+(1.0-r)*gptrunc(itrunc)
         else
            gp(i) = r*gptrunc(1)+(1.0-r)*gptrunc(itrunc)
         end if
      end do
   end do
end if

end subroutine interp_gp

!----------------------------------------------------------------------
! Subroutine: interp_sp_from_gp
! Purpose: spectral interpolation from grid-point interpolation
!----------------------------------------------------------------------
subroutine interp_sp_from_gp(nntrunc,sptrunc,nn,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: sp(nn)

! Local variables
real(8) :: gptrunc(nntrunc),gp(nn)

! Inverse FFT
call sp2gp(nn,sptrunc,gptrunc)

! Grid-point interpolation
call interp_gp(nntrunc,gptrunc,nn,gp)

! FFT
call gp2sp(nn,gp,sp)

end subroutine interp_sp_from_gp

!----------------------------------------------------------------------
! Subroutine: interp_gp_from_sp
! Purpose: grid-point interpolation from spectral interpolation
!----------------------------------------------------------------------
subroutine interp_gp_from_sp(nntrunc,gptrunc,nn,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: gp(nn)

! Local variables
real(8) :: sptrunc(nntrunc),sp(nn)

! FFT
call gp2sp(nntrunc,gptrunc,sptrunc)

! Spectral interpolation
call interp_sp(nntrunc,sptrunc,nn,sp)

! Inverse FFT
call sp2gp(nn,sp,gp)

end subroutine interp_gp_from_sp

!----------------------------------------------------------------------
! Subroutine: interp_planczos
! Purpose: PLanczos interpolation
!----------------------------------------------------------------------
subroutine interp_planczos(nntrunc,sigmabtrunc,spvartrunc,sptrunc,nn,sigmab,spvar,planczos_from_planczosif,interp_b,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: planczos_from_planczosif
logical,intent(in) :: interp_b
real(8),intent(out) :: sp(nn)

! Local variables

if (planczos_from_planczosif) then
   ! Grid-point interpolation
   if (interp_b) then
      ! Reference is B space
   else
      ! Reference is inverse B space
   end if
else
   ! Spectral interpolation 
   call interp_sp(nntrunc,sptrunc,nn,sp)
end if

end subroutine interp_planczos

!----------------------------------------------------------------------
! Subroutine: interp_planczosif
! Purpose: PLanczosIF interpolation
!----------------------------------------------------------------------
subroutine interp_planczosif(nntrunc,sigmabtrunc,spvartrunc,gptrunc,nn,sigmab,spvar,planczosif_from_planczos,interp_b,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: planczosif_from_planczos
logical,intent(in) :: interp_b
real(8),intent(out) :: gp(nn)

! Local variables
real(8) :: sptrunc(nntrunc),sp(nn)
real(8) :: bgptrunc(nntrunc),b(nn),guess(nn)

if (planczosif_from_planczos) then
   ! Spectral interpolation
   call apply_u_ad(nntrunc,sigmabtrunc,spvartrunc,gptrunc,sptrunc)  
   call interp_sp(nntrunc,sptrunc,nn,sp)
   call apply_u(nn,sigmab,spvar,sp,b)
   call interp_gp(nntrunc,gptrunc,nn,guess)
   call inverse_b(nn,sigmab,spvar,b,guess,gp)  
else
   ! Grid-point interpolation
   if (interp_b) then
      ! Reference is B space
      call apply_b(nntrunc,sigmabtrunc,spvartrunc,gptrunc,bgptrunc)
      call interp_gp(nntrunc,bgptrunc,nn,b)
      call interp_gp(nntrunc,gptrunc,nn,guess)
      call inverse_b(nn,sigmab,spvar,b,guess,gp)  
   else
      ! Reference is inverse B space
      call interp_gp(nntrunc,gptrunc,nn,gp)
   end if
end if

end subroutine interp_planczosif

end module interp
