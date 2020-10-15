!----------------------------------------------------------------------
! Module: fft
! Purpose: FFT-related routines
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module fft

implicit none

INCLUDE 'fftw3.f'

contains
!----------------------------------------------------------------------
! Subroutine: gp2sp
! Purpose: grid-point to spectral
!----------------------------------------------------------------------
subroutine gp2sp(nn,gp,sp)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: gp(nn)
real(8),intent(out) :: sp(nn)

! Local variables
integer :: i
integer(8) :: plan
real(8) :: gptmp(nn),norm
complex(8) :: cp(nn/2+1)

! Setup plan
gptmp = gp
call dfftw_plan_dft_r2c_1d(plan,nn,gptmp,cp,fftw_estimate)

! Apply DFT
gptmp = gp
call dfftw_execute_dft_r2c(plan,gptmp,cp)
cp(2:nn/2) = sqrt(2.0)*cp(2:nn/2)

! Destroy plan
call dfftw_destroy_plan(plan)

! Normalize
norm = 1.0/sqrt(real(nn,8))
cp = cp*norm

! Complex to real
sp(1) = real(cp(1),8)
do i=2,nn/2
   sp(2*(i-1)) = real(cp(i),8)
   sp(2*(i-1)+1) = aimag(cp(i))
end do
sp(nn) = real(cp(nn/2+1),8)

end subroutine gp2sp

!----------------------------------------------------------------------
! Subroutine: sp2gp
! Purpose: spectral to grid-point
!----------------------------------------------------------------------
subroutine sp2gp(nn,sp,gp)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sp(nn)
real(8),intent(out) :: gp(nn)

! Local variables
integer :: i
integer(8) :: plan
real(8) :: norm
complex(8) :: cp(nn/2+1),cptmp(nn/2+1)

! Real to complex
cp(1) = cmplx(sp(1),0.0)
do i=2,nn/2
   cp(i) = cmplx(sp(2*(i-1)),sp(2*(i-1)+1))
end do
cp(nn/2+1) = cmplx(sp(nn),0.0)

! Setup plan
cptmp = cp
call dfftw_plan_dft_c2r_1d(plan,nn,cptmp,gp,fftw_estimate)

! Apply DFT
cptmp = cp
cptmp(2:nn/2) = cptmp(2:nn/2)/sqrt(2.0)
call dfftw_execute_dft_c2r(plan,cptmp,gp)

! Destroy plan
call dfftw_destroy_plan(plan)

! Normalize
norm = 1.0/sqrt(real(nn,8))
gp = gp*norm

end subroutine sp2gp

!----------------------------------------------------------------------
! Subroutine: fft_test
! Purpose: test FFT direct, inverse and adjoint
!----------------------------------------------------------------------
subroutine fft_test(nn)

implicit none

! Passed variables
integer,intent(in) :: nn

! Local variables
real(8) :: gpsave(nn),gp(nn),gp1(nn),gp2(nn),sumgp
real(8) :: spsave(nn),sp(nn),sp1(nn),sp2(nn),sumsp

! Initialization
call random_number(gpsave)
call gp2sp(nn,gpsave,spsave)

! Direct + inverse test
gp = gpsave
call gp2sp(nn,gp,sp)
call sp2gp(nn,sp,gp)
write(*,'(a,e15.8)') 'Direct + inverse test:           ',maxval(abs(gp-gpsave))

! Inverse + direct test
sp = spsave
call sp2gp(nn,sp,gp)
call gp2sp(nn,gp,sp)
write(*,'(a,e15.8)') 'Inverse + direct test:           ',maxval(abs(sp-spsave))

! Adjoint test
call random_number(gp1)
call random_number(gp2)
call gp2sp(nn,gp2,sp2)
call gp2sp(nn,gp1,sp1)
call sp2gp(nn,sp2,gp2)
sumgp = sum(gp1*gp2)
sumsp = sum(sp1*sp2)
write(*,'(a,e15.8)') 'Adjoint test:                    ',sumgp-sumsp

end subroutine fft_test

end module fft
