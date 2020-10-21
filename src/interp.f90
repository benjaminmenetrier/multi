!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use fft
use type_bmatrix
use type_lmp

implicit none

logical,parameter :: gp_from_sp = .true. ! Use spectral interpolation for model space

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

if (nntrunc==nn) then
   sp = sptrunc
elseif (nntrunc<nn) then
   ! Initialize
   sp = 0.0

   ! Copy
   sp(1:nntrunc) = sptrunc
end if

end subroutine interp_sp

!----------------------------------------------------------------------
! Subroutine: interp_gp
! Purpose: grid-point interpolation (model)
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
real(8) :: r,sptrunc(nntrunc),sp(nn)

if (gp_from_sp) then
   ! Direct Fourier transform
   call gp2sp(nntrunc,gptrunc,sptrunc)

   ! Grid-point interpolation
   call interp_sp(nntrunc,sptrunc,nn,sp)

   ! Inverse Fourier transform
   call sp2gp(nn,sp,gp)
else
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
end if

end subroutine interp_gp

!----------------------------------------------------------------------
! Subroutine: interp_incr_control
! Purpose: increment interpolation in control space
!----------------------------------------------------------------------
subroutine interp_incr_control(nntrunc,sptrunc,nn,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: sp(nn)

! Spectral interpolation
call interp_sp(nntrunc,sptrunc,nn,sp)

end subroutine interp_incr_control

!----------------------------------------------------------------------
! Subroutine: interp_lmp_control
! Purpose: LMP interpolation in control space
!----------------------------------------------------------------------
subroutine interp_lmp_control(ni,io,nn,lmp_lanczos)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(lmp_type),intent(inout) :: lmp_lanczos(io)

! Local variables
integer :: jo,ii,ji
real(8) :: proj,norm

if (io==1) then
   write(*,'(a)') 'ERROR: interpolation of LMP should not be called for io=1'
   stop
end if

select case (trim(lmp_lanczos(1)%mode))
case ('spectral','ritz')
   do jo=2,io
      ! Copy eigenpairs and omega
      lmp_lanczos(io)%outer(jo)%eigenval = lmp_lanczos(jo)%outer(jo)%eigenval
      lmp_lanczos(io)%outer(jo)%eigenvec = lmp_lanczos(jo)%outer(jo)%eigenvec
      lmp_lanczos(io)%outer(jo)%omega = lmp_lanczos(jo)%outer(jo)%omega

      do ii=1,ni+1
         ! Interpolation of Lanczos vectors
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%lancvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii))

         ! Orthogonality test
         do ji=1,ii-1
            proj = sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji)*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii))/sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji)*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji))
            ! if (abs(proj)>1.0e-6) then
            !    write(*,'(a)') 'ERROR: orthogonality lost in interp_lmp_control'
            !    stop
            ! end if
         end do
         norm = sqrt(sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)**2))
         ! if (abs(norm-1.0)>1.0e-6) then
         !    write(*,'(a)') 'ERROR: orthogonality lost in interp_lmp_control'
         !    stop
         ! end if
      end do

      ! Ritz vectors
      lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_lanczos(io)%outer(jo)%eigenvec)
   end do
end select

end subroutine interp_lmp_control

!----------------------------------------------------------------------
! Subroutine: interp_incr_model
! Purpose: increment interpolation in model space
!----------------------------------------------------------------------
subroutine interp_incr_model(nntrunc,gptrunc,nn,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: gp(nn)

! Grid-point interpolation
call interp_gp(nntrunc,gptrunc,nn,gp)

end subroutine interp_incr_model

!----------------------------------------------------------------------
! Subroutine: interp_lmp_model
! Purpose: LMP interpolation in model space
!----------------------------------------------------------------------
subroutine interp_lmp_model(ni,io,nn,bmatrix,lmp_planczosif)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(bmatrix_type),intent(in) :: bmatrix(io)
type(lmp_type),intent(inout) :: lmp_planczosif(io)

! Local variables
integer :: jo,ii,ji
real(8) :: proj,norm

if (io==1) then
   write(*,'(a)') 'ERROR: interpolation of LMP should not be called for io=1'
   stop
end if

select case (trim(lmp_planczosif(1)%mode))
case ('spectral','ritz')
   do jo=2,io
      ! Copy eigenpairs and omega
      lmp_planczosif(io)%outer(jo)%eigenval = lmp_planczosif(jo)%outer(jo)%eigenval
      lmp_planczosif(io)%outer(jo)%eigenvec = lmp_planczosif(jo)%outer(jo)%eigenvec
      lmp_planczosif(io)%outer(jo)%omega = lmp_planczosif(jo)%outer(jo)%omega

      do ii=1,ni+1
         ! Interpolation of Lanczos vectors
         call interp_gp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%lancvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii))

         ! Regeneration of other Lanczos vectors
         call lmp_apply(lmp_planczosif(io),nn(io),ni,jo-1,lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii),lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii))
         call bmatrix_apply(bmatrix(io),nn(io),lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii),lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii))

         ! Orthogonality test
         do ji=1,ii-1
            proj = sum(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ji)*lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii))/sum(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ji)*lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ji))
            ! if (abs(proj)>1.0e-6) then
            !    write(*,'(a)') 'ERROR: orthogonality lost in interp_lmp_model'
            !    stop
            ! end if
         end do
         norm = sqrt(sum(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)*lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii)))
         ! if (abs(norm-1.0)>1.0e-6) then
         !    write(*,'(a)') 'ERROR: orthogonality lost in interp_lmp_model'
         !    stop
         ! end if
      end do

      ! Ritz vectors
      lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
      lmp_planczosif(io)%outer(jo)%ritzvec1(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
      lmp_planczosif(io)%outer(jo)%ritzvec2(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
   end do
end select

end subroutine interp_lmp_model

!----------------------------------------------------------------------
! Subroutine: interp_test
! Purpose: test interpolation for the multi-resolution case
!----------------------------------------------------------------------
subroutine interp_test(no,nn,bmatrix)

implicit none

! Passed variables
integer,intent(in) :: no
integer,intent(in) :: nn(no)
type(bmatrix_type),intent(in) :: bmatrix(no)

! Local variables
integer :: io,jo
real(8),allocatable :: gptrunc(:),gp(:)
real(8),allocatable :: sptrunc(:),sp1(:),sp2(:)

do io=1,no
   do jo=1,io-1
      ! Allocation
      allocate(gptrunc(nn(jo)))
      allocate(gp(nn(io)))
      allocate(sptrunc(nn(jo)))
      allocate(sp1(nn(io)))
      allocate(sp2(nn(io)))

      ! Generate random vector
      call random_number(gptrunc)

      ! Apply adjoint of the square-root of B and then interpolate
      call bmatrix_apply_sqrt_ad(bmatrix(jo),nn(jo),gptrunc,sptrunc)
      call interp_incr_control(nn(jo),sptrunc,nn(io),sp1)

      ! Interpolate and then apply adjoint of the square-root of B
      call interp_incr_model(nn(jo),gptrunc,nn(io),gp)
      call bmatrix_apply_sqrt_ad(bmatrix(io),nn(io),gp,sp2)

      ! Compute difference
      write(*,'(a,i4,a,i4,a,e15.8)') 'Interpolation test for ',nn(jo),' ~> ',nn(io),': ',maxval(abs(sp2-sp1))

      ! Release memory
      deallocate(gptrunc)
      deallocate(gp)
      deallocate(sptrunc)
      deallocate(sp1)
      deallocate(sp2)
   end do
end do

end subroutine interp_test

end module interp
