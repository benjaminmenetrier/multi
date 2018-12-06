!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use type_bmatrix
use type_lmp
use fft

implicit none

logical,parameter :: sp_from_gp_ref = .false.
logical,parameter :: gp_from_sp_ref = .false.
logical,parameter :: orthotest = .false.

contains

!----------------------------------------------------------------------
! Subroutine: interp_sp
! Purpose: spectral interpolation (zero padding)
!----------------------------------------------------------------------
subroutine interp_sp(nntrunc,sptrunc,nn,sp,sp_from_gp_in)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: sp(nn)
logical,intent(in),optional :: sp_from_gp_in

! Local variables
real(8) :: norm,gptrunc(nntrunc),gp(nn)
logical :: sp_from_gp

if (present(sp_from_gp_in)) then
   sp_from_gp = sp_from_gp_in
else
   sp_from_gp = sp_from_gp_ref
end if

if (sp_from_gp) then
   ! Inverse Fourier transform
   call sp2gp(nntrunc,sptrunc,gptrunc)

   ! Grid-point interpolation
   call interp_gp(nntrunc,gptrunc,nn,gp,.false.)

   ! Direct Fourier transform
   call gp2sp(nn,gp,sp)
else
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
end if

end subroutine interp_sp

!----------------------------------------------------------------------
! Subroutine: interp_gp
! Purpose: grid-point interpolation (linear)
!----------------------------------------------------------------------
subroutine interp_gp(nntrunc,gptrunc,nn,gp,gp_from_sp_in)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(out) :: gp(nn)
logical,intent(in),optional :: gp_from_sp_in

! Local variables
integer :: rat,itrunc,j,i
real(8) :: r,sptrunc(nntrunc),sp(nn)
logical :: gp_from_sp

if (present(gp_from_sp_in)) then
   gp_from_sp = gp_from_sp_in
else
   gp_from_sp = gp_from_sp_ref
end if

if (gp_from_sp) then
   ! Direct Fourier transform
   call gp2sp(nntrunc,gptrunc,sptrunc)

   ! Grid-point interpolation
   call interp_sp(nntrunc,sptrunc,nn,sp,.false.)

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
! Subroutine: interp_lanczos_incr
! Purpose: lanczos interpolation for increment
!----------------------------------------------------------------------
subroutine interp_lanczos_incr(nntrunc,bmatrix_trunc,sptrunc,nn,bmatrix,lanczos_from_planczosif,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
type(bmatrix_type),intent(in) :: bmatrix_trunc
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
type(bmatrix_type),intent(in) :: bmatrix
logical,intent(in) :: lanczos_from_planczosif
real(8),intent(out) :: sp(nn)

! Local variables
real(8) :: gptrunc(nntrunc),gp(nn)
real(8) :: b(nntrunc),guess(nntrunc)

if (lanczos_from_planczosif) then
   ! Interpolation from PLanczosIF
   call bmatrix_apply_sqrt(bmatrix_trunc,nntrunc,sptrunc,b)
   guess = 0.0
   call bmatrix_apply_inv(bmatrix_trunc,nntrunc,b,guess,gptrunc)
   call interp_planczosif_incr(nntrunc,bmatrix_trunc,gptrunc,nn,bmatrix,.false.,gp)
   call bmatrix_apply_sqrt_ad(bmatrix,nn,gp,sp)  
else
   ! Spectral interpolation 
   call interp_sp(nntrunc,sptrunc,nn,sp)
end if

end subroutine interp_lanczos_incr

!----------------------------------------------------------------------
! Subroutine: interp_lanczos_lmp
! Purpose: lanczos interpolation for LMP
!----------------------------------------------------------------------
subroutine interp_lanczos_lmp(ni,io,nn,bmatrix,lmp_lanczos,lmp_planczosif)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(bmatrix_type),intent(in) :: bmatrix(io)
type(lmp_type),intent(inout) :: lmp_lanczos(io)
type(lmp_type),intent(inout),optional :: lmp_planczosif(io)

! Local variables
integer :: jo,ip,jp
real(8) :: proj,norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)
logical :: ortho_loss,norm_loss

do jo=1,io
   ! Copy number of vectors and values
   lmp_lanczos(io)%outer(jo)%np = lmp_lanczos(jo)%outer(jo)%np
   lmp_lanczos(io)%outer(jo)%val(1:lmp_lanczos(io)%outer(jo)%np) = lmp_lanczos(jo)%outer(jo)%val(1:lmp_lanczos(io)%outer(jo)%np)
end do

if (present(lmp_planczosif)) then
   ! Copy number of vectors and values
   lmp_planczosif(io)%outer(io)%np = lmp_lanczos(io)%outer(io)%np
   lmp_planczosif(io)%outer(io)%val = lmp_lanczos(io)%outer(io)%val

   if (io>1) then
      ! Allocation
      allocate(avec(nn(io-1)))
      allocate(bvec(nn(io-1)))
      allocate(cvec(nn(io-1)))
      allocate(dvec(nn(io-1)))
      allocate(evec(nn(io-1)))
      allocate(guess(nn(io-1)))

      ! Transform Ritz vectors from Lanczos space to PLanczosIF space
      do ip=1,lmp_lanczos(io)%outer(io)%np
         ! Initial vector
         avec = lmp_lanczos(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip)

         ! Apply Lanczos LMP square-root
         call lmp_apply_sqrt(lmp_lanczos(io-1),nn(io-1),io-1,avec,bvec)

         ! Apply square-root of B
         call bmatrix_apply_sqrt(bmatrix(io-1),nn(io-1),bvec,dvec)

         ! Apply B inverse
         guess = 0.0
         call bmatrix_apply_inv(bmatrix(io-1),nn(io-1),dvec,guess,cvec)

         ! Apply PLanczosIF LMP inverse
         call lmp_apply_inv(lmp_planczosif(io-1),nn(io-1),io-1,cvec,evec)

         ! Final vector
         lmp_planczosif(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip) = evec
      end do

      ! Release memory
      deallocate(avec)
      deallocate(bvec)
      deallocate(cvec)
      deallocate(dvec)
      deallocate(evec)
      deallocate(guess)
   end if

   ! Interpolate preconditioning vectors
   call interp_planczosif_lmp(ni,io,nn,bmatrix,lmp_planczosif)

   ! Allocation
   allocate(avec(nn(io)))
   allocate(bvec(nn(io)))
   allocate(cvec(nn(io)))
   allocate(evec(nn(io)))

   ! Transform Ritz vectors from PLanczosIF space to Lanczos space
   do jo=1,io
      do ip=1,lmp_lanczos(io)%outer(jo)%np
         ! Initial vector
         evec = lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)

         ! Apply PLanczosIF LMP
         call lmp_apply(lmp_planczosif(io),nn(io),jo-1,evec,cvec)

         ! Apply adjoint square-root of B
         call bmatrix_apply_sqrt_ad(bmatrix(io),nn(io),cvec,bvec)

         ! Apply Lanczos LMP square-root inverse
         call lmp_apply_sqrt_inv(lmp_lanczos(io),nn(io),jo-1,bvec,avec)

         ! Final vectors
         lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip) = avec
         lmp_lanczos(io)%outer(jo)%vec1(1:nn(io),ip) = avec
      end do
   end do
 
   ! Release memory
   deallocate(avec)
   deallocate(bvec)
   deallocate(cvec)
   deallocate(evec)
else
   ! Initialization
   ortho_loss = .false.
   norm_loss = .false.

   do jo=1,io
      ! Interpolation of Ritz vectors
      do ip=1,lmp_lanczos(io)%outer(jo)%np
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%ritzvec_trunc(1:nn(jo-1),ip),nn(io),lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip))
      end do
      
      ! Gram-Schmidt orthogonalization
      do ip=1,lmp_lanczos(io)%outer(jo)%np
         do jp=1,ip-1
            proj = sum(lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),jp)*lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip))/sum(lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),jp)*lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),jp))
            lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip) = lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)-proj*lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),jp)
            if (abs(proj)>1.0e-6) ortho_loss = .true.
         end do
         norm = sqrt(sum(lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)**2))
         lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip) = lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)/norm
         if (abs(norm-1.0)>1.0e-6) norm_loss = .true.
      end do
   end do

   if (orthotest) then
      ! Print orthogonality/norm loss
      if (ortho_loss) write(*,*) '     Orthogonalization loss in Lanczos'
      if (norm_loss) write(*,*) '     Normalization loss in Lanczos'
   end if
end if

if (orthotest) then
   ! Test orthogonality
   do jo=1,io
      do ip=1,lmp_lanczos(io)%outer(jo)%np
         do jp=1,lmp_lanczos(io)%outer(jo)%np
            write(*,*) '     Lanczos ortho test',io,jo,ip,jp,sum(lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)*lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),jp))
         end do
      end do
   end do
end if

do jo=1,io
   ! Copy Ritz vector to vector 1
   do ip=1,lmp_lanczos(io)%outer(jo)%np
      lmp_lanczos(io)%outer(jo)%vec1(1:nn(io),ip) = lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)
   end do
end do

end subroutine interp_lanczos_lmp

!----------------------------------------------------------------------
! Subroutine: interp_planczosif_incr
! Purpose: planczosif interpolation for increment
!----------------------------------------------------------------------
subroutine interp_planczosif_incr(nntrunc,bmatrix_trunc,gptrunc,nn,bmatrix,planczosif_from_lanczos,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
type(bmatrix_type),intent(in) :: bmatrix_trunc
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
type(bmatrix_type),intent(in) :: bmatrix
logical,intent(in) :: planczosif_from_lanczos
real(8),intent(out) :: gp(nn)

! Local variables
real(8) :: sptrunc(nntrunc),sp(nn)
real(8) :: b(nn),guess(nn)

if (planczosif_from_lanczos) then
   ! Interpolation from Lanczos
   call bmatrix_apply_sqrt_ad(bmatrix_trunc,nntrunc,gptrunc,sptrunc)  
   call interp_lanczos_incr(nntrunc,bmatrix_trunc,sptrunc,nn,bmatrix,.false.,sp)
   call bmatrix_apply_sqrt(bmatrix,nn,sp,b)
   call interp_gp(nntrunc,gptrunc,nn,guess)
   call bmatrix_apply_inv(bmatrix,nn,b,guess,gp)  
else
   ! Grid-point interpolation
   call interp_gp(nntrunc,gptrunc,nn,gp)
end if

end subroutine interp_planczosif_incr

!----------------------------------------------------------------------
! Subroutine: interp_planczosif_lmp
! Purpose: planczosif interpolation for LMP
!----------------------------------------------------------------------
subroutine interp_planczosif_lmp(ni,io,nn,bmatrix,lmp_planczosif,lmp_lanczos)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(bmatrix_type),intent(in) :: bmatrix(io)
type(lmp_type),intent(inout) :: lmp_planczosif(io)
type(lmp_type),intent(inout),optional :: lmp_lanczos(io)

! Local variables
integer :: jo,ip,jp
real(8) :: proj,norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)
logical :: ortho_loss,norm_loss

do jo=1,io
   ! Copy number of vectors and values
   lmp_planczosif(io)%outer(jo)%np = lmp_planczosif(jo)%outer(jo)%np
   lmp_planczosif(io)%outer(jo)%val(1:lmp_planczosif(io)%outer(jo)%np) = lmp_planczosif(jo)%outer(jo)%val(1:lmp_planczosif(io)%outer(jo)%np)
end do

if (present(lmp_lanczos)) then
   ! Copy number of vectors and values
   lmp_lanczos(io)%outer(io)%np = lmp_planczosif(io)%outer(io)%np
   lmp_lanczos(io)%outer(io)%val = lmp_planczosif(io)%outer(io)%val

   if (io>1) then
      ! Allocation
      allocate(avec(nn(io-1)))
      allocate(bvec(nn(io-1)))
      allocate(cvec(nn(io-1)))
      allocate(evec(nn(io-1)))

      ! Transform Ritz vectors from PLanczosIF space to Lanczos space
      do ip=1,lmp_planczosif(io)%outer(io)%np
         ! Initial vector
         evec = lmp_planczosif(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip)

         ! Apply PLanczosIF LMP
         call lmp_apply(lmp_planczosif(io-1),nn(io-1),io-1,evec,cvec)

         ! Apply adjoint square-root of B
         call bmatrix_apply_sqrt_ad(bmatrix(io-1),nn(io-1),cvec,bvec)

         ! Apply Lanczos LMP square-root inverse
         call lmp_apply_sqrt_inv(lmp_lanczos(io-1),nn(io-1),io-1,bvec,avec)

         ! Final vector
         lmp_lanczos(io)%outer(io)%ritzvec_trunc(1:nn(io-1),ip) = avec
      end do

      ! Release memory
      deallocate(avec)
      deallocate(bvec)
      deallocate(cvec)
      deallocate(evec)
   end if

   ! Interpolate preconditioning vectors
   call interp_lanczos_lmp(ni,io,nn,bmatrix,lmp_lanczos)

   ! Allocation
   allocate(avec(nn(io)))
   allocate(bvec(nn(io)))
   allocate(cvec(nn(io)))
   allocate(dvec(nn(io)))
   allocate(evec(nn(io)))
   allocate(guess(nn(io)))

   ! Transform Ritz vectors from Lanczos space to PLanczosIF space
   do jo=1,io
      do ip=1,lmp_planczosif(io)%outer(jo)%np
         ! Initial vector
         avec = lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),ip)

         ! Apply Lanczos LMP square-root
         call lmp_apply_sqrt(lmp_lanczos(io),nn(io),jo-1,avec,bvec)

         ! Apply square-root of B
         call bmatrix_apply_sqrt(bmatrix(io),nn(io),bvec,dvec)

         ! Apply B inverse
         guess = 0.0
         call bmatrix_apply_inv(bmatrix(io),nn(io),dvec,guess,cvec)

         ! Apply PLanczosIF LMP inverse
         call lmp_apply_inv(lmp_planczosif(io),nn(io),jo-1,cvec,evec)

         ! Final vectors
         lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip) = evec
         lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),ip) = dvec
         lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip) = cvec
      end do
   end do
 
   ! Release memory
   deallocate(avec)
   deallocate(bvec)
   deallocate(cvec)
   deallocate(dvec)
   deallocate(evec)
   deallocate(guess)
else
   ! Initialization
   ortho_loss = .false.
   norm_loss = .false.

   do jo=1,io
      ! Interpolation of Ritz vectors
      do ip=1,lmp_planczosif(io)%outer(jo)%np
         call interp_gp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%ritzvec_trunc(1:nn(jo-1),ip),nn(io),lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip))
      end do
   
      ! Gram-Schmidt orthogonalization
      do ip=1,lmp_planczosif(io)%outer(jo)%np
         do jp=1,ip-1
            proj = sum(lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),jp)*lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip))/sum(lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),jp)*lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),jp))
            lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip) = lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)-proj*lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),jp)
            if (abs(proj)>1.0e-6) ortho_loss = .true.
         end do
         if (jo==1) then
            lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip) = lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)
         else
            call lmp_apply(lmp_planczosif(io),nn(io),jo-1,lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip),lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip))
         end if
         call bmatrix_apply(bmatrix(io),nn(io),lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip),lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),ip))
         norm = sqrt(sum(lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)*lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),ip)))
         lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip) = lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)/norm
         lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),ip) = lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),ip)/norm
         lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip) = lmp_planczosif(io)%outer(jo)%vec2(1:nn(io),ip)/norm
         if (abs(norm-1.0)>1.0e-6) norm_loss = .true.
      end do
   end do

   if (orthotest) then
      ! Print orthogonality/norm loss
      if (ortho_loss) write(*,*) '     Orthogonalization loss in PLanczosIF'
      if (norm_loss) write(*,*) '     Normalization loss in PLanczosIF'
   end if
end if

if (orthotest) then
   ! Test orthogonality
   do jo=1,io
      do ip=1,lmp_planczosif(io)%outer(jo)%np
         do jp=1,lmp_planczosif(io)%outer(jo)%np
            write(*,*) '     PLanczosIF ortho test',io,jo,ip,jp,sum(lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),ip)*lmp_planczosif(io)%outer(jo)%vec1(1:nn(io),jp))
         end do
      end do
   end do
end if

end subroutine interp_planczosif_lmp

end module interp
