!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use bmatrix
use type_lmp
use fft

implicit none

logical,parameter :: orthotest = .false.

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
! Subroutine: interp_lanczos_incr
! Purpose: lanczos interpolation for increment
!----------------------------------------------------------------------
subroutine interp_lanczos_incr(nntrunc,sigmabtrunc,spvartrunc,sptrunc,nn,sigmab,spvar,lanczos_from_planczosif,sp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: sptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: lanczos_from_planczosif
real(8),intent(out) :: sp(nn)

! Local variables
real(8) :: gptrunc(nntrunc),gp(nn)
real(8) :: b(nntrunc),guess(nntrunc)

if (lanczos_from_planczosif) then
   ! Grid-point interpolation
   call apply_u(nntrunc,sigmabtrunc,spvartrunc,sptrunc,b)
   guess = 0.0
   call inverse_b(nntrunc,sigmabtrunc,spvartrunc,b,guess,gptrunc)
   call interp_planczosif_incr(nntrunc,sigmabtrunc,spvartrunc,gptrunc,nn,sigmab,spvar,.false.,gp)
   call apply_u_ad(nn,sigmab,spvar,gp,sp)  
else
   ! Spectral interpolation 
   call interp_sp(nntrunc,sptrunc,nn,sp)
end if

end subroutine interp_lanczos_incr

!----------------------------------------------------------------------
! Subroutine: interp_lanczos_lmp
! Purpose: lanczos interpolation for LMP
!----------------------------------------------------------------------
subroutine interp_lanczos_lmp(ni,io,nn,sigmab,spvar,lmp_lanczos,lmp_planczosif)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
real(8),intent(in) :: sigmab(nn(io),io)
real(8),intent(in) :: spvar(nn(io),io)
type(lmp_type),intent(inout) :: lmp_lanczos(io,io)
type(lmp_type),intent(inout),optional :: lmp_planczosif(io,io)

! Local variables
integer :: jo,ip,jp
real(8) :: norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)

do jo=1,io
   ! Copy number of vectors and values
   lmp_lanczos(jo,io)%np = lmp_lanczos(jo,jo)%np
   lmp_lanczos(jo,io)%val(1:lmp_lanczos(jo,io)%np) = lmp_lanczos(jo,jo)%val(1:lmp_lanczos(jo,io)%np)
end do

if (present(lmp_planczosif)) then
   ! Copy number of vectors and values
   lmp_planczosif(io,io)%np = lmp_lanczos(io,io)%np
   lmp_planczosif(io,io)%val = lmp_lanczos(io,io)%val

   if (io>1) then
      ! Allocation
      allocate(avec(nn(io-1)))
      allocate(bvec(nn(io-1)))
      allocate(cvec(nn(io-1)))
      allocate(dvec(nn(io-1)))
      allocate(evec(nn(io-1)))
      allocate(guess(nn(io-1)))

      ! Transform Ritz vectors from Lanczos space to PLanczosIF space
      do ip=1,lmp_lanczos(io,io)%np
         ! Initial vector
         avec = lmp_lanczos(io,io)%ritzvec_trunc(1:nn(io-1),ip)

         ! Apply Lanczos LMP square-root
         call lmp_apply_sqrt(nn(io-1),io-1,lmp_lanczos(1:io-1,io-1),avec,bvec)

         ! Apply square-root of B
         call apply_u(nn(io-1),sigmab(1:nn(io-1),io-1),spvar(1:nn(io-1),io-1),bvec,dvec)

         ! Apply B inverse
         guess = 0.0
         call inverse_b(nn(io-1),sigmab(1:nn(io-1),io-1),spvar(1:nn(io-1),io-1),dvec,guess,cvec)

         ! Apply PLanczosIF LMP inverse
         call lmp_apply_inv(nn(io-1),io-1,lmp_planczosif(1:io-1,io-1),cvec,evec)

         ! Final vector
Write(*,*) 'test',lmp_planczosif(io,io)%ritzvec_trunc(1:5,ip)
Write(*,*) 'tost',evec(1:5)
         lmp_planczosif(io,io)%ritzvec_trunc(1:nn(io-1),ip) = evec
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
   call interp_planczosif_lmp(ni,io,nn,sigmab,spvar,lmp_planczosif)

   ! Allocation
   allocate(avec(nn(io)))
   allocate(bvec(nn(io)))
   allocate(cvec(nn(io)))
   allocate(evec(nn(io)))

   ! Transform Ritz vectors from PLanczosIF space to Lanczos space
   do jo=1,io
      do ip=1,lmp_lanczos(jo,io)%np
         ! Initial vector
         evec = lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)

         ! Apply PLanczosIF LMP
         call lmp_apply(nn(io),jo-1,lmp_planczosif(1:jo-1,io),evec,cvec)

         ! Apply adjoint square-root of B
         call apply_u_ad(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),cvec,bvec)

         ! Apply Lanczos LMP square-root inverse
         call lmp_apply_sqrt_inv(nn(io),jo-1,lmp_lanczos(1:jo-1,io),bvec,avec)

         ! Final vector
         lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip) = avec
      end do
   end do
 
   ! Release memory
   deallocate(avec)
   deallocate(bvec)
   deallocate(cvec)
   deallocate(evec)
else
   allocate(evec(nn(io)))
   allocate(dvec(nn(io)))
   allocate(cvec(nn(io)))
   allocate(bvec(nn(io)))
   allocate(avec(nn(io)))
   allocate(guess(nn(io)))

   do jo=1,io
      ! Interpolation of Ritz vectors
      do ip=1,lmp_lanczos(jo,io)%np
         call interp_sp(nn(jo-1),lmp_lanczos(jo,jo)%ritzvec_trunc(1:nn(jo-1),ip),nn(io),lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip))
      end do
      
      ! Gram-Schmidt orthogonalization
      do ip=1,lmp_lanczos(jo,io)%np
         do jp=1,ip-1
            lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip) = lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)-sum(lmp_lanczos(jo,io)%ritzvec(1:nn(io),jp)*lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip))/sum(lmp_lanczos(jo,io)%ritzvec(1:nn(io),jp)*lmp_lanczos(jo,io)%ritzvec(1:nn(io),jp))*lmp_lanczos(jo,io)%ritzvec(1:nn(io),jp)
         end do
         norm = sqrt(sum(lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)**2))
         lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip) = lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)/norm
      end do
      
      if (orthotest) then
         ! Test orthogonality
         do ip=1,lmp_lanczos(jo,io)%np
            do jp=1,lmp_lanczos(jo,io)%np
               write(*,*) 'Ortho test',io,jo,ip,jp,sum(lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)*lmp_lanczos(jo,io)%ritzvec(1:nn(io),jp))
            end do
         end do
      end if
   end do
end if

do jo=1,io
   ! Copy Ritz vector to vector 1
   do ip=1,lmp_lanczos(jo,io)%np
      lmp_lanczos(jo,io)%vec1(1:nn(io),ip) = lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)
   end do
end do

end subroutine interp_lanczos_lmp

!----------------------------------------------------------------------
! Subroutine: interp_planczosif_incr
! Purpose: planczosif interpolation for increment
!----------------------------------------------------------------------
subroutine interp_planczosif_incr(nntrunc,sigmabtrunc,spvartrunc,gptrunc,nn,sigmab,spvar,planczosif_from_lanczos,gp)

implicit none

! Passed variables
integer,intent(in) :: nntrunc
real(8),intent(in) :: sigmabtrunc(nntrunc)
real(8),intent(in) :: spvartrunc(nntrunc)
real(8),intent(in) :: gptrunc(nntrunc)
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
logical,intent(in) :: planczosif_from_lanczos
real(8),intent(out) :: gp(nn)

! Local variables
real(8) :: sptrunc(nntrunc),sp(nn)
real(8) :: b(nn),guess(nn)

if (planczosif_from_lanczos) then
   ! Spectral interpolation
   call apply_u_ad(nntrunc,sigmabtrunc,spvartrunc,gptrunc,sptrunc)  
   call interp_lanczos_incr(nntrunc,sigmabtrunc,spvartrunc,sptrunc,nn,sigmab,spvar,.false.,sp)
   call apply_u(nn,sigmab,spvar,sp,b)
   call interp_gp(nntrunc,gptrunc,nn,guess)
   call inverse_b(nn,sigmab,spvar,b,guess,gp)  
else
   ! Grid-point interpolation
   call interp_gp(nntrunc,gptrunc,nn,gp)
end if

end subroutine interp_planczosif_incr

!----------------------------------------------------------------------
! Subroutine: interp_planczosif_lmp
! Purpose: planczosif interpolation for LMP
!----------------------------------------------------------------------
subroutine interp_planczosif_lmp(ni,io,nn,sigmab,spvar,lmp_planczosif,lmp_lanczos)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
real(8),intent(in) :: sigmab(nn(io),io)
real(8),intent(in) :: spvar(nn(io),io)
type(lmp_type),intent(inout) :: lmp_planczosif(io,io)
type(lmp_type),intent(inout),optional :: lmp_lanczos(io,io)

! Local variables
integer :: jo,ip,jp
real(8) :: norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)

do jo=1,io
   ! Copy number of vectors and values
   lmp_planczosif(jo,io)%np = lmp_planczosif(jo,jo)%np
   lmp_planczosif(jo,io)%val(1:lmp_planczosif(jo,io)%np) = lmp_planczosif(jo,jo)%val(1:lmp_planczosif(jo,io)%np)
end do

if (present(lmp_lanczos)) then
   ! Copy number of vectors and values
   lmp_lanczos(io,io)%np = lmp_planczosif(io,io)%np
   lmp_lanczos(io,io)%val = lmp_planczosif(io,io)%val

   if (io>1) then
      ! Allocation
      allocate(avec(nn(io-1)))
      allocate(bvec(nn(io-1)))
      allocate(cvec(nn(io-1)))
      allocate(evec(nn(io-1)))

      ! Transform Ritz vectors from PLanczosIF space to Lanczos space
      do ip=1,lmp_planczosif(io,io)%np
         ! Initial vector
         evec = lmp_planczosif(io,io)%ritzvec_trunc(1:nn(io-1),ip)

         ! Apply PLanczosIF LMP
         call lmp_apply(nn(io-1),io-1,lmp_planczosif(1:io-1,io-1),evec,cvec)

         ! Apply adjoint square-root of B
         call apply_u_ad(nn(io-1),sigmab(1:nn(io-1),io-1),spvar(1:nn(io-1),io-1),cvec,bvec)

         ! Apply Lanczos LMP square-root inverse
         call lmp_apply_sqrt_inv(nn(io-1),io-1,lmp_lanczos(1:io-1,io-1),bvec,avec)

         ! Final vector
         lmp_lanczos(io,io)%ritzvec_trunc(1:nn(io-1),ip) = avec
      end do

      ! Release memory
      deallocate(avec)
      deallocate(bvec)
      deallocate(cvec)
      deallocate(evec)
   end if

   ! Interpolate preconditioning vectors
   call interp_lanczos_lmp(ni,io,nn,sigmab,spvar,lmp_lanczos)

   ! Allocation
   allocate(avec(nn(io)))
   allocate(bvec(nn(io)))
   allocate(cvec(nn(io)))
   allocate(dvec(nn(io)))
   allocate(evec(nn(io)))
   allocate(guess(nn(io)))

   ! Transform Ritz vectors from Lanczos space to PLanczosIF space
   do jo=1,io
      do ip=1,lmp_planczosif(jo,io)%np
         ! Initial vector
         avec = lmp_lanczos(jo,io)%ritzvec(1:nn(io),ip)

         ! Apply Lanczos LMP square-root
         call lmp_apply_sqrt(nn(io),jo-1,lmp_lanczos(1:jo-1,io),avec,bvec)

         ! Apply square-root of B
         call apply_u(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),bvec,dvec)

         ! Apply B inverse
         guess = 0.0
         call inverse_b(nn(io),sigmab(1:nn(io),io),spvar(1:nn(io),io),dvec,guess,cvec)

         ! Apply PLanczosIF LMP inverse
         call lmp_apply_inv(nn(io),jo-1,lmp_planczosif(1:jo-1,io),cvec,evec)

         ! Final vectors
         lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip) = evec
         lmp_planczosif(jo,io)%vec1(1:nn(io),ip) = dvec
         lmp_planczosif(jo,io)%vec2(1:nn(io),ip) = cvec
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
   do jo=1,io
      ! Interpolation of Ritz vectors
      do ip=1,lmp_planczosif(jo,io)%np
         call interp_gp(nn(jo-1),lmp_planczosif(jo,jo)%ritzvec_trunc(1:nn(jo-1),ip),nn(io),lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip))
      end do
   
      ! Gram-Schmidt orthogonalization
      do ip=1,lmp_planczosif(jo,io)%np
         do jp=1,ip-1
            lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip) = lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)-sum(lmp_planczosif(jo,io)%vec1(1:nn(io),jp)*lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip))/sum(lmp_planczosif(jo,io)%vec1(1:nn(io),jp)*lmp_planczosif(jo,io)%ritzvec(1:nn(io),jp))*lmp_planczosif(jo,io)%ritzvec(1:nn(io),jp)
         end do
         if (jo==1) then
            lmp_planczosif(jo,io)%vec2(1:nn(io),ip) = lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)
         else
            call lmp_apply(nn(io),jo-1,lmp_planczosif(1:jo-1,io),lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip),lmp_planczosif(jo,io)%vec2(1:nn(io),ip))
         end if
         call apply_b(nn(io),sigmab(:,io),spvar(:,io),lmp_planczosif(jo,io)%vec2(1:nn(io),ip),lmp_planczosif(jo,io)%vec1(1:nn(io),ip))
         norm = sqrt(sum(lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)*lmp_planczosif(jo,io)%vec1(1:nn(io),ip)))
         lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip) = lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)/norm
         lmp_planczosif(jo,io)%vec1(1:nn(io),ip) = lmp_planczosif(jo,io)%vec1(1:nn(io),ip)/norm
         lmp_planczosif(jo,io)%vec2(1:nn(io),ip) = lmp_planczosif(jo,io)%vec2(1:nn(io),ip)/norm
      end do
   
      if (orthotest) then
         ! Test orthogonality
         do ip=1,lmp_planczosif(jo,io)%np
            do jp=1,lmp_planczosif(jo,io)%np
               write(*,*) 'Ortho test',io,jo,ip,jp,sum(lmp_planczosif(jo,io)%ritzvec(1:nn(io),ip)*lmp_planczosif(jo,io)%vec1(1:nn(io),jp))
            end do
         end do
      end if
   end do
end if

end subroutine interp_planczosif_lmp

end module interp
