!----------------------------------------------------------------------
! Module: interp
! Purpose: interpolations
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module interp

use type_bmatrix
use type_hmatrix
use type_lmp
use type_rmatrix
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
! Subroutine: interp_incr_control
! Purpose: increment interpolation in control space
!----------------------------------------------------------------------
subroutine interp_incr_control(nntrunc,bmatrix_trunc,sptrunc,nn,bmatrix,lanczos_from_planczosif,sp)

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
   call interp_incr_linear(nntrunc,bmatrix_trunc,gptrunc,nn,bmatrix,.false.,gp)
   call bmatrix_apply_sqrt_ad(bmatrix,nn,gp,sp)  
else
   ! Spectral interpolation 
   call interp_sp(nntrunc,sptrunc,nn,sp)
end if

end subroutine interp_incr_control

!----------------------------------------------------------------------
! Subroutine: interp_lmp_control
! Purpose: LMP interpolation in control space
!----------------------------------------------------------------------
subroutine interp_lmp_control(ni,io,nn,bmatrix,hmatrix,rmatrix,nobs,lmp_lanczos,lmp_planczosif)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(bmatrix_type),intent(in) :: bmatrix(io)
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: nobs
type(lmp_type),intent(inout) :: lmp_lanczos(io)
type(lmp_type),intent(inout),optional :: lmp_planczosif(io)

! Local variables
integer :: jo,ii,ji,info
real(8) :: proj,norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)
real(8),allocatable :: vtmp(:),vtmp2(:),xtmp(:),ytmp(:),ytmp2(:)
real(8),allocatable :: G(:,:),X(:,:)
logical :: ortho_loss,norm_loss

if (io==1) then
   write(*,*) 'ERROR: interpolation of LMP should not be called for io=1'
   stop
end if
   
select case (trim(lmp_lanczos(1)%mode))
case ('spectral','ritz')
   do jo=2,io
      ! Copy eigenpairs and omega
      lmp_lanczos(io)%outer(jo)%eigenval = lmp_lanczos(jo)%outer(jo)%eigenval
      lmp_lanczos(io)%outer(jo)%eigenvec = lmp_lanczos(jo)%outer(jo)%eigenvec
      lmp_lanczos(io)%outer(jo)%omega = lmp_lanczos(jo)%outer(jo)%omega
   end do

   if (present(lmp_planczosif)) then
      ! Copy eigenpairs and omega
      lmp_planczosif(io)%outer(io)%eigenval = lmp_lanczos(io)%outer(io)%eigenval
      lmp_planczosif(io)%outer(io)%eigenvec = lmp_lanczos(io)%outer(io)%eigenvec
      lmp_planczosif(io)%outer(io)%omega = lmp_lanczos(io)%outer(io)%omega

      ! Allocation
      allocate(avec(nn(io-1)))
      allocate(bvec(nn(io-1)))
      allocate(cvec(nn(io-1)))
      allocate(dvec(nn(io-1)))
      allocate(evec(nn(io-1)))
      allocate(guess(nn(io-1)))

      ! Transform Lanczos vectors from control space to linear space
      do ii=1,ni+1
         ! Initial vector
         avec = lmp_lanczos(io)%outer(io)%lancvec_trunc(1:nn(io-1),ii)

         ! Apply Lanczos LMP square-root
         call lmp_apply_sqrt(lmp_lanczos(io-1),nn(io-1),ni,io-1,avec,bvec)

         ! Apply square-root of B
         call bmatrix_apply_sqrt(bmatrix(io-1),nn(io-1),bvec,dvec)

         ! Apply B inverse
         guess = 0.0
         call bmatrix_apply_inv(bmatrix(io-1),nn(io-1),dvec,guess,cvec)

         ! Apply PLanczosIF LMP inverse
         call lmp_apply_inv(lmp_planczosif(io-1),nn(io-1),ni,io-1,cvec,evec)

         ! Final vector
         lmp_planczosif(io)%outer(io)%lancvec_trunc(1:nn(io-1),ii) = evec
      end do

      ! Release memory
      deallocate(avec)
      deallocate(bvec)
      deallocate(cvec)
      deallocate(dvec)
      deallocate(evec)
      deallocate(guess)
   
      ! Interpolate preconditioning vectors
      call interp_lmp_linear(ni,io,nn,bmatrix,hmatrix,rmatrix,nobs,lmp_planczosif)
   
      ! Allocation
      allocate(avec(nn(io)))
      allocate(bvec(nn(io)))
      allocate(cvec(nn(io)))
      allocate(evec(nn(io)))
   
      ! Transform Lanczos vectors from linear space to control space
      do jo=2,io
         do ii=1,ni+1
            ! Initial vector
            evec = lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)
   
            ! Apply PLanczosIF LMP
            call lmp_apply(lmp_planczosif(io),nn(io),ni,jo-1,evec,cvec)
   
            ! Apply adjoint square-root of B
            call bmatrix_apply_sqrt_ad(bmatrix(io),nn(io),cvec,bvec)
   
            ! Apply Lanczos LMP square-root inverse
            call lmp_apply_sqrt_inv(lmp_lanczos(io),nn(io),ni,jo-1,bvec,avec)
   
            ! Final vectors
            lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii) = avec
         end do

         ! Ritz vectors
         lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_lanczos(io)%outer(jo)%eigenvec)
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
   
      do jo=2,io
         ! Interpolation of Lanczos vectors
         do ii=1,ni+1
            call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%lancvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii))
         end do
         
         ! Gram-Schmidt orthogonalization
         do ii=1,ni+1
            do ji=1,ii-1
               proj = sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji)*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii))/sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji)*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji))
               lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii) = lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)-proj*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji)
               if (abs(proj)>1.0e-6) ortho_loss = .true.
            end do
            norm = sqrt(sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)**2))
            lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii) = lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)/norm
            if (abs(norm-1.0)>1.0e-6) norm_loss = .true.
         end do

         ! Ritz vectors
         lmp_lanczos(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_lanczos(io)%outer(jo)%eigenvec)
      end do
   
      if (orthotest) then
         ! Print orthogonality/norm loss
         if (ortho_loss) write(*,*) '     Orthogonalization loss in Lanczos'
         if (norm_loss) write(*,*) '     Normalization loss in Lanczos'
      end if
   end if
   
   if (orthotest) then
      ! Test orthogonality
      do jo=2,io
         do ii=1,ni+1
            do ji=1,ni+1
               write(*,*) '     Lanczos ortho test',io,jo,ii,ji,sum(lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)*lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ji))
            end do
         end do
      end do
   end if
case ('quasi-newton')
   do jo=2,io
      ! Interpolation of descent vectors
      do ii=1,ni
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%pvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%pvec(1:nn(io),ii))
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%qvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%qvec(1:nn(io),ii))
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%rvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%rvec(1:nn(io),ii))
      end do
   end do
case ('full')
   ! Allocation
   allocate(vtmp(nn(io)))
   allocate(vtmp2(nn(io)))
   allocate(xtmp(nn(io)))
   allocate(ytmp(nobs))
   allocate(ytmp2(nobs))
   allocate(G(ni,ni))
   allocate(X(ni,ni))

   do jo=2,io
      ! Interpolation of descent vector
      do ii=1,ni
         call interp_sp(nn(jo-1),lmp_lanczos(jo)%outer(jo)%pvec_trunc(1:nn(jo-1),ii),nn(io),lmp_lanczos(io)%outer(jo)%pvec(1:nn(io),ii))
      end do
   end do

   do jo=2,io
      ! Compute the other vector
      do ii=1,ni
         call lmp_apply_sqrt(lmp_lanczos(io),nn(io),ni,jo-1,lmp_lanczos(io)%outer(jo)%pvec(:,ii),vtmp)
         call bmatrix_apply_sqrt(bmatrix(io),nn(io),vtmp,xtmp)
         call hmatrix_apply(hmatrix,nn(io),xtmp,nobs,ytmp)
         call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
         call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn(io),xtmp)
         call bmatrix_apply_sqrt_ad(bmatrix(io),nn(io),xtmp,vtmp2)
         vtmp = vtmp+vtmp2
         call lmp_apply_sqrt_ad(lmp_lanczos(io),nn(io),ni,jo-1,vtmp,lmp_lanczos(io)%outer(jo)%qvec(:,ii))
      end do

      ! Compute Cholesky decompositions
      G = matmul(transpose(lmp_lanczos(io)%outer(jo)%qvec),lmp_lanczos(io)%outer(jo)%pvec)
      call dpotrf('U',ni,G,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dpotrf: ',info
         stop
      end if
      do ii=2,ni
         G(ii,1:ii-1) = 0.0
      end do
      X = matmul(transpose(lmp_lanczos(io)%outer(jo)%pvec),lmp_lanczos(io)%outer(jo)%pvec)
      call dpotrf('U',ni,X,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dpotrf: ',info
         stop
      end if
      do ii=2,ni
         X(ii,1:ii-1) = 0.0
      end do

      ! Inverse upper triangular matrices
      lmp_lanczos(io)%outer(jo)%Ginv = G
      call dtrtri('U','N',ni,lmp_lanczos(io)%outer(jo)%Ginv,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dtrtri: ',info
         stop
      end if
      lmp_lanczos(io)%outer(jo)%Xinv = X
      call dtrtri('U','N',ni,lmp_lanczos(io)%outer(jo)%Xinv,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dtrtri: ',info
         stop
      end if
   end do
end select

end subroutine interp_lmp_control

!----------------------------------------------------------------------
! Subroutine: interp_incr_linear
! Purpose: increment interpolation in linear space
!----------------------------------------------------------------------
subroutine interp_incr_linear(nntrunc,bmatrix_trunc,gptrunc,nn,bmatrix,planczosif_from_lanczos,gp)

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
   call interp_incr_control(nntrunc,bmatrix_trunc,sptrunc,nn,bmatrix,.false.,sp)
   call bmatrix_apply_sqrt(bmatrix,nn,sp,b)
   call interp_gp(nntrunc,gptrunc,nn,guess)
   call bmatrix_apply_inv(bmatrix,nn,b,guess,gp)  
else
   ! Grid-point interpolation
   call interp_gp(nntrunc,gptrunc,nn,gp)
end if

end subroutine interp_incr_linear

!----------------------------------------------------------------------
! Subroutine: interp_lmp_linear
! Purpose: LMP interpolation in linear space
!----------------------------------------------------------------------
subroutine interp_lmp_linear(ni,io,nn,bmatrix,hmatrix,rmatrix,nobs,lmp_planczosif,lmp_lanczos)

implicit none

! Passed variables
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nn(io)
type(bmatrix_type),intent(in) :: bmatrix(io)
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: nobs
type(lmp_type),intent(inout) :: lmp_planczosif(io)
type(lmp_type),intent(inout),optional :: lmp_lanczos(io)

! Local variables
integer :: jo,ii,ji,info
real(8) :: proj,norm
real(8),allocatable :: avec(:),bvec(:),cvec(:),dvec(:),evec(:)
real(8),allocatable :: guess(:)
real(8),allocatable :: xtmp(:),ytmp(:),ytmp2(:)
real(8),allocatable :: G(:,:)
logical :: ortho_loss,norm_loss

if (io==1) then
   write(*,*) 'ERROR: interpolation of LMP should not be called for io=1'
   stop
end if

select case (trim(lmp_planczosif(1)%mode))
case ('spectral','ritz')
   do jo=2,io
      ! Copy eigenpairs and omega
      lmp_planczosif(io)%outer(jo)%eigenval = lmp_planczosif(jo)%outer(jo)%eigenval
      lmp_planczosif(io)%outer(jo)%eigenvec = lmp_planczosif(jo)%outer(jo)%eigenvec
      lmp_planczosif(io)%outer(jo)%omega = lmp_planczosif(jo)%outer(jo)%omega
   end do

   if (present(lmp_lanczos)) then
      ! Copy eigenpairs and last beta
      lmp_lanczos(io)%outer(io)%eigenval = lmp_planczosif(io)%outer(io)%eigenval
      lmp_lanczos(io)%outer(io)%eigenvec = lmp_planczosif(io)%outer(io)%eigenvec
      lmp_lanczos(io)%outer(io)%omega = lmp_planczosif(io)%outer(io)%omega
  
      if (io>1) then
         ! Allocation
         allocate(avec(nn(io-1)))
         allocate(bvec(nn(io-1)))
         allocate(cvec(nn(io-1)))
         allocate(evec(nn(io-1)))
   
         ! Transform Lanczos vectors from linear space to control space
         do ii=1,ni+1
            ! Initial vector
            evec = lmp_planczosif(io)%outer(io)%lancvec_trunc(1:nn(io-1),ii)
   
            ! Apply PLanczosIF LMP
            call lmp_apply(lmp_planczosif(io-1),nn(io-1),ni,io-1,evec,cvec)
   
            ! Apply adjoint square-root of B
            call bmatrix_apply_sqrt_ad(bmatrix(io-1),nn(io-1),cvec,bvec)
   
            ! Apply Lanczos LMP square-root inverse
            call lmp_apply_sqrt_inv(lmp_lanczos(io-1),nn(io-1),ni,io-1,bvec,avec)
   
            ! Final vector
            lmp_lanczos(io)%outer(io)%lancvec_trunc(1:nn(io-1),ii) = avec
         end do
   
         ! Release memory
         deallocate(avec)
         deallocate(bvec)
         deallocate(cvec)
         deallocate(evec)
      end if
   
      ! Interpolate preconditioning vectors
      call interp_lmp_control(ni,io,nn,bmatrix,hmatrix,rmatrix,nobs,lmp_lanczos)
   
      ! Allocation
      allocate(avec(nn(io)))
      allocate(bvec(nn(io)))
      allocate(cvec(nn(io)))
      allocate(dvec(nn(io)))
      allocate(evec(nn(io)))
      allocate(guess(nn(io)))
   
      ! Transform Lanczos vectors from control space to linear space
      do jo=2,io
         do ii=1,ni+1
            ! Initial vector
            avec = lmp_lanczos(io)%outer(jo)%lancvec(1:nn(io),ii)
   
            ! Apply Lanczos LMP square-root
            call lmp_apply_sqrt(lmp_lanczos(io),nn(io),ni,jo-1,avec,bvec)
   
            ! Apply square-root of B
            call bmatrix_apply_sqrt(bmatrix(io),nn(io),bvec,dvec)
   
            ! Apply B inverse
            guess = 0.0
            call bmatrix_apply_inv(bmatrix(io),nn(io),dvec,guess,cvec)
   
            ! Apply PLanczosIF LMP inverse
            call lmp_apply_inv(lmp_planczosif(io),nn(io),ni,jo-1,cvec,evec)
   
            ! Final vectors
            lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii) = evec
            lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii) = dvec
            lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii) = cvec
         end do

         ! Ritz vectors
         lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
         lmp_planczosif(io)%outer(jo)%ritzvec1(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
         lmp_planczosif(io)%outer(jo)%ritzvec2(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
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
   
      do jo=2,io
         ! Interpolation of Lanczos vectors
         do ii=1,ni+1
            call interp_gp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%lancvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii))
         end do
      
         ! Gram-Schmidt orthogonalization
         do ii=1,ni+1
            do ji=1,ii-1
               proj = sum(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ji)*lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii))/sum(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ji)*lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ji))
               lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii) = lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)-proj*lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ji)
               if (abs(proj)>1.0e-6) ortho_loss = .true.
            end do
            call lmp_apply(lmp_planczosif(io),nn(io),ni,jo-1,lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii),lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii))
            call bmatrix_apply(bmatrix(io),nn(io),lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii),lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii))
            norm = sqrt(sum(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)*lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii)))
            lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii) = lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)/norm
            lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii) = lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ii)/norm
            lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii) = lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),ii)/norm
            if (abs(norm-1.0)>1.0e-6) norm_loss = .true.
         end do

         ! Ritz vectors
         lmp_planczosif(io)%outer(jo)%ritzvec(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
         lmp_planczosif(io)%outer(jo)%ritzvec1(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
         lmp_planczosif(io)%outer(jo)%ritzvec2(1:nn(io),1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec2(1:nn(io),1:ni),lmp_planczosif(io)%outer(jo)%eigenvec)
      end do
   
      if (orthotest) then
         ! Print orthogonality/norm loss
         if (ortho_loss) write(*,*) '     Orthogonalization loss in PLanczosIF'
         if (norm_loss) write(*,*) '     Normalization loss in PLanczosIF'
      end if
   end if
   
   if (orthotest) then
      ! Test orthogonality
      do jo=2,io
         do ii=1,ni+1
            do ji=1,ni
               write(*,*) '     PLanczosIF ortho test',io,jo,ii,ji,sum(lmp_planczosif(io)%outer(jo)%lancvec(1:nn(io),ii)*lmp_planczosif(io)%outer(jo)%lancvec1(1:nn(io),ji))
            end do
         end do
      end do
   end if
case ('quasi-newton')
   do jo=2,io
      ! Interpolation of descent vectors
      do ii=1,ni
         call interp_sp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%pbarvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%pbarvec(1:nn(io),ii))
         call interp_sp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%pvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%pvec(1:nn(io),ii))
         call interp_sp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%qvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%qvec(1:nn(io),ii))
         call interp_sp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%tvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%tvec(1:nn(io),ii))
      end do
   end do
case ('full')
   ! Allocation
   allocate(xtmp(nn(io)))
   allocate(ytmp(nobs))
   allocate(ytmp2(nobs))
   allocate(G(ni,ni))

   do jo=2,io
      ! Interpolation of descent vector
      do ii=1,ni
         call interp_sp(nn(jo-1),lmp_planczosif(jo)%outer(jo)%pbarvec_trunc(1:nn(jo-1),ii),nn(io),lmp_planczosif(io)%outer(jo)%pbarvec(1:nn(io),ii))
      end do

      ! Compute the other vectors
      do ii=1,ni
         call bmatrix_apply(bmatrix(io),nn(io),lmp_planczosif(io)%outer(jo)%pbarvec(:,ii),lmp_planczosif(io)%outer(jo)%pvec(:,ii))
         call hmatrix_apply(hmatrix,nn(io),lmp_planczosif(io)%outer(jo)%pvec(:,ii),nobs,ytmp)
         call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
         call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn(io),xtmp)
         lmp_planczosif(io)%outer(jo)%qvec(:,ii) = lmp_planczosif(io)%outer(jo)%pbarvec(:,ii)+xtmp

         call hmatrix_apply(hmatrix,nn(io),lmp_planczosif(io)%outer(jo)%pvec(:,ii),nobs,ytmp)
         call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
         call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn(io),xtmp)
         call bmatrix_apply(bmatrix(io),nn(io),xtmp,lmp_planczosif(io)%outer(jo)%tvec(:,ii))
         lmp_planczosif(io)%outer(jo)%tvec(:,ii) = lmp_planczosif(io)%outer(jo)%pvec(:,ii)+lmp_planczosif(io)%outer(jo)%tvec(:,ii)
      end do

      ! Compute Cholesky decomposition
      G = matmul(transpose(lmp_planczosif(io)%outer(jo)%qvec),lmp_planczosif(io)%outer(jo)%pvec)
      call dpotrf('U',ni,G,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dpotrf: ',info
         stop
      end if
      do ii=2,ni
         G(ii,1:ii-1) = 0.0
      end do

      ! Inverse upper triangular matrice
      call dtrtri('U','N',ni,G,ni,info)
      if (info/=0) then
         write(*,*) 'Error in dtrtri: ',info
         stop
      end if

      ! Compute final inverse
      lmp_planczosif(io)%outer(jo)%Ginv = matmul(G,transpose(G))
   end do
end select

end subroutine interp_lmp_linear

end module interp
