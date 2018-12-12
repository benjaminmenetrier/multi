!----------------------------------------------------------------------
! Module: type_lmp
! Purpose: LMP preconditioner
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_lmp

implicit none

type outer_type
   ! Spectral / Ritz LMP
   real(8),allocatable :: eigenval(:)
   real(8),allocatable :: eigenvec(:,:)
   real(8),allocatable :: omega(:)
   real(8),allocatable :: lancvec_trunc(:,:)
   real(8),allocatable :: lancvec(:,:)
   real(8),allocatable :: lancvec1(:,:)
   real(8),allocatable :: lancvec2(:,:)
   real(8),allocatable :: ritzvec(:,:)
   real(8),allocatable :: ritzvec1(:,:)
   real(8),allocatable :: ritzvec2(:,:)

   ! Quasi-Newton LMP
   real(8),allocatable :: pvec_trunc(:,:)
   real(8),allocatable :: pbarvec_trunc(:,:)
   real(8),allocatable :: qvec_trunc(:,:)
   real(8),allocatable :: rvec_trunc(:,:)
   real(8),allocatable :: tvec_trunc(:,:)
   real(8),allocatable :: pvec(:,:)
   real(8),allocatable :: pbarvec(:,:)
   real(8),allocatable :: qvec(:,:)
   real(8),allocatable :: rvec(:,:)
   real(8),allocatable :: tvec(:,:)

   ! Full LMP
   real(8),allocatable :: Ginv(:,:)
   real(8),allocatable :: Xinv(:,:)
end type outer_type

type lmp_type
   character(len=1024) :: mode
   character(len=1024) :: space
   integer :: io
   type(outer_type),allocatable :: outer(:)
end type lmp_type

contains

!----------------------------------------------------------------------
! Subroutine: lmp_alloc
! Purpose: allocate LMP
!----------------------------------------------------------------------
subroutine lmp_alloc(lmp,nn,ni,io,mode,space)

implicit none

! Passed variables
type(lmp_type),intent(inout) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
character(len=*),intent(in) :: mode
character(len=*),intent(in) :: space

! Local variables
integer :: jo

! Release memory
if (allocated(lmp%outer)) then
   do jo=1,lmp%io
      if (allocated(lmp%outer(jo)%eigenval)) deallocate(lmp%outer(jo)%eigenval)
      if (allocated(lmp%outer(jo)%eigenvec)) deallocate(lmp%outer(jo)%eigenvec)
      if (allocated(lmp%outer(jo)%omega)) deallocate(lmp%outer(jo)%omega)
      if (allocated(lmp%outer(jo)%lancvec_trunc)) deallocate(lmp%outer(jo)%lancvec_trunc)
      if (allocated(lmp%outer(jo)%lancvec)) deallocate(lmp%outer(jo)%lancvec)
      if (allocated(lmp%outer(jo)%lancvec1)) deallocate(lmp%outer(jo)%lancvec1)
      if (allocated(lmp%outer(jo)%lancvec2)) deallocate(lmp%outer(jo)%lancvec2)
      if (allocated(lmp%outer(jo)%ritzvec)) deallocate(lmp%outer(jo)%ritzvec)
      if (allocated(lmp%outer(jo)%ritzvec1)) deallocate(lmp%outer(jo)%ritzvec1)
      if (allocated(lmp%outer(jo)%ritzvec2)) deallocate(lmp%outer(jo)%ritzvec2)
      if (allocated(lmp%outer(jo)%pvec_trunc)) deallocate(lmp%outer(jo)%pvec_trunc)
      if (allocated(lmp%outer(jo)%pbarvec_trunc)) deallocate(lmp%outer(jo)%pbarvec_trunc)
      if (allocated(lmp%outer(jo)%qvec_trunc)) deallocate(lmp%outer(jo)%qvec_trunc)
      if (allocated(lmp%outer(jo)%rvec_trunc)) deallocate(lmp%outer(jo)%rvec_trunc)
      if (allocated(lmp%outer(jo)%tvec_trunc)) deallocate(lmp%outer(jo)%tvec_trunc)
      if (allocated(lmp%outer(jo)%pvec)) deallocate(lmp%outer(jo)%pvec)
      if (allocated(lmp%outer(jo)%pbarvec)) deallocate(lmp%outer(jo)%pbarvec)
      if (allocated(lmp%outer(jo)%qvec)) deallocate(lmp%outer(jo)%qvec)
      if (allocated(lmp%outer(jo)%rvec)) deallocate(lmp%outer(jo)%rvec)
      if (allocated(lmp%outer(jo)%tvec)) deallocate(lmp%outer(jo)%tvec)
      if (allocated(lmp%outer(jo)%Ginv)) deallocate(lmp%outer(jo)%Ginv)
      if (allocated(lmp%outer(jo)%Xinv)) deallocate(lmp%outer(jo)%Xinv)
   end do
   deallocate(lmp%outer)
end if

! Allocation
lmp%mode = trim(mode)
lmp%space = trim(space)
lmp%io = io
allocate(lmp%outer(io))
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      allocate(lmp%outer(jo)%eigenval(ni))
      allocate(lmp%outer(jo)%eigenvec(ni,ni))
      allocate(lmp%outer(jo)%lancvec_trunc(nn,ni+1))
      allocate(lmp%outer(jo)%lancvec(nn,ni+1))
      allocate(lmp%outer(jo)%ritzvec(nn,ni))
      if (trim(lmp%space)=='linear') then
         allocate(lmp%outer(jo)%lancvec1(nn,ni+1))
         allocate(lmp%outer(jo)%lancvec2(nn,ni+1))
         allocate(lmp%outer(jo)%ritzvec1(nn,ni))
         allocate(lmp%outer(jo)%ritzvec2(nn,ni))
      end if
      allocate(lmp%outer(jo)%omega(ni))
   case ('quasi-newton')
      allocate(lmp%outer(jo)%pvec_trunc(nn,ni))
      allocate(lmp%outer(jo)%pvec(nn,ni))
      allocate(lmp%outer(jo)%qvec_trunc(nn,ni))
      allocate(lmp%outer(jo)%qvec(nn,ni))
      if (trim(lmp%space)=='control') then
         allocate(lmp%outer(jo)%rvec_trunc(nn,ni))
         allocate(lmp%outer(jo)%rvec(nn,ni))
      end if
      if (trim(lmp%space)=='linear') then
         allocate(lmp%outer(jo)%pbarvec_trunc(nn,ni))
         allocate(lmp%outer(jo)%pbarvec(nn,ni))
         allocate(lmp%outer(jo)%tvec_trunc(nn,ni))
         allocate(lmp%outer(jo)%tvec(nn,ni))
      end if
   case ('full')
      allocate(lmp%outer(jo)%pvec(nn,ni))
      allocate(lmp%outer(jo)%qvec(nn,ni))
      allocate(lmp%outer(jo)%Ginv(ni,ni))
      if (trim(lmp%space)=='control') then
         allocate(lmp%outer(jo)%pvec_trunc(nn,ni))
         allocate(lmp%outer(jo)%Xinv(ni,ni))
      end if
      if (trim(lmp%space)=='linear') then
         allocate(lmp%outer(jo)%pbarvec_trunc(nn,ni))
         allocate(lmp%outer(jo)%pbarvec(nn,ni))
         allocate(lmp%outer(jo)%tvec(nn,ni))
      end if
   end select
end do

end subroutine lmp_alloc

!----------------------------------------------------------------------
! Subroutine: lmp_apply
! Purpose: apply LMP
!----------------------------------------------------------------------
subroutine lmp_apply(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: tmp,tau,eta(ni),zeta,pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Spectral LMP
      do ii=1,ni
         px = px+lmp%outer(jo)%ritzvec2(:,ii)*(1.0/lmp%outer(jo)%eigenval(ii)-1.0)*sum(lmp%outer(jo)%ritzvec1(:,ii)*x)
      end do

      ! Last Ritz term
      if (trim(lmp%mode)=='ritz') then
         tmp = 0.0
         do ii=1,ni
            tmp = tmp+lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec1(:,ii)*x)
         end do
         do ii=1,ni
            px = px-lmp%outer(jo)%ritzvec2(:,ii)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%lancvec1(:,ni+1)*x)-lmp%outer(jo)%lancvec2(:,ni+1)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec1(:,ii)*x)+lmp%outer(jo)%ritzvec2(:,ii)*lmp%outer(jo)%omega(ii)*tmp
         end do
      end if
   case ('quasi-newton')
      do ii=ni,1,-1
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         eta(ii) = tau*sum(lmp%outer(jo)%pvec(:,ii)*px)
         px = px-eta(ii)*lmp%outer(jo)%qvec(:,ii)
      end do
      do ii=1,ni
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         zeta = tau*sum(lmp%outer(jo)%tvec(:,ii)*px)
         px = px+(eta(ii)-zeta)*lmp%outer(jo)%pbarvec(:,ii)
      end do
   case ('full')
      pxtmp = px
      px = px-matmul(lmp%outer(jo)%qvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%pvec),px)))
      px = px-matmul(lmp%outer(jo)%pbarvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%tvec),px)))
      px = px+matmul(lmp%outer(jo)%pbarvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%pvec),pxtmp)))
   end select
end do

end subroutine lmp_apply

!----------------------------------------------------------------------
! Subroutine: lmp_apply_ad
! Purpose: apply LMP adjoint
!----------------------------------------------------------------------
subroutine lmp_apply_ad(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: tmp,tau,eta(ni),zeta,pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=io,2,-1
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Spectral LMP adjoint
      do ii=1,ni
         px = px+lmp%outer(jo)%ritzvec1(:,ii)*(1.0/lmp%outer(jo)%eigenval(ii)-1.0)*sum(lmp%outer(jo)%ritzvec2(:,ii)*x)
      end do

      ! Last Ritz term adjoint
      if (trim(lmp%mode)=='ritz') then
         tmp = 0.0
         do ii=1,ni
            tmp = tmp+lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec2(:,ii)*x)
         end do
         do ii=1,ni
            px = px-lmp%outer(jo)%lancvec1(:,ni+1)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec2(:,ii)*x)-lmp%outer(jo)%ritzvec1(:,ii)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%lancvec2(:,ni+1)*x)+lmp%outer(jo)%ritzvec1(:,ii)*lmp%outer(jo)%omega(ii)*tmp
         end do
      end if
   case ('quasi-newton')
      do ii=ni,1,-1
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         eta(ii) = tau*sum(lmp%outer(jo)%pbarvec(:,ii)*px)
         px = px-eta(ii)*lmp%outer(jo)%tvec(:,ii)
      end do
      do ii=1,ni
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         zeta = tau*sum(lmp%outer(jo)%qvec(:,ii)*px)
         px = px+(eta(ii)-zeta)*lmp%outer(jo)%pvec(:,ii)
      end do
   case ('full')
      pxtmp = px
      px = px-matmul(lmp%outer(jo)%tvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%pbarvec),px)))
      px = px-matmul(lmp%outer(jo)%pvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%qvec),px)))
      px = px+matmul(lmp%outer(jo)%pvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%pbarvec),pxtmp)))
   end select
end do

end subroutine lmp_apply_ad

!----------------------------------------------------------------------
! Subroutine: lmp_apply_inv
! Purpose: apply LMP, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_inv(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Initialization
      pxtmp = px

      ! Spectral LMP inverse
      do ii=1,ni
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(lmp%outer(jo)%eigenval(ii)-1.0)*sum(lmp%outer(jo)%ritzvec1(:,ii)*px)
      end do

      ! Last Ritz term inverse
      if (trim(lmp%mode)=='ritz') then
         do ii=1,ni
! TODO
         end do
      end if

      ! Update
      px = pxtmp
   case ('quasi-newton')
! TODO
   case ('full')
! TODO
   end select
end do

end subroutine lmp_apply_inv

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt
! Purpose: apply LMP square-root
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: pxtmp(nn),tau

! Initialization
px = x

! Previous outer loops
do jo=io,2,-1
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Initialization
      pxtmp = px

      ! Spectral LMP square-root
      do ii=1,ni
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(1.0/sqrt(lmp%outer(jo)%eigenval(ii))-1.0)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
      end do

      ! Last Ritz term square-root
      if (trim(lmp%mode)=='ritz') then
         do ii=1,ni
            pxtmp = pxtmp-lmp%outer(jo)%ritzvec(:,ii)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%lancvec(:,ni+1)*px)
         end do
      end if

      ! Update
      px = pxtmp
   case ('quasi-newton')
      do ii=1,ni
         ! Quasi-Newton LMP square-root
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         px = px-lmp%outer(jo)%pvec(:,ii)*sum((tau*lmp%outer(jo)%qvec(:,ii)+sqrt(tau/sum(lmp%outer(jo)%rvec(:,ii)**2))*lmp%outer(jo)%rvec(:,ii))*px)
      end do
   case ('full')
      px = px-matmul(lmp%outer(jo)%pvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%Ginv),matmul(transpose(lmp%outer(jo)%qvec),px))))+matmul(lmp%outer(jo)%pvec,matmul(lmp%outer(jo)%Ginv,matmul(lmp%outer(jo)%Xinv,matmul(transpose(lmp%outer(jo)%pvec),px))))
   end select
end do

end subroutine lmp_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_ad
! Purpose: apply LMP square-root adjoint
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_ad(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: pxtmp(nn),tau

! Initialization
px = x

! Previous outer loops
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Initialization
      pxtmp = px

      ! Spectral LMP square-root adjoint
      do ii=1,ni
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(1.0/sqrt(lmp%outer(jo)%eigenval(ii))-1.0)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
      end do

      ! Last Ritz term square-root adjoint
      if (trim(lmp%mode)=='ritz') then
         do ii=1,ni
            pxtmp = pxtmp-lmp%outer(jo)%lancvec(:,ni+1)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
         end do
      end if

      ! Update
      px = pxtmp
   case ('quasi-newton')
      do ii=ni,1,-1
         ! Quasi-Newton LMP square-root
         tau = 1.0/sum(lmp%outer(jo)%pvec(:,ii)*lmp%outer(jo)%qvec(:,ii))
         px = px-(tau*lmp%outer(jo)%qvec(:,ii)+sqrt(tau/sum(lmp%outer(jo)%rvec(:,ii)**2))*lmp%outer(jo)%rvec(:,ii))*sum(lmp%outer(jo)%pvec(:,ii)*px)
      end do
   case ('full')
      px = px-matmul(lmp%outer(jo)%qvec,matmul(lmp%outer(jo)%Ginv,matmul(transpose(lmp%outer(jo)%Ginv),matmul(transpose(lmp%outer(jo)%pvec),px))))+matmul(lmp%outer(jo)%pvec,matmul(transpose(lmp%outer(jo)%Xinv),matmul(transpose(lmp%outer(jo)%Ginv),matmul(transpose(lmp%outer(jo)%pvec),px))))
   end select
end do

end subroutine lmp_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_inv
! Purpose: apply LMP square-root, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_inv(lmp,nn,ni,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ii
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Initialization
      pxtmp = px

      ! Spectral LMP square-root inverse
      do ii=1,ni
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(sqrt(lmp%outer(jo)%eigenval(ii))-1.0)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
      end do

      ! Last Ritz term square-root inverse
      if (trim(lmp%mode)=='ritz') then
         do ii=1,ni
            pxtmp = pxtmp-lmp%outer(jo)%ritzvec(:,ii)*lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%lancvec(:,ni+1)*px)
         end do
      end if

      ! Update
      px = pxtmp
   case ('quasi-newton')
! TODO
   case ('full')
! TODO
   end select
end do

end subroutine lmp_apply_sqrt_inv

end module type_lmp
