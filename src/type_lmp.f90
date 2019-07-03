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
      if (trim(lmp%space)=='model') then
         allocate(lmp%outer(jo)%lancvec1(nn,ni+1))
         allocate(lmp%outer(jo)%lancvec2(nn,ni+1))
         allocate(lmp%outer(jo)%ritzvec1(nn,ni))
         allocate(lmp%outer(jo)%ritzvec2(nn,ni))
      end if
      allocate(lmp%outer(jo)%omega(ni))
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
real(8) :: tmp

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
real(8) :: tmp

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
   end select
end do

end subroutine lmp_apply_ad

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
real(8) :: pxtmp(nn)

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
real(8) :: pxtmp(nn)

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
   end select
end do

end subroutine lmp_apply_sqrt_ad

end module type_lmp
