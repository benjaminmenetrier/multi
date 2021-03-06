!----------------------------------------------------------------------
! Module: type_lmp
! Purpose: LMP preconditioner
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_lmp

use tools_kinds
use type_geom

implicit none

type outer_type
   ! Spectral / Ritz LMP
   real(kind_real),allocatable :: eigenval(:)
   real(kind_real),allocatable :: eigenvec(:,:)
   real(kind_real),allocatable :: omega(:)
   real(kind_real),allocatable :: lancvec(:,:)
   real(kind_real),allocatable :: lancvec1(:,:)
   real(kind_real),allocatable :: lancvec2(:,:)
   real(kind_real),allocatable :: ritzvec(:,:)
   real(kind_real),allocatable :: ritzvec1(:,:)
   real(kind_real),allocatable :: ritzvec2(:,:)
end type outer_type

type lmp_type
   character(len=1024) :: mode
   integer :: io
   real(kind_real),allocatable :: lancvec_trunc(:,:)
   type(outer_type),allocatable :: outer(:)
contains
   procedure :: alloc => lmp_alloc
   procedure :: apply => lmp_apply
   procedure :: apply_ad => lmp_apply_ad
   procedure :: apply_sqrt => lmp_apply_sqrt
   procedure :: apply_sqrt_ad => lmp_apply_sqrt_ad
end type lmp_type

contains

!----------------------------------------------------------------------
! Subroutine: lmp_alloc
! Purpose: allocate LMP
!----------------------------------------------------------------------
subroutine lmp_alloc(lmp,no,geom,ni,io,mode,algo_name)

implicit none

! Passed variables
class(lmp_type),intent(inout) :: lmp
integer,intent(in) :: no
type(geom_type),intent(in) :: geom(no)
integer,intent(in) :: ni
integer,intent(in) :: io
character(len=*),intent(in) :: mode
character(len=*),intent(in) :: algo_name

! Local variables
integer :: jo

! Release memory
if (allocated(lmp%outer)) then
   if (allocated(lmp%lancvec_trunc)) deallocate(lmp%lancvec_trunc)
   do jo=2,lmp%io
      if (allocated(lmp%outer(jo)%eigenval)) deallocate(lmp%outer(jo)%eigenval)
      if (allocated(lmp%outer(jo)%eigenvec)) deallocate(lmp%outer(jo)%eigenvec)
      if (allocated(lmp%outer(jo)%omega)) deallocate(lmp%outer(jo)%omega)
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
lmp%io = io
if (io>1) then
   allocate(lmp%lancvec_trunc(geom(io-1)%nh,ni+1))
   allocate(lmp%outer(2:io))
   do jo=2,io
      select case (trim(lmp%mode))
      case ('spectral','ritz')
         allocate(lmp%outer(jo)%eigenval(ni))
         allocate(lmp%outer(jo)%eigenvec(ni,ni))
         allocate(lmp%outer(jo)%lancvec(geom(io)%nh,ni+1))
         allocate(lmp%outer(jo)%ritzvec(geom(io)%nh,ni))
         if (trim(algo_name)=='planczosif') then
            allocate(lmp%outer(jo)%lancvec1(geom(io)%nh,ni+1))
            allocate(lmp%outer(jo)%lancvec2(geom(io)%nh,ni+1))
            allocate(lmp%outer(jo)%ritzvec1(geom(io)%nh,ni))
            allocate(lmp%outer(jo)%ritzvec2(geom(io)%nh,ni))
         end if
         allocate(lmp%outer(jo)%omega(ni))
      end select
   end do
end if

! Initialization at missing value
if (io>1) then
   lmp%lancvec_trunc = msv_real
   do jo=2,io
      select case (trim(lmp%mode))
      case ('spectral','ritz')
         lmp%outer(jo)%eigenval = msv_real
         lmp%outer(jo)%eigenvec = msv_real
         lmp%outer(jo)%lancvec = msv_real
         lmp%outer(jo)%ritzvec = msv_real
         if (trim(algo_name)=='planczosif') then
            lmp%outer(jo)%lancvec1 = msv_real
            lmp%outer(jo)%lancvec2 = msv_real
            lmp%outer(jo)%ritzvec1 = msv_real
            lmp%outer(jo)%ritzvec2 = msv_real
         end if
         lmp%outer(jo)%omega = msv_real
      end select
   end do
end if

end subroutine lmp_alloc

!----------------------------------------------------------------------
! Subroutine: lmp_apply
! Purpose: apply LMP
!----------------------------------------------------------------------
subroutine lmp_apply(lmp,geom,ni,io,x,px)

implicit none

! Passed variables
class(lmp_type),intent(in) :: lmp
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni
integer,intent(in) :: io
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: px(geom%nh)

! Local variables
integer :: jo,ii
real(kind_real) :: tmp

! Initialization
px = x

! Previous outer loops
do jo=2,io
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Spectral LMP
      do ii=1,ni
         px = px+lmp%outer(jo)%ritzvec2(:,ii)*(one/lmp%outer(jo)%eigenval(ii)-one)*sum(lmp%outer(jo)%ritzvec1(:,ii)*x)
      end do

      ! Last Ritz term
      if (trim(lmp%mode)=='ritz') then
         tmp = zero
         do ii=1,ni
            tmp = tmp+lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec1(:,ii)*x)
         end do
         do ii=1,ni
            px = px-lmp%outer(jo)%ritzvec2(:,ii)*lmp%outer(jo)%omega(ii) &
 & *sum(lmp%outer(jo)%lancvec1(:,ni+1)*x)-lmp%outer(jo)%lancvec2(:,ni+1)*lmp%outer(jo)%omega(ii) &
 & *sum(lmp%outer(jo)%ritzvec1(:,ii)*x)+lmp%outer(jo)%ritzvec2(:,ii)*lmp%outer(jo)%omega(ii)*tmp
         end do
      end if
   end select
end do

end subroutine lmp_apply

!----------------------------------------------------------------------
! Subroutine: lmp_apply_ad
! Purpose: apply LMP adjoint
!----------------------------------------------------------------------
subroutine lmp_apply_ad(lmp,geom,ni,io,x,px)

implicit none

! Passed variables
class(lmp_type),intent(in) :: lmp
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni
integer,intent(in) :: io
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: px(geom%nh)

! Local variables
integer :: jo,ii
real(kind_real) :: tmp

! Initialization
px = x

! Previous outer loops
do jo=io,2,-1
   select case (trim(lmp%mode))
   case ('spectral','ritz')
      ! Spectral LMP adjoint
      do ii=1,ni
         px = px+lmp%outer(jo)%ritzvec1(:,ii)*(one/lmp%outer(jo)%eigenval(ii)-one)*sum(lmp%outer(jo)%ritzvec2(:,ii)*x)
      end do

      ! Last Ritz term adjoint
      if (trim(lmp%mode)=='ritz') then
         tmp = zero
         do ii=1,ni
            tmp = tmp+lmp%outer(jo)%omega(ii)*sum(lmp%outer(jo)%ritzvec2(:,ii)*x)
         end do
         do ii=1,ni
            px = px-lmp%outer(jo)%lancvec1(:,ni+1)*lmp%outer(jo)%omega(ii) &
 & *sum(lmp%outer(jo)%ritzvec2(:,ii)*x)-lmp%outer(jo)%ritzvec1(:,ii)*lmp%outer(jo)%omega(ii) &
 & *sum(lmp%outer(jo)%lancvec2(:,ni+1)*x)+lmp%outer(jo)%ritzvec1(:,ii)*lmp%outer(jo)%omega(ii)*tmp
         end do
      end if
   end select
end do

end subroutine lmp_apply_ad

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt
! Purpose: apply LMP square-root
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt(lmp,geom,ni,io,x,px)

implicit none

! Passed variables
class(lmp_type),intent(in) :: lmp
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni
integer,intent(in) :: io
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: px(geom%nh)

! Local variables
integer :: jo,ii
real(kind_real) :: pxtmp(geom%nh)

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
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(one/sqrt(lmp%outer(jo)%eigenval(ii))-one)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
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
subroutine lmp_apply_sqrt_ad(lmp,geom,ni,io,x,px)

implicit none

! Passed variables
class(lmp_type),intent(in) :: lmp
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni
integer,intent(in) :: io
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: px(geom%nh)

! Local variables
integer :: jo,ii
real(kind_real) :: pxtmp(geom%nh)

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
         pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ii)*(one/sqrt(lmp%outer(jo)%eigenval(ii))-one)*sum(lmp%outer(jo)%ritzvec(:,ii)*px)
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
