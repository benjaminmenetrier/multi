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
   integer :: np
   real(8),allocatable :: ritzvec_trunc(:,:)
   real(8),allocatable :: ritzvec(:,:)
   real(8),allocatable :: vec1(:,:)
   real(8),allocatable :: vec2(:,:)
   real(8),allocatable :: val(:)
end type outer_type

type lmp_type
   integer :: io
   type(outer_type),allocatable :: outer(:)
end type lmp_type

contains

!----------------------------------------------------------------------
! Subroutine: lmp_alloc
! Purpose: allocate LMP
!----------------------------------------------------------------------
subroutine lmp_alloc(lmp,nn,ni,io)

implicit none

! Passed variables
type(lmp_type),intent(inout) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io

! Local variables
integer :: jo

! Release memory
if (allocated(lmp%outer)) then
   do jo=1,lmp%io
      if (allocated(lmp%outer(jo)%ritzvec_trunc)) deallocate(lmp%outer(jo)%ritzvec_trunc)
      if (allocated(lmp%outer(jo)%ritzvec)) deallocate(lmp%outer(jo)%ritzvec)
      if (allocated(lmp%outer(jo)%vec1)) deallocate(lmp%outer(jo)%vec1)
      if (allocated(lmp%outer(jo)%vec2)) deallocate(lmp%outer(jo)%vec2)
      if (allocated(lmp%outer(jo)%val)) deallocate(lmp%outer(jo)%val)
   end do
   deallocate(lmp%outer)
end if

! Allocation
lmp%io = io
allocate(lmp%outer(io))
do jo=1,io
   allocate(lmp%outer(jo)%ritzvec_trunc(nn,ni))
   allocate(lmp%outer(jo)%ritzvec(nn,ni))
   allocate(lmp%outer(jo)%vec1(nn,ni))
   allocate(lmp%outer(jo)%vec2(nn,ni))
   allocate(lmp%outer(jo)%val(ni))
end do

end subroutine lmp_alloc

!----------------------------------------------------------------------
! Subroutine: lmp_apply
! Purpose: apply LMP
!----------------------------------------------------------------------
subroutine lmp_apply(lmp,nn,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip

! Initialization
px = x

! Previous outer loops
do jo=1,io
   do ip=1,lmp%outer(jo)%np
      px = px+lmp%outer(jo)%vec2(:,ip)*(1.0/lmp%outer(jo)%val(ip)-1.0)*sum(lmp%outer(jo)%vec1(:,ip)*x)
   end do
end do

end subroutine lmp_apply

!----------------------------------------------------------------------
! Subroutine: lmp_apply_inv
! Purpose: apply LMP, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_inv(lmp,nn,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=1,io
   ! Initialization
   pxtmp = px

   ! LMP inverse
   do ip=1,lmp%outer(jo)%np
      pxtmp = pxtmp+lmp%outer(jo)%ritzvec(:,ip)*(lmp%outer(jo)%val(ip)-1.0)*sum(lmp%outer(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_inv

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt
! Purpose: apply LMP square-root
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt(lmp,nn,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=io,1,-1
   ! Initialization
   pxtmp = px

   ! LMP square-root
   do ip=1,lmp%outer(jo)%np
      pxtmp = pxtmp+lmp%outer(jo)%vec1(:,ip)*(1.0/sqrt(lmp%outer(jo)%val(ip))-1.0)*sum(lmp%outer(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_ad
! Purpose: apply LMP square-root adjoint
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_ad(lmp,nn,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=1,io
   ! Initialization
   pxtmp = px

   ! LMP square-root
   do ip=1,lmp%outer(jo)%np
      pxtmp = pxtmp+lmp%outer(jo)%vec1(:,ip)*(1.0/sqrt(lmp%outer(jo)%val(ip))-1.0)*sum(lmp%outer(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_inv
! Purpose: apply LMP square-root, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_inv(lmp,nn,io,x,px)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: io
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=1,io
   ! Initialization
   pxtmp = px

   ! LMP square-root inverse
   do ip=1,lmp%outer(jo)%np
      pxtmp = pxtmp+lmp%outer(jo)%vec1(:,ip)*(sqrt(lmp%outer(jo)%val(ip))-1.0)*sum(lmp%outer(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_sqrt_inv

end module type_lmp
