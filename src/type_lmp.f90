!----------------------------------------------------------------------
! Module: type_lmp
! Purpose: LMP preconditioner
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_lmp

implicit none

type lmp_type
   integer :: np
   real(8),allocatable :: ritzvec_trunc(:,:)
   real(8),allocatable :: ritzvec(:,:)
   real(8),allocatable :: vec1(:,:)
   real(8),allocatable :: vec2(:,:)
   real(8),allocatable :: val(:)
end type lmp_type

real(8),parameter :: tol = 1.0e-6
integer,parameter :: niter = 50

contains

!----------------------------------------------------------------------
! Subroutine: lmp_alloc
! Purpose: allocate LMP
!----------------------------------------------------------------------
subroutine lmp_alloc(lmp,nn,ni)

implicit none

! Passed variables
type(lmp_type),intent(inout) :: lmp
integer,intent(in) :: nn
integer,intent(in) :: ni

! Release memory
if (allocated(lmp%ritzvec_trunc)) deallocate(lmp%ritzvec_trunc)
if (allocated(lmp%ritzvec)) deallocate(lmp%ritzvec)
if (allocated(lmp%vec1)) deallocate(lmp%vec1)
if (allocated(lmp%vec2)) deallocate(lmp%vec2)
if (allocated(lmp%val)) deallocate(lmp%val)

! Allocation
allocate(lmp%ritzvec_trunc(nn,ni))
allocate(lmp%ritzvec(nn,ni))
allocate(lmp%vec1(nn,ni))
allocate(lmp%vec2(nn,ni))
allocate(lmp%val(ni))

end subroutine lmp_alloc

!----------------------------------------------------------------------
! Subroutine: lmp_copy
! Purpose: copy LMP
!----------------------------------------------------------------------
subroutine lmp_copy(lmp_in,lmp_out)

implicit none

! Passed variables
type(lmp_type),intent(in) :: lmp_in
type(lmp_type),intent(out) :: lmp_out

! Local variables
integer :: nn,ni

! Size
nn = size(lmp_in%ritzvec,1)
ni = size(lmp_in%ritzvec,2)

! Allocation
call lmp_alloc(lmp_out,nn,ni)

! Copy
lmp_out%np = lmp_in%np
lmp_out%ritzvec_trunc = lmp_in%ritzvec_trunc
lmp_out%ritzvec = lmp_in%ritzvec
lmp_out%vec1 = lmp_in%vec1
lmp_out%vec2 = lmp_in%vec2
lmp_out%val = lmp_in%val

end subroutine lmp_copy

!----------------------------------------------------------------------
! Subroutine: lmp_apply
! Purpose: apply LMP
!----------------------------------------------------------------------
subroutine lmp_apply(nn,io,lmp,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: io
type(lmp_type),intent(in) :: lmp(io)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip

! Initialization
px = x

! Previous outer loops
do jo=1,io
   do ip=1,lmp(jo)%np
      px = px+lmp(jo)%vec2(:,ip)*(1.0/lmp(jo)%val(ip)-1.0)*sum(lmp(jo)%vec1(:,ip)*x)
   end do
end do

end subroutine lmp_apply

!----------------------------------------------------------------------
! Subroutine: lmp_apply_inv
! Purpose: apply LMP, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_inv(nn,io,lmp,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: io
type(lmp_type),intent(in) :: lmp(io)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)
real(8) :: rmse_rat,xtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=1,io
   ! Initialization
   pxtmp = px

   ! LMP inverse
   do ip=1,lmp(jo)%np
      pxtmp = pxtmp+lmp(jo)%ritzvec(:,ip)*(lmp(jo)%val(ip)-1.0)*sum(lmp(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

! Check RMSE ratio
call lmp_apply(nn,io,lmp,px,xtmp)
rmse_rat = sqrt(sum((x-xtmp)**2)/sum(x**2))
if (rmse_rat>tol) then
   write(*,*) '      ',x(1:5)
   write(*,*) '      ',xtmp(1:5)
   write(*,*) '      LMP inversion RMSE ratio:',rmse_rat
end if

end subroutine lmp_apply_inv

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt
! Purpose: apply LMP square-root
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt(nn,io,lmp,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: io
type(lmp_type),intent(in) :: lmp(io)
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
   do ip=1,lmp(jo)%np
      pxtmp = pxtmp+lmp(jo)%vec1(:,ip)*(1.0/sqrt(lmp(jo)%val(ip))-1.0)*sum(lmp(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_ad
! Purpose: apply LMP square-root adjoint
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_ad(nn,io,lmp,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: io
type(lmp_type),intent(in) :: lmp(io)
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
   do ip=1,lmp(jo)%np
      pxtmp = pxtmp+lmp(jo)%vec1(:,ip)*(1.0/sqrt(lmp(jo)%val(ip))-1.0)*sum(lmp(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

end subroutine lmp_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: lmp_apply_sqrt_inv
! Purpose: apply LMP square-root, inverse
!----------------------------------------------------------------------
subroutine lmp_apply_sqrt_inv(nn,io,lmp,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: io
type(lmp_type),intent(in) :: lmp(io)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,ip
real(8) :: pxtmp(nn)
real(8) :: rmse_rat,xtmp(nn)

! Initialization
px = x

! Previous outer loops
do jo=1,io
   ! Initialization
   pxtmp = px

   ! LMP square-root inverse
   do ip=1,lmp(jo)%np
      pxtmp = pxtmp+lmp(jo)%vec1(:,ip)*(sqrt(lmp(jo)%val(ip))-1.0)*sum(lmp(jo)%vec1(:,ip)*px)
   end do

   ! Update
   px = pxtmp
end do

! Check RMSE ratio
call lmp_apply_sqrt(nn,io,lmp,px,xtmp)
rmse_rat = sqrt(sum((x-xtmp)**2)/sum(x**2))
if (rmse_rat>tol) then
   write(*,*) x(1:5)
   write(*,*) px(1:5)
   write(*,*) xtmp(1:5)
   write(*,*) '      LMP sqrt inversion RMSE ratio:',rmse_rat
   stop
end if

end subroutine lmp_apply_sqrt_inv

end module type_lmp
