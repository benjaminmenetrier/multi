!----------------------------------------------------------------------
! Module: tools_rand
! Purpose: random number generator
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module tools_rand

use tools_kinds
use tools_const

implicit none

contains

!----------------------------------------------------------------------
! Subroutine: set_seed
! Purpose: set random seed
!----------------------------------------------------------------------
subroutine set_seed(rand_seed)

implicit none

! Passed variables
logical,intent(in) :: rand_seed

! Local variable
integer :: nseed,seed_offset,iseed
integer,allocatable :: seed(:)

! Get seed size
call random_seed(size=nseed)

! Allocation
allocate(seed(nseed))

! Define seed offset
if (rand_seed) then
   call system_clock(count=seed_offset)
else
   seed_offset = 0
end if

! Define seed
do iseed=1,nseed
   seed(iseed) = seed_offset+iseed
end do

! Set seed
call random_seed(put=seed)

end subroutine set_seed

!----------------------------------------------------------------------
! Subroutine: rand_normal
! Purpose: random number with a Normal distribution
!----------------------------------------------------------------------
subroutine rand_normal(n,rand)

implicit none

! Passed variables
integer,intent(in) :: n
real(kind_real),intent(out) :: rand(n)

! Local variables
integer :: iset,i
real(kind_real) :: rsq,v1,v2,fac,gset,gasdev

! Normal distribution
iset = 0
do i=1,n
   if(iset==0) then
      rsq = zero
      do while((rsq>=one).or.(rsq<=zero))
         call random_number(v1)
         v1=two*v1-one
         call random_number(v2)
         v2=two*v2-one
         rsq=v1**2+v2**2
      end do
      fac = sqrt(-two*log(rsq)/rsq)
      gset = v1*fac
      gasdev = v2*fac
      iset = 1
   else
      gasdev = gset
      iset = 0
   end if
   rand(i) = gasdev
end do

end subroutine rand_normal

end module tools_rand
