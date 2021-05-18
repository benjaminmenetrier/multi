!----------------------------------------------------------------------
! Module: tools_netcdf
! Purpose: NetCDF tools
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module tools_netcdf

use netcdf

implicit none

contains

!----------------------------------------------------------------------
! Subroutine: ncerr
! Purpose: handle NetCDF error
!----------------------------------------------------------------------
subroutine ncerr(subr,info)

implicit none

! Passed variables
character(len=*),intent(in) :: subr !< Calling subroutine
integer,intent(in) :: info          !< Info index

! Check status
if (info/=nf90_noerr) then
   write(*,'(a)') 'NetCDF error in '//trim(subr)//': '//nf90_strerror(info)
   stop
end if

end subroutine ncerr

end module tools_netcdf
