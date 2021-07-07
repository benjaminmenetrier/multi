!----------------------------------------------------------------------
! Module: tools_const
! Purpose: constants
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module tools_const

use tools_kinds

implicit none

real(kind_real),parameter :: zero = 0.0_kind_real
real(kind_real),parameter :: one = 1.0_kind_real
real(kind_real),parameter :: two = 2.0_kind_real
real(kind_real),parameter :: pi = acos(-one)

! Missing values
integer,parameter :: msv_int = -999
real(kind_real),parameter :: msv_real = -999.0_kind_real

end module tools_const
