!----------------------------------------------------------------------
! Module: tools_kinds
!> Kinds definition
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module tools_kinds

use, intrinsic :: iso_c_binding
use netcdf

implicit none

! C kinds
integer,parameter :: kind_int = c_int                        !< Integer kind
integer,parameter :: kind_float = c_float                    !< Float kind
integer,parameter :: kind_double = c_double                  !< Double kind

! Real kind alias for the whole code
integer,parameter :: kind_real = c_double                    !< Real kind alias for the whole code

! NetCDF kind alias for the whole code
integer,parameter :: nc_kind_int = nf90_int                  !< NetCDF integer kind alias
integer,parameter :: nc_kind_real = nf90_double              !< NetCDF real kind alias

! Huge values
integer,parameter :: huge_int = huge(0_kind_int)             !< Integer huge
real(kind_real),parameter :: huge_real = huge(0.0_kind_real) !< Real huge

private
public :: kind_float,kind_double,kind_real
public :: nc_kind_int,nc_kind_real
public :: huge_int,huge_real

end module tools_kinds
