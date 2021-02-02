!----------------------------------------------------------------------
! Module: type_geom
! Purpose: geometry
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_geom

use, intrinsic :: iso_c_binding
use netcdf
use tools_const
use tools_netcdf

implicit none

include 'fftw3.f03'


type geom_type
   integer             :: nx
   integer             :: ny
   integer             :: nh
   integer             :: kmax
   integer             :: lmax
   real(8),allocatable :: x(:)
   real(8),allocatable :: y(:)
   character(len=1024) :: interp_method
contains
   procedure :: setup => geom_setup
   procedure :: write => geom_write
   procedure :: complex_to_real => geom_complex_to_real
   procedure :: real_to_complex => geom_real_to_complex
   procedure :: complex_to_full => geom_complex_to_full
   procedure :: full_to_complex => geom_full_to_complex
   procedure :: gp2sp => geom_gp2sp
   procedure :: sp2gp => geom_sp2gp
   procedure :: spdotprod => geom_spdotprod
   procedure :: interp => geom_interp
   procedure :: interp_spectral => geom_interp_spectral
   procedure :: interp_bilinear => geom_interp_bilinear
   procedure :: interp_bilinear_ad => geom_interp_bilinear_ad
   procedure :: interp_nearest => geom_interp_nearest
   procedure :: transitive_interp_test => geom_transitive_interp_test
end type geom_type

contains

!----------------------------------------------------------------------
! Subroutine: geom_setup
! Purpose: setup geometry
!----------------------------------------------------------------------
subroutine geom_setup(geom,nx,ny,interp_method)
    
implicit none

! Passed variables
class(geom_type),intent(inout) :: geom
integer,intent(in) :: nx
integer,intent(in) :: ny
character(len=1024):: interp_method

! Local variables
integer :: ix,iy
real(8) :: dp1,dp2
real(8),allocatable :: gpsave(:),gp(:),gp1(:),gp2(:)
real(8),allocatable :: spsave(:),sp(:),sp1(:),sp2(:)
complex(8),allocatable :: cpsave(:,:),cp(:,:),cp_fullsave(:,:),cp_full(:,:)

! Check dimensions
if (mod(nx,2)==0) then
!   write(*,'(a)') '      Error: nx should be odd'
!   stop
end if
if (mod(ny,2)==0) then
!   write(*,'(a)') '      Error: ny should be odd'
!  stop
end if

! Copy dimensions and attributes
geom%nx = nx
geom%ny = ny
geom%interp_method = interp_method

! Derived dimensions
geom%nh = geom%nx*geom%ny
geom%kmax = geom%nx/2
geom%lmax = geom%ny/2

! Allocation
allocate(geom%x(geom%nx))
allocate(geom%y(geom%ny))
allocate(gpsave(geom%nh))
allocate(gp(geom%nh))
allocate(gp1(geom%nh))
allocate(gp2(geom%nh))
allocate(spsave(geom%nh))
allocate(sp(geom%nh))
allocate(sp1(geom%nh))
allocate(sp2(geom%nh))
allocate(cpsave(geom%kmax+1,geom%ny))
allocate(cp(geom%kmax+1,geom%ny))
allocate(cp_fullsave(0:geom%kmax,-geom%lmax:geom%lmax))
allocate(cp_full(0:geom%kmax,-geom%lmax:geom%lmax))

! Grid coordinates
do ix=1,geom%nx
   geom%x(ix) = real(ix-1,8)/real(geom%nx,8)
end do
do iy=1,geom%ny
   geom%y(iy) = real(iy-1,8)/real(geom%ny,8)
end do

write(*,'(a,i4,a,i4)') '      Grid: ',geom%nx,' x ',geom%ny
write(*,'(a,i4,a,i4)') '      Spectral: ',geom%kmax,' x ',geom%lmax

! Initialization
call random_number(gpsave)
call geom%gp2sp(gpsave,spsave)
call geom%real_to_complex(spsave,cpsave)
call geom%complex_to_full(cp,cp_fullsave)

! Complex to real to complex
cp = cpsave
call geom%complex_to_real(cp,sp)
call geom%real_to_complex(sp,cp)
write(*,'(a,e15.8)') '      Complex to real to complex: ',maxval(abs(cp-cpsave))

! Real to complex to real
sp = spsave
call geom%real_to_complex(sp,cp)
call geom%complex_to_real(cp,sp)
write(*,'(a,e15.8)') '      Real to complex to real:    ',maxval(abs(sp-spsave))

! Complex to full to complex
cp = cpsave
call geom%complex_to_full(cp,cp_full)
call geom%full_to_complex(cp_full,cp)
write(*,'(a,e15.8)') '      Complex to full to complex: ',maxval(abs(cp-cpsave))

! Full to complex to full
cp_full = cp_fullsave
call geom%full_to_complex(cp_full,cp)
call geom%complex_to_full(cp,cp_full)
write(*,'(a,e15.8)') '      Full to complex to full:    ',maxval(abs(cp_full-cp_fullsave))

! GP to SP to GP
gp = gpsave
call geom%gp2sp(gp,sp)
call geom%sp2gp(sp,gp)
write(*,'(a,e15.8)') '      GP to SP to GP:             ',maxval(abs(gp-gpsave))

! SP to GP to SP
sp = spsave
call geom%sp2gp(sp,gp)
call geom%gp2sp(gp,sp)
write(*,'(a,e15.8)') '      SP to GP to SP:             ',maxval(abs(sp-spsave))

! Dot product
call random_number(gp1)
call random_number(gp2)
call geom%gp2sp(gp1,sp1)
call geom%gp2sp(gp2,sp2)
dp1 = sum(gp1*gp1)
dp2 = geom%spdotprod(sp1,sp1)
write(*,'(a,e15.8)') '      Dot product:                ',2.0*abs(dp1-dp2)/abs(dp1+dp2)

! SP to GP adjoint
call geom%sp2gp(sp2,gp2)
dp1 = sum(gp1*gp2)
dp2 = geom%spdotprod(sp1,sp2)
write(*,'(a,e15.8)') '      SP to GP adjoint:           ',2.0*abs(dp1-dp2)/abs(dp1+dp2)

! Null interpolation
call geom%interp(geom,gp1,gp2)
write(*,'(a,e15.8)') '      Null interpolation:         ',maxval(abs(gp1-gp2))

end subroutine geom_setup

!----------------------------------------------------------------------
! Subroutine: geom_write
! Purpose: write geometry
!----------------------------------------------------------------------
subroutine geom_write(geom,grpid)

implicit none

! Passed variables
class(geom_type),intent(inout) :: geom
integer,intent(in) :: grpid

! Local variables
integer :: nx_id,ny_id,x_id,y_id

! Create dimensions
call ncerr('geom_write',nf90_def_dim(grpid,'nx',geom%nx,nx_id))
call ncerr('geom_write',nf90_def_dim(grpid,'ny',geom%ny,ny_id))

! Create variables
call ncerr('geom_write',nf90_def_var(grpid,'x_coord',nf90_double,(/nx_id/),x_id))
call ncerr('geom_write',nf90_def_var(grpid,'y_coord',nf90_double,(/ny_id/),y_id))

! Write variables
call ncerr('geom_write',nf90_put_var(grpid,x_id,geom%x))
call ncerr('geom_write',nf90_put_var(grpid,y_id,geom%y))

end subroutine geom_write

!----------------------------------------------------------------------
! Subroutine: geom_complex_to_real
! Purpose: complex to real mapping
!----------------------------------------------------------------------
subroutine geom_complex_to_real(geom,cp,sp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
complex(8),intent(in) :: cp(geom%kmax+1,geom%ny)
real(8),intent(out) :: sp(geom%nh)

! Local variables
integer :: i,j,isp

! Complex to real
isp = 0
i = 1
j = 1
isp = isp+1
sp(isp) = real(cp(i,j),8)
do j=2,geom%ny/2+1
   isp = isp+1
   sp(isp) = real(cp(i,j),8)
   isp = isp+1
   sp(isp) = aimag(cp(i,j))
end do
do i=2,geom%nx/2+1
   do j=1,geom%ny
      isp = isp+1
      sp(isp) = real(cp(i,j),8)
      isp = isp+1
      sp(isp) = aimag(cp(i,j))
   end do
end do

end subroutine geom_complex_to_real

!----------------------------------------------------------------------
! Subroutine: geom_real_to_complex
! Purpose: real to complex mapping
!----------------------------------------------------------------------
subroutine geom_real_to_complex(geom,sp,cp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
real(8),intent(in) :: sp(geom%nh)
complex(8),intent(out) :: cp(geom%kmax+1,geom%ny)

! Local variables
integer :: i,j,isp

! Real to complex
isp = 0
i = 1
j = 1
isp = isp+1
cp(i,j) = dcmplx(sp(isp),0.0)
do j=2,geom%ny/2+1
   isp = isp+1
   cp(i,j) = dcmplx(sp(isp),sp(isp+1))
   cp(i,geom%ny-j+2) = dconjg(cp(i,j))
   isp = isp+1
end do
do i=2,geom%nx/2+1
   do j=1,geom%ny
      isp = isp+1
      cp(i,j) = dcmplx(sp(isp),sp(isp+1))
      isp = isp+1
   end do
end do

end subroutine geom_real_to_complex

!----------------------------------------------------------------------
! Subroutine: geom_complex_to_full
! Purpose: complex array mapping
!----------------------------------------------------------------------
subroutine geom_complex_to_full(geom,cp,cp_full)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
complex(8),intent(in) :: cp(geom%kmax+1,geom%ny)
complex(8),intent(out) :: cp_full(0:geom%kmax,-geom%lmax:geom%lmax)

! Complex array mapping
cp_full(0:geom%kmax,0:geom%lmax) = cp(1:geom%kmax+1,1:geom%lmax+1)
cp_full(1:geom%kmax,-geom%lmax:-1) = cp(2:geom%kmax+1,geom%ny-geom%lmax+1:geom%ny)
cp_full(0,-1:-geom%lmax:-1) = dcmplx(0.0,0.0)

end subroutine geom_complex_to_full

!----------------------------------------------------------------------
! Subroutine: geom_full_to_complex
! Purpose: complex array mapping, inverse
!----------------------------------------------------------------------
subroutine geom_full_to_complex(geom,cp_full,cp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
complex(8),intent(in) :: cp_full(0:geom%kmax,-geom%lmax:geom%lmax)
complex(8),intent(out) :: cp(geom%kmax+1,geom%ny)

! Complex array mapping, inverse
cp(1:geom%kmax+1,1:geom%lmax+1) = cp_full(0:geom%kmax,0:geom%lmax)
cp(2:geom%kmax+1,geom%ny-geom%lmax+1:geom%ny) = cp_full(1:geom%kmax,-geom%lmax:-1)

end subroutine geom_full_to_complex

!----------------------------------------------------------------------
! Subroutine: geom_gp2sp
! Purpose: grid-point to spectral
!----------------------------------------------------------------------
subroutine geom_gp2sp(geom,gp,sp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
real(8),intent(in) :: gp(geom%nh)
real(8),intent(out) :: sp(geom%nh)

! Local variables
integer(8) :: plan
real(8) :: norm
real(8) :: gp_2d(geom%nx,geom%ny),gp_2d_save(geom%nx,geom%ny)
complex(8) :: cp(geom%kmax+1,geom%ny)

! Reshape vector
gp_2d = reshape(gp,(/geom%nx,geom%ny/))

! Create plan
call dfftw_plan_dft_r2c_2d(plan,geom%nx,geom%ny,gp_2d,cp,fftw_estimate)

! Compute FFT
cp = dcmplx(0.0,0.0)
call dfftw_execute_dft_r2c(plan,gp_2d,cp)

! Auto-adjoint factor
cp(1,2:geom%ny) = cp(1,2:geom%ny)*sqrt(2.0_8)
cp(2:geom%kmax+1,:) = cp(2:geom%kmax+1,:)*sqrt(2.0_8)

! Normalize
norm = 1.0/sqrt(real(geom%nh,8))
cp = cp*norm

! Complex to real
call geom%complex_to_real(cp,sp)

! Destroy plan
call dfftw_destroy_plan(plan)

end subroutine geom_gp2sp

!----------------------------------------------------------------------
! Subroutine: geom_sp2gp
! Purpose: spectral to grid-point
!----------------------------------------------------------------------
subroutine geom_sp2gp(geom,sp,gp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
real(8),intent(in) :: sp(geom%nh)
real(8),intent(out) :: gp(geom%nh)

! Local variables
integer(8) :: plan
real(8) :: gp_2d(geom%nx,geom%ny),norm
complex(8) :: cp(geom%kmax+1,geom%ny)

! Real to complex
call geom%real_to_complex(sp,cp)

! Auto-adjoint factor
cp(1,2:geom%ny) = cp(1,2:geom%ny)/sqrt(2.0_8)
cp(2:geom%kmax+1,:) = cp(2:geom%kmax+1,:)/sqrt(2.0_8)

! Create plan
call dfftw_plan_dft_c2r_2d(plan,geom%nx,geom%ny,cp,gp_2d,fftw_estimate)

! Compute FFT
gp_2d = 0.0
call dfftw_execute_dft_c2r(plan,cp,gp_2d)

! Normalize
norm = 1.0/sqrt(real(geom%nh,8))
gp_2d = gp_2d*norm

! Reshape vector
gp = pack(gp_2d,.true.)

! Destroy plan
call dfftw_destroy_plan(plan)

end subroutine geom_sp2gp

!----------------------------------------------------------------------
! Function: geom_spdotprod
! Purpose: spectral dot product
!----------------------------------------------------------------------
function geom_spdotprod(geom,sp1,sp2) result(dp)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom
real(8),intent(in) :: sp1(geom%nh)
real(8),intent(in) :: sp2(geom%nh)

! Returned variable
real(8) :: dp

! Local variables
complex(8) :: cp1(geom%kmax+1,geom%ny),cp2(geom%kmax+1,geom%ny)
complex(8) :: cp1_full(0:geom%kmax,-geom%lmax:geom%lmax),cp2_full(0:geom%kmax,-geom%lmax:geom%lmax)

! Real to complex
call geom%real_to_complex(sp1,cp1)
call geom%real_to_complex(sp2,cp2)

! Remap complex array
call geom%complex_to_full(cp1,cp1_full)
call geom%complex_to_full(cp2,cp2_full)

! Dot product
dp = real(sum(cp1_full*dconjg(cp2_full)),8)

end function geom_spdotprod

!----------------------------------------------------------------------
! Subroutine: geom_interp
! Purpose: field interpolation interface
!----------------------------------------------------------------------
subroutine geom_interp(geom_in,geom_out,gp_in,gp_out)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom_in
type(geom_type),intent(in) :: geom_out
real(8),intent(in) :: gp_in(geom_in%nh)
real(8),intent(out) :: gp_out(geom_out%nh)

! Local variables
integer :: ix,iy,ih
real(8),allocatable :: sp_in(:),sp_out(:)

select case (trim(geom_in%interp_method))
case('spectral')
   ! Allocation
   allocate(sp_in(geom_in%nh))
   allocate(sp_out(geom_out%nh))

   ! Direct Fourier transform
   call geom_in%gp2sp(gp_in,sp_in)

   ! Spectral interpolation
   call geom_in%interp_spectral(geom_out,sp_in,sp_out)

   ! Inverse Fourier transform
   call geom_out%sp2gp(sp_out,gp_out)

   ! Release memory
   deallocate(sp_in)
   deallocate(sp_out)
case ('bilinear')
   ! Initialization
   ih = 0
   do iy=1,geom_out%ny
      do ix=1,geom_out%nx
         ! Apply scalar bilinear interpolation
         ih = ih+1
         call geom_in%interp_bilinear(geom_out%x(ix),geom_out%y(iy),gp_in,gp_out(ih))
      end do
   end do
case ('nearest')
   ih=0
   do iy=1,geom_out%ny
      do ix=1,geom_out%nx
         ih = ih+1
         ! Apply scalar nearest neighbor interpolation
         call geom_in%interp_nearest(geom_out%x(ix),geom_out%y(iy),gp_in,gp_out(ih))
      end do
   end do
case default
   write(*,'(a)') 'Error: wrong interpolation method'
  stop
end select

end subroutine geom_interp

!----------------------------------------------------------------------
! Subroutine: geom_interp_spectral
! Purpose: spectral interpolation (zero padding)
!----------------------------------------------------------------------
subroutine geom_interp_spectral(geom_in,geom_out,sp_in,sp_out)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom_in
type(geom_type),intent(in) :: geom_out
real(8),intent(in) :: sp_in(geom_in%nh)
real(8),intent(out) :: sp_out(geom_out%nh)

! Local variables
integer :: kmaxmin,lmaxmin
complex(8),allocatable :: cp_in(:,:),cp_in_full(:,:)
complex(8),allocatable :: cp_out(:,:),cp_out_full(:,:)

if ((geom_in%nx==geom_out%nx).and.(geom_in%ny==geom_out%ny)) then
   sp_out = sp_in
else
   ! Allocation
   allocate(cp_in(geom_in%kmax+1,geom_in%ny))
   allocate(cp_in_full(0:geom_in%kmax,-geom_in%lmax:geom_in%lmax))
   allocate(cp_out(geom_out%kmax+1,geom_out%ny))
   allocate(cp_out_full(0:geom_out%kmax,-geom_out%lmax:geom_out%lmax))

   ! Real to complex
   call geom_in%real_to_complex(sp_in,cp_in)

   ! Complex to full
   call geom_in%complex_to_full(cp_in,cp_in_full)

   ! Initialize
   cp_out_full = dcmplx(0.0,0.0)

   ! Copy
   kmaxmin = min(geom_in%kmax,geom_out%kmax)
   lmaxmin = min(geom_in%lmax,geom_out%lmax)
   cp_out_full(0:kmaxmin,-lmaxmin:lmaxmin) = cp_in_full(0:kmaxmin,-lmaxmin:lmaxmin)

   ! Full to complex
   call geom_out%full_to_complex(cp_out_full,cp_out)

   ! Complex to real
   call geom_out%complex_to_real(cp_out,sp_out)

   ! Release memory
   deallocate(cp_in)
   deallocate(cp_in_full)
   deallocate(cp_out)
   deallocate(cp_out_full)
end if

end subroutine geom_interp_spectral

!----------------------------------------------------------------------
! Subroutine: geom_interp_bilinear
! Purpose: bilinar interpolation
!----------------------------------------------------------------------
subroutine geom_interp_bilinear(geom_in,x_out,y_out,gp_in,gp_out)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom_in
real(8),intent(in) :: x_out
real(8),intent(in) :: y_out
real(8),intent(in) :: gp_in(geom_in%nh)
real(8),intent(out) :: gp_out

! Local variables
integer :: ix_inf,ix_sup,iy_inf,iy_sup
integer :: ih_inf_inf,ih_inf_sup,ih_sup_inf,ih_sup_sup
real(8) :: rx_inf,rx_sup,ry_inf,ry_sup

! Bilinear interpolation indices and coefficients
ix_inf = floor(x_out*real(geom_in%nx,8))+1
iy_inf = floor(y_out*real(geom_in%ny,8))+1
if (ix_inf<geom_in%nx) then
   ix_sup = ix_inf+1
   rx_inf = (geom_in%x(ix_sup)-x_out)/(geom_in%x(ix_sup)-geom_in%x(ix_inf))
   rx_sup = (x_out-geom_in%x(ix_inf))/(geom_in%x(ix_sup)-geom_in%x(ix_inf))
else
   ix_sup = 1
   rx_inf = (1.0-x_out)/(1.0-geom_in%x(ix_inf))
   rx_sup = (x_out-geom_in%x(ix_inf))/(1.0-geom_in%x(ix_inf))
end if
if (iy_inf<geom_in%ny) then
   iy_sup = iy_inf+1
   ry_inf = (geom_in%y(iy_sup)-y_out)/(geom_in%y(iy_sup)-geom_in%y(iy_inf))
   ry_sup = (y_out-geom_in%y(iy_inf))/(geom_in%y(iy_sup)-geom_in%y(iy_inf))
else
   iy_sup = 1
   ry_inf = (1.0-y_out)/(1.0-geom_in%y(iy_inf))
   ry_sup = (y_out-geom_in%y(iy_inf))/(1.0-geom_in%y(iy_inf))
end if

! Apply bilinear interpolation
ih_inf_inf = ix_inf+(iy_inf-1)*geom_in%nx
ih_inf_sup = ix_inf+(iy_sup-1)*geom_in%nx
ih_sup_inf = ix_sup+(iy_inf-1)*geom_in%nx
ih_sup_sup = ix_sup+(iy_sup-1)*geom_in%nx
gp_out = rx_inf*ry_inf*gp_in(ih_inf_inf) &
     & + rx_inf*ry_sup*gp_in(ih_inf_sup) &
     & + rx_sup*ry_inf*gp_in(ih_sup_inf) &
     & + rx_sup*ry_sup*gp_in(ih_sup_sup)

end subroutine geom_interp_bilinear

!----------------------------------------------------------------------
! Subroutine: geom_interp_bilinear_ad
! Purpose: bilinar interpolation, adjoint
!----------------------------------------------------------------------
subroutine geom_interp_bilinear_ad(geom_in,x_out,y_out,gp_out,gp_in)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom_in
real(8),intent(in) :: x_out
real(8),intent(in) :: y_out
real(8),intent(in) :: gp_out
real(8),intent(inout) :: gp_in(geom_in%nh)

! Local variables
integer :: ix_inf,ix_sup,iy_inf,iy_sup
integer :: ih_inf_inf,ih_inf_sup,ih_sup_inf,ih_sup_sup
real(8) :: rx_inf,rx_sup,ry_inf,ry_sup

! Bilinear interpolation indices and coefficients
ix_inf = floor(x_out*real(geom_in%nx,8))+1
iy_inf = floor(y_out*real(geom_in%ny,8))+1
if (ix_inf<geom_in%nx) then
   ix_sup = ix_inf+1
   rx_inf = (geom_in%x(ix_sup)-x_out)/(geom_in%x(ix_sup)-geom_in%x(ix_inf))
   rx_sup = (x_out-geom_in%x(ix_inf))/(geom_in%x(ix_sup)-geom_in%x(ix_inf))
else
   ix_sup = 1
   rx_inf = (1.0-x_out)/(1.0-geom_in%x(ix_inf))
   rx_sup = (x_out-geom_in%x(ix_inf))/(1.0-geom_in%x(ix_inf))
end if
if (iy_inf<geom_in%ny) then
   iy_sup = iy_inf+1
   ry_inf = (geom_in%y(iy_sup)-y_out)/(geom_in%y(iy_sup)-geom_in%y(iy_inf))
   ry_sup = (y_out-geom_in%y(iy_inf))/(geom_in%y(iy_sup)-geom_in%y(iy_inf))
else
   iy_sup = 1
   ry_inf = (1.0-y_out)/(1.0-geom_in%y(iy_inf))
   ry_sup = (y_out-geom_in%y(iy_inf))/(1.0-geom_in%y(iy_inf))
end if

! Apply bilinear interpolation adjoint
ih_inf_inf = ix_inf+(iy_inf-1)*geom_in%nx
ih_inf_sup = ix_inf+(iy_sup-1)*geom_in%nx
ih_sup_inf = ix_sup+(iy_inf-1)*geom_in%nx
ih_sup_sup = ix_sup+(iy_sup-1)*geom_in%nx
gp_in(ih_inf_inf) = gp_in(ih_inf_inf)+rx_inf*ry_inf*gp_out
gp_in(ih_inf_sup) = gp_in(ih_inf_sup)+rx_inf*ry_sup*gp_out
gp_in(ih_sup_inf) = gp_in(ih_sup_inf)+rx_sup*ry_inf*gp_out
gp_in(ih_sup_sup) = gp_in(ih_sup_sup)+rx_sup*ry_sup*gp_out

end subroutine geom_interp_bilinear_ad

!----------------------------------------------------------------------
! Subroutine: geom_interp_nearest
! Purpose: nearest neighbor interpolation
!----------------------------------------------------------------------
subroutine geom_interp_nearest(geom_in,x_out,y_out,gp_in,gp_out)

implicit none

! Passed variables
class(geom_type),intent(in) :: geom_in
real(8),intent(in) :: x_out
real(8),intent(in) :: y_out
real(8),intent(in) :: gp_in(geom_in%nh)
real(8),intent(out) :: gp_out

! Local variables
integer :: ix_inf,ix_sup,iy_inf,iy_sup
integer :: ix,iy,ih
real(8) :: dist_min,distx_inf,distx_sup,disty_inf,disty_sup

! Interpolation indices and coefficients
ix_inf = floor(x_out*real(geom_in%nx,8))+1
iy_inf = floor(y_out*real(geom_in%ny,8))+1
if (ix_inf<geom_in%nx) then
   ix_sup = ix_inf+1
   distx_inf = (x_out-geom_in%x(ix_inf))**2
   distx_sup = (geom_in%x(ix_sup)-x_out)**2
else
   ix_sup = 1
   distx_inf = (x_out-geom_in%x(ix_inf))**2
   distx_sup = (1.0-x_out)**2
end if
if (iy_inf<geom_in%ny) then
   iy_sup = iy_inf+1
   disty_inf = (y_out-geom_in%y(iy_inf))**2
   disty_sup = (geom_in%y(iy_sup)-y_out)**2
else
   iy_sup = 1
   disty_inf = (y_out-geom_in%y(iy_inf))**2
   disty_sup = (1.0-y_out)**2
end if

! Find minimum (squared) distance
dist_min = huge(1.0)
if (distx_inf+disty_inf<dist_min) then
   dist_min = distx_inf+disty_inf
   ix = ix_inf
   iy = iy_inf
end if
if (distx_inf+disty_sup<dist_min) then
   dist_min = distx_inf+disty_sup
   ix = ix_inf
   iy = iy_sup
end if
if (distx_sup+disty_inf<dist_min) then
   dist_min = distx_sup+disty_inf
   ix = ix_sup
   iy = iy_inf
end if
if (distx_sup+disty_sup<dist_min) then
   dist_min = distx_sup+disty_sup
   ix = ix_sup
   iy = iy_sup
end if

! Copy nearest neighbor value
ih = ix+(iy-1)*geom_in%nx
gp_out = gp_in(ih)

end subroutine geom_interp_nearest

!----------------------------------------------------------------------
! Subroutine: geom_transitive_interp_test
! Purpose: test interpolation transitivity
!----------------------------------------------------------------------
subroutine geom_transitive_interp_test(geom,geom_full)

implicit none
  
! Passed variables
class(geom_type),intent(in) :: geom
type(geom_type),intent(in) :: geom_full

! Local variables
integer :: ix,iy,ih
integer :: nx_inter,ny_inter
real(8) :: test_upscaling,test_downscaling,test_right_inverse
real(8),allocatable :: x(:),x_save(:),x_full(:),x_full_save(:),x_inter(:)
logical :: is_transitive
type(geom_type)  :: geom_inter

write(*,'(a,i1)') '   Geometry setup for intermediate resolution (transitivity test)'

! Setup intermediate geometry
nx_inter = (geom%nx+geom_full%nx)/2
if (mod(nx_inter,2)==0) nx_inter = nx_inter+1
ny_inter = (geom%ny+geom_full%ny)/2
if (mod(ny_inter,2)==0) ny_inter = ny_inter+1
call geom_inter%setup(nx_inter,ny_inter,geom%interp_method)

! Allocation
allocate(x(geom%nh))
allocate(x_save(geom%nh))
allocate(x_full(geom_full%nh))
allocate(x_full_save(geom_full%nh))
allocate(x_inter(geom_inter%nh))

!  Upscaling test
do iy=1,geom%ny
   do ix=1,geom%nx
      ih = ix+(iy-1)*geom%nx
      x(ih) = geom%x(ix)
   end do
end do
!call random_number(x)
call geom%interp(geom_full,x,x_full_save)
call geom%interp(geom_inter,x,x_inter)
call geom_inter%interp(geom_full,x_inter,x_full)
test_upscaling = maxval(abs(x_full-x_full_save))/maxval(abs(x))

!  Downscaling test
call random_number(x_full)
call geom_full%interp(geom,x_full,x_save)
call geom_full%interp(geom_inter,x_full,x_inter)
call geom_inter%interp(geom,x_inter,x)
test_downscaling = maxval(abs(x-x_save))/maxval(abs(x_full))

! Right-inverse test
call random_number(x)
x_save = x
call geom%interp(geom_full,x,x_full)
call geom_full%interp(geom,x_full,x)
test_right_inverse = maxval(abs(x-x_save))/maxval(abs(x_save))
!print*, x
!print*, x_save
!print*, x-x_save
!print*

! Print results
is_transitive = (test_upscaling<1.0e-12).and.(test_downscaling<1.0e-12).and.(test_right_inverse<1.0e-12)
write(*,'(a,l)') '      Transitive interpolation: ', is_transitive
write(*,'(a,e15.8)') '       - Upscaling:               ',test_upscaling
write(*,'(a,e15.8)') '       - Downscaling:             ',test_downscaling
write(*,'(a,e15.8)') '       - Right-inverse:           ',test_right_inverse

! Release memory
deallocate(x)
deallocate(x_save)
deallocate(x_full)
deallocate(x_full_save)
deallocate(x_inter)

end subroutine geom_transitive_interp_test

end module type_geom
