!----------------------------------------------------------------------
! Module: type_bmatrix
! Purpose: B matrix
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_bmatrix

use netcdf
use tools_kinds
use tools_netcdf
use tools_rand
use type_geom

implicit none

type bmatrix_type
   real(kind_real),allocatable :: spvar(:,:)
   real(kind_real),allocatable :: sigmab(:)
contains
   procedure :: setup => bmatrix_setup
   procedure :: write => bmatrix_write
   procedure :: apply_spvar => bmatrix_apply_spvar
   procedure :: apply_spvar_sqrt => bmatrix_apply_spvar_sqrt
   procedure :: apply_spvar_sqrt_inv => bmatrix_apply_spvar_sqrt_inv
   procedure :: apply_sqrt => bmatrix_apply_sqrt
   procedure :: apply_sqrt_ad => bmatrix_apply_sqrt_ad
   procedure :: apply_sqrt_inv => bmatrix_apply_sqrt_inv
   procedure :: apply_sqrt_ad_inv => bmatrix_apply_sqrt_ad_inv
   procedure :: apply => bmatrix_apply
   procedure :: apply_inv => bmatrix_apply_inv
   procedure :: randomize => bmatrix_randomize
end type bmatrix_type

contains

!----------------------------------------------------------------------
! Subroutine: bmatrix_setup
! Purpose: setup B matrix
!----------------------------------------------------------------------
subroutine bmatrix_setup(bmatrix,geom,geom_full,sigmabvar,Lb,spvarmin)

implicit none

! Passed variables
class(bmatrix_type),intent(inout) :: bmatrix
type(geom_type),intent(in) :: geom
type(geom_type),intent(in) :: geom_full
real(kind_real),intent(in) :: sigmabvar
real(kind_real),intent(in) :: Lb
real(kind_real),intent(in) :: spvarmin

! Local variables
integer :: k,l,ih,ix,iy
real(kind_real) :: kstarsq,arg,dp1,dp2,dist
real(kind_real) :: x(geom_full%nh),v(geom_full%nh)
real(kind_real) :: x1(geom%nh),x2(geom%nh),x1out(geom%nh),x2out(geom%nh),v1(geom%nh),v2(geom%nh)
real(kind_real) :: x1_2d(geom%nx,geom%ny),x2_2d(geom%nx,geom%ny),sigmab_2d(geom%nx,geom%ny)

! Release memory
if (allocated(bmatrix%sigmab)) deallocate(bmatrix%sigmab)

! Allocation
allocate(bmatrix%spvar(0:geom_full%kmax,-geom_full%lmax:geom_full%lmax))
allocate(bmatrix%sigmab(geom%nh))

! Compute spectral variance
do k=0,geom_full%kmax
   do l=-geom_full%lmax,geom_full%lmax
      ! Compute k*^2
      kstarsq = real(k**2+l**2,kind_real)

      ! Apply Gaussian spectral variance
      arg = two*kstarsq*(pi*Lb)**2
      if (arg<-log(spvarmin)) then
         bmatrix%spvar(k,l) = exp(-arg)
      else
         bmatrix%spvar(k,l) = spvarmin
      end if
   end do
end do

! Compute grid-point variance
ih = 0
do iy=1,geom%ny
   do ix=1,geom%nx
      ih = ih+1
      bmatrix%sigmab(ih) = one+sigmabvar*sin(2*pi*geom%x(ix))*sin(2*pi*geom%y(iy))
   end do
end do

! Compute normalization
x = zero
x(1) = one
call geom_full%gp2sp(x,v)
call bmatrix%apply_spvar(geom_full,v)
call geom_full%sp2gp(v,x)
bmatrix%spvar = bmatrix%spvar/x(1)

! Direct + inverse test on B
call random_number(x1)
call bmatrix%apply(geom,x1,x1out)
call bmatrix%apply_inv(geom,x1out,x2)
write(*,'(a,e15.8)') '      Direct + inverse test on B:    ',maxval(abs(x1-x2))

! Inverse + direct test on B
call random_number(x1)
call bmatrix%apply_inv(geom,x1,x1out)
call bmatrix%apply(geom,x1out,x2)
write(*,'(a,e15.8)') '      Inverse + direct test on B:    ',maxval(abs(x1-x2))

! Test on U being the square-root of B
call random_number(x1)
call bmatrix%apply(geom,x1,x1out)
call bmatrix%apply_sqrt_ad(geom,x1,v1)
call bmatrix%apply_sqrt(geom,v1,x2)
write(*,'(a,e15.8)') '      U = square-root of B:          ',maxval(abs(x1out-x2))

! Adjoint test on U
call random_number(x1)
call random_number(x2)
call geom%gp2sp(x2,v2)
call bmatrix%apply_sqrt_ad(geom,x1,v1)
call bmatrix%apply_sqrt(geom,v2,x2)
dp1 = sum(x1*x2)
dp2 = geom%spdotprod(v1,v2)
write(*,'(a,e15.8)') '      Adjoint test on U:             ',two*abs(dp1-dp2)/abs(dp1+dp2)

! Adjoint test on B
call random_number(x1)
call random_number(x2)
call bmatrix%apply(geom,x1,x1out)
call bmatrix%apply(geom,x2,x2out)
dp1 = sum(x1*x2out)
dp2 = sum(x2*x1out)
write(*,'(a,e15.8)') '      Auto-adjoint test on B:        ',two*abs(dp1-dp2)/abs(dp1+dp2)

! Print correlation shape
x1_2d = zero
x1_2d(1,1) = one
x1 = pack(x1_2d,.true.)
call bmatrix%apply(geom,x1,x2)
x2_2d = reshape(x2,(/geom%nx,geom%ny/))
sigmab_2d = reshape(bmatrix%sigmab,(/geom%nx,geom%ny/))
x2_2d = x2_2d/(sigmab_2d(1,1)*sigmab_2d)

iy = 1
write(*,'(a)',advance='no') '      X-coordinate:                 '
do ix=1,geom%nx
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   write(*,'(f6.2)',advance='no') geom%x(ix)
end do
write(*,'(a)')
write(*,'(a)',advance='no') '      Theoretical correlation shape:'
do ix=1,geom%nx
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   if (ix<geom%nx/2) then
      dist = geom%x(ix)
   else
      dist = one-geom%x(ix)
   end if
   write(*,'(f6.2)',advance='no') exp(-0.5*(dist/Lb)**2)
end do
write(*,'(a)')
write(*,'(a)',advance='no') '      Effective correlation shape:  '
do ix=1,geom%nx
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   write(*,'(f6.2)',advance='no') x2_2d(ix,iy)
end do
write(*,'(a)')
write(*,'(a)',advance='no') '      Y-coordinate:                 '
ix = 1
do iy=1,geom%ny
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   write(*,'(f6.2)',advance='no') geom%y(iy)
end do
write(*,'(a)')
write(*,'(a)',advance='no') '      Theoretical correlation shape:'
do iy=1,geom%ny
   ih = 1+(iy-1)*geom%nx
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   if (iy<geom%ny/2) then
      dist = geom%y(iy)
   else
      dist = one-geom%y(iy)
   end if
   write(*,'(f6.2)',advance='no') exp(-0.5*(dist/Lb)**2)
end do
write(*,'(a)')
write(*,'(a)',advance='no') '      Effective correlation shape:  '
do iy=1,geom%ny
   ih = 1+(iy-1)*geom%nx
   if (x2_2d(ix,iy)<1.0e-3_kind_real) exit
   write(*,'(f6.2)',advance='no') x2_2d(ix,iy)
end do
write(*,'(a)')

end subroutine bmatrix_setup

!----------------------------------------------------------------------
! Subroutine: bmatrix_write
! Purpose: write B matrix
!----------------------------------------------------------------------
subroutine bmatrix_write(bmatrix,geom,grpid)

implicit none

! Passed variables
class(bmatrix_type),intent(inout) :: bmatrix
type(geom_type),intent(in) :: geom
integer,intent(in) :: grpid

! Local variables
integer :: nx_id,ny_id,sigmab_id,dirac_cov_id,dirac_cor_id
real(kind_real) :: x(geom%nh),Bx(geom%nh)
real(kind_real) :: x1_2d(geom%nx,geom%ny),Bx_2d(geom%nx,geom%ny),sigmab_2d(geom%nx,geom%ny),Cx_2d(geom%nx,geom%ny)

! Dirac test at coordinates (1,1) 
x1_2d = zero
x1_2d(1,1) = one
x = pack(x1_2d,.true.)
call bmatrix%apply(geom,x,Bx)
Bx_2d = reshape(Bx,(/geom%nx,geom%ny/))
sigmab_2d = reshape(bmatrix%sigmab,(/geom%nx,geom%ny/))
Cx_2d = Bx_2d/(sigmab_2d(1,1)*sigmab_2d)

! Get dimensions
call ncerr('bmatrix_write',nf90_inq_dimid(grpid,'nx',nx_id))
call ncerr('bmatrix_write',nf90_inq_dimid(grpid,'ny',ny_id))

! Create variables
call ncerr('bmatrix_write',nf90_def_var(grpid,'sigmab',nc_kind_real,(/nx_id,ny_id/),sigmab_id))
call ncerr('bmatrix_write',nf90_def_var(grpid,'dirac_cov',nc_kind_real,(/nx_id,ny_id/),dirac_cov_id))
call ncerr('bmatrix_write',nf90_def_var(grpid,'dirac_cor',nc_kind_real,(/nx_id,ny_id/),dirac_cor_id))

! Write variables
call ncerr('bmatrix_write',nf90_put_var(grpid,sigmab_id,sigmab_2d))
call ncerr('bmatrix_write',nf90_put_var(grpid,dirac_cov_id,Bx_2d))
call ncerr('bmatrix_write',nf90_put_var(grpid,dirac_cor_id,Cx_2d))

end subroutine bmatrix_write

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_spvar
! Purpose: apply spectral variances
!----------------------------------------------------------------------
subroutine bmatrix_apply_spvar(bmatrix,geom,v)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(inout) :: v(geom%nh)

! Local variables
integer :: k,l
real(kind_real) :: kstarsq,spvar
complex(8) :: cp(geom%kmax+1,geom%ny),cp_full(0:geom%kmax,-geom%lmax:geom%lmax)

! Real to complex
call geom%real_to_complex(v,cp)

! Complex to full
call geom%complex_to_full(cp,cp_full)

! Apply spectral variance
do k=0,geom%kmax
   do l=-geom%lmax,geom%lmax
      cp_full(k,l) = cp_full(k,l)*bmatrix%spvar(k,l)
   end do
end do

! Full to complex
call geom%full_to_complex(cp_full,cp)

! Complex to real
call geom%complex_to_real(cp,v)

end subroutine bmatrix_apply_spvar

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_spvar_sqrt
! Purpose: apply spectral variances, square-root
!----------------------------------------------------------------------
subroutine bmatrix_apply_spvar_sqrt(bmatrix,geom,v)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(inout) :: v(geom%nh)

! Local variables
integer :: k,l
real(kind_real) :: kstarsq,spvar
complex(8) :: cp(geom%kmax+1,geom%ny),cp_full(0:geom%kmax,-geom%lmax:geom%lmax)

! Real to complex
call geom%real_to_complex(v,cp)

! Complex to full
call geom%complex_to_full(cp,cp_full)

! Apply spectral variance square-root
do k=0,geom%kmax
   do l=-geom%lmax,geom%lmax
      cp_full(k,l) = cp_full(k,l)*sqrt(bmatrix%spvar(k,l))
   end do
end do

! Full to complex
call geom%full_to_complex(cp_full,cp)

! Complex to real
call geom%complex_to_real(cp,v)

end subroutine bmatrix_apply_spvar_sqrt

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_spvar_sqrt_inv
! Purpose: apply spectral variances, square-root inverse
!----------------------------------------------------------------------
subroutine bmatrix_apply_spvar_sqrt_inv(bmatrix,geom,v)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(inout) :: v(geom%nh)

! Local variables
integer :: k,l
real(kind_real) :: kstarsq,spvar
complex(8) :: cp(geom%kmax+1,geom%ny),cp_full(0:geom%kmax,-geom%lmax:geom%lmax)

! Real to complex
call geom%real_to_complex(v,cp)

! Complex to full
call geom%complex_to_full(cp,cp_full)

! Apply spectral variance square-root
do k=0,geom%kmax
   do l=-geom%lmax,geom%lmax
      cp_full(k,l) = cp_full(k,l)/sqrt(bmatrix%spvar(k,l))
   end do
end do

! Full to complex
call geom%full_to_complex(cp_full,cp)

! Complex to real
call geom%complex_to_real(cp,v)

end subroutine bmatrix_apply_spvar_sqrt_inv

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt
! Purpose: apply square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt(bmatrix,geom,v,x)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: v(geom%nh)
real(kind_real),intent(out) :: x(geom%nh)

! Local variables
real(kind_real) :: vtmp(geom%nh)

! Adjoint FFT inverse
call geom%gp2sp(v,vtmp)

! Apply spectral standard-deviation
call bmatrix%apply_spvar_sqrt(geom,vtmp)

! Adjoint FFT
call geom%sp2gp(vtmp,x)

! Apply grid-point standard-deviation
x = x*bmatrix%sigmab

end subroutine bmatrix_apply_sqrt

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt_ad
! Purpose: apply adjoint of the square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt_ad(bmatrix,geom,x,v)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: v(geom%nh)

! Local variables
real(kind_real) :: xtmp(geom%nh),vtmp(geom%nh)

! Apply grid-point standard-deviation
xtmp = x*bmatrix%sigmab

! FFT
call geom%gp2sp(xtmp,vtmp)

! Apply spectral standard-deviation
call bmatrix%apply_spvar_sqrt(geom,vtmp)

! FFT inverse
call geom%sp2gp(vtmp,v)

end subroutine bmatrix_apply_sqrt_ad

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt_inv
! Purpose: apply inverse square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt_inv(bmatrix,geom,x,v)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: v(geom%nh)

! Local variables
real(kind_real) :: xtmp(geom%nh),vtmp(geom%nh)

! Apply grid-point standard-deviation inverse
xtmp = x/bmatrix%sigmab

! Adjoint FFT inverse
call geom%gp2sp(xtmp,vtmp)

! Apply spectral standard-deviation inverse
call bmatrix%apply_spvar_sqrt_inv(geom,vtmp)

! Adjoint FFT
call geom%sp2gp(vtmp,v)

end subroutine bmatrix_apply_sqrt_inv

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_sqrt_ad_inv
! Purpose: apply inverse adjoint of the square-root of the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply_sqrt_ad_inv(bmatrix,geom,v,x)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: v(geom%nh)
real(kind_real),intent(out) :: x(geom%nh)

! Local variables
real(kind_real) :: vtmp(geom%nh)

! FFT
call geom%gp2sp(v,vtmp)

! Apply spectral standard-deviation inverse
call bmatrix%apply_spvar_sqrt_inv(geom,vtmp)

! FFT inverse
call geom%sp2gp(vtmp,x)

! Apply grid-point standard-deviation inverse
x = x/bmatrix%sigmab

end subroutine bmatrix_apply_sqrt_ad_inv

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply
! Purpose: apply B matrix
!----------------------------------------------------------------------
subroutine bmatrix_apply(bmatrix,geom,x,bx)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: bx(geom%nh)

! Local variables
real(kind_real) :: v(geom%nh)

! Apply adjoint of the square-root of the B matrix
call bmatrix%apply_sqrt_ad(geom,x,v)

! Apply square-root of the B matrix
call bmatrix%apply_sqrt(geom,v,bx)

end subroutine bmatrix_apply

!----------------------------------------------------------------------
! Subroutine: bmatrix_apply_inv
! Purpose: apply B matrix inverse
!----------------------------------------------------------------------
subroutine bmatrix_apply_inv(bmatrix,geom,x,binvx)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(in) :: x(geom%nh)
real(kind_real),intent(out) :: binvx(geom%nh)

! Local variables
real(kind_real) :: v(geom%nh)

! Apply inverse of the square-root of the B matrix
call bmatrix%apply_sqrt_inv(geom,x,v)

! Apply inverse adjoint of the square-root of the B matrix
call bmatrix%apply_sqrt_ad_inv(geom,v,binvx)

end subroutine bmatrix_apply_inv

!----------------------------------------------------------------------
! Subroutine: bmatrix_randomize
! Purpose: randomize the B matrix
!----------------------------------------------------------------------
subroutine bmatrix_randomize(bmatrix,geom,x)

implicit none

! Passed variables
class(bmatrix_type),intent(in) :: bmatrix
type(geom_type),intent(in) :: geom
real(kind_real),intent(out) :: x(geom%nh)

! Local variables
real(kind_real) :: nu(geom%nh)

! Gaussian random vector
call rand_normal(geom%nh,nu)

! Apply B matrix square-root
call bmatrix%apply_sqrt(geom,nu,x)

end subroutine bmatrix_randomize

end module type_bmatrix
