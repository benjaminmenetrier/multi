!----------------------------------------------------------------------
! Module: type_algo
! Purpose: lanczos and planczosif algorithms
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_algo

use tools_const
use tools_kinds
use type_bmatrix
use type_geom
use type_hoperator
use type_lmp
use type_rmatrix

implicit none

type algo_type
   integer :: nimax
   real(kind_real),allocatable :: xg(:)
   real(kind_real),allocatable :: dvb(:)
   real(kind_real),allocatable :: dxbbar(:)
   real(kind_real),allocatable :: d(:)

   real(kind_real),allocatable :: dx(:,:)
   real(kind_real),allocatable :: dva(:)
   real(kind_real),allocatable :: dxabar(:)
   real(kind_real),allocatable :: jb_quad(:)
   real(kind_real),allocatable :: jo_quad(:)
   real(kind_real),allocatable :: j_quad(:)
   real(kind_real),allocatable :: jb_full(:)
   real(kind_real),allocatable :: jo_full(:)
   real(kind_real),allocatable :: j_full(:)
   real(kind_real),allocatable :: eigenval(:)
   real(kind_real),allocatable :: eigenvec(:,:)
   real(kind_real),allocatable :: lancvec(:,:)
   real(kind_real) :: lastbeta
   real(kind_real),allocatable :: rho_sqrt(:)
   real(kind_real),allocatable :: beta(:)
contains
   procedure :: alloc => algo_alloc
   procedure :: write => algo_write
   procedure :: apply_lanczos => algo_apply_lanczos
   procedure :: apply_planczosif => algo_apply_planczosif
   procedure :: compute_cost => algo_compute_cost
end type algo_type

contains

!----------------------------------------------------------------------
! Subroutine: algo_alloc
! Purpose: allocate algo
!----------------------------------------------------------------------
subroutine algo_alloc(algo,geom,ni,nobs)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni
integer,intent(in) :: nobs

! Release memory
if (allocated(algo%xg)) deallocate(algo%xg)
if (allocated(algo%dvb)) deallocate(algo%dvb)
if (allocated(algo%dxbbar)) deallocate(algo%dxbbar)
if (allocated(algo%d)) deallocate(algo%d)
if (allocated(algo%dx)) deallocate(algo%dx)
if (allocated(algo%dva)) deallocate(algo%dva)
if (allocated(algo%dxabar)) deallocate(algo%dxabar)
if (allocated(algo%jb_quad)) deallocate(algo%jb_quad)
if (allocated(algo%jo_quad)) deallocate(algo%jo_quad)
if (allocated(algo%j_quad)) deallocate(algo%j_quad)
if (allocated(algo%jb_full)) deallocate(algo%jb_full)
if (allocated(algo%jo_full)) deallocate(algo%jo_full)
if (allocated(algo%j_full)) deallocate(algo%j_full)
if (allocated(algo%eigenval)) deallocate(algo%eigenval)
if (allocated(algo%eigenvec)) deallocate(algo%eigenvec)
if (allocated(algo%lancvec)) deallocate(algo%lancvec)
if (allocated(algo%rho_sqrt)) deallocate(algo%rho_sqrt)
if (allocated(algo%beta)) deallocate(algo%beta)

! Allocation
allocate(algo%xg(geom%nh))
allocate(algo%dvb(geom%nh))
allocate(algo%dxbbar(geom%nh))
allocate(algo%d(nobs))
allocate(algo%dx(geom%nh,ni))
allocate(algo%dva(geom%nh))
allocate(algo%dxabar(geom%nh))
allocate(algo%jb_quad(0:ni))
allocate(algo%jo_quad(0:ni))
allocate(algo%j_quad(0:ni))
allocate(algo%jb_full(0:ni))
allocate(algo%jo_full(0:ni))
allocate(algo%j_full(0:ni))
allocate(algo%eigenval(ni))
allocate(algo%eigenvec(ni,ni))
allocate(algo%lancvec(geom%nh,ni+1))
allocate(algo%rho_sqrt(0:ni))
allocate(algo%beta(0:ni))

! Initialization at missing value
algo%nimax = msv_int
algo%xg = msv_real
algo%dvb = msv_real
algo%dxbbar = msv_real
algo%d = msv_real
algo%dx = msv_real
algo%dva = msv_real
algo%dxabar = msv_real
algo%jb_quad = msv_real
algo%jo_quad = msv_real
algo%j_quad = msv_real
algo%jb_full = msv_real
algo%jo_full = msv_real
algo%j_full = msv_real
algo%eigenval = msv_real
algo%eigenvec = msv_real
algo%lancvec = msv_real
algo%lastbeta = msv_real
algo%rho_sqrt = msv_real
algo%beta = msv_real

end subroutine algo_alloc

!----------------------------------------------------------------------
! Subroutine: algo_write
! Purpose: write algo
!----------------------------------------------------------------------
subroutine algo_write(algo,geom,grpid,subgrpid)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
integer,intent(in) :: grpid
integer,intent(in) :: subgrpid

! Local variables
integer :: ii
integer :: nx_id,ny_id,nimax_id,nimax_plus_one_id
integer :: dx_id,jb_quad_id,jo_quad_id,j_quad_id,jb_full_id,jo_full_id,j_full_id,eigenval_id,rho_sqrt_id,beta_id
real(kind_real) :: dx_2d(geom%nx,geom%ny)

! Create dimensions
call ncerr('algo_write',nf90_inq_dimid(grpid,'nx',nx_id))
call ncerr('algo_write',nf90_inq_dimid(grpid,'ny',ny_id))
call ncerr('algo_write',nf90_def_dim(subgrpid,'nimax',algo%nimax,nimax_id))
call ncerr('algo_write',nf90_def_dim(subgrpid,'nimax_plus_one',algo%nimax+1,nimax_plus_one_id))

! Create variables
call ncerr('algo_write',nf90_def_var(subgrpid,'dx',nc_kind_real,(/nx_id,ny_id,nimax_id/),dx_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jb_quad',nc_kind_real,(/nimax_plus_one_id/),jb_quad_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jo_quad',nc_kind_real,(/nimax_plus_one_id/),jo_quad_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'j_quad',nc_kind_real,(/nimax_plus_one_id/),j_quad_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jb_full',nc_kind_real,(/nimax_plus_one_id/),jb_full_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jo_full',nc_kind_real,(/nimax_plus_one_id/),jo_full_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'j_full',nc_kind_real,(/nimax_plus_one_id/),j_full_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'eigenval',nc_kind_real,(/nimax_id/),eigenval_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'rho_sqrt',nc_kind_real,(/nimax_plus_one_id/),rho_sqrt_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'beta',nc_kind_real,(/nimax_plus_one_id/),beta_id))

! Write variables
do ii=1,algo%nimax
   dx_2d = reshape(algo%dx(:,ii),(/geom%nx,geom%ny/))
   call ncerr('algo_write',nf90_put_var(subgrpid,dx_id,dx_2d,(/1,1,ii/),(/geom%nx,geom%ny,1/)))
end do
call ncerr('algo_write',nf90_put_var(subgrpid,jb_quad_id,algo%jb_quad(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jo_quad_id,algo%jo_quad(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,j_quad_id,algo%j_quad(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jb_full_id,algo%jb_full(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jo_full_id,algo%jo_full(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,j_full_id,algo%j_full(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,eigenval_id,algo%eigenval(1:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,rho_sqrt_id,algo%rho_sqrt(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,beta_id,algo%beta(0:algo%nimax)))

end subroutine algo_write

!----------------------------------------------------------------------
! Subroutine: algo_apply_lanczos
! Purpose: Lanczos algorithm in control space
!----------------------------------------------------------------------
subroutine algo_apply_lanczos(algo,geom,bmatrix,hoperator,rmatrix,ni,lmp,shutoff_type,shutoff_value)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hoperator_type),intent(in) :: hoperator
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: shutoff_type
real(kind_real),intent(in) :: shutoff_value

! Local variables
integer :: ii,ji,info
real(kind_real) :: accuracy
real(kind_real) :: xtmp(geom%nh),ytmp(hoperator%nobs),ytmp2(hoperator%nobs),y(ni,ni),subdiag(ni-1),work(max(1,2*ni-2))
real(kind_real) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(kind_real) :: vtmp(geom%nh),vtmp2(geom%nh)
real(kind_real) :: p(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni),rhs(ni),u(geom%nh,ni),v(geom%nh,0:ni+1)
real(kind_real) :: w(geom%nh,0:ni)
real(kind_real),allocatable :: eigenvec(:,:)
logical :: convergence

! Initialization
s(:,0) = zero
v(:,0) = zero
call rmatrix%apply_inv(algo%d,ytmp)
call hoperator%apply_ad(geom,algo%xg,ytmp,xtmp)
call bmatrix%apply_sqrt_ad(geom,xtmp,vtmp)
vtmp = algo%dvb+vtmp
call lmp%apply_sqrt_ad(geom,ni,lmp%io,vtmp,r(:,0))
beta(0) = sqrt(sum(r(:,0)**2))
v(:,1) = r(:,0)/beta(0)
beta(1) = zero
rhs = zero
rhs(1) = beta(0)
convergence = .false.
ii = 0

! Initialize quadratic cost function
algo%jb_quad(0) = 0.5*sum((zero-algo%dvb)**2)
call rmatrix%apply_inv(algo%d,ytmp)
algo%jo_quad(0) = 0.5*sum(algo%d*ytmp)
algo%j_quad(0) = algo%jb_quad(0)+algo%jo_quad(0)

do while ((.not.convergence).and.(ii<ni))
   ! Update
   ii = ii+1
   call lmp%apply_sqrt(geom,ni,lmp%io,v(:,ii),vtmp)
   call bmatrix%apply_sqrt(geom,vtmp,xtmp)
   call hoperator%apply_tl(geom,algo%xg,xtmp,ytmp)
   call rmatrix%apply_inv(ytmp,ytmp2)
   call hoperator%apply_ad(geom,algo%xg,ytmp2,xtmp)
   call bmatrix%apply_sqrt_ad(geom,xtmp,vtmp2)
   vtmp = vtmp+vtmp2
   call lmp%apply_sqrt_ad(geom,ni,lmp%io,vtmp,q(:,ii))
   q(:,ii) = q(:,ii)-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*v(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   beta(ii+1) = sqrt(sum(w(:,ii)**2))
   v(:,ii+1) = w(:,ii)/beta(ii+1)

   ! Stopping criterion based on beta
   if (shutoff_type==2) then
      if (beta(ii+1)<shutoff_value) then 
         write(*,'(a)') '         Convergence reached, based on beta'
         convergence = .true.
         algo%nimax = ii-1
      end if
   end if

   ! Compute eigenpairs
   if (ii==1) then
      algo%eigenval(1) = alpha(1)
      algo%eigenvec(1,1) = one
   else
      allocate(eigenvec(ii,ii))
      algo%eigenval(1:ii) = alpha(1:ii)
      subdiag(1:ii-1) = beta(2:ii)
      call dsteqr('I',ii,algo%eigenval(1:ii),subdiag(1:ii-1),eigenvec,ii,work(1:2*ii-2),info)
      algo%eigenvec(1:ii,1:ii) = eigenvec
      if (info/=0) then
         write(*,'(a,i2)') '         Error in dsteqr: ',info
         stop
      end if
      algo%eigenval(1:ii) = max(algo%eigenval(1:ii),one)
      deallocate(eigenvec)
   end if

   ! Update increment
   y(1:ii,ii) = matmul(algo%eigenvec(1:ii,1:ii),matmul(transpose(algo%eigenvec(1:ii,1:ii)),rhs(1:ii))/algo%eigenval(1:ii))
   s(:,ii) = matmul(v(:,1:ii),y(1:ii,ii))
   call lmp%apply_sqrt(geom,ni,lmp%io,s(:,ii),u(:,ii))
   call bmatrix%apply_sqrt(geom,u(:,ii),algo%dx(:,ii))

   ! Stopping criterion based on the Ritz pairs approximation
   if (shutoff_type==3) then
      do ji=1,ii
         accuracy = abs(beta(ii+1)*s(ji,ii))/algo%eigenval(ji)
         if (accuracy>shutoff_value) then
            write(*,'(a)') '         Convergence reached, based on the Ritz pairs approximation'
            convergence = .true.
            algo%nimax = ii-1
         end if
      end do
   end if

   ! Compute quadratic cost function
   algo%jb_quad(ii) = 0.5*sum((u(:,ii)-algo%dvb)**2)
   call hoperator%apply_tl(geom,algo%xg,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv((algo%d-ytmp),ytmp2)
   algo%jo_quad(ii) = 0.5*sum((algo%d-ytmp)*ytmp2)
   algo%j_quad(ii) = algo%jb_quad(ii)+algo%jo_quad(ii)

   ! Stopping criterion based on the Jb increase
   if (shutoff_type==1) then
      if (algo%jb_quad(ii)/algo%jb_quad(ii-1)<shutoff_value) then
         write(*,'(a)') '         Convergence reached, based on the Jb increase'
         convergence = .true.
         algo%nimax = ii-1
       end if
   end if

   ! Stopping criterion based on the J stagnation
   if (abs(algo%j_quad(ii)-algo%j_quad(ii-1))<1.0e-12_kind_real*algo%j_quad(ii-1)) then
      write(*,'(a)') '         Convergence reached, based on the J stagnation'
      convergence = .true.
      algo%nimax = ii-1
   end if

   ! Check cost function
   if (algo%j_quad(ii)>algo%j_quad(ii-1)) then  
      write(*,'(a)') '         Error: increasing cost function in Lanczos'
      stop
   end if
end do

if (.not.convergence) then
   ! No convergence
   write(*,'(a)') '         Convergence not reached'
   algo%nimax = ni
end if

! Final update
if (algo%nimax>0) algo%dva = u(:,algo%nimax)

! Lanczos vectors
if (algo%nimax>0) algo%lancvec = v(:,1:algo%nimax+1)

! Last beta value
if (algo%nimax>0) algo%lastbeta = beta(algo%nimax+1)

! Lanczos to CG
beta(0) = sqrt(sum(r(:,0)**2))
p(:,0) = r(:,0)
rho(0) = sum(r(:,0)**2)
algo%rho_sqrt(0) = sqrt(rho(0))
algo%beta(0) = beta(0)
do ii=1,algo%nimax
   r(:,ii) = -beta(ii+1)*y(ii,ii)*v(:,ii+1)
   rho(ii) = sum(r(:,ii)**2)
   algo%rho_sqrt(ii)=sqrt(rho(ii))
   beta(ii) = rho(ii)/rho(ii-1)
   algo%beta(ii) = beta(ii)
   p(:,ii) = r(:,ii)+beta(ii)*p(:,ii-1)
   alpha(ii-1) = sqrt(sum((s(:,ii)-s(:,ii-1))**2)/sum(p(:,ii-1)**2))
   q(:,ii-1) = -(r(:,ii)-r(:,ii-1))/alpha(ii-1)
end do

end subroutine algo_apply_lanczos

!----------------------------------------------------------------------
! Subroutine: algo_apply_planczosif
! Purpose: PLanczosIF algorithm in linear space
!----------------------------------------------------------------------
subroutine algo_apply_planczosif(algo,geom,bmatrix,hoperator,rmatrix,ni,lmp,shutoff_type,shutoff_value)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hoperator_type),intent(in) :: hoperator
type(rmatrix_type),intent(in) :: rmatrix
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: shutoff_type
real(kind_real),intent(in) :: shutoff_value

! Local variables
integer :: ii,ji,info
real(kind_real) :: accuracy
real(kind_real) :: xtmp(geom%nh),ytmp(hoperator%nobs),ytmp2(hoperator%nobs),y(ni,ni),rhs(ni),subdiag(ni-1),work(max(1,2*ni-2))
real(kind_real) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(kind_real) :: dxb(geom%nh),dxbar(geom%nh,ni)
real(kind_real) :: l(geom%nh,0:ni),p(geom%nh,0:ni),pbar(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni)
real(kind_real) :: tbar(geom%nh,0:ni),t(geom%nh,0:ni),v(geom%nh,0:ni+1),w(geom%nh,0:ni),zbar(geom%nh,0:ni+1),z(geom%nh,0:ni+1)
real(kind_real),allocatable :: eigenvec(:,:)
logical :: convergence

! Initialization
s(:,0) = zero
v(:,0) = zero
call rmatrix%apply_inv(algo%d,ytmp)
call hoperator%apply_ad(geom,algo%xg,ytmp,xtmp)
r(:,0) = algo%dxbbar+xtmp
call lmp%apply(geom,ni,lmp%io,r(:,0),tbar(:,0))
call bmatrix%apply(geom,tbar(:,0),t(:,0))
beta(0) = sqrt(sum(r(:,0)*t(:,0)))
v(:,1) = r(:,0)/beta(0)
zbar(:,1) = tbar(:,0)/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = zero
rhs = zero
rhs(1) = beta(0)
convergence = .false.
ii = 0

! Initialize cost function
call bmatrix%apply(geom,algo%dxbbar,dxb)
algo%jb_quad(0) = 0.5*sum((zero-dxb)*(zero-algo%dxbbar))
call rmatrix%apply_inv(algo%d,ytmp)
algo%jo_quad(0) = 0.5*sum(algo%d*ytmp)
algo%j_quad(0) = algo%jb_quad(0)+algo%jo_quad(0)

do while ((.not.convergence).and.(ii<ni))
   ! Update
   ii = ii+1
   call hoperator%apply_tl(geom,algo%xg,z(:,ii),ytmp)
   call rmatrix%apply_inv(ytmp,ytmp2)
   call hoperator%apply_ad(geom,algo%xg,ytmp2,xtmp)
   q(:,ii) = zbar(:,ii)+xtmp-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call lmp%apply(geom,ni,lmp%io,w(:,ii),tbar(:,ii))
   call bmatrix%apply(geom,tbar(:,ii),t(:,ii))
   beta(ii+1) = sqrt(sum(w(:,ii)*t(:,ii)))
   v(:,ii+1) = w(:,ii)/beta(ii+1)
   zbar(:,ii+1) = tbar(:,ii)/beta(ii+1)
   z(:,ii+1) = t(:,ii)/beta(ii+1)

   ! Stopping criterion based on beta
   if (shutoff_type==2) then
      if (beta(ii+1)<shutoff_value) then 
         write(*,'(a)') '         Convergence reached, based on beta'
         convergence = .true.
         algo%nimax = ii-1
      end if
   end if

   ! Compute eigenpairs
   if (ii==1) then
      algo%eigenval(1) = alpha(1)
      algo%eigenvec(1,1) = one
   else
      allocate(eigenvec(ii,ii))
      algo%eigenval(1:ii) = alpha(1:ii)
      subdiag(1:ii-1) = beta(2:ii)
      call dsteqr('I',ii,algo%eigenval(1:ii),subdiag(1:ii-1),eigenvec,ii,work(1:2*ii-2),info)
      algo%eigenvec(1:ii,1:ii) = eigenvec
      if (info/=0) then
         write(*,'(a,i2)') '         Error in dsteqr: ',info
         stop
      end if
      algo%eigenval(1:ii) = max(algo%eigenval(1:ii),one)
      deallocate(eigenvec)
   end if

   ! Update increment
   y(1:ii,ii) = matmul(algo%eigenvec(1:ii,1:ii),matmul(transpose(algo%eigenvec(1:ii,1:ii)),rhs(1:ii))/algo%eigenval(1:ii))
   s(:,ii) = matmul(z(:,1:ii),y(1:ii,ii))
   dxbar(:,ii) = matmul(zbar(:,1:ii),y(1:ii,ii))
   algo%dx(:,ii) = s(:,ii)

   ! Stopping criterion based on the Ritz pairs approximation
   if (shutoff_type==3) then
      do ji=1,ii
         accuracy = abs(beta(ii+1)*s(ji,ii))/algo%eigenval(ji)
         if (accuracy>shutoff_value) then
            write(*,'(a)') '         Convergence reached, based on the Ritz pairs approximation'
            convergence = .true.
            algo%nimax = ii-1
         end if
      end do
   end if

   ! Compute quadratic cost function
   algo%jb_quad(ii) = 0.5*sum((algo%dx(:,ii)-dxb)*(dxbar(:,ii)-algo%dxbbar))
   call hoperator%apply_tl(geom,algo%xg,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv(algo%d-ytmp,ytmp2)
   algo%jo_quad(ii) = 0.5*sum((algo%d-ytmp)*ytmp2)
   algo%j_quad(ii) = algo%jb_quad(ii)+algo%jo_quad(ii)

   ! Stopping criterion based on the Jb increase
   if (shutoff_type==1) then
      if (algo%jb_quad(ii)/algo%jb_quad(ii-1)<shutoff_value) then
         write(*,'(a)') '         Convergence reached, based on the Jb increase'
         convergence = .true.
         algo%nimax = ii-1
       end if
   end if

   ! Stopping criterion based on the J stagnation
   if (abs(algo%j_quad(ii)-algo%j_quad(ii-1))<1.0e-12_kind_real*algo%j_quad(ii-1)) then
      write(*,'(a)') '         Convergence reached, based on the J stagnation'
      convergence = .true.
      algo%nimax = ii-1
   end if

   ! Check cost function
   if (algo%j_quad(ii)>algo%j_quad(ii-1)) then  
      write(*,'(a)') '         Error: increasing cost function in Lanczos'
      stop
   end if
end do

if (.not.convergence) then
   ! No convergence
   write(*,'(a)') '         Convergence not reached'
   algo%nimax = ni
end if

! Final update
if (algo%nimax>0) algo%dxabar = dxbar(:,algo%nimax)

! Lanczos vectors
if (algo%nimax>0) algo%lancvec = v(:,1:algo%nimax+1)

! Last beta value
if (algo%nimax>0) algo%lastbeta = beta(algo%nimax+1)

! PLanczosIF to PCGIF
call bmatrix%apply(geom,r(:,0),l(:,0))
call lmp%apply(geom,algo%nimax,lmp%io,r(:,0),zbar(:,0))
call lmp%apply_ad(geom,algo%nimax,lmp%io,l(:,0),z(:,0))
rho(0) = sum(r(:,0)*z(:,0))
beta(0) = sqrt(sum(r(:,0)*z(:,0)))
pbar(:,0) = zbar(:,0)
p(:,0) = z(:,0)
algo%rho_sqrt(0)=sqrt(rho(0))
algo%beta(0)=beta(0)
do ii=1,algo%nimax
   r(:,ii) = -beta(ii+1)*y(ii,ii)*v(:,ii+1)
   call bmatrix%apply(geom,r(:,ii),l(:,ii))
   call lmp%apply(geom,algo%nimax,lmp%io,r(:,ii),zbar(:,ii))
   call lmp%apply_ad(geom,algo%nimax,lmp%io,l(:,ii),z(:,ii))
   rho(ii) = sum(r(:,ii)*z(:,ii))
   algo%rho_sqrt(ii)=sqrt(rho(ii))
   beta(ii) = rho(ii)/rho(ii-1)
   algo%beta(ii) = beta(ii)
   pbar(:,ii) = zbar(:,ii)+beta(ii)*pbar(:,ii-1)
   p(:,ii) = z(:,ii)+beta(ii)*p(:,ii-1)
   alpha(ii-1) = sqrt(sum((s(:,ii)-s(:,ii-1))**2)/sum(p(:,ii-1)**2))
   q(:,ii-1) = -(r(:,ii)-r(:,ii-1))/alpha(ii-1)
   t(:,ii-1) = -(l(:,ii)-l(:,ii-1))/alpha(ii-1)
end do

end subroutine algo_apply_planczosif

!----------------------------------------------------------------------
! Subroutine: algo_compute_cost
! Purpose: compute full cost function
!----------------------------------------------------------------------
subroutine algo_compute_cost(algo,geom,bmatrix,hoperator,rmatrix,xb)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hoperator_type),intent(in) :: hoperator
type(rmatrix_type),intent(in) :: rmatrix
real(kind_real),intent(in) :: xb(geom%nh)

! Local variables
integer :: ii
real(kind_real) :: xtmp(geom%nh),ytmp(hoperator%nobs),ytmp2(hoperator%nobs)

! Initialize full cost function
call bmatrix%apply_inv(geom,algo%xg-xb,xtmp)
algo%jb_full(0) = 0.5*sum((algo%xg-xb)*xtmp)
call hoperator%apply(geom,algo%xg,ytmp)
ytmp = ytmp-hoperator%yo
call rmatrix%apply_inv(ytmp,ytmp2)
algo%jo_full(0) = 0.5*sum(ytmp*ytmp2)
algo%j_full(0) = algo%jb_full(0)+algo%jo_full(0)

do ii=1,algo%nimax
   ! Compute full cost function
   call bmatrix%apply_inv(geom,algo%xg+algo%dx(:,ii)-xb,xtmp)
   algo%jb_full(ii) = 0.5*sum((algo%xg+algo%dx(:,ii)-xb)*xtmp)
   call hoperator%apply(geom,(algo%xg+algo%dx(:,ii)),ytmp)
   ytmp = ytmp-hoperator%yo
   call rmatrix%apply_inv(ytmp,ytmp2)
   algo%jo_full(ii) = 0.5*sum(ytmp*ytmp2)
   algo%j_full(ii) = algo%jb_full(ii)+algo%jo_full(ii)
end do

end subroutine algo_compute_cost

end module type_algo
