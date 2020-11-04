!----------------------------------------------------------------------
! Module: type_algo
! Purpose: lanczos and planczosif algorithms
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_algo

use type_bmatrix
use type_geom
use type_hmatrix
use type_lmp
use type_rmatrix

implicit none

type algo_type
   integer :: nimax
   real(8),allocatable :: dx(:,:)
   real(8),allocatable :: dva(:)
   real(8),allocatable :: dxabar(:)
   real(8),allocatable :: jb(:)
   real(8),allocatable :: jo(:)
   real(8),allocatable :: j(:)
   real(8),allocatable :: jb_nl(:)
   real(8),allocatable :: jo_nl(:)
   real(8),allocatable :: j_nl(:)
   real(8),allocatable :: eigenval(:)
   real(8),allocatable :: eigenvec(:,:)
   real(8),allocatable :: lancvec(:,:)
   real(8) :: lastbeta
   real(8),allocatable :: rho_sqrt(:)
   real(8),allocatable :: beta(:)
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
subroutine algo_alloc(algo,geom,ni)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
integer,intent(in) :: ni

! Release memory
if (allocated(algo%dx)) deallocate(algo%dx)
if (allocated(algo%dva)) deallocate(algo%dva)
if (allocated(algo%dxabar)) deallocate(algo%dxabar)
if (allocated(algo%jb)) deallocate(algo%jb)
if (allocated(algo%jo)) deallocate(algo%jo)
if (allocated(algo%j)) deallocate(algo%j)
if (allocated(algo%jb_nl)) deallocate(algo%jb_nl)
if (allocated(algo%jo_nl)) deallocate(algo%jo_nl)
if (allocated(algo%j_nl)) deallocate(algo%j_nl)
if (allocated(algo%eigenval)) deallocate(algo%eigenval)
if (allocated(algo%eigenvec)) deallocate(algo%eigenvec)
if (allocated(algo%lancvec)) deallocate(algo%lancvec)
if (allocated(algo%rho_sqrt)) deallocate(algo%rho_sqrt)
if (allocated(algo%beta)) deallocate(algo%beta)

! Allocation
allocate(algo%dx(geom%nh,ni))
allocate(algo%dva(geom%nh))
allocate(algo%dxabar(geom%nh))
allocate(algo%jb(0:ni))
allocate(algo%jo(0:ni))
allocate(algo%j(0:ni))
allocate(algo%jb_nl(0:ni))
allocate(algo%jo_nl(0:ni))
allocate(algo%j_nl(0:ni))
allocate(algo%eigenval(ni))
allocate(algo%eigenvec(ni,ni))
allocate(algo%lancvec(geom%nh,ni+1))
allocate(algo%rho_sqrt(0:ni))
allocate(algo%beta(0:ni))

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
integer :: dx_id,jb_id,jo_id,j_id,jb_nl_id,jo_nl_id,j_nl_id,eigenval_id,rho_sqrt_id,beta_id
real(8) :: dx_2d(geom%nx,geom%ny)

! Create dimensions
call ncerr('algo_write',nf90_inq_dimid(grpid,'nx',nx_id))
call ncerr('algo_write',nf90_inq_dimid(grpid,'ny',ny_id))
call ncerr('algo_write',nf90_def_dim(subgrpid,'nimax',algo%nimax,nimax_id))
call ncerr('algo_write',nf90_def_dim(subgrpid,'nimax_plus_one',algo%nimax+1,nimax_plus_one_id))

! Create variables
call ncerr('algo_write',nf90_def_var(subgrpid,'dx',nf90_double,(/nx_id,ny_id,nimax_id/),dx_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jb',nf90_double,(/nimax_plus_one_id/),jb_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jo',nf90_double,(/nimax_plus_one_id/),jo_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'j',nf90_double,(/nimax_plus_one_id/),j_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jb_nl',nf90_double,(/nimax_plus_one_id/),jb_nl_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'jo_nl',nf90_double,(/nimax_plus_one_id/),jo_nl_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'j_nl',nf90_double,(/nimax_plus_one_id/),j_nl_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'eigenval',nf90_double,(/nimax_id/),eigenval_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'rho_sqrt',nf90_double,(/nimax_plus_one_id/),rho_sqrt_id))
call ncerr('algo_write',nf90_def_var(subgrpid,'beta',nf90_double,(/nimax_plus_one_id/),beta_id))

! Write variables
do ii=1,algo%nimax
   dx_2d = reshape(algo%dx(:,ii),(/geom%nx,geom%ny/))
   call ncerr('algo_write',nf90_put_var(subgrpid,dx_id,dx_2d))
end do
call ncerr('algo_write',nf90_put_var(subgrpid,jb_id,algo%jb(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jo_id,algo%jo(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,j_id,algo%j(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jb_nl_id,algo%jb_nl(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,jo_nl_id,algo%jo_nl(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,j_nl_id,algo%j_nl(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,eigenval_id,algo%eigenval(1:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,rho_sqrt_id,algo%rho_sqrt(0:algo%nimax)))
call ncerr('algo_write',nf90_put_var(subgrpid,beta_id,algo%beta(0:algo%nimax)))

end subroutine algo_write

!----------------------------------------------------------------------
! Subroutine: algo_apply_lanczos
! Purpose: Lanczos algorithm in control space
!----------------------------------------------------------------------
subroutine algo_apply_lanczos(algo,geom,bmatrix,hmatrix,rmatrix,dvbbar,d,ni,lmp,shutoff_type,shutoff_value)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: dvbbar(geom%nh)
real(8),intent(in) :: d(hmatrix%nobs)
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
integer,intent(in)  :: shutoff_type
real(8),intent(in)  :: shutoff_value

! Local variables
integer :: ii,ji,info
real(8) :: accuracy
real(8) :: xtmp(geom%nh),ytmp(hmatrix%nobs),ytmp2(hmatrix%nobs),y(ni,ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: vtmp(geom%nh),vtmp2(geom%nh)
real(8) :: p(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni),rhs(ni),u(geom%nh,ni),v(geom%nh,0:ni+1),w(geom%nh,0:ni)
real(8),allocatable :: eigenvec(:,:)
logical :: convergence

! Initialization
s(:,0) = 0.0
v(:,0) = 0.0
call rmatrix%apply_inv(d,ytmp)
call hmatrix%apply_ad(geom,ytmp,xtmp)
call bmatrix%apply_sqrt_ad(geom,xtmp,vtmp)
vtmp = dvbbar+vtmp
call lmp%apply_sqrt_ad(geom,ni,lmp%io,vtmp,r(:,0))
beta(0) = sqrt(sum(r(:,0)**2))
v(:,1) = r(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)
convergence = .false.
ii = 0

! Initialize linear cost function
algo%jb(0) = 0.5*sum((0.0-dvbbar)**2)
call rmatrix%apply_inv(d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)
algo%j(0) = algo%jb(0)+algo%jo(0)

do while ((.not.convergence).and.(ii<ni))
   ! Update
   ii = ii+1
   call lmp%apply_sqrt(geom,ni,lmp%io,v(:,ii),vtmp)
   call bmatrix%apply_sqrt(geom,vtmp,xtmp)
   call hmatrix%apply(geom,xtmp,ytmp)
   call rmatrix%apply_inv(ytmp,ytmp2)
   call hmatrix%apply_ad(geom,ytmp2,xtmp)
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
      algo%eigenvec(1,1) = 1.0
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
      algo%eigenval(1:ii) = max(algo%eigenval(1:ii),1.0)
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

   ! Compute linear cost function
   algo%jb(ii) = 0.5*sum((u(:,ii)-dvbbar)**2)
   call hmatrix%apply(geom,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv((d-ytmp),ytmp2)
   algo%jo(ii) = 0.5*sum((d-ytmp)*ytmp2)
   algo%j(ii) = algo%jb(ii)+algo%jo(ii)

   ! Stopping criterion based on the Jb increase
   if (shutoff_type==1) then
      if (algo%jb(ii)/algo%jb(ii-1)<shutoff_value) then
         write(*,'(a)') '         Convergence reached, based on the Jb increase'
         convergence = .true.
         algo%nimax = ii-1
       end if
   end if

   ! Check cost function
   if (algo%j(ii)>algo%j(ii-1)) then
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
algo%dva = u(:,algo%nimax)

! Lanczos vectors
algo%lancvec = v(:,1:algo%nimax+1)

! Last beta value
algo%lastbeta = beta(algo%nimax+1)

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
subroutine algo_apply_planczosif(algo,geom,bmatrix,hmatrix,rmatrix,dxbbar,d,ni,lmp,shutoff_type,shutoff_value)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: dxbbar(geom%nh)
real(8),intent(in) :: d(hmatrix%nobs)
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
integer,intent(in) :: shutoff_type
real(8),intent(in) :: shutoff_value

! Local variables
integer :: ii,ji,info
real(8) :: accuracy
real(8) :: xtmp(geom%nh),ytmp(hmatrix%nobs),ytmp2(hmatrix%nobs),y(ni,ni),rhs(ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: dxb(geom%nh),dxbar(geom%nh,ni)
real(8) :: l(geom%nh,0:ni),p(geom%nh,0:ni),pbar(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni)
real(8) :: tbar(geom%nh,0:ni),t(geom%nh,0:ni),v(geom%nh,0:ni+1),w(geom%nh,0:ni),zbar(geom%nh,0:ni+1),z(geom%nh,0:ni+1)
real(8),allocatable :: eigenvec(:,:)
logical :: convergence

! Initialization
s(:,0) = 0.0
v(:,0) = 0.0
call rmatrix%apply_inv(d,ytmp)
call hmatrix%apply_ad(geom,ytmp,xtmp)
r(:,0) = dxbbar+xtmp
call lmp%apply(geom,ni,lmp%io,r(:,0),tbar(:,0))
call bmatrix%apply(geom,tbar(:,0),t(:,0))
beta(0) = sqrt(sum(r(:,0)*t(:,0)))
v(:,1) = r(:,0)/beta(0)
zbar(:,1) = tbar(:,0)/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)
convergence = .false.
ii = 0

! Initialize cost function
call bmatrix%apply(geom,dxbbar,dxb)
algo%jb(0) = 0.5*sum((0.0-dxb)*(0.0-dxbbar))
call rmatrix%apply_inv(d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)
algo%j(0) = algo%jb(0)+algo%jo(0)

do while ((.not.convergence).and.(ii<ni))
   ! Update
   ii = ii+1
   call hmatrix%apply(geom,z(:,ii),ytmp)
   call rmatrix%apply_inv(ytmp,ytmp2)
   call hmatrix%apply_ad(geom,ytmp2,xtmp)
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
      algo%eigenvec(1,1) = 1.0
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
      algo%eigenval(1:ii) = max(algo%eigenval(1:ii),1.0)
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

   ! Compute cost function
   algo%jb(ii) = 0.5*sum((algo%dx(:,ii)-dxb)*(dxbar(:,ii)-dxbbar))
   call hmatrix%apply(geom,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv(d-ytmp,ytmp2)
   algo%jo(ii) = 0.5*sum((d-ytmp)*ytmp2)
   algo%j(ii) = algo%jb(ii)+algo%jo(ii)

   ! Stopping criterion based on the Jb increase
   if (shutoff_type==1) then
      if (algo%jb(ii)/algo%jb(ii-1)<shutoff_value) then
         write(*,'(a)') '         Convergence reached, based on the Jb increase'
         convergence = .true.
         algo%nimax = ii-1
       end if
   end if

   ! Check cost function
   if (algo%j(ii)>algo%j(ii-1)) then
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
algo%dxabar = dxbar(:,algo%nimax)

! Lanczos vectors
algo%lancvec = v(:,1:algo%nimax+1)

! Last beta value
algo%lastbeta = beta(algo%nimax+1)

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
! Purpose: compute nonlinear cost function
!----------------------------------------------------------------------
subroutine algo_compute_cost(algo,geom,bmatrix,hmatrix,rmatrix,xb,xg,ni)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: xb(geom%nh)
real(8),intent(in) :: xg(geom%nh)
integer,intent(in) :: ni

! Local variables
integer :: ii
real(8) :: xtmp(geom%nh),ytmp(hmatrix%nobs),ytmp2(hmatrix%nobs)

! Initialize nonlinear cost function
call bmatrix%apply_inv(geom,xg-xb,xtmp)
algo%jb_nl(0) = 0.5*sum((xg-xb)*xtmp)
call hmatrix%apply(geom,xg,ytmp)
ytmp = ytmp-hmatrix%yo
call rmatrix%apply_inv(ytmp,ytmp2)
algo%jo_nl(0) = 0.5*sum(ytmp*ytmp2)
algo%j_nl(0) = algo%jb_nl(0)+algo%jo_nl(0)

do ii=1,algo%nimax
   ! Compute nonlinear cost function
   call bmatrix%apply_inv(geom,xg+algo%dx(:,ii)-xb,xtmp)
   algo%jb_nl(ii) = 0.5*sum((xg+algo%dx(:,ii)-xb)*xtmp)
   call hmatrix%apply(geom,(xg+algo%dx(:,ii)),ytmp)
   ytmp = ytmp-hmatrix%yo
   call rmatrix%apply_inv(ytmp,ytmp2)
   algo%jo_nl(ii) = 0.5*sum(ytmp*ytmp2)
   algo%j_nl(ii) = algo%jb_nl(ii)+algo%jo_nl(ii)
end do

end subroutine algo_compute_cost

end module type_algo
