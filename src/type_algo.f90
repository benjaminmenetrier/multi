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
use type_obsloc
use type_rmatrix

implicit none

type algo_type
   real(8),allocatable :: dx(:,:)
   real(8),allocatable :: dva(:)
   real(8),allocatable :: jb(:)
   real(8),allocatable :: jo(:)
   real(8),allocatable :: jb_nl(:)
   real(8),allocatable :: jo_nl(:)
   real(8),allocatable :: eigenval(:)
   real(8),allocatable :: eigenvec(:,:)
   real(8),allocatable :: lancvec(:,:)
   real(8) :: lastbeta
   real(8),allocatable :: rho_sqrt(:)
   real(8),allocatable :: beta(:)
contains
   procedure :: alloc => algo_alloc
   procedure :: apply_lanczos => algo_apply_lanczos
   procedure :: apply_planczosif => algo_apply_planczosif
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
if (allocated(algo%jb)) deallocate(algo%jb)
if (allocated(algo%jo)) deallocate(algo%jo)
if (allocated(algo%jb_nl)) deallocate(algo%jb_nl)
if (allocated(algo%jo_nl)) deallocate(algo%jo_nl)
if (allocated(algo%eigenval)) deallocate(algo%eigenval)
if (allocated(algo%eigenvec)) deallocate(algo%eigenvec)
if (allocated(algo%lancvec)) deallocate(algo%lancvec)
if (allocated(algo%rho_sqrt)) deallocate(algo%rho_sqrt)
if (allocated(algo%beta)) deallocate(algo%beta)

! Allocation
allocate(algo%dx(geom%nh,ni))
allocate(algo%dva(geom%nh))
allocate(algo%jb(0:ni))
allocate(algo%jo(0:ni))
allocate(algo%jb_nl(0:ni))
allocate(algo%jo_nl(0:ni))
allocate(algo%eigenval(ni))
allocate(algo%eigenvec(ni,ni))
allocate(algo%lancvec(geom%nh,ni+1))
allocate(algo%rho_sqrt(0:ni))
allocate(algo%beta(0:ni))

end subroutine algo_alloc

!----------------------------------------------------------------------
! Subroutine: algo_apply_lanczos
! Purpose: Lanczos algorithm in control space
!----------------------------------------------------------------------
subroutine algo_apply_lanczos(algo,geom,bmatrix,hmatrix,rmatrix,xb,xg,dvb_bar,d,ni,lmp,shutoff_type,shutoff_value)

implicit none

! Passed variables
class(algo_type),intent(inout) :: algo
type(geom_type),intent(in) :: geom
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: xb(geom%nh)
real(8),intent(in) :: xg(geom%nh)
real(8),intent(in) :: dvb_bar(geom%nh)
real(8),intent(in) :: d(hmatrix%nobs)
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
integer,intent(in)  :: shutoff_type
real(8),intent(in)  :: shutoff_value

! Local variables
integer :: ii,iii,info
real(8) :: xtmp(geom%nh),ytmp(hmatrix%nobs),ytmp2(hmatrix%nobs),y(ni,ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: vtmp(geom%nh),vtmp2(geom%nh)
real(8) :: p(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni),rhs(ni),u(geom%nh,ni),v(geom%nh,0:ni+1),w(geom%nh,0:ni)
real(8) :: accuracy
real(8),allocatable :: eigenvec(:,:)

! Initialization
s(:,0) = 0.0
v(:,0) = 0.0
call rmatrix%apply_inv(d,ytmp)
call hmatrix%apply_ad(geom,ytmp,xtmp)
call bmatrix%apply_sqrt_ad(geom,xtmp,vtmp)
vtmp = dvb_bar+vtmp
call lmp%apply_sqrt_ad(geom,ni,lmp%io,vtmp,r(:,0))
beta(0) = sqrt(sum(r(:,0)**2))
v(:,1) = r(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize linear cost function
algo%jb(0) = 0.5*sum((0.0-dvb_bar)**2)
call rmatrix%apply_inv(d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)

! Initialize nonlinear cost function
call bmatrix%apply_inv(geom,xg-xb,xtmp)
algo%jb_nl(0) = 0.5*sum((xg-xb)*xtmp)
call hmatrix%apply(geom,xg,ytmp)
ytmp = ytmp-hmatrix%yo
call rmatrix%apply_inv(ytmp,ytmp2)
algo%jo_nl(0) = 0.5*sum(ytmp*ytmp2)

do ii=1,ni
   ! Update
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

   ! Stop criterion:
   if (shutoff_type==2) then
      if (beta(ii+1)<shutoff_value) then
         write(*,'(a)') 'Convergence reached'
         stop
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
         write(*,'(a,i2)') 'Error in dsteqr: ',info
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

   ! rhs check:

   ! Stop criterion on the Ritz pairs approximation:
   if (shutoff_type==3) then
      do iii=1,ii
         accuracy=beta(ii+1)*s(iii,ii)/algo%eigenval(iii)
         write(*,'(a,e15.8,a,e15.8,a,e15.8,a,e15.8)') 'acc check:',accuracy,' ',beta(ii+1),' ',s(iii,ii),' ',algo%eigenval(iii)
         if (accuracy > shutoff_value) then
            !write(*,'(a,e15.8)') 'Unacceptable Ritz pair, with accuracy: ',
            !stop
         end if
      end do
   end if

   ! Compute linear cost function
   algo%jb(ii) = 0.5*sum((u(:,ii)-dvb_bar)**2)
   call hmatrix%apply(geom,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv((d-ytmp),ytmp2)
   algo%jo(ii) = 0.5*sum((d-ytmp)*ytmp2)

   ! Compute nonlinear cost function
   call bmatrix%apply_inv(geom,xg+algo%dx(:,ii)-xb,xtmp)
   algo%jb_nl(ii) = 0.5*sum((xg+algo%dx(:,ii)-xb)*xtmp)
   call hmatrix%apply(geom,(xg+algo%dx(:,ii)),ytmp)
   ytmp = ytmp-hmatrix%yo
   call rmatrix%apply_inv(ytmp,ytmp2)
   algo%jo_nl(ii) = 0.5*sum(ytmp*ytmp2)

   ! Stop criterion:
   if (shutoff_type==1) then
      if (algo%jb(ii)/algo%jb(ii-1)<shutoff_value) then
         write(*,'(a)') 'Convergence reached'
         stop
      end if
   end if
   ! Check cost function
   ! if (algo%jb(ii)+algo%jo(ii)>algo%jb(ii-1)+algo%jo(ii-1)) then
   !    write(*,'(a)') 'ERROR: increasing cost function in Lanczos'
   !    stop
   ! end if
end do

! Final update
algo%dva = u(:,ni)

! Lanczos vectors
algo%lancvec = v(:,1:ni+1)

! Last beta value
algo%lastbeta = beta(ni+1)

! Lanczos to CG
beta(0) = sqrt(sum(r(:,0)**2))
p(:,0) = r(:,0)
rho(0) = sum(r(:,0)**2)
algo%rho_sqrt(0) = sqrt(rho(0))
algo%beta(0) = beta(0)

do ii=1,ni
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
subroutine algo_apply_planczosif(algo,geom,bmatrix,hmatrix,rmatrix,dxbbar,d,ni,lmp,dxabar,shutoff_type,shutoff_value)

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
real(8),intent(out) :: dxabar(geom%nh)
integer,intent(in)  :: shutoff_type
real(8),intent(in)  ::shutoff_value

! Local variables
integer :: ii,info
real(8) :: xtmp(geom%nh),ytmp(hmatrix%nobs),ytmp2(hmatrix%nobs),y(ni,ni),rhs(ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: dxb(geom%nh),dxbar(geom%nh,ni),l(geom%nh,0:ni),p(geom%nh,0:ni),pbar(geom%nh,0:ni),q(geom%nh,0:ni),r(geom%nh,0:ni),s(geom%nh,0:ni),tbar(geom%nh,0:ni),t(geom%nh,0:ni),v(geom%nh,0:ni+1),w(geom%nh,0:ni),zbar(geom%nh,0:ni+1),z(geom%nh,0:ni+1)

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

! Initialize cost function
call bmatrix%apply(geom,dxbbar,dxb)
algo%jb(0) = 0.5*sum((0.0-dxb)*(0.0-dxbbar))
call rmatrix%apply_inv(d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)

do ii=1,ni
   ! Update
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

   ! Stop criterion:
   if (shutoff_type==2) then
      if (beta(ii+1)<shutoff_value) then
         write(*,'(a)') 'Convergence reached'
         stop
      end if
   end if

   ! Compute eigenpairs
   if (ii==1) then
      algo%eigenval(1) = alpha(1)
      algo%eigenvec(1,1) = 1.0
   else
      algo%eigenval(1:ii) = alpha(1:ii)
      subdiag(1:ii-1) = beta(2:ii)
      call dsteqr('I',ii,algo%eigenval(1:ii),subdiag(1:ii-1),algo%eigenvec(1:ii,1:ii),ii,work(1:2*ii-2),info)
      if (info/=0) then
         write(*,'(a,i2)') 'Error in dsteqr: ',info
         stop
      end if
      algo%eigenval(1:ii) = max(algo%eigenval(1:ii),1.0)
   end if

   ! Update increment
   y(1:ii,ii) = matmul(algo%eigenvec(1:ii,1:ii),matmul(transpose(algo%eigenvec(1:ii,1:ii)),rhs(1:ii))/algo%eigenval(1:ii))
   s(:,ii) = matmul(z(:,1:ii),y(1:ii,ii))
   dxbar(:,ii) = matmul(zbar(:,1:ii),y(1:ii,ii))
   algo%dx(:,ii) = s(:,ii)

   ! Compute cost function
   algo%jb(ii) = 0.5*sum((algo%dx(:,ii)-dxb)*(dxbar(:,ii)-dxbbar))
   call hmatrix%apply(geom,algo%dx(:,ii),ytmp)
   call rmatrix%apply_inv(d-ytmp,ytmp2)
   algo%jo(ii) = 0.5*sum((d-ytmp)*ytmp2)

   ! Stop criterion:
   if (shutoff_type==1) then
      if (algo%jb(ii)/algo%jb(ii-1)<shutoff_value) then
         write(*,'(a)') 'Convergence reached'
         stop
      end if
   end if

   ! Check cost function
   ! if (algo%jb(ii)+algo%jo(ii)>algo%jb(ii-1)+algo%jo(ii-1)) then
   !    write(*,'(a)') 'ERROR: increasing cost function in PLanczosIF'
   !    stop
   ! end if
end do

! Final update
dxabar = dxbar(:,ni)

! Lanczos vectors
algo%lancvec = v(:,1:ni+1)

! Last beta value
algo%lastbeta = beta(ni+1)

! PLanczosIF to PCGIF
call bmatrix%apply(geom,r(:,0),l(:,0))
call lmp%apply(geom,ni,lmp%io,r(:,0),zbar(:,0))
call lmp%apply_ad(geom,ni,lmp%io,l(:,0),z(:,0))
rho(0) = sum(r(:,0)*z(:,0))
beta(0) = sqrt(sum(r(:,0)*z(:,0)))
pbar(:,0) = zbar(:,0)
p(:,0) = z(:,0)

algo%rho_sqrt(0)=sqrt(rho(0))
algo%beta(0)=beta(0)

do ii=1,ni
   r(:,ii) = -beta(ii+1)*y(ii,ii)*v(:,ii+1)
   call bmatrix%apply(geom,r(:,ii),l(:,ii))
   call lmp%apply(geom,ni,lmp%io,r(:,ii),zbar(:,ii))
   call lmp%apply_ad(geom,ni,lmp%io,l(:,ii),z(:,ii))
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

end module type_algo
