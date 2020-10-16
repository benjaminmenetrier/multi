!----------------------------------------------------------------------
! Module: type_algo
! Purpose: lanczos and planczosif algorithms
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module type_algo

use type_bmatrix
use type_hmatrix
use type_lmp
use type_rmatrix

implicit none

type algo_type
   real(8),allocatable :: dx(:,:)
   real(8),allocatable :: jb(:)
   real(8),allocatable :: jo(:)
   real(8),allocatable :: eigenval(:)
   real(8),allocatable :: eigenvec(:,:)
   real(8),allocatable :: lancvec(:,:)
   real(8) :: lastbeta
   real(8),allocatable :: rho_sqrt(:)
   real(8),allocatable :: beta(:)
end type algo_type

contains

!----------------------------------------------------------------------
! Subroutine: algo_alloc
! Purpose: allocate algo
!----------------------------------------------------------------------
subroutine algo_alloc(algo,nn,ni)

implicit none

! Passed variables
type(algo_type),intent(inout) :: algo
integer,intent(in) :: nn
integer,intent(in) :: ni

! Release memory
if (allocated(algo%dx)) deallocate(algo%dx)
if (allocated(algo%jb)) deallocate(algo%jb)
if (allocated(algo%jo)) deallocate(algo%jo)
if (allocated(algo%eigenval)) deallocate(algo%eigenval)
if (allocated(algo%eigenvec)) deallocate(algo%eigenvec)
if (allocated(algo%lancvec)) deallocate(algo%lancvec)

if (allocated(algo%rho_sqrt)) deallocate(algo%rho_sqrt)
if (allocated(algo%beta)) deallocate(algo%beta)

! Allocation
allocate(algo%dx(nn,ni))
allocate(algo%jb(0:ni))
allocate(algo%jo(0:ni))
allocate(algo%eigenval(ni))
allocate(algo%eigenvec(ni,ni))
allocate(algo%lancvec(nn,ni+1))

allocate(algo%rho_sqrt(0:ni))
allocate(algo%beta(0:ni))

end subroutine algo_alloc

        
!----------------------------------------------------------------------
! Subroutine: algo_apply_lanczos
! Purpose: Lanczos algorithm in control space
!----------------------------------------------------------------------
subroutine algo_apply_lanczos(algo,nn,bmatrix,hmatrix,rmatrix,dvb,nobs,d,ni,lmp,dva,shutoff_type,shutoff_value)

implicit none

! Passed variables
type(algo_type),intent(inout) :: algo
integer,intent(in) :: nn
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: dvb(nn)
integer,intent(in) :: nobs
real(8),intent(in) :: d(nobs)
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
real(8),intent(out) :: dva(nn)
integer,intent(in)  :: shutoff_type
real(8),intent(in)  :: shutoff_value

! Local variables
integer :: ii,info
real(8) :: xtmp(nn),ytmp(nobs),ytmp2(nobs),y(ni,ni),rhs(ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: vtmp(nn),vtmp2(nn)
real(8) :: p(nn,0:ni),q(nn,0:ni),r(nn,0:ni),s(nn,0:ni),u(nn,ni),v(nn,0:ni+1),w(nn,0:ni)

! Initialization
s(:,0) = 0.0
v(:,0) = 0.0
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
call hmatrix_apply_ad(hmatrix,nobs,ytmp,nn,xtmp)
call bmatrix_apply_sqrt_ad(bmatrix,nn,xtmp,vtmp)
vtmp = dvb+vtmp
call lmp_apply_sqrt_ad(lmp,nn,ni,lmp%io,vtmp,r(:,0))
beta(0) = sqrt(sum(r(:,0)**2))
v(:,1) = r(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize cost function
algo%jb(0) = 0.5*sum((0.0-dvb)**2)
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)

do ii=1,ni
   ! Update
   call lmp_apply_sqrt(lmp,nn,ni,lmp%io,v(:,ii),vtmp)
   call bmatrix_apply_sqrt(bmatrix,nn,vtmp,xtmp)
   call hmatrix_apply(hmatrix,nn,xtmp,nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
   call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn,xtmp)
   call bmatrix_apply_sqrt_ad(bmatrix,nn,xtmp,vtmp2)
   vtmp = vtmp+vtmp2
   call lmp_apply_sqrt_ad(lmp,nn,ni,lmp%io,vtmp,q(:,ii))
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
   s(:,ii) = matmul(v(:,1:ii),y(1:ii,ii))
   call lmp_apply_sqrt(lmp,nn,ni,lmp%io,s(:,ii),u(:,ii))
   call bmatrix_apply_sqrt(bmatrix,nn,u(:,ii),algo%dx(:,ii))

   ! Compute cost function
   algo%jb(ii) = 0.5*sum((u(:,ii)-dvb)**2)
   call hmatrix_apply(hmatrix,nn,algo%dx(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,(d-ytmp),ytmp2)
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
   !    write(*,'(a)') 'ERROR: increasing cost function in Lanczos'
   !    stop
   ! end if
end do

! Final update
dva = u(:,ni)

! Lanczos vectors
algo%lancvec = v(:,1:ni+1)

! Last beta value
algo%lastbeta = beta(ni+1)

! ancienne version:
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
subroutine algo_apply_planczosif(algo,nn,bmatrix,hmatrix,rmatrix,dxbbar,nobs,d,ni,lmp,dxabar,shutoff_type,shutoff_value)

implicit none

! Passed variables
type(algo_type),intent(inout) :: algo
integer,intent(in) :: nn
type(bmatrix_type),intent(in) :: bmatrix
type(hmatrix_type),intent(in) :: hmatrix
type(rmatrix_type),intent(in) :: rmatrix
real(8),intent(in) :: dxbbar(nn)
integer,intent(in) :: nobs
real(8),intent(in) :: d(nobs)
integer,intent(in) :: ni
type(lmp_type),intent(in) :: lmp
real(8),intent(out) :: dxabar(nn)
integer,intent(in)  :: shutoff_type
real(8),intent(in)  ::shutoff_value

! Local variables
integer :: ii,info
real(8) :: xtmp(nn),ytmp(nobs),ytmp2(nobs),y(ni,ni),rhs(ni),subdiag(ni-1),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1),rho(0:ni)
real(8) :: dxb(nn),dxbar(nn,ni),l(nn,0:ni),p(nn,0:ni),pbar(nn,0:ni),q(nn,0:ni),r(nn,0:ni),s(nn,0:ni),tbar(nn,0:ni),t(nn,0:ni),v(nn,0:ni+1),w(nn,0:ni),zbar(nn,0:ni+1),z(nn,0:ni+1)

! Initialization
s(:,0) = 0.0
v(:,0) = 0.0
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
call hmatrix_apply_ad(hmatrix,nobs,ytmp,nn,xtmp)
r(:,0) = dxbbar+xtmp
call lmp_apply(lmp,nn,ni,lmp%io,r(:,0),tbar(:,0))
call bmatrix_apply(bmatrix,nn,tbar(:,0),t(:,0))
beta(0) = sqrt(sum(r(:,0)*t(:,0)))
v(:,1) = r(:,0)/beta(0)
zbar(:,1) = tbar(:,0)/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize cost function
call bmatrix_apply(bmatrix,nn,dxbbar,dxb)
algo%jb(0) = 0.5*sum((0.0-dxb)*(0.0-dxbbar))
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
algo%jo(0) = 0.5*sum(d*ytmp)

do ii=1,ni
   ! Update
   call hmatrix_apply(hmatrix,nn,z(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
   call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn,xtmp)
   q(:,ii) = zbar(:,ii)+xtmp-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call lmp_apply(lmp,nn,ni,lmp%io,w(:,ii),tbar(:,ii))
   call bmatrix_apply(bmatrix,nn,tbar(:,ii),t(:,ii))
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
   call hmatrix_apply(hmatrix,nn,algo%dx(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,d-ytmp,ytmp2)
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
call bmatrix_apply(bmatrix,nn,r(:,0),l(:,0))
call lmp_apply(lmp,nn,ni,lmp%io,r(:,0),zbar(:,0))
call lmp_apply_ad(lmp,nn,ni,lmp%io,l(:,0),z(:,0))
rho(0) = sum(r(:,0)*z(:,0))
beta(0) = sqrt(sum(r(:,0)*z(:,0)))
pbar(:,0) = zbar(:,0)
p(:,0) = z(:,0)

algo%rho_sqrt(0)=sqrt(rho(0))
algo%beta(0)=beta(0)

do ii=1,ni
   r(:,ii) = -beta(ii+1)*y(ii,ii)*v(:,ii+1)
   call bmatrix_apply(bmatrix,nn,r(:,ii),l(:,ii))
   call lmp_apply(lmp,nn,ni,lmp%io,r(:,ii),zbar(:,ii))
   call lmp_apply_ad(lmp,nn,ni,lmp%io,l(:,ii),z(:,ii))
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
