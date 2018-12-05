!----------------------------------------------------------------------
! Module: algo
! Purpose: lanczos and planczosif algorithms, preconditioner application
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module algo

use type_bmatrix
use type_hmatrix
use type_lmp
use type_rmatrix

implicit none

contains

!----------------------------------------------------------------------
! Subroutine: lanczos
! Purpose: preconditioned Lanczos
!----------------------------------------------------------------------
subroutine lanczos(nn,bmatrix,hmatrix,rmatrix,dvb,nobs,d,ni,lmp,dva,dx,jb,jo,ritzval,ritzvec)

implicit none

! Passed variables
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
real(8),intent(out) :: dx(nn,ni)
real(8),intent(out) :: jb(0:ni)
real(8),intent(out) :: jo(0:ni)
real(8),intent(out) :: ritzval(ni)
real(8),intent(out) :: ritzvec(nn,ni)

! Local variables
integer :: ii,ji,ierr
real(8) :: xtmp(nn),ytmp(nobs),ytmp2(nobs),y(ni),rhs(ni),subdiag(ni-1),eigenval(ni),eigenvec(ni,ni),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1)
real(8) :: vtmp(nn),vtmp2(nn)
real(8) :: q(nn,ni),r0(nn),s(nn,ni),u(nn,ni),v(nn,0:ni+1),w(nn,0:ni)

! Initialization
v(:,0) = 0.0
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
call hmatrix_apply_ad(hmatrix,nobs,ytmp,nn,xtmp)
call bmatrix_apply_sqrt_ad(bmatrix,nn,xtmp,vtmp)
vtmp = dvb+vtmp
call lmp_apply_sqrt_ad(lmp,nn,lmp%io,vtmp,r0)
beta(0) = sqrt(sum(r0**2))
v(:,1) = r0/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize cost function
jb(0) = 0.5*sum((0.0-dvb)**2)
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
jo(0) = 0.5*sum(d*ytmp)

do ii=1,ni
   ! Update
   call lmp_apply_sqrt(lmp,nn,lmp%io,v(:,ii),vtmp)
   call bmatrix_apply_sqrt(bmatrix,nn,vtmp,xtmp)
   call hmatrix_apply(hmatrix,nn,xtmp,nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
   call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn,xtmp)
   call bmatrix_apply_sqrt_ad(bmatrix,nn,xtmp,vtmp2)
   vtmp = vtmp+vtmp2
   call lmp_apply_sqrt_ad(lmp,nn,lmp%io,vtmp,q(:,ii))
   q(:,ii) = q(:,ii)-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*v(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   beta(ii+1) = sqrt(sum(w(:,ii)**2))
   v(:,ii+1) = w(:,ii)/beta(ii+1)

   ! Compute eigenpairs
   if (ii==1) then
      eigenval(1) = alpha(1)
      eigenvec(1,1) = 1.0
   else
      eigenval(1:ii) = alpha(1:ii)
      subdiag(1:ii-1) = beta(2:ii)
      call dsteqr('I',ii,eigenval(1:ii),subdiag(1:ii-1),eigenvec(1:ii,1:ii),ii,work(1:2*ii-2),ierr)
      if (ierr/=0) then
         write(*,'(a)') 'Error in dsteqr'
         stop
      end if
      eigenval(1:ii) = max(eigenval(1:ii),1.0)
   end if

   ! Update increment
   y(1:ii) = matmul(eigenvec(1:ii,1:ii),matmul(transpose(eigenvec(1:ii,1:ii)),rhs(1:ii))/eigenval(1:ii))
   s(:,ii) = matmul(v(:,1:ii),y(1:ii))
   call lmp_apply_sqrt(lmp,nn,lmp%io,s(:,ii),u(:,ii))
   call bmatrix_apply_sqrt(bmatrix,nn,u(:,ii),dx(:,ii))

   ! Compute cost function
   jb(ii) = 0.5*sum((u(:,ii)-dvb)**2)
   call hmatrix_apply(hmatrix,nn,dx(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,(d-ytmp),ytmp2)
   jo(ii) = 0.5*sum((d-ytmp)*ytmp2)

   ! Check cost function
   if (jb(ii)+jo(ii)<jb(ii-1)+jo(ii-1)) then
      ! Final update
      dva = u(:,ii)

      ! Ritz pairs
      ritzval(1:ii) = eigenval(1:ii)
      ritzvec(:,1:ii) = matmul(v(:,1:ii),eigenvec(1:ii,1:ii))
   else
      ! Set remaining values
      do ji=ii,ni
         dx(:,ji) = dx(:,ii)
         jb(ji) = jb(ii)
         jo(ji) = jo(ii)
         ritzval(ji) = 1.0
         ritzvec(:,ji) = 0.0
      end do
      exit
   end if
end do

end subroutine lanczos

!----------------------------------------------------------------------
! Subroutine: planczosif
! Purpose: preconditioned Lanczos "inverse-free"
!----------------------------------------------------------------------
subroutine planczosif(nn,bmatrix,hmatrix,rmatrix,dxbbar,nobs,d,ni,lmp,dxabar,dx,jb,jo,ritzval,ritzvec)

implicit none

! Passed variables
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
real(8),intent(out) :: dx(nn,ni)
real(8),intent(out) :: jb(0:ni)
real(8),intent(out) :: jo(0:ni)
real(8),intent(out) :: ritzval(ni)
real(8),intent(out) :: ritzvec(nn,ni)

! Local variables
integer :: ii,ji,ierr
real(8) :: xtmp(nn),ytmp(nobs),ytmp2(nobs),y(ni),rhs(ni),subdiag(ni-1),eigenval(ni),eigenvec(ni,ni),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1)
real(8) :: dxb(nn),dxbar(nn,ni),q(nn,ni),r0(nn),tbar(nn,0:ni),t(nn,0:ni),v(nn,0:ni+1),w(nn,0:ni),zbar(nn,0:ni+1),z(nn,0:ni+1)

! Initialization
v(:,0) = 0.0
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
call hmatrix_apply_ad(hmatrix,nobs,ytmp,nn,xtmp)
r0 = dxbbar+xtmp
call lmp_apply(lmp,nn,lmp%io,r0,tbar(:,0))
call bmatrix_apply(bmatrix,nn,tbar(:,0),t(:,0))
beta(0) = sqrt(sum(r0*t(:,0)))
v(:,1) = r0/beta(0)
zbar(:,1) = tbar(:,0)/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize cost function
call bmatrix_apply(bmatrix,nn,dxbbar,dxb)
jb(0) = 0.5*sum((0.0-dxb)*(0.0-dxbbar))
call rmatrix_apply_inv(rmatrix,nobs,d,ytmp)
jo(0) = 0.5*sum(d*ytmp)

do ii=1,ni
   ! Update
   call hmatrix_apply(hmatrix,nn,z(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,ytmp,ytmp2)
   call hmatrix_apply_ad(hmatrix,nobs,ytmp2,nn,xtmp)
   q(:,ii) = zbar(:,ii)+xtmp-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call lmp_apply(lmp,nn,lmp%io,w(:,ii),tbar(:,ii))
   call bmatrix_apply(bmatrix,nn,tbar(:,ii),t(:,ii))
   beta(ii+1) = sqrt(sum(w(:,ii)*t(:,ii)))
   v(:,ii+1) = w(:,ii)/beta(ii+1)
   zbar(:,ii+1) = tbar(:,ii)/beta(ii+1)
   z(:,ii+1) = t(:,ii)/beta(ii+1)

   ! Compute eigenpairs
   if (ii==1) then
      eigenval(1) = alpha(1)
      eigenvec(1,1) = 1.0
   else
      eigenval(1:ii) = alpha(1:ii)
      subdiag(1:ii-1) = beta(2:ii)
      call dsteqr('I',ii,eigenval(1:ii),subdiag(1:ii-1),eigenvec(1:ii,1:ii),ii,work(1:2*ii-2),ierr)
      if (ierr/=0) then
         write(*,'(a)') 'Error in dsteqr'
         stop
      end if
      eigenval(1:ii) = max(eigenval(1:ii),1.0)
   end if

   ! Update increment
   y(1:ii) = matmul(eigenvec(1:ii,1:ii),matmul(transpose(eigenvec(1:ii,1:ii)),rhs(1:ii))/eigenval(1:ii))
   dxbar(:,ii) = matmul(zbar(:,1:ii),y(1:ii))
   dx(:,ii) = matmul(z(:,1:ii),y(1:ii))

   ! Compute cost function
   jb(ii) = 0.5*sum((dx(:,ii)-dxb)*(dxbar(:,ii)-dxbbar))
   call hmatrix_apply(hmatrix,nn,dx(:,ii),nobs,ytmp)
   call rmatrix_apply_inv(rmatrix,nobs,d-ytmp,ytmp2)
   jo(ii) = 0.5*sum((d-ytmp)*ytmp2)

   ! Check cost function
   if (jb(ii)+jo(ii)<jb(ii-1)+jo(ii-1)) then
      ! Final update
      dxabar = dxbar(:,ii)

      ! Ritz pairs
      ritzval(1:ii) = eigenval(1:ii)
      ritzvec(:,1:ii) = matmul(v(:,1:ii),eigenvec(1:ii,1:ii))
   else
      ! Set remaining values
      do ji=ii,ni
         dx(:,ji) = dx(:,ii)
         jb(ji) = jb(ii)
         jo(ji) = jo(ii)
         ritzval(ji) = 1.0
         ritzvec(:,ji) = 0.0
      end do
      exit
   end if
end do

end subroutine planczosif

end module algo
