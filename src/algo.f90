!----------------------------------------------------------------------
! Module: algo
! Purpose: PLanczos and PLanczosIF algorithms, preconditioner application
! Author: Benjamin Menetrier
! Licensing: this code is distributed under the CeCILL-C license
! Copyright © 2018 IRIT
!----------------------------------------------------------------------
module algo

use bmatrix
use hmatrix

implicit none

private
public :: planczos,planczosif

contains

!----------------------------------------------------------------------
! Subroutine: planczos
! Purpose: preconditioned Lanczos
!----------------------------------------------------------------------
subroutine planczos(nn,sigmab,spvar,sigmao,dvb,nobs,d,ni,io,nprec,precvec,dva,dx,jb,jo,ritzval,vritzvec)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: sigmao
real(8),intent(in) :: dvb(nn)
integer,intent(in) :: nobs
real(8),intent(in) :: d(nobs)
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nprec(io)
real(8),intent(in) :: precvec(nn,ni,io,2)
real(8),intent(out) :: dva(nn)
real(8),intent(out) :: dx(nn,ni)
real(8),intent(out) :: jb(0:ni)
real(8),intent(out) :: jo(0:ni)
real(8),intent(out) :: ritzval(ni)
real(8),intent(out) :: vritzvec(nn,ni)

! Local variables
integer :: ii,ji,ierr
real(8) :: xtmp(nn),ytmp(nobs),y(ni),rhs(ni),subdiag(ni-1),eigenval(ni),eigenvec(ni,ni),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1)
real(8) :: vtmp(nn)
real(8) :: q(nn,ni),r0(nn),t(nn,0:ni),u(nn,ni),v(nn,0:ni+1),w(nn,0:ni),z(nn,0:ni+1)

! Initialization
v(:,0) = 0.0
ytmp = d/sigmao**2
call apply_h_ad(nobs,ytmp,nn,xtmp)
call apply_u_ad(nn,sigmab,spvar,xtmp,vtmp)
r0 = dvb+vtmp
call apply_lmp(nn,ni,io,nprec,precvec,r0,t(:,0))
beta(0) = sqrt(sum(t(:,0)*r0))
v(:,1) = r0/beta(0)
z(:,1) = t(:,0)/beta(0)
beta(1) = 0.0
rhs = 0.0
rhs(1) = beta(0)

! Initialize cost function
jb(0) = 0.5*sum((0.0-dvb)*(0.0-dvb))
ytmp = d/sigmao
jo(0) = 0.5*sum(ytmp*ytmp)

do ii=1,ni
   ! Update
   call apply_u(nn,sigmab,spvar,z(:,ii),xtmp)
   call apply_h(nn,xtmp,nobs,ytmp)
   ytmp = ytmp/sigmao**2
   call apply_h_ad(nobs,ytmp,nn,xtmp)
   call apply_u_ad(nn,sigmab,spvar,xtmp,vtmp)
   q(:,ii) = z(:,ii)+vtmp-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call apply_lmp(nn,ni,io,nprec,precvec,w(:,ii),t(:,ii))
   beta(ii+1) = sqrt(sum(t(:,ii)*w(:,ii)))
   v(:,ii+1) = w(:,ii)/beta(ii+1)
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
   u(:,ii) = matmul(z(:,1:ii),y(1:ii))
   call apply_u(nn,sigmab,spvar,u(:,ii),dx(:,ii))

   ! Compute cost function
   jb(ii) = 0.5*sum((u(:,ii)-dvb)*(u(:,ii)-dvb))
   call apply_h(nn,dx(:,ii),nobs,ytmp)
   ytmp = (d-ytmp)/sigmao
   jo(ii) = 0.5*sum(ytmp*ytmp)

   ! Check cost function
   if (jb(ii)+jo(ii)<jb(ii-1)+jo(ii-1)) then
      ! Final update
      dva = u(:,ii)

      ! Ritz pairs
      ritzval(1:ii) = eigenval(1:ii)
      vritzvec(:,1:ii) = matmul(z(:,1:ii),eigenvec(1:ii,1:ii))
   else
      ! Set remaining values
      do ji=ii,ni
         dx(:,ji) = dx(:,ii)
         jb(ji) = jb(ii)
         jo(ji) = jo(ii)
         ritzval(ji) = 1.0
         vritzvec(:,ji) = 0.0
      end do
      exit
   end if
end do

end subroutine planczos

!----------------------------------------------------------------------
! Subroutine: planczosif
! Purpose: preconditioned Lanczos "inverse-free"
!----------------------------------------------------------------------
subroutine planczosif(nn,sigmab,spvar,sigmao,dxbbar,nobs,d,ni,io,nprec,precvec,dxabar,dx,jb,jo,ritzval,dxbarritzvec,dxritzvec)

implicit none

! Passed variables
integer,intent(in) :: nn
real(8),intent(in) :: sigmab(nn)
real(8),intent(in) :: spvar(nn)
real(8),intent(in) :: sigmao
real(8),intent(in) :: dxbbar(nn)
integer,intent(in) :: nobs
real(8),intent(in) :: d(nobs)
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nprec(io)
real(8),intent(in) :: precvec(nn,ni,io,2)
real(8),intent(out) :: dxabar(nn)
real(8),intent(out) :: dx(nn,ni)
real(8),intent(out) :: jb(0:ni)
real(8),intent(out) :: jo(0:ni)
real(8),intent(out) :: ritzval(ni)
real(8),intent(out) :: dxbarritzvec(nn,ni)
real(8),intent(out) :: dxritzvec(nn,ni)

! Local variables
integer :: ii,ji,ierr
real(8) :: xtmp(nn),ytmp(nobs),y(ni),rhs(ni),subdiag(ni-1),eigenval(ni),eigenvec(ni,ni),work(max(1,2*ni-2))
real(8) :: alpha(0:ni),beta(0:ni+1)
real(8) :: dxb(nn),dxbar(nn,ni),q(nn,ni),r0(nn),tbar(nn,0:ni),t(nn,0:ni),v(nn,0:ni+1),w(nn,0:ni),zbar(nn,0:ni+1),z(nn,0:ni+1)

! Initialization
ytmp = d/sigmao**2
call apply_h_ad(nobs,ytmp,nn,xtmp)
r0 = dxbbar+xtmp
call apply_lmp(nn,ni,io,nprec,precvec,r0,tbar(:,0))
call apply_b(nn,sigmab,spvar,tbar(:,0),t(:,0))
beta(0) = sqrt(sum(r0*t(:,0)))
beta(1) = 0.0
v(:,0) = 0.0
v(:,1) = r0/beta(0)
zbar(:,1) = tbar(:,0)/beta(0)
z(:,1) = t(:,0)/beta(0)
rhs = 0.0
rhs(1) = beta(0)
call apply_b(nn,sigmab,spvar,dxbbar,dxb)

! Initialize cost function
jb = 0.0
jo = 0.0
jb(0) = 0.5*sum((0.0-dxb)*(0.0-dxbbar))
ytmp = d/sigmao
jo(0) = 0.5*sum(ytmp*ytmp)

do ii=1,ni
   ! Update
   call apply_h(nn,z(:,ii),nobs,ytmp)
   ytmp = ytmp/sigmao**2
   call apply_h_ad(nobs,ytmp,nn,xtmp)
   q(:,ii) = zbar(:,ii)+xtmp-beta(ii)*v(:,ii-1)
   alpha(ii) = sum(q(:,ii)*z(:,ii))
   w(:,ii) = q(:,ii)-alpha(ii)*v(:,ii)
   call apply_lmp(nn,ni,io,nprec,precvec,w(:,ii),tbar(:,ii))
   call apply_b(nn,sigmab,spvar,tbar(:,ii),t(:,ii))
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
   call apply_h(nn,dx(:,ii),nobs,ytmp)
   ytmp = (d-ytmp)/sigmao
   jo(ii) = 0.5*sum(ytmp*ytmp)

   ! Check cost function
   if (jb(ii)+jo(ii)<jb(ii-1)+jo(ii-1)) then
      ! Final update
      dxabar = dxbar(:,ii)

      ! Ritz pairs
      ritzval(1:ii) = eigenval(1:ii)
      dxbarritzvec(:,1:ii) = matmul(zbar(:,1:ii),eigenvec(1:ii,1:ii))
      dxritzvec(:,1:ii) = matmul(z(:,1:ii),eigenvec(1:ii,1:ii))
   else
      ! Set remaining values
      do ji=ii,ni
         dx(:,ji) = dx(:,ii)
         jb(ji) = jb(ii)
         jo(ji) = jo(ii)
         ritzval(ji) = 1.0
         dxbarritzvec(:,ji) = 0.0
         dxritzvec(:,ji) = 0.0
      end do
      exit
   end if
end do

end subroutine planczosif

!----------------------------------------------------------------------
! Subroutine: apply_lmp
! Purpose: apply LMP
!----------------------------------------------------------------------
subroutine apply_lmp(nn,ni,io,nprec,precvec,x,px)

implicit none

! Passed variables
integer,intent(in) :: nn
integer,intent(in) :: ni
integer,intent(in) :: io
integer,intent(in) :: nprec(io)
real(8),intent(in) :: precvec(nn,ni,io,2)
real(8),intent(in) :: x(nn)
real(8),intent(out) :: px(nn)

! Local variables
integer :: jo,iprec

! Initialization
px = x

! Previous outer loops
do jo=1,io-1
   do iprec=1,nprec(jo)
      px = px-sum(precvec(:,iprec,jo,1)*x)*precvec(:,iprec,jo,2)
   end do
end do

end subroutine apply_lmp

end module algo
