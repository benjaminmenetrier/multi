program main

use netcdf
use tools_netcdf
use tools_rand
use type_algo
use type_bmatrix
use type_geom
use type_hmatrix
use type_lmp
use type_rmatrix

implicit none

! Parameters
integer,parameter   :: nomax = 9           ! Maximum number of outer iterations
integer,parameter   :: namelist_unit = 9   ! Namelist unit

! Namelist parameters
integer             :: no                  ! Number of outer iterations
integer             :: ni                  ! Number of inner iterations
character(len=8)    :: lmp_mode            ! LMP mode ('none', 'spectral', 'ritz')
logical             :: test_ortho          ! Test orthogonality
integer             :: shutoff_type        ! Stopping criterion according: 1-Jb, 2-beta, 3-Ritz convergence (else: no criterion)
real(8)             :: shutoff_value       ! Stopping criterion threshold
logical             :: consistent_strategy ! Use the consistent strategy detailed in the note
logical             :: use_binv            ! Use B inverse in the inconsistent strategy
logical             :: equiv_precond       ! Equivalent preconditionings (conditions on B and interpolations are valid)
integer             :: nx(nomax)           ! X direction sizes
integer             :: ny(nomax)           ! Y direction sizes
integer             :: nobs                ! Number of observations
real(8)             :: sigma_obs           ! Observation error standard deviation
real(8)             :: sigmabvar           ! Grid-point standard deviation variations amplitude
real(8)             :: Lb                  ! Correlation length-scale
real(8)             :: spvarmin            ! Minimum spectral variance (inverse fails if this is too small)
logical             :: new_seed            ! New random seed
character(len=1024) :: filename            ! Filename

! Local variables
integer                        :: ncid,subgrpid,nx_id,ny_id,nobs_id,xb_id,xg_id,hxg_id,d_id
integer                        :: io,jo,ii,ji
integer,allocatable            :: grpid(:)
real(8)                        :: proj,norm
real(8),allocatable            :: xb_full(:),xg_full(:),dxb_full(:),dxa_full(:)
real(8),allocatable            :: xb(:),xg(:),dxb(:),dxbbar(:),dvbbar(:),dxa_prev(:),dva_interp(:),dxabar_interp(:)
real(8),allocatable            :: xb_2d(:,:),xg_2d(:,:)
real(8),allocatable            :: d(:),hxg(:)
character(len=1024)            :: grpname
type(algo_type),allocatable    :: algo_lanczos(:),algo_planczosif(:)
type(bmatrix_type),allocatable :: bmatrix(:)
type(geom_type),allocatable    :: geom(:)
type(hmatrix_type)             :: hmatrix
type(lmp_type),allocatable     :: lmp_lanczos(:),lmp_planczosif(:)
type(rmatrix_type)             :: rmatrix

! Namelist blocks
namelist/solver/ &
 & no, &
 & ni, &
 & lmp_mode, &
 & test_ortho, &
 & shutoff_type, &
 & shutoff_value, &
 & consistent_strategy, &
 & use_binv, &
 & equiv_precond
namelist/resolutions/ &
 & nx, &
 & ny
namelist/observations/ &
 & nobs, &
 & sigma_obs
namelist/background/ &
 & sigmabvar, &
 & Lb, &
 & spvarmin
namelist/miscellanous/ &
 & new_seed, &
 & filename

! Initialize namelist (default values)
no = 1
ni = 10
lmp_mode = 'none'
test_ortho = .false.
shutoff_type = 0
shutoff_value = 0.0
consistent_strategy = .true.
use_binv = .false.
equiv_precond = .true.
nx = 101
ny = 101
nobs = 2000
sigma_obs = 0.1
sigmabvar = 0.0
Lb = 0.12
spvarmin = 1.0e-5
new_seed = .false.
filename = 'output'

! Read namelist
open(unit=namelist_unit,file='namelist',status='old',action='read')
read(namelist_unit,nml=solver)
read(namelist_unit,nml=resolutions)
read(namelist_unit,nml=observations)
read(namelist_unit,nml=background)
read(namelist_unit,nml=miscellanous)
close(unit=namelist_unit)

! Print namelist
write(*,'(a)') 'Namelist'
write(*,nml=solver)
write(*,nml=resolutions)
write(*,nml=observations)
write(*,nml=background)
write(*,nml=miscellanous)
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! General initialization
!--------------------------------------------------------------------------------
write(*,'(a)') 'General initialization'

! Set seed
call set_seed(new_seed)

! Allocation (number of outer iterations)
allocate(grpid(no))
allocate(geom(no))
allocate(bmatrix(no))
allocate(algo_lanczos(no))
allocate(algo_planczosif(no))
allocate(lmp_lanczos(no))
allocate(lmp_planczosif(no))

! Create NetCDF file
call ncerr('main',nf90_create(trim(filename)//'.nc',ior(nf90_clobber,nf90_netcdf4),ncid))

! Set missing value
call ncerr('main',nf90_put_att(ncid,nf90_global,'_FillValue',-999.0))

! Create groups for each outer iteration
do io=1,no
   write(grpname,'(a,i1)') 'outer_',io
   call ncerr('main',nf90_def_grp(ncid,grpname,grpid(io)))
end do

! Setup geometries
do io=1,no
   write(*,'(a,i1)') '   Geometry setup for outer iteration ',io
   call geom(io)%setup(nx(io),ny(io),equiv_precond)
   call geom(io)%write(grpid(io))
end do

! Allocation (model/control space)
allocate(xb_full(geom(no)%nh))
allocate(xg_full(geom(no)%nh))

! Setup B matrices
do io=1,no
   write(*,'(a,i1)') '   B matrix setup for outer iteration ',io
   if (geom(io)%equiv_precond) then
      call bmatrix(io)%setup(geom(io),geom(no),sigmabvar,Lb,spvarmin)
   else
      call bmatrix(io)%setup(geom(io),geom(io),sigmabvar,Lb,spvarmin)
   end if
   call bmatrix(io)%write(geom(io),grpid(io))
end do

! Generate background state
call bmatrix(no)%randomize(geom(no),xb_full)
do io=1,no
   ! Allocation
   allocate(xb(geom(io)%nh))
   allocate(xb_2d(geom(io)%nx,geom(io)%ny))

   ! Interpolate background at current resolution
   call geom(no)%interp_gp(geom(io),xb_full,xb)

   ! Reshape background
   xb_2d = reshape(xb,(/geom(io)%nx,geom(io)%ny/))

   ! Get dimensions
   call ncerr('main',nf90_inq_dimid(grpid(io),'nx',nx_id))
   call ncerr('main',nf90_inq_dimid(grpid(io),'ny',ny_id))

   ! Create variable
   call ncerr('main',nf90_def_var(grpid(io),'xb',nf90_double,(/nx_id,ny_id/),xb_id))

   ! Write variable
   call ncerr('main',nf90_put_var(grpid(io),xb_id,xb_2d))

   ! Release memory
   deallocate(xb)
   deallocate(xb_2d)
end do

! Allocation (obs space)
allocate(d(nobs))
allocate(hxg(nobs))

! Setup obserations locations
call hmatrix%setup(nobs)
call hmatrix%write(ncid)

! Setup R matrix
call rmatrix%setup(nobs,sigma_obs)

! Setup observations
call rmatrix%randomize(hmatrix%yo)

write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Multi-incremental Lanczos in control space
!--------------------------------------------------------------------------------
write(*,'(a)') 'Multi-incremental Lanczos in control space'
do io=1,no
   write(*,'(a,i1,a,i4,a,i4)') '   Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny

   ! Allocation
   allocate(xb(geom(io)%nh))
   allocate(xg(geom(io)%nh))
   allocate(xg_2d(geom(io)%nx,geom(io)%ny))
   allocate(dvbbar(geom(io)%nh))
   call algo_lanczos(io)%alloc(geom(io),ni)
   call lmp_lanczos(io)%alloc(no,geom,ni,io,lmp_mode,'control')

   ! Test H matrix
   write(*,'(a)') '      Test H matrix'
   call hmatrix%test(geom(io))

   ! Interpolate background at current resolution
   call geom(no)%interp_gp(geom(io),xb_full,xb)

   if (consistent_strategy) then
      ! Allocation
      allocate(dva_interp(geom(io)%nh))
      allocate(dxb(geom(io)%nh))
      allocate(dxb_full(geom(no)%nh))

      ! Initialization
      dvbbar = 0.0

      do jo=1,io-1
         ! Interpolate analysis increment of previous iterations at current resolution
         call geom(jo)%interp_sp(geom(io),algo_lanczos(jo)%dva,dva_interp)

         ! Add contributions
         dvbbar = dvbbar-dva_interp
      end do

      ! Compute background increment at current resolution
      call bmatrix(io)%apply_sqrt(geom(io),dvbbar,dxb)

      ! Compute guess at current resolution
      xg = xb-dxb

      ! Interpolate background increment at full resolution
      call geom(io)%interp_gp(geom(no),dxb,dxb_full)

      ! Compute guess at full resolution
      xg_full = xb_full-dxb_full

      ! Release memory
      deallocate(dva_interp)
      deallocate(dxb)
      deallocate(dxb_full)
   else
      ! Compute guess at full resolution
      if (io==1) then
         ! Background at full resolution
         xg_full = xb_full
      else
         ! Allocation
         allocate(dxa_prev(geom(io-1)%nh))
         allocate(dxa_full(geom(no)%nh))

         ! Compute analysis increment of the previous outer iteration
         call bmatrix(io-1)%apply_sqrt(geom(io-1),algo_lanczos(io-1)%dva,dxa_prev)

         ! Interpolate analysis increment of the previous outer iteration at full resolution
         call geom(io-1)%interp_gp(geom(no),dxa_prev,dxa_full)

         ! Add analysis increment of the previous outer iteration at full resolution
         xg_full = xg_full+dxa_full
 
         ! Release memory
         deallocate(dxa_prev)
         deallocate(dxa_full)
      end if

      ! Interpolate guess at current resolution
      call geom(no)%interp_gp(geom(io),xg_full,xg)

      ! Right-hand side
      if (use_binv) then
         ! Allocation
         allocate(dxb(geom(io)%nh))
         allocate(dxbbar(geom(io)%nh))

         ! Compute background increment at current resolution
         dxb = xg-xb

         ! Apply B inverse
         call bmatrix(io)%apply_inv(geom(io),dxb,dxbbar)

         ! Apply B square-root
         call bmatrix(io)%apply_sqrt_ad(geom(io),dxbbar,dvbbar)

         ! Release memory
         deallocate(dxb)
         deallocate(dxbbar)
      else
         ! Allocation
         allocate(dva_interp(geom(io)%nh))

         ! Initialization
         dvbbar = 0.0

         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_sp(geom(io),algo_lanczos(jo)%dva,dva_interp)

            ! Add contributions
            dvbbar = dvbbar-dva_interp
         end do

         ! Release memory
         deallocate(dva_interp)
      end if
   end if

   ! Compute innovation
   call hmatrix%apply(geom(io),xg,hxg)
   d = hmatrix%yo-hxg

   ! Copy and interpolate LMP vectors
   if (io>1) then
      select case (trim(lmp_mode))
      case ('spectral','ritz')
         ! Copy eigenpairs
         lmp_lanczos(io)%outer(io)%eigenval = algo_lanczos(io-1)%eigenval
         lmp_lanczos(io)%outer(io)%eigenvec = algo_lanczos(io-1)%eigenvec

         ! Copy Lanczos vectors
         lmp_lanczos(io)%lancvec_trunc = algo_lanczos(io-1)%lancvec

         ! Compute omega
         lmp_lanczos(io)%outer(io)%omega = algo_lanczos(io-1)%eigenvec(ni,:) &
 & *algo_lanczos(io-1)%lastbeta/algo_lanczos(io-1)%eigenval

         do jo=2,io
            ! Copy eigenpairs and omega
            lmp_lanczos(io)%outer(jo)%eigenval = lmp_lanczos(jo)%outer(jo)%eigenval
            lmp_lanczos(io)%outer(jo)%eigenvec = lmp_lanczos(jo)%outer(jo)%eigenvec
            lmp_lanczos(io)%outer(jo)%omega = lmp_lanczos(jo)%outer(jo)%omega

            do ii=1,ni+1
               ! Interpolation of Lanczos vectors
               call geom(jo-1)%interp_sp(geom(io),lmp_lanczos(jo)%lancvec_trunc(:,ii),lmp_lanczos(io)%outer(jo)%lancvec(:,ii))

               if (test_ortho) then
                  ! Orthogonality test
                  do ji=1,ii-1
                     proj = sum(lmp_lanczos(io)%outer(jo)%lancvec(:,ji)*lmp_lanczos(io)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp_lanczos(io)%outer(jo)%lancvec(:,ji)*lmp_lanczos(io)%outer(jo)%lancvec(:,ji))
                     if (abs(proj)>1.0e-6) then
                        write(*,'(a)') 'ERROR: orthogonality lost'
                        stop
                     end if
                  end do
                  norm = sqrt(sum(lmp_lanczos(io)%outer(jo)%lancvec(:,ii)**2))
                  if (abs(norm-1.0)>1.0e-6) then
                     write(*,'(a)') 'ERROR: orthogonality lost '
                     stop
                  end if
               end if
            end do

            ! Ritz vectors
            lmp_lanczos(io)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp_lanczos(io)%outer(jo)%lancvec(:,1:ni), &
 & lmp_lanczos(io)%outer(jo)%eigenvec)
         end do
      end select
   end if

   ! Minimization
   write(*,'(a)') '      Minimization'
   call algo_lanczos(io)%apply_lanczos(geom(io),bmatrix(io),hmatrix,rmatrix,dvbbar,d,ni,lmp_lanczos(io), &
 & shutoff_type,shutoff_value)

   ! Compute nonlinear cost function
   call algo_lanczos(io)%compute_cost(geom(io),bmatrix(io),hmatrix,rmatrix,xb,xg,ni)

   ! Result
   write(*,'(a)') '      Cost function: J = Jb+Jo / J_nl = Jb_nl+Jo_nl'
   do ii=0,ni
      write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '         Iteration ',ii,':', &
 & algo_lanczos(io)%j(ii),' = ',algo_lanczos(io)%jb(ii),' + ',algo_lanczos(io)%jo(ii),' / ', &
 & algo_lanczos(io)%j_nl(ii),' = ',algo_lanczos(io)%jb_nl(ii),' + ',algo_lanczos(io)%jo_nl(ii)
   end do

   ! Create subgroup for this algo
   call ncerr('main',nf90_def_grp(grpid(io),'lanczos',subgrpid))

   ! Reshape guess
   xb_2d = reshape(xb,(/geom(io)%nx,geom(io)%ny/))

   ! Get dimensions
   call ncerr('main',nf90_inq_dimid(grpid(io),'nx',nx_id))
   call ncerr('main',nf90_inq_dimid(grpid(io),'ny',ny_id))
   call ncerr('main',nf90_inq_dimid(ncid,'nobs',nobs_id))

   ! Create variables
   call ncerr('main',nf90_def_var(subgrpid,'xg',nf90_double,(/nx_id,ny_id/),xg_id))
   call ncerr('main',nf90_def_var(subgrpid,'hxg',nf90_double,(/nobs_id/),hxg_id))
   call ncerr('main',nf90_def_var(subgrpid,'d',nf90_double,(/nobs_id/),d_id))

   ! Write variable
   call ncerr('main',nf90_put_var(subgrpid,xg_id,xg_2d))
   call ncerr('main',nf90_put_var(subgrpid,hxg_id,hxg))
   call ncerr('main',nf90_put_var(subgrpid,d_id,d))

   ! Write algo data
   call algo_lanczos(io)%write(geom(io),grpid(io),subgrpid)

   ! Release memory
   deallocate(xb)
   deallocate(xg)
   deallocate(xg_2d)
   deallocate(dvbbar)
end do
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Multi-incremental PLanczosIF in model space
!--------------------------------------------------------------------------------
write(*,'(a)') 'Multi-incremental PLanczosIF in model space'
do io=1,no
   write(*,'(a,i1,a,i4,a,i4)') '   Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny

   ! Allocation
   allocate(xb(geom(io)%nh))
   allocate(xg(geom(io)%nh))
   allocate(xg_2d(geom(io)%nx,geom(io)%ny))
   allocate(dxbbar(geom(io)%nh))
   call algo_planczosif(io)%alloc(geom(io),ni)
   call lmp_planczosif(io)%alloc(no,geom,ni,io,lmp_mode,'model')

   ! Test H matrix
   write(*,'(a)') '      Test H matrix'
   call hmatrix%test(geom(io))

   ! Interpolate background at current resolution
   call geom(no)%interp_gp(geom(io),xb_full,xb)

   if (consistent_strategy) then
      ! Allocation
      allocate(dxabar_interp(geom(io)%nh))
      allocate(dxb(geom(io)%nh))
      allocate(dxb_full(geom(no)%nh))

      ! Initialization
      dxbbar = 0.0

      do jo=1,io-1
         ! Interpolate analysis increment of previous iterations at current resolution
         call geom(jo)%interp_gp(geom(io),algo_planczosif(jo)%dxabar,dxabar_interp)

         ! Add contributions
         dxbbar = dxbbar-dxabar_interp
      end do

      ! Compute background increment at current resolution
      call bmatrix(io)%apply(geom(io),dxbbar,dxb)

      ! Compute guess at current resolution
      xg = xb-dxb

      ! Interpolate background increment at full resolution
      call geom(io)%interp_gp(geom(no),dxb,dxb_full)

      ! Compute guess at full resolution
      xg_full = xb_full-dxb_full

      ! Release memory
      deallocate(dxabar_interp)
      deallocate(dxb)
      deallocate(dxb_full)
   else
      ! Compute guess at full resolution
      if (io==1) then
         ! Background at full resolution
         xg_full = xb_full
      else
         ! Allocation
         allocate(dxa_prev(geom(io-1)%nh))
         allocate(dxa_full(geom(no)%nh))

         ! Compute analysis increment of the previous outer iteration
         call bmatrix(io-1)%apply(geom(io-1),algo_planczosif(io-1)%dxabar,dxa_prev)

         ! Interpolate analysis increment of the previous outer iteration at full resolution
         call geom(io-1)%interp_gp(geom(no),dxa_prev,dxa_full)

         ! Add analysis increment of the previous outer iteration at full resolution
         xg_full = xg_full+dxa_full
 
         ! Release memory
         deallocate(dxa_prev)
         deallocate(dxa_full)
      end if

      ! Interpolate guess at current resolution
      call geom(no)%interp_gp(geom(io),xg_full,xg)

      ! Right-hand side
      if (use_binv) then
         ! Allocation
         allocate(dxb(geom(io)%nh))

         ! Compute background increment at current resolution
         dxb = xg-xb

         ! Apply B inverse
         call bmatrix(io)%apply_inv(geom(io),dxb,dxbbar)

         ! Release memory
         deallocate(dxb)
      else
         ! Allocation
         allocate(dxabar_interp(geom(io)%nh))

         ! Initialization
         dxbbar = 0.0

         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_gp(geom(io),algo_planczosif(jo)%dxabar,dxabar_interp)

            ! Add contributions
            dxbbar = dxbbar-dxabar_interp
         end do

         ! Release memory
         deallocate(dxabar_interp)
      end if
   end if

   ! Compute innovation
   call hmatrix%apply(geom(io),xg,hxg)
   d = hmatrix%yo-hxg

   ! Copy and interpolate LMP vectors
   if (io>1) then
      select case (trim(lmp_mode))
      case ('spectral','ritz')
         ! Copy eigenpairs
         lmp_planczosif(io)%outer(io)%eigenval = algo_planczosif(io-1)%eigenval
         lmp_planczosif(io)%outer(io)%eigenvec = algo_planczosif(io-1)%eigenvec

         ! Copy Lanczos vectors
         lmp_planczosif(io)%lancvec_trunc = algo_planczosif(io-1)%lancvec

         ! Compute omega
         lmp_planczosif(io)%outer(io)%omega = algo_planczosif(io-1)%eigenvec(ni,:) &
 & *algo_planczosif(io-1)%lastbeta/algo_planczosif(io-1)%eigenval

         do jo=2,io
            ! Copy eigenpairs and omega
            lmp_planczosif(io)%outer(jo)%eigenval = lmp_planczosif(jo)%outer(jo)%eigenval
            lmp_planczosif(io)%outer(jo)%eigenvec = lmp_planczosif(jo)%outer(jo)%eigenvec
            lmp_planczosif(io)%outer(jo)%omega = lmp_planczosif(jo)%outer(jo)%omega

            do ii=1,ni+1
               ! Interpolation of Lanczos vectors
               call geom(jo-1)%interp_gp(geom(io),lmp_planczosif(jo)%lancvec_trunc(:,ii), &
 & lmp_planczosif(io)%outer(jo)%lancvec(:,ii))

               ! Regeneration of other Lanczos vectors
               call lmp_planczosif(io)%apply(geom(io),ni,jo-1,lmp_planczosif(io)%outer(jo)%lancvec(:,ii), &
 & lmp_planczosif(io)%outer(jo)%lancvec2(:,ii))
               call bmatrix(io)%apply(geom(io),lmp_planczosif(io)%outer(jo)%lancvec2(:,ii), &
 & lmp_planczosif(io)%outer(jo)%lancvec1(:,ii))

               if (test_ortho) then
                  ! Orthogonality test
                  do ji=1,ii-1
                     proj = sum(lmp_planczosif(io)%outer(jo)%lancvec1(:,ji)*lmp_planczosif(io)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp_planczosif(io)%outer(jo)%lancvec1(:,ji)*lmp_planczosif(io)%outer(jo)%lancvec(:,ji))
                     if (abs(proj)>1.0e-6) then
                        write(*,'(a)') 'ERROR: orthogonality lost'
                        stop
                     end if
                  end do
                  norm = sqrt(sum(lmp_planczosif(io)%outer(jo)%lancvec(:,ii)*lmp_planczosif(io)%outer(jo)%lancvec1(:,ii)))
                  if (abs(norm-1.0)>1.0e-6) then
                     write(*,'(a)') 'ERROR: orthogonality lost'
                     stop
                  end if
               end if
            end do

            ! Ritz vectors
            lmp_planczosif(io)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec(:,1:ni), &
 & lmp_planczosif(io)%outer(jo)%eigenvec)
            lmp_planczosif(io)%outer(jo)%ritzvec1(:,1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec1(:,1:ni) &
 & ,lmp_planczosif(io)%outer(jo)%eigenvec)
            lmp_planczosif(io)%outer(jo)%ritzvec2(:,1:ni) = matmul(lmp_planczosif(io)%outer(jo)%lancvec2(:,1:ni), &
 & lmp_planczosif(io)%outer(jo)%eigenvec)
         end do
      end select
   end if

   ! Minimization
   write(*,'(a)') '      Minimization'
   call algo_planczosif(io)%apply_planczosif(geom(io),bmatrix(io),hmatrix,rmatrix,dxbbar,d,ni,lmp_planczosif(io), &
 & shutoff_type,shutoff_value)

   ! Compute nonlinear cost function
   call algo_planczosif(io)%compute_cost(geom(io),bmatrix(io),hmatrix,rmatrix,xb,xg,ni)

   ! Result
   write(*,'(a)') '      Cost function: J = Jb+Jo / J_nl = Jb_nl+Jo_nl'
   do ii=0,ni
      write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '         Iteration ',ii,':', &
 & algo_planczosif(io)%j(ii),' = ',algo_planczosif(io)%jb(ii),' + ',algo_planczosif(io)%jo(ii),' / ', &
 & algo_planczosif(io)%j(ii),' = ',algo_planczosif(io)%jb_nl(ii),' + ',algo_planczosif(io)%jo_nl(ii)
   end do

   ! Create subgroup for this algo
   call ncerr('main',nf90_def_grp(grpid(io),'planczosif',subgrpid))

   ! Reshape guess
   xb_2d = reshape(xb,(/geom(io)%nx,geom(io)%ny/))

   ! Get dimensions
   call ncerr('main',nf90_inq_dimid(grpid(io),'nx',nx_id))
   call ncerr('main',nf90_inq_dimid(grpid(io),'ny',ny_id))
   call ncerr('main',nf90_inq_dimid(ncid,'nobs',nobs_id))

   ! Create variables
   call ncerr('main',nf90_def_var(subgrpid,'xg',nf90_double,(/nx_id,ny_id/),xg_id))
   call ncerr('main',nf90_def_var(subgrpid,'hxg',nf90_double,(/nobs_id/),hxg_id))
   call ncerr('main',nf90_def_var(subgrpid,'d',nf90_double,(/nobs_id/),d_id))

   ! Write variable
   call ncerr('main',nf90_put_var(subgrpid,xg_id,xg_2d))
   call ncerr('main',nf90_put_var(subgrpid,hxg_id,hxg))
   call ncerr('main',nf90_put_var(subgrpid,d_id,d))

   ! Write algo data
   call algo_planczosif(io)%write(geom(io),grpid(io),subgrpid)

   ! Release memory
   deallocate(xb)
   deallocate(xg)
   deallocate(xg_2d)
   deallocate(dxbbar)
end do
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Lanczos-PLanczosIF comparison
!--------------------------------------------------------------------------------
write(*,'(a)') 'Lanczos-PLanczosIF comparison'
do io=1,no
   write(*,'(a,i1,a,i4,a,i4)') '   Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny
   do ii=0,ni
      write(*,'(a,i3,a,e15.8,a,e15.8,a,e15.8)') '      Iteration ',ii,':', algo_lanczos(io)%j(ii)-algo_planczosif(io)%j(ii), &
 & ' = ',algo_lanczos(io)%jb(ii)-algo_planczosif(io)%jb(ii),' + ',algo_lanczos(io)%jo(ii)-algo_planczosif(io)%jo(ii)
   end do
end do

! Close NetCDF file
call ncerr('main',nf90_close(ncid))

end program main
