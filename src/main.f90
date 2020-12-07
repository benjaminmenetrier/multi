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
integer,parameter   :: nm = 3              ! Number of methods (theoretical, standard, alternative)
integer,parameter   :: na = 2              ! Number of algorithms (lanczos, planczosif)
integer,parameter   :: nomax = 9           ! Maximum number of outer iterations
integer,parameter   :: namelist_unit = 9   ! Namelist unit

! Namelist parameters
character(len=1024) :: namelist_name       ! name of the namelist file to use (default="namelist")

integer             :: no                  ! Number of outer iterations
integer             :: ni                  ! Number of inner iterations
character(len=8)    :: lmp_mode            ! LMP mode ('none', 'spectral', 'ritz')
logical             :: test_ortho          ! Test orthogonality
integer             :: shutoff_type        ! Stopping criterion according: 1-Jb, 2-beta, 3-Ritz convergence (else: no criterion)
real(8)             :: shutoff_value       ! Stopping criterion threshold
logical             :: transitive_interp   ! Transitive interpolator (grid-point interpolator defined with FFTs and spetral interpolator)
logical             :: projective_Bmatrix  ! Projective B matrix (low-resolution B matrix is a a projection of the high-resolution B matrix)
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
integer                        :: ncid,subgrpid,nx_id,ny_id,xb_id,xg_id,hxg_id,d_id
integer                        :: x_obs_id,y_obs_id,nobs_id,obs_val_id
integer                        :: io,jo,ii,ji,im,jm,ia,ja
integer,allocatable            :: grpid(:,:)
real(8)                        :: proj,norm
real(8)                        :: diff(2,2,nm,nm)
real(8),allocatable            :: xb_full(:),xg_full(:),dxb_full(:),dxa_full(:)
real(8),allocatable            :: xb(:),xg(:),dxb(:),dxbbar(:),dvb(:),dva(:),dxa(:),dxabar(:),dxa_prev(:)
real(8),allocatable            :: xb_2d(:,:),xg_2d(:,:)
real(8),allocatable            :: d(:),hxg(:)
character(len=1024)            :: method(3)
character(len=1024)            :: grpname
type(algo_type),allocatable    :: algo(:,:,:)
type(bmatrix_type),allocatable :: bmatrix(:)
type(geom_type),allocatable    :: geom(:)
type(hmatrix_type)             :: hmatrix
type(lmp_type),allocatable     :: lmp_lanczos(:,:),lmp_planczosif(:,:)
type(rmatrix_type)             :: rmatrix

! Namelist blocks
namelist/solver/ &
 & no, &
 & ni, &
 & lmp_mode, &
 & test_ortho, &
 & shutoff_type, &
 & shutoff_value, &
 & transitive_interp, &
 & projective_Bmatrix
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
transitive_interp = .true.
nx = 101
ny = 101
nobs = 2000
sigma_obs = 0.1
sigmabvar = 0.0
Lb = 0.12
spvarmin = 1.0e-5
new_seed = .false.
filename = 'output.nc'

! Read the name of the namelist to use
call get_command_argument(1,namelist_name)
write(*,'(a)') 'Reading namelist: '//trim(namelist_name)

! Read namelist
!open(unit=namelist_unit,file='namelist',status='old',action='read')
open(unit=namelist_unit,file=namelist_name,status='old',action='read')
read(namelist_unit,nml=solver)
read(namelist_unit,nml=resolutions)
read(namelist_unit,nml=observations)
read(namelist_unit,nml=background)
read(namelist_unit,nml=miscellanous)
close(unit=namelist_unit)

! Print namelist
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
   
! Methods
method(1) = 'theoretical'
method(2) = 'standard'
method(3) = 'alternative'

! Allocation (number of outer iterations)
allocate(grpid(0:no,nm))
allocate(geom(no))
allocate(bmatrix(no))
allocate(algo(no,2,nm))
allocate(lmp_lanczos(no,nm))
allocate(lmp_planczosif(no,nm))

! Create NetCDF file
call ncerr('main',nf90_create(trim(filename),ior(nf90_clobber,nf90_netcdf4),ncid))

! Set missing value
call ncerr('main',nf90_put_att(ncid,nf90_global,'_FillValue',-999.0))

! Create groups for each method and each outer iteration
do im=1,nm
   call ncerr('main',nf90_def_grp(ncid,method(im),grpid(0,im)))
   do io=1,no
      write(grpname,'(a,i1)') 'outer_',io
      call ncerr('main',nf90_def_grp(grpid(0,im),grpname,grpid(io,im)))
   end do
end do

! Setup geometries
do io=1,no
   write(*,'(a,i1)') '   Geometry setup for outer iteration ',io
   call geom(io)%setup(nx(io),ny(io),transitive_interp)
   do im=1,nm
      call geom(io)%write(grpid(io,im))
   end do
end do

! Allocation (model/control space)
allocate(xb_full(geom(no)%nh))
allocate(xg_full(geom(no)%nh))

! Setup B matrices
do io=1,no
   write(*,'(a,i1)') '   B matrix setup for outer iteration ',io
   if (projective_Bmatrix) then
      call bmatrix(io)%setup(geom(io),geom(no),sigmabvar,Lb,spvarmin)
   else
      call bmatrix(io)%setup(geom(io),geom(io),sigmabvar,Lb,spvarmin)
   end if
   do im=1,nm
      call bmatrix(io)%write(geom(io),grpid(io,im))
   end do
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

   do im=1,nm
      ! Get dimensions
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'nx',nx_id))
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'ny',ny_id))

      ! Create variable
      call ncerr('main',nf90_def_var(grpid(io,im),'xb',nf90_double,(/nx_id,ny_id/),xb_id))

      ! Write variable
      call ncerr('main',nf90_put_var(grpid(io,im),xb_id,xb_2d))
   end do

   ! Release memory
   deallocate(xb)
   deallocate(xb_2d)
end do

! Allocation (obs space)
allocate(d(nobs))
allocate(hxg(nobs))

! Setup obserations locations
call hmatrix%setup(nobs)
call hmatrix%write(ncid,x_obs_id,y_obs_id,nobs_id)

! Setup R matrix
call rmatrix%setup(nobs,sigma_obs)

! Setup observations
call rmatrix%randomize(hmatrix%yo)
! Write observations:
! Create obs variable
call ncerr('main',nf90_def_var(ncid,'obs_val',nf90_double,(/nobs_id/),obs_val_id))
! Write obs variable
call ncerr('main',nf90_put_var(ncid,obs_val_id,hmatrix%yo))


write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Methods (theoretical, standard and alternative)
!--------------------------------------------------------------------------------
do im=1,nm
   write(*,'(a,i1,a,i4,a,i4)') '   Method: '//trim(method(im))

   !--------------------------------------------------------------------------------
   ! Multi-incremental Lanczos in control space
   !--------------------------------------------------------------------------------
   write(*,'(a)') 'Multi-incremental Lanczos in control space'
   do io=1,no
      write(*,'(a,i1,a,i4,a,i4)') '      Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny

      ! Allocation
      allocate(xb(geom(io)%nh))
      allocate(xg(geom(io)%nh))
      allocate(xg_2d(geom(io)%nx,geom(io)%ny))
      allocate(dvb(geom(io)%nh))
      call algo(io,1,im)%alloc(geom(io),ni)
      call lmp_lanczos(io,im)%alloc(no,geom,ni,io,lmp_mode,'control')

      ! Test H matrix
      write(*,'(a)') '      Test H matrix'
      call hmatrix%test(geom(io))

      ! Interpolate background at current resolution
      call geom(no)%interp_gp(geom(io),xb_full,xb)

      ! Compte guess and first term of the right-hand side
      if (im==1) then
         ! Theoretical method

         ! Compute guess at full resolution
         if (io==1) then
            ! Background at full resolution
            xg_full = xb_full
         else
            ! Allocation
            allocate(dxa_prev(geom(io-1)%nh))
            allocate(dxa_full(geom(no)%nh))

            ! Compute analysis increment of the previous outer iteration at full resolution
            call bmatrix(io-1)%apply_sqrt(geom(io-1),algo(io-1,1,im)%dva,dxa_prev)
            call geom(io-1)%interp_gp(geom(no),dxa_prev,dxa_full)

            ! Add analysis increment of the previous outer iteration at full resolution
            xg_full = xg_full+dxa_full

            ! Release memory
            deallocate(dxa_prev)
            deallocate(dxa_full)
         end if

         ! Interpolate guess at current resolution
         call geom(no)%interp_gp(geom(io),xg_full,xg)

         ! Allocation
         allocate(dxb(geom(io)%nh))
         allocate(dxbbar(geom(io)%nh))

         ! Compute background increment at current resolution
         dxb = xb-xg

         ! Apply B inverse
         call bmatrix(io)%apply_inv(geom(io),dxb,dxbbar)

         ! Apply B square-root adjoint
         call bmatrix(io)%apply_sqrt_ad(geom(io),dxbbar,dvb)

         ! Release memory
         deallocate(dxb)
         deallocate(dxbbar)
      elseif (im==2) then
         ! Standard method

         ! Compute guess at full resolution
         if (io==1) then
            ! Background at full resolution
            xg_full = xb_full
         else
            ! Allocation
            allocate(dxa_full(geom(no)%nh))
            allocate(dva(geom(io)%nh))
            allocate(dxa(geom(io)%nh))

            ! Compute analysis increment of the previous outer iteration at full resolution
            call geom(io-1)%interp_sp(geom(io),algo(io-1,1,im)%dva,dva)
            call bmatrix(io)%apply_sqrt(geom(io),dva,dxa)
            call geom(io)%interp_gp(geom(no),dxa,dxa_full)

            ! Add analysis increment of the previous outer iteration at full resolution
            xg_full = xg_full+dxa_full
 
            ! Release memory
            deallocate(dva)
            deallocate(dxa)
            deallocate(dxa_full)
         end if

         ! Interpolate guess at current resolution
         call geom(no)%interp_gp(geom(io),xg_full,xg)

         ! Allocation
         allocate(dva(geom(io)%nh))

         ! Initialization
         dvb = 0.0

         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_sp(geom(io),algo(jo,1,im)%dva,dva)

             ! Add contributions
            dvb = dvb-dva
         end do

         ! Release memory
         deallocate(dva)
      elseif (im==3) then
         ! Alternative method

         ! Allocation
         allocate(dva(geom(io)%nh))
         allocate(dxb(geom(io)%nh))
         allocate(dxb_full(geom(no)%nh))

         ! Initialization
         dvb = 0.0

         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_sp(geom(io),algo(jo,1,im)%dva,dva)

            ! Add contributions
            dvb = dvb-dva
         end do

         ! Compute background increment at current resolution
         call bmatrix(io)%apply_sqrt(geom(io),dvb,dxb)

         ! Compute guess at current resolution
         xg = xb-dxb

         ! Interpolate background increment at full resolution
         call geom(io)%interp_gp(geom(no),dxb,dxb_full)

         ! Compute guess at full resolution
         xg_full = xb_full-dxb_full

         ! Release memory
         deallocate(dva)
         deallocate(dxb)
         deallocate(dxb_full)
      end if

      ! Compute innovation
      call hmatrix%apply(geom(io),xg,hxg)
      d = hmatrix%yo-hxg
   
      ! Copy and interpolate LMP vectors
      if (io>1) then
         select case (trim(lmp_mode))
         case ('spectral','ritz')
            ! Copy eigenpairs
            lmp_lanczos(io,im)%outer(io)%eigenval = algo(io-1,1,im)%eigenval
            lmp_lanczos(io,im)%outer(io)%eigenvec = algo(io-1,1,im)%eigenvec
   
            ! Copy Lanczos vectors
            lmp_lanczos(io,im)%lancvec_trunc = algo(io-1,1,im)%lancvec
   
            ! Compute omega
            lmp_lanczos(io,im)%outer(io)%omega = algo(io-1,1,im)%eigenvec(ni,:) &
 & *algo(io-1,1,im)%lastbeta/algo(io-1,1,im)%eigenval
   
            do jo=2,io
               ! Copy eigenpairs and omega
               lmp_lanczos(io,im)%outer(jo)%eigenval = lmp_lanczos(jo,im)%outer(jo)%eigenval
               lmp_lanczos(io,im)%outer(jo)%eigenvec = lmp_lanczos(jo,im)%outer(jo)%eigenvec
               lmp_lanczos(io,im)%outer(jo)%omega = lmp_lanczos(jo,im)%outer(jo)%omega
   
               do ii=1,ni+1
                  ! Interpolation of Lanczos vectors
                  call geom(jo-1)%interp_sp(geom(io),lmp_lanczos(jo,im)%lancvec_trunc(:,ii), &
 & lmp_lanczos(io,im)%outer(jo)%lancvec(:,ii))
   
                  if (test_ortho) then
                     ! Orthogonality test
                     do ji=1,ii-1
                        proj = sum(lmp_lanczos(io,im)%outer(jo)%lancvec(:,ji)*lmp_lanczos(io,im)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp_lanczos(io,im)%outer(jo)%lancvec(:,ji)*lmp_lanczos(io,im)%outer(jo)%lancvec(:,ji))
                        if (abs(proj)>1.0e-6) then
                           write(*,'(a)') 'ERROR: orthogonality lost'
                           stop
                        end if
                     end do
                     norm = sqrt(sum(lmp_lanczos(io,im)%outer(jo)%lancvec(:,ii)**2))
                     if (abs(norm-1.0)>1.0e-6) then
                        write(*,'(a)') 'ERROR: orthogonality lost '
                        stop
                     end if
                  end if
               end do
   
               ! Ritz vectors
               lmp_lanczos(io,im)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp_lanczos(io,im)%outer(jo)%lancvec(:,1:ni), &
 & lmp_lanczos(io,im)%outer(jo)%eigenvec)
            end do
         end select
      end if
   
      ! Minimization
      write(*,'(a)') '      Minimization'
      call algo(io,1,im)%apply_lanczos(geom(io),bmatrix(io),hmatrix,rmatrix,dvb,d,ni,lmp_lanczos(io,im), &
 & shutoff_type,shutoff_value)
   
      ! Compute nonlinear cost function
      call algo(io,1,im)%compute_cost(geom(io),bmatrix(io),hmatrix,rmatrix,xb,xg,ni)
   
      ! Result
      write(*,'(a)') '      Cost function: J = Jb+Jo / J_nl = Jb_nl+Jo_nl'
      do ii=0,ni
         write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '         Iteration ',ii,':', &
 & algo(io,1,im)%j(ii),' = ',algo(io,1,im)%jb(ii),' + ',algo(io,1,im)%jo(ii),' / ', &
 & algo(io,1,im)%j_nl(ii),' = ',algo(io,1,im)%jb_nl(ii),' + ',algo(io,1,im)%jo_nl(ii)
      end do
   
      ! Create subgroup for this algo
      call ncerr('main',nf90_def_grp(grpid(io,im),'lanczos',subgrpid))
   
      ! Reshape guess
      xg_2d = reshape(xg,(/geom(io)%nx,geom(io)%ny/))
   
      ! Get dimensions
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'nx',nx_id))
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'ny',ny_id))
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
      call algo(io,1,im)%write(geom(io),grpid(io,im),subgrpid)
   
      ! Release memory
      deallocate(xb)
      deallocate(xg)
      deallocate(xg_2d)
      deallocate(dvb)
   end do
   write(*,'(a)') ''
end do
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Methods (theoretical, standard and alternative)
!--------------------------------------------------------------------------------
do im=1,nm
   write(*,'(a,i1,a,i4,a,i4)') '      Method: '//trim(method(im))
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
      call algo(io,2,im)%alloc(geom(io),ni)
      call lmp_planczosif(io,im)%alloc(no,geom,ni,io,lmp_mode,'model')
   
      ! Test H matrix
      write(*,'(a)') '      Test H matrix'
      call hmatrix%test(geom(io))
   
      ! Interpolate background at current resolution
      call geom(no)%interp_gp(geom(io),xb_full,xb)
   
      ! Compte guess and first term of the right-hand side
      if (im==1) then
         ! Theoretical method

         ! Compute guess at full resolution
         if (io==1) then
            ! Background at full resolution
            xg_full = xb_full
         else
            ! Allocation
            allocate(dxa_prev(geom(io-1)%nh))
            allocate(dxa_full(geom(no)%nh))
   
            ! Compute analysis increment of the previous outer iteration at full resolution
            call bmatrix(io-1)%apply(geom(io-1),algo(io-1,2,im)%dxabar,dxa_prev)
            call geom(io-1)%interp_gp(geom(no),dxa_prev,dxa_full)
   
            ! Add analysis increment of the previous outer iteration at full resolution
            xg_full = xg_full+dxa_full
    
            ! Release memory
            deallocate(dxa_prev)
            deallocate(dxa_full)
         end if
   
         ! Interpolate guess at current resolution
         call geom(no)%interp_gp(geom(io),xg_full,xg)
   
         ! Allocation
         allocate(dxb(geom(io)%nh))
   
         ! Compute background increment at current resolution
         dxb = xb-xg
   
         ! Apply B inverse
         call bmatrix(io)%apply_inv(geom(io),dxb,dxbbar)
   
         ! Release memory
         deallocate(dxb)
      elseif (im==2) then
         ! Standard method

         ! Compute guess at full resolution
         if (io==1) then
            ! Background at full resolution
            xg_full = xb_full
         else
            ! Allocation
            allocate(dxa_full(geom(no)%nh))
            allocate(dxabar(geom(io)%nh))
            allocate(dxa(geom(io)%nh))
   
            ! Compute analysis increment of the previous outer iteration at full resolution
            call geom(io-1)%interp_gp(geom(io),algo(io-1,2,im)%dxabar,dxabar)
            call bmatrix(io)%apply(geom(io),dxabar,dxa)
            call geom(io)%interp_gp(geom(no),dxa,dxa_full)
   
            ! Add analysis increment of the previous outer iteration at full resolution
            xg_full = xg_full+dxa_full
    
            ! Release memory
            deallocate(dxabar)
            deallocate(dxa)
            deallocate(dxa_full)
         end if
   
         ! Interpolate guess at current resolution
         call geom(no)%interp_gp(geom(io),xg_full,xg)
   
         ! Allocation
         allocate(dxabar(geom(io)%nh))
   
         ! Initialization
         dxbbar = 0.0
   
         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_gp(geom(io),algo(jo,2,im)%dxabar,dxabar)
   
            ! Add contributions
            dxbbar = dxbbar-dxabar
         end do
   
         ! Release memory
         deallocate(dxabar)
      elseif (im==3) then
         ! Alternative method

         ! Allocation
         allocate(dxabar(geom(io)%nh))
         allocate(dxb(geom(io)%nh))
         allocate(dxb_full(geom(no)%nh))
   
         ! Initialization
         dxbbar = 0.0
   
         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_gp(geom(io),algo(jo,2,im)%dxabar,dxabar)
   
            ! Add contributions
            dxbbar = dxbbar-dxabar
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
         deallocate(dxabar)
         deallocate(dxb)
         deallocate(dxb_full)
      end if
   
      ! Compute innovation
      call hmatrix%apply(geom(io),xg,hxg)
      d = hmatrix%yo-hxg
   
      ! Copy and interpolate LMP vectors
      if (io>1) then
         select case (trim(lmp_mode))
         case ('spectral','ritz')
            ! Copy eigenpairs
            lmp_planczosif(io,im)%outer(io)%eigenval = algo(io-1,2,im)%eigenval
            lmp_planczosif(io,im)%outer(io)%eigenvec = algo(io-1,2,im)%eigenvec
   
            ! Copy Lanczos vectors
            lmp_planczosif(io,im)%lancvec_trunc = algo(io-1,2,im)%lancvec
   
            ! Compute omega
            lmp_planczosif(io,im)%outer(io)%omega = algo(io-1,2,im)%eigenvec(ni,:) &
 & *algo(io-1,2,im)%lastbeta/algo(io-1,2,im)%eigenval
   
            do jo=2,io
               ! Copy eigenpairs and omega
               lmp_planczosif(io,im)%outer(jo)%eigenval = lmp_planczosif(jo,im)%outer(jo)%eigenval
               lmp_planczosif(io,im)%outer(jo)%eigenvec = lmp_planczosif(jo,im)%outer(jo)%eigenvec
               lmp_planczosif(io,im)%outer(jo)%omega = lmp_planczosif(jo,im)%outer(jo)%omega
   
               do ii=1,ni+1
                  ! Interpolation of Lanczos vectors
                  call geom(jo-1)%interp_gp(geom(io),lmp_planczosif(jo,im)%lancvec_trunc(:,ii), &
 & lmp_planczosif(io,im)%outer(jo)%lancvec(:,ii))
   
                  ! Regeneration of other Lanczos vectors
                  call lmp_planczosif(io,im)%apply(geom(io),ni,jo-1,lmp_planczosif(io,im)%outer(jo)%lancvec(:,ii), &
 & lmp_planczosif(io,im)%outer(jo)%lancvec2(:,ii))
                  call bmatrix(io)%apply(geom(io),lmp_planczosif(io,im)%outer(jo)%lancvec2(:,ii), &
 & lmp_planczosif(io,im)%outer(jo)%lancvec1(:,ii))
   
                  if (test_ortho) then
                     ! Orthogonality test
                     do ji=1,ii-1
                        proj = sum(lmp_planczosif(io,im)%outer(jo)%lancvec1(:,ji)*lmp_planczosif(io,im)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp_planczosif(io,im)%outer(jo)%lancvec1(:,ji)*lmp_planczosif(io,im)%outer(jo)%lancvec(:,ji))
                        if (abs(proj)>1.0e-6) then
                           write(*,'(a)') 'ERROR: orthogonality lost'
                           stop
                        end if
                     end do
                     norm = sqrt(sum(lmp_planczosif(io,im)%outer(jo)%lancvec(:,ii)*lmp_planczosif(io,im)%outer(jo)%lancvec1(:,ii)))
                     if (abs(norm-1.0)>1.0e-6) then
                        write(*,'(a)') 'ERROR: orthogonality lost'
                        stop
                     end if
                  end if
               end do
   
               ! Ritz vectors
               lmp_planczosif(io,im)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp_planczosif(io,im)%outer(jo)%lancvec(:,1:ni), &
 & lmp_planczosif(io,im)%outer(jo)%eigenvec)
               lmp_planczosif(io,im)%outer(jo)%ritzvec1(:,1:ni) = matmul(lmp_planczosif(io,im)%outer(jo)%lancvec1(:,1:ni) &
 & ,lmp_planczosif(io,im)%outer(jo)%eigenvec)
               lmp_planczosif(io,im)%outer(jo)%ritzvec2(:,1:ni) = matmul(lmp_planczosif(io,im)%outer(jo)%lancvec2(:,1:ni), &
 & lmp_planczosif(io,im)%outer(jo)%eigenvec)
            end do
         end select
      end if
   
      ! Minimization
      write(*,'(a)') '      Minimization'
      call algo(io,2,im)%apply_planczosif(geom(io),bmatrix(io),hmatrix,rmatrix,dxbbar,d,ni,lmp_planczosif(io,im), &
 & shutoff_type,shutoff_value)
   
      ! Compute nonlinear cost function
      call algo(io,2,im)%compute_cost(geom(io),bmatrix(io),hmatrix,rmatrix,xb,xg,ni)
   
      ! Result
      write(*,'(a)') '      Cost function: J = Jb+Jo / J_nl = Jb_nl+Jo_nl'
      do ii=0,ni
         write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '         Iteration ',ii,':', &
 & algo(io,2,im)%j(ii),' = ',algo(io,2,im)%jb(ii),' + ',algo(io,2,im)%jo(ii),' / ', &
 & algo(io,2,im)%j(ii),' = ',algo(io,2,im)%jb_nl(ii),' + ',algo(io,2,im)%jo_nl(ii)
      end do

      ! Create subgroup for this algo
      call ncerr('main',nf90_def_grp(grpid(io,im),'planczosif',subgrpid))
   
      ! Reshape guess
      xg_2d = reshape(xg,(/geom(io)%nx,geom(io)%ny/))
   
      ! Get dimensions
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'nx',nx_id))
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'ny',ny_id))
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
      call algo(io,2,im)%write(geom(io),grpid(io,im),subgrpid)
   
      ! Release memory
      deallocate(xb)
      deallocate(xg)
      deallocate(xg_2d)
      deallocate(dxbbar)
   end do
   write(*,'(a)') ''
end do
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! General comparison
!--------------------------------------------------------------------------------
write(*,'(a)') 'General comparison'
diff = 0.0
do im=1,nm
   do jm=1,nm
      do ia=1,na
         do ja=1,na
            ! Compute max relative difference between methods/preconditionings
            do io=1,no
               do ii=0,ni
                  diff(ia,ja,im,jm) = max(diff(ia,ja,im,jm), &
 & 2.0*abs(algo(io,ia,im)%j(ii)-algo(io,ja,jm)%j(ii))/abs(algo(io,ia,im)%j(ii)+algo(io,ja,jm)%j(ii)))
               end do
           end do
           write(*,'(a,i1,a,i1,a,i1,a,i1,a,e11.4)') 'Methods ',im,'/',jm,' :: algos ',ia,'/',ja,' => ',diff(ia,ja,im,jm)
        end do
      end do
   end do
end do

! Close NetCDF file
call ncerr('main',nf90_close(ncid))

end program main
