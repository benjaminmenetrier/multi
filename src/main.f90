program main

use netcdf
use tools_kinds
use tools_netcdf
use tools_rand
use type_algo
use type_bmatrix
use type_geom
use type_hoperator
use type_lmp
use type_rmatrix

implicit none

! Parameters
integer,parameter :: nmmax = 3           ! Maximum number of methods (theoretical, standard, alternative)
integer,parameter :: namax = 2           ! Maximum number of algorithms (lanczos, planczosif)
integer,parameter :: nomax = 99          ! Maximum number of outer iterations
integer,parameter :: namelist_unit = 9   ! Namelist unit
real(kind_real),parameter :: threshold = 1.0e-10_kind_real ! Threshold for similar results

! Namelist parameters
character(len=1024) :: namelist_name    ! Namelist file name (default value is 'namelist')
integer :: nm                           ! Number of methods
character(len=1024) :: method(nmmax)    ! Solver methods
integer :: na                           ! Number of algorithms
character(len=1024) :: algorithm(namax) ! Solver algorithms
integer :: no                           ! Number of outer iterations
integer :: ni                           ! Number of inner iterations
character(len=1024) :: lmp_mode         ! LMP mode ('none', 'spectral', 'ritz')
logical :: test_ortho                   ! Test orthogonality
integer :: shutoff_type                 ! Stopping criterion according: 1-Jb, 2-beta, 3-Ritz convergence (else: no criterion)
real(kind_real) :: shutoff_value        ! Stopping criterion threshold
character(len=1024) :: interp_method    ! Interpolation method ('spectral',  or not)
logical :: projective_Bmatrix           ! Projective B matrix (low-resolution B matrix is a a projection of the high-resolution B matrix)
integer :: nx(nomax)                    ! X direction sizes
integer :: ny(nomax)                    ! Y direction sizes
integer :: nobs                         ! Number of observations
real(kind_real) :: sigma_obs            ! Observation error standard deviation
real(kind_real) :: Hnl_coeff            ! Noninear coefficient for H operator
real(kind_real) :: sigmabvar            ! Grid-point standard deviation variations amplitude
real(kind_real) :: Lb                   ! Correlation length-scale
real(kind_real) :: spvarmin             ! Minimum spectral variance (inverse fails if this is too small)
logical :: new_seed                     ! New random seed
character(len=1024) :: filename         ! Output file name

! Local variables
integer :: ncid,subgrpid,nx_id,ny_id,xb_id,xg_id,hxg_id,d_id,iobs,i
integer :: x_obs_id,y_obs_id,nobs_id,yt_id,yo_id,nx_full_id,ny_full_id,x_full_id,y_full_id,xt_id
integer :: io,jo,ii,ji,im,jm,ia,ja
integer,allocatable :: grpid(:,:)
real(kind_real),allocatable :: x_obs(:),y_obs(:),obs_err(:)
real(kind_real) :: proj,norm
real(kind_real),allocatable :: maxdiff(:,:,:,:)
real(kind_real),allocatable :: xb_full(:),xg_full(:),dxb_full(:),dxa_full(:,:),xt(:),yt(:)
real(kind_real),allocatable :: xb_full_io(:,:),xg_full_io(:,:)
real(kind_real),allocatable :: xb(:),dxb(:),dxbbar(:),dva(:),dxa(:),dxa_tmp(:),dxabar(:),dxa_prev(:)
real(kind_real),allocatable :: xb_2d(:,:),xg_2d(:,:),xt_2d(:,:)
real(kind_real),allocatable :: hxg(:)
character(len=1024) :: grpname
type(algo_type),allocatable :: algo(:,:,:)
type(bmatrix_type),allocatable :: bmatrix(:)
type(geom_type),allocatable :: geom(:)
type(hoperator_type) :: hoperator
type(lmp_type),allocatable :: lmp(:,:,:)
type(rmatrix_type) :: rmatrix

! Namelist blocks
namelist/solver/ &
 & nm, &
 & method, &
 & na, &
 & algorithm, &
 & no, &
 & ni, &
 & lmp_mode, &
 & test_ortho, &
 & shutoff_type, &
 & shutoff_value, &
 & interp_method, &
 & projective_Bmatrix
namelist/resolutions/ &
 & nx, &
 & ny
namelist/observations/ &
 & nobs, &
 & sigma_obs, &
 & Hnl_coeff
namelist/background/ &
 & sigmabvar, &
 & Lb, &
 & spvarmin
namelist/miscellanous/ &
 & new_seed, &
 & filename

! Initialize namelist (default values)
nm = 3
method(1) = 'theoretical'
method(2) = 'standard'
method(3) = 'alternative'
na = 2
algorithm(1) = 'lanczos'
algorithm(2) = 'planczosif'
no = 1
ni = 10
lmp_mode = 'none'
test_ortho = .false.
shutoff_type = 0
shutoff_value = zero
interp_method = 'spectral'
nx = 101
ny = 101
nobs = 2000
sigma_obs = 0.1
Hnl_coeff = 0.
sigmabvar = zero
Lb = 0.12
spvarmin = 1.0e-5_kind_real
new_seed = .false.
filename = 'output.nc'

! Read the name of the namelist to use
if (command_argument_count()==0) then
   namelist_name = 'namelist'
else
   call get_command_argument(1,namelist_name)
end if
write(*,'(a)') 'Reading namelist: '//trim(namelist_name)

! Read namelist
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

! Allocation (number of outer iterations)
allocate(grpid(0:no,nm))
allocate(geom(no))
allocate(bmatrix(no))
allocate(algo(no,2,nm))
allocate(lmp(no,na,nm))

! Create NetCDF file
call ncerr('main',nf90_create(trim(filename),ior(nf90_clobber,nf90_netcdf4),ncid))

! Set missing value
call ncerr('main',nf90_put_att(ncid,nf90_global,'_FillValue',msv_real))

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
   call geom(io)%setup(nx(io),ny(io),interp_method)
   do im=1,nm
      call geom(io)%write(grpid(io,im))
   end do
end do

! Test interpolation transitivity
if ((no>1).and.(geom(1)%nh<geom(no)%nh)) call geom(1)%transitive_interp_test(geom(no))

! Setup B matrices
do io=1,no
   write(*,'(a,i1)') '   B matrix setup for outer iteration ',io
   if (projective_Bmatrix) then
      if (abs(sigmabvar)>zero) then
         write(*,'(a)') 'Error: sigmabvar should be 0 for a projective B matrix family'
         stop
      end if
      call bmatrix(io)%setup(geom(io),geom(no),sigmabvar,Lb,spvarmin)
   else
      call bmatrix(io)%setup(geom(io),geom(io),sigmabvar,Lb,spvarmin)
   end if
   do im=1,nm
      call bmatrix(io)%write(geom(io),grpid(io,im))
   end do
end do

! Allocation (model/control space)

! Background
allocate(xb_full(geom(no)%nh))

! Guess
allocate(xg_full(geom(no)%nh))

! Truth
allocate(xt(geom(no)%nh))
allocate(xt_2d(geom(no)%nx,geom(no)%ny))

! Initialize the states (xt, xb and obs):

! Truth generation
call bmatrix(no)%randomize(geom(no),xt)
xt = one+xt
xt_2d = reshape(xt,(/geom(no)%nx,geom(no)%ny/))

! Generate background state
call bmatrix(no)%randomize(geom(no),xb_full)
xb_full = xt+xb_full

! Write initial variables

! Create dimensions
call ncerr('main',nf90_def_dim(ncid,'nx_full',geom(no)%nx,nx_full_id))
call ncerr('main',nf90_def_dim(ncid,'ny_full',geom(no)%ny,ny_full_id))

! Create variables
call ncerr('main',nf90_def_var(ncid,'x_full',nc_kind_real,(/nx_full_id/),x_full_id))
call ncerr('main',nf90_def_var(ncid,'y_full',nc_kind_real,(/ny_full_id/),y_full_id))
call ncerr('main',nf90_def_var(ncid,'xt',nc_kind_real,(/nx_full_id,ny_full_id/),xt_id))

! Write variables
call ncerr('main',nf90_put_var(ncid,x_full_id,geom(no)%x))
call ncerr('main',nf90_put_var(ncid,y_full_id,geom(no)%y))
call ncerr('main',nf90_put_var(ncid,xt_id,xt_2d))

! Allocation (obs space)
allocate(yt(nobs))
allocate(hxg(nobs))

! Setup obserations locations
call hoperator%setup(nobs,Hnl_coeff)

! Setup true observations
call hoperator%apply(geom(no),xt,yt)

! Setup R matrix
call rmatrix%setup(nobs,sigma_obs)

! Setup observations perturbations
call rmatrix%randomize(hoperator%yo)

! Setup observations
hoperator%yo = yt+hoperator%yo

! Write observations

! Create dimensions
call ncerr('main',nf90_def_dim(ncid,'nobs',hoperator%nobs,nobs_id))

! Create variables
call ncerr('main',nf90_def_var(ncid,'x_obs',nc_kind_real,(/nobs_id/),x_obs_id))
call ncerr('main',nf90_def_var(ncid,'y_obs',nc_kind_real,(/nobs_id/),y_obs_id))
call ncerr('main',nf90_def_var(ncid,'yt',nc_kind_real,(/nobs_id/),yt_id))
call ncerr('main',nf90_def_var(ncid,'yo',nc_kind_real,(/nobs_id/),yo_id))

! Write variables
call ncerr('main',nf90_put_var(ncid,x_obs_id,hoperator%x_obs))
call ncerr('main',nf90_put_var(ncid,y_obs_id,hoperator%y_obs))
call ncerr('main',nf90_put_var(ncid,yt_id,yt))
call ncerr('main',nf90_put_var(ncid,yo_id,hoperator%yo))

! Release memory
deallocate(yt)

do io=1,no
   ! Allocation
   allocate(xb(geom(io)%nh))
   allocate(xb_2d(geom(io)%nx,geom(io)%ny))

   ! Interpolate background at current resolution
   call geom(no)%interp(geom(io),xb_full,xb)

   ! Reshape background
   xb_2d = reshape(xb,(/geom(io)%nx,geom(io)%ny/))

   do im=1,nm
      ! Get dimensions
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'nx',nx_id))
      call ncerr('main',nf90_inq_dimid(grpid(io,im),'ny',ny_id))

      ! Create variable
      call ncerr('main',nf90_def_var(grpid(io,im),'xb',nc_kind_real,(/nx_id,ny_id/),xb_id))

      ! Write variable
      call ncerr('main',nf90_put_var(grpid(io,im),xb_id,xb_2d))
   end do

   ! Release memory
   deallocate(xb)
   deallocate(xb_2d)
end do

write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Loop on methods
!--------------------------------------------------------------------------------
do im=1,nm
   write(*,'(a,i1,a,i4,a,i4)') 'Method: '//trim(method(im))

   !--------------------------------------------------------------------------------
   ! Loop on algorithms
   !--------------------------------------------------------------------------------
   do ia=1,na
      write(*,'(a)') '   Algorithm: '//trim(algorithm(ia))

      ! Allocation
      if ((trim(method(im))=='theoretical').or.(trim(method(im))=='standard')) allocate(dxa_full(geom(no)%nh,no))

      do io=1,no
         write(*,'(a,i1,a,i4,a,i4)') '      Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny

         ! Allocation
         allocate(xb(geom(io)%nh))
         allocate(xg_2d(geom(io)%nx,geom(io)%ny))
         call algo(io,ia,im)%alloc(geom(io),ni,nobs)
         call lmp(io,ia,im)%alloc(no,geom,ni,io,lmp_mode,algorithm(ia))

         ! Test H matrix
         write(*,'(a)') '         Test H matrix'
         call hoperator%test(geom(io))

         ! Interpolate background at current resolution
         call geom(no)%interp(geom(io),xb_full,xb)
         
         ! Compte guess and first term of the right-hand side
         if (trim(method(im))=='theoretical') then
            ! Theoretical method

            ! Compute guess at full resolution
            if (io==1) then
               ! Background at full resolution
               xg_full = xb_full
            else
               ! Allocation
               allocate(dxa_prev(geom(io-1)%nh))

               ! Compute analysis increment from the previous outer iteration at full resolution
               if (trim(algorithm(ia))=='lanczos') then
                  call bmatrix(io-1)%apply_sqrt(geom(io-1),algo(io-1,ia,im)%dva,dxa_prev)
                  call geom(io-1)%interp(geom(no),dxa_prev,dxa_full(:,io-1))
               elseif (trim(algorithm(ia))=='planczosif') then
                  call bmatrix(io-1)%apply(geom(io-1),algo(io-1,ia,im)%dxabar,dxa_prev)
                  call geom(io-1)%interp(geom(no),dxa_prev,dxa_full(:,io-1))
               end if

               ! Add analysis increment of the previous outer iteration at full resolution
               xg_full = xg_full+dxa_full(:,io-1)

               ! Release memory
               deallocate(dxa_prev)
            end if
   
            ! Interpolate guess at current resolution
            call geom(no)%interp(geom(io),xg_full,algo(io,ia,im)%xg)
            
            ! Allocation
            allocate(dxb(geom(io)%nh))

            ! Compute background increment at current resolution
            dxb = xb-algo(io,ia,im)%xg

            ! Apply B inverse
            call bmatrix(io)%apply_inv(geom(io),dxb,algo(io,ia,im)%dxbbar)

            ! Apply B square-root adjoint
            if (trim(algorithm(ia))=='lanczos') call bmatrix(io)%apply_sqrt_ad(geom(io),algo(io,ia,im)%dxbbar,algo(io,ia,im)%dvb)

            ! Release memory
            deallocate(dxb)
         elseif (trim(method(im))=='standard') then
            ! Standard method

            ! Compute guess at full resolution
            if (io==1) then
               ! Background at full resolution
               xg_full = xb_full
            else
               if (projective_Bmatrix) then
                  ! Allocation
                  allocate(dxa(geom(io-1)%nh))
   
                  ! Compute analysis increment from the previous outer iteration at full resolution
                  if (trim(algorithm(ia))=='lanczos') then
                     call bmatrix(io-1)%apply_sqrt(geom(io-1),algo(io-1,ia,im)%dva,dxa)
                  elseif (trim(algorithm(ia))=='planczosif') then
                     call bmatrix(io-1)%apply(geom(io-1),algo(io-1,ia,im)%dxabar,dxa)
                  end if
                  call geom(io-1)%interp(geom(no),dxa,dxa_full(:,io-1))
                  
                  ! Add analysis increment of the previous outer iteration at full resolution
                  xg_full = xg_full+dxa_full(:,io-1)
 
                  ! Release memory
                  deallocate(dxa)
               else
                  ! Allocation
                  allocate(dxa(geom(io)%nh))
                  allocate(dxa_tmp(geom(io)%nh))   

                  ! Initialization
                  dxa = zero

                  ! Compute analysis increment from the previous outer iterations at full resolution
                  if (trim(algorithm(ia))=='lanczos') then
                     ! Allocation
                     allocate(dva(geom(io)%nh))

                     ! Loop over previous iterations
                     do jo=1,io-1
                        call geom(jo)%interp(geom(io),algo(jo,ia,im)%dva,dva)
                        call bmatrix(io)%apply_sqrt(geom(io),dva,dxa_tmp)
                        dxa = dxa+dxa_tmp
                     end do
                     do jo=1,io-2
                        call geom(no)%interp(geom(io),dxa_full(:,jo),dxa_tmp)
                        dxa = dxa-dxa_tmp
                     end do

                     ! Release memory
                     deallocate(dva)
                  elseif (trim(algorithm(ia))=='planczosif') then
                     ! Allocation
                     allocate(dxabar(geom(io)%nh))

                     ! Loop over previous iterations
                     do jo=1,io-1
                        call geom(jo)%interp(geom(io),algo(jo,ia,im)%dxabar,dxabar)
                        call bmatrix(io)%apply(geom(io),dxabar,dxa_tmp)
                        dxa = dxa+dxa_tmp
                     end do
                     do jo=1,io-2
                        call geom(no)%interp(geom(io),dxa_full(:,jo),dxa_tmp)
                        dxa = dxa-dxa_tmp
                     end do

                     ! Release memory
                     deallocate(dxabar)
                  end if                   
                  call geom(io)%interp(geom(no),dxa,dxa_full(:,io-1))
                  
                  ! Add analysis increment of the previous outer iteration at full resolution
                  xg_full = xg_full+dxa_full(:,io-1)
 
                  ! Release memory
                  deallocate(dxa)
                  deallocate(dxa_tmp)
               end if
            end if
   
            ! Interpolate guess at current resolution
            call geom(no)%interp(geom(io),xg_full,algo(io,ia,im)%xg)
            
            if (trim(algorithm(ia))=='lanczos') then 
               ! Allocation
               allocate(dva(geom(io)%nh))

               ! Initialization
               algo(io,ia,im)%dvb = zero
   
               do jo=1,io-1
                  ! Interpolate analysis increment of previous iterations at current resolution
                  call geom(jo)%interp(geom(io),algo(jo,ia,im)%dva,dva)
   
                  ! Add contributions
                  algo(io,ia,im)%dvb = algo(io,ia,im)%dvb-dva
               end do
   
               ! Release memory
               deallocate(dva)
            elseif (trim(algorithm(ia))=='planczosif') then  
               ! Allocation
               allocate(dxabar(geom(io)%nh))
   
               ! Initialization
               algo(io,ia,im)%dxbbar = zero
   
               do jo=1,io-1
                  ! Interpolate analysis increment of previous iterations at current resolution
                  call geom(jo)%interp(geom(io),algo(jo,ia,im)%dxabar,dxabar)
                  
                  ! Add contributions
                  algo(io,ia,im)%dxbbar = algo(io,ia,im)%dxbbar-dxabar
               end do
   
               ! Release memory
               deallocate(dxabar)
            end if
         elseif (trim(method(im))=='alternative') then
            ! Alternative method
            allocate(dxb(geom(io)%nh))
            allocate(dxb_full(geom(no)%nh))
   
            if (trim(algorithm(ia))=='lanczos') then
               ! Allocation
               allocate(dva(geom(io)%nh))

               ! Initialization
               algo(io,ia,im)%dvb = zero
   
               do jo=1,io-1
                  ! Interpolate analysis increment of previous iterations at current resolution
                  call geom(jo)%interp(geom(io),algo(jo,ia,im)%dva,dva)
   
                  ! Add contributions
                  algo(io,ia,im)%dvb = algo(io,ia,im)%dvb-dva
               end do
   
               ! Compute background increment at current resolution
               call bmatrix(io)%apply_sqrt(geom(io),algo(io,ia,im)%dvb,dxb)

               ! Release memory
               deallocate(dva)
            elseif (trim(algorithm(ia))=='planczosif') then 
               ! Allocation
               allocate(dxabar(geom(io)%nh))

               ! Initialization
               algo(io,ia,im)%dxbbar = zero
   
               do jo=1,io-1
                  ! Interpolate analysis increment of previous iterations at current resolution
                  call geom(jo)%interp(geom(io),algo(jo,ia,im)%dxabar,dxabar)
                  
                  ! Add contributions
                  algo(io,ia,im)%dxbbar = algo(io,ia,im)%dxbbar-dxabar
               end do
   
               ! Compute background increment at current resolution
               call bmatrix(io)%apply(geom(io),algo(io,ia,im)%dxbbar,dxb)

               ! Release memory
               deallocate(dxabar)
            end if
   
            ! Compute guess at current resolution
            algo(io,ia,im)%xg = xb-dxb
   
            ! Interpolate background increment at full resolution
            call geom(io)%interp(geom(no),dxb,dxb_full)
            
            ! Compute guess at full resolution
            xg_full = xb_full-dxb_full
   
            ! Release memory
            deallocate(dxb)
            deallocate(dxb_full)
         end if
   
         ! Compute innovation
         call hoperator%apply(geom(io),algo(io,ia,im)%xg,hxg)
         algo(io,ia,im)%d = hoperator%yo-hxg
      
         ! Copy and interpolate LMP vectors
         if (io>1) then
            select case (trim(lmp_mode))
            case ('spectral','ritz')
               ! Copy eigenpairs
               lmp(io,ia,im)%outer(io)%eigenval = algo(io-1,ia,im)%eigenval
               lmp(io,ia,im)%outer(io)%eigenvec = algo(io-1,ia,im)%eigenvec
      
               ! Copy Lanczos vectors
               lmp(io,ia,im)%lancvec_trunc = algo(io-1,ia,im)%lancvec
      
               ! Compute omega
               lmp(io,ia,im)%outer(io)%omega = algo(io-1,ia,im)%eigenvec(ni,:)*algo(io-1,ia,im)%lastbeta/algo(io-1,ia,im)%eigenval

               do jo=2,io
                  ! Copy eigenpairs and omega
                  lmp(io,ia,im)%outer(jo)%eigenval = lmp(jo,ia,im)%outer(jo)%eigenval
                  lmp(io,ia,im)%outer(jo)%eigenvec = lmp(jo,ia,im)%outer(jo)%eigenvec
                  lmp(io,ia,im)%outer(jo)%omega = lmp(jo,ia,im)%outer(jo)%omega
      
                  do ii=1,ni+1
                     if (trim(algorithm(ia))=='lanczos') then
                        ! Interpolation of Lanczos vectors
                        call geom(jo-1)%interp(geom(io),lmp(jo,ia,im)%lancvec_trunc(:,ii),lmp(io,ia,im)%outer(jo)%lancvec(:,ii))
                     elseif (trim(algorithm(ia))=='planczosif') then
                        ! Interpolation of Lanczos vectors
                        call geom(jo-1)%interp(geom(io),lmp(jo,ia,im)%lancvec_trunc(:,ii),lmp(io,ia,im)%outer(jo)%lancvec(:,ii))
   
                        ! Regeneration of other Lanczos vectors
                        call lmp(io,ia,im)%apply(geom(io),ni,jo-1,lmp(io,ia,im)%outer(jo)%lancvec(:,ii), &
 & lmp(io,ia,im)%outer(jo)%lancvec2(:,ii))
                        call bmatrix(io)%apply(geom(io),lmp(io,ia,im)%outer(jo)%lancvec2(:,ii), &
 & lmp(io,ia,im)%outer(jo)%lancvec1(:,ii))
                     end if

                     if (test_ortho) then
                        ! Orthogonality test
                        do ji=1,ii-1
                           if (trim(algorithm(ia))=='lanczos') then
                              proj = sum(lmp(io,ia,im)%outer(jo)%lancvec(:,ji)*lmp(io,ia,im)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp(io,ia,im)%outer(jo)%lancvec(:,ji)*lmp(io,ia,im)%outer(jo)%lancvec(:,ji))
                           elseif (trim(algorithm(ia))=='planczosif') then
                              proj = sum(lmp(io,ia,im)%outer(jo)%lancvec1(:,ji)*lmp(io,ia,im)%outer(jo)%lancvec(:,ii)) &
 & /sum(lmp(io,ia,im)%outer(jo)%lancvec1(:,ji)*lmp(io,ia,im)%outer(jo)%lancvec(:,ji))
                           end if
                           if (abs(proj)>1.0e-6_kind_real) then
                              write(*,'(a)') 'Error: orthogonality lost'
                              stop
                           end if
                        end do
                        if (trim(algorithm(ia))=='lanczos') then
                           norm = sqrt(sum(lmp(io,ia,im)%outer(jo)%lancvec(:,ii)**2))
                        elseif (trim(algorithm(ia))=='planczosif') then
                           norm = sqrt(sum(lmp(io,ia,im)%outer(jo)%lancvec(:,ii)*lmp(io,ia,im)%outer(jo)%lancvec1(:,ii)))
                        end if
                        if (abs(norm-one)>1.0e-6_kind_real) then
                           write(*,'(a)') 'Error: orthogonality lost '
                           stop
                        end if
                     end if
                  end do
      
                  ! Ritz vectors
                  if (trim(algorithm(ia))=='lanczos') then
                     lmp(io,ia,im)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp(io,ia,im)%outer(jo)%lancvec(:,1:ni), &
 & lmp(io,ia,im)%outer(jo)%eigenvec)
                  elseif (trim(algorithm(ia))=='planczosif') then
                     lmp(io,ia,im)%outer(jo)%ritzvec(:,1:ni) = matmul(lmp(io,ia,im)%outer(jo)%lancvec(:,1:ni), &
 & lmp(io,ia,im)%outer(jo)%eigenvec)
                     lmp(io,ia,im)%outer(jo)%ritzvec1(:,1:ni) = matmul(lmp(io,ia,im)%outer(jo)%lancvec1(:,1:ni), &
 & lmp(io,ia,im)%outer(jo)%eigenvec)
                     lmp(io,ia,im)%outer(jo)%ritzvec2(:,1:ni) = matmul(lmp(io,ia,im)%outer(jo)%lancvec2(:,1:ni), &
 & lmp(io,ia,im)%outer(jo)%eigenvec)

                  end if
               end do
            end select
         end if
      
         ! Minimization
         write(*,'(a)') '         Minimization'
         if (trim(algorithm(ia))=='lanczos') then
            call algo(io,ia,im)%apply_lanczos(geom(io),bmatrix(io),hoperator,rmatrix,ni,lmp(io,ia,im), &
 & shutoff_type,shutoff_value)
         elseif (trim(algorithm(ia))=='planczosif') then
            call algo(io,ia,im)%apply_planczosif(geom(io),bmatrix(io),hoperator,rmatrix,ni,lmp(io,ia,im), &
 & shutoff_type,shutoff_value)
         end if

         ! Compute nonlinear cost function
         call algo(io,ia,im)%compute_cost(geom(io),bmatrix(io),hoperator,rmatrix,xb)
         
         ! Result
         write(*,'(a)') '         Cost function: J = Jb+Jo / j_full = Jb_full+Jo_full'
         do ii=0,ni
            write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '            Iteration ',ii,':', &
    & algo(io,ia,im)%j_quad(ii),' = ',algo(io,ia,im)%jb_quad(ii),' + ',algo(io,ia,im)%jo_quad(ii),' / ', &
    & algo(io,ia,im)%j_full(ii),' = ',algo(io,ia,im)%jb_full(ii),' + ',algo(io,ia,im)%jo_full(ii)
         end do
      
         ! Create subgroup for this algo
         call ncerr('main',nf90_def_grp(grpid(io,im),algorithm(ia),subgrpid))
      
         ! Reshape guess
         xg_2d = reshape(algo(io,ia,im)%xg,(/geom(io)%nx,geom(io)%ny/))
      
         ! Get dimensions
         call ncerr('main',nf90_inq_dimid(grpid(io,im),'nx',nx_id))
         call ncerr('main',nf90_inq_dimid(grpid(io,im),'ny',ny_id))
         call ncerr('main',nf90_inq_dimid(ncid,'nobs',nobs_id))
      
         ! Create variables
         call ncerr('main',nf90_def_var(subgrpid,'xg',nc_kind_real,(/nx_id,ny_id/),xg_id))
         call ncerr('main',nf90_def_var(subgrpid,'hxg',nc_kind_real,(/nobs_id/),hxg_id))
         call ncerr('main',nf90_def_var(subgrpid,'d',nc_kind_real,(/nobs_id/),d_id))
      
         ! Write variable
         call ncerr('main',nf90_put_var(subgrpid,xg_id,xg_2d))
         call ncerr('main',nf90_put_var(subgrpid,hxg_id,hxg))
         call ncerr('main',nf90_put_var(subgrpid,d_id,algo(io,ia,im)%d))
      
         ! Write algo data
         call algo(io,ia,im)%write(geom(io),grpid(io,im),subgrpid)
      
         ! Release memory
         deallocate(xb)
         deallocate(xg_2d)
      end do

      ! Release memory
      if ((trim(method(im))=='theoretical').or.(trim(method(im))=='standard')) deallocate(dxa_full)
   end do
end do
write(*,'(a)') ''

!--------------------------------------------------------------------------------
! General comparison
!--------------------------------------------------------------------------------
write(*,'(a)') 'General comparison'
write(*,'(a)') ''

! Allocation
allocate(maxdiff(na,na,nm,nm))

do io=1,no
   write(*,'(a,i1,a,i4,a,i4)') '   Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny
   write(*,'(a)') ''

   ! Initialization
   maxdiff = zero

   ! Compute maximum difference
   do im=1,nm
      do jm=1,nm
         do ia=1,na
            do ja=1,na
               ! Compute max relative difference between methods/preconditionings
               do ii=0,ni
                  maxdiff(ia,ja,im,jm) = max(maxdiff(ia,ja,im,jm), &
 & two*abs(algo(io,ia,im)%j_quad(ii)-algo(io,ja,jm)%j_quad(ii))/abs(algo(io,ia,im)%j_quad(ii)+algo(io,ja,jm)%j_quad(ii)))
               end do
            end do
         end do
      end do
   end do

   write(*,'(a16)',advance='no') ''
   do ja=1,na
      do jm=1,nm
         write(*,'(a16)',advance='no') '################'
      end do
   end do
   write(*,'(a)') '#'
   write(*,'(a16)',advance='no') ''
   do ja=1,na
      do jm=1,nm
         write(*,'(a,a13,a)',advance='no') '# ',algorithm(ja),' '
      end do
   end do
   write(*,'(a)') '#'
   write(*,'(a16)',advance='no') ''
   do ja=1,na
      do jm=1,nm
         write(*,'(a,a13,a)',advance='no') '# ',method(jm),' '
      end do
   end do
   write(*,'(a)') '#'
   write(*,'(a16)',advance='no') '################'
   do ja=1,na
      do jm=1,nm
         write(*,'(a16)',advance='no') '################'
      end do
   end do
   write(*,'(a)') '#'
   do ia=1,na
      do im=1,nm
         write(*,'(a,a13,a)',advance='no') '# ',algorithm(ia),' '
         do ja=1,na
            do jm=1,nm
               write(*,'(a,e13.6,a)',advance='no') '# ',maxdiff(ia,ja,im,jm),' '
            end do
         end do
         write(*,'(a)') '#'
         write(*,'(a,a13,a)',advance='no') '# ',method(im),' '
         do ja=1,na
            do jm=1,nm
               if (maxdiff(ia,ja,im,jm)<threshold) then
                  write(*,'(a,a,a13,a,a)',advance='no') '# ',char(27)//'[0;32m','    similar     ',char(27)//'[0;0m',' '
               else
                  write(*,'(a,a,a13,a,a)',advance='no') '# ',char(27)//'[0;91m','   different    ',char(27)//'[0;0m',' '
               end if
            end do
         end do
         write(*,'(a)') '#'
         write(*,'(a16)',advance='no') '################'
         do ja=1,na
            do jm=1,nm
               write(*,'(a16)',advance='no') '################'
            end do
         end do
         write(*,'(a)') '#'
      end do
   end do
   write(*,'(a)') ''
end do

! Close NetCDF file
call ncerr('main',nf90_close(ncid))

end program main
