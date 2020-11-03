program main

!use interp
use tools_rand
use type_algo
use type_bmatrix
use type_geom
use type_hmatrix
use type_lmp
use type_rmatrix


implicit none

! Maximum parameters
integer,parameter   :: nomax = 10                    ! Maximum number of outer iterations

! Namelist parameters
integer          :: no                            ! Number of outer iterations
integer          :: ni                            ! Number of inner iterations
character(len=8) :: lmp_mode                      ! LMP mode ('none', 'spectral', 'ritz')
integer          :: shutoff_type                  ! Stop criterion according: 1-Jb, 2-beta, (else: no criterion)
real(8)          :: shutoff_value                 ! Stop criterion threshold
logical          :: consistent_strategy           ! Use the consistent strategy detailed in the note
logical          :: use_binv                      ! Use B inverse in the inconsistent strategy
integer          :: nx(nomax)                     ! X direction size
integer          :: ny(nomax)                     ! Y direction size
integer          :: nobs                          ! Number of observations
real(8)          :: sigma_obs                     ! Observation error standard deviation
real(8)          :: sigmabvar                     ! Grid-point standard deviation variations amplitude
real(8)          :: Lb                            ! Correlation length-scale
logical          :: new_seed                      ! New random seed

! Output files units
integer,parameter :: namelist_unit = 9
integer,parameter :: delta_test_unit = 10
integer,parameter :: Bdelta_test_unit = 11
integer,parameter :: Hdelta_test_unit = 12
integer,parameter :: lanczos_control_space_unit = 13
integer,parameter :: lanczos_control_space_outer_grid_unit = 14
integer,parameter :: lanczos_control_space_outer_obs_unit = 15
integer,parameter :: PlanczosIF_model_space_unit = 16
integer,parameter :: PlanczosIF_model_space_outer_grid_unit = 17
integer,parameter :: PlanczosIF_model_space_outer_obs_unit = 18
integer,parameter :: lanczos_control_vs_PlanczosIF_model_unit = 19

! Local variables
integer                        :: io,jo,ii
real(8),allocatable            :: xb_full(:),xg_full(:),dxb_full(:),dxa_full(:)
real(8),allocatable            :: xb(:),xg(:),dxb(:),dxb_bar(:),dvb_bar(:),dxa_prev(:),dva_interp(:)
real(8),allocatable            :: d(:),hxg(:)
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
 & shutoff_type, &
 & shutoff_value, &
 & consistent_strategy, &
 & use_binv
namelist/resolutions/ &
 & nx, &
 & ny
namelist/observations/ &
 & nobs, &
 & sigma_obs
namelist/background/ &
 & sigmabvar, &
 & Lb
namelist/miscellanous/ &
 & new_seed

! Read namelist
open(unit=namelist_unit,file='namelist',status='old',action='read')
read(namelist_unit,nml=solver)
read(namelist_unit,nml=resolutions)
read(namelist_unit,nml=observations)
read(namelist_unit,nml=background)
read(namelist_unit,nml=miscellanous)

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
allocate(geom(no))
allocate(bmatrix(no))
allocate(algo_lanczos(no))
allocate(algo_planczosif(no))
allocate(lmp_lanczos(no))
allocate(lmp_planczosif(no))

! Setup geometries
do io=1,no
   write(*,'(a,i2)') '   Geometry setup for outer iteration ',io
   call geom(io)%setup(nx(io),ny(io))
end do

! Allocation (model/control space)
allocate(xb_full(geom(no)%nh))
allocate(xg_full(geom(no)%nh))

! Setup B matrices
do io=1,no
   write(*,'(a,i2)') '   B matrix setup for outer iteration ',io
   call bmatrix(io)%setup(geom(io),geom(no),sigmabvar,Lb)
end do

! Generate background state
call bmatrix(no)%randomize(geom(no),xb_full)

! Allocation (obs space)
allocate(d(nobs))
allocate(hxg(nobs))

! Setup obserations locations
call hmatrix%setup(nobs)

! Setup R matrix
call rmatrix%setup(nobs,sigma_obs)

! Setup observations
call rmatrix%randomize(hmatrix%yo)

! Open files
open(delta_test_unit,file='results/delta_test.dat')
open(Bdelta_test_unit,file='results/Bdelta_test.dat')
open(Hdelta_test_unit,file='results/Hdelta_test.dat')
open(lanczos_control_space_unit,file='results/lanczos_control_space.dat')
open(lanczos_control_space_outer_grid_unit,file='results/lanczos_control_space_outer_grid.dat')
open(lanczos_control_space_outer_obs_unit,file='results/lanczos_control_space_outer_obs.dat')
open(PlanczosIF_model_space_unit,file='results/PlanczosIF_model_space.dat')
open(PlanczosIF_model_space_outer_grid_unit,file='results/PlanczosIF_model_space_outer_grid.dat')
open(PlanczosIF_model_space_outer_obs_unit,file='results/PlanczosIF_model_space_outer_obs.dat')
open(lanczos_control_vs_PlanczosIF_model_unit,file='results/lanczos_control_vs_PlanczosIF_model.dat')

! Write headers
write(Bdelta_test_unit,'(a)') '# B matrix: 1:outer_loop 2:line 3:column 4:element value'
write(Hdelta_test_unit,'(a)') '# H matrix: 1:outer_loop 2:line 3:column 4:element value'
write(lanczos_control_space_unit,'(a)') '# Outer iteration , resolution , Inner iteration , J=Jb+Jo , Jb , Jo, sqrt(rho), beta'
write(lanczos_control_space_outer_grid_unit,'(a)') '# outer iteration, indices, coord, dva_interp, dvb, dxb, xb, xg'
write(lanczos_control_space_outer_obs_unit,'(a)') '# outer iteration, indices, coord, hxg, yo, d'
write(PlanczosIF_model_space_unit,'(a)') '# Outer iteration , resolution , Inner iteration , J=Jb+Jo , Jb , Jo, sqrt(rho), beta'
write(PlanczosIF_model_space_outer_grid_unit,'(a)') '# outer iteration, indices, coord, dxabar_interp, dxbbar, dxb, xb, xg'
write(PlanczosIF_model_space_outer_obs_unit,'(a)') '# outer iteration, indices, coord, hxg, yo, d'
write(lanczos_control_vs_PlanczosIF_model_unit,'(a)') '# Outer iteration , resolution , Inner iteration , delta_J , delta_Jb , delta_Jo, sqrt(rho), beta'

write(*,'(a)') ''

!--------------------------------------------------------------------------------
! Multi-incremental Lanczos in control space
!--------------------------------------------------------------------------------
write(*,'(a)') 'Multi-incremental Lanczos in control space'
do io=1,no
   write(*,'(a,i2,a,i4,a,i4)') '   Outer iteration ',io,', resolution: ',geom(io)%nx,' x ',geom(io)%ny

   ! Allocation
   allocate(xb(geom(io)%nh))
   allocate(xg(geom(io)%nh))
   allocate(dvb_bar(geom(io)%nh))
   call algo_lanczos(io)%alloc(geom(io),ni)
   call lmp_lanczos(io)%alloc(geom(io),ni,io,lmp_mode,'control')

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
      dvb_bar = 0.0

      do jo=1,io-1
         ! Interpolate analysis increment of previous iterations at current resolution
         call geom(jo)%interp_sp(geom(io),algo_lanczos(jo)%dva,dva_interp)

         ! Add contributions
         dvb_bar = dvb_bar-dva_interp
      end do

      ! Compute background increment at current resolution
      call bmatrix(io)%apply_sqrt(geom(io),dvb_bar,dxb)

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
         allocate(dxb_bar(geom(io)%nh))

         ! Compute background increment at current resolution
         dxb = xg-xb

         ! Apply B inverse
         call bmatrix(io)%apply_inv(geom(io),dxb,dxb_bar)

         ! Apply B square-root
         call bmatrix(io)%apply_sqrt_ad(geom(io),dxb_bar,dvb_bar)

         ! Release memory
         deallocate(dxb)
         deallocate(dxb_bar)
      else
         ! Allocation
         allocate(dva_interp(geom(io)%nh))

         ! Initialization
         dvb_bar = 0.0

         do jo=1,io-1
            ! Interpolate analysis increment of previous iterations at current resolution
            call geom(jo)%interp_sp(geom(io),algo_lanczos(jo)%dva,dva_interp)

            ! Add contributions
            dvb_bar = dvb_bar-dva_interp
         end do

         ! Release memory
         deallocate(dva_interp)
      end if
   end if

   ! Compute innovation
   call hmatrix%apply(geom(io),xg,hxg)
   d = hmatrix%yo-hxg

!   if (io>1) then
!      ! Copy preconditioning vectors
!      select case (trim(lmp_mode))
!      case ('spectral','ritz')
!         ! Copy eigenpairs
!         lmp_lanczos(io)%outer(io)%eigenval = algo_lanczos(io-1)%eigenval
!         lmp_lanczos(io)%outer(io)%eigenvec = algo_lanczos(io-1)%eigenvec

         ! Copy Lanczos vectors
!         lmp_lanczos(io)%outer(io)%lancvec_trunc(1:nn(io-1),:) = algo_lanczos(io-1)%lancvec(1:nn(io-1),:)

         ! Compute omega
!         lmp_lanczos(io)%outer(io)%omega = algo_lanczos(io-1)%eigenvec(ni,:)*algo_lanczos(io-1)%lastbeta/algo_lanczos(io-1)%eigenval
!      end select

      ! Interpolate preconditioning vectors
!      call interp_lmp_control(ni,io,nn(1:io),lmp_lanczos(1:io))
!   end if

   ! Minimization
   write(*,'(a)') '      Minimization (J = Jb+Jo / J_nl = Jb_nl+Jo_nl)'
   call algo_lanczos(io)%apply_lanczos(geom(io),bmatrix(io),hmatrix,rmatrix,xb,xg,dvb_bar,d,ni,lmp_lanczos(io),shutoff_type,shutoff_value)


   ! Result
   do ii=0,ni
      write(*,'(a,i3,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4,a,e11.4)') '         Iteration ',ii,':',algo_lanczos(io)%jb(ii)+algo_lanczos(io)%jo(ii),' = ',algo_lanczos(io)%jb(ii),' + ',algo_lanczos(io)%jo(ii),' / ',algo_lanczos(io)%jb_nl(ii)+algo_lanczos(io)%jo_nl(ii),' = ',algo_lanczos(io)%jb_nl(ii),' + ',algo_lanczos(io)%jo_nl(ii)
   end do

   ! Release memory
   deallocate(xb)
   deallocate(xg)
   deallocate(dvb_bar)
end do
write(*,'(a)') ''

! Close files
close(delta_test_unit)
close(Bdelta_test_unit)
close(Hdelta_test_unit)
close(lanczos_control_space_unit)
close(lanczos_control_space_outer_grid_unit)
close(lanczos_control_space_outer_obs_unit)
close(PlanczosIF_model_space_unit)
close(PlanczosIF_model_space_outer_grid_unit)
close(PlanczosIF_model_space_outer_obs_unit)
close(lanczos_control_vs_PlanczosIF_model_unit)

end program main
