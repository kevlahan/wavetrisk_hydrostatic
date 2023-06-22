program Tsunami
  ! Simple tsunami test case
  !
  ! The permeability penalization parameter eta = dt_cfl and the porosity parameter alpha is set in the input file.
  use main_mod
  use test_case_mod
  use io_mod
  implicit none
  integer                                :: N
  integer, dimension(2)                  :: Nx, Ny
  real(8), dimension(:,:),   allocatable :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export
  real(8), dimension(2)                  :: lon_lat_range
  real(8), dimension(:,:,:), allocatable :: field2d_save
  logical                                :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371.229d0 * KM              ! mean radius of the Earth
  grav_accel     = 9.80616d0  * METRE/SECOND**2 ! gravitational acceleration 
  omega          = 7.29211d-5 * RAD/SECOND      ! Earth’s angular velocity 
  ref_density    = 1027d0     * KG/METRE**3     ! reference density (seawater)

  ! Bathymetry parameters
  min_depth      = -50d0  * METRE               ! minimum allowed depth (must be negative)
  max_depth      = -4d0   * KM                  ! maximum allowed depth (must be negative)
  etopo_res      = 4                            ! resolution of etopo data in arcminutes (if used)
  penalize       = .true.                       ! penalize land regions 
  etopo_coast    = .true.                       ! use etopo data for coastlines (i.e. penalization)
  etopo_bathy    = .true.                       ! use etopo data for bathymetry
  npts_penal     = 4d0                          ! smooth penalization mask over npts_penal grid points
  npts_topo      = 4d0                          ! smooth over topography over npts_topo grid points	

  dH             =  7d0   * METRE               ! initial perturbation to the free surface
  pert_radius    =  1d3 * KM                    ! radius of Gaussian free surface perturbation
  lon_c          = -50d0  * DEG                 ! longitude location of perturbation
  lat_c          =  25d0  * DEG                 ! latitude  location of perturbation

  ! Parameters for 2D projection
  N              = 1024                         ! size of lat-lon grid in 2D projection
  lon_lat_range  = (/2d0*MATH_PI, MATH_PI/)     ! region to save in 2D projection

  ! Numerical parameters
  mode_split         = .false.                  ! use explicit time step for accuracy
  timeint_type       = "RK33"                   ! time scheme
  cfl_num            = 1d0                      ! CFL number
  
  compressible       = .false.                  ! incompressible
  adapt_dt           = .true.                   ! adapt time step
  remap              = .false.                  ! remap vertical layers

  Laplace_order_init = 2                        ! viscosity type
  C_visc(S_MASS)     = 0d-3                     ! dimensionless viscosity of S_MASS
  C_visc(S_TEMP)     = 0d-3                     ! dimensionless viscosity of S_TEMP
  C_visc(S_VELO)     = 1d-3                     ! dimensionless viscosity of S_VELO (rotu, divu)
  
  log_mass           = .true.                   ! mass diagnostics

  ! Dimensional scaling
  wave_speed         = sqrt (grav_accel*abs(max_depth)) ! inertia-gravity wave speed based on maximum allowed depth
  Udim               = wave_speed                   ! velocity scale
  Ldim               = 2d0 * pert_radius              ! length scale (free surface perturbation width)
  Tdim               = Ldim/Udim                    ! time scale
  Hdim               = abs (max_depth)              ! vertical length scale

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0d0
  do while (time < time_end)
     call start_timing;  call time_step (dt_write, aligned); call stop_timing

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)

        ! Save fields
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program Tsunami





