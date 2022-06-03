program upwelling
  ! Simple wind-driven upwelling/downwelling in a periodic zonal channel.
  ! An extension to the sphere of a test case used in ROMS and CROCO.
  use main_mod
  use test_case_mod
  use io_mod
  implicit none
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius             = 240d0     * KM
  omega              = 6d-5      * RAD/SECOND
  grav_accel         = 9.80616d0 * METRE/SECOND**2   ! gravitational acceleration 
  p_top              = 0d0       * hPa               ! pressure at free surface
  ref_density        = 1027d0    * KG/METRE**3       ! reference density at depth (seawater)

  ! Numerical method parameters
  default_thresholds = .true.                        ! use default threshold
  adapt_dt           = .true.                        ! adapt time step
  match_time         = .true.                        ! avoid very small time steps when saving (if false) 
  mode_split         = .true.                        ! split barotropic mode if true
  penalize           = .true.                        ! penalize land regions
  alpha              = 1d-6                          ! porosity used in penalization
  npts_penal         = 2.5                           ! number of points to smooth over in penalization
  timeint_type       = "RK4"                         ! always use RK4
  compressible       = .false.                       ! always run with incompressible equations
  remapscalar_type   = "PPR"                         ! optimal remapping scheme
  remapvelo_type     = "PPR"                         ! optimal remapping scheme
  
  Laplace_order_init = 1                              
  Laplace_order = Laplace_order_init
  implicit_diff_sclr = .false.
  implicit_diff_divu = .false.
  vert_diffuse       = .true.                        ! include vertical diffusion
  tke_closure        = .false.                       ! use analytic scheme for eddy viscosity and eddy diffusion
         
  ! Depth and layer parameters
  sigma_z            = .true.                        ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
  coords             = "croco"                       ! grid type for pure sigma grid ("croco" or "uniform")
  max_depth          = -150d0 * METRE                ! total depth
  min_depth          =  -25d0 * METRE                ! minimum depth
  Tcline             =  -50d0 * METRE                ! position of thermocline

  ! Land mass parameters
  width              = 80d0 * KM                     ! width of zonal channel
  lat_width          = (width/radius)/DEG            ! width of zonal channel (in degrees)
  lat_c              = 45d0                          ! centre of zonal channel (in degrees)
  
  ! Bottom friction
  bottom_friction_case = 3d-4 / SECOND          

  ! Wind stress
  tau_0              = 0.1d0

  ! Equation of state variables
  a_0                = 0.28d0 / CELSIUS
  b_0                = 0d0
  mu_1               = 0d0
  mu_2               = 0d0
  T_ref              = 14d0   * CELSIUS
  
  ! Vertical level to save
  save_zlev          = 8

  ! Characteristic scales
  wave_speed         = sqrt (grav_accel*abs(max_depth))   ! inertia-gravity wave speed 
  f0                 = 2d0*omega*sin(45d0*DEG)            ! representative Coriolis parameter
  beta               = 2d0*omega*cos(45d0*DEG) / radius   ! beta parameter at 30 degrees latitude
  Rd                 = wave_speed / f0                    ! barotropic Rossby radius of deformation                   

  ! Dimensional scaling
  drho               = 2.5d0 * KG/METRE**3                ! magnitude of density perturbation
  Udim               = 0.35d0 * METRE/SECOND              ! velocity scale
  Ldim               = lat_width*DEG * radius             ! length scale 
  Tdim               = Ldim/Udim                          ! time scale
  Hdim               = abs (max_depth)                    ! vertical length scale
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions
  
  ! Initialize variables
  call initialize (run_id)

  ! Initialize diagnostic variables
  call init_diagnostics

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
     call stop_timing

     call update_diagnostics

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (run_id, rebalance)
           call init_diagnostics !! resets diagnostics !!
        end if

        ! Save fields
        call vertical_velocity
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program upwelling
