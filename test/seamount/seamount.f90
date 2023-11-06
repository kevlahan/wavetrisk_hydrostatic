program Seamount
  ! Seamount case evaluates pressure gradient error (Beckmann and Haidvogel 1993). An isolated Gaussian bathymetry with flat isopycnals
  ! and hybrid sigma coordinates in z (vertical levels not aligned with isopycnals). The exact solution is a state of rest.
  ! The density difference is -3 kg/m^3 and decreases exponentially to zero with increasing depth. Parameters are set as in
  ! Debreu and Kevlahan (2020).
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  integer :: l
  real(8) :: total_eta
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius_earth   = 6371.229                * KM
  omega_earth    = 7.29211d-5              * RAD/SECOND
  grav_accel     = 9.80616                 * METRE/SECOND**2 ! gravitational acceleration 
  p_top          = 0.0_8                   * hPa             ! pressure at free surface
  ref_density    = 1000                    * KG/METRE**3     ! reference density at depth (seawater)

  scale          = 41.75d0                                   ! scaling factor for small planet
  radius         = radius_earth/scale                        ! mean radius of the small planet
  omega          = omega_earth                              ! angular velocity (scaled for small planet to keep beta constant)
  !omega          = 1d-4/2d0                * RAD/SECOND      ! angular velocity for seamount at north pole

  ! Depth and layer parameters
  min_depth      =   -50 * METRE                             ! minimum allowed depth (must be negative)
  max_depth      = -5000 * METRE                             ! total depth
  tau_0          =     0 * NEWTON/METRE**2                   ! maximum wind stress
  
  ! Seamount
  lat_c          = 43.29 * DEG                               ! latitude of seamount
  !lat_c          = 90d0 * DEG
  lon_c          =     0 * DEG                               ! longitude
  h0             =  4500 * METRE                             ! height of seamount
  width          =    40 * KM                                ! radius of seamount
  delta          =   500 * METRE                             ! vertical decay of density
  
  ! Numerical method parameters
  match_time         = .false.                               ! avoid very small time steps when saving 
  mode_split         = .true.                                ! split barotropic mode if true
  penalize           = .false.                               ! no penalization
  timeint_type       = "RK4"                                 ! always use RK4
  compressible       = .false.                               ! always run with incompressible equations
  default_thresholds = .true.                                ! always use default threshold
  adapt_dt           = .true.                                ! always adapt time steps
  remapscalar_type   = "PPR"                                 ! optimal remapping scheme
  remapvelo_type     = "PPR"                                 ! optimal remapping scheme

  Laplace_order_init = 1                             
  Laplace_order      = Laplace_order_init
  visc               = 50d0 * METRE**2/SECOND                ! viscosity
  
  theta1             = 0.8d0
  theta2             = 0.8d0

  vert_diffuse       = .false.                               ! no vertical diffusion

  ! Characteristic scales
  wave_speed     = sqrt (grav_accel*abs(max_depth))          ! inertia-gravity wave speed 
  f0             = 2d0*omega*sin(lat_c)                      ! representative Coriolis parameter
  beta           = 2d0*omega*cos(lat_c) / radius             ! beta parameter at 30 degrees latitude
  Rd             = wave_speed / f0                           ! barotropic Rossby radius of deformation                   
  Ro             = 0d0                                       ! Rossby number

  ! Bottom friction
  bottom_friction_case = 0 / SECOND

  ! Dimensional scaling
  Udim           = 1d0                                       ! velocity scale (arbitrary)
  Ldim           = width                                     ! length scale 
  Tdim           = Ldim/Udim                                 ! time scale
  Hdim           = abs (max_depth)                           ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read test case parameters
  call read_test_case_parameters
  
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
  total_cpu_time = 0d0
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
     call stop_timing

     call update_diagnostics

     tke_sea = total_ke ("adaptive") / (4*MATH_PI*radius**2)
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (run_id, rebalance)
           call init_diagnostics !! resets diagnostics !!
        end if

        ! Save fields
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program Seamount
