program turb
  ! Homogeneous turbulence test case.
  ! Initials conditions are a set of 17 baroclinic jets in geostrophic balance.
  use main_mod
  use test_case_mod
  use io_vtk_mod
  implicit none
  logical :: aligned
  
  ! Initialize mpi, shared variables and domainss
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  ! Parameters defining domain (based on beta-plane values)
  radius             = 1000d0 * KM                     ! meridional width of zonal channel
  omega              = 6d-5                            ! planet rotation
  grav_accel         = 9.80616d0 * METRE/SECOND**2     ! gravitational acceleration 
  ref_density        = 1028d0 * KG/METRE**3         ! reference density at depth (maximum density)

  ! Numerical method parameters
  default_thresholds = .true.                          ! use default threshold
  
  match_time         = .true.                          ! avoid very small time steps when saving (if false) 
  penalize           = .false.                          ! penalize land regions
  compressible       = .false.                         ! always run with incompressible equations
  remapscalar_type   = "PPR"                           ! optimal remapping scheme
  remapvelo_type     = "PPR"                           ! optimal remapping scheme

  ! Time stepping parameters
  timeint_type       = "RK3"                           ! time integration scheme
  adapt_dt           = .true.                          ! adapt time step
  mode_split         = .true.                          ! split barotropic mode if true
  theta1             = 0.8d0                           ! external pressure gradient (1 = fully implicit, 0.5 = Crank-Nicolson) stable if > 0.75
  theta2             = 0.8d0                           ! barotropic flow divergence (1 = fully implicit, 0.5 = Crank-Nicolson) stable if > 0.75
  
  ! Horizontal diffusion
  Laplace_order_init = 2                         
  Laplace_order      = Laplace_order_init

  ! Vertical diffusion
  vert_diffuse       = .true.                       
  tke_closure        = .true.
         
  ! Depth and layer parameters
  sigma_z            = .true.                          ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
  coords             = "croco"                         ! grid type for pure sigma grid ("croco" or "uniform")
  max_depth          = -1000d0 * METRE                 ! total depth
  min_depth          = max_depth                       ! minimum depth
  Tcline             =  -100d0 * METRE                 ! thermocline

  ! Bottom friction
  bottom_friction_case = 5d-3 / SECOND

  ! Relaxation to initial zonal flow (nudging)
  tau_nudge          = 50d0 * DAY

  ! Wind stress
  tau_0              = 0d0

  ! Equation of state variables
  a_0                = 0.28d0 / CELSIUS
  b_0                = 0d0
  mu_1               = 0d0
  mu_2               = 0d0
  mu_1               = 0d0
  mu_2               = 0d0
  T_ref              = 14d0   * CELSIUS
  
  ! Vertical level to save
  save_zlev          = zlevels

  ! Characteristic scales
  drho               = -3d0 * KG/METRE**3                                         ! magnitude of density perturbation (linear stratification)
  grav_reduced       = grav_accel * abs(drho)/ref_density                         ! reduced gravity
  wave_speed         = sqrt (grav_accel*abs(max_depth))                           ! inertia-gravity wave speed
  drho_dz            = abs (drho / max_depth)                                     ! density gradient
  bv                 = sqrt (grav_accel * abs(drho_dz)/ref_density)               ! Brunt-Vaisala frequency
  c1                 = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification
  
  ! Dimensional scaling
  Udim               = 0.5d0 * METRE/SECOND              ! velocity scale (initial maximum initial zonal jet velocity)
  Ldim               = 10d0 * DEG * radius               ! length scale (initial width of zonal jets)
  Tdim               = Ldim/Udim                         ! time scale
  Hdim               = abs (max_depth)                   ! vertical length scale

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions
  
  ! Initialize variables
  call initialize (run_id)

  ! Set interval for adapting grid based on the horizontal advective velocity scale (i.e. advect no more than one grid point before adapting)
  iadapt = CFL_adv * nint ((dx_min/Udim) / dt_init)

  ! Initialize 2D projections and zonal averages
  Nproj = sqrt (20d0 * 4**min_level) ! size of 2D projection: Nproj x Nproj/2
  call initialize_projection (Nproj)
  allocate (y2(Ny(1):Ny(2),1:zlevels,1:4),       y2_0(Ny(1):Ny(2),1:zlevels,1:4))
  allocate (zonal(Ny(1):Ny(2),1:zlevels,1:4), zonal_0(Ny(1):Ny(2),1:zlevels,1:4))
  call zonal_mean (zonal_0, y2_0) ; zonal = zonal_0 ; y2 = y2_0

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

     if (tau_nudge /= 0d0) then
        if (modulo (istep, 10) == 0) call zonal_mean (zonal, y2)
        call euler (sol, wav_coeff, trend_nudge, dt)
     end if

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
end program turb
