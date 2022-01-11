program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  !
  ! The permeability penalization parameter eta = dt_cfl and the porosity parameter alpha is set in the input file.
  ! The initial condition has a 500 m thick less dense layer over a 1500 m thick reference density layer (only the upper
  ! layer is included in the single vertical layer shallow water case).
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  integer :: l
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius_earth   = 6371.229d0                * KM
  omega_earth    = 7.29211d-5                * RAD/SECOND
  grav_accel     = 9.80616d0                 * METRE/SECOND**2 ! gravitational acceleration 
  p_top          = 0d0                       * hPa             ! pressure at free surface
  ref_density    = 1028d0                    * KG/METRE**3     ! reference density at depth (seawater)

  radius         = radius_earth/scale                          ! mean radius of the small planet
  omega          = omega_earth/scale_omega                     ! angular velocity (scaled for small planet to keep beta constant)

  ! Numerical method parameters
  match_time         = .true.                        ! avoid very small time steps when saving 
  mode_split         = .true.                         ! split barotropic mode if true
  timeint_type       = "RK3"                          
  rebalance          = .true.
  penalize           = .true.                         ! penalize land regions
  compressible       = .false.                        ! always run with incompressible equations
  remapscalar_type   = "PPR"                          ! optimal remapping scheme
  remapvelo_type     = "PPR"                          ! optimal remapping scheme
  Laplace_order_init = 1                              
  Laplace_order      = Laplace_order_init
  nstep_init         = 5                               ! take nstep_init small steps on restart
  log_mass           = .true.

  ! Depth and layer parameters
  etopo_res      = 4                                    ! resolution of etopo data in arcminutes (if used) 
  etopo_coast    = .false.                              ! use etopo data for coastlines (i.e. penalization)
  min_depth      =   -50d0 * METRE                      ! minimum allowed depth (must be negative)
  if (zlevels == 1) then                                ! maximum allowed depth (must be negative)
     max_depth   = -4000d0 * METRE                      ! total depth
     halocline   = -4000d0 * METRE                      ! location of top (less dense) layer in two layer case
     mixed_layer = -4000d0 * METRE                      ! location of layer forced by surface wind stress
     drho        =     0d0 * KG/METRE**3                ! density perturbation at free surface 
     tau_0       =   0.8d0 * NEWTON/METRE**2            ! maximum wind stress
     u_wbc       =   1.5d0 * METRE/SECOND               ! estimated western boundary current speed
  elseif (zlevels == 2) then
     max_depth   = -4000d0 * METRE                      ! total depth
     halocline   = -1000d0 * METRE                      ! location of top (less dense) layer in two layer case
     mixed_layer = -1000d0 * METRE                      ! location of layer forced by surface wind stress
     drho        =    -8d0 * KG/METRE**3                ! density perturbation at free surface (density of top layer is rho0 + drho/2)
     tau_0       =   0.4d0 * NEWTON/METRE**2            ! maximum wind stress
     u_wbc       =   1.7d0 * METRE/SECOND               ! estimated western boundary current speed
  elseif (zlevels >= 3) then
     max_depth   = -4000d0 * METRE                      ! total depth
     halocline   = -1000d0 * METRE                      ! location of top (less dense) layer in two layer case
     mixed_layer =  -100d0 * METRE                      ! location of layer forced by surface wind stress
     drho        =    -8d0 * KG/METRE**3                ! density perturbation at free surface (linear stratification from free surface to halocline)
     tau_0       =   0.4d0 * NEWTON/METRE**2            ! maximum wind stress
     u_wbc       =   1.7d0 * METRE/SECOND               ! estimated western boundary current speed
  end if

  ! Vertical level to save
  save_zlev = zlevels 

  ! Characteristic scales
  wave_speed     = sqrt (grav_accel*abs(max_depth))        ! inertia-gravity wave speed 
  f0             = 2d0*omega*sin(45*DEG)                     ! representative Coriolis parameter
  beta           = 2d0*omega*cos(45*DEG) / radius            ! beta parameter at 30 degrees latitude
  Rd             = wave_speed / f0                         ! barotropic Rossby radius of deformation                   
  drho_dz        = drho / halocline                        ! density gradient
  bv             = sqrt (grav_accel * abs(drho_dz)/ref_density) ! Brunt-Vaisala frequency

  if (zlevels == 2) then
     c1 = sqrt (grav_accel*abs(drho)/2d0/ref_density*halocline*(max_depth-halocline)/abs(max_depth)) ! two-layer internal wave speed
  elseif (zlevels >= 3) then
     c1 = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification
  endif

  ! First baroclinic Rossby radius of deformation
  if (zlevels == 1) then
     Rb = 0d0
  elseif (zlevels == 2) then
     Rb = c1 / f0                                 
  else
     Rb = bv * abs(max_depth) / (MATH_PI*f0)
  end if

  ! Bottom friction
  if (drag) then
     bottom_friction_case = 4d-4 / abs(max_depth) ! nemo value
  else
     bottom_friction_case = 0d0 / SECOND
  end if

  ! Internal wave friction 
  if (drho == 0d0 .or. remap) then
     wave_friction = 0d0
  else
     wave_friction = 0d0*u_wbc / Rb / 3d0                   ! three e-folding growth times of internal wave (requires accurate u_wbc estimate)
!!$     wave_friction = 1d0/ (200d0 * HOUR)                 ! fixed
  end if

  ! Relaxation of buoyancy to mean profile 
  if (remap) then
     k_T = 1d0 / (30d0 * DAY)               
  else
     k_T = 0d0
  end if

  delta_I        = sqrt (u_wbc/beta)                  ! inertial layer
  delta_sm       = u_wbc/f0                           ! barotropic submesoscale  

  ! Dimensional scaling
  Udim           = u_wbc                              ! velocity scale
  Ldim           = delta_I                            ! length scale 
  Tdim           = Ldim/Udim                          ! time scale
  Hdim           = abs (max_depth)                    ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions
  
  ! Initialize variables
  call initialize (run_id)

  ! Set interval for adapting grid based on the horizontal advective velocity scale (i.e. advect no more than one grid point before adapting)
  iadapt    = CFL_adv * nint ((dx_min/Udim) / dt_init)
  irebalance = iadapt ! rebalance interval using charm++/AMPI

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  total_cpu_time = 0d0
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
     if (k_T /= 0.0_8) call euler (sol, wav_coeff, trend_relax, dt)
     call stop_timing

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

  if (rank == 0)  write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  call finalize
end program Drake

