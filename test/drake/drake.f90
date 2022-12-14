program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius_earth   = 6371.229d0 * KM                      ! radius of Earth
  omega_earth    = 7.29211d-5 * RAD/SECOND              ! rotation rate of Earth
  grav_accel     = 9.80616d0  * METRE/SECOND**2         ! gravitational acceleration 
  ref_density    = 1028d0     * KG/METRE**3             ! reference density at depth (seawater)
  
  ! Numerical method parameters
  timeint_type       = "RK3"                            ! time scheme
  match_time         = .true.                           ! avoid very small time steps when saving 
  mode_split         = .true.                           ! split barotropic mode if true
  penalize           = .true.                           ! penalize land regions
  compressible       = .false.                          ! always run with incompressible equations
  remapscalar_type   = "PPR"                           ! remapping scheme for scalars
  remapvelo_type     = "PPR"                           ! remapping scheme for velocity

  coarse_iter         = 40                              ! maximum number of coarse scale bicgstab iterations for elliptic solver
  fine_iter           = 100                             ! maximum number of fine scale jacobi iterations for elliptic solver
  tol_elliptic        = 1d-6                            ! tolerance for coarse scale bicgstab elliptic solver
  tol_jacobi          = 1d-3                            ! tolerance for fine scale jacobi iterations
  
  Laplace_order_init = 1                                ! use Laplacian viscosity
  nstep_init         = 10                               ! take nstep_init small steps on restart
  cfl_num            = 25d0                             ! cfl number
  save_zlev          = zlevels                          ! vertical layer to save

  npts_penal         = 4.5d0                            ! smooth mask over this many grid points 
  etopo_coast        = .false.                          ! etopo data for coastlines (i.e. penalization)
  etopo_res          = 4                                ! resolution of etopo data in arcminutes

  if (zlevels == 1) then                                ! maximum allowed depth (must be negative)
     max_depth            = -4000d0 * METRE             ! total depth
     halocline            = -4000d0 * METRE             ! location of top (less dense) layer in two layer case
     thermocline          = -4000d0 * METRE             ! location of layer forced by surface wind stress
     drho                 =     0d0 * KG/METRE**3       ! density perturbation at free surface 
     tau_0                =   0.4d0 * NEWTON/METRE**2   ! maximum wind stress
     bottom_friction_case =    5d-3 * METRE/SECOND      ! bottom friction   
     u_wbc                =   1.5d0 * METRE/SECOND      ! estimated western boundary current speed
  elseif (zlevels == 2) then
     remap                = .true.
     iremap               = 10
     max_depth            = -4000d0 * METRE             ! total depth
     thermocline          = -1000d0 * METRE             ! location of layer forced by surface wind stress
     halocline            = -1000d0 * METRE             ! location of top (less dense) layer in two layer case
     drho                 =    -8d0 * KG/METRE**3       ! density perturbation at free surface (density of top layer is rho0 + drho/2)
     tau_0                =   0.4d0 * NEWTON/METRE**2   ! maximum wind stress
     bottom_friction_case =    5d-3 * METRE/SECOND      ! bottom friction
     Ku                   =     4d0 * METRE**2/SECOND   ! viscosity for vertical diffusion
     u_wbc                =   1.5d0 * METRE/SECOND      ! estimated western boundary current speed
  elseif (zlevels >= 3) then
     sigma_z              = .true.                      ! sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     vert_diffuse         = .true.
     tke_closure          = .false.
     remap                = .true.
     iremap               =   5
     thermocline          =    -100d0 * METRE           ! location of surface mixed layer
     halocline            =   -1100d0 * METRE           ! location less dense layer 
     drho                 = -10.481d0 * KG/METRE**3     ! density perturbation at free surface (density of top layer is rho0 + drho/2)
     tau_0                =     0.4d0 * NEWTON/METRE**2 ! maximum wind stress
     max_depth            =   -4000d0 * METRE           ! total depth
     ref_density          = 1027.75d0 * KG/METRE**3     ! reference density
     bottom_friction_case =      5d-3 * METRE/SECOND    ! bottom friction                      
     u_wbc                =       1d0 * METRE/SECOND    ! estimated western boundary current speed
  end if

  ! Characteristic scales
  radius         = radius_earth/scale                           ! mean radius of the small planet
  omega          = omega_earth/scale_omega                      ! angular velocity (scaled for small planet to keep beta constant)
  wave_speed     = sqrt (grav_accel*abs(max_depth))             ! inertia-gravity wave speed 
  f0             = 2d0*omega*sin(45d0*DEG)                      ! representative Coriolis parameter
  beta           = 2d0*omega*cos(45d0*DEG) / radius             ! beta parameter at 45 degrees latitude
  Rd             = wave_speed / f0                              ! barotropic Rossby radius of deformation                   
  drho_dz        = drho / halocline                             ! density gradient
  bv             = sqrt (grav_accel * abs(drho_dz)/ref_density) ! Brunt-Vaisala frequency
  delta_I        = sqrt (u_wbc/beta)                            ! inertial layer
  delta_sm       = u_wbc / f0                                   ! barotropic submesoscale

  ! Baroclinic wave speed
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

  ! Relaxation of buoyancy to mean profile 
  if (zlevels == 2 .and. remap) then
     k_T = 1d0 / (30d0 * DAY)               
  else
     k_T = 0d0
  end if

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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Parameters that require viscosity
  delta_M = (visc_rotu/beta)**(1d0/(2*Laplace_order_init+1))  ! Munk layer scale
  Rey     = u_wbc * delta_I / visc_rotu                       ! Reynolds number of western boundary current
  Ro      = u_wbc / (delta_M*f0)                              ! Rossby number (based on boundary current)

  if (zlevels <= 2) bottom_friction_case = beta * delta_M/4d0 ! ensure that delta_S = delta_M/4
  
  delta_S = bottom_friction_case / (abs(max_depth) * beta)    ! Stommel layer scale

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set interval for adapting grid based on the horizontal advective velocity scale (i.e. advect no more than one grid point before adapting)
  iadapt     = 1        ! Drake unstable with trend computation over entire grid if iadapt > 1
  irebalance = 4*iadapt ! rebalance interval using charm++/AMPI

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
     if (k_T /= 0d0) call euler (sol, wav_coeff, trend_relax, dt)
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

