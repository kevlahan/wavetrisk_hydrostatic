program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  use io_vtk_mod
  use main_mod
  use test_case_mod
  implicit none
  real(dp) :: Area_min, dx_min, dz, visc
  logical  :: relax = .false.

  call init_arch_mod 
  call init_comm_mpi_mod
  call read_test_case_parameters

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Numerical method parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  timeint_type            = "RK3" 
  adapt_dt                = .true.
  compressible            = .false.
  default_thresholds      = .false. ! needed because do not have good initial estimate of layer dependence on variable norms
  log_min_mass            = .false.
  mode_split              = .true.
  penalize                = .true.                
  split_mean_perturbation = .true.

  if (mode_split) then
     cfl_num              = 15.0_dp
  else
     cfl_num              =  0.9_dp                             
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Earth parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  radius_earth   =   6371.229 * KM                      
  omega_earth    = 7.29211e-5 * RAD/SECOND              
  H_earth        =          4 * KM                     
  g_earth        =    9.80616 * METRE/SECOND**2        
  ref_density    =       1030 * KG/METRE**3
  c_p            =    3991.87 * JOULE/(KG*KELVIN) ! specific heat at constant pressure for seawater

  porosity       = 1.0e-2_dp

  ! Earth scaling factors
  L_norm         = radius_earth
  H_norm         = H_earth
  U_norm         = sqrt (g_earth * H_earth)
  T_norm         = L_norm / U_norm

  ! Use normalized equation
  normalized = .false.
  if (normalized) then
     grav_accel  = g_earth * H_norm / U_norm**2
     omega       = omega_earth * T_norm
     radius      = radius_earth / L_norm                     
  else
     grav_accel  = g_earth
     omega       = omega_earth  / scale_omega              ! angular velocity (scaled for small planet to keep beta constant)
     radius      = radius_earth / scale                    ! mean radius of the small planet
  end if

  f0             = 2*omega * sin (40 * DEG)                ! representative Coriolis parameter
  beta           = 2*omega * cos (40 * DEG) / radius       ! beta parameter at 45 degrees latitude

  min_depth      = -50 * METRE / H_norm                    ! minimum allowed depth (must be negative)
  k_T            =       1 / (30 * DAY)                    ! relaxation time to mean buoyancy profile (if relax = .true.)
  
  if (zlevels == 1) then
     relax                = .false.
     sigma_z              = .false.
     vert_diffuse         = .false.

     max_depth            = -H_earth
     coords               = "uniform"
     z_mixed              = max_depth                      ! location of top (less dense) layer in two layer case
     z_linear             = max_depth                      ! location of layer forced by surface wind stress
     drho                 =       0 * KG/METRE**3          ! density perturbation at free surface
     tau_0                =     0.4 * NEWTON/METRE**2      ! maximum wind stress
     u_wbc                =       1 * METRE/SECOND         ! estimated western boundary current speed
     
     bottom_friction_case =    rb_0                        ! constant bottom friction
  elseif (zlevels >= 2) then
     relax                = .true.                         ! relax to mean vertical stratification
     remap                = .true.                         ! remap vertical coordinates
     sigma_z              = .true.                         ! sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     tke_closure          = .false.                        ! use analytic profiles for eddy viscosity/diffusivity
     vert_diffuse         = .true.                         ! use vertical diffusion model

     stratification       = "tanh"                         ! type of stratification (tanh, linear)
     z_mixed              =    -200 * METRE                ! bottom of constant density surface mixed
     z_linear             =    -500 * METRE                ! bottom of linear stratification layer below mixed layer
     max_depth            =   -4000 * METRE                ! total depth
     
     bottom_friction_case = rb_0                           ! constant bottom friction equal to NEMO value 4e-4

     drho                 =      -4 * KG/METRE**3          ! density perturbation at free surface at poles
     tau_0                =     0.1 * NEWTON/METRE**2      ! maximum wind stress
     u_wbc                =       1 * METRE/SECOND         ! estimated western boundary current speed (tanh)

     ! Solar flux
     Q_sr                 =       0 * WATT/METRE**2        ! incoming solar radiation heat flux (set to zero to turn off solar forcing)
     R_lw                 =     1.0_dp                     ! proportion of flux in longwave
     xi_lw                =     200 * METRE                ! penetration depth of solar flux
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Characteristic scales
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Area_min       = 4*MATH_PI * radius**2 / number_hex (max_level)
  dx_min         = sqrt (2 / sqrt(3.0_dp) * Area_min)
  dz             = abs (max_depth) / dble (zlevels)                           ! layer depth scale
  
  wave_speed     = sqrt (grav_accel * abs (max_depth))                        ! inertia-gravity wave speed
  dt_init        = cfl_num * 0.85 * dx_min / wave_speed                       ! average time step
  visc           = C_Drake * 2.51/dt_init * (dx_min**2/24/1.65)**Laplace_rotu ! viscosity for RK3
  delta_I        = sqrt (u_wbc/beta)                                          ! inertial layer
  delta_M        = (visc/beta)**(1.0_dp/(2*Laplace_rotu + 1))                 ! Munk layer scale
  delta_sm       = u_wbc / f0                                                 ! barotropic submesoscale
  delta_S        = bottom_friction_case / (abs(max_depth) * beta)             ! Stommel layer scale
  Rey            = u_wbc * delta_sm**(2*Laplace_rotu - 1) / visc              ! Reynolds number of western boundary current
  Ro             = u_wbc / (delta_M*f0)                                       ! Rossby number (based on boundary current)

  ! First baroclinic mode speed for linear stratification
  if (zlevels == 1) then
     c1 = 0.0_dp
  else
     h_linear       = z_mixed - z_linear                                      ! approximate thickness of linear stratification region
     bv             = sqrt (- grav_accel/ref_density * drho/h_linear)         ! Brunt-Vaisala frequency
     if (zlevels == 2) then                                                      
        c1 = sqrt (grav_accel * abs (drho) /ref_density * abs(z_mixed) * (z_mixed - max_depth) / abs (max_depth)) 
     elseif (zlevels >= 3) then                                                  
        c1 = bv * h_linear / MATH_PI                                   
     endif

  end if
  lambda0        = wave_speed / f0                                            ! external scale
  lambda1        = c1 / f0                                                    ! mesoscale

  ! Dimensional scaling
  Hdim           = abs (max_depth)                                            ! vertical length scale
  Ldim           = delta_I                                                    ! length scale
  Udim           = u_wbc                                                      ! velocity scale
  Tdim           = Ldim / Udim                                                ! time scale
  Mudim          = ref_density * dz                                           ! rho_dz scale
  Thetadim       =        drho * dz                                           ! buoyancy scale

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Initialization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (etopo_bathy .or. etopo_coast) call read_etopo_data
  call assign_functions
  call initialize (run_id)
  call print_test_case_parameters
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(a,/)') &
       '----------------------------------------------------- Start simulation run &
       &-------------------------------------------------------'
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  total_cpu_time = 0.0_dp
  do while (time < time_end)
     call start_timing
     call time_step 
     if (relax) call euler (sol, wav_coeff, trend_relax, dt)
     call stop_timing

     call print_log
  end do
  if (rank == 0) write (6,'(a,es11.4)') "Total cpu time = ", total_cpu_time
  call finalize
end program Drake
