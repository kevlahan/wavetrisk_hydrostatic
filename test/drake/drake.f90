program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  use main_mod
  use test_case_mod
  use io_vtk_mod
  use vert_diffusion_mod
  implicit none
  real(8) :: dz, visc

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  call read_test_case_parameters

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Numerical method parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  default_thresholds      = .true.
  scale_aware             = .true.                    
  mode_split              = .true.                     
  split_mean_perturbation = .true.
  adapt_dt                = .false.
  nstep_init              = 10
  if (mode_split) then
     cfl_num              = 20.0_dp
     timeint_type         = "RK3"                         
  else
     cfl_num              = 0.3_dp                             
     timeint_type         = "RK45"                         
  end if
  match_time              = .true.                 
  compressible            = .false.                
  log_min_mass            = .false.                
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Earth parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  radius_earth   =   6371.229 * KM                      
  omega_earth    = 7.29211e-5 * RAD/SECOND              
  H_earth        =          4 * KM                     
  g_earth        =    9.80616 * METRE/SECOND**2        
  ref_density    =       1030 * KG/METRE**3            

  ! Earth scaling factors
  L_norm         = radius_earth
  H_norm         = H_earth
  U_norm         = sqrt (g_earth * H_earth)
  T_norm         = L_norm / U_norm

  ! Use normalized equation
  normalized = .false.
  if (normalized) then
     radius      = radius_earth / L_norm                     
     omega       = omega_earth * T_norm
     grav_accel  = g_earth * H_norm / U_norm**2
  else
     radius      = radius_earth / scale                    ! mean radius of the small planet
     omega       = omega_earth  / scale_omega              ! angular velocity (scaled for small planet to keep beta constant)
     grav_accel  = g_earth
  end if

  f0             = 2*omega * sin (40 * DEG)                ! representative Coriolis parameter
  beta           = 2*omega * cos (40 * DEG) / radius       ! beta parameter at 45 degrees latitude

  min_depth      = -50 * METRE / H_norm                    ! minimum allowed depth (must be negative)
  
  ! Topography (etopo smoothing not yet implemented)
  alpha          = 1e-1_dp
  penalize       = .true.                                  ! penalize land regions
  etopo_bathy    = .false.                                 ! etopo data for bathymetry
  etopo_coast    = .false.                                 ! etopo data for coastlines (i.e. penalization)
  etopo_res      = 4                                       ! resolution of etopo or analytic data in arcminutes
 
  if (zlevels == 1) then
     sigma_z              = .false.
     vert_diffuse         = .false.

     max_depth            = -H_earth
     coords               = "uniform"
     mixed_layer          = max_depth                      ! location of top (less dense) layer in two layer case
     thermocline          = max_depth                      ! location of layer forced by surface wind stress
     drho                 =       0 * KG/METRE**3          ! density perturbation at free surface
     tau_0                =     0.4 * NEWTON/METRE**2      ! maximum wind stress
     u_wbc                =       1 * METRE/SECOND         ! estimated western boundary current speed
     
     bottom_friction_case =    rb_0                        ! constant bottom friction
     k_T                  =     0.0_dp                     ! relaxation to mean buoyancy profile
  elseif (zlevels >= 2) then
     coords               = "uniform"
     sigma_z              = .true.                         ! sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     max_depth            =   -4000 * METRE                ! total depth
     thermocline          =   -4000 * METRE                ! linear stratification region between thermocline and mixed_layer
     mixed_layer          =    -200 * METRE                ! constant density at depth < mixed_layer

     remap                = .true.
     min_mass_remap       = 0.7_dp

     bottom_friction_case = rb_0                           ! constant bottom friction equal to NEMO value 4e-4
     vert_diffuse         = .true.
     tke_closure          = .false.

     drho                 =      -4 * KG/METRE**3          ! density perturbation at free surface at poles
     tau_0                =     0.1 * NEWTON/METRE**2      ! maximum wind stress
     u_wbc                =       1 * METRE/SECOND         ! estimated Western boundary current speed
     k_T                  =       1 / (30 * DAY)           ! relaxation to mean buoyancy profile

     Kt_max               = 5e-3_dp
     Kv_max               = 5e-3_dp
  end if

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Characteristic scales
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Area_min       = hex_area_avg (max_level)
  dx_min         = dx_avg (max_level)
  wave_speed     = sqrt (grav_accel * abs(max_depth))                 ! inertia-gravity wave speed
  dt_init        = cfl_num * 0.85 * dx_min / (wave_speed + u_wbc)     ! average time step
  visc           = C_Drake * Area_min**Laplace_rotu / dt_init         ! viscosity
  Rd             = wave_speed / f0                                    ! barotropic Rossby radius of deformation             
  drho_dz        = drho / (mixed_layer-thermocline)                   ! density gradient
  bv             = sqrt (grav_accel * abs(drho_dz)/ref_density)       ! Brunt-Vaisala frequency
  delta_I        = sqrt (u_wbc/beta)                                  ! inertial layer
  delta_M        = (visc/beta)**(1d0/(2d0*Laplace_rotu+1))            ! Munk layer scale
  delta_sm       = u_wbc / f0                                         ! barotropic submesoscale
  delta_S        = bottom_friction_case / (abs(max_depth) * beta)     ! Stommel layer scale
  Fr             = u_wbc / (bv*abs(max_depth))                        ! Froude number
  Rey            = u_wbc * delta_sm**(2*Laplace_rotu-1) / visc        ! Reynolds number of western boundary current
  Ro             = u_wbc / (delta_M*f0)                               ! Rossby number (based on boundary current)

  ! Baroclinic wave speed
  if (zlevels == 2) then     ! two-layer internal wave speed
     c1 = sqrt (grav_accel * abs(drho) /ref_density * mixed_layer * (max_depth-mixed_layer) / abs(max_depth)) 
  elseif (zlevels >= 3) then ! first baroclinic mode speed for linear stratification
     c1 = bv * sqrt (abs(max_depth) / grav_accel) / MATH_PI * wave_speed                                      
  endif
  lambda0        = wave_speed / f0                                    ! external scale
  lambda1        = c1 / f0                                            ! mesoscale
 
  ! First baroclinic Rossby radius of deformation
  if (zlevels == 1) then
     Rb = 0d0
  elseif (zlevels == 2) then
     Rb = c1 / f0                                 
  else
     Rb = bv * abs(max_depth) / (MATH_PI*f0)
  end if

  dz = abs (max_depth) / dble(zlevels) ! layer depth scale

  ! Dimensional scaling
  Ldim           = delta_I             ! length scale 
  Hdim           = abs (max_depth)     ! vertical length scale
  Tdim           = Ldim / Udim         ! time scale
  
  Mudim          = ref_density * dz    ! rho_dz scale
  Thetadim       =        drho * dz    ! buoyancy scale
  Udim           = u_wbc               ! velocity scale

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Initialization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (etopo_bathy .or. etopo_coast) call read_etopo_data
  call assign_functions
  call initialize (run_id)
  call print_test_case_parameters
  call write_and_export (iwrite) ! save initial conditions
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(a,/)') &
       '----------------------------------------------------- Start simulation run &
       -------------------------------------------------------'
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  total_cpu_time = 0d0
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
     if (k_T /= 0d0) call euler (sol, wav_coeff, trend_relax, dt)
     call stop_timing

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)
        call write_and_export (iwrite) ! save fields
     end if
  end do

  if (rank == 0)  write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  call finalize
end program

