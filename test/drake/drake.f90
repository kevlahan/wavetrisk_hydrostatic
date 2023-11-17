program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  use main_mod
  use test_case_mod
  use io_mod
  use vert_diffusion_mod
  implicit none

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Earth parameters
  radius_earth   = 6371.229d0 * KM                      ! radius of Earth
  omega_earth    = 7.29211d-5 * RAD/SECOND              ! rotation rate of Earth
  H_earth        =        4d0 * KM                      ! mean ocean depth of Earth
  g_earth        = 9.80616d0  * METRE/SECOND**2         ! gravitational acceleration 
  ref_density    =     1030d0 * KG/METRE**3             ! reference density at depth (seawater)

  ! Earth scaling factors
  L_norm         = radius_earth
  H_norm         = H_earth
  U_norm         = sqrt (g_earth * H_earth)
  T_norm         = L_norm / U_norm

  ! Use normalized equation
  normalized = .false.

  if (normalized) then
     radius      = radius_earth / L_norm                     
     omega       =  omega_earth * T_norm
     grav_accel  =  g_earth * H_norm / U_norm**2
  else
     radius      = radius_earth / scale                  ! mean radius of the small planet
     omega       = omega_earth  / scale_omega            ! angular velocity (scaled for small planet to keep beta constant)
     grav_accel  = g_earth
  end if

  f0             = 2d0*omega*sin(45d0*DEG)               ! representative Coriolis parameter
  beta           = 2d0*omega*cos(45d0*DEG) / radius      ! beta parameter at 45 degrees latitude

  ! Free surface perturbation parameters
  dH             =   0d0  * METRE / H_norm               ! initial perturbation to the free surface
  pert_radius    =   1d3  * KM    / L_norm               ! radius of Gaussian free surface perturbation
  lon_c          = -50d0  * DEG                          ! longitude location of perturbation
  lat_c          =  25d0  * DEG                          ! latitude  location of perturbation

  min_depth      = -50d0  * METRE / H_norm               ! minimum allowed depth (must be negative)
  
  ! Numerical method parameters
  default_thresholds = .true.

  mode_split         = .true.                            ! split barotropic mode if true
  if (mode_split) then
     cfl_num         = 15d0
     timeint_type    = "RK3"                         
  else
     cfl_num         = 0.3d0                             
     timeint_type    = "RK45"                         
  end if
  
  match_time         = .true.                           ! avoid very small time steps when saving 
  compressible       = .false.                          ! always run with incompressible equations

  remapscalar_type   = "PPR"                            ! remapping scheme for scalars
  remapvelo_type     = "PPR"                            ! remapping scheme for velocity

  nstep_init         = 10                               ! take nstep_init small steps on restart
  adapt_dt           = .false.

  ! Topography (etopo smoothing not yet implemented)
  alpha              = 1d-1
  penalize           = .true.                           ! penalize land regions
  etopo_bathy        = .false.                          ! etopo data for bathymetry
  etopo_coast        = .false.                          ! etopo data for coastlines (i.e. penalization)
  etopo_res          = 4                                ! resolution of etopo or analytic data in arcminutes

  dx_min             = sqrt (4d0/sqrt(3d0) * 4d0*MATH_PI*radius**2/(20d0*4d0**max_level))              
  dx_max             = sqrt (4d0/sqrt(3d0) * 4d0*MATH_PI*radius**2/(20d0*4d0**min_level))

  Laplace_order_init = 2                                ! Laplacian if 1, bi-Laplacian if 2. No diffusion if 0.
  C_visc(S_MASS)     = 0d-3                             ! dimensionless viscosity of S_MASS
  C_visc(S_TEMP)     = 0d-3                             ! dimensionless viscosity of S_TEMP
  C_visc(S_VELO)     = 5d-4                             ! dimensionless viscosity of S_VELO (rotu, divu)

  if (zlevels == 1) then
     sigma_z              = .false.
     vert_diffuse         = .false.
     max_depth            = -H_earth
     coords               = "uniform"
     mixed_layer          = max_depth                     ! location of top (less dense) layer in two layer case
     thermocline          = max_depth                     ! location of layer forced by surface wind stress
     drho                 =     0d0 * KG/METRE**3         ! density perturbation at free surface
     tau_0                =   0.4d0 * NEWTON/METRE**2     ! maximum wind stress
     u_wbc                =     1d0 * METRE/SECOND        ! estimated western boundary current speed
     k_T                  =     0d0                       ! relaxation to mean buoyancy profile
  elseif (zlevels >= 2) then
     coords               = "uniform"
     sigma_z              = .true.                        ! sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     max_depth            =   -4000d0 * METRE             ! total depth, constant density at depth > thermocline
     thermocline          =   -4000d0 * METRE             ! linear stratification region between thermocline and mixed_layer
     mixed_layer          =    -200d0 * METRE             ! constant density at depth < mixed_layer

     remap                = .true.
     iremap               =   5
     
     vert_diffuse         = .true.
     tke_closure          = .false.
     if (zlevels == 2) then
        Kt_const         =  0d0 * METRE**2 / SECOND       ! eddy diffusion
        Kv_bottom         = 4d0 * METRE**2 / SECOND       ! eddy viscosity
     end if
     e_min                =  0d0                          ! minimum TKE
     patankar             = .false.                       ! avoid noise with zero initial velocity
     enhance_diff         = .false.
     
     drho                 =      -2d0 * KG/METRE**3       ! density perturbation at free surface at poles
     tau_0                =     0.1d0 * NEWTON/METRE**2   ! maximum wind stress
     u_wbc                =       1d0 * METRE/SECOND      ! estimated western boundary current speed
     k_T                  =       1d0 / (30d0 * DAY)      ! relaxation to mean buoyancy profile
  end if
  
  ! Characteristic scales
  wave_speed         = sqrt (grav_accel*abs(max_depth))   ! inertia-gravity wave speed
  dt_cfl             = cfl_num * dx_min / wave_speed
  dt_init            = dt_cfl

  visc_sclr(S_MASS)  = C_visc(S_MASS) * dx_min**(2d0*Laplace_order_init)/dt_cfl 
  visc_sclr(S_TEMP)  = C_visc(S_TEMP) * dx_min**(2d0*Laplace_order_init)/dt_cfl 
  visc_divu          = C_visc(S_VELO) * dx_min**(2d0*Laplace_order_init)/dt_cfl 
  visc_rotu          = C_visc(S_VELO) * dx_min**(2d0*Laplace_order_init)/dt_cfl 

  ! Set bottom friction so lateral viscosity dominates
  bottom_friction_case = 0.25d0 * (visc_rotu * beta**(2d0*Laplace_order_init))**(1d0/(2d0*Laplace_order_init+1d0))

  Rd             = wave_speed / f0                                    ! barotropic Rossby radius of deformation             
  drho_dz        = drho / (mixed_layer-thermocline)                   ! density gradient
  bv             = sqrt (grav_accel * abs(drho_dz)/ref_density)       ! Brunt-Vaisala frequency
  delta_I        = sqrt (u_wbc/beta)                                  ! inertial layer
  delta_M        = (visc_rotu/beta)**(1d0/(2d0*Laplace_order_init+1)) ! Munk layer scale
  delta_sm       = u_wbc / f0                                         ! barotropic submesoscale
  delta_S        = bottom_friction_case / (abs(max_depth) * beta)     ! Stommel layer scale
  Fr             = u_wbc / (bv*abs(max_depth))                        ! Froude number
  Rey            = u_wbc * delta_sm**(2d0*Laplace_order_init-1d0) / visc_rotu ! Reynolds number of western boundary current
  Ro             = u_wbc / (delta_M*f0)                               ! Rossby number (based on boundary current)

  ! Baroclinic wave speed
  if (zlevels == 2) then
     c1 = sqrt (grav_accel*abs(drho)/ref_density*mixed_layer*(max_depth-mixed_layer)/abs(max_depth)) ! two-layer internal wave speed
  elseif (zlevels >= 3) then
     c1 = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification
  endif
  lambda0        = wave_speed / f0                              ! external scale
  lambda1        = c1 / f0                                      ! mesoscale

  ! First baroclinic Rossby radius of deformation
  if (zlevels == 1) then
     Rb = 0d0
  elseif (zlevels == 2) then
     Rb = c1 / f0                                 
  else
     Rb = bv * abs(max_depth) / (MATH_PI*f0)
  end if
  
  ! Dimensional scaling
  Udim           = u_wbc                              ! velocity scale
  Ldim           = delta_I                            ! length scale 
  Tdim           = Ldim/Udim                          ! time scale
  Hdim           = abs (max_depth)                    ! vertical length scale
  
  if (etopo_bathy .or. etopo_coast) call read_etopo_data

  s_test = radius/4d0  ! parameter for elliptic solver test
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)

  ! Initialize random numbers
  call random_seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Set interval for adapting grid based on the horizontal advective velocity scale (i.e. advect no more than one grid point before adapting)
  iadapt     = 1        ! Drake unstable with trend computation over entire grid if iadapt > 1
  irebalance = 4*iadapt ! rebalance interval using charm++/AMPI

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)

  elliptic_solver => SRJ
  
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
end program

