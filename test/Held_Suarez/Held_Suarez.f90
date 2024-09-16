program Held_Suarez
  ! Held & Suarez (1994) test case
  ! Bulletin of the American Meteorological Society 75 (10), 1825-1830
  use main_mod
  use test_case_mod
  use lnorms_mod
  use physics_Held_Suarez_mod
  implicit none

  integer        :: l
  real(8)        :: area_sphere, dx_scaling, nu_CAM, nu, nu_scaling, res_scaling, rx0_max, rx1_max
  logical        :: aligned
  character(256) :: input_file

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Initialize random number generator
  call initialize_seed

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371d0    * KM                   ! mean radius of the Earth
  grav_accel     = 9.8d0     * METRE/SECOND**2      ! gravitational acceleration 
  omega          = 7.292d-5  * RAD/SECOND           ! Earth's angular velocity in radians per second
  ref_density    = ref_density_air                  ! reference density of air
  call std_surf_pres (0d0, p_0)                     ! surface pressure from USA standard atmosphere model
  c_p            = 1004d0  * JOULE/(KG*KELVIN)      ! specific heat at constant pressure in joules per kilogram Kelvin
  kappa          = 2d0/7d0                          ! kappa
  R_d            = kappa * c_p                      ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_v            = c_p - R_d                        ! specific heat at constant volume c_v = c_p - R_d
  gamma          = c_p / c_v                        ! heat capacity ratio

  ! Local test case parameters
  T_0            = 300d0      * KELVIN              ! reference temperature
  T_mean         = 315d0      * KELVIN              ! mean temperature
  T_tropo        = 200d0      * KELVIN              ! tropopause temperature
  k_a            = 1d0/40d0   / DAY                 ! cooling at free surface of atmosphere
  k_f            = 1d0        / DAY                 ! Rayleigh friction
  k_s            = 1d0/4d0    / DAY                 ! cooling at surface
  delta_T        = 65d0       * KELVIN/METRE        ! meridional temperature gradient
  delta_theta    = 10d0       * KELVIN/METRE        ! vertical temperature gradient
  sigma_b        = 0.7d0                            ! normalized tropopause pressure height
  
  ! Local test case parameters (Jablonowski and Williamson 2006 zonally symmetric initial conditions)
  u_0            = 70d0     * METRE/SECOND          ! maximum velocity of zonal wind
  gamma_T        = 5d-3     * KELVIN/METRE          ! temperature lapse rate
  delta_T2       = 4.8d5    * KELVIN                ! empirical temperature difference
  
  sigma_0        = 0.252d0                          ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2d0                            ! value of sigma at the tropopause

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                              ! temperature scale (both theta and T from DYNAMICO)
  Pdim           = p_0                              ! pressure scale

  ! Dimensional scaling
  specvoldim     = (R_d * Tempdim) / Pdim           ! specific volume scale
  wave_speed     = sqrt (gamma * Pdim * specvoldim) ! acoustic wave speed

  Udim           = u_0                              ! velocity scale
  Tdim           = 1d0  * DAY                       ! time scale
  Ldim           = Udim * Tdim                      ! length scale
  Hdim           = wave_speed**2 / grav_accel       ! vertical length scale

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Numerical method parameters
  compressible             = .true.                 ! compressible equations
  split_mean_perturbation  = .false.                ! split prognostic variables into mean and fluctuations
  uniform                  = .false.                ! hybrid vertical grid (based on A, B coefficients)
  adapt_dt                 = .false.                ! fixed time step
  default_thresholds       = .false.                ! thresholding type
  remap                    = .true.                 ! use vertical remapping
  iremap                   = 4                      ! remap every 4 dt
  timeint_type             = "RK3"                  ! time integration scheme (use RK34, RK45 or RK4)
  Laplace_order_init       = 2                      ! bi-Laplacian horizontal diffusion
  scale_aware              = .false.                ! do not use scale-aware viscosity
  analytic_topo            = "none"                 ! type of analytic topography (mountains or none if NCAR_topo = .false.)
  log_min_mass             = .false.                ! compute minimum mass at each dt (for checking stability issues)
  log_total_mass           = .false.                ! check whether total mass is conserved (for debugging)

  ! Average hexagon areas and horizontal resolution
  area_sphere = 4d0*MATH_PI * radius**2 
  Area_min    = area_sphere / (10d0 * 4d0**max_level)
  Area_max    = area_sphere / (10d0 * 4d0**min_level)

  res_scaling = (120d0*KM) / (100d0*KM)                        ! ratio between TRiSK  and CAM dx, approximately (120/100)^4
  dx_min      = sqrt (2d0 / sqrt(3d0) * Area_min)              
  dx_max      = sqrt (2d0 / sqrt(3d0) * Area_max)
  dx_scaling  = 2d0 ** (dble (6 - max_level))                  ! scaling factor compared to approximately J6 base CAM value

  ! Time adaptivity parameters
  if (adapt_dt) then
     cfl_max = 1d0                                             ! maximum cfl number
     cfl_min = cfl_max                                         ! minimum cfl number
     T_cfl   = 5d-1 * DAY                                      ! time over which to increase cfl number from cfl_min to cfl_max
     cfl_num = cfl_min                                         ! initialize cfl number
     dt_init = cfl_min * 0.85d0 * dx_min / (wave_speed + Udim) ! initial time step     (0.85 factor corrects for minimum dx)
     dt_max  = cfl_max * 0.85d0 * dx_min / (wave_speed + Udim) ! equilibrium time step (0.85 factor corrects for minimum dx)
  else
     dt_init = 300d0 * SECOND * dx_scaling                     ! Lauritzen value
     dt_max  = dt_init
  end if

  ! CAM values for viscosity for 1 degree are 1e15 except 2.5e15 for divu
  ! note that wavetrisk J6 resolution is about 120 km while 1 degree CAM resolution is about 110 km
  !
  ! Rescale by difference in grid sizes to get equivalent viscosity
  nu_scaling         = (res_scaling * dx_scaling)**(2*Laplace_order_init)
  nu_CAM             = 1d15 * METRE**4/SECOND
  nu                 = nu_CAM * nu_scaling
  
  nu_sclr            = nu         ! CAM value (scaled)
  nu_rotu            = nu / 4d0   ! smaller to respect smaller stability limit for rotu Laplacian
  nu_divu            = nu * 2.5d0 ! increase by CAM ratio

  ! Equivalent non-dimensional viscosities
  C_visc(S_MASS)     = nu_sclr * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_visc(S_TEMP)     = nu_sclr * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_visc(S_VELO)     = nu_rotu * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_div              = nu_divu * dt_max / (3d0 * Area_min**Laplace_order_init)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)
  call print_test_case_parameters

  ! Save initial conditions
  call omega_velocity
  call write_and_export (iwrite)

  ! Compute hydrostatic error factors for topography
  if (NCAR_topo .or. analytic_topo=="mountains") then
     if (rank == 0) write (6,'(/,a)') "Level      rx0_max   rx1_max"
     do l = min_level, max_level
        call cal_rx0_max (l, rx0_max)
        call cal_rx1_max (l, rx1_max)
        if (rank == 0) write (6,'(i2,8x,2(es8.2,3x))') l, rx0_max, rx1_max
     end do
  end if
   
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  
  total_cpu_time = 0d0; time_start = time
  do while (time < time_end)
     cfl_num = cfl (time) ! gradually increase cfl number
     
     call start_timing
     call time_step (dt_write, aligned)
     call euler (sol, wav_coeff, trend_physics_Held_Suarez, dt)
     call stop_timing
     
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)

        ! Save fields (after reloading checkpoint)
        call omega_velocity
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program Held_Suarez
