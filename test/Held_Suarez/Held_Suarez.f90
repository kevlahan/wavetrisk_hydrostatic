program Held_Suarez
  ! Held & Suarez (1994) test case
  ! Bulletin of the American Meteorological Society 75 (10), 1825-1830
  use main_mod
  use ops_mod
  use test_case_mod
  use io_mod
  use topo_grid_descriptor_mod
  implicit none

  integer        :: l
  real(8)        :: fine_mass
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
  ! Reference pressure (mean surface pressure) in Pascals
  if (NCAR_topo) then
     call std_surf_pres (0d0, p_0)
  else
     p_0         = 1000d0  * hPa                
  end if
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
  sigma_c        = 1d0 - sigma_b
  
  ! Local test case parameters (Jablonowski and Williamson 2006 zonally symmetric initial conditions)
  u_0            = 70d0     * METRE/SECOND          ! maximum velocity of zonal wind
  gamma_T        = 5d-3     * KELVIN/METRE          ! temperature lapse rate
  delta_T2       = 4.8d5    * KELVIN                ! empirical temperature difference
  
  sigma_0        = 0.252d0                          ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2d0                            ! value of sigma at the tropopause

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                              ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 70d0     * KELVIN                ! temperature scale for tolerances
  Pdim           = p_0                              ! pressure scale
  dPdim          = 80d0     * hPa                   ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d * Tempdim) / Pdim           ! specific volume scale
  wave_speed     = sqrt (gamma * Pdim * specvoldim) ! acoustic wave speed

  Udim           = u_0                              ! velocity scale
  Tdim           = 1d0  * DAY                       ! time scale
  Ldim           = Udim * Tdim                      ! length scale
  Hdim           = wave_speed**2 / grav_accel       ! vertical length scale

  ! Numerical method parameters
  Area_min           = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**max_level)
  Area_max           = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**min_level)

  dx_min             = sqrt (2d0 / sqrt(3d0) * Area_min)              
  dx_max             = sqrt (2d0 / sqrt(3d0) * Area_max)

  timeint_type       = "RK4"
  iremap             = 4

  default_thresholds = .false.
  compressible       = .true.
  remap              = .true.
  uniform            = .false.

  adapt_dt           = .false.
  if (adapt_dt) then
     cfl_min = 0.1d0                                           ! minimum cfl number
     cfl_max = 1d0                                             ! maximum cfl number
     T_cfl   = 1d0 * DAY                                       ! time over which to increase cfl number from cfl_min to cfl_max
     cfl_num = cfl_min                                         ! initialize cfl number
     dt_init = cfl_min * 0.85d0 * dx_min / (wave_speed + Udim) ! initial time step     (0.85 factor corrects for minimum dx)
     dt_max  = cfl_max * 0.85d0 * dx_min / (wave_speed + Udim) ! equilibrium time step (0.85 factor corrects for minimum dx)
  else
     dt_init = 300d0 * SECOND * 2**(6 - max_level)             ! Lauritzen value
     dt_max  = dt_init
  end if

  ! Diffusion parameters: always use hyperdiffusion
  Laplace_order_init = 2 ! Laplacian if 1, bi-Laplacian if 2. No diffusion if 0

  ! Use scale-aware viscosity (if .false. viscosity depends only on dt)
  scale_aware        = .false.                      

  ! CAM values for viscosity
  nu_sclr            = 1.0d15
  nu_rotu            = 1.0d15
  nu_divu            = 2.5d15

  ! Equivalent non-dimensional viscosities
  C_visc(S_MASS)     = nu_sclr * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_visc(S_TEMP)     = nu_sclr * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_visc(S_VELO)     = nu_rotu * dt_max / (3d0 * Area_min**Laplace_order_init) 
  C_div              = nu_divu * dt_max / (3d0 * Area_min**Laplace_order_init) 
  
  ! Adapt on mean variables (fluctuations are initially zero)
  log_min_mass       = .false.
  log_total_mass     = .false.

  ! Topography data levels (need to use same DOMAIN_LEVEL as used to generate topography data)
  topo_min_level     = 4
  topo_max_level     = 6

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)
  call print_test_case_parameters

  ! Save initial conditions
  call omega_velocity
  call write_and_export (iwrite)
  
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0d0; time_start = time

  do while (time < time_end)
     cfl_num = cfl (time) ! gradually increase cfl number
     call start_timing
     call time_step (dt_write, aligned)
     !if (time >= 200*DAY .and. modulo (istep, 100) == 0) call statistics
     call euler (sol, wav_coeff, trend_cooling, dt)
     call stop_timing
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        
        if (remap) call remap_vertical_coordinates

        ! Save fields
        call omega_velocity
        call write_and_export (iwrite)
        
        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call write_checkpoint (run_id, rebalance) 

           ! Save statistics
           if (time >= 200*DAY .and. modulo (istep, 100) == 0) then
              call combine_stats
              if (rank == 0) call write_out_stats
           end if
        end if
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(a,eS11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program Held_Suarez
