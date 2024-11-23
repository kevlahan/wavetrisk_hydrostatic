program climate
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    Climate simulation using Held and Suarez (1994) or Simple Physics (Hourdin 1993) subgrid scale model
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use main_mod
  use test_case_mod
  use lnorms_mod
  implicit none

  integer        :: l
  real(8)        :: rx0_max, rx1_max
  logical        :: aligned
  character(256) :: input_file

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  
   ! Read test case parameters
  call read_test_case_parameters
  
  ! Initialize random number generator 
  call initialize_seed

  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Numerical method parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  scale_aware              = .true.                            ! scale-aware viscosity
  adapt_dt                 = .true.                            ! adapt time step
  default_thresholds       = .true.                            ! thresholding type
  cfl_num                  = 1d0                               ! CFL number

  dt_phys                  = 0 * MINUTE                        ! interval for physics split step
  remap                    = .true.                            ! use vertical remapping
  min_mass_remap           = 0.8d0                             ! minimum mass at which to remap
 
  split_mean_perturbation  = .false.                           ! split prognostic variables into mean and fluctuations
  compressible             = .true.                            ! compressible equations
  uniform                  = .false.                           ! hybrid vertical grid (based on A, B coefficients)
  timeint_type             = "RK3"                             ! time integration scheme (use RK34, RK45, RK3 or RK4)
  Laplace_order_init       = 2                                 ! bi-Laplacian horizontal diffusion
  analytic_topo            = "none"                            ! type of analytic topography (mountains or none if NCAR_topo = .false.)

  log_min_mass             = .true.                            ! compute minimum mass at each dt (for checking stability issues)
  log_total_mass           = .false.                           ! check whether total mass is conserved (for debugging)


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Local test case parameters (default values for many parameters set in physics module)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  physics_model            = .true.                            ! use physics model sub-step (type is determined in input)
  Ekman_ic                 = .true.                            ! Ekman (T) or zero (F) velocity initial conditions
  u_0                      =  30d0     * METRE/SECOND          ! geostrophic velocity
  T_0                      = 285d0     * KELVIN                ! reference temperature (simple physics)

  ! Simple physics submodel switches
  turbulence_model         = .true.                            ! vertical diffusion module
  radiation_model          = .true.                            ! radiation module
  convecAdj_model          = .true.                            ! convective adjustment module
  diurnal                  = .true.                            ! diurnal cycle 

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Standard (shared) parameter values for the simulation
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  radius                   = 6371d0    * KM                    ! mean radius of the Earth
  grav_accel               = 9.8d0     * METRE/SECOND**2       ! gravitational acceleration 
  omega                    = 7.292d-5  * RAD/SECOND            ! Earth's angular velocity in radians per second
  ref_density              = ref_density_air                   ! reference density of air
  p_top                    = 0.01d0    * Pa                    ! pressure at the top in Pascals
  call std_surf_pres (0d0, p_0)                                ! surface pressure from USA standard atmosphere model
  c_p                      = 1004d0    * JOULE/(KG*KELVIN)     ! specific heat at constant pressure in joules per kilogram Kelvin
  R_d                      = 287d0     * JOULE/(KG*KELVIN)     ! set to a whole number
  ref_density              = 1.204d0   * KG/METRE**3           ! reference density 
  kappa                    = R_d / c_p                         ! kappa
  c_v                      = c_p - R_d                         ! specific heat at constant volume c_v = c_p - R_d
  gamma                    = c_p / c_v                         ! heat capacity ratio
  wave_speed               = sqrt (gamma * (R_d * T_0) )       ! acoustic wave speed
  max_depth                = wave_speed**2 / grav_accel        ! depth of atmosphere
  dz                       = max_depth / dble (zlevels)        ! representative layer height

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Dimensional scaling
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Tdim                     = 1d0  * DAY                         ! time scale
  Ldim                     = u_0 * Tdim                         ! length scale
  Hdim                     = max_depth                          ! vertical length scale
  Pdim                     = p_0                                ! pressure scale
  Tempdim                  = T_0                                ! temperature scale (both theta and T from DYNAMICO)
  
  Mudim                    = ref_density * dz                   ! mu scale
  Thetadim                 = ref_density * dz * Tempdim         ! Theta scale
  Udim                     = u_0                                ! velocity scale
  

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Initialization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  call assign_functions

  ! Initialize variables
  call initialize (run_id)

  call print_test_case_parameters

  ! Save initial conditions
  call omega_velocity
  !call write_and_export (iwrite)
  !if (physics_type == "Simple") call mean_values (0) ! processing for the physics package mean values

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
     call start_timing
     call time_step (dt_write, aligned) 
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
!        if (physics_type == "Simple") call mean_values (iwrite) 
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program climate
