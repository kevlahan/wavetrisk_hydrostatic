program climate
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !    Climate simulation using Held and Suarez (1994) or Simple Physics (Hourdin 1993) subgrid scale model
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use main_mod
  use test_case_mod
  use io_vtk_mod
  implicit none
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  call read_test_case_parameters  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Numerical method parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  timeint_type             = "RK3"                            ! time integration scheme 
  adapt_dt                 = .false.                          ! adapt time step
  compressible             = .true.                           ! compressible equations
  default_thresholds       = .false.                          ! thresholding type
  log_min_mass             = .true.                           ! compute minimum mass at each dt (for checking stability issues)
  log_total_mass           = .false.                          ! check whether total mass is conserved (for debugging)
  remap                    = .true.                           ! use vertical remapping
  split_mean_perturbation  = .true.                           ! split prognostic variables into mean and fluctuations
  uniform                  = .false.                          ! hybrid vertical grid (based on A, B coefficients)

  Laplace_sclr             = 2                                ! scalars
  Laplace_divu             = 2                                ! div u
  Laplace_rotu             = 2                                ! rot u 
  min_mass_remap           = 0.7d0                            ! minimum mass at which to remap
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Local test case parameters (default values for many parameters set in physics module)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  physics_model            = .true.                           ! use physics model sub-step (type is determined in input)

  ! Simple physics sub-models
  convecAdj_model          = .true.                           ! convective adjustment module
  diurnal                  = .true.                           ! diurnal cycle 
  radiation_model          = .true.                           ! radiation module
  turbulence_model         = .true.                           ! vertical diffusion module
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Standard (shared) paraggmeter values for the simulation
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  c_p                      = 1004     * JOULE/(KG*KELVIN)     ! specific heat at constant pressure in joules per kilogram Kelvin
  grav_accel               = 9.8      * METRE/SECOND**2       ! gravitational acceleration 
  omega                    = 7.292d-5 * RAD/SECOND            ! Earth's angular velocity in radians per second
  R_d                      = 287      * JOULE/(KG*KELVIN)     ! set to a whole number
  radius                   = 6371     * KM                    ! mean radius of the Earth
  ref_density              = ref_density_air                  ! reference density of air
  T_0                      = 285      * KELVIN                ! reference temperature (simple physics)
  u_0                      =  30      * METRE/SECOND          ! geostrophic velocity
  
  ! Derived quantities
  c_v                      = c_p - R_d                        ! specific heat at constant volume c_v = c_p - R_d
  gamma                    = c_p / c_v                        ! heat capacity ratio
  kappa                    = R_d / c_p                        ! kappa

  wave_speed               = sqrt (gamma * (R_d * T_0) )      ! acoustic wave speed
  max_depth                = wave_speed**2 / grav_accel       ! depth of atmosphere
  dz                       = max_depth / dble (zlevels)       ! representative layer height
  
  call std_surf_pres (0d0, p_0)                               ! reference pressure (USA standard atmosphere model)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Dimensional scaling
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  Tdim                     = 1  * DAY                         ! time scale
  Ldim                     = u_0 * Tdim                       ! length scale
  Hdim                     = max_depth                        ! vertical length scale
  Pdim                     = p_0                              ! pressure scale
  Tempdim                  = T_0                              ! temperature scale (both theta and T from DYNAMICO)
  
  Mudim                    = ref_density * dz                 ! mu scale
  Thetadim                 = ref_density * dz * Tempdim       ! Theta scale
  Udim                     = u_0                              ! velocity scale 

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Initialization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call assign_functions
  call initialize_seed
  call initialize (run_id)
  call print_test_case_parameters
  
  open (unit = 12, file = trim(run_id)//'_log', action = 'WRITE', form = 'FORMATTED', position = 'APPEND')
  !call write_and_export (iwrite)
  
  total_cpu_time = 0d0; time_start = time
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Run simulation
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  if (rank == 0) write (6,'(a,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
 
  do while (time < time_end)
     call start_timing ; call time_step (dt_write, aligned) ; call stop_timing
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)
        call write_and_export (iwrite) 
     end if
  end do
  close (12)
  if (rank == 0) write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
  call finalize
end program climate
