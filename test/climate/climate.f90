program climate
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Climate simulation using Held and Suarez (1994) or Simple Physics (Hourdin 1993) subgrid scale model
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use audit316_mod, only: audit316_init,audit316_finish,audit316_begin,audit316_end,audit316_i
  use boundary310_mod, only: b310_init,b310_begin,b310_end,b310_flush
  use parallel_block_mpi_mod, only: snapshot_groups312,audit_transport_footer322
  use runtime322_mod
  use campaign312_mod
  use mpi_f08, only: MPI_Wtime
  use shared_mod, only: dt,n_active,tol,min_level,max_level
  use comm_mpi_mod, only: get_timing
  use main_mod, only: checkpoint312
  use kind_mod,   only : dp
  use shared_mod, only : DAY, JOULE, KELVIN, KG, KM, METRE, N_VARIABLE, RAD, SECOND, &
       Hdim, Ldim, Mudim, Pdim, Tdim, Tempdim, Thetadim, Udim, &
       adapt_dt, c_g, c_p, c_s, c_v, compressible, default_thresholds, gamma, grav_accel, &
       kappa, level_start, log_min_mass, max_depth, omega, p_0, physics_model, physics_type, R_d, radius, &
       ref_density, ref_density_air, run_id, time, time_end, uniform, vtk_grid, wave_speed, zlevels

  use adapt_mod,           only : WT_after_step 
  use arch_mod,            only : finalize, init_arch_mod, rank
  use comm_mpi_mod,        only : init_comm_mpi_mod, start_timing, stop_timing, write_load_conn
  use domain_mod,          only : sol, trend, wav_coeff
  use init_physics_mod,    only : convecAdj_model, diurnal, obliquity, radiation_model, soil_model, turbulence_model
  use checkpoint_mod,      only : dump_adapt_mpi 
  use io_vtk_mod,          only : write_and_export
  use main_mod,            only : initialize, time_step
  use parallel_block_mpi_mod, only : report_parallel_block_profile
  use multi_level_mod,     only : trend_ml
  use std_atm_profile_mod, only : std_surf_pres

  use test_case_mod,       only : assign_functions, dz, initialize_seed, print_log, print_test_case_parameters,  &
       read_test_case_parameters, T_0, topo_test, total_cpu_time, u_0

  implicit none
  real(dp) :: init_start312,startup312
  integer :: j312

  ! Initialize mpi, shared variables and domains
  call init_arch_mod(defer_domain_limit=.true.)
  call init_runtime322()
  call audit316_init(rank)
  call init_comm_mpi_mod
  call read_test_case_parameters   

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Numerical method parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Laplace_divu = 1
  adapt_dt                 = .false.                          ! adapt time step
  if (physics_type == "Held_Suarez") adapt_dt = .false.
  compressible             = .true.                           ! compressible equations
  default_thresholds       = .false.                          ! thresholding type
  log_min_mass             = .true.                           ! compute minimum mass at each dt (for checking stability issues)
  topo_test                = .false.                          ! no physics model stationary flow test
  uniform                  = .false.                          ! hybrid vertical grid (based on A, B coefficients)
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Local test case parameters (default values for many parameters set in physics module)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (.not. topo_test) physics_model = .true.                 ! use physics model sub-step (type is determined in input)

  ! Simple physics sub-models
  obliquity                = 23.5_dp                          ! seasons
  convecAdj_model          = .true.                           ! convective adjustment module
  diurnal                  = .true.                           ! diurnal cycle
  soil_model               = .true.                           ! include soil model
  radiation_model          = .true.                           ! radiation module
  turbulence_model         = .true.                           ! vertical diffusion module
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Standard (shared) parameters
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  c_p                      = 1004     * JOULE/(KG*KELVIN)     ! specific heat at constant pressure
  grav_accel               = 9.8      * METRE/SECOND**2       ! gravitational acceleration 
  omega                    = 7.292e-5 * RAD/SECOND            ! Earth's angular velocity
  R_d                      = 287      * JOULE/(KG*KELVIN)     ! set to a whole number
  radius                   = 6371     * KM                    ! mean radius of the Earth
  ref_density              = ref_density_air                  ! reference density of air
  T_0                      = 285      * KELVIN                ! reference temperature (simple physics)
  u_0                      =  30      * METRE/SECOND          ! geostrophic velocity
  
  ! Derived quantities
  c_v                      = c_p - R_d                        ! specific heat at constant volume c_v = c_p - R_d
  gamma                    = c_p / c_v                        ! heat capacity ratio
  kappa                    = R_d / c_p                        ! kappa

  max_depth                = R_d * T_0 / grav_accel           ! scale height for dry air
  dz                       = max_depth / dble (zlevels)       ! representative layer height

  c_s                      = sqrt (gamma * (R_d * T_0) )      ! acoustic wave speed
  c_g                      = sqrt (grav_accel * max_depth)    ! gravity wave speed
  wave_speed               = max (c_g, c_s)                   ! wave speed used for CFL number 
  
  call std_surf_pres (0.0_dp, p_0)                            ! reference pressure (USA standard atmosphere model)
  
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
  call configure312
  ! Stage321 permits the bounded four-step single-level J5 control only.
  if(tol<=0.0_dp.and.steps312>1) then
     if(tol<0.0_dp.or.min_level/=5.or.max_level/=5) &
          error stop "Stage312 tol=0 allows one timestep except single-level J5"
  end if
  init_start312=MPI_Wtime()
  call audit316_begin("initialization",0_int64)
  call initialize_seed
  call initialize (run_id)
  call audit316_end()
  startup312=MPI_Wtime()-init_start312
  call print_test_case_parameters
  call checkpoint312(9900)
  call snapshot_groups312("initial")

  ! Save initial tendencies in prognostic variables when testing topography
  if (topo_test .and. time_end < 1e-16_dp) then
     if (rank == 0) write (6,'(a)') "Saving trend data"

     call trend_ml (sol(1:N_VARIABLE,1:zlevels), trend(1:N_VARIABLE,1:zlevels))
     sol = trend
     call write_and_export (vtk_grid)
     
     ! Save checkpoint with trend data for spectrum computation
     call WT_after_step (sol, wav_coeff, level_start-1)
     call write_load_conn (0)
     call dump_adapt_mpi  (0)
  end if
  
  total_cpu_time = 0.0_dp
  open (unit = 12, file = trim(run_id)//'_log', action = 'WRITE', form = 'FORMATTED', position = 'APPEND')
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Run simulation
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
  if (rank == 0) write (6,'(a)') &
       '----------------------------------------------------- Start simulation run &
       & ------------------------------------------------------'
  call b310_init
  call observe_runtime322(1)
  do j312=1,steps312
     step312=j312
     if(time>=time_end) error stop "Stage312 time_end reached before step limit"
     if(j312==2) call start_window322()
     call audit316_begin("step",int(j312,int64))
     call start_timing
     call b310_begin
     call time_step
     call b310_end
     call stop_timing
     call audit316_end()
     timing312(j312)=get_timing()
     physical312(j312)=time
     dt312(j312)=dt
     active312(j312)=int(sum(n_active),int64)
     call print_log
     if(j312==4) call end_window322()
     call observe_runtime322(j312+1)
  end do
  ! All numerical steps and their outer timing/logging have completed. In lean
  ! mode this is the sole flush; output events inside time_step retain records.
  call report_parallel_block_profile(reset=.true.,after_trajectory=.true.)
  call b310_flush(rank)
  call flush312(rank,startup312)
  if(steps312>0) call checkpoint312(9901)
  ! Defer the map-override snapshot/checks until the existing checkpoint barriers.
  call snapshot_groups312("final")
  ! Existing checkpoint ends in all-rank barriers; no new output before it.
  call finish_runtime322(rank,steps312)
  close (12)
  if (rank == 0) write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
  call audit316_i("driver.steps",[steps312])
  call audit_transport_footer322()
  call audit316_finish()
  call finalize
end program climate
