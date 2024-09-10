!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: Simple_Physics.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Program file to run the Simple physics with the test case used for the original driver by Thomas Dubos
!!!!!!!!!!!!!!!!!!!!!!!!
program Simple_Physics
   use main_mod
   use ops_mod
   use test_case_mod
   use io_mod
   use init_physics_mod
   use physics_call_mod
   use phys_processing_mod

   implicit none
   logical        :: aligned
   character(256) :: input_file

   ! Initialize mpi, shared variables and domains
   call init_arch_mod
   call init_comm_mpi_mod

   ! Initialize random number generator
   call initialize_seed

   ! Numerical parameters
   adapt_dt       = .true.                           ! adapt time step
   compressible   = .true.                           ! compressible physics
   uniform        = .false.                          ! hybrid vertical pressure grid
   cfl_num        = 1d0                              ! cfl number
   Laplace_order_init = 2                            ! hyperdiffusion
   timeint_type   = "RK3"                            ! time integration scheme (use RK34, RK45 or RK4)
   iremap         = 10                               ! remap interval

   ! Standard (shared) parameter values for the simulation
   radius         = 6400d0    * KM                   ! mean radius of the Earth
   grav_accel     = 9.8d0     * METRE/SECOND**2      ! gravitational acceleration
   p_0            = 1000d0    * hPa                  ! reference pressure (mean surface pressure) in Pascals
   p_top          = 0.01d0    * Pa                   ! pressure at the top in Pascals
   c_p            = 1004d0    * JOULE/(KG*KELVIN)    ! specific heat at constant pressure in joules per kilogram Kelvin
   R_d            = 287d0     * JOULE/(KG*KELVIN)    ! set to a whole number
   c_v            = c_p - R_d * JOULE/(KG*KELVIN)    ! specific heat at constant volume c_v = c_p - R_d
   ref_density    = 1.204d0   * KG/METRE**3          ! reference density (kg/m^3)
   kappa          = R_d / c_p                        ! kappa
   gamma          = c_p / c_v                        ! heat capacity ratio

   ! Local initial conditions test case parameters
   T_0            = 250d0      * KELVIN              ! reference temperature
   u_0            = 30d0       * METRE/SECOND        ! geostrophic wind speed
   e_thick        = 10d0       * KM                  ! Eckman Layer Thickness in meters
  
   ! Dimensions for scaling tendencies
   Tempdim        = T_0                              ! temperature scale (both theta and T from DYNAMICO)
   dTempdim       = 70d0       * KELVIN              ! temperature scale for tolerances
   Pdim           = p_0                              ! pressure scale
   dPdim          = 80d0       * hPa                 ! scale of surface pressure variation determining mass tolerance scale

   ! Dimensional scaling
   specvoldim     = (R_d * Tempdim) / Pdim           ! specific volume scale
   wave_speed     = sqrt (gamma * Pdim * specvoldim) ! acoustic wave speed

   Udim           = u_0                              ! velocity scale
   Tdim           = 1d0  * DAY                       ! time scale - a solar day
   Ldim           = Udim * Tdim                      ! length scale
   Hdim           = wave_speed**2 / grav_accel       ! vertical length scale

   ! Read test case parameters
   call read_test_case_parameters
   
   ! Initialize functions
   call assign_functions

   ! Initialize physics grid parameters (it sets zmin if soil is being used)
   call init_soil_grid

   ! Initialize variables
   call initialize (run_id)

   ! Initialize the physics and physics function pointers
   call init_physics

   ! Physics call initializations if checkpointing
   if (cp_idx > 0) call physics_checkpoint_restart

   ! Save initial conditions
   call print_test_case_parameters
   
   !call write_and_export (iwrite)
   call mean_values (0) ! processing for the physics package mean values

   if (rank == 0) write (6,'(A,/)') &
      '----------------------------------------------------- Start simulation run &
      ------------------------------------------------------'
   open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
   
   total_cpu_time = 0d0
   do while (time < time_end)
      call start_timing
      call time_step (dt_write, aligned) ! dynamics step
      !When no dynamics called use : call timestep_placeholder(dt_write, aligned)
      call euler (sol(1:N_VARIABLE,1:ZLEVELS), wav_coeff(1:N_VARIABLE,1:ZLEVELS), trend_physics, dt) ! physics step
      call stop_timing
      call print_log

      if (aligned) then
         iwrite = iwrite+1
         if (remap) call remap_vertical_coordinates

         ! Save checkpoint (and rebalance)
         if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance) 

         ! Save fields
         call write_and_export (iwrite)
         call mean_values      (iwrite) 
      end if
   end do

   call write_checkpoint (run_id, rebalance) ! checkpoint after final day
   call mean_values(INT(time_end))           ! save means of final day
   
   if (rank == 0) then
      close (12)
      write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
   end if
   call finalize
end program Simple_Physics
