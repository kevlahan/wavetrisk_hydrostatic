!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: Simple_Physics.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Jan 6th
! Description: Program file to run the Simple physics with the test case used for the original driver by Thomas Dubos
!!!!!!!!!!!!!!!!!!!!!!!!
program Simple_Physics
   use main_mod
   use ops_mod
   use test_case_mod
   use io_mod
   implicit none

   logical        :: aligned
   character(256) :: input_file

   ! Initialize mpi, shared variables and domains
   call init_arch_mod
   call init_comm_mpi_mod

   ! Initialize random number generator
   call initialize_seed

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Standard (shared) parameter values for the simulation
   radius         = 6400      * KM                 ! mean radius of the Earth
   grav_accel     = 9.8       * METRE/SECOND**2    ! gravitational acceleration
   p_0            = 1000      * hPa                ! reference pressure (mean surface pressure) in Pascals
   p_top          = 0.01       * Pa                 ! pressure at the top in Pascals
   c_p            = 1004.0_8  * JOULE/(KG*KELVIN)  ! specific heat at constant pressure in joules per kilogram Kelvin
   !R_d            = 296.945007 * JOULE/(KG*KELVIN)  ! ideal gas constant for dry air in joules per kilogram Kelvin
   ! Set Rd to 287 whole number
   R_d            = 287 * JOULE/(KG*KELVIN)
   c_v            = c_p - R_d * JOULE/(KG*KELVIN)  ! specific heat at constant volume c_v = c_p - R_d
   ref_density    = 1.204     * KM                  ! Reference density (km/m^3)

   kappa          = R_d/c_p                    ! kappa
   gamma          = c_p/c_v                        ! heat capacity ratio

   ! Local test case parameters
   T_0            = 250      * KELVIN              ! reference temperature
   u_0            = 30       * METRE/SECOND        ! geostrophic wind speed
   e_thick        = 10       * KM                  ! Eckman Layer Thickness in meters
   soil_mod       = .TRUE.                         ! (T) soil module is on
   Nsoil          = 10                             ! Number of soil layers

   ! Dimensions for scaling tendencies
   Tempdim        = T_0                            ! temperature scale (both theta and T from DYNAMICO)
   dTempdim       = 70       * KELVIN              ! temperature scale for tolerances
   Pdim           = p_0                            ! pressure scale
   dPdim          = 80       * hPa                 ! scale of surface pressure variation determining mass tolerance scale

   ! Dimensional scaling
   specvoldim     = (R_d*Tempdim)/Pdim          ! specific volume scale
   wave_speed     = sqrt(gamma*Pdim*specvoldim) ! acoustic wave speed

   Udim           = 30 * METRE/SECOND           ! velocity scale
   Tdim           = 1  * DAY                    ! time scale - a solar day
   Ldim           = Udim*Tdim                   ! length scale
   Hdim           = wave_speed**2/grav_accel    ! vertical length scale
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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
   call write_and_export (iwrite)
   call mean_values(0)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (rank == 0) write (6,'(A,/)') &
      '----------------------------------------------------- Start simulation run &
      ------------------------------------------------------'
   open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
   total_cpu_time = 0.0_8
   do while (time < time_end)
      call start_timing
      call time_step (dt_write, aligned) ! dynamics step
      !When no dynamics called use : call timestep_placeholder(dt_write, aligned)
      if (time >= 200*DAY .and. modulo (istep, 100) == 0) call statistics
      call euler (sol(1:N_VARIABLE,1:ZLEVELS), wav_coeff(1:N_VARIABLE,1:ZLEVELS), trend_physics, dt) ! physics step
      call stop_timing

      call sum_total_mass (.false.)
      call print_log

      if (aligned) then
         iwrite = iwrite+1
         if (remap) call remap_vertical_coordinates

         if (modulo (iwrite, CP_EVERY) == 0) then
            call write_checkpoint (run_id, rebalance) ! save checkpoint (and rebalance)

            ! Save statistics
            call combine_stats
            if (rank == 0) call write_out_stats
         end if

         ! Save fields
         call write_and_export (iwrite)
         call mean_values(iwrite)
      end if
   end do

   call write_checkpoint (run_id, rebalance) ! checkpoint after final day
   call mean_values(INT(time_end)) ! save means of final day
   if (rank == 0) then
      close (12)
      write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
   end if
   call finalize
end program Simple_Physics
