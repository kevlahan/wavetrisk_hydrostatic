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

  character(255) :: topo_operation

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Initialize random number generator
  call initialize_seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371      * KM                 ! mean radius of the Earth
  grav_accel     = 9.8       * METRE/SECOND**2    ! gravitational acceleration 
  omega          = 7.292d-5  * RAD/SECOND         ! Earth's angular velocity in radians per second
  p_0            = 1000      * hPa                ! reference pressure (mean surface pressure) in Pascals
  c_p            = 1004.0_8  * JOULE/(KG*KELVIN)  ! specific heat at constant pressure in joules per kilogram Kelvin
  kappa          = 2.0_8/7.0_8                    ! kappa
  R_d            = kappa*c_p * JOULE/(KG*KELVIN)  ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_v            = c_p - R_d * JOULE/(KG*KELVIN)  ! specific heat at constant volume c_v = c_p - R_d
  gamma          = c_p/c_v                        ! heat capacity ratio

  ! Local test case parameters
  T_0            = 300      * KELVIN              ! reference temperature
  T_mean         = 315      * KELVIN              ! mean temperature
  T_tropo        = 200      * KELVIN              ! tropopause temperature
  k_a            = 1.0_8/40 / DAY                 ! cooling at free surface of atmosphere
  k_f            = 1.0_8    / DAY                 ! Rayleigh friction
  k_s            = 1.0_8/4  / DAY                 ! cooling at surface
  delta_T        = 60       * KELVIN/METRE        ! meridional temperature gradient
  delta_theta    = 10       * KELVIN/METRE        ! vertical temperature gradient

  sigma_b        = 0.7_8                          ! normalized tropopause pressure height
  sigma_c        = 1.0_8-sigma_b
  
  ! Local test case parameters (Jablonowski and Williamson 2006 zonally symmetric initial conditions)
  u_0            = 35       * METRE/SECOND        ! maximum velocity of zonal wind
  gamma_T        = 5.0d-3   * KELVIN/METRE        ! temperature lapse rate
  delta_T2       = 4.8d5    * KELVIN              ! empirical temperature difference
  
  sigma_0        = 0.252_8                        ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2_8                          ! value of sigma at the tropopause

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                            ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 70       * KELVIN              ! temperature scale for tolerances
  Pdim           = p_0                            ! pressure scale
  dPdim          = 80       * hPa                 ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim          ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim) ! acoustic wave speed

  Udim           = 30 * METRE/SECOND           ! velocity scale
  Tdim           = 1  * DAY                    ! time scale
  Ldim           = Udim*Tdim                   ! length scale
  Hdim           = wave_speed**2/grav_accel    ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)
  call print_test_case_parameters

  ! Generate topography data
  topo_operation = "read"
  select case (topo_operation)
  case ("write") ! write out file descriptor for topography
     call write_grid_coords
  case ("read") ! read in geopotential
     call read_geopotential ("J06_gmted2010_modis_bedmachine_nc3000_NoAniso_Laplace0120_20231003.nc")
     call forward_topo_transform (topography, wav_topography)
     call inverse_topo_transform (wav_topography, topography)

     ! Check mass conservation
     fine_mass = integrate_hex (topo, z_null, level_end)
     do l = level_end-1, level_start, -1
        write (6,'(a,i2,a,es10.4)') &
             "Relative mass error at level ", l, " = ", abs (integrate_hex (topo, z_null, l) - fine_mass)/fine_mass
     end do
  end select

  ! Save initial conditions
  call write_and_export (iwrite)

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
     if (time >= 200*DAY .and. modulo (istep, 100) == 0) call statistics
     call euler (sol, wav_coeff, trend_cooling, dt)
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
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program Held_Suarez
