program DCMIP2012c4
  ! DCMIP (2012) test case 4: baroclinic instability
  ! Jablonowski and Williamson (2006) QJR Meteorol Soc (2006), 132, 2943–2975
  use main_mod
  use test_case_mod
  use io_mod
  use io_vtk_mod
  implicit none

  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371.229    * KM                ! mean radius of the Earth 
  grav_accel     = 9.80616     * METRE/SECOND**2   ! gravitational acceleration 
  omega          = 7.29212d-5  * RAD/SECOND        ! Earth’s angular velocity 
  p_0            = 1000        * hPa               ! reference pressure (mean surface pressure) 
  R_d            = 287         * JOULE/(KG*KELVIN) ! ideal gas constant for dry air
  c_p            = 1004.64     * JOULE/(KG*KELVIN) ! specific heat at constant pressure 
  c_v            = 717.6       * JOULE/(KG*KELVIN) ! specfic heat at constant volume c_v = R_d - c_p
  
  gamma          = c_p/c_v                         ! heat capacity ratio
  kappa          = 2.0_8/7.0_8                     ! kappa=R_d/c_p

  ! Local test case parameters
  u_0            = 35          * METRE/SECOND      ! maximum velocity of zonal wind
  u_p            = 1           * METRE/SECOND      ! maximum perturbation to zonal wind
  R_pert         = radius/10   * METRE             ! radius of perturbation to zonal wind
  T_0            = 288         * KELVIN            ! temperature in Kelvin
  gamma_T        = 5d-3        * KELVIN/METRE      ! temperature lapse rate
  delta_T        = 4.8d5       * KELVIN            ! empirical temperature difference
  lon_c          = MATH_PI/9   * RAD               ! longitude location of perturbation to zonal wind
  lat_c          = 2*MATH_PI/9 * RAD               ! latitude location of perturbation to zonal wind

  sigma_0        = 0.252_8                         ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2_8                           ! value of sigma at the tropopause
  
  ! Dimensions for scaling tendencies
  Tempdim        = T_0                             ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 80          * KELVIN            ! temperature scale for tolerances
  Pdim           = p_0                             ! pressure scale
  dPdim          = 80          * hPa               ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim              ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim)     ! acoustic wave speed

  Udim           = 3*u_0                           ! velocity scale (factor 3 from adaptive threshold runs)
  Tdim           = DAY                             ! time scale
  Ldim           = Udim*Tdim                       ! length scale
  Hdim           = wave_speed**2/grav_accel        ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize functions
  call assign_functions
 
  ! Initialize variables
  call initialize (run_id)

  ! Save initial conditions
  call print_test_case_parameters
  call write_and_export (iwrite)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim(run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     ! Time step
     call start_timing; call time_step (dt_write, aligned); call stop_timing
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        !if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)

        ! Save fields
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program DCMIP2012c4
