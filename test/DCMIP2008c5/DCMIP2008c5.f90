program DCMIP2008c5
  ! DCMIP2008c5 test case 5: Mountain-induced Rossby wave
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none

  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371.229   * KM                 ! mean radius of the Earth
  grav_accel     = 9.80616    * METRE/SECOND**2    ! gravitational acceleration
  omega          = 7.29211d-5 * RAD/SECOND         ! Earth’s angular velocity
  p_0            = 1001.456   * hPa                ! reference pressure (mean surface pressure)
  ref_surf_press = 930        * hPa                ! reference surface pressure
  R_d            = 287.04     * JOULE/(KG*KELVIN)  ! ideal gas constant for dry air
  c_p            = 1004.64    * JOULE/(KG*KELVIN)  ! specific heat at constant pressure
  c_v            = 717.6_8    * JOULE/(KG*KELVIN)  ! specfic heat at constant volume c_v = R_d - c_p
  
  gamma          = c_p/c_v                         ! heat capacity ratio
  kappa          = 2.0_8/7.0_8                     ! kappa=R_d/c_p

  ! Local test case parameters
  d2             = (1500*KM)**2                    ! square of half width of Gaussian mountain profile
  h_0            = 2          * KM                 ! mountain height
  lon_c          = MATH_PI/2  * RAD                ! longitude location of mountain
  lat_c          = MATH_PI/6  * RAD                ! latitude location of mountain
  T_0            = 288        * KELVIN             ! temperature 
  u_0            = 20         * METRE/SECOND       ! velocity
  N_freq         = sqrt(grav_accel**2/(c_p*T_0))   ! Brunt-Vaisala buoyancy frequency

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                             ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 30         * KELVIN             ! temperature scale for tolerances
  Pdim           = ref_surf_press                  ! pressure scale
  dPdim          = 250        * hPa                ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim              ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim)     ! acoustic wave speed

  Udim           = u_0                             ! velocity scale
  Ldim           = 2*sqrt(d2)                      ! length scale (mountain width)
  Tdim           = Ldim/Udim                       ! time scale (advection past mountain)
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
  open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned)
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
  
  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program DCMIP2008c5
