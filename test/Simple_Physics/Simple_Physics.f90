!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: Simple_Physics.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Program file to run the Simple physics with the test case used for the original driver by Thomas Dubos
!!!!!!!!!!!!!!!!!!!!!!!!
program Simple_Physics
   use main_mod
   use test_case_mod
   use io_mod
   use init_physics_mod
   use physics_simple_mod
   use phys_processing_mod
   use callkeys, only : lverbose

   implicit none
   real(8)        :: area_sphere, dx_scaling, nu_CAM, nu, nu_scaling, res_scaling
   logical        :: aligned
   character(256) :: input_file

   ! Initialize mpi, shared variables and domains
   call init_arch_mod
   call init_comm_mpi_mod

   ! Initialize random number generator
   call initialize_seed

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Numerical parameters
   split_mean_perturbation = .false.                 ! split mean and perturbation (T)
   adapt_dt                = .false.                 ! adapt time step
   default_thresholds      = .false.                 ! thresholding type
   compressible            = .true.                  ! compressible physics
   uniform                 = .false.                 ! hybrid vertical pressure grid
   cfl_num                 = 1d0                     ! cfl number
   Laplace_order_init      = 2                       ! hyperdiffusion
   timeint_type            = "RK3"                   ! time integration scheme (use RK34, RK45 or RK4)
   iremap                  = 10                      ! remap interval
   log_min_mass            = .false.                 ! compute minimum mass at each dt (for checking stability issues)
   log_total_mass          = .false.                 ! check whether total mass is conserved (for debugging)

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

   ! Dimensions for scaling tendencies
   Tempdim        = T_0                              ! temperature scale (both theta and T from DYNAMICO)
   dTempdim       = 70d0       * KELVIN              ! temperature scale for tolerances
   Pdim           = p_0                              ! pressure scale
   dPdim          = 80d0       * hPa                 ! scale of surface pressure variation determining mass tolerance scale
   specvoldim     = (R_d * Tempdim) / Pdim           ! specific volume scale

   wave_speed     = sqrt (gamma * Pdim * specvoldim) ! acoustic wave speed

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
      cfl_min = cfl_max/2d0                                         ! minimum cfl number
      T_cfl   = DAY                                       ! time over which to increase cfl number from cfl_min to cfl_max
      cfl_num = cfl_min                                         ! initialize cfl number
      dt_init = cfl_min * 0.85d0 * dx_min / (wave_speed + u_0)  ! initial time step     (0.85 factor corrects for minimum dx)
      dt_max  = dt_max
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
   C_visc(S_MASS)     = nu_sclr * dt_init / (3d0 * Area_min**Laplace_order_init) 
   C_visc(S_TEMP)     = nu_sclr * dt_init / (3d0 * Area_min**Laplace_order_init) 
   C_visc(S_VELO)     = nu_rotu * dt_init / (3d0 * Area_min**Laplace_order_init) 
   C_div              = nu_divu * dt_init / (3d0 * Area_min**Laplace_order_init)

   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Physics Model Parameters !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   radiation_mod  = .true.                           ! (T) radiation module is on
   turbulence_mod = .true.                           ! (T) vertical diffusion module is on
   convecAdj_mod  = .true.                           ! (T) convective adjustment module is on
   soil_mod       = .false.                          ! (T) soil module is on (if F use surface flux only)

   ! Physics Package planet test case parameters
   gas_molarmass  = 28.9702532d0                     ! molar mass of main gain (used to set ideal gas const in pacakage)
   perihelion     = 150d0                            ! planet perihelion distance [MMkm]
   aphelion       = 150d0                            ! planet aphelion distance   [MMkm]
   perihelion_day = 0d0                              ! perihelion day
   obliquity      = 0d0                              ! planet axial tilt/obliquity
   sea_surf       = 0.01d0                           ! sea surface roughness length scale  [m]
   soil_surf      = 0.01d0                           ! soil surface roughness length scale [m]
   sea_inertia    = 3000d0                           ! sea thermal  inertia [J/m^3K]
   soil_inertia   = 3000d0                           ! soil thermal inertia [J/m^3K]
   sea_emissive   = 1d0                              ! sea emissivity
   soil_emmisive  = 1d0                              ! soil emissivity
   min_turbmix    = 100d0                            ! minimum turbulent mixing length [m]
   sw_atten       = 0.99d0                           ! attenuation of shortwave radiation coefficient
   lw_atten       = 0.08d0                           ! attenuation of longwave radiation coefficient
   seasons        = .false.                          ! seasons flag **** Does not do anything ****
   diurnal        = .true.                           ! diurnal cycle flag

   ! Single precision parameters
   sea_albedo     = 0.112e0                          ! sea albedo
   soil_albedo    = 0.112e0                          ! soil albedo
   emin_turb      = 1e-16                            ! minimum turbulent kinetic energy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
   ! Dimensional scaling
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
   lverbose = .false.

   ! Physics call initializations if checkpointing
   if (cp_idx > 0) call physics_checkpoint_restart

   ! Save initial conditions
   call print_test_case_parameters
   
   call write_and_export (iwrite)
   !call mean_values (0) ! processing for the physics package mean values

   if (rank == 0) write (6,'(a,/)') &
      '----------------------------------------------------- Start simulation run &
      ------------------------------------------------------'
   open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')

   total_cpu_time = 0d0; ; time_start = time
   do while (time < time_end)
      cfl_num = cfl (time) ! gradually increase cfl number
      
      call start_timing
      call time_step (dt_write, aligned)
      !call euler (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_physics_simple, dt) ! physics step
      !call WT_after_step (sol(:,zmin:0), wav_coeff(:,zmin:0), level_start-1)
      call stop_timing
      call print_log

      if (aligned) then
         iwrite = iwrite+1
         if (remap) call remap_vertical_coordinates

         ! Save checkpoint (and rebalance)
         if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance) 

         ! Save fields
         call write_and_export (iwrite) ; lverbose = .false.
         !call mean_values      (iwrite)
      end if
   end do

   call write_checkpoint (run_id, rebalance) ! checkpoint after final day
   !call mean_values (int (time_end))         ! save means of final day
   
   if (rank == 0) then
      close (12)
      write (6,'(a,es11.4)') 'Total cpu time = ', total_cpu_time
   end if
   call finalize
end program Simple_Physics
