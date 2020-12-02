program Held_Suarez
  ! Held & Suarez (1994) test case
  ! Bulletin of the American Meteorological Society 75 (10), 1825-1830
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
  radius         = 6371      * KM                 ! mean radius of the Earth
  grav_accel     = 9.8       * METRE/SECOND**2    ! gravitational acceleration 
  omega          = 7.292d-5  * RAD/SECOND         ! Earth's angular velocity in radians per second
  p_0            = 1000      * hPa                ! reference pressure (mean surface pressure) in Pascals
  c_p            = 1004.0_8  * JOULE/(KG*KELVIN)  ! specific heat at constant pressure in joules per kilogram Kelvin
  R_d            = kappa*c_p * JOULE/(KG*KELVIN)  ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_v            = c_p - R_d * JOULE/(KG*KELVIN)  ! specific heat at constant volume c_v = c_p - R_d
  
  kappa          = 2.0_8/7.0_8                    ! kappa
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

  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, run_id)

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
     call time_step (dt_write, aligned, set_thresholds)
     if (time >= 200*DAY .and. modulo (istep, 100) == 0) call statistics
     call euler (sol, wav_coeff, trend_cooling, dt)
     call stop_timing

     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        if (remap) call remap_vertical_coordinates

        if (modulo (iwrite, CP_EVERY) == 0) then
           call write_checkpoint (dump, load, run_id, rebalance) ! save checkpoint (and rebalance)

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

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (q, dom, id, idE, idNE, idN, v, zlev, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)                           :: physics_scalar_flux
  type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
  type(domain)                                         :: dom
  integer                                              :: d, id, idE, idNE, idN, v, zlev
  logical, optional                                    :: type

  integer                    :: id_i
  real(8), dimension(1:EDGE) :: d_e, grad, l_e
  logical                    :: local_type

  if (present(type)) then
     local_type = type
  else
     local_type = .false.
  end if

  id_i = id + 1
  d = dom%id + 1

  if (Laplace_order == 0) then
     physics_scalar_flux = 0.0_8
  else
     if (.not.local_type) then ! usual flux at edges E, NE, N
        l_e =  dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
        d_e =  dom%len%elts(EDGE*id+1:EDGE*id_i)
     else ! flux at SW corner
        l_e(RT+1) = dom%pedlen%elts(EDGE*idE+RT+1)
        l_e(DG+1) = dom%pedlen%elts(EDGE*idNE+DG+1)
        l_e(UP+1) = dom%pedlen%elts(EDGE*idN+UP+1)
        d_e(RT+1) = -dom%len%elts(EDGE*idE+RT+1)
        d_e(DG+1) = -dom%len%elts(EDGE*idNE+DG+1)
        d_e(UP+1) = -dom%len%elts(EDGE*idN+UP+1)
     end if
     
     ! Calculate gradients
     if (Laplace_order == 1) then
        grad = grad_physics (q(v,zlev)%data(d)%elts)
     elseif (Laplace_order == 2) then
        grad = grad_physics (Laplacian_scalar(v)%data(d)%elts)
     end if

     ! Complete scalar diffusion
     physics_scalar_flux = (-1)**Laplace_order * visc_sclr(v) * grad * l_e
  end if
contains
  function grad_physics (scalar)
    implicit none
    real(8), dimension(1:EDGE) :: grad_physics
    real(8), dimension(:)      :: scalar

    grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))   / d_e(RT+1)
    grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1)) / d_e(DG+1)
    grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))   / d_e(UP+1)
  end function grad_physics
end function physics_scalar_flux

function physics_scalar_source (q, id, zlev)
  ! Additional physics for the source term of the scalar trend
  use domain_mod
  implicit none
  real(8), dimension(scalars(1):scalars(2))            :: physics_scalar_source
  integer                                              :: id, zlev
  type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q

  physics_scalar_source = 0.0_8
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: id, visc_scale
  real(8), dimension(1:EDGE) :: diffusion

  id = idx (i, j, offs, dims)

  visc_scale = 1.0_8!(max_level/dom%level%elts(id+1))**(2*Laplace_order_init-1)

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu()) * visc_scale
  end if

  ! Total physics for source term of velocity trend
  physics_velo_source =  diffusion
contains
  function grad_divu()
    implicit none
    real(8), dimension(3) :: grad_divu

    integer :: idE, idN, idNE

    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    grad_divu(RT+1) = (divu(idE+1) - divu(id+1))  /dom%len%elts(EDGE*id+RT+1)
    grad_divu(DG+1) = (divu(id+1)  - divu(idNE+1))/dom%len%elts(EDGE*id+DG+1)
    grad_divu(UP+1) = (divu(idN+1) - divu(id+1))  /dom%len%elts(EDGE*id+UP+1)
  end function grad_divu

  function curl_rotu()
    implicit none
    real(8), dimension(3) :: curl_rotu

    integer :: idS, idW

    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    
    curl_rotu(RT+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1))/dom%pedlen%elts(EDGE*id+RT+1)
    curl_rotu(DG+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+DG+1)
    curl_rotu(UP+1) = (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1))/dom%pedlen%elts(EDGE*id+UP+1)
  end function curl_rotu
end function physics_velo_source




