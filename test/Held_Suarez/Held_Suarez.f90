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
  external       :: trend_cooling

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Initialize random number generator
  call initialize_seed

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6.371d6                     ! mean radius of the Earth in meters
  grav_accel     = 9.8_8                       ! gravitational acceleration in meters per second squared
  omega          = 7.292d-5                    ! Earth's angular velocity in radians per second
  p_0            = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
  c_p            = 1004.0_8                    ! specific heat at constant pressure in joules per kilogram Kelvin
  kappa          = 2.0_8/7.0_8                 ! kappa
  R_d            = kappa*c_p                   ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_v            = c_p - R_d                   ! specific heat at constant volume c_v = c_p - R_d
  gamma          = c_p/c_v                     ! heat capacity ratio
  gk             = grav_accel*kappa

  ! Local test case parameters
  T_0            = 300.0_8                     ! reference temperature
  T_mean         = 315.0_8                     ! mean temperature
  T_tropo        = 200.0_8                     ! tropopause temperature
  sigma_b        = 0.7_8                       ! normalized tropopause pressure height
  sigma_c        = 1.0_8-sigma_b                     
  k_a            = 1.0_8/(40*DAY)              ! cooling at free surface of atmosphere
  k_f            = 1.0_8/DAY                   ! Rayleigh friction
  k_s            = 1.0_8/(4*DAY)               ! cooling at surface
  delta_T        = 60.0_8                      ! meridional temperature gradient
  delta_theta    = 10.0_8                      ! vertical temperature gradient

  ! Local test case parameters (Jablonowski and Williamson 2006 zonally symmetric initial conditions)
  u_0            = 35.0_8                      ! maximum velocity of zonal wind
  gamma_T        = 5.0d-3                      ! temperature lapse rate
  delta_T2       = 4.8d5                       ! empirical temperature difference
  sigma_0        = 0.252_8                     ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2_8                       ! value of sigma at the tropopause

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                         ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 7.0d1                       ! temperature scale for tolerances
  Pdim           = p_0                         ! pressure scale
  dPdim          = 8.0d3                       ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim          ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim) ! acoustic wave speed

  Udim           = 30.0_8                      ! velocity scale
  Tdim           = DAY                         ! time scale
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
     if (time >= 100*DAY .and. modulo (istep, 20) == 0) call statistics
     call euler (trend_cooling, dt)
     call stop_timing

     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        if (remap .and. min_allowed_mass /= 1.0_8) call remap_vertical_coordinates

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
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  integer                                  :: id_i, v
  real(8)                                  :: dx, visc
  real(8), parameter                       :: C = 1d-2
  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: grad
  logical                                  :: local_type

  if (present(type)) then
     local_type = type
  else
     local_type = .false.
  end if

  id_i = id + 1  

  if (Laplace_order == 0) then
     physics_scalar_flux = 0.0_8
  else
     ! Calculate gradients
     if (Laplace_order == 1) then
        grad(S_MASS,:) = grad_physics (mass)
        grad(S_TEMP,:) = grad_physics (temp)
     elseif (Laplace_order == 2) then
        do v = S_MASS, S_TEMP
           grad(v,:) = grad_physics (Laplacian_scalar(v)%data(dom%id+1)%elts)
        end do
     end if

     ! Fluxes of physics
     if (.not.local_type) then ! Usual flux at edges E, NE, N
        do v = S_MASS, S_TEMP
           physics_scalar_flux(v,:) = visc_sclr(v) * grad(v,:) * dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
        end do
     else ! Flux at edges W, SW, S
        physics_scalar_flux(:,RT+1) = visc_sclr * grad(:,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(:,DG+1) = visc_sclr * grad(:,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(:,UP+1) = visc_sclr * grad(:,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)
     end if

     ! Find correct sign for diffusion on left hand side of the equation
     physics_scalar_flux = (-1)**Laplace_order * physics_scalar_flux
  end if
contains
  function grad_physics (scalar)
    implicit none
    real(8), dimension(1:EDGE) :: grad_physics
    real(8), dimension(:)      :: scalar

    if (.not.local_type) then ! Usual gradient at edges of hexagon E, NE, N
       grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1) 
       grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1) 
       grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1) 
    else ! Gradient for southwest edges of hexagon W, SW, S
       grad_physics(RT+1) = -(scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*idE+RT+1) 
       grad_physics(DG+1) = -(scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*idNE+DG+1)
       grad_physics(UP+1) = -(scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*idN+UP+1) 
    end if
  end function grad_physics
end function physics_scalar_flux

function physics_scalar_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the scalar trend
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

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

  integer                    :: id
  real(8), dimension(1:EDGE) :: diffusion

  id = idx (i, j, offs, dims)

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
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

subroutine trend_cooling (q, dq)
  ! Trend for Held-Suarez cooling
  use domain_mod
  use ops_mod
  use time_integr_mod
  implicit none
  type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, dq

  integer :: d, k, p

  call update_array_bdry (sol, NONE)

  ! Current surface pressure
  call cal_surf_press (sol)

  do k = 1, zlevels
     do d = 1, size(grid)
        mass  =>  q(S_MASS,k)%data(d)%elts
        temp  =>  q(S_TEMP,k)%data(d)%elts
        velo  =>  q(S_VELO,k)%data(d)%elts
        dmass => dq(S_MASS,k)%data(d)%elts
        dtemp => dq(S_TEMP,k)%data(d)%elts
        dvelo => dq(S_VELO,k)%data(d)%elts
        do p = 2, grid(d)%patch%length
           call apply_onescale_to_patch (cal_pressure,  grid(d), p-1, k, 0, 1)
           call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
           call apply_onescale_to_patch (trend_velo,    grid(d), p-1, k, 0, 0)
        end do
        nullify (mass, temp, velo, dmass, dtemp, dvelo)
     end do
  end do
  dq%bdry_uptodate = .false.
contains
  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Trend for cooling step (relaxation to equilibrium temperature)
    use test_case_mod
    use main_mod
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer     :: id_i
    real(8)     :: k_T, lat, lon, theta_equil

    id_i = idx (i, j, offs, dims) + 1
    call cart2sph (dom%node%elts(id_i), lon, lat)
    call cal_theta_eq (dom%press%elts(id_i), dom%surf_press%elts(id_i), lat, theta_equil, k_T)

    dmass(id_i) = 0.0_8
    dtemp(id_i) = - k_T * (temp(id_i) - theta_equil*mass(id_i))
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    ! Velocity trend for cooling step (Rayleigh friction)
    use test_case_mod
    use main_mod
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i
    real(8) :: k_v, sigma

    id = idx (i, j, offs, dims)
    id_i = id+1

    sigma = (dom%press%elts(id_i) - p_top) / (dom%surf_press%elts(id_i) - p_top)
    k_v = k_f * max (0.0_8, (sigma-sigma_b)/sigma_c)

    dvelo(EDGE*id+1:EDGE*id_i) = - k_v * velo(EDGE*id+1:EDGE*id_i)
  end subroutine trend_velo
end subroutine trend_cooling




