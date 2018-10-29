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

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6.371d6                  ! mean radius of the Earth in meters
  grav_accel     = 9.8_8                       ! gravitational acceleration in meters per second squared
  omega          = 7.292d-5                    ! Earth’s angular velocity in radians per second
  ref_press      = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
  ref_surf_press = ref_press                   ! reference surface pressure
  R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_p            = 1004.0_8                    ! specific heat at constant pressure in joules per kilogram Kelvin
  c_v            = 717.6_8                     ! specific heat at constant volume c_v = R_d - c_p
  gamma          = c_p/c_v                     ! heat capacity ratio
  kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

  ! Local test case parameters
  T_0            = 300.0_8                     ! reference temperature
  T_mean         = 315.0_8                     ! mean temperature
  T_tropo        = 200.0_8                     ! tropopause temperature
  eta_b          = 0.7_8                       ! normalized tropopause pressure height
  k_a            = 1.0_8/(40*DAY)              ! cooling at free surface of atmosphere
  k_f            = 1.0_8/DAY                   ! Rayleigh friction
  k_s            = 1.0_8/(4*DAY)               ! cooling at surface
  delta_T        = 60.0_8                      ! meridional temperature gradient
  delta_theta    = 10.0_8                      ! vertical temperature gradient

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                         ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 7.0d1                       ! temperature scale for tolerances
  Pdim           = ref_surf_press              ! pressure scale
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
     call stop_timing

     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        if (remap) call remap_vertical_coordinates (set_thresholds)

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (dump, load, run_id)

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
  ! NOTE: call with arguments (dom, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  integer                                  :: v
  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: grad
  real(8), dimension(:), pointer           :: flx
  logical                                  :: local_type

  if (present(type)) then
     local_type = type
  else
     local_type = .false.
  end if

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
        physics_scalar_flux(S_MASS,RT+1) = viscosity_mass * grad(S_MASS,RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
        physics_scalar_flux(S_MASS,DG+1) = viscosity_mass * grad(S_MASS,DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
        physics_scalar_flux(S_MASS,UP+1) = viscosity_mass * grad(S_MASS,UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

        physics_scalar_flux(S_TEMP,RT+1) = viscosity_temp * grad(S_TEMP,RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
        physics_scalar_flux(S_TEMP,DG+1) = viscosity_temp * grad(S_TEMP,DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
        physics_scalar_flux(S_TEMP,UP+1) = viscosity_temp * grad(S_TEMP,UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
     else ! Flux at edges W, SW, S
        physics_scalar_flux(S_MASS,RT+1) = viscosity_mass * grad(S_MASS,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(S_MASS,DG+1) = viscosity_mass * grad(S_MASS,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(S_MASS,UP+1) = viscosity_mass * grad(S_MASS,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)

        physics_scalar_flux(S_TEMP,RT+1) = viscosity_temp * grad(S_TEMP,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(S_TEMP,DG+1) = viscosity_temp * grad(S_TEMP,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(S_TEMP,UP+1) = viscosity_temp * grad(S_TEMP,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)
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
  ! Newton cooling to equilibrium potential temperature theta_equil
  use domain_mod
  use test_case_mod
  use ops_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  integer                           :: id_i
  real(8)                           :: cooling, k_T, lat, lon, theta_equil

  id_i = idx (i, j, offs, dims) + 1

  ! Newton cooling
  call cart2sph (dom%node%elts(id_i), lon, lat)
  call cal_theta_eq (dom%press%elts(id_i)/ref_press, dom%press%elts(id_i)/dom%surf_press%elts(id_i), lat, theta_equil, k_T)
  cooling = -k_T * (temp(id_i) - theta_equil*mass(id_i))

  physics_scalar_source(S_MASS) = 0.0_8
  physics_scalar_source(S_TEMP) = cooling
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  !
  ! In this test case we add Rayleigh friction and Laplacian diffusion
  use domain_mod
  use ops_mod
  use test_case_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(Domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer :: e, id, id_i
  real(8) :: eta, k_v
  real(8), dimension(1:EDGE) :: diffusion, curl_rotu, grad_divu, Rayleigh

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     grad_divu = gradi_e (divu, dom, i, j, offs, dims)
     curl_rotu = curlv_e (vort, dom, i, j, offs, dims)
     diffusion =  (-1)**(Laplace_order-1) * (viscosity_divu(zlev) * grad_divu - viscosity_rotu * curl_rotu)
  end if
  
  ! Total physics for source term of velocity trend
  id = idx (i, j, offs, dims)
  id_i = id+1
  eta = dom%press%elts(id_i)/dom%surf_press%elts(id_i)
  k_v = k_f * max (0.0_8, (eta - eta_b)/(1.0_8-eta_b))
  do e = 1, EDGE
     Rayleigh(e) = -k_v * velo(EDGE*id+e)
  end do
  
  physics_velo_source = diffusion + Rayleigh
end function physics_velo_source

subroutine trend_cooling (q, dq)
  ! Trend for Held-Suarez cooling
  ! (not used currently since cooling and Rayleigh friction are handled in main trend)
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
           call apply_onescale_to_patch (cal_pressure, grid(d), p-1, k, 0, 1)
           call apply_onescale_to_patch (trend_mass,   grid(d), p-1, k, 0, 1)
           call apply_onescale_to_patch (trend_temp,   grid(d), p-1, k, 0, 1)
           call apply_onescale_to_patch (trend_velo,   grid(d), p-1, k, 0, 0)
        end do
        nullify (mass, temp, velo, dmass, dtemp, dvelo)
     end do
  end do
  dq%bdry_uptodate = .false.
contains
  subroutine trend_mass (dom, i, j, zlev, offs, dims)
    ! Mass trend for cooling step
    use test_case_mod
    use main_mod
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer     :: id_i

    id_i = idx(i, j, offs, dims) + 1
    dmass(id_i) = 0.0_8
  end subroutine trend_mass

  subroutine trend_temp (dom, i, j, zlev, offs, dims)
    ! Temperature trend for cooling step (relaxation to equilibrium temperature)
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

    id_i = idx (i, j, offs, dims)+1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       call cart2sph (dom%node%elts(id_i), lon, lat)
       call cal_theta_eq (dom%press%elts(id_i)/ref_press, dom%press%elts(id_i)/dom%surf_press%elts(id_i), lat, theta_equil, k_T)
       dtemp(id_i) = - k_T * (temp(id_i) - theta_equil*mass(id_i))
    else
       dtemp(id_i) = 0.0_8
    end if
  end subroutine trend_temp

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

    integer :: e, id, id_e, id_i
    real(8) :: k_v

    id = idx (i, j, offs, dims)
    id_i = id+1

    do e = 1, EDGE
       id_e = EDGE*id + e
       if (dom%mask_e%elts(EDGE*id+e) >= ADJZONE) then
          k_v = k_f * max (0.0_8, (dom%press%elts(id_i)/dom%surf_press%elts(id_i) - eta_b)/(1.0_8-eta_b))
          dvelo(id_e) = - k_v * velo(id_e)
       else
          dvelo(id_e) = 0.0_8
       end if
    end do
  end subroutine trend_velo
end subroutine trend_cooling




