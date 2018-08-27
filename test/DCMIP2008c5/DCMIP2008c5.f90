program DCMIP2008c5
  ! DCMIP2008c5 test case 5: Mountain-induced Rossby wave
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none

  logical :: aligned

  ! Basic initialization of structures (grid, geometry etc)
  call init_main_mod 
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, qe, vort)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters ("test_case.in")
 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6.371229d6                   ! mean radius of the Earth in meters
  grav_accel     = 9.80616_8                    ! gravitational acceleration in meters per second squared
  omega          = 7.29211d-5                   ! Earth’s angular velocity in radians per second
  ref_press      = 100145.6_8                   ! reference pressure (mean surface pressure) in Pascals
  ref_surf_press = 930.0d2                      ! reference surface pressure
  R_d            = 287.04_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_p            = 1004.64_8                    ! specific heat at constant pressure in joules per kilogram Kelvin
  c_v            = 717.6_8                      ! specfic heat at constant volume c_v = R_d - c_p
  gamma          = c_p/c_v                      ! heat capacity ratio
  kappa          = 2.0_8/7.0_8                  ! kappa=R_d/c_p

  ! Local test case parameters
  d2             = 1.5d6**2                     ! square of half width of Gaussian mountain profile in meters
  h_0            = 2.0d3                        ! mountain height in meters
  lon_c          = MATH_PI/2                    ! longitude location of mountain
  lat_c          = MATH_PI/6                    ! latitude location of mountain
  T_0            = 288.0_8                      ! temperature in Kelvin
  u_0            = 20.0_8                       ! velocity in meters per second
  N_freq         = sqrt(grav_accel**2/(c_p*T_0)) ! Brunt-Vaisala buoyancy frequency

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                          ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 3.0d1                        ! temperature scale for tolerances
  Pdim           = ref_surf_press               ! pressure scale
  dPdim          = 2.5d3                        ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim           ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim)  ! acoustic wave speed

  Udim           = u_0                          ! velocity scale
  Ldim           = 2*sqrt(d2)                   ! length scale (mountain width)
  Tdim           = Ldim/Udim                    ! time scale (advection past mountain)
  Hdim           = wave_speed**2/grav_accel     ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  ! Initialize vertical grid
  call initialize_a_b_vert

  ! Determine vertical level to save
  call set_save_level

  ! Initialize thresholds to default values 
  call initialize_thresholds

  ! Initialize time step and viscosities
  call initialize_dt_viscosity
    
  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, test_case)
  call sum_total_mass (.true.)
  call barrier

  ! Save initial conditions
  call write_and_export (iwrite)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (rank == 0) write (6,'(A,/)') &
       '----------------------------------------------------- Start simulation run &
       ------------------------------------------------------'
  open (unit=12, file=trim (test_case)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  if (resume <= 0) iwrite = 0
  total_cpu_time = 0.0_8
  do while (time < time_end)
     call start_timing
     call time_step (dt_write, aligned, set_thresholds)
     call stop_timing
     
     call sum_total_mass (.false.)
     call print_log

     ! Save fields
     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates (set_thresholds)
        call write_and_export (iwrite)

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite,CP_EVERY) == 0) call write_checkpoint (dump, load, test_case, .false.)
     end if
  end do
  
  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program DCMIP2008c5

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

  if (max(viscosity_mass, viscosity_temp) == 0.0_8) then
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
  ! Additional physics for the source term of the scalar trend
  ! Newton cooling to equilibrium potential temperature theta_equil
  use domain_mod
  use test_case_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source(S_MASS) = 0.0_8
  physics_scalar_source(S_TEMP) = 0.0_8
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

  integer                      :: e
  real(8), dimension(1:EDGE) :: diffusion,  curl_rotu, grad_divu

  if (max(maxval(viscosity_divu), viscosity_rotu)==0.0_8) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     grad_divu = gradi_e (divu, dom, i, j, offs, dims)
     curl_rotu = curlv_e (vort, dom, i, j, offs, dims)
     do e = 1, EDGE 
        diffusion(e) = viscosity_divu(zlev) * grad_divu(e) - viscosity_rotu * curl_rotu(e)
     end do
  end if

  ! Find correct sign of diffusion on right hand side of equation
  diffusion = (-1)**(Laplace_order-1) * diffusion

  ! Total physics for source term of velocity trend
  do e = 1, EDGE
     physics_velo_source(e) =  diffusion(e)
  end do
end function physics_velo_source


