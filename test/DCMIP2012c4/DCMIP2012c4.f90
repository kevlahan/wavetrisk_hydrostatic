program DCMIP2012c4
  ! DCMIP (2012) test case 4: baroclinic instability
  ! Jablonowski and Williamson (2006) QJR Meteorol Soc (2006), 132, 2943–2975
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
  radius         = 6.371229d6                  ! mean radius of the Earth in meters
  grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
  omega          = 7.29212d-5                  ! Earth’s angular velocity in radians per second
  p_0            = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
  R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_p            = 1004.64_8                   ! specific heat at constant pressure in joules per kilogram Kelvin
  c_v            = 717.6_8                     ! specfic heat at constant volume c_v = R_d - c_p
  gamma          = c_p/c_v                     ! heat capacity ratio
  kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

  ! Local test case parameters
  u_0            = 35.0_8                      ! maximum velocity of zonal wind
  u_p            = 1.0_8                       ! maximum perturbation to zonal wind
  R_pert         = radius/10.0_8               ! radius of perturbation to zonal wind
  T_0            = 288.0_8                     ! temperature in Kelvin
  gamma_T        = 5.0d-3                      ! temperature lapse rate
  delta_T        = 4.8d5                       ! empirical temperature difference
  sigma_0        = 0.252_8                     ! value of sigma at reference level (level of the jet)
  sigma_t        = 0.2_8                       ! value of sigma at the tropopause
  lon_c          = MATH_PI/9                   ! longitude location of perturbation to zonal wind
  lat_c          = 2*MATH_PI/9                 ! latitude location of perturbation to zonal wind

  ! Dimensions for scaling tendencies
  Tempdim        = T_0                         ! temperature scale (both theta and T from DYNAMICO)
  dTempdim       = 8.0d1                       ! temperature scale for tolerances
  Pdim           = p_0                         ! pressure scale
  dPdim          = 8.0d3                       ! scale of surface pressure variation determining mass tolerance scale

  ! Dimensional scaling
  specvoldim     = (R_d*Tempdim)/Pdim          ! specific volume scale
  wave_speed     = sqrt(gamma*Pdim*specvoldim) ! acoustic wave speed

  Udim           = 3*u_0                       ! velocity scale (factor 3 from adaptive threshold runs)
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
  open (unit=12, file=trim(run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
  total_cpu_time = 0.0_8
  do while (time < time_end)
     ! Time step
     call start_timing; call time_step (dt_write, aligned, set_thresholds); call stop_timing

     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite+1
        if (remap) call remap_vertical_coordinates

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
end program DCMIP2012c4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (dom, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  use ops_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  integer                                  :: id_i, v
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
           physics_scalar_flux(v,:) = viscosity_mass * grad(v,:) * dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
        end do
     else ! Flux at edges W, SW, S
        physics_scalar_flux(:,RT+1) = viscosity_mass * grad(:,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(:,DG+1) = viscosity_mass * grad(:,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(:,UP+1) = viscosity_mass * grad(:,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)
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
  use test_case_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source = 0.0_8
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

  real(8), dimension(1:EDGE) :: diffusion,  curl_rotu, grad_divu

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     grad_divu = gradi_e (divu, dom, i, j, offs, dims)
     curl_rotu = curlv_e (vort, dom, i, j, offs, dims)
     diffusion =  (-1)**(Laplace_order-1) * (viscosity_divu(zlev) * grad_divu - viscosity_rotu * curl_rotu)
  end if
  
  ! Total physics for source term of velocity trend
  physics_velo_source =  diffusion
end function physics_velo_source


