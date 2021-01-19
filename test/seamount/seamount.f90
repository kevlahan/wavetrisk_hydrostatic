program Seamount
  ! Seamount case evaluates pressure gradient error (Beckmann and Haidvogel 1993). An isolated Gaussian bathymetry with flat isopycnals
  ! and hybrid sigma coordinates in z (vertical levels not aligned with isopycnals). The exact solution is a state of rest.
  ! The density difference is -3 kg/m^3 and decreases exponentially to zero with increasing depth. Parameters are set as in
  ! Debreu and Kevlahan (2020).
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  integer :: l
  real(8) :: total_eta
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius_earth   = 6371.229                * KM
  omega_earth    = 7.29211d-5              * RAD/SECOND
  grav_accel     = 9.80616                 * METRE/SECOND**2 ! gravitational acceleration 
  p_top          = 0.0_8                   * hPa             ! pressure at free surface
  ref_density    = 1000                    * KG/METRE**3     ! reference density at depth (seawater)

  scale          = 41.75d0                                     ! scaling factor for small planet
  radius         = radius_earth/scale                        ! mean radius of the small planet
  omega          = omega_earth                               ! angular velocity (scaled for small planet to keep beta constant)

  ! Depth and layer parameters
  min_depth   =   -50 * METRE                                ! minimum allowed depth (must be negative)
  max_depth   = -5000 * METRE                                ! total depth
  tau_0       =     0 * NEWTON/METRE**2                      ! maximum wind stress
  
  ! Seamount
  lat_c          = 43.29 * DEG                               ! latitude of seamount
  lon_c          =     0 * DEG                               ! longitude
  h0             =  4500 * METRE                             ! height of seamount
  width          =    40 * KM                                ! radius of seamount
  delta          =   500 * METRE                             ! vertical decay of density
  visc           =    50 * METRE**2/SECOND                   ! viscosity for rotu
  drag           = .false.                                   ! no bottom friction
  
  ! Numerical method parameters
  match_time         = .false.                               ! avoid very small time steps when saving 
  mode_split         = .true.                                ! split barotropic mode if true
  penalize           = .false.                               ! no penalization
  timeint_type       = "RK4"                                 ! always use RK4
  compressible       = .false.                               ! always run with incompressible equations
  mean_split         = .true.                                ! always split into mean and fluctuation (solve for fluctuation)
  default_thresholds = .true.                                ! always use default threshold
  adapt_dt           = .true.                                ! always adapt time steps
  remapscalar_type   = "PPR"                                 ! optimal remapping scheme
  remapvelo_type     = "PPR"                                 ! optimal remapping scheme
  Laplace_order_init = 1                              
  Laplace_order = Laplace_order_init

  ! Vertical level to save
  save_zlev = zlevels 

  ! Characteristic scales
  wave_speed     = sqrt (grav_accel*abs(max_depth))        ! inertia-gravity wave speed 
  f0             = 2*omega*sin(lat_c)                      ! representative Coriolis parameter
  beta           = 2*omega*cos(lat_c) / radius             ! beta parameter at 30 degrees latitude
  Rd             = wave_speed / f0                         ! barotropic Rossby radius of deformation                   
  drho_dz        = drho/max_depth                          ! approximate density gradient
  bv             = sqrt (grav_accel * drho_dz/ref_density) ! Brunt-Vaisala frequency
  Ro             = 0.0_8                                   ! Rossby number
  bu             = bv * abs(max_depth) / (f0 * width)      ! Burger number

  c1 = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification

  ! First baroclinic Rossby radius of deformation
  Rb = bv * abs(max_depth) / (MATH_PI*f0)

  ! Bottom friction
  bottom_friction = 0.0_8

  ! Internal wave friction 
  wave_friction = 0.0_8
  
  ! Relaxation of buoyancy to mean profile 
  k_T           = 0.0_8

  ! Dimensional scaling
  Udim           = 1.0_8                              ! velocity scale (arbitrary)
  Ldim           = width                              ! length scale 
  Tdim           = Ldim/Udim                          ! time scale
  Hdim           = abs (max_depth)                    ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize variables
  call initialize (run_id)

  ! Initialize diagnostic variables
  call init_diagnostics

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
     if (k_T /= 0.0_8) call euler (sol, wav_coeff, trend_relax, dt)
     call stop_timing

     call update_diagnostics

     tke = total_ke ("adaptive") / (4*MATH_PI*radius**2)
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (run_id, rebalance)
           call init_diagnostics !! resets diagnostics !!
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
end program Seamount

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

  if (Laplace_order == 0 .or. maxval (visc_sclr) == 0.0_8) then
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

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  use domain_mod
  use test_case_mod
  use ops_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: d, id, id_i
  real(8), dimension(1:EDGE) :: diffusion

  d = dom%id + 1
  id = idx (i, j, offs, dims)
  id_i = id + 1

  diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())

  physics_velo_source = diffusion
contains
  function grad_divu()
    implicit none
    real(8), dimension(3) :: grad_divu

    integer :: idE, idN, idNE

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    grad_divu(RT+1) = (divu(idE+1) - divu(id+1))   / dom%len%elts(EDGE*id+RT+1)
    grad_divu(DG+1) = (divu(id+1)  - divu(idNE+1)) / dom%len%elts(EDGE*id+DG+1)
    grad_divu(UP+1) = (divu(idN+1) - divu(id+1))   / dom%len%elts(EDGE*id+UP+1)
  end function grad_divu

  function curl_rotu()
    implicit none
    real(8), dimension(3) :: curl_rotu

    integer :: idS, idW

    idS = idx (i,   j-1, offs, dims)
    idW = idx (i-1, j,   offs, dims)

    curl_rotu(RT+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)
    curl_rotu(DG+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
    curl_rotu(UP+1) = (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
  end function curl_rotu
end function physics_velo_source
