program jet
  ! Baroclinic jet test case based on the beta plane configuration in Soufflet et al (Ocean Modelling 98, 36-50, 2016)
  use main_mod
  use test_case_mod
  use io_mod
  implicit none
  logical :: aligned
  
  ! Initialize mpi, shared variables and domainss
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  ! Parameters defining domain (based on beta-plane values)
  soufflet           = .false.                         ! set radius to exactly match Soufflet domain
  if (soufflet) then
     lat_c           = 30d0                            ! centre of zonal channel (in degrees)
     width           = 2000d0 * KM                     ! zonal channel width
     L_jet           = 0.8d0 * width                   ! width of jet transition region
     beta            = 1.6d-11 / (METRE * SECOND)      ! beta parameter
     f0              = 1d-4    / SECOND                ! Coriolis parameter
     omega           = f0 / (2d0*sin(lat_c*DEG))       ! planet rotation
     radius          = f0 / (beta * tan (lat_c*DEG))   ! planet radius to exactly match Soufflet beta plane
  else
     lat_c           = 30d0                            ! centre of zonal channel (in degrees)
     radius          = 1000d0 * KM                     ! meridional width of zonal channel
     width           = radius                          ! zonal channel width
     L_jet           = 0.4d0 * width                   ! width of jet transition region
     f0              = 1d-4  / SECOND                  ! Coriolis parameter
     omega           = f0 / (2d0*sin(lat_c*DEG))       ! planet rotation
     beta            = 2d0*omega*cos(lat_c*DEG)/radius ! beta parameter
  end if
  
  grav_accel         = 9.80616d0 * METRE/SECOND**2     ! gravitational acceleration 
  ref_density        = 1027.75d0 * KG/METRE**3         ! reference density at depth (maximum density)

  ! Numerical method parameters
  default_thresholds = .true.                          ! use default threshold
  
  match_time         = .true.                          ! avoid very small time steps when saving (if false) 
  penalize           = .true.                          ! penalize land regions
  alpha              = 1d-2                            ! porosity used in penalizatio
  npts_penal         = 4.5d0                           ! number of points to smooth over in penalization
  coarse_iter        = 20                              ! number of coarse scale iterations of elliptic solver
  fine_iter          = 20                              ! number of fine scale iterations of elliptic solver
  tol_elliptic       = 1d-8                            ! coarse scale tolerance of elliptic solver
  compressible       = .false.                         ! always run with incompressible equations
  remapscalar_type   = "PPR"                           ! optimal remapping scheme
  remapvelo_type     = "PPR"                           ! optimal remapping scheme

  ! Time stepping parameters
  timeint_type       = "RK3"                           ! time integration scheme
  adapt_dt           = .true.                          ! adapt time step
  mode_split         = .true.                          ! split barotropic mode if true
  theta1             = 0.8d0                           ! external pressure gradient (1 = fully implicit, 0.5 = Crank-Nicolson) stable if > 0.75
  theta2             = 0.8d0                           ! barotropic flow divergence (1 = fully implicit, 0.5 = Crank-Nicolson) stable if > 0.75
  
  ! Horizontal diffusion
  Laplace_order_init = 2                         
  Laplace_order      = Laplace_order_init

  ! Vertical diffusion
  vert_diffuse       = .true.                       
  tke_closure        = .true.
         
  ! Depth and layer parameters
  sigma_z            = .true.                          ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
  coords             = "croco"                         ! grid type for pure sigma grid ("croco" or "uniform")
  max_depth          = -4000d0 * METRE                 ! total depth
  min_depth          = max_depth                       ! minimum depth
  Tcline             =  -100d0 * METRE                 ! thermocline

  ! Land mass parameter
  lat_width          = (width/radius)/DEG              ! width of zonal channel (in degrees)
  
  ! Bottom friction
  bottom_friction_case = 5d-3 / SECOND

  ! Relaxation to initial zonal flow (nudging)
  tau_nudge          = 50d0 * DAY

  ! Wind stress
  tau_0              = 0d0

  ! Equation of state variables
  a_0                = 0.28d0 / CELSIUS
  b_0                = 0d0
  mu_1               = 0d0
  mu_2               = 0d0
  mu_1               = 0d0
  mu_2               = 0d0
  T_ref              = 14d0   * CELSIUS
  
  ! Vertical level to save
  save_zlev          = zlevels
!!$  save_zlev          = zlevels-5

  ! Characteristic scales
  L_pyc              = dz_b(2)                                      ! pycnocline
  drho               = drho_surf(2)                                 ! magnitude of density perturbation
  grav_reduced       = grav_accel * drho/ref_density                ! reduced gravity
  wave_speed         = sqrt (grav_accel*abs(max_depth))             ! inertia-gravity wave speed
  drho_dz            = drho / abs (L_pyc)                           ! approximate density gradient
  bv                 = sqrt (grav_accel * abs(drho_dz)/ref_density) ! Brunt-Vaisala frequency
  
!!$  c1                 = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification
!!$  Rb                 = bv * abs(max_depth) / (MATH_PI*f0) ! first baroclinic Rossby radius of deformation
  c1                 = sqrt (grav_reduced *  L_pyc)                 ! first baroclinic mode speed for linear stratification
  Rb                 = c1 / f0                                      ! alternate definition of first baroclinic Rossby radius of deformation
  
  Rd                 = wave_speed / f0                              ! barotropic Rossby radius of deformation                    
  
  ! Dimensional scaling
  if (soufflet) then ! velocity scale
     Udim = 0.2d0 * METRE/SECOND              
  else
     Udim = 1d0 * METRE/SECOND
  end if
  Ldim               = L_jet                             ! length scale 
  Tdim               = Ldim/Udim                         ! time scale
  Hdim               = abs (max_depth)                   ! vertical length scale

  iadapt             =  dt_init                           ! adapt horizontal grid every iadapt time steps
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize functions
  call assign_functions
  
  ! Initialize variables
  call initialize (run_id)

  ! Set interval for adapting grid based on the horizontal advective velocity scale (i.e. advect no more than one grid point before adapting)
  iadapt = CFL_adv * nint ((dx_min/Udim) / dt_init)

  ! Initialize 2D projections and zonal averages
  Nproj = sqrt (20d0 * 4**min_level) ! size of 2D projection: Nproj x Nproj/2
  call initialize_projection (Nproj)
  allocate (y2(Ny(1):Ny(2),1:zlevels,1:4),       y2_0(Ny(1):Ny(2),1:zlevels,1:4))
  allocate (zonal(Ny(1):Ny(2),1:zlevels,1:4), zonal_0(Ny(1):Ny(2),1:zlevels,1:4))
  call zonal_mean (zonal_0, y2_0) ; zonal = zonal_0 ; y2 = y2_0

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
  total_cpu_time = 0d0
  do while (time < time_end)
     call start_timing
     
     call time_step (dt_write, aligned)

     if (tau_nudge /= 0d0) then
        if (modulo (istep, 10) == 0) call zonal_mean (zonal, y2)
        call euler (sol, wav_coeff, trend_nudge, dt)
     end if

     call stop_timing

     call update_diagnostics

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (run_id, rebalance)
           call init_diagnostics !! resets diagnostics !!
        end if

        ! Save fields
        call vertical_velocity
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program jet

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

  if (implicit_diff_sclr .or. Laplace_order == 0 .or. maxval (visc_sclr) == 0d0) then
     physics_scalar_flux = 0.0_8
  else
     id_i = id + 1
     d = dom%id + 1
     
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

  physics_scalar_source = 0d0
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  ! wind stress and bottom friction are included as surface fluxes in the split eddy viscosity split step
  use test_case_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: id
  real(8), dimension(1:EDGE) :: diffusion

  id = idx (i, j, offs, dims)
  
  ! Horizontal diffusion
  if (implicit_diff_divu) then
     diffusion = - visc_rotu * curl_rotu()
  else
     diffusion = (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
  end if

  physics_velo_source = diffusion 
contains
  function grad_divu()
    implicit none
    real(8), dimension(3) :: grad_divu

    integer :: idE, idN, idNE

    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

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
