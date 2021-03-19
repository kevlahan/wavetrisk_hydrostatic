program upwelling
  ! Simple wind-driven upwelling/downwelling in a periodic zonal channel.
  ! An extension to the sphere of a test case used in ROMS and CROCO.
  use main_mod
  use test_case_mod
  use io_mod
  implicit none
  integer :: it
  logical :: aligned
  real(8) :: radius_earth

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard parameters
  radius             = 240     * KM
  omega              = 6d-5    * RAD/SECOND
  grav_accel         = 9.80616 * METRE/SECOND**2     ! gravitational acceleration 
  p_top              = 0.0_8   * hPa                 ! pressure at free surface
  ref_density        = 1027    * KG/METRE**3         ! reference density at depth (seawater)

  ! Numerical method parameters
  default_thresholds = .true.                        ! use default threshold
  adapt_dt           = .true.                        ! adapt time step
  match_time         = .true.                        ! avoid very small time steps when saving (if false) 
  mode_split         = .true.                        ! split barotropic mode if true
  penalize           = .true.                        ! penalize land regions
  alpha              = 1d-6                         ! porosity used in penalization
  npts_penal         = 4.5                           ! number of points to smooth over in penalization
  coarse_iter        = 20                            ! number of coarse scale iterations of elliptic solver
  fine_iter          = 20                            ! number of fine scale iterations of elliptic solver
  tol_elliptic       = 1d-8                          ! coarse scale tolerance of elliptic solver
  timeint_type       = "RK4"                         ! always use RK4
  compressible       = .false.                       ! always run with incompressible equations
  remapscalar_type   = "PPR"                         ! optimal remapping scheme
  remapvelo_type     = "PPR"                         ! optimal remapping scheme
  
  Laplace_order_init = 1                              
  Laplace_order = Laplace_order_init
  implicit_diff_sclr = .false.
  implicit_diff_divu = .false.
  vert_diffuse       = .true.                        ! include vertical diffusion
  
  coords             = "croco"                       ! chebyshev, croco, roms, or uniform
  
  ! Depth and layer parameters
  max_depth          = -150 * METRE                  ! total depth

  ! Land mass parameters
  width              = 80 * KM                       ! width of zonal channel
  lat_width          = (width/radius)/DEG            ! width of zonal channel (in degrees)
  lat_c              = 45                            ! centre of zonal channel (in degrees)
  
  ! Bottom friction
  friction_upwelling = 3d-4

  ! Vertical diffusion type
  rich_diff          = .false.                       ! richardson if true, roms if false               

  ! Wind stress
  tau_0              = 0.1_8

  ! Vertical level to save
  save_zlev          = 3

  ! Characteristic scales
  wave_speed         = sqrt (grav_accel*abs(max_depth))  ! inertia-gravity wave speed 
  f0                 = 2*omega*sin(45*DEG)               ! representative Coriolis parameter
  beta               = 2*omega*cos(45*DEG) / radius      ! beta parameter at 30 degrees latitude
  Rd                 = wave_speed / f0                   ! barotropic Rossby radius of deformation                   

  ! Dimensional scaling
  drho               = 2.5 * KG/METRE**3                  ! magnitude of density perturbation
  Udim               = 0.5_8                              ! velocity scale
  Ldim               = lat_width*DEG * radius             ! length scale 
  Tdim               = Ldim/Udim                          ! time scale
  Hdim               = abs (max_depth)                    ! vertical length scale
  
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
     call time_step (dt_write, aligned, &
          bottom_friction=friction_upwelling, wind_d=upwelling_drag, source_b=upwelling_bottom, source_t=upwelling_top)
     call stop_timing

     call update_diagnostics

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
        call vertical_velocity
        call write_and_export (iwrite)
     end if
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
  end if
  call finalize
end program upwelling

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

  if (implicit_diff_sclr .or. Laplace_order == 0 .or. maxval (visc_sclr) == 0.0_8) then
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

  physics_scalar_source = 0.0_8
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
