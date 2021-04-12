program tke
  ! 1D test cases for TKE closure eddy viscosity model
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
  omega              = 0      * RAD/SECOND
  grav_accel         = 9.81 * METRE/SECOND**2     ! gravitational acceleration 
  p_top              = 0.0_8   * hPa                 ! pressure at free surface
  ref_density        = 1024    * KG/METRE**3         ! reference density at depth (seawater)

  ! Numerical method parameters
  default_thresholds = .true.                        ! use default threshold
  adapt_dt           = .true.                        ! adapt time step
  match_time         = .true.                        ! avoid very small time steps when saving (if false) 
  mode_split         = .true.                        ! split barotropic mode if true
  penalize           = .false.                        ! penalize land regions
  alpha              = 1d-6                         ! porosity used in penalization
  npts_penal         = 2.5                           ! number of points to smooth over in penalization
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
         
  ! Depth and layer parameters
  sigma_z            = .true.                        ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
  coords             = "croco"                       ! grid type for pure sigma grid ("croco" or "uniform")
  max_depth          = -50 * METRE                  ! total depth
  min_depth          = -50 * METRE                  ! minimum depth
  Tcline             = -50 * METRE                  ! position of thermocline

  ! Land mass parameters
  width              = 80 * KM                       ! width of zonal channel
  lat_width          = (width/radius)/DEG            ! width of zonal channel (in degrees)
  lat_c              = 45                            ! centre of zonal channel (in degrees)
  
  ! Bottom friction
  friction_upwelling = 0.0_8        

  ! Wind stress
  u_star             = 0.01 * METRE/SEC
  tau_0              = ref_density * u_star**2

  ! Other parameters
  alpha_0             = 2d-4
  N_0                 = 0.01 / SEC
  T_0                 = 16.0

  ! Vertical level to save
  save_zlev          = 8

  ! Characteristic scales
  wave_speed         = sqrt (grav_accel*abs(max_depth))  ! inertia-gravity wave speed 
  f0                 = 2*omega*sin(45*DEG)               ! representative Coriolis parameter
  beta               = 2*omega*cos(45*DEG) / radius      ! beta parameter at 30 degrees latitude
  Rd                 = wave_speed / f0                   ! barotropic Rossby radius of deformation                   

  ! Dimensional scaling
  drho               = 2.5 * KG/METRE**3                  ! magnitude of density perturbation
  Udim               = 0.35_8                              ! velocity scale
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
   
     call implicit_vertical_diffusion (bottom_friction, wind_tau, wind_d, source_b, source_t)

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        call remap_vertical_coordinates
        if (remap .and. modulo (istep, iremap) == 0) call remap_vertical_coordinates

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

  physics_scalar_flux = 0.0_8
  end if
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

  physics_velo_source = 0.0_8 
end function physics_velo_source
