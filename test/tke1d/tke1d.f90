program tke1d
  ! 1D test cases for TKE closure eddy viscosity model
  use main_mod
  use test_case_mod
  use io_mod
  use vert_diffusion_mod
  implicit none
  integer              :: ialign, idt, it
  real(8)              :: align_time
  logical              :: aligned
  
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
  compressible       = .false.                       ! always run with incompressible equations
  adapt_dt           = .false.                        ! adapt time step
  match_time         = .true.                        ! avoid very small time steps when saving (if false) 
  penalize           = .false.                        ! penalize land regions
  remapscalar_type   = "PPR"                         ! optimal remapping scheme
  remapvelo_type     = "PPR"                         ! optimal remapping scheme
  vert_diffuse       = .true.                        ! include vertical diffusion
  tke_closure       = .true.                        ! use TKE closure to determine eddy viscosity and eddy diffusion
  
  ! Depth and layer parameters
  sigma_z            = .true.                        ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
  coords             = "croco"                       ! grid type for pure sigma grid ("croco" or "uniform")
  max_depth          = -50 * METRE                  ! total depth
  min_depth          = -50 * METRE                  ! minimum depth

  ! Bottom friction
  friction_tke       = 0.0_8        

  ! Wind stress  
  u_0                = 0.01 * METRE/SECOND
  tau_0              = ref_density * u_0**2

  ! Other parameters
  alpha_0             = 2d-4
  N_0                 = 0.01 / SECOND
  T_0                 = 16.0
  beta_eos            = 0.28_8

  ! Vertical level to save
  save_zlev          = 8

  ! Characteristic scales
  wave_speed         = sqrt (grav_accel*abs(max_depth))  ! inertia-gravity wave speed 

  ! Dimensional scaling
  Udim               = 0.35_8                              ! velocity scale
  Ldim               = radius             ! length scale 
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
     istep       = istep+1
     istep_cumul = istep_cumul+1
     dt = dt_new
     n_patch_old = grid%patch%length
     n_node_old  = grid%node%length
     idt    = nint (dt, 8)
     ialign = nint (align_time, 8)
     if (ialign > 0 .and. istep > 20) then
        aligned = (modulo (itime+idt,ialign) < modulo (itime,ialign))
     else
        aligned = .false.
     end if

     ! Modify time step
     if (aligned .and. match_time) then
        idt = ialign - modulo (itime,ialign)
        dt = idt
     end if

     call implicit_vertical_diffusion (r_d, tau, tke_drag, tke_bottom, tke_top)

     min_mass = cpt_min_mass ()
     call sum_total_mass (.false.)
     itime = itime + idt
     if (match_time) then
        time  = itime
     else
        time = time + dt
     end if
     dt_new = cpt_dt ()

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
end program tke1d

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
