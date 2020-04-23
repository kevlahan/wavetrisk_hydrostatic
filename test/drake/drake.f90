program Drake
  ! Simplified Drake passage test case on small planet
  ! (inspired by Ferreira, Marshall and Rose 2011, J Climate 24, 992-1012)
  !
  ! The permeability penalization parameter eta = dt_cfl and the porosity parameter alpha is set in the input file.
  ! The initial condition has a 500 m thick less dense layer over a 1500 m thick reference density layer (only the upper
  ! layer is included in the single vertical layer shallow water case).
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  logical :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  scale          = 6                                  ! scale factor for small planet (1/6 Earth radius)
  radius         = 6371.229/scale   * KM              ! mean radius of the small planet
  grav_accel     = 9.80616          * METRE/SECOND**2 ! gravitational acceleration 
  omega          = 7.29211d-5/scale * RAD/SECOND      ! angular velocity (scaled for small planet to keep beta constant)
  p_top          = 0.0_8            * hPa             ! pressure at free surface
  ref_density    = 1028             * KG/METRE**3     ! reference density at depth (seawater)
  drho           = -1               * KG/METRE**3     ! density difference in upper layer

  ! Numerical method parameters
  timeint_type       = "RK4"                          ! always use RK4
  compressible       = .false.                        ! always run with incompressible equations
  mean_split         = .true.                         ! always split into mean and fluctuation (solve for fluctuation)
  remapscalar_type   = "2PPM"                         ! optimal remapping scheme
  remapvelo_type     = "2PPM"                         ! optimal remapping scheme
  Laplace_order_init = 1                              ! always run with Laplacian diffusion
  Laplace_order = Laplace_order_init
  
  ! Local test case parameters
  min_depth      = -50 * METRE                        ! minimum allowed depth (must be negative)
  if (zlevels == 1) then                              ! maximum allowed depth (must be negative)
     max_depth   =  -500 * METRE                     
  else
     max_depth   = -2000 * METRE
!!$     max_depth   = -4000 * METRE
  end if
  
  top_layer      = -500 * METRE                       ! location of top (less dense) layer in two layer case
  
  ! Characteristic scales
  wave_speed     = sqrt (grav_accel*abs(max_depth))                            ! inertia-gravity wave speed 
  f0             = 2*omega*sin(30*DEG)          * RAD/SECOND                   ! representative Coriolis parameter
  beta           = 2*omega*cos(30*DEG) / radius * RAD/(SECOND*METRE)           ! beta parameter at 30 degrees latitude
  L_R            = wave_speed / f0              * METRE                        ! Rossby radius
  bv             = sqrt (grav_accel/ref_density*abs(drho)/abs(top_layer))      ! Brunt-Vaisala frequency         
  c1             = sqrt (bv**2*abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed 
!!$  u_wbc          = 4                            * METRE/SECOND                 ! western boundary current (H = 1 km, top = 500 m)
  u_wbc          = 2                            * METRE/SECOND                 ! western boundary current (H = 2 km, top = 500 m)
!!$  u_wbc          = 1                            * METRE/SECOND                 ! western boundary current (H = 4 km, top = 1000 m)
  delta_I        = sqrt (u_wbc/beta)            * METRE                        ! inertial layer
  delta_sm       = u_wbc/f0                     * METRE                        ! barotropic submesoscale  

  ! Dimensional scaling
  Udim           = u_wbc                              ! velocity scale
  Ldim           = delta_I                            ! length scale 
  Tdim           = Ldim/Udim                          ! time scale
  Hdim           = abs (max_depth)                    ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, run_id)

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
     call start_timing;  call time_step (dt_write, aligned, set_thresholds); call stop_timing

     call update_diagnostics
     
     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call deallocate_diagnostics
           call write_checkpoint (dump, load, run_id, rebalance)
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
end program

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
  use test_case_mod
  use ops_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                         :: d, id, id_i, idE, idN, idNE
  real(8)                         :: r
  real(8), dimension(0:NORTHEAST) :: full_mass
  real(8), dimension(1:EDGE)      :: diffusion, drag_force, mass_e, permeability, stress, tau_wind, u_mag, wind_force

  d = dom%id + 1
  id = idx (i, j, offs, dims)
  id_i = id + 1

  idE  = idx (i+1, j,   offs, dims) + 1
  idNE = idx (i+1, j+1, offs, dims) + 1
  idN  = idx (i,   j+1, offs, dims) + 1

  ! Interpolate mass to edges
  full_mass(0:NORTHEAST) = mass((/id,idN,idE,id,id,idNE/)+1) + mean_m((/id,idN,idE,id,id,idNE/)+1)
  mass_e(RT+1) = interp (full_mass(0), full_mass(EAST))
  mass_e(DG+1) = interp (full_mass(0), full_mass(NORTHEAST))
  mass_e(UP+1) = interp (full_mass(0), full_mass(NORTH))

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
  end if

  ! Wind stress per unit length in top layer only
  if (zlev == zlevels) then
     idE  = idx (i+1, j,   offs, dims) + 1
     idNE = idx (i+1, j+1, offs, dims) + 1
     idN  = idx (i,   j+1, offs, dims) + 1

     tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id_i), dom%node%elts(idE))
     tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE), dom%node%elts(id_i))
     tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id_i), dom%node%elts(idN))
  else
     tau_wind = 0.0_8
  end if
  wind_force = tau_wind / mass_e
  
  ! Bottom stress applied in lowest layer only (as in NEMO)
  if (zlev == 1) then
     ! Quadratic 
!!$     u_mag(RT+1) = sqrt (2 * interp (ke(id_i), ke(idE)))
!!$     u_mag(DG+1) = sqrt (2 * interp (ke(id_i), ke(idNE)))
!!$     u_mag(UP+1) = sqrt (2 * interp (ke(id_i), ke(idN)))
!!$     drag_force = - friction_coeff * u_mag * velo(EDGE*id+RT+1:EDGE*id+UP+1)
     ! Linear
     drag_force = - friction_coeff * velo(EDGE*id+RT+1:EDGE*id+UP+1)
  else
     drag_force = 0.0_8
  end if
  drag_force = drag_force / (mass_e/ref_density)

  ! Internal wave drag to reduce oscillation amplitude (energy neutral)
  if (zlevels == 2) then
     r = friction_coeff / abs(max_depth)
     if (zlev == 1) then
        drag_force = drag_force - r * (velo(EDGE*id+RT+1:EDGE*id+UP+1) - sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
     elseif (zlev == 2) then
        drag_force = drag_force - r * (velo(EDGE*id+RT+1:EDGE*id+UP+1) - sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
     end if
  else
     drag_force = 0.0_8
  end if
  
  ! Permeability
  permeability = - penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)/eta * velo(EDGE*id+RT+1:EDGE*id+UP+1)
  
  ! Complete source term for velocity trend
  ! (do not include drag and wind stress in solid regions)
  if (penal_node(zlevels)%data(d)%elts(id_i) == 0.0_8) then ! sea
     physics_velo_source = diffusion + permeability + drag_force + wind_force
  else ! land
     physics_velo_source = diffusion + permeability 
  end if
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





