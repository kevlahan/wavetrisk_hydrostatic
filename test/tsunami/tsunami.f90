program Tsunami
  ! Simple tsunami test case
  !
  ! The permeability penalization parameter eta = dt_cfl and the porosity parameter alpha is set in the input file.
  use main_mod
  use test_case_mod
  use io_mod
  implicit none
  integer                                :: N
  integer, dimension(2)                  :: Nx, Ny
  real(8), dimension(:,:),   allocatable :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export
  real(8), dimension(2)                  :: lon_lat_range
  real(8), dimension(:,:,:), allocatable :: field2d_save
  logical                                :: aligned

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  radius         = 6371.229   * KM              ! mean radius of the Earth
  grav_accel     = 9.80616    * METRE/SECOND**2 ! gravitational acceleration 
  omega          = 7.29211d-5 * RAD/SECOND      ! Earth’s angular velocity 
  p_top          = 0.0_8      * hPa             ! pressure at free surface
  ref_density    = 1027       * KG/METRE**3     ! reference density (seawater)

  ! Local test case parameters
  min_depth      = -50  * METRE                 ! minimum allowed depth (must be negative)
  max_depth      = -4   * KM                    ! maximum allowed depth (must be negative)

  dH             =  7   * METRE                 ! initial perturbation to the free surface
  pert_radius    =  1e3 * KM                    ! radius of Gaussian free surface perturbation
  lon_c          = -50  * DEG                   ! longitude location of perturbation
  lat_c          =  25  * DEG                   ! latitude  location of perturbation

  ! Dimensional scaling
  wave_speed     = sqrt (grav_accel*abs(max_depth)) ! inertia-gravity wave speed based on maximum allowed depth
  Udim           = wave_speed                   ! velocity scale
  Ldim           = 2 * pert_radius              ! length scale (free surface perturbation width)
  Tdim           = Ldim/Udim                    ! time scale
  Hdim           = abs (max_depth)              ! vertical length scale

  ! Parameters for 2D projection
  N              = 1024                         ! size of lat-lon grid in 2D projection
  lon_lat_range  = (/2*MATH_PI, MATH_PI/)       ! region to save in 2D projection

  mode_split     = .false.                      ! use explicit time step for accuracy
  compressible   = .false.                      ! incompressible
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)

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
     call start_timing;  call time_step (dt_write, aligned); call stop_timing

     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint (run_id, rebalance)

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

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: d, id, id_i, visc_scale
  real(8), dimension(1:EDGE) :: diffusion

  d = dom%id + 1
  id = idx (i, j, offs, dims)
  id_i = id + 1

  visc_scale = 1.0_8!(max_level/dom%level%elts(id+1))**(2*Laplace_order_init-1)

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu()) * visc_scale
  end if

  ! Total physics for source term of velocity trend including volume penalization
  physics_velo_source = diffusion 
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





