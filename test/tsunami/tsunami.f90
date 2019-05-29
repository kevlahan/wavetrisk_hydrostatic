program Tsunami
  ! Simple tsunami test case
  !
  ! The permeability penalization parameter eta = 4 dt_cfl and the porosity parameter alpha is set in the input file.
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
  radius         = 6371.229*KM                  ! mean radius of the Earth [m]
  grav_accel     = 9.80616_8                    ! gravitational acceleration [m/s^2]
  omega          = 7.29211d-5                   ! Earth’s angular velocity [rad/s]
  p_top          = 0.0_8                        ! pressure at free surface
  ref_density    = 1.0d3                        ! reference density (water) [kg/m^3]

  ! Local test case parameters
  mean_depth     = -3*KM                        ! mean depth coordinate [m] (must be negative)
  dH             = 1d-3 * abs (mean_depth)      ! perturbation to the free surface [m]
  pert_radius    = 1500*KM                      ! radius of Gaussian free surface perturbation [m]
  lon_c          = MATH_PI/2                    ! longitude location of perturbation
  lat_c          = MATH_PI/6                    ! latitude location of perturbation

  ! Dimensional scaling
  wave_speed     = sqrt (grav_accel*abs(mean_depth)) ! inertia-gravity wave speed
  Udim           = grav_accel*dH/wave_speed/2                   ! velocity scale
  Ldim           = 2 * pert_radius              ! length scale (free surface perturbation width)
  Tdim           = Ldim/Udim                    ! time scale (advection past mountain)
  Hdim           = abs (mean_depth)             ! vertical length scale
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Read test case parameters
  call read_test_case_parameters

  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, run_id)

  allocate (n_patch_old(size(grid)))

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
     n_patch_old = grid%patch%length
     call start_timing
     call time_step (dt_write, aligned, set_thresholds)
     call stop_timing

     call update_depth_penalization
     
     call sum_total_mass (.false.)
     call print_log

     if (aligned) then
        iwrite = iwrite + 1
        if (remap) call remap_vertical_coordinates

        ! Save checkpoint (and rebalance)
        if (modulo (iwrite, CP_EVERY) == 0) then
           call write_checkpoint (dump, load, run_id, rebalance)
           if (rebalance) call update_depth_penalization
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
end program Tsunami

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  integer                                  :: id_i, v
  real(8)                                  :: dx, visc
  real(8), parameter                       :: C = 1d-2
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
           physics_scalar_flux(v,:) = visc_sclr(v) * grad(v,:) * dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
        end do
     else ! Flux at edges W, SW, S
        physics_scalar_flux(:,RT+1) = visc_sclr * grad(:,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(:,DG+1) = visc_sclr * grad(:,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(:,UP+1) = visc_sclr * grad(:,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)
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
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source = 0.0_8
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  use domain_mod
  use time_integr_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: id, id_i
  real(8), dimension(1:EDGE) :: diffusion

  id   = idx (i, j, offs, dims)
  id_i = id + 1

  if (Laplace_order == 0) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
  end if

  ! Total physics for source term of velocity trend
  physics_velo_source =  diffusion - penal(zlev)%data(dom%id+1)%elts(id+1)/eta * velo(EDGE*id+1:EDGE*id_i)
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


