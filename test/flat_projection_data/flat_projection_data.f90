program flat_projection_data
  ! Post-processing of checkpoint data to calculate flat projection
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none

  integer        :: k
  real(8)        :: timing, total_cpu_time
  character(255) :: command
  logical        :: aligned, write_init

  ! Basic initialization of structures (grid, geometry etc)
  call init_main_mod 
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, qe, vort)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters ("flat_projection_data.in")

  ! DCMIP 2012c4
  if (trim(test_case) == 'DCMIP2012c4') then
     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
     omega          = 7.29212d-5                  ! Earth’s angular velocity in radians per second
     ref_press      = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = ref_press                   ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p
     u_0            = 35.0_8                      ! maximum velocity of zonal wind
     eta_0          = 0.252_8                     ! value of eta at reference level (level of the jet)
     compressible   = .true.                      ! Compressible equations
     uniform        = .false.                     ! Type of vertical grid
  else
     write(6,'(A)') "Test case not supported"
     stop
  end if
  resume = check_start
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  ! Initialize vertical grid
  call initialize_a_b_vert

  ! Initialize thresholds to default values 
  call initialize_thresholds

  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, test_case)
  call barrier

  call export_2d (cart2sph2, 3000000+100*cp_idx, (/-N/2, N/2/), (/-N/4, N/4/), (/2.0_8*MATH_PI, MATH_PI/), test_case)

  do cp_idx = resume+1, check_end
     resume = NONE
     call restart_full (set_thresholds, load, test_case)
     call export_2d (cart2sph2, 3000000+100*cp_idx, (/-N/2, N/2/), (/-N/4, N/4/), (/2.0_8*MATH_PI, MATH_PI/), test_case)
  end do
  call finalize
end program flat_projection_data

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type
  
  physics_scalar_flux(S_MASS,:) = 0.0_8
end function physics_scalar_flux

function grad_physics (scalar, dom, id, idE, idNE, idN, local_type)
  use domain_mod
  use test_case_mod
  implicit none

  real(8), dimension(1:EDGE)               :: grad_physics
  real(8), dimension(:), pointer           :: scalar
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical                                  :: local_type

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

function physics_scalar_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the scalar trend
  ! Newton cooling to equilibrium potential temperature theta_equil
  use domain_mod
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
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(Domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  physics_velo_source = 0.0_8
end function physics_velo_source


