module flat_projection_data_mod
  ! 2d projection from saved data
  use main_mod
  use remap_mod
  implicit none

  character(255)                     :: test_case 
  integer                            :: check_start, check_end, iwrite
  logical                            :: uniform
contains
   subroutine read_test_case_parameters (filename)
    implicit none
    character(*)   :: filename
    integer        :: fid = 500
    real(8)        :: pressures
    character(255) :: varname

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, check_start
    read (fid,*) varname, check_end
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, level_save
    read (fid,*) varname, pressures
    close(fid)

    if (rank==0) then
       write (6,'(A,A)')      "test_case          = ", test_case
       write (6,'(A,i12)')    "first file         = ", check_start
       write (6,'(A,i12)')    "first file         = ", check_end
       write (6,'(A,i3)')     "min_level          = ", min_level
       write (6,'(A,i3)')     "max_level          = ", max_level
       write (6,'(A,i3)')     "zlevels            = ", zlevels
       write (6,'(A,i3)')     "level_save         = ", level_save
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", pressures
       write (6,*) ' '
    end if
    pressure_save = (/pressures/) * 1.0d2 ! Convert to Pascals
  end subroutine read_test_case_parameters

  subroutine set_surf_geopot
    implicit none
    integer ::  d, p

    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (set_surfgeopot, grid(d), p-1, z_null, 0, 1)
       end do
    end do
  end subroutine set_surf_geopot
  
  subroutine set_surfgeopot (dom, i, j, zlev, offs, dims)
    ! Initialize surface geopotential after restart
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type (Coord) :: x_i
    integer      :: id

    id   = idx(i, j, offs, dims)
    x_i  = dom%node%elts(id+1)

    ! Surfaced geopotential
    dom%surf_geopot%elts(id+1) = surf_geopot_fun(x_i)
  end subroutine set_surfgeopot

   real(8) function surf_geopot_fun (x_i)
    ! Surface geopotential
    implicit none
    Type(Coord) :: x_i

    surf_geopot_fun = 0.0_8
  end function surf_geopot_fun
 
  subroutine initialize_a_b_vert
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    if (allocated(a_vert)) deallocate(a_vert)
    if (allocated(b_vert)) deallocate(b_vert)
    if (allocated(a_vert_mass)) deallocate(a_vert_mass)
    if (allocated(b_vert_mass)) deallocate(b_vert_mass)
    allocate (a_vert(1:zlevels), b_vert(1:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          press_infty = 0.0_8
          a_vert(k) = real(k-1)/real(zlevels) * press_infty/ref_press
          b_vert(k) = 1.0_8 - real(k-1)/real(zlevels)
       end do
    else
       if (zlevels==18) then
          a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1.0_8 /)
       elseif (zlevels==26) then
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
               0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
               0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1.0_8 /)
       elseif (zlevels==30) then
          a_vert = (/ 0.00225523952394724, 0.00503169186413288, 0.0101579474285245, 0.0185553170740604, 0.0306691229343414, &
               0.0458674766123295, 0.0633234828710556, 0.0807014182209969, 0.0949410423636436, 0.11169321089983, & 
               0.131401270627975, 0.154586806893349, 0.181863352656364, 0.17459799349308, 0.166050657629967, &
               0.155995160341263, 0.14416541159153, 0.130248308181763, 0.113875567913055, 0.0946138575673103, &
               0.0753444507718086, 0.0576589405536652, 0.0427346378564835, 0.0316426791250706, 0.0252212174236774, &
               0.0191967375576496, 0.0136180268600583, 0.00853108894079924, 0.00397881818935275, 0.0, 0.0 /)
          b_vert = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0393548272550106, &
               0.0856537595391273, 0.140122056007385, 0.204201176762581, 0.279586911201477, 0.368274360895157,  &
               0.47261056303978, 0.576988518238068, 0.672786951065063, 0.753628432750702, 0.813710987567902, &
               0.848494648933411, 0.881127893924713, 0.911346435546875, 0.938901245594025, 0.963559806346893, &
               0.985112190246582, 1.0 /)
       elseif (zlevels==49) then
          a_vert=(/0.002251865_8, 0.003983890_8, 0.006704364_8, 0.01073231_8, 0.01634233_8, 0.02367119_8, &
               0.03261456_8, 0.04274527_8, 0.05382610_8, 0.06512175_8, 0.07569850_8, 0.08454283_8, &
               0.08396310_8, 0.08334103_8, 0.08267352_8, 0.08195725_8, 0.08118866_8, 0.08036393_8, &
               0.07947895_8, 0.07852934_8, 0.07751036_8, 0.07641695_8, 0.07524368_8, 0.07398470_8, &
               0.07263375_8, 0.07118414_8, 0.06962863_8, 0.06795950_8, 0.06616846_8, 0.06424658_8, &
               0.06218433_8, 0.05997144_8, 0.05759690_8, 0.05504892_8, 0.05231483_8, 0.04938102_8, &
               0.04623292_8, 0.04285487_8, 0.03923006_8, 0.03534049_8, 0.03116681_8, 0.02668825_8, &
               0.02188257_8, 0.01676371_8, 0.01208171_8, 0.007959612_8, 0.004510297_8, 0.001831215_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
               0.006755112_8, 0.01400364_8, 0.02178164_8, 0.03012778_8, 0.03908356_8, 0.04869352_8, &
               0.05900542_8, 0.07007056_8, 0.08194394_8, 0.09468459_8, 0.1083559_8, 0.1230258_8, &
               0.1387673_8, 0.1556586_8, 0.1737837_8, 0.1932327_8, 0.2141024_8, 0.2364965_8, &
               0.2605264_8, 0.2863115_8, 0.3139801_8, 0.3436697_8, 0.3755280_8, 0.4097133_8, &
               0.4463958_8, 0.4857576_8, 0.5279946_8, 0.5733168_8, 0.6219495_8, 0.6741346_8, &
               0.7301315_8, 0.7897776_8, 0.8443334_8, 0.8923650_8, 0.9325572_8, 0.9637744_8, &
               0.9851122_8, 1.0_8/)
       else
          write(0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
          stop
       end if

       ! DCMIP order is opposite to ours
       if (.not. uniform) then
          a_vert = a_vert(zlevels+1:1:-1)
          b_vert = b_vert(zlevels+1:1:-1)
       end if

       ! Set pressure at infinity
       press_infty = a_vert(zlevels+1)*ref_press ! note that b_vert at top level is 0, a_vert is small but non-zero

       ! Set mass coefficients
       b_vert_mass = b_vert(1:zlevels)-b_vert(2:zlevels+1)
       a_vert_mass = ((a_vert(1:zlevels)-a_vert(2:zlevels+1))*ref_press + b_vert_mass*press_infty)/grav_accel
    end if
  end subroutine initialize_a_b_vert
 
  subroutine dump (fid)
    implicit none
    integer :: fid

    write(fid) itime
    write(fid) iwrite
    write(fid) tol_mass, tol_temp, tol_velo
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid

    read(fid) itime
    read(fid) iwrite
    read(fid) tol_mass, tol_temp, tol_velo
  end subroutine load

  subroutine set_thresholds
    implicit none
  end subroutine set_thresholds
    
  subroutine apply_initial_conditions
    implicit none
  end subroutine apply_initial_conditions
end module flat_projection_data_mod

program flat_projection_data
  use main_mod
  use flat_projection_data_mod
  implicit none

  character(255) :: command

  ! Initialize grid etc
  call init_main_mod 

  ! Nullify all pointers initially
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, qe, vort)

  save_levels = 1; allocate(pressure_save(1:save_levels))  ! number of vertical levels to save

  ! Read test case parameters
  call read_test_case_parameters ("flat_projection_data.in")

  allocate (tol_mass(1:zlevels), tol_temp(1:zlevels), tol_velo(1:zlevels))

  ! Parameters for the simulation
  radius         = 6.371229d6   ! mean radius of the Earth in meters
  grav_accel     = 9.80616_8    ! gravitational acceleration in meters per second squared
  ref_press      = 1.0d5        ! reference pressure (mean surface pressure) in Pascals
  ref_surf_press = ref_press    ! reference surface pressure
  R_d            = 287.0_8      ! ideal gas constant for dry air in joules per kilogram Kelvin
  kappa          = 2.0_8/7.0_8  ! kappa=R_d/c_p
  compressible   = .true.       ! Compressible equations
  uniform        = .false.      ! Type of vertical grid

   ! Initialize vertical grid
  call initialize_a_b_vert

  ! Initialize variables
  resume = check_start
  call initialize (apply_initial_conditions, set_thresholds, dump, load, test_case)

  call set_surf_geopot 

  call barrier

  call export_2d (cart2sph2, 3000000+100*cp_idx, (/-256, 256/), (/-128, 128/), (/2.0_8*MATH_PI, MATH_PI/), test_case)

  do cp_idx = resume+1, check_end
     resume = NONE
     call restart_full (set_thresholds, load, test_case)
     call export_2d (cart2sph2, 3000000+100*cp_idx, (/-256, 256/), (/-128, 128/), (/2.0_8*MATH_PI, MATH_PI/), test_case)
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
  use flat_projection_data_mod
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
  use flat_projection_data_mod
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


