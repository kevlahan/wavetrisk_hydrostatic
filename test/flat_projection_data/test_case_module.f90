module test_case_mod
  ! Module file for flat_projection_data
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use equation_of_state_mod
  use ops_mod
  use std_atm_profile_mod
  use io_mod
  implicit none
  integer :: mean_beg, mean_end, cp_2d, N
  real(8) :: npts_penal, ref_surf_press, scale
  real(8) :: dPdim, R_ddim, specvoldim, dTempdim
  logical :: mean_split, zonal

  ! DCMIP2012c4
  real(8) :: eta_0, u_0 
  ! DCMIP2008c5
  real(8) :: d2, h_0, lat_c, lon_c
  ! Drake
  real(8)               :: drho, halocline, mixed_layer
  real(8), dimension(2) :: density_drake, height
  ! Seamount
  real(8) :: delta, h0, lat_val, lon_val, width
  character(255) :: coords, stratification
  ! Upwelling
  real(8)            :: lat_width, Tcline, T_0
  real(8), parameter :: slope = 1.75d-4 ! slope parameter (larger value -> steeper slope)
  real(8), parameter :: shift = 8
  ! Jet
  real(8) :: beta, f0, L_jet
  logical :: soufflet
  ! Simple Physics
  logical :: climatology
  type(Float_Field), dimension(:), allocatable, target :: simple_phys_temp, simple_phys_zonal, simple_phys_merid, simple_phys_vels
  ! Held Suarez
  real(8) :: delta_T, delta_theta, sigma_b, sigma_c, k_a, k_f, k_s, T_mean, T_tropo
  real(8) :: delta_T2, sigma_t, sigma_v, sigma_0, gamma_T
  real(8) :: cfl_max, cfl_min, T_cfl
contains
  subroutine assign_functions
    ! Assigns generic pointer functions to functions defined in test cases
    implicit none

    ! Standard functions
    apply_initial_conditions => apply_initial_conditions_case
    dump                     => dump_case
    load                     => load_case
    initialize_a_b_vert      => initialize_a_b_vert_case
    initialize_dt_viscosity  => initialize_dt_viscosity_case
    initialize_thresholds    => initialize_thresholds_case
    physics_scalar_flux      => physics_scalar_flux_case
    physics_velo_source      => physics_velo_source_case
    set_save_level           => set_save_level_case
    set_thresholds           => set_thresholds_case
    surf_geopot              => surf_geopot_case
    update                   => update_case
    z_coords                 => z_coords_case
  end subroutine assign_functions

  function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type)
    use domain_mod
    implicit none
    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
    type(domain)                                         :: dom
    integer                                              :: d, id, idE, idNE, idN, v, zlev
    logical, optional                                    :: type

    physics_scalar_flux_case = 0.0_8
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none
    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    physics_velo_source_case = 0.0_8
  end function physics_velo_source_case

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential
    implicit none
    integer :: d, id
    
    type(Coord) :: x_i
    real(8)     :: amp, b_max, c1, cs2, sn2, lon, lat, rgrc, y

    x_i = grid(d)%node%elts(id)

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2

    if (trim (test_case) == "DCMIP2012c4") then
       c1 = u_0*cos((1.0_8-eta_0)*MATH_PI/2)**1.5

       surf_geopot_case = c1 * (c1 * (-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8)  + &
            radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
    elseif (trim (test_case) == "DCMIP2008c5") then
       rgrc = radius*acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

       surf_geopot_case = grav_accel*h_0*exp__flush (-rgrc**2/d2)
    elseif (trim (test_case) == "Held_Suarez") then
       if (NCAR_topo) then ! add non-zero surface geopotential
          surf_geopot_case = grav_accel * topography%data(d)%elts(id)
       else ! surface geopotential from Jablonowski and Williamson (2006)
          call cart2sph (grid(d)%node%elts(id), lon, lat)

          c1 = u_0 * cos((1d0 - sigma_0) * MATH_PI/2d0)**1.5
          cs2 = cos (lat)**2; sn2 = sin (lat)**2

          surf_geopot_case = c1 * (c1 * (-2d0 * sn2**3 * (cs2 + 1d0/3d0) + 10d0/63d0) &
               + radius * omega * (8d0/5d0 * cs2**1.5 * (sn2 + 2d0/3d0) - MATH_PI/4d0)) &
               + grav_accel * topography%data(d)%elts(id)
       end if
    elseif (trim (test_case) == "seamount") then
       rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))
       surf_geopot_case = grav_accel*h0 * exp__flush (-(rgrc/width)**2)
    elseif (trim (test_case) == "upwelling") then
       b_max = abs (max_depth - min_depth)
       lat = lat / DEG
       if (abs(lat-lat_c) <= lat_width/2) then
          amp = b_max / (1.0_8 - tanh (-slope*width/shift))
          y = (lat - (lat_c - lat_width/2))/180 * MATH_PI*radius ! y = 0 at low latitude boundary of channel                                               
          surf_geopot_case = amp * (1.0_8 - tanh (slope * (f(y) - width/shift)))
       else
          surf_geopot_case = b_max
       end if
       surf_geopot_case = grav_accel * surf_geopot_case
    else
       surf_geopot_case = grav_accel * 0d0
    end if
  end function surf_geopot_case

  real(8) function surf_geopot_latlon (lat, lon)
    ! Surface geopotential
    implicit none
    real(8)            :: lon, lat
    real(8)            :: amp, b_max, c1, cs2, sn2, rgrc, y

    ! Find latitude and longitude from Cartesian coordinates
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2

    if (trim (test_case) == "DCMIP2012c4") then
       c1 = u_0*cos((1.0_8-eta_0)*MATH_PI/2)**1.5

       surf_geopot_latlon = c1 * (c1 * (-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8)  + &
            radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
    elseif (trim (test_case) == "DCMIP2008c5") then
       rgrc = radius*acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

       surf_geopot_latlon = grav_accel*h_0*exp__flush (-rgrc**2/d2)
    elseif (trim (test_case) == "Held_Suarez") then
       c1 = u_0*cos((1.0_8-eta_0)*MATH_PI/2)**1.5
       surf_geopot_latlon = c1*(c1*(-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8) &
            + radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
    elseif (trim (test_case) == "seamount") then
       rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))
       surf_geopot_latlon = grav_accel*h0 * exp__flush (-(rgrc/width)**2)
    elseif (trim (test_case) == "upwelling") then
       b_max = abs (max_depth - min_depth)
       lat = lat / DEG
       if (abs(lat-lat_c) <= lat_width/2) then
          amp = b_max / (1.0_8 - tanh (-slope*width/shift))
          y = (lat - (lat_c - lat_width/2))/180 * MATH_PI*radius ! y = 0 at low latitude boundary of channel                                               
          surf_geopot_latlon = amp * (1.0_8 - tanh (slope * (f(y) - width/shift)))
       else
          surf_geopot_latlon = b_max
       end if
       surf_geopot_latlon = grav_accel * surf_geopot_latlon
    else
       surf_geopot_latlon = grav_accel * 0d0
    end if
  end function surf_geopot_latlon

  real(8) function f (y)
    implicit none
    real(8) :: y

    if (y <= width/2) then
       f = y
    else
       f = width - y
    end if
  end function f

  subroutine initialize_a_b_vert_case
    implicit none
    integer :: k
    real(8) :: z
    real(8), dimension(6) :: p

    ! Allocate vertical grid parameters
    if (trim(test_case) == "drake" .or. trim(test_case) == "seamount" .or. trim(test_case) == "upwelling" &
         .or. trim(test_case) == "jet") then
       allocate (a_vert(0:zlevels), b_vert(0:zlevels))
       allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
    else
       allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
       allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
    end if

    if (trim (test_case) == 'DCMIP2008c5'.or. trim (test_case) == 'DCMIP2012c4' .or. trim (test_case) == 'Held_Suarez') then
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
       elseif (zlevels==32) then
          a_vert = (/  0.00225523952394724d0, 0.00503169186413288d0, 0.0101579474285245d0, &
               0.0185553170740604d0, 0.0297346755951211d0, 0.0392730012536049d0, &
               0.0471144989132881d0, 0.0562404990196228d0, 0.0668004974722862d0, &
               0.0807014182209969d0, 0.0949410423636436d0, 0.11169321089983d0, &
               0.131401270627975d0, 0.154586806893349d0, 0.181863352656364d0, &
               0.17459799349308d0, 0.166050657629967d0, 0.155995160341263d0, 0.14416541159153d0, &
               0.130248308181763d0, 0.113875567913055d0, 0.0946138575673103d0, &
               0.0753444507718086d0, 0.0576589405536652d0, 0.0427346378564835d0, &
               0.0316426791250706d0, 0.0252212174236774d0, 0.0191967375576496d0, &
               0.0136180268600583d0, 0.00853108894079924d0, 0.00397881818935275d0, 0d0, 0d0 /)

          b_vert = (/ 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.0393548272550106d0, &
               0.0856537595391273d0, 0.140122056007385d0, 0.204201176762581d0, &
               0.279586911201477d0, 0.368274360895157d0, 0.47261056303978d0, &
               0.576988518238068d0, 0.672786951065063d0, 0.753628432750702d0, &
               0.813710987567902d0, 0.848494648933411d0, 0.881127893924713d0, &
               0.911346435546875d0, 0.938901245594025d0, 0.963559806346893d0, &
               0.985112190246582d0, 1d0 /)
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
       a_vert = a_vert(zlevels+1:1:-1) * p_0
       b_vert = b_vert(zlevels+1:1:-1)

       ! Set pressure at infinity
       p_top = a_vert(zlevels+1) ! note that b_vert at top level is 0, a_vert is small but non-zero

       ! Set mass coefficients
       a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
       b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
    elseif (trim (test_case) == "seamount") then
       if (trim (coords) == "uniform") then 
          do k = 0, zlevels
             a_vert(k) = dble(k)/dble(zlevels)
             b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
          end do
       elseif (trim (coords) == "chebyshev") then
          a_vert(0) = 0.0_8; b_vert(0) = 1.0_8
          do k = 1, zlevels-1
             a_vert(k) = (1.0_8 + cos (dble(2*k-1)/dble(2*(zlevels-1)) * MATH_PI)) / 2
             b_vert(k) = (1.0_8 + cos (dble(2*k-1)/dble(2*(zlevels-1)) * MATH_PI)) / 2
          end do
          a_vert(zlevels) = 1.0_8; b_vert(zlevels) = 0.0_8
       end if

       ! Vertical grid spacing
       a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
       b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
    elseif (trim (test_case) == "upwelling" .or. trim (test_case) == "jet") then
       b_vert(0) = 1.0_8 ; b_vert(zlevels) = 0.0_8

       if (trim (coords) == "uniform") then
          do k = 1, zlevels-1
             b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
          end do
       elseif (trim(coords) == "croco" .or. trim(coords) == "roms") then
          if (trim(coords) == "croco") then
             p = (/ -5.7831,  18.9754, -24.6521,  16.1698, -5.7092, 0.9972 /)
          elseif (trim(coords) == "roms") then
             p = (/ 0.0, 0.85230,  -3.2106,  4.4745, -3.1587,  1.0343 /)
          end if
          do k = 1, zlevels-1
             z = dble(k)/dble(zlevels)
             b_vert(k) = p(1)*z**5 + p(2)*z**4 + p(3)*z**3 + p(4)*z**2 + p(5)*z + p(6)
          end do
       end if
       a_vert = 1.0_8 - b_vert

       ! Vertical grid spacing
       a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
       b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
    else 
         do k = 1, zlevels+1
            a_vert(k) = dble(k-1)/dble(zlevels) * p_top
            b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
         end do
       a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
       b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1) 
    end if
  end subroutine initialize_a_b_vert_case

  function z_coords_case (eta_surf, z_s)
    ! Hybrid sigma-z vertical coordinates to minimize inclination of layers to geopotential
    ! near the free surface over strong bathymetry gradients.
    ! Reference: similar to Shchepetkin and McWilliams (JCP vol 228, 8985-9000, 2009)
    !
    ! Sets the a_vert parameter that depends on eta_surf (but not b_vert).
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    integer                       :: k
    real(8)                       :: cff, cff1, cff2, hc, z_0
    real(8), parameter            :: theta_b = 0d0, theta_s = 7d0
    real(8), dimension(0:zlevels) :: Cs, sc

    hc = min (abs(min_depth), abs(Tcline))

    cff1 = 1.0_8 / sinh (theta_s)
    cff2 = 0.5d0 / tanh (0.50 * theta_s)

    sc(0) = -1.0_8
    Cs(0) = -1.0_8
    cff = 1d0 / dble(zlevels)
    do k = 1, zlevels
       sc(k) = cff * dble (k - zlevels)
       Cs(k) = (1.0_8 - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5d0)) - 0.5d0)
    end do

    z_coords_case(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1.0_8 - z_0 / z_s
       z_coords_case(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords_case

  subroutine cal_AB
    ! Computes A and B coefficients for hybrid vertical grid as in LMDZ
    implicit none

    integer                         :: l
    real(8)                         :: snorm
    real(8), dimension(1:zlevels)   :: dsig
    real(8), dimension(1:zlevels+1) :: sig

    snorm  = 0.0_8
    do l = 1, zlevels
       dsig(l) = 1.0_8 + 7 * sin (MATH_PI*(l-0.5_8)/(zlevels+1))**2 ! LMDZ standard (concentrated near top and surface)
       snorm = snorm + dsig(l)
    end do

    do l = 1, zlevels
       dsig(l) = dsig(l)/snorm
    end do

    sig(zlevels+1) = 0.0_8
    do l = zlevels, 1, -1
       sig(l) = sig(l+1) + dsig(l)
    end do

    b_vert(zlevels+1) = 0.0_8
    do  l = 1, zlevels
       b_vert(l) = exp (1.0_8 - 1/sig(l)**2)
       a_vert(l) = (sig(l) - b_vert(l)) * p_0
    end do
    b_vert(1) = 1.0_8
    a_vert(1) = 0.0_8
    a_vert(zlevels+1) = (sig(zlevels+1) - b_vert(zlevels+1)) * p_0
  end subroutine cal_AB

  subroutine read_test_case_parameters
    implicit none
    integer        :: fid = 500
    real(8)        :: press_save
    character(255) :: filename, varname

    ! Find input parameters file name
    if (iargc() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'flat_projection_data.in'
    end if
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, mean_beg
    read (fid,*) varname, mean_end
    read (fid,*) varname, cp_2d
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, uniform
    read (fid,*) varname, level_save
    read (fid,*) varname, N
    read (fid,*) varname, press_save
    read (fid,*) varname, lat_val
    read (fid,*) varname, lon_val
    read (fid,*) varname, zonal
    read (fid,*) varname, NCAR_topo
    read (fid,*) varname, topo_file
    read (fid,*) varname, topo_min_level
    read (fid,*) varname, topo_max_level  
    close(fid)

    if (rank==0) then
       write (6,'(A,A)')      "test_case              = ", trim (test_case)
       write (6,'(A,A)')      "run_id                 = ", trim (run_id)
       write (6,'(A,i5)')     "number of domains      = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors   = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL           = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL            = ", PATCH_LEVEL
       write (6,'(A,i4)')     "first mean file        = ", mean_beg
       write (6,'(A,i4)')     "last mean file         = ", mean_end
       write (6,'(A,i4)')     "file to save as 2d     = ", cp_2d
       write (6,'(A,i3)')     "min_level              = ", min_level
       write (6,'(A,i3)')     "max_level              = ", max_level
       write (6,'(A,i3)')     "zlevels                = ", zlevels
       write (6,'(A,L1)')     "uniform                = ", uniform
       write (6,'(A,i3)')     "level_save             = ", level_save
       write (6,'(A,i5)')     "N                      = ", N
       write (6,'(A,es10.4)') "pressure_save (hPa)    = ", press_save
       write (6,'(A,es10.4)') "lat_val = ", lat_val
       write (6,'(A,es10.4)') "lon_val = ", lon_val
       write( 6,'(A,L1)')     "zonal average          = ", zonal
       write( 6,'(A,L1)')     "NCAR_topo              = ", NCAR_topo
       write (6,'(A,A)')      "topo_file              = ", trim (topo_file)
       write (6,'(a,i3)')     "topo_min_level       = ", topo_min_level
       write (6,'(a,i3)')     "topo_max_level       = ", topo_max_level
       write (6,*) ' '
    end if

    allocate (pressure_save(1))
    pressure_save(1) = 100*press_save

    if (climatology) call check_climatology_mean_beg_2D
  end subroutine read_test_case_parameters

  subroutine topo_flat (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer  :: d, e, id, id_e, id_i

    type(Coord)                    :: p
    type(Coord), dimension(1:EDGE) :: q

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1
    p = dom%node%elts(id_i)

    select case (itype)
    case ("bathymetry")
       topography%data(d)%elts(id_i) = max_depth + surf_geopot (d, id_i) / grav_accel
    case ("penalize")
       if (trim (test_case) == "upwelling" .or. trim (test_case) == "jet") then
          penal_node(zlev)%data(d)%elts(id_i) = mask (p)

          q(RT+1) = dom%node%elts(idx(i+1, j,   offs, dims)+1)
          q(DG+1) = dom%node%elts(idx(i+1, j+1, offs, dims)+1)
          q(UP+1) = dom%node%elts(idx(i,   j+1, offs, dims)+1)
          do e = 1, EDGE
             id_e = EDGE*id + e
             penal_edge(zlev)%data(d)%elts(id_e) = interp (mask(p), mask(q(e)))
          end do
       end if
    end select
  end subroutine topo_flat

  real(8) function mask (p)
    ! Land mass mask
    implicit none
    type(Coord) :: p

    real(8) :: dlat, lat, lon
    real(8) :: n_smth_N, n_smth_S, width_N, width_S

    dlat = 0.5d0*npts_penal * (dx_max/radius) / DEG ! widen channel to account for boundary smoothing
    width_S = 90d0 + (lat_c - (lat_width/2d0 + dlat))
    width_N = 90d0 - (lat_c + (lat_width/2d0 + dlat))

    n_smth_S = 4d0*radius * width_S*DEG / (dx_max * npts_penal)
    n_smth_N = 4d0*radius * width_N*DEG / (dx_max * npts_penal)

    call cart2sph (p, lon, lat)

    mask = exp__flush (- abs((lat/DEG+90d0)/width_S)**n_smth_S) + exp__flush (- abs((lat/DEG-90d0)/width_N)**n_smth_N)
  end function mask

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_flat (dom, i, j, zlev, offs, dims, 'bathymetry')
  end subroutine set_bathymetry

  subroutine set_penal (dom, i, j, zlev, offs, dims)
    ! Set penalization mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (penalize) then
       call topo_flat (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
    end if
  end subroutine set_penal

  subroutine apply_initial_conditions_case
    implicit none
    integer :: d, k, l

    do l = level_start, level_end
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_start, level_end
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    implicit none
    integer :: l
    
    do l = level_start, level_end
       if (NCAR_topo) call apply_onescale (assign_topo, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine update_case

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Dummy routine
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values
    use utils_mod
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k
    real(8)                       :: eta, rho, z_s
    real(8)                       :: k_T, lat, lon, p, p_s, pot_temp
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord)                   :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    x_i = dom%node%elts(id_i)

    eta = 0d0
    z_s = topography%data(d)%elts(id_i)

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if

    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zmax
       rho = porous_density (d, id_i, k)
       if (trim (test_case) == "drake") then
          if (k == zlevels+1) then
             sol_mean(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
             sol_mean(S_TEMP,k)%data(d)%elts(id_i) = 0.0_8
          else
             sol_mean(S_MASS,k)%data(d)%elts(id_i) = ref_density * height(k) 
             sol_mean(S_TEMP,k)%data(d)%elts(id_i) = (density_drake(k) - ref_density) * height(k)
          end if
       elseif (trim (test_case) == "seamount") then
          if (k == zlevels+1) then
             sol_mean(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
             sol_mean(S_TEMP,k)%data(d)%elts(id_i) = 0.0_8
          else
             if (mean_split) then
                sol_mean(S_MASS,k)%data(d)%elts(id_i) = rho * dz(k)
                sol_mean(S_TEMP,k)%data(d)%elts(id_i) = rho * dz(k) * buoy_flat (eta, z_s, k)
             else
                sol_mean(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
                sol_mean(S_TEMP,k)%data(d)%elts(id_i) = 0.0_8
             end if
          end if
       elseif (trim (test_case) == "upwelling" .or. trim (test_case) == "jet") then
          if (k == zlevels+1) then
             sol_mean(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
             sol_mean(S_TEMP,k)%data(d)%elts(id_i) = 0.0_8
          else
             sol_mean(S_MASS,k)%data(d)%elts(id_i) = rho * dz(k)
             sol_mean(S_TEMP,k)%data(d)%elts(id_i) = 0.0_8
          end if
          sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8   
       elseif (trim (test_case) == "Held_Suarez") then
          call cart2sph (x_i, lon, lat)

          dom%surf_press%elts(id+1) = surf_pressure (d, id+1)
          p_s = dom%surf_press%elts(id+1)

          p = 0.5d0 * (a_vert(k) + a_vert(k+1) + (b_vert(k) + b_vert(k+1)) * p_s) ! pressure at level k

          call cal_theta_eq (p, p_s, lat, pot_temp, k_T)

          sol_mean(S_MASS,k)%data(d)%elts(id+1) = a_vert_mass(k) + b_vert_mass(k) * p_s / grav_accel
          sol_mean(S_TEMP,k)%data(d)%elts(id+1) = sol_mean(S_MASS,k)%data(d)%elts(id+1) * pot_temp
          sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
       end if
    end do
  end subroutine init_mean

  subroutine cal_theta_eq (p, p_s, lat, theta_equil, k_T)
    ! Returns equilibrium potential temperature theta_equil and Newton cooling constant k_T
    use domain_mod
    implicit none
    real(8) :: p, p_s, lat, theta_equil, k_T

    real(8) :: cs2, sigma, theta_force, theta_tropo

    cs2 = cos (lat)**2

    sigma = (p - p_top) / (p_s - p_top)

    k_T = k_a + (k_s - k_a) * max (0d0, (sigma - sigma_b) / sigma_c) * cs2**2

    theta_tropo = T_tropo * (p / p_0)**(-kappa)  ! potential temperature at tropopause

    theta_force = T_mean - delta_T * (1d0 - cs2) - delta_theta * cs2 * log (p / p_0)

    theta_equil = max (theta_tropo, theta_force) ! equilibrium temperature
  end subroutine cal_theta_eq

  real(8) function surf_pressure (d, id) 
    ! Surface pressure
    implicit none
    integer :: d, id

    real(8) :: z_s

    if (NCAR_topo) then ! use standard atmosphere
       z_s = surf_geopot_case (d, id) / grav_accel
       call std_surf_pres (z_s, surf_pressure)
    else
       surf_pressure = p_0
    end if
  end function surf_pressure

  real(8) function buoy_flat (eta, z_s, zlev)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    integer :: zlev
    real(8) :: eta, z_s

    real(8)            :: rho, z, z1, z2
    logical, parameter :: roms = .false.

    if (trim (test_case) == "drake") then
       z = a_vert(zlev) * eta + b_vert(zlev) * z_s
       if (zlevels /= 1 .and. z >= halocline) then
          buoy_flat = - (1.0_8 - z/halocline) * drho/ref_density
       else
          buoy_flat = 0.0_8
       end if
    elseif (trim (test_case) == "seamount") then
       z1 = a_vert(zlev-1) * eta + b_vert(zlev-1) * z_s
       z2 = a_vert(zlev)   * eta + b_vert(zlev)   * z_s

       rho = 0.5 * (density_flat (z1) + density_flat (z2))

       buoy_flat = (ref_density - rho) / ref_density 
    else
       buoy_flat = 0.0_8 ! no mean component
    end if
  end function buoy_flat

  real(8) function density_flat (z)
    implicit none
    real(8) :: z

    if (trim (test_case) == "seamount") then
       if (trim(stratification) == "linear") then 
          density_flat = ref_density + drho * (max_depth-z)/max_depth
       elseif (trim(stratification) == "exponential") then
          density_flat = ref_density + drho * exp__flush (z/delta)
       end if
    elseif (trim (test_case) == "upwelling") then
       density_flat = density_eos (S_ref, temp_init (z),  z)
    else
       if (rank == 0) write(6,'(A)') "Test case not supported"
       stop
    end if
  end function density_flat

  real(8) function temp_init (z)
    implicit none
    real(8) :: z

    real(8)            :: hz, strat, z0, z1
    real(8), parameter :: h0 = 150_8, hz_0 = 6.5_8, z0_0 = -35_8, z1_0 = -75_8

    strat = abs(max_depth)
    hz = hz_0 * abs(max_depth/h0)
    z0 = z0_0 * abs(max_depth/h0)
    z1 = z1_0 * abs(max_depth/h0)

    temp_init = T_ref + 4*tanh ((z - z0) / hz) + (z - z1) / strat
  end function temp_init

  subroutine set_thresholds_case
    ! Dummy routine
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    implicit none
    real(8) :: area

    area = 4d0*MATH_PI*radius**2/(20d0*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

    area = 4d0*MATH_PI*radius**2/(20d0*4**min_level)
    dx_max = sqrt (4d0/sqrt(3d0) * area)
  end subroutine initialize_dt_viscosity_case

  real(8) function tau (p)
    ! Magnitude of wind stress at node p (dummy routine)
    implicit none
    type(Coord) :: p

    tau = 0.1_8
  end function tau

  subroutine set_save_level_case
    implicit none
    save_zlev = zlevels
  end subroutine set_save_level_case

  subroutine dump_case (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump_case

  subroutine load_case (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load_case

  subroutine init_physics_climatology
   implicit none
   integer :: k

   ! allocate climatology summation arrays
   allocate(simple_phys_temp(0:zlevels))
   allocate(simple_phys_vels(0:zlevels))
   allocate(simple_phys_zonal(0:zlevels))
   allocate(simple_phys_merid(0:zlevels))

   !Initialize arrays to zero ! Column 0 is for the sol_save, since simple_phys always has atlest zlevels 0:zmax in S_Temp
   simple_phys_temp = sol(S_TEMP, 0:zlevels)
   simple_phys_vels = sol(S_VELO, 0:zlevels)
   simple_phys_zonal = sol(S_TEMP, 0:zlevels)
   simple_phys_merid = sol(S_TEMP, 0:zlevels)
   do k = 0, zlevels
      call zero_float_field (simple_phys_temp(k), S_TEMP)
      call zero_float_field (simple_phys_vels(k), S_VELO)
      call zero_float_field (simple_phys_zonal(k), S_TEMP)
      call zero_float_field (simple_phys_merid(k), S_TEMP)
   end do
   
 end subroutine init_physics_climatology

 subroutine climatology_add_temp (dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id_i, d
   real(8) :: full_mass, full_theta, potential_temp, temperature
   ! Get id of column and domain
   id_i = idx (i, j, offs, dims) + 1
   d = dom%id + 1

   ! calulate the pressure
   call cal_pressure(dom, i, j, zlev, offs, dims)

   ! calculate the temperature
   full_mass = mass(id_i) + mean_m(id_i)
   full_theta = temp(id_i) + mean_t(id_i)
   potential_temp = full_theta/full_mass
   temperature = potential_temp * ((dom%press%elts(id_i)/dom%surf_press%elts(id_i))**kappa)

   ! sum the temperature
   temp1(id_i) = temp1(id_i) + temperature

   !calculate the density -> for the KEs
   dom%ke%elts(id_i) = dom%press%elts(id_i)/(temperature*R_d)

 end subroutine climatology_add_temp

 subroutine climatology_add_velocities (dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id, d, e

   ! Get id of column and domain
   id = idx (i, j, offs, dims)
   d = dom%id + 1

   !update the velocities
   do e = RT,UP
      velo_2d(EDGE*id+e+1) = velo_2d(EDGE*id+e+1) + velo(EDGE*id+e+1)
   end do
 end subroutine climatology_add_velocities

 subroutine climatology_add_KE(dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id_i, d, e

   ! Get id of column and domain
   id_i = idx (i, j, offs, dims) + 1
   d = dom%id + 1

   !Calculate the zonal and meridional velocities
   call interp_UVW_latlon(dom, i, j, zlev, offs, dims)


   ! !update the KE's
   simple_phys_zonal(zlev)%data(d)%elts(id_i) = simple_phys_zonal(zlev)%data(d)%elts(id_i)+(0.5*dom%ke%elts(id_i)*(velo1(id_i)**2))
   simple_phys_merid(zlev)%data(d)%elts(id_i) = simple_phys_merid(zlev)%data(d)%elts(id_i)+(0.5*dom%ke%elts(id_i)*(velo2(id_i)**2))
 end subroutine climatology_add_KE

 subroutine climatology_temp_mean(dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id_i, d

   ! Get id of column and domain
   id_i = idx (i, j, offs, dims) + 1
   d = dom%id + 1

   temp1(id_i) = temp1(id_i)/(mean_end-mean_beg+1)

 end subroutine climatology_temp_mean

 subroutine climatology_velocity_mean(dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id, d, e

   ! Get id of column and domain
   id = idx (i, j, offs, dims) 
   d = dom%id + 1

   do e = RT,UP
      velo_2d(EDGE*id+e+1) = velo_2d(EDGE*id+e+1)/(mean_end-mean_beg+1)
   end do

 end subroutine climatology_velocity_mean

 subroutine climatology_KE_mean(dom, i, j, zlev, offs, dims)
   implicit none
   type (Domain)                  :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id_i, d, e

   ! Get id of column and domain
   id_i = idx (i, j, offs, dims) + 1
   d = dom%id + 1

   simple_phys_zonal(zlev)%data(d)%elts(id_i) = simple_phys_zonal(zlev)%data(d)%elts(id_i)/(mean_end-mean_beg+1)
   simple_phys_merid(zlev)%data(d)%elts(id_i) = simple_phys_merid(zlev)%data(d)%elts(id_i)/(mean_end-mean_beg+1)

  end subroutine climatology_KE_mean

  subroutine check_climatology_mean_beg_2D
   implicit none

   if (cp_2d .ne. mean_end) then
      if (rank == 0) write(6,'(A)') "Attempting to get Simple physics climatology, but only works if mean_beg = cp_2d in input file"
      stop
   end if
  end subroutine check_climatology_mean_beg_2D

  subroutine cal_temp_dens (dom, i, j, zlev, offs, dims)
   ! Compute temperature in compressible case
   implicit none
   type(Domain)                   :: dom
   integer                        :: i, j, zlev
   integer, dimension(N_BDRY+1)   :: offs
   integer, dimension(2,N_BDRY+1) :: dims

   integer :: id, d, k
   real(8), dimension(1:zlevels) :: p

   d = dom%id + 1
   id = idx(i, j, offs, dims)

   ! Integrate the pressure upwards
   p(1) = dom%surf_press%elts(id+1) - 0.5*grav_accel*sol(S_MASS,1)%data(d)%elts(id+1)
   do k = 2, zlevels
      p(k) = p(k-1) - grav_accel*interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k-1)%data(d)%elts(id+1))
   end do

   ! Temperature at all vertical levels (saved in exner_fun) and density (saved in penal_node)
   do k = 1, zlevels
      exner_fun(k)%data(d)%elts(id+1) = sol(S_TEMP,k)%data(d)%elts(id+1)/sol(S_MASS,k)%data(d)%elts(id+1) * (p(k)/p_0)**kappa
      penal_node(k)%data(d)%elts(id+1) = p(k) / (exner_fun(k)%data(d)%elts(id+1) * R_d)
   end do

   ! Temperature at save levels (saved in trend)
   do k = 1, save_levels
      trend(1,k)%data(d)%elts(id+1) = sol_save(S_TEMP,k)%data(d)%elts(id+1)/sol_save(S_MASS,k)%data(d)%elts(id+1) * &
           (pressure_save(k)/p_0)**kappa
   end do

 end subroutine cal_temp_dens

 subroutine deallocate_climatology
   ! Deallocated extra structures
   implicit NONE
   deallocate(simple_phys_temp)
   deallocate(simple_phys_vels)
   deallocate(simple_phys_zonal)
   deallocate(simple_phys_merid)
 end subroutine deallocate_climatology

end module test_case_mod
