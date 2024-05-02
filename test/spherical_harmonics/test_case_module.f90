module test_case_mod
  ! Module file for spherical_harmonics
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none
  integer        :: angular_order, cp_beg, cp_end, lwin, N, ntaper, k_min, k_max
  real(8)        :: concentration, lat0, lon0, ref_surf_press, theta0
  real(8)        :: dPdim, R_ddim, specvoldim, dTempdim
  character(255) :: spec_type
  logical        :: local_spec
  character(255) :: coords

  ! DCMIP2012c4
  real(8) :: eta_0, u_0 
  ! DCMIP2008c5
  real(8) :: d2, h_0, lat_c, lon_c
  ! Drake
  integer               :: npts_penal
  real(8)               :: drho, halocline, radius_earth, scale
  real(8), dimension(2) :: density_drake, height
  ! Jet
  real(8) :: beta, f0, Tcline
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
    ! Additional physics for the flux term of the scalar trend
    ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
    !
    ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at southwest edges W, SW, S
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
    type(domain)                                         :: dom
    integer                                              :: d, id, idE, idNE, idN, v, zlev
    logical, optional                                    :: type

    physics_scalar_flux_case = 0d0
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    physics_velo_source_case = 0d0
  end function physics_velo_source_case

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential
    implicit none
    integer      :: d, id
    
    Type(Coord) :: x_i
    real(8)     :: c1, cs2, sn2, lon, lat, rgrc

    x_i = grid(d)%node%elts(id)

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2

    if (trim (test_case) == "DCMIP2012c4") then
       c1 = u_0*cos((1d0-eta_0)*MATH_PI/2)**1.5

       surf_geopot_case = c1 * (c1 * (-2d0*sn2**3*(cs2 + 1d0/3d0) + 10d0/63d0)  + &
            radius*omega*(8d0/5d0*cs2**1.5d0*(sn2 + 2d0/3d0) - MATH_PI/4d0))
    elseif (trim (test_case) == "DCMIP2008c5") then
       rgrc = radius*acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

       surf_geopot_case = grav_accel*h_0*exp__flush (-rgrc**2/d2)
    elseif (trim (test_case) == "Held_Suarez") then
       c1 = u_0*cos((1d0-eta_0)*MATH_PI/2)**1.5
       surf_geopot_case = c1*(c1*(-2*sn2**3*(cs2 + 1d0/3d0) + 10d0/63d0) &
            + radius*omega*(8d0/5d0*cs2**1.5*(sn2 + 2d0/3d0) - MATH_PI/4))
       ! surf_geopot_case = 0d0
    elseif (trim (test_case) == "drake" .or. trim (test_case) == "jet") then
       surf_geopot_case = 0d0
    else
       write(6,'(A)') "Test case not supported"
       call abort
    end if
  end function surf_geopot_case

  subroutine initialize_a_b_vert_case
    implicit none
    integer               :: k
    real(8)               :: z
    real(8), dimension(6) :: p

    ! Allocate vertical grid parameters
    allocate (a_vert(0:zlevels),      b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (trim (test_case) == 'DCMIP2008c5'.or. trim (test_case) == 'DCMIP2012c4' .or. trim (test_case) == 'Held_Suarez') then
       if (zlevels==18) then
          a_vert=(/0.00251499d0, 0.00710361d0, 0.01904260d0, 0.04607560d0, 0.08181860d0, &
               0.07869805d0, 0.07463175d0, 0.06955308d0, 0.06339061d0, 0.05621774d0, 0.04815296d0, &
               0.03949230d0, 0.03058456d0, 0.02193336d0, 0.01403670d0, 0.007458598d0, 0.002646866d0, &
               0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0.03756984d0, 0.08652625d0, 0.1476709d0, 0.221864d0, &
               0.308222d0, 0.4053179d0, 0.509588d0, 0.6168328d0, 0.7209891d0, 0.816061d0, 0.8952581d0, &
               0.953189d0, 0.985056d0, 1d0 /)
       elseif (zlevels==26) then
          a_vert=(/0.002194067d0, 0.004895209d0, 0.009882418d0, 0.01805201d0, 0.02983724d0, 0.04462334d0, 0.06160587d0, &
               0.07851243d0, 0.07731271d0, 0.07590131d0, 0.07424086d0, 0.07228744d0, 0.06998933d0, 0.06728574d0, 0.06410509d0, &
               0.06036322d0, 0.05596111d0, 0.05078225d0, 0.04468960d0, 0.03752191d0, 0.02908949d0, 0.02084739d0, 0.01334443d0, &
               0.00708499d0, 0.00252136d0, 0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.01505309d0, 0.03276228d0, 0.05359622d0, &
               0.07810627d0, 0.1069411d0, 0.1408637d0, 0.1807720d0, 0.2277220d0, 0.2829562d0, 0.3479364d0, 0.4243822d0, &
               0.5143168d0, 0.6201202d0, 0.7235355d0, 0.8176768d0, 0.8962153d0, 0.9534761d0, 0.9851122d0, 1d0 /)
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
          a_vert=(/0.002251865d0, 0.003983890d0, 0.006704364d0, 0.01073231d0, 0.01634233d0, 0.02367119d0, &
               0.03261456d0, 0.04274527d0, 0.05382610d0, 0.06512175d0, 0.07569850d0, 0.08454283d0, &
               0.08396310d0, 0.08334103d0, 0.08267352d0, 0.08195725d0, 0.08118866d0, 0.08036393d0, &
               0.07947895d0, 0.07852934d0, 0.07751036d0, 0.07641695d0, 0.07524368d0, 0.07398470d0, &
               0.07263375d0, 0.07118414d0, 0.06962863d0, 0.06795950d0, 0.06616846d0, 0.06424658d0, &
               0.06218433d0, 0.05997144d0, 0.05759690d0, 0.05504892d0, 0.05231483d0, 0.04938102d0, &
               0.04623292d0, 0.04285487d0, 0.03923006d0, 0.03534049d0, 0.03116681d0, 0.02668825d0, &
               0.02188257d0, 0.01676371d0, 0.01208171d0, 0.007959612d0, 0.004510297d0, 0.001831215d0, &
               0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
               0.006755112d0, 0.01400364d0, 0.02178164d0, 0.03012778d0, 0.03908356d0, 0.04869352d0, &
               0.05900542d0, 0.07007056d0, 0.08194394d0, 0.09468459d0, 0.1083559d0, 0.1230258d0, &
               0.1387673d0, 0.1556586d0, 0.1737837d0, 0.1932327d0, 0.2141024d0, 0.2364965d0, &
               0.2605264d0, 0.2863115d0, 0.3139801d0, 0.3436697d0, 0.3755280d0, 0.4097133d0, &
               0.4463958d0, 0.4857576d0, 0.5279946d0, 0.5733168d0, 0.6219495d0, 0.6741346d0, &
               0.7301315d0, 0.7897776d0, 0.8443334d0, 0.8923650d0, 0.9325572d0, 0.9637744d0, &
               0.9851122d0, 1d0/)
       else
          write(0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
          stop
       end if
       ! DCMIP order is opposite to ours
       a_vert = a_vert(zlevels:0:-1)
       b_vert = b_vert(zlevels:0:-1)
    elseif (trim (test_case) == 'drake') then
       if (zlevels == 2) then 
          a_vert(0) = 0d0; a_vert(1) = 0d0;                 a_vert(2) = 1d0
          b_vert(0) = 1d0; b_vert(1) = halocline/max_depth; b_vert(2) = 0d0
       else
          if (trim (coords) == "uniform") then 
             do k = 0, zlevels
                b_vert(k) = 1d0 - dble(k)/dble(zlevels)
             end do
          elseif (trim (coords) == "chebyshev") then
             do k = 0, zlevels
                b_vert(k) = (1d0 + cos (dble(k)/dble(zlevels) * MATH_PI)) / 2d0
             end do
          elseif (trim (coords) == "chebyshev_half") then
             do k = 0, zlevels
                b_vert(k) = 1d0 - sin (dble(k)/dble(zlevels) * MATH_PI/2d0)
             end do
          else
             write (6,*) "Selected vertical coordinate type not supported ..."
             call abort
          end if
          a_vert = 1d0 - b_vert
       end if
    elseif (trim (test_case) == 'jet') then
       a_vert = 0d0 ; b_vert = 0d0
       if (trim (coords) == "uniform") then
          do k = 2, zlevels-1
             b_vert(k) = 1d0 - dble(k)/dble(zlevels)
          end do
       elseif (trim(coords) == "croco") then
          p = (/ -5.7831,  18.9754, -24.6521,  16.1698, -5.7092, 0.9972 /)
          do k = 1, zlevels-1
             z = dble(k)/dble(zlevels)
             b_vert(k) = p(1)*z**5 + p(2)*z**4 + p(3)*z**3 + p(4)*z**2 + p(5)*z + p(6)
          end do
       end if
       a_vert = 1d0 - b_vert
    end if
    
    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

  subroutine read_test_case_parameters
    implicit none
    integer        :: fid = 500
    character(255) :: filename, varname

    ! Find input parameters file name
    if (iargc() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'spherical_harmonics.in'
    end if
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, spec_type
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, cp_beg
    read (fid,*) varname, cp_end
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, vert_diffuse
    read (fid,*) varname, save_zlev
    read (fid,*) varname, uniform
    read (fid,*) varname, level_fill
    read (fid,*) varname, N
    read (fid,*) varname, local_spec
    read (fid,*) varname, lat0
    read (fid,*) varname, lon0
    read (fid,*) varname, theta0
    read (fid,*) varname, concentration
    read (fid,*) varname, ntaper
    read (fid,*) varname, angular_order
    close(fid)

    if (save_zlev > zlevels) save_zlev = zlevels ! avoid incorrect choice of save_zlev

    if (save_zlev == -1) then
       k_min = 1 ; k_max = zlevels
    else
       k_min = save_zlev; k_max = save_zlev
    end if
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A,A)')      "spec_type              = ", trim (spec_type)
       write (6,'(A,A)')      "test_case              = ", trim (test_case)
       write (6,'(A,A)')      "run_id                 = ", trim (run_id)
       write (6,'(A,i5)')     "number of domains      = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors   = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL           = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL            = ", PATCH_LEVEL
       write (6,'(A,i4)')     "First checkpoint       = ", cp_beg
       write (6,'(A,i4)')     "Last checkpoint        = ", cp_end
       write (6,'(A,i3)')     "min_level              = ", min_level
       write (6,'(A,i3)')     "max_level              = ", max_level
       write (6,'(A,i3)')     "zlevels                = ", zlevels
       write (6,'(A,l1)')     "vert_diffuse           = ", vert_diffuse
       write (6,'(A,i3)')     "save_zlev              = ", save_zlev
       write(6,'(A,l1)')      "uniform                = ", uniform
       write (6,'(A,i3)')     "level_fill             = ", level_fill
       write (6,'(A,i5)')     "N                      = ", N
       write (6,'(a,l1)')     "local_spec             = ", local_spec
       write (6,'(A,es11.4)') "lat0                   = ", lat0
       write (6,'(A,es11.4)') "lon0                   = ", lon0
       write (6,'(A,es10.4)') "theta0                 = ", theta0
       write (6,'(A,es10.4)') "concentration          = ", concentration
       write (6,'(A,i3)')     "ntaper                 = ", ntaper
       write (6,'(A,i3)')     "angular_order          = ", angular_order
       write (6,*) ' '
       call print_density_pert
    end if
  end subroutine print_test_case_parameters
  
  subroutine print_density_pert
    implicit none
    integer     :: k
    real(8)     :: z
    type(Coord) :: x_i 

    x_i = Coord (radius, 0d0, 0d0)

    write (6,'(a)') " Layer    z"      
    do k = 1, zlevels
       z = 0.5d0 * (b_vert(k)+b_vert(k-1)) * max_depth
       write (6, '(2x,i2, 3x, es9.2)') k, z
    end do
    write (6,'(A)') &
         '*********************************************************************&
         *************************************************************'
  end subroutine print_density_pert

  subroutine topo_sphere (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: npts
    character(*)                   :: itype

    integer            :: d, e, id, id_e, id_i, ii, is0, it0, jj, s, t
    real(8)            :: dx, lat, lat0, lat_width, lon, mask, M_topo, n_lat, n_lon, r, s0, t0, sw_topo, topo_sum, wgt
    real(8), parameter :: lat_max = 60, lat_min = -35, lon_width = 15
    type(Coord)        :: p, q

    id = idx (i, j, offs, dims)
    id_i = id + 1

    select case (itype)
    case ("bathymetry")
       topography%data(d)%elts(id_i) = max_depth + surf_geopot_case (d, id_i) / grav_accel
    case ("penalize")
       call cart2sph (dom%node%elts(id_i), lon, lat)
       dx = dx_max

       ! Analytic land mass with smoothing
       lat_width = (lat_max - lat_min) / 2
       lat0 = lat_max - lat_width

       n_lat = 4*radius * lat_width*DEG / (dx * npts_penal)
       n_lon = 4*radius * lon_width*DEG / (dx * npts_penal)

       mask = exp__flush (- abs((lat/DEG-lat0)/lat_width)**n_lat - abs(lon/DEG/(lon_width))**n_lon) ! constant longitude width

       d = dom%id + 1
       penal_node(zlev)%data(d)%elts(id_i) = mask
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = max (penal_edge(zlev)%data(d)%elts(id_e), mask)
       end do
    end select
  end subroutine topo_sphere

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_sphere (dom, i, j, zlev, offs, dims, 'bathymetry')
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
       call topo_sphere (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0       
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
       do k = 1, zmax
          call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! dummy routine

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
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id, id_i
    real (8)    :: dz, eta_surf, rho, z
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    if (trim (test_case) == "drake") then
       x_i  = dom%node%elts(id_i)
       eta_surf = 0d0

       if (zlev == zlevels+1) then
          sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
       else
          dz =  b_vert_mass(zlev) * max_depth
          z = 0.5d0 * (b_vert(zlev) + b_vert(zlev-1)) * topography%data(d)%elts(id_i)

          rho = porous_density (d, id_i, zlev)
          sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = rho * dz
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy_init (x_i, z)
       end if
    else
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end if
  end subroutine init_mean

  real(8) function buoyancy_init (x_i, z)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    if (zlevels /= 1 .and. z >= halocline) then
       buoyancy_init = - (1d0 - z/halocline) * drho/ref_density
    else
       buoyancy_init = 0d0
    end if
  end function buoyancy_init

  subroutine set_thresholds_case
    ! Dummy routine
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    allocate (threshold(1:N_VARIABLE,1:zlevels));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,1:zlevels)); threshold_def = 0d0
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    implicit none
  end subroutine initialize_dt_viscosity_case

  subroutine set_save_level_case
    implicit none
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

    if (trim (test_case) == "jet") then
       hc = min (abs(min_depth), abs(Tcline))

       cff1 = 1d0 / sinh (theta_s)
       cff2 = 0.5d0 / tanh (0.5d0 * theta_s)

       sc(0) = -1d0
       Cs(0) = -1d0
       cff = 1d0 / dble(zlevels)
       do k = 1, zlevels
          sc(k) = cff * dble (k - zlevels)
          Cs(k) = (1d0 - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5d0)) - 0.5d0)
       end do

       z_coords_case(0) = z_s
       do k = 1, zlevels
          cff = hc * (sc(k) - Cs(k))
          z_0 = cff - Cs(k) * z_s
          a_vert(k) = 1d0 - z_0 / z_s
          z_coords_case(k) = eta_surf * a_vert(k) + z_0
       end do
    else
       z_coords_case = 0d0
    end if
  end function z_coords_case
end module test_case_mod
