module test_case_mod
  ! Module file for climate test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use std_atm_profile_mod
  use io_mod
  use sso_mod
  use physics_trend_mod
  implicit none

  ! Standard variables
  integer :: CP_EVERY, resume_init
  real(8) :: time_start, total_cpu_time

  ! Test case variables
  real(8) :: Area_max, Area_min, C_div, dt_max
  real(8) :: topo_Area_min, topo_dx_min
  real(8) :: cfl_max, cfl_min, T_cfl, nu_sclr, nu_rotu, nu_divu, T_0, u_0

  ! Model parameters
  real(8)        :: e_thick       = 10d0 * KM ! Ekman initial conditions
  character(255) :: analytic_topo = "none"    ! mountains or none (used if NCAR_topo = .false.)

  ! From simple physics
  logical        :: Ekman_ic      = .false.   
  logical        :: scale_aware   = .false.
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
    ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
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
       physics_scalar_flux_case = 0d0
    else
       if (.not.local_type) then ! usual flux at edges E, NE, N
          l_e =  dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
          d_e =  dom%len%elts(EDGE*id+1:EDGE*id_i)
       else ! flux at SW corner
          l_e(RT+1) = dom%pedlen%elts(EDGE*idE+RT+1)
          l_e(DG+1) = dom%pedlen%elts(EDGE*idNE+DG+1)
          l_e(UP+1) = dom%pedlen%elts(EDGE*idN+UP+1)
          d_e(RT+1) =  - dom%len%elts(EDGE*idE+RT+1)
          d_e(DG+1) =  - dom%len%elts(EDGE*idNE+DG+1)
          d_e(UP+1) =  - dom%len%elts(EDGE*idN+UP+1)
       end if

       ! Calculate gradients
       if (Laplace_order == 1) then
          grad = grad_physics (q(v,zlev)%data(d)%elts)
       elseif (Laplace_order == 2) then
          grad = grad_physics (Laplacian_scalar(v)%data(d)%elts)
       end if
       physics_scalar_flux_case = (-1d0)**Laplace_order * grad * l_e
    end if
  contains
    function grad_physics (scalar)
      implicit none
      real(8), dimension(1:EDGE) :: grad_physics
      real(8), dimension(:)      :: scalar

      grad_physics(RT+1) = (visc(idE) * scalar(idE+1) - visc(id)   * scalar(id+1))   / d_e(RT+1)
      grad_physics(DG+1) = (visc(id)  * scalar(id+1)  - visc(idNE) * scalar(idNE+1)) / d_e(DG+1)
      grad_physics(UP+1) = (visc(idN) * scalar(idN+1) - visc(id)   * scalar(id+1))   / d_e(UP+1)
    end function grad_physics
    
    real(8) function visc (id)
      ! Scale aware viscosity
      ! factor 1.5 ensures that maximum stable C_visc matches theoretical estimate of 1/6^Laplace_order
      implicit none
      integer :: id
      real(8) :: Area

      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         Area = 1d0 / dom%areas%elts(id+1)%hex_inv
      else
         Area = Area_min
      end if
      
      visc = C_visc(v) * nu_scale (Area, dt)
    end function visc
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    ! Additional physics for the source term of the velocity trend
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8), dimension(1:EDGE) :: blocking_drag, wave_drag

    id = idx (i, j, offs, dims)

    ! Scale aware viscosity
    if (Laplace_order == 0) then
       physics_velo_source_case = 0d0
    else
       physics_velo_source_case = (-1d0)**(Laplace_order-1) * (grad_divu () - curl_rotu ()) 
    end if
  contains
    function grad_divu ()
      implicit none
      real(8), dimension(1:EDGE) :: grad_divu

      integer :: idE, idN, idNE
      
      idE  = idx (i+1, j,   offs, dims)
      idN  = idx (i,   j+1, offs, dims)
      idNE = idx (i+1, j+1, offs, dims)

      grad_divu(RT+1) = (visc_div(idE) * divu(idE+1) - visc_div(id)   * divu(id+1))   / dom%len%elts(EDGE*id+RT+1)
      grad_divu(DG+1) = (visc_div(id)  * divu(id+1)  - visc_div(idNE) * divu(idNE+1)) / dom%len%elts(EDGE*id+DG+1)
      grad_divu(UP+1) = (visc_div(idN) * divu(idN+1) - visc_div(id)   * divu(id+1))   / dom%len%elts(EDGE*id+UP+1)
    end function grad_divu

    function curl_rotu ()
      implicit none
      real(8), dimension(1:EDGE) :: curl_rotu

      integer :: idS, idW

      idS  = idx (i,   j-1, offs, dims)
      idW  = idx (i-1, j,   offs, dims)

      curl_rotu(RT+1) = (visc_rot(id)  * vort(TRIAG*id +LORT+1) - visc_rot(idS) * vort(TRIAG*idS+UPLT+1)) &
           / dom%pedlen%elts(EDGE*id+RT+1)
      curl_rotu(DG+1) = (visc_rot(id)  * vort(TRIAG*id +LORT+1) - visc_rot(id)  * vort(TRIAG*id +UPLT+1)) &
           / dom%pedlen%elts(EDGE*id+DG+1)
      curl_rotu(UP+1) = (visc_rot(idW) * vort(TRIAG*idW+LORT+1) - visc_rot(id)  * vort(TRIAG*id +UPLT+1)) &
           / dom%pedlen%elts(EDGE*id+UP+1)
    end function curl_rotu
    
    real(8) function visc_div (id)
      ! Scale aware viscosity
      implicit none
      integer :: id
      real(8) :: Area

      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         Area = 1d0 / dom%areas%elts(id+1)%hex_inv
      else
         Area = Area_min
      end if
      
      visc_div = C_div * nu_scale (Area, dt)
    end function visc_div

    real(8) function visc_rot (id)
      ! Scale aware viscosity
      implicit none
      integer :: id
      real(8) :: Area
      
      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         Area = 1d0 / dom%areas%elts(id+1)%hex_inv
      else
         Area = Area_min
      end if
      
      visc_rot = C_visc(S_VELO) * nu_scale (Area, dt)
    end function visc_rot
  end function physics_velo_source_case

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8) :: k_T, lat, lon, p, p_s
    
    d   = dom%id+1
    id  = idx (i, j, offs, dims)

    call cart2sph (dom%node%elts(id+1), lon, lat)

    if (NCAR_topo) then ! surface pressure from multilevel topography
       p_s = dom%surf_press%elts(id+1)
    else                ! surface pressure from standard atmosphere
       call std_surf_pres (topography%data(d)%elts(id+1), p_s)
    end if
    
    do k = 1, zlevels
       p = 0.5d0 * (a_vert(k) + a_vert(k+1) + (b_vert(k) + b_vert(k+1)) * p_s) ! pressure at level k

       if (split_mean_perturbation) then
          sol(S_MASS,k)%data(d)%elts(id+1) = 0d0
          sol(S_TEMP,k)%data(d)%elts(id+1) = 0d0
       else
          sol(S_MASS,k)%data(d)%elts(id+1) = a_vert_mass(k) + b_vert_mass(k) * p_s / grav_accel
          sol(S_TEMP,k)%data(d)%elts(id+1) = sol(S_MASS,k)%data(d)%elts(id+1) * theta_init (p)
       end if

       ! Initial velocity with Ekman layer velocity
       if (ekman_ic) then
          call vel2uvw (dom, i, j, k, offs, dims)
       else
          sol(S_VELO,k)%data(d)%elts(id_edge(id)) = 0d0
       end if
    end do
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: id, d, k
    real(8) :: k_T, lat, lon, p, p_s
    
    d  = dom%id+1
    id = idx (i, j, offs, dims)

    call cart2sph (dom%node%elts(id+1), lon, lat)

    if (NCAR_topo) then ! surface pressure from multilevel topography
       p_s = dom%surf_press%elts(id+1)
    else                ! surface pressure from standard atmosphere
       call std_surf_pres (topography%data(d)%elts(id+1), p_s)
    end if

    do k = 1, zlevels
       p = 0.5d0 * (a_vert(k) + a_vert(k+1) + (b_vert(k) + b_vert(k+1)) * p_s) ! pressure at level k

       if (split_mean_perturbation) then
          sol_mean(S_MASS,k)%data(d)%elts(id+1) = a_vert_mass(k) + b_vert_mass(k) * p_s / grav_accel
          sol_mean(S_TEMP,k)%data(d)%elts(id+1) = sol_mean(S_MASS,k)%data(d)%elts(id+1) * theta_init (p)
       else
          sol_mean(S_MASS,k)%data(d)%elts(id+1) = 0d0
          sol_mean(S_TEMP,k)%data(d)%elts(id+1) = 0d0
       end if
       sol_mean(S_VELO,k)%data(d)%elts(id_edge(id)) = 0d0
    end do
  end subroutine init_mean
  
  real(8) function theta_init (p)
    ! Initial potential temperature profile
    implicit none
    real(8) :: p

    theta_init = T_0 * (p/p_0)**(-kappa)
  end function theta_init

  subroutine vel2uvw (dom, i, j, zlev, offs, dims)
    ! Sets the velocities on the computational grid given a function vel_fun that provides zonal and meridional velocities
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    external                        :: vel_fun

    integer      :: d, id, idE, idN, idNE
    type (Coord) :: vel, x_i, x_E, x_N, x_NE

    d = dom%id+1

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id  +1)
    x_E  = dom%node%elts(idE +1)
    x_N  = dom%node%elts(idN +1)
    x_NE = dom%node%elts(idNE+1)

    vel = vel_init ()
    
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = inner (direction (x_i,  x_E), vel)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = inner (direction (x_NE, x_i), vel)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = inner (direction (x_i,  x_N), vel)
  contains    
    function vel_init ()
      ! Zonal latitude-dependent wind
      implicit none
      type(Coord) :: vel_init

      real(8)     :: D_e, lon, lat, p, p_s, phi, u, v
      type(Coord) :: e_zonal, e_merid

      call cart2sph (x_i, lon, lat)

      if (NCAR_topo) then ! surface pressure from multilevel topography
         p_s = dom%surf_press%elts(id+1)
      else                ! surface pressure from standard atmosphere
         call std_surf_pres (topography%data(d)%elts(id+1), p_s)
      end if

      p = 0.5d0 * (a_vert(zlev) + a_vert(zlev+1) + (b_vert(zlev) + b_vert(zlev+1)) * p_s) ! pressure at level k

      phi = - R_d * T_0 * log (p / p_s)
      
      D_e = e_thick * grav_accel

      u = u_0 * (1d0 - exp (-phi/D_e) * cos (phi/D_e)) * cos (lat) ! zonal velocity 
      v = u_0 * (1d0 - exp (-phi/D_e) * sin (phi/D_e)) * cos (lat) ! meridional velocity

      e_zonal = Coord (-sin(lon),           cos(lon),               0d0) 
      e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) 

      ! Velocity vector
      vel_init = u * e_zonal + v * e_merid
    end function vel_init
  end subroutine vel2uvw

  real(8) function surf_geopot_case (d, id)
    ! Set geopotential and topography
    implicit none
    integer :: d, id
    
    surf_geopot_case = grav_accel * topography%data(d)%elts(id)
  end function surf_geopot_case

  subroutine init_topo (dom, i, j, zlev, offs, dims)
    ! Assigns analytic topography
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: lat, lon, width

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    select case (analytic_topo)
    case ("mountains")
       width = 8d0 * dx_min
       call cart2sph (grid(d)%node%elts(id), lon, lat)
       
       topography%data(d)%elts(id) = &
            tanh_profile (4d3*METRE,  65*DEG,   95*DEG,  25*DEG,  40*DEG) + & ! Himalayas
            tanh_profile (4d3*METRE, -60*DEG,  -50*DEG, -60*DEG, -10*DEG)     ! Andes
    case ("none")
       topography%data(d)%elts(id) = 0d0
    end select
  contains
    real(8) function tanh_profile (height, lon_min, lon_max, lat_min, lat_max)
      implicit none
      real(8), intent(in) :: height, lon_min, lon_max, lat_min, lat_max

      tanh_profile = height * (profile1d (lat, lat_min, lat_max) * profile1d (lon, lon_min, lon_max))
    end function tanh_profile

    real(8) function profile1d (x, xmin, xmax)
      implicit none
      real(8) :: x, xmin, xmax

      profile1d = prof (x, xmax) - prof (x, xmin)
    end function profile1d

    real(8) function prof (x, x0)
      implicit none
      real(8) :: x, x0

      prof = 0.5d0 * (1d0 - tanh ((x - x0)/((width/5d0)/radius)))
    end function prof
    
    real(8) function ellipse_profile (lon_0, lat_0, e, height, sigma, theta)
      ! Elliptical smoothed mountain, ellipticity 0 < e <= 1
      ! Minimum resolution of gradient = N
      implicit none
      real(8), intent(in) :: lon_0, lat_0 ! centre of ellipse (in radians)
      real(8), intent(in) :: e            ! ellipticity
      real(8), intent(in) :: height       ! height of ellipse (in metres)
      real(8), intent(in) :: sigma        ! size of ellipse (in radians)
      real(8), intent(in) :: theta        ! orientation (in radians)

      real(8)            :: dtheta, p, lat_loc, lat_rot, lon_loc, lon_rot, rsq, sigma_x, sigma_y

      real(8), parameter :: npts_slope = 5d0 ! resolve slope with this many cells

      dtheta = dx_min / radius

      sigma_x = sigma
      sigma_y = sigma_x * sqrt (1d0 - e**2)

      ! Transform coordinates (shift and rotate)
      lon_loc = lon - lon_0; lat_loc = lat - lat_0

      lon_rot = lon_loc * cos (theta) - lat_loc * sin (theta)
      lat_rot = lon_loc * sin (theta) + lat_loc * cos (theta)

      rsq = (lon_rot/sigma_x)**2 + (lat_rot/sigma_y)**2

      ! Order of hyper Gaussian is 2 p
      p = (log (0.01d0) / log (1d0 - npts_slope * dtheta/(2d0*sigma_y))) / 2d0

      ellipse_profile = height * exp__flush (-rsq**p)
    end function ellipse_profile
  end subroutine init_topo

  subroutine vel_fun (lon, lat, u, v)
    ! Random initial wind
    implicit none
    real(8) :: lon, lat, u, v

    real(8) :: rgrc
    real(8) :: lat_c, lon_c, r
    real(8) :: amp = 1d0 ! amplitude of random noise

    ! Zonal velocity component
    call random_number (r)
    u = amp * 2d0 * (r - 0.5d0)

    ! Meridional velocity component
    call random_number (r)
    v = amp * 2d0 * (r - 0.5d0)
  end subroutine vel_fun

  subroutine initialize_a_b_vert_case
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = dble(k-1)/dble(zlevels) * p_top
          b_vert(k) = 1d0 - dble(k-1)/dble(zlevels)
       end do
    else
       if (zlevels == 18) then
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
          if (rank == 0) then
             write (6,'(/,a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             write (6,'(a)'  ) "!                                                  !"
             write (6,'(a)'  ) "!    zlevels choice not supported ... aborting     !"
             write (6,'(a)'  ) "!                                                  !"
             write (6,'(a,/)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
             call abort
          end if
       end if
       a_vert = a_vert(zlevels+1:1:-1) * p_0
       b_vert = b_vert(zlevels+1:1:-1)
    end if

    p_top = a_vert(zlevels+1) ! assumes b_vert(zlevels+1) = 0

    ! Set mass coefficients
    a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
    b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
  end subroutine initialize_a_b_vert_case

  subroutine read_test_case_parameters
    implicit none
    integer            :: k, v
    integer, parameter :: fid = 500, funit = 400
    character(255)     :: bash_cmd, command, filename, varname
    character(2)       :: var_file
    logical            :: file_exists

    ! Find input parameters file name
    if (command_argument_count() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'test_case.in'
    end if
    if (rank == 0) write (6,'(/,a,a,/)') "Input file = ", trim (filename)

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, physics_type
    read (fid,*) varname, run_id
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, Nsoil
    read (fid,*) varname, NCAR_topo
    read (fid,*) varname, sso
    read (fid,*) varname, topo_file
    read (fid,*) varname, topo_min_level
    read (fid,*) varname, topo_max_level   
    read (fid,*) varname, tol
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)

    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init

    ! Bins for zonal statistics
    !nbins = sqrt (10*4**max_level/2) ! consistent with maximum resolution
    nbins = 300
    allocate (Nstats(zlevels,nbins), Nstats_glo(zlevels,nbins)) ; Nstats = 0 ; Nstats_glo = 0
    allocate (zonal_avg(zlevels,nbins,nvar_zonal), zonal_avg_glo(zlevels,nbins,nvar_zonal))
    zonal_avg = 0d0; zonal_avg_glo = 0d0
    allocate (bounds(1:nbins-1))
    dbin = 1.8d2/nbins
    bounds = -90+dbin + dbin*(/ (ibin, ibin = 0, nbins-1) /)

    ! Initialize rank 0 with saved statistics data if present
    if (rank == 0) then
       inquire (file = trim(run_id)//'.3.tgz', exist = file_exists)
       if (file_exists) then
          command = 'gtar xzf '//trim(run_id)//'.3.tgz'
          write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
          call system (trim(bash_cmd))

          write (var_file, '(i2.2)') 00
          open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="UNFORMATTED", action='READ')
          read (funit) Nstats
          close (funit)

          do v = 1, nvar_zonal
             write (var_file, '(i2)') v+10
             open (unit=funit, file=trim(run_id)//'.3.'//var_file, form="FORMATTED", action='READ')
             do k = zlevels, 1, -1
                read (funit,*) zonal_avg(k,:,v)
             end do
             close (funit)
          end do

          ! Convert variances to sums of squares for statistics computation
          zonal_avg(:,:,2) = zonal_avg(:,:,2) * (Nstats - 1)
          do v = 6, nvar_zonal
             zonal_avg(:,:,v) = zonal_avg(:,:,v) * (Nstats - 1)
          end do
       end if
    end if
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    if (rank==0) then
       write (6,'(a)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(a)')        "RUN PARAMETERS"
       write (6,'(a,a)')      "test_case               = ", trim (test_case)
       write (6,'(a,a)')      "physics_type            = ", trim (physics_type)
       write (6,'(a,a)')      "run_id                  = ", trim (run_id)
       write (6,'(a,l1)')     "compressible            = ", compressible
       write (6,'(a,l1)')     "split_mean_perturbation = ", split_mean_perturbation 
       write (6,'(a,i3)')     "min_level               = ", min_level
       write (6,'(a,i3)')     "max_level               = ", max_level
       write (6,'(a,i5)')     "number of domains       = ", N_GLO_DOMAIN
       write (6,'(a,i5)')     "number of processors    = ", n_process
       write (6,'(a,i5)')     "DOMAIN_LEVEL            = ", DOMAIN_LEVEL
       write (6,'(a,i5)')     "PATCH_LEVEL             = ", PATCH_LEVEL
       write (6,'(a,i3)')     "zlevels                 = ", zlevels
       write (6,'(a,i3)')     "Nsoil                   = ", Nsoil
       write (6,'(a,l1)')     "uniform                 = ", uniform
       write (6,'(a,l1)')     "remap                   = ", remap
       write (6,'(a,i3)')     "iremap                  = ", iremap
       write (6,'(a,l1)')     "default_thresholds      = ", default_thresholds
       write (6,'(a,es10.4)') "tolerance               = ", tol
       write (6,'(a,i1)')     "optimize_grid           = ", optimize_grid
       write (6,'(a,l1)')     "adapt_dt                = ", adapt_dt
       write (6,'(a,es10.4)') "dt_init                 = ", dt_init
       write (6,'(a,es10.4)') "cfl_num                 = ", cfl_num
       write (6,'(a,a)')      "timeint_type            = ", trim (timeint_type)
       
       write (6,'(a,i1,/)')     "Laplace_order           = ", Laplace_order
       write (6,'(a,/,a,/,/,a,es8.2,/,a,es8.2,/)') "Stability limits:", &
            "[Klemp 2017 Damping Characteristics of Horizontal Laplacian Diffusion Filters Mon Weather Rev 145, 4365-4379.]", &
            "C_visc(S_MASS) and C_visc(S_TEMP) <  (1/6)**Laplace_order = ", (1d0/6d0)**Laplace_order, &
            "                   C_visc(S_VELO) < (1/24)**Laplace_order = ", (1d0/24d0)**Laplace_order
       if (scale_aware) then
          write (6,'(a,/)') "Scale-aware horizontal viscosity"
       else
          write (6,'(a,/)') "Horizontal viscosity based on dx_min"
       end if

       write (6,'(a)') "Non-dimensional viscosities"
       write (6,'(3(a,es8.2/))') "C_visc(S_MASS)           = ", C_visc(S_MASS), "C_visc(S_TEMP)           = ", C_visc(S_TEMP), &
            "C_visc(S_VELO)           = ", C_visc(S_VELO)

       write (6,'(a)')        "Approximate viscosities on finest grid"
       write (6,'(a,es8.2)') "nu_scalar                = ", nu_sclr
       write (6,'(a,es8.2)') "nu_rot                   = ", nu_rotu
       write (6,'(a,es8.2,/)') "nu_div                   = ", nu_divu

       write (6,'(a,es10.4)') "dt_max           [s]     = ", dt_max
       write (6,'(a,es10.4)') "dt_write         [d]     = ", dt_write / DAY
       write (6,'(a,i3)')     "CP_EVERY                 = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance                = ", rebalance
       write (6,'(a,es10.4)') "time_end         [d]     = ", time_end / DAY
       write (6,'(a,i6)')     "resume                   = ", resume_init

       write (6,'(/,a)')      "STANDARD PARAMETERS"
       write (6,'(a,es10.4)') "radius          [km]     = ", radius / KM
       write (6,'(a,es10.4)') "omega        [rad/s]     = ", omega
       write (6,'(a,es10.4)') "ref_density [kg/m^3]     = ", ref_density
       write (6,'(a,es10.4)') "p_0           [hPa]      = ", p_0/100d0
       write (6,'(a,es10.4)') "p_top         [hPa]      = ", p_top/100d0
       write (6,'(a,es10.4)') "R_d      [J/(kg K)]      = ", R_d
       write (6,'(a,es10.4)') "c_p      [J/(kg K)]      = ", c_p
       write (6,'(a,es10.4)') "c_v      [J/(kg K)]      = ", c_v
       write (6,'(a,es10.4)') "gamma                    = ", gamma
       write (6,'(a,es10.4)') "kappa                    = ", kappa
       write (6,'(a,f10.1)')  "dx_max         [km]      = ", dx_max / KM
       write (6,'(a,f10.1)')  "dx_min         [km]      = ", dx_min / KM

       write (6,'(/,a)')      "TEST CASE PARAMETERS"
       write (6,'(a,f5.1)')   "T_0             [K]      = ", T_0
       write (6,'(a,es10.4)') "T_mean          [K]      = ", T_mean
       write (6,'(a,es10.4)') "T_tropo         [K]      = ", T_tropo
       write (6,'(a,es10.4)') "sigma_b                  = ", sigma_b
       write (6,'(a,es10.4)') "k_a           [1/d]      = ", k_a / DAY
       write (6,'(a,es10.4)') "k_f           [1/d]      = ", k_f / DAY
       write (6,'(a,es10.4)') "k_s           [1/d]      = ", k_s / DAY
       write (6,'(a,es10.4)') "delta_T       [K/m]      = ", delta_T
       write (6,'(a,es10.4)') "delta_theta   [K/m]      = ", delta_theta
       write (6,'(a,es10.4,/)') "wave_speed    [m/s]      = ", wave_speed
       write (6,'(a,l)')      "NCAR_topo                = ", NCAR_topo
       if (NCAR_topo) then
          write (6,'(a,l)')      "sso                      = ", sso
          write (6,'(a,l)')      "wave_drag                = ", wave_drag
          write (6,'(a,l)')      "blocking_drag            = ", blocking_drag
          write (6,'(a,a)')      "topo_file                = ", trim (topo_file)
          write (6,'(a,i3)')     "topo_min_level           = ", topo_min_level
          write (6,'(a,i3)')     "topo_max_level           = ", topo_max_level
       else
          write (6,'(a,a)')      "analytic_topo            = ", analytic_topo
       end if
       write (6,'(a)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load, total_layers
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

    total_layers = size (threshold, 2)

    if (rank == 0) then
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i12,2(a,es8.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', sum (threshold(S_MASS,:))/total_layers, &
            ' temp tol = ', sum (threshold(S_TEMP,:))/total_layers, &
            ' velo tol = ', sum (threshold(S_VELO,:))/total_layers, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(5(es15.9,1x),i2,1x,i12,1x,2(es15.9,1x))')  &
            time/DAY, dt, sum (threshold(S_MASS,:))/total_layers, sum (threshold(S_TEMP,:))/total_layers, &
            sum (threshold(S_VELO,:))/total_layers, level_end, sum (n_active), rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    use lnorms_mod
    implicit none

    call cal_lnorm ("2")
    lnorm(S_VELO,:) = Udim

    threshold_def = tol * lnorm
  end subroutine initialize_thresholds_case

  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    implicit none

    if (default_thresholds) then
       threshold = threshold_def
    else
       call cal_lnorm ("2")
       threshold = tol * lnorm
       threshold(S_VELO,:) = max (threshold(S_VELO,:), threshold_def(S_VELO,:))
    end if
  end subroutine set_thresholds_case

  subroutine initialize_dt_viscosity_case
    ! Evaluate viscosity (for finest grid), find equivalent non-dimensional viscosities C_visc and set time step
    ! (based on CAM 120 km resolution value)
    implicit none
    real(8) :: area_sphere, dx_scaling, nu, nu_dim, C_max, nu_scaling
    
    real(8), parameter :: nu_CAM       = 1d15 * METRE**4/SECOND     ! CAM value for horizontal resolution dx = 120 km
    real(8), parameter :: res_scaling  = (136d0*KM) / (100d0*KM)    ! ratio between TRiSK and CAM grid resolutions
                                                                    ! (TRiSK grid: 107 km <= dx <= 136 km at J6)

    dx_scaling  = 2d0 ** (dble (6 - max_level))                     ! scaling factor compared to approximately J6 base CAM value
    C_max       = 1d0/6d0**Laplace_order                            ! maximum stable non-dimensional viscosity for scalars and div u

    ! Average hexagon areas and horizontal resolution
    area_sphere = 4d0*MATH_PI * radius**2 
    Area_min    = area_sphere / (10d0 * 4d0**max_level)
    Area_max    = area_sphere / (10d0 * 4d0**min_level)
    dx_min      = sqrt (2d0 / sqrt(3d0) * Area_min)              
    dx_max      = sqrt (2d0 / sqrt(3d0) * Area_max)

    ! Time step parameters
    if (adapt_dt) then
       cfl_max = 1d0                                                ! maximum cfl number
       cfl_min = cfl_max                                            ! minimum cfl number
       T_cfl   = 1d-1 * DAY                                         ! time over which to increase cfl number from cfl_min to cfl_max
       cfl_num = cfl_min                                            ! initialize cfl number
       dt_init = cfl_min * 0.85d0 * dx_min / (wave_speed + u_0)     ! initial time step     (0.85 factor corrects for minimum dx)
       dt_max  = cfl_max * 0.85d0 * dx_min / (wave_speed + u_0)     ! equilibrium time step (0.85 factor corrects for minimum dx)
    else
       dt_init = 300d0 * SECOND * dx_scaling                        ! CAM value
       dt_max  = dt_init
    end if

    ! Set viscosities
    nu_dim         = nu_scale (Area_min, dt_max)                    ! viscosity scale on finest grid
    nu_scaling     = (res_scaling * dx_scaling)**(2*Laplace_order)
    nu             = nu_CAM * nu_scaling                            ! scaled CAM viscosity

    ! Limit viscosity to stable values
    nu_sclr        = min (nu,          nu_dim * C_max                     )  
    nu_rotu        = min (nu,          nu_dim * C_max / 4d0**Laplace_order)
    nu_divu        = min (nu * 2.5d0,  nu_dim * C_max                     )
    
    ! Equivalent non-dimensional viscosities
    ! (CAM non-dimensional viscosity is C =  (1e15 m^4/s) (300 s) / (120 km)^4 = 1.45e-3)
    C_visc(S_MASS) = nu_sclr / nu_dim 
    C_visc(S_TEMP) = nu_sclr / nu_dim 
    C_visc(S_VELO) = nu_rotu / nu_dim
    C_div          = nu_divu / nu_dim
  end subroutine initialize_dt_viscosity_case

  real(8) function nu_scale (Area, dt)
    ! Viscosity scale
    ! (factor 1.5 ensures stability limit matches theoretical value)
    implicit none
    real(8) :: Area, dt

    nu_scale = 1.5d0 * Area**Laplace_order / dt
  end function nu_scale

  subroutine apply_initial_conditions_case
    implicit none
    
    if (NCAR_topo) then
       call apply_bdry (assign_NCAR_topo, z_null, 0, 1)
    else
       call apply_bdry (init_topo, z_null, 0, 1)
    end if
    topography%bdry_uptodate = .false.
    
    call apply_bdry (init_mean, z_null, 0, 1)
    call apply_bdry (init_sol,  z_null, 0, 1)
    sol%bdry_uptodate        = .false.
    sol_mean%bdry_uptodate   = .false.

    call update_bdry (topography, NONE)
    call update_bdry (sol,        NONE) 
    call update_bdry (sol_mean,   NONE)

    ! SSO parameters (requires topography)
    if (NCAR_TOPO .and. sso) then
       call apply (cal_sso_param, z_null)
       
       sso_param%bdry_uptodate = .false.
       call update_bdry (sso_param, NONE)
    end if
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update sol_mean and topography on new grid
    use wavelet_mod
    implicit none
    integer :: d, p

    if (istep /= 0) then
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             if (NCAR_topo) then
                call apply_onescale_to_patch (assign_NCAR_topo, grid(d), p-1, z_null, 0, 1)
             else
                call apply_onescale_to_patch (init_topo, grid(d), p-1, z_null, 0, 1)
             end if
             call apply_onescale_to_patch (init_mean, grid(d), p-1, z_null, 0, 1)
          end do
       end do
    else ! need to set values over entire grid on restart
       if (NCAR_topo) then
          call apply_bdry (assign_NCAR_topo, z_null, 0, 1)
       else
          call apply_bdry (init_topo, z_null, 0, 1)
       end if
       call apply_bdry (init_mean, z_null, 0, 1)
    end if
    topography%bdry_uptodate = .false.
    sol_mean%bdry_uptodate   = .false.
    
    call update_bdry (topography, NONE)
    call update_bdry (sol_mean,   NONE)
    
    ! SSO parameters (requires topography)
    if (NCAR_topo .and. sso) then
       if (istep /= 0) then
          do d = 1, size(grid)
             do p = n_patch_old(d)+1, grid(d)%patch%length
                if (sso) call apply_onescale_to_patch (cal_sso_param, grid(d), p-1, z_null, 0, 1)
             end do
          end do
       else ! need to set values over entire grid on restart
          call apply_bdry (cal_sso_param, z_null, 0, 1)
       end if
       sso_param%bdry_uptodate = .false.
       
       call update_bdry (sso_param, NONE)
    end if
  end subroutine update_case

  subroutine set_save_level_case
    implicit none
  end subroutine set_save_level_case

  subroutine initialize_seed
    implicit none

    integer                            :: k
    integer, dimension(1:8)            :: values
    integer, dimension(:), allocatable :: seed

    call date_and_time (values=values)
    call random_seed(size=k)
    allocate (seed(1:k))
    seed = values(8)
    call random_seed (put=seed)
  end subroutine initialize_seed

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
    ! Dummy routine
    ! (see upwelling test case for example)
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    z_coords_case = 0d0
  end function z_coords_case

  real(8) function cfl (t)
    ! Gradually increase cfl number from cfl_min to cfl_max over T_cfl 
    implicit none
    real(8) :: t

    if (t - time_start <= T_cfl) then
       cfl = cfl_min + (cfl_max - cfl_min) * sin (MATH_PI/2d0 * (t - time_start) / T_cfl)
    else
       cfl = cfl_max
    end if
  end function cfl
end module test_case_mod
