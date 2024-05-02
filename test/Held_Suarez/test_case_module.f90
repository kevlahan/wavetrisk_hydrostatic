module test_case_mod
  ! Module file for Held & Suarez (1994) test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use std_atm_profile_mod
  use io_mod
  implicit none

  ! Standard variables
  integer :: CP_EVERY, resume_init
  real(8) :: time_start, total_cpu_time

  ! Test case variables
  real(8) :: Area_max, Area_min, C_div, delta_T, delta_theta, dt_max, k_a, k_f, k_s, specvoldim, T_0, T_mean, T_tropo
  real(8) :: delta_T2, sigma_t, sigma_v, sigma_0, gamma_T, sigma_b, sigma_c, u_0
  real(8) :: cfl_max, cfl_min, T_cfl, nu_sclr, nu_rotu, nu_divu
  logical :: scale_aware = .true.
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
      ! factor 3 ensures that maximum stable C_visc matches theoretical estimate of 1/6^Laplace_order
      implicit none
      integer :: id

      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         visc = C_visc(v) * 3d0 * dom%areas%elts(id+1)%hex_inv**(-Laplace_order) / dt
      else
         visc = C_visc(v) * 3d0 * Area_min**Laplace_order / dt
      end if
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
      real(8), dimension(3) :: grad_divu

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
      real(8), dimension(3) :: curl_rotu

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
      ! factor 3 ensures that maximum stable C_visc matches theoretical estimate of 1/24^Laplace_order
      implicit none
      integer :: id

      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         visc_div = C_div * 3d0 * dom%areas%elts(id+1)%hex_inv**(-Laplace_order) / dt
      else
         visc_div = C_div * 3d0 * Area_min**Laplace_order / dt
      end if
    end function visc_div

    real(8) function visc_rot (id)
      ! Scale aware viscosity
      ! factor 3 ensures that maximum stable C_visc matches theoretical estimate of 1/24^Laplace_order
      implicit none
      integer :: id

      if (scale_aware .and. dom%areas%elts(id+1)%hex_inv /= 0d0) then
         visc_rot = C_visc(S_VELO) * 3d0 * dom%areas%elts(id+1)%hex_inv**(-Laplace_order) / dt
      else
         visc_rot = C_visc(S_VELO) * 3d0 * Area_min**Laplace_order / dt
      end if
    end function visc_rot
  end function physics_velo_source_case
  
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: id, d, k
    
    d  = dom%id+1
    id = idx (i, j, offs, dims)

    do k = 1, zlevels
       sol(S_MASS,k)%data(d)%elts(id+1) = 0d0 
       sol(S_TEMP,k)%data(d)%elts(id+1) = 0d0
       if (NCAR_topo) then
          sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
       else
          call vel2uvw (dom, i, j, k, offs, dims, vel_fun)
       end if
    end do
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: id, d, k, l
    real(8)     :: k_T, lat, lon, p, p_s, pot_temp
    type(Coord) :: x_i
    
    d   = dom%id+1
    id  = idx (i, j, offs, dims)
    l   = dom%level%elts(id+1)
    x_i = dom%node%elts(id+1)
    
    call cart2sph (x_i, lon, lat)

    ! Surface pressure
    dom%surf_press%elts(id+1) = surf_pressure (d, id+1)
    p_s = dom%surf_press%elts(id+1)

    do k = 1, zlevels
       p = 0.5d0 * (a_vert(k) + a_vert(k+1) + (b_vert(k) + b_vert(k+1)) * p_s) ! pressure at level k

       ! Potential temperature
       call cal_theta_eq (p, p_s, lat, pot_temp, k_T)

       sol_mean(S_MASS,k)%data(d)%elts(id+1) = a_vert_mass(k) + b_vert_mass(k) * p_s / grav_accel
       sol_mean(S_TEMP,k)%data(d)%elts(id+1) = sol_mean(S_MASS,k)%data(d)%elts(id+1) * pot_temp
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
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
    !theta_equil = Tempdim
  end subroutine cal_theta_eq

  real(8) function set_temp (x_i, sigma)
    ! From Jablonowski and Williamson (2006)
    implicit none
    type(Coord) :: x_i
    real(8)     :: sigma

    real(8) :: cs2, lon, lat, sn2, Tmean

    call cart2sph (x_i, lon, lat)
    sn2 = sin (lat)**2
    cs2 = cos (lat)**2

    if (sigma >= sigma_t) then
       Tmean = T_0 * sigma**(R_d * Gamma_T / grav_accel)
    else
       Tmean = T_0 * sigma**(R_d * Gamma_T / grav_accel) + delta_T * (sigma_t - sigma)**5
    end if

    set_temp = Tmean + 0.75d0 * sigma * MATH_PI * u_0 / R_d * sin(sigma_v) * sqrt(cos(sigma_v)) * &
         (2d0*u_0*cos(sigma_v)**1.5 * (-2d0*sn2**3 * (cs2 + 1d0/3d0) + 10d0/63d0) &
         + radius * omega * (8d0/5d0*cs2**1.5 * (sn2+2/3d0) - MATH_PI/4d0))
  end function set_temp

  real(8) function surf_geopot_case (d, id)
    implicit none
    integer :: d, id
    
    real(8) :: c1, cs2, lon, lat, sn2

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
  end function surf_geopot_case

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

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    real(8) :: rgrc
    real(8) :: lat_c, lon_c, r
    real(8) :: amp = 1d0 ! amplitude of random noise

    ! Zonal velocity component
    call random_number (r)
    u = u_0 * cos (sigma_v)**1.5 * sin (2d0*lat)**2  + amp * 2d0 * (r - 0.5d0)

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

  subroutine cal_AB
    ! Computes A and B coefficients for hybrid vertical grid as in LMDZ
    implicit none

    integer                         :: l
    real(8)                         :: snorm
    real(8), dimension(1:zlevels)   :: dsig
    real(8), dimension(1:zlevels+1) :: sig

    snorm  = 0d0
    do l = 1, zlevels
       dsig(l) = 1d0 + 7 * sin (MATH_PI*(l-0.5d0)/(zlevels+1))**2 ! LMDZ standard (concentrated near top and surface)
       !dsig(l) = 1d0 + 7 * cos (MATH_PI/2*(l-0.5d0)/(zlevels+1))**2 ! Concentrated at top
       !dsig(l) = 1d0 + 7 * sin (MATH_PI/2*(l-0.5d0)/(zlevels+1))**2 ! Concentrated at surface
       snorm = snorm + dsig(l)
    end do

    do l = 1, zlevels
       dsig(l) = dsig(l)/snorm
    end do

    sig(zlevels+1) = 0d0
    do l = zlevels, 1, -1
       sig(l) = sig(l+1) + dsig(l)
    end do

    b_vert(zlevels+1) = 0d0
    do  l = 1, zlevels
       b_vert(l) = exp (1d0 - 1/sig(l)**2)
       a_vert(l) = (sig(l) - b_vert(l)) * p_0
    end do
    b_vert(1) = 1d0
    a_vert(1) = 0d0
    a_vert(zlevels+1) = (sig(zlevels+1) - b_vert(zlevels+1)) * p_0
  end subroutine cal_AB

  subroutine read_test_case_parameters
    implicit none
    integer            :: k, v
    integer, parameter :: fid = 500, funit = 400
    real(8)            :: press_save
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
    read (fid,*) varname, run_id
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, save_zlev
    read (fid,*) varname, NCAR_topo
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
    Laplace_order = Laplace_order_init

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
    integer :: k
    real(8) :: p_save

    if (rank==0) then
       k = save_zlev
       p_save = 0.5d0 * (a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1)) * p_0)
       
       write (6,'(a)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(a)')        "RUN PARAMETERS"
       write (6,'(a,a)')      "test_case           = ", trim (test_case)
       write (6,'(a,a)')      "run_id              = ", trim (run_id)
       write (6,'(a,l1)')     "compressible        = ", compressible
       write (6,'(a,i3)')     "min_level           = ", min_level
       write (6,'(a,i3)')     "max_level           = ", max_level
       write (6,'(a,i5)')     "number of domains   = ", N_GLO_DOMAIN
       write (6,'(a,i5)')     "number of processors = ", n_process
       write (6,'(a,i5)')     "DOMAIN_LEVEL        = ", DOMAIN_LEVEL
       write (6,'(a,i5)')     "PATCH_LEVEL         = ", PATCH_LEVEL
       write (6,'(a,i3)')     "zlevels             = ", zlevels
       write (6,'(a,i3,1x,f6.1,a)')     "save_zlev           = ", save_zlev, p_save/1d2, ' [hPa]'
       write (6,'(a,l1)')     "uniform             = ", uniform
       write (6,'(a,l1)')     "remap               = ", remap
       write (6,'(a,i3)')     "iremap              = ", iremap
       write (6,'(a,l1)')     "default_thresholds  = ", default_thresholds
       write (6,'(a,es10.4)') "tolerance           = ", tol
       write (6,'(a,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(a,l1)')     "adapt_dt            = ", adapt_dt
       write (6,'(a,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(a,a)')      "timeint_type        = ", trim (timeint_type)
       
       write (6,'(a,i1,/)')     "Laplace_order       = ", Laplace_order_init
       write (6,'(a,/,a,/,/,a,es8.2,/,a,es8.2,/)') "Stability limits:", &
            "[Klemp 2017 Damping Characteristics of Horizontal Laplacian Diffusion Filters Mon Weather Rev 145, 4365-4379.]", &
            "C_visc(S_MASS) and C_visc(S_TEMP) <  (1/6)**Laplace_order = ", (1d0/6d0)**Laplace_order_init, &
            "                   C_visc(S_VELO) < (1/24)**Laplace_order = ", (1d0/24d0)**Laplace_order_init
       if (scale_aware) then
          write (6,'(a,/)') "Scale-aware horizontal viscosity"
       else
          write (6,'(a,/)') "Horizontal viscosity based on dx_min"
       end if

       write (6,'(a)') "Non-dimensional viscosities"
       write (6,'(3(a,es8.2/))') "C_visc(S_MASS)     = ", C_visc(S_MASS), "C_visc(S_TEMP)     = ", C_visc(S_TEMP), &
            "C_visc(S_VELO)     = ", C_visc(S_VELO)

       write (6,'(a)')        "Approximate viscosities on finest grid"
       write (6,'(a,es8.2)') "nu_scalar           = ", nu_sclr
       write (6,'(a,es8.2)') "nu_rot              = ", nu_rotu
       write (6,'(a,es8.2,/)') "nu_div              = ", nu_divu

       write (6,'(a,es10.4)') "dt_max          [s] = ", dt_max
       write (6,'(a,es10.4)') "dt_write        [d] = ", dt_write / DAY
       write (6,'(a,i3)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance           = ", rebalance
       write (6,'(a,es10.4)') "time_end        [d] = ", time_end / DAY
       write (6,'(a,i6)')     "resume              = ", resume_init

       write (6,'(/,a)')      "STANDARD PARAMETERS"
       write (6,'(a,es10.4)') "radius         [km] = ", radius / KM
       write (6,'(a,es10.4)') "omega       [rad/s] = ", omega
       write (6,'(a,es10.4)') "p_0           [hPa] = ", p_0/100d0
       write (6,'(a,es10.4)') "p_top         [hPa] = ", p_top/100d0
       write (6,'(a,es10.4)') "R_d      [J/(kg K)] = ", R_d
       write (6,'(a,es10.4)') "c_p      [J/(kg K)] = ", c_p
       write (6,'(a,es10.4)') "c_v      [J/(kg K)] = ", c_v
       write (6,'(a,es10.4)') "gamma               = ", gamma
       write (6,'(a,es10.4)') "kappa               = ", kappa
       write (6,'(a,f10.1)') "dx_max         [km] = ", dx_max / KM
       write (6,'(a,f10.1)') "dx_min         [km] = ", dx_min / KM

       write (6,'(/,a)')      "TEST CASE PARAMETERS"
       write (6,'(a,f5.1)')   "T_0             [K] = ", T_0
       write (6,'(a,es10.4)') "T_mean          [K] = ", T_mean
       write (6,'(a,es10.4)') "T_tropo         [K] = ", T_tropo
       write (6,'(a,es10.4)') "sigma_b             = ", sigma_b
       write (6,'(a,es10.4)') "k_a           [1/d] = ", k_a / DAY
       write (6,'(a,es10.4)') "k_f           [1/d] = ", k_f / DAY
       write (6,'(a,es10.4)') "k_s           [1/d] = ", k_s / DAY
       write (6,'(a,es10.4)') "delta_T       [K/m] = ", delta_T
       write (6,'(a,es10.4)') "delta_theta   [K/m] = ", delta_theta
       write (6,'(a,es10.4,/)') "wave_speed    [m/s] = ", wave_speed
       write (6,'(a,l)')      "NCAR_topo           = ", NCAR_topo
       write (6,'(a,a)')      "topo_file           = ", trim (topo_file)
       write (6,'(a,i3)')     "topo_min_level      = ", topo_min_level
       write (6,'(a,i3)')     "topo_max_level      = ", topo_max_level
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

    allocate (threshold(1:N_VARIABLE,zmin:zlevels));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,zmin:zlevels)); threshold_def = 0d0

    call cal_lnorm ("inf")
    lnorm(S_VELO,:) = Udim

    threshold_def = tol * lnorm

    threshold_def(S_VELO,:) = threshold_def(S_VELO,:) 
  end subroutine initialize_thresholds_case

  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    implicit none
    real(8), dimension(1:N_VARIABLE,zmin:zlevels) :: lnorm_mean
    
    call cal_lnorm ("inf")
    lnorm(S_VELO,:) = max (Udim,  (lnorm(S_VELO,:)))
    
    threshold = tol * lnorm

    threshold(S_VELO,:) = threshold(S_VELO,:) 
  end subroutine set_thresholds_case

  subroutine initialize_dt_viscosity_case
  end subroutine initialize_dt_viscosity_case

  subroutine apply_initial_conditions_case
    implicit none
    integer :: l
    
    do l = level_start, level_end
       if (NCAR_topo) call apply_onescale (assign_topo, l, z_null, 0, 1)
       call apply_onescale (init_mean, l, z_null, 0, 1)
       call apply_onescale (init_sol,  l, z_null, 0, 1)
    end do
    topography%bdry_uptodate = .false.
    sol%bdry_uptodate        = .false.
    sol_mean%bdry_uptodate   = .false.
    
    call update_bdry       (topography, NONE)
    call update_array_bdry (sol,        NONE)
    call update_array_bdry (sol_mean,   NONE)
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update means, bathymetry and penalization mask
    ! not needed in this test case
    use wavelet_mod
    implicit none
    integer :: d, l, p

    if (istep /= 0) then
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             if (NCAR_topo) call apply_onescale_to_patch (assign_topo, grid(d), p-1, z_null, 0, 1)
             call apply_onescale_to_patch (init_mean, grid(d), p-1, z_null, 0, 1)
          end do
       end do
    else ! need to set values over entire grid on restart
       do l = level_start, level_end
          if (NCAR_topo) call apply_onescale (assign_topo, l, z_null, 0, 1)
          call apply_onescale (init_mean, l, z_null, 0, 1)
       end do
    end if

    topography%bdry_uptodate = .false.
    sol%bdry_uptodate        = .false.
    sol_mean%bdry_uptodate   = .false.
    
    call update_bdry       (topography, NONE)
    call update_array_bdry (sol,        NONE)
    call update_array_bdry (sol_mean,   NONE)
  end subroutine update_case

  subroutine vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
    ! Sets the velocities on the computational grid given a function vel_fun that provides zonal and meridional velocities
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    external                        :: vel_fun

    integer      :: d, id, idE, idN, idNE
    type (Coord) :: x_i, x_E, x_N, x_NE

    d = dom%id+1

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel (vel_fun, x_i,  x_E)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel (vel_fun, x_NE, x_i)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel (vel_fun, x_i,  x_N)
  end subroutine vel2uvw

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

  subroutine trend_cooling (q, dq)
    ! Trend for Held-Suarez cooling
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (sol, NONE, 27)

    ! Current surface pressure
    call cal_surf_press_HS (sol)

    do k = 1, zlevels
       do d = 1, size(grid)
          mean_m =>   sol_mean(S_MASS,k)%data(d)%elts
          mass   =>  q(S_MASS,k)%data(d)%elts
          temp   =>  q(S_TEMP,k)%data(d)%elts
          velo   =>  q(S_VELO,k)%data(d)%elts
          dmass  => dq(S_MASS,k)%data(d)%elts
          dtemp  => dq(S_TEMP,k)%data(d)%elts
          dvelo  => dq(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (cal_press_HS,  grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_velo,    grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dtemp, dvelo, mass, temp, velo)
       end do
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_cooling

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Trend for cooling step (relaxation to equilibrium temperature)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: k_T, lat, lon, sigma, theta_equil

    id = idx (i, j, offs, dims) + 1

    dmass(id) = 0d0

    call cart2sph (dom%node%elts(id), lon, lat)
    
    call cal_theta_eq (dom%press%elts(id), dom%surf_press%elts(id), lat, theta_equil, k_T)

    if (NCAR_topo) then 
       sigma = (dom%press%elts(id) - p_top) / (dom%surf_press%elts(id) - p_top)

       if (sigma > 0.7d0) then ! no temperature relaxation in lower part of atmosphere
          dtemp(id) = 0d0
       else
          dtemp(id) = - k_T * temp(id)       
       end if
    else
       dtemp(id) = - k_T * temp(id)
    end if
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    ! Velocity trend for cooling step (Rayleigh friction)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i
    real(8) :: k_v, sigma

    id = idx (i, j, offs, dims)
    id_i = id + 1

    sigma = (dom%press%elts(id_i) - p_top) / (dom%surf_press%elts(id_i) - p_top)
    k_v = k_f * max (0d0, (sigma - sigma_b) / sigma_c)
    dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = - k_v * velo(EDGE*id+RT+1:EDGE*id+UP+1)
  end subroutine trend_velo

  subroutine cal_press_HS (dom, i, j, zlev, offs, dims)
    ! Integrate pressure up from surface to top layer
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: full_mass, p_upper

    id_i = idx (i, j, offs, dims) + 1

    full_mass = mass(id_i) + mean_m(id_i)
    p_upper = dom%press_lower%elts(id_i) - grav_accel * full_mass

    dom%press%elts(id_i) = 0.5d0 * (dom%press_lower%elts(id_i) + p_upper)
    dom%press_lower%elts(id_i) = p_upper
  end subroutine cal_press_HS

  subroutine cal_surf_press_HS (q)
    implicit none
    ! Compute surface pressure and save in press_lower for upward integration
    ! Set geopotential to surface geopotential for upward integration
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q

    integer :: d, k, mass_type, p

    call apply (set_surf_geopot_HS, z_null)

    do d = 1, size(grid)
       grid(d)%surf_press%elts = 0d0
       do k = 1, zlevels
          mass   =>        q(S_MASS,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass_HS, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, mean_m)
       end do
       grid(d)%surf_press%elts = grav_accel * grid(d)%surf_press%elts + p_top

       grid(d)%press_lower%elts = grid(d)%surf_press%elts
    end do
  end subroutine cal_surf_press_HS

  subroutine column_mass_HS (dom, i, j, zlev, offs, dims)
    ! Sum up mass
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: full_mass

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    full_mass = mean_m(id) + mass(id)

    dom%surf_press%elts(id) = dom%surf_press%elts(id) + full_mass
  end subroutine column_mass_HS

  subroutine set_surf_geopot_HS (dom, i, j, zlev, offs, dims)
    ! Set initial geopotential to surface geopotential
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    dom%geopot%elts(id) = surf_geopot_case (d, id)
  end subroutine set_surf_geopot_HS

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
