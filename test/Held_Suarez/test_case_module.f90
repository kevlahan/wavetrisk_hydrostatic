module test_case_mod
  ! Module file for Held & Suarez (1994) test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none

  ! Standard variables
  integer :: CP_EVERY, resume_init
  real(8) :: dt_cfl, total_cpu_time, dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim

  ! Test case variables
  real(8) :: delta_T, delta_theta, sigma_b, sigma_c, k_a, k_f, k_s, T_0, T_mean, T_tropo
  real(8) :: delta_T2, sigma_t, sigma_v, sigma_0, gamma_T, u_0
  logical :: NCAR_topo
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
    real(8)                    :: visc
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

       if (Laplace_order == 0) then
          visc = 0d0
       else ! scale aware viscosity
          visc =  C_visc(v) * dom%len%elts(EDGE*id+RT+1)**(2d0*Laplace_order_init)/dt
       end if
       physics_scalar_flux_case = visc * (-1)**Laplace_order * grad * l_e
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

    integer                    :: id
    real(8)                    :: visc
    
    id = idx (i, j, offs, dims)

    ! Scale aware viscosity
    if (Laplace_order == 0) then
       visc = 0d0
    else
       visc = C_visc(S_VELO) * dom%len%elts(EDGE*id+RT+1)**(2d0*Laplace_order)/dt
    end if
    physics_velo_source_case =  visc * (-1)**(Laplace_order-1) * (grad_divu() - curl_rotu())
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
  end function physics_velo_source_case

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! From Jablonowski and Williamson (2006) without perturbation
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: p, p_s, pot_temp, sigma

    d = dom%id+1

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    ! Surface pressure
    dom%surf_press%elts(id+1) = surf_pressure (x_i)
    p_s = dom%surf_press%elts(id+1)

    ! Pressure at level zlev
    p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

    ! Normalized pressure
    sigma = (p - p_top) / (p_s - p_top)
    sigma_v = (sigma - sigma_0) * MATH_PI/2

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*p_s/grav_accel

    ! Potential temperature
    pot_temp = set_temp (x_i, sigma) * (p/p_0)**(-kappa)
    !     call cal_theta_eq (p, p_s, lat, theta_equil, k_T)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)

    ! Means are zero
    sol_mean(S_MASS,zlev)%data(d)%elts(id+1)                      = 0d0
    sol_mean(S_TEMP,zlev)%data(d)%elts(id+1)                      = 0d0
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
  end subroutine init_sol

  subroutine cal_theta_eq (p, p_s, lat, theta_equil, k_T)
    ! Returns equilibrium potential temperature theta_equil and Newton cooling constant k_T
    use domain_mod
    implicit none
    real(8) :: p, p_s, lat, theta_equil, k_T

    real(8) :: cs2, sigma, theta_force, theta_tropo

    cs2 = cos (lat)**2

    sigma = (p - p_top) / (p_s - p_top)

    k_T = k_a + (k_s-k_a) * max (0d0, (sigma-sigma_b)/sigma_c) * cs2**2

    theta_tropo = T_tropo * (p/p_0)**(-kappa) ! Potential temperature at tropopause

    theta_force = T_mean - delta_T*(1d0-cs2) - delta_theta*cs2 * log (p/p_0)

    theta_equil = max (theta_tropo, theta_force) ! Equilibrium temperature
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
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel)
    else
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel) + delta_T * (sigma_t - sigma)**5
    end if

    set_temp = Tmean + 0.75_8 * sigma*MATH_PI*u_0/R_d * sin(sigma_v) * sqrt(cos(sigma_v)) * &
         (2*u_0*cos(sigma_v)**1.5*(-2*sn2**3*(cs2+1/3d0) + 10/63d0) &
         + radius*omega*(8/5d0*cs2**1.5*(sn2+2/3d0) - MATH_PI/4))
  end function set_temp

  real(8) function surf_geopot_case (dom, id)
    ! Surface geopotential from Jablonowski and Williamson (2006)
    implicit none
    integer       :: id
    type (Domain) :: dom
    
    Type(Coord) :: x_i
    real(8)     :: c1, cs2, lon, lat, sn2

    if (NCAR_topo) then 
       surf_geopot_case = topography%data(d)%elts(id)   
    else ! surface geopotential from Jablonowski and Williamson (2006)
       x_i = dom%node%elts(id)
       call cart2sph (x_i, lon, lat)
       cs2 = cos (lat)**2; sn2 = sin (lat)**2

       c1 = u_0 * cos((1d0 - sigma_0) * MATH_PI/2d0)**1.5

       surf_geopot_case =  c1 * (c1 * (-2d0 * sn2**3 * (cs2 + 1d0/3d0) + 10d0/63d0) &
            + radius * omega * (8d0/5d0 * cs2**1.5 * (sn2 + 2d0/3d0) - MATH_PI/4d0))
    end if
  end function surf_geopot_case

  real(8) function surf_pressure (x_i)
    ! Surface pressure
    implicit none
    type(Coord) :: x_i

    surf_pressure = p_0
  end function surf_pressure

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    real(8) :: r

    call random_number (r)

    u = u_0 * cos (sigma_v)**1.5 * sin (2*lat)**2 + r ! Zonal velocity component
    !    u = 0d0 ! Uniform
    v = 0d0                                 ! Meridional velocity component
  end subroutine vel_fun

  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                     :: k
    real(8), dimension(1:N_VARIABLE,1:zlevels) :: threshold_new
    character(3), parameter                     :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       call cal_lnorm_sol (sol, order)
       !threshold_new = max (tol*lnorm, threshold_def) ! Avoid very small thresholds before instability develops
       threshold_new = tol*lnorm
    end if

    if (istep >= 10) then
       threshold = 0.01*threshold_new + 0.99*threshold
    else
       threshold = threshold_new
    end if
  end subroutine set_thresholds_case

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
          a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1d0 /)
       elseif (zlevels==26) then
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
               0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
               0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1d0 /)
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
               0d0, 0d0 /)
          b_vert=(/0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, &
               0.006755112_8, 0.01400364_8, 0.02178164_8, 0.03012778_8, 0.03908356_8, 0.04869352_8, &
               0.05900542_8, 0.07007056_8, 0.08194394_8, 0.09468459_8, 0.1083559_8, 0.1230258_8, &
               0.1387673_8, 0.1556586_8, 0.1737837_8, 0.1932327_8, 0.2141024_8, 0.2364965_8, &
               0.2605264_8, 0.2863115_8, 0.3139801_8, 0.3436697_8, 0.3755280_8, 0.4097133_8, &
               0.4463958_8, 0.4857576_8, 0.5279946_8, 0.5733168_8, 0.6219495_8, 0.6741346_8, &
               0.7301315_8, 0.7897776_8, 0.8443334_8, 0.8923650_8, 0.9325572_8, 0.9637744_8, &
               0.9851122_8, 1d0/)
       end if
       a_vert = a_vert(zlevels+1:1:-1) * p_0
       b_vert = b_vert(zlevels+1:1:-1)
    end if

    ! LMDZ grid
    !call cal_AB

    p_top = a_vert(zlevels+1)

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
       dsig(l) = 1d0 + 7 * sin (MATH_PI*(l-0.5_8)/(zlevels+1))**2 ! LMDZ standard (concentrated near top and surface)
       !dsig(l) = 1d0 + 7 * cos (MATH_PI/2*(l-0.5_8)/(zlevels+1))**2 ! Concentrated at top
       !dsig(l) = 1d0 + 7 * sin (MATH_PI/2*(l-0.5_8)/(zlevels+1))**2 ! Concentrated at surface
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
    if (rank == 0) write (6,'(a,A)') "Input file = ", trim (filename)

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, NCAR_topo
    read (fid,*) varname, tol
    read (fid,*) varname, press_save
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)

    allocate (pressure_save(1))
    pressure_save(1) = 1.0d2*press_save
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
          command = 'tar xzf '//trim(run_id)//'.3.tgz'
          write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
          call system (bash_cmd)

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
       write (6,'(a,l1)')     "uniform             = ", uniform
       write (6,'(a,l1)')     "remap               = ", remap
       write (6,'(a,a)')      "remapscalar_type    = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type      = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap              = ", iremap
       write (6,'(a,l1)')     "default_thresholds  = ", default_thresholds
       write (6,'(a,es10.4)') "tolerance           = ", tol
       write (6,'(a,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(a,l1)')     "adapt_dt            = ", adapt_dt
       write (6,'(a,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(a,a)')      "timeint_type        = ", trim (timeint_type)       
       write (6,'(a,es10.4)') "pressure_save (hPa) = ", pressure_save(1)/100
       write (6,'(a,i1)')     "Laplace_order       = ", Laplace_order_init
       write (6,'(a,i2)')     "n_diffuse           = ", n_diffuse
       write (6,'(a,es10.4)') "dt_write (day)      = ", dt_write/DAY
       write (6,'(a,i6)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance           = ", rebalance
       write (6,'(a,es10.4)') "time_end (day)      = ", time_end/DAY
       write (6,'(a,i6)')     "resume              = ", resume_init

       write (6,'(/,a)')      "STANDARD PARAMETERS"
       write (6,'(a,es10.4)') "radius              = ", radius
       write (6,'(a,es10.4)') "omega               = ", omega
       write (6,'(a,es10.4)') "p_0   (hPa)         = ", p_0/100
       write (6,'(a,es10.4)') "p_top (hPa)         = ", p_top/100
       write (6,'(a,es10.4)') "R_d                 = ", R_d
       write (6,'(a,es10.4)') "c_p                 = ", c_p
       write (6,'(a,es10.4)') "c_v                 = ", c_v
       write (6,'(a,es10.4)') "gamma               = ", gamma
       write (6,'(a,es10.4)') "kappa               = ", kappa

       write (6,'(/,a)')      "TEST CASE PARAMETERS"
       write (6,'(a,es10.4)') "T_0                 = ", T_0
       write (6,'(a,es10.4)') "T_mean              = ", T_mean
       write (6,'(a,es10.4)') "T_tropo             = ", T_tropo
       write (6,'(a,es10.4)') "sigma_b             = ", sigma_b
       write (6,'(a,es10.4)') "k_a                 = ", k_a
       write (6,'(a,es10.4)') "k_f                 = ", k_f
       write (6,'(a,es10.4)') "k_s                 = ", k_s
       write (6,'(a,es10.4)') "delta_T             = ", delta_T
       write (6,'(a,es10.4)') "delta_theta         = ", delta_theta
       write (6,'(a)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

    if (rank == 0) then
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i12,4(a,es8.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', sum (threshold(S_MASS,:))/zlevels, &
            ' temp tol = ', sum (threshold(S_TEMP,:))/zlevels, &
            ' velo tol = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(5(es15.9,1x),i2,1x,i12,1x,4(es15.9,1x))')  &
            time/DAY, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    integer :: k

    allocate (threshold(1:N_VARIABLE,1:zlevels));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,1:zlevels)); threshold_def = 0d0

    lnorm(S_MASS,:) = dPdim/grav_accel
    do k = 1, zlevels
       lnorm(S_TEMP,k) = (a_vert_mass(k) + b_vert_mass(k)*Pdim/grav_accel)*dTempdim
    end do
    lnorm(S_TEMP,:) = lnorm(S_TEMP,:) + Tempdim*lnorm(S_MASS,:) ! Add component due to tendency in mass 
    lnorm(S_VELO,:) = Udim
    threshold_def = tol * lnorm
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case
  end subroutine initialize_dt_viscosity_case

  subroutine apply_initial_conditions_case
    implicit none
    integer :: k, l

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update means, bathymetry and penalization mask
    ! not needed in this test case
    use wavelet_mod
    implicit none
    integer :: d, k, l, p

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

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, x_i,  x_E)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel(vel_fun, x_NE, x_i)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, x_i,  x_N)
  end subroutine vel2uvw

  subroutine set_save_level_case
    ! Determines closest vertical level to desired pressure
    implicit none
    integer :: k
    real(8) :: dpress, p, save_press

    dpress = 1d16; save_zlev = 0
    do k = 1, zlevels
       p = 0.5 * (a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1))*p_0)
       if (abs(p-pressure_save(1)) < dpress) then
          dpress = abs (p-pressure_save(1))
          save_zlev = k
          save_press = p
       end if
    end do
    if (rank==0) write (6,'(/,A,i2,A,f5.1,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate pressure = ", save_press/100, " hPa)"
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
          mean_m =>  sol_mean(S_MASS,k)%data(d)%elts
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

    integer :: id_i
    real(8) :: k_T, lat, lon, theta_equil

    id_i = idx (i, j, offs, dims) + 1
    call cart2sph (dom%node%elts(id_i), lon, lat)
    call cal_theta_eq (dom%press%elts(id_i), dom%surf_press%elts(id_i), lat, theta_equil, k_T)

    dmass(id_i) = 0d0
    dtemp(id_i) = - k_T * (temp(id_i) - theta_equil*mass(id_i))
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
    id_i = id+1

    sigma = (dom%press%elts(id_i) - p_top) / (dom%surf_press%elts(id_i) - p_top)
    k_v = k_f * max (0d0, (sigma-sigma_b)/sigma_c)

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

    dom%press%elts(id_i) = 0.5 * (dom%press_lower%elts(id_i) + p_upper)
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
          mass   => q(S_MASS,k)%data(d)%elts
          temp   => q(S_TEMP,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass_HS, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, mean_m, mean_t, temp)
       end do
       grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + p_top

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

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    dom%surf_press%elts(id_i) = dom%surf_press%elts(id_i) + mass(id_i)
  end subroutine column_mass_HS

  subroutine set_surf_geopot_HS (dom, i, j, zlev, offs, dims)
    ! Set initial geopotential to surface geopotential
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    dom%geopot%elts(id_i) = surf_geopot_case (dom%node%elts(id_i))
  end subroutine set_surf_geopot_HS

  function z_coords_case (eta_surf, z_s)
    ! Dummy routine
    ! (see upwelling test case for example)
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    z_coords_case = 0d0
  end function z_coords_case
end module test_case_mod
