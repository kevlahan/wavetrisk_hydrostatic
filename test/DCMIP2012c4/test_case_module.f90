module test_case_mod
  ! Module file for DCMIP2012c4
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none

  ! Standard variables
  integer :: CP_EVERY, resume_init
  real(8) :: dt_cfl, total_cpu_time, dPdim, R_ddim, specvoldim, dTempdim

  ! Test case variables
  real(8) :: delta_T, sigma, sigma_t, sigma_v, sigma_0, gamma_T, lat_c, lon_c, R_pert, T_0, u_p, u_0
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


    physics_scalar_flux_case = 0.0_8
    if (Laplace_sclr /= 0) then
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
       if (Laplace_sclr == 1) then
          grad = grad_physics (q(v,zlev)%data(d)%elts)
       elseif (Laplace_sclr == 2) then
          grad = grad_physics (Laplacian_scalar(v)%data(d)%elts)
       end if

       ! Complete scalar diffusion
       physics_scalar_flux_case = (-1)**Laplace_sclr * visc_sclr(v) * grad * l_e
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
    real(8), dimension(1:EDGE) :: diffusion

    id = idx (i, j, offs, dims)

    if (Laplace_rotu == 0) then
       diffusion = 0.0_8
    else
       ! Calculate Laplacian of velocity
       diffusion =  (-1)**(Laplace_rotu-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
    end if

    ! Total physics for source term of velocity trend
    physics_velo_source_case =  diffusion
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
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: p, p_s, pot_temp

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
    dom%surf_press%elts(id+1) = surf_pressure (d, id+1)
    p_s = dom%surf_press%elts(id+1)

    ! Pressure at level zlev
    p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

    ! Normalized pressure
    sigma = p / p_s
    sigma_v = (sigma - sigma_0) * MATH_PI/2

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*p_s/grav_accel

    ! Potential temperature
    pot_temp =  set_temp(x_i) * (p/p_0)**(-kappa)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)

    ! Means are zero
    sol_mean(S_MASS,zlev)%data(d)%elts(id+1)                      = 0d0
    sol_mean(S_TEMP,zlev)%data(d)%elts(id+1)                      = 0d0
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
  end subroutine init_sol

  real(8) function set_temp (x_i)
    implicit none
    type(Coord) :: x_i

    real(8) :: cs2, lon, lat, sn2, Tmean

    call cart2sph (x_i, lon, lat)
    sn2 = sin(lat)**2
    cs2 = cos(lat)**2

    if (sigma >= sigma_t) then
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel)
    else
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel) + delta_T * (sigma_t - sigma)**5
    end if

    set_temp = Tmean + 0.75_8 * sigma*MATH_PI*u_0/R_d * sin(sigma_v) * sqrt(cos(sigma_v)) * &
         (2*u_0*cos(sigma_v)**1.5*(-2*sn2**3*(cs2+1/3.0_8) + 10/63.0_8) &
         + radius*omega*(8/5.0_8*cs2**1.5*(sn2+2/3.0_8) - MATH_PI/4))
  end function set_temp

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential
    implicit none
    integer :: d, id
    
    type(Coord) :: x_i
    real(8)     :: c1, cs2, lon, lat, sn2

    x_i = grid(d)%node%elts(id)

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2

    c1 = u_0*cos((1.0_8-sigma_0)*MATH_PI/2)**1.5
    surf_geopot_case =  c1*(c1*(-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8) &
         + radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
  end function surf_geopot_case
  
  real(8) function surf_pressure (d, id)
    ! Surface pressure
    implicit none
    integer :: d, id

    surf_pressure = p_0
  end function surf_pressure

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    real(8) :: rgrc

    ! Great circle distance
    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))

    u = u_0*cos(sigma_v)**1.5*sin(2*lat)**2 + u_p*exp__flush(-(rgrc/R_pert)**2)  ! Zonal velocity component
    v = 0.0_8         ! Meridional velocity component
  end subroutine vel_fun

  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    real(8), dimension(1:N_VARIABLE,1:zlevels) :: threshold_new
    character(3), parameter                    :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       call cal_lnorm (order)
       threshold_new = max (tol*lnorm, threshold_def) ! Avoid very small thresholds before instability develops
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
    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = dble(k-1)/dble(zlevels) * p_top
          b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
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
          write (0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
          stop
       end if

       ! DCMIP order is opposite to ours
       a_vert = a_vert(zlevels+1:1:-1) * p_0
       b_vert = b_vert(zlevels+1:1:-1)
    end if

    ! LMDZ grid
    !call cal_AB

    ! Set pressure at infinity
    p_top = a_vert(zlevels+1) ! note that b_vert at top level is 0, a_vert is small but non-zero

    ! Set mass coefficients 
    a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1)) / grav_accel
    b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
  end subroutine initialize_a_b_vert_case

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
       !dsig(l) = 1.0_8 + 7 * cos (MATH_PI/2*(l-0.5_8)/(zlevels+1))**2 ! Concentrated at top
       !dsig(l) = 1.0_8 + 7 * sin (MATH_PI/2*(l-0.5_8)/(zlevels+1))**2 ! Concentrated at surface
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
    integer, parameter :: fid = 500
    character(255)     :: filename, varname

    ! Find input parameters file name
    if (command_argument_count() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'test_case.in'
    end if
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

    open(unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, compressible
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, uniform
    read (fid,*) varname, remap
    read (fid,*) varname, remapscalar_type
    read (fid,*) varname, remapvelo_type
    read (fid,*) varname, iremap
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, tol
    read (fid,*) varname, optimize_grid
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, cfl_num
    read (fid,*) varname, timeint_type
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, rebalance
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)

    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case           = ", trim (test_case)
       write (6,'(A,A)')      "run_id              = ", trim (run_id)
       write (6,'(A,L1)')     "compressible        = ", compressible
       write (6,'(A,i3)')     "min_level           = ", min_level
       write (6,'(A,i3)')     "max_level           = ", max_level
       write (6,'(A,i5)')     "number of domains   = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL        = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL         = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels             = ", zlevels
       write (6,'(A,L1)')     "uniform             = ", uniform
       write (6,'(A,L1)')     "remap               = ", remap
       write (6,'(a,a)')      "remapscalar_type    = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type      = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap              = ", iremap
       write (6,'(A,L1)')     "default_thresholds  = ", default_thresholds
       write (6,'(A,es10.4)') "tolerance           = ", tol
       write (6,'(A,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt            = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(a,a)')      "timeint_type        = ", trim (timeint_type)
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", pressure_save(1)/100
       write (6,'(A,es10.4)') "dt_write (days)     = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance           = ", rebalance
       write (6,'(A,es10.4)') "time_end (days)     = ", time_end/DAY
       write (6,'(A,i6)')     "resume              = ", resume_init

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius              = ", radius
       write (6,'(A,es10.4)') "omega               = ", omega
       write (6,'(A,es10.4)') "p_0   (hPa)         = ", p_0/100
       write (6,'(A,es10.4)') "p_top (hPa)         = ", p_top/100
       write (6,'(A,es10.4)') "R_d                 = ", R_d
       write (6,'(A,es10.4)') "c_p                 = ", c_p
       write (6,'(A,es10.4)') "c_v                 = ", c_v
       write (6,'(A,es10.4)') "gamma               = ", gamma
       write (6,'(A,es10.4)') "kappa               = ", kappa

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es10.4)') "u_0                  = ", u_0
       write (6,'(A,es10.4)') "u_p                  = ", u_p
       write (6,'(A,es10.4)') "R_pert               = ", R_pert
       write (6,'(A,es10.4)') "T_0                  = ", T_0
       write (6,'(A,es10.4)') "T_0                  = ", R_pert
       write (6,'(A,es10.4)') "gamma_T              = ", gamma_T
       write (6,'(A,es10.4)') "delta_T              = ", delta_T
       write (6,'(A,es10.4)') "sigma_0              = ", sigma_0
       write (6,'(A,es10.4)') "sigma_t              = ", sigma_t
       write (6,'(A,es10.4)') "lon_c                = ", lon_c
       write (6,'(A,es10.4)') "lat_c                = ", lat_c

       write (6,'(A)') &
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
            time/HOUR, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer :: k

    lnorm(S_MASS,:) = dPdim/grav_accel
    do k = 1, zlevels
       lnorm(S_TEMP,k) = (a_vert_mass(k) + b_vert_mass(k)*Pdim/grav_accel)*dTempdim
    end do
    lnorm(S_TEMP,:) = lnorm(S_TEMP,:) + Tempdim*lnorm(S_MASS,:) ! Add component due to tendency in mass 
    lnorm(S_VELO,:) = Udim
    threshold_def = tol * lnorm
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, tau_divu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

    ! Diffusion constants
    C_sclr = 2d-3       ! <= 1e-3 for hyperdiffusion 
    C_divu = C_sclr   ! from eigenvalues of second order Laplacian
    C_rotu = C_sclr / 4**Laplace_rotu ! <= 1.09e-3 for hyperdiffusion (lower than exact limit 1/24^2 = 1.7e-3 due to non-uniform grid)

    ! CFL limit for time step
    dt_cfl = cfl_num*dx_min/(wave_speed+Udim) * 0.85 ! corrected for dynamic value
    dt_init = dt_cfl

    tau_sclr = dt_cfl / C_sclr
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_rotu == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_rotu == 1 .or. Laplace_rotu == 2) then
       visc_sclr = dx_min**(2*Laplace_rotu) / tau_sclr
       visc_divu = dx_min**(2*Laplace_rotu) / tau_divu
       visc_rotu = dx_min**(2*Laplace_rotu) / tau_rotu
    elseif (Laplace_rotu > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,3(a,es8.2),a,/)') "dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2))') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity_case

  subroutine apply_initial_conditions_case
    implicit none
    integer :: k, l

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -1, 1)
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

    z_coords_case = 0.0_8
  end function z_coords_case
end module test_case_mod
