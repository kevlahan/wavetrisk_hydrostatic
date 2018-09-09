module test_case_mod
  ! Module file for DCMIP2012c4
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  ! Standard variables
  integer                              :: iwrite, CP_EVERY, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Test case variables
  real(8) :: delta_T, eta, eta_t, eta_v, eta_0, gamma_T, lat_c, lon_c,  R_pert, T_0, u_p, u_0
contains
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: column_mass, lev_press, pot_temp, p_top, p_bot

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
    column_mass = dom%surf_press%elts(id+1)/grav_accel

    ! Pressure at level zlev
    lev_press = 0.5*(a_vert(zlev)+a_vert(zlev+1))*ref_press + 0.5*(b_vert(zlev)+b_vert(zlev+1))*dom%surf_press%elts(id+1)

    ! Normalized pressure
    eta = lev_press/dom%surf_press%elts(id+1)
    eta_v = (eta - eta_0) * MATH_PI/2

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*column_mass
    
    ! Potential temperature
    pot_temp =  set_temp(x_i) * (lev_press/ref_press)**(-kappa)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  real(8) function set_temp (x_i)
    implicit none
    type(Coord) :: x_i

    real(8) :: cs2, lon, lat, sn2, Tmean

    call cart2sph (x_i, lon, lat)
    sn2 = sin(lat)**2
    cs2 = cos(lat)**2

    if (eta>=eta_t) then
       Tmean = T_0*eta**(R_d*Gamma_T/grav_accel)
    else
       Tmean = T_0*eta**(R_d*Gamma_T/grav_accel) + delta_T * (eta_t - eta)**5
    end if

    set_temp = Tmean + 0.75_8 * eta*MATH_PI*u_0/R_d * sin(eta_v) * sqrt(cos(eta_v)) * &
         (2*u_0*cos(eta_v)**1.5*(-2*sn2**3*(cs2+1/3.0_8) + 10/63.0_8) + radius*omega*(8/5.0_8*cs2**1.5*(sn2+2/3.0_8) - MATH_PI/4))
  end function set_temp

  real(8) function surf_geopot (x_i)
    ! Surface geopotential
    implicit none
    Type(Coord) :: x_i
    real(8)     :: c1, cs2, lon, lat, sn2

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos(lat)**2
    sn2 = sin(lat)**2

    c1 = u_0*cos((1.0_8-eta_0)*MATH_PI/2)**1.5
    surf_geopot =  c1*(c1*(-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8)  + radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
  end function surf_geopot

  real(8) function surf_pressure (x_i)
    ! Surface pressure
    implicit none
    type(Coord) :: x_i

    surf_pressure = ref_press
  end function surf_pressure

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v
    
    real(8) :: rgrc

    ! Great circle distance
    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))

    u = u_0*cos(eta_v)**1.5*sin(2*lat)**2 + u_p*exp__flush(-(rgrc/R_pert)**2)  ! Zonal velocity component
    v = 0.0_8         ! Meridional velocity component
  end subroutine vel_fun

 subroutine set_thresholds
    ! Set thresholds dynamically (trend or sol must be known)
    use wavelet_mod
    implicit none
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: lnorm, threshold_new

    character(3), parameter :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       if (adapt_trend) then
          call cal_lnorm (trend, order, lnorm)
       else
          call cal_lnorm (sol,   order, lnorm)
       end if
       threshold_new = tol * lnorm
    end if
    
    if (istep /= 0) then
       threshold = 0.9_8*threshold + 0.1_8*threshold_new
    else
       threshold = threshold_new
    end if
  end subroutine set_thresholds

  subroutine initialize_a_b_vert
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = dble(k-1)/dble(zlevels) * press_infty/ref_press
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
       a_vert = a_vert(zlevels+1:1:-1)
       b_vert = b_vert(zlevels+1:1:-1)
    end if
    
    ! Set pressure at infinity
    press_infty = a_vert(zlevels+1)*ref_press ! note that b_vert at top level is 0, a_vert is small but non-zero

    ! Set mass coefficients
    b_vert_mass = b_vert(1:zlevels)-b_vert(2:zlevels+1)
    a_vert_mass = ((a_vert(1:zlevels)-a_vert(2:zlevels+1))*ref_press + b_vert_mass*press_infty)/grav_accel
  end subroutine initialize_a_b_vert

  subroutine read_test_case_parameters (filename)
    implicit none
    character(*)   :: filename

    integer, parameter :: fid = 500
    real(8)            :: press_save
    character(255)     :: varname

    open(unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, compressible
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, uniform
    read (fid,*) varname, remap
    read (fid,*) varname, min_allowed_mass
    read (fid,*) varname, adapt_trend
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, perfect
    read (fid,*) varname, tol
    read (fid,*) varname, optimize_grid
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, cfl_num
    read (fid,*) varname, press_save
    read (fid,*) varname, Laplace_order
    read (fid,*) varname, tau_diffusion
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume

    allocate (pressure_save(1))
    pressure_save(1) = 1.0d2*press_save

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ***********************************************************'
       write (6,'(A,A)')      "test_case           = ", trim (test_case)
       write (6,'(A,A)')      "run_id              = ", trim (run_id)
       write (6,'(A,L1)')     "compressible        = ", compressible
       write (6,'(A,i3)')     "min_level           = ", min_level
       write (6,'(A,i3)')     "max_level           = ", max_level
       write (6,'(A,i5)')     "number of domains   = ", N_GLO_DOMAIN
       write (6,'(A,i3)')     "zlevels             = ", zlevels
       write (6,'(A,L1)')     "uniform             = ", uniform
       write (6,'(A,L1)')     "remap               = ", remap
       write (6,'(A,es10.4)') "min_allowed_mass    = ", min_allowed_mass
       write (6,'(A,L1)')     "adapt_trend         = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds  = ", default_thresholds
       write (6,'(A,L1)')     "perfect             = ", perfect
       write (6,'(A,es10.4)') "tolerance           = ", tol
       write (6,'(A,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt            = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", press_save
       write (6,'(A,i1)')     "Laplace_order       = ", Laplace_order
       write (6,'(A,es10.4)') "tau_diffusion (h)   = ", tau_diffusion
       write (6,'(A,es10.4)') "dt_write            = ", dt_write
       write (6,'(A,i6)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(A,es10.4)') "time_end            = ", time_end 
       write (6,'(A,i6)')     "resume              = ", resume
       write (6,'(A)') ' '
    end if
    close(fid)
    dt_write = dt_write * MINUTE
    tau_diffusion = tau_diffusion * HOUR
    time_end = time_end * HOUR
  end subroutine read_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    
    if (rank == 0) then
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i9,4(a,es8.2,1x))') &
            'time [h] = ', time/HOUR, &
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

       write (12,'(5(es15.9,1x),i2,1x,i9,1x,4(es15.9,1x))')  &
            time/HOUR, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    integer :: k
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: lnorm
    
    allocate (threshold(S_MASS:S_VELO,1:zlevels));     threshold     = 0.0_8
    allocate (threshold_def(S_MASS:S_VELO,1:zlevels)); threshold_def = 0.0_8

    lnorm(S_MASS,:) = dPdim/grav_accel
    do k = 1, zlevels
       lnorm(S_TEMP,k) = (a_vert_mass(k) + b_vert_mass(k)*Pdim/grav_accel)*dTempdim
    end do
    lnorm(S_TEMP,:) = lnorm(S_TEMP,:) + Tempdim*lnorm(S_MASS,:) ! Add component due to tendency in mass 
    lnorm(S_VELO,:) = Udim
    if (adapt_trend) lnorm = lnorm/Tdim
    threshold_def = tol * lnorm
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity
    use wavelet_mod
    implicit none
    
    integer :: k
    real(8) :: Area_lozenge, k_max, visc

    allocate (viscosity_divu(1:zlevels))
    
    ! Average area of smallest lozenges
    Area_lozenge = 4*MATH_PI*radius**2/(10*4**max_level + 2)

    ! Smallest triangle edge length
    dx_min = sqrt (Area_lozenge/(sqrt(3.0_8)/2))

    ! Largest wavenumber on regular lozenge grid
    k_max = MATH_PI/(sqrt(3.0_8)*dx_min)

    ! CFL limit for time step
    dt_cfl = cfl_num*dx_min/(wave_speed+Udim)
    dt_init = dt_cfl

    if (fresh_start) then
       ! Viscosity constant from eigenvalues of Laplacian
       if (Laplace_order == 0) then
          viscosity_mass = 0.0_8
          viscosity_temp = 0.0_8
          viscosity_divu = 0.0_8
          viscosity_rotu = 0.0_8
       elseif (Laplace_order == 1 .or. Laplace_order == 2) then
          L_diffusion = L_diffusion / 2**(max_level-min_level) ! Correct length scales for finest grid
 
          viscosity_mass = L_diffusion(1)**(2*Laplace_order) / tau_diffusion
          viscosity_temp = L_diffusion(1)**(2*Laplace_order) / tau_diffusion
          viscosity_divu = L_diffusion(2)**(2*Laplace_order) / tau_diffusion
          viscosity_rotu = L_diffusion(3)**(2*Laplace_order) / tau_diffusion
       elseif (Laplace_order > 2) then
          if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
          stop
       end if
    end if
    visc = max (viscosity_mass, viscosity_temp, maxval (viscosity_divu), viscosity_rotu)
    
    if (rank == 0) then
       write (6,'(4(A,es8.2), A)')   'dx_min  = ', dx_min, ' k_max  = ', k_max, ' dt_cfl = ', dt_cfl, &
            ' diffusion stability constant = ', dt_cfl/dx_min**(2*Laplace_order)*visc, " (should be < 0.25)"
       write (6,'(4(A,es8.2),/)') 'Viscosity_mass = ', viscosity_mass, ' Viscosity_temp = ', viscosity_temp, &
            ' Viscosity_divu = ', sum (viscosity_divu)/zlevels, ' Viscosity_rotu = ', viscosity_rotu
    end if
  end subroutine initialize_dt_viscosity

  subroutine apply_initial_conditions
    implicit none
    integer :: k, l

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -1, 1)
       end do
    end do
  end subroutine apply_initial_conditions

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

  subroutine set_save_level
    ! Determines closest vertical level to desired pressure
    implicit none
    integer :: k
    real(8) :: dpress, lev_press, save_press

    dpress = 1d16; save_zlev = 0
    do k = 1, zlevels
       lev_press = 0.5*ref_press*(a_vert(k)+a_vert(k+1) + b_vert(k)+b_vert(k+1))
       if (abs(lev_press-pressure_save(1)) < dpress) then
          dpress = abs(lev_press-pressure_save(1))
          save_zlev = k
          save_press = lev_press
       end if
    end do
    if (rank==0) write (6,'(/,A,i2,A,f5.1,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate pressure = ", save_press/100, " hPa)"
  end subroutine set_save_level

  subroutine dump (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load
end module test_case_mod
