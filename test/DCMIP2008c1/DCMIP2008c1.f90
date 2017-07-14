!TO DO: swap out kinetic energy (wind.f90 DYNAMICO)
! add coriolis correction to timestep calculation
! fix ARK method

module DCMIP2008c1_mod ! DCMIP2008 test case 5 parameters
  use main_mod
  implicit none

  ! all dimensional, subscript _t denotes test case values in contrast to those set in shared.f90
  real(8), parameter :: grav_accel_t = 9.80616_8 !gravitational acceleration in meters per second squared
  real(8), parameter :: omega_t      = 7.29212e-5 !Earth’s angular velocity in radians per second
  real(8), parameter :: f0_t         = 2.0_8*omega_t !Coriolis parameter
  real(8), parameter :: u_0_t        = 35.0_8 !velocity in meters per second
  real(8), parameter :: T_0_t        = 288.0_8 !temperature in Kelvin
  real(8), parameter :: DeltaT_t     = 4.8e+05 !empirical temperature difference in the atmosphere in Kelvin
  real(8), parameter :: d2_t         = (1500.0e+03)**2 !square of half width of Gaussian mountain profile in meters
  real(8), parameter :: h_0_t        = 2000.0_8 !mountain height in meters
  real(8), parameter :: a_t          = 6.371229e6 !mean radius of the Earth in meters
  real(8), parameter :: ref_press_t  = 100000.6_8 !reference pressure (mean surface pressure) in Pascals
  real(8), parameter :: R_d_t        = 287.04_8 !ideal gas constant for dry air in joules per kilogram Kelvin
  real(8), parameter :: c_p_t        = 1004.64_8 !specific heat at constant pressure in joules per kilogram Kelvin
  real(8), parameter :: kappa_t      = R_d_t/c_p_t !kappa=R_d/c_p
  real(8), parameter :: N_t          = sqrt(grav_accel_t*grav_accel_t/(c_p_t*T_0_t)) !Brunt-Vaisala buoyancy frequency
  real(8), parameter :: cst_density_t= 100.0_8 !density for the incompressible case in kilogram per metres cubed
  real(8)            :: press_infty_t !pressure at the top of the model in Pascals, is set later using a's and b's
  real(8), parameter :: eta_t        = 0.2_8 !non-dimensional location of tropopause
  real(8), parameter :: eta_0        = 0.252_8 !non-dimensional constant to do with eta
  real(8), parameter :: Gamma_T      = 0.005_8 !temperature lapse rate in Kelvin per metres

  ! Dimensional scaling
  real(8), parameter :: Ldim = sqrt(d2_t)                             ! horizontal length scale
  real(8), parameter :: Hdim = h_0_t                                  ! vertical length scale
  real(8), parameter :: Udim = u_0_t                                  ! velocity scale
  real(8), parameter :: acceldim = Udim*Udim/Hdim                     ! acceleration scale
  real(8), parameter :: Tdim = Ldim/Udim                              ! time scale
  real(8), parameter :: Tempdim = T_0_t                               ! temperature scale (both theta and T from DYNAMICO)
  real(8), parameter :: pdim = ref_press_t                            ! pressure scale
  real(8), parameter :: R_ddim = R_d_t                                ! R_d scale
  real(8), parameter :: massdim = pdim*Hdim/(Tempdim*R_d_t)           ! mass (=rho*dz following DYNAMICO) scale
  real(8), parameter :: specvoldim = (R_d_t*Tempdim)/pdim             ! specific volume scale
  real(8), parameter :: geopotdim = acceldim*massdim*specvoldim/Hdim  ! geopotential scale

  real(8) :: csq

  real(8) :: VELO_SCALE

  real(8), parameter :: LAND = 1
  real(8), parameter :: SEA  = 0
  character(255) IC_file

  integer :: CP_EVERY 

  real(8) :: Hmin, eta, alpha, dh_min, dh_max, dx_min, dx_max, kmin
  real(8) :: initotalmass, totalmass, timing, total_time, dh, mean_mass(1:18), mean_temp(1:18), nr_nodes

  logical const_bathymetry, wasprinted

  real(8) max_dh

  integer iwrite, j

contains
  subroutine apply_initial_conditions()
    integer l, d, p, k

    wasprinted=.false.
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale(init_sol, l, k, 0, 1)
          wasprinted=.false.
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine sum_total_mass(initialgo)
    integer k
    logical initialgo

    k=1 !select vertical level
    if (initialgo) then
       initotalmass=integrate_hex(mass_pert, level_start, k)
    else
       totalmass=integrate_hex(mass_pert, level_start, k)
       !        if (rank.eq.0) write(*,'(A,ES23.14)') 'integr_hex relative change in mass', abs(totalmass-initotalmass)/initotalmass
    end if
  end subroutine sum_total_mass

  subroutine write_and_print_step()
    real(4) timing
    timing = get_timing()
    if (rank .eq. 0) write(1011,'(3(ES13.4,1X), I3, 2(1X, I9), 1(1X,ES13.4))') &
         time, dt, timing, level_end, n_active, VELO_SCALE
  end subroutine write_and_print_step

  subroutine cpt_max_dh(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .gt. 0) then
       if (dom%level%elts(id+1) .eq. level_end .or. dom%mask_n%elts(id+1) .eq. ADJZONE) &
            max_dh = max(max_dh, abs(sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)))
    end if
  end subroutine cpt_max_dh

  subroutine initialize_a_b_vert()
    allocate(a_vert(zlevels+1), b_vert(zlevels+1))

    if (zlevels.eq.18) then
       a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
            0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
            0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
            0.0_8, 0.0_8 /)
       a_vert=a_vert(19:1:-1) ! DCMIP order is opposite ours
       b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
            0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
            0.953189_8, 0.985056_8, 1.0_8/)
       b_vert=b_vert(19:1:-1)
    elseif (zlevels.eq.49) then
       a_vert=(/0.002251865_8, 0.003983890_8, 0.006704364_8, 0.01073231_8, 0.01634233_8, 0.02367119_8, &
            0.03261456_8, 0.04274527_8, 0.05382610_8, 0.06512175_8, 0.07569850_8, 0.08454283_8, &
            0.08396310_8, 0.08334103_8, 0.08267352_8, 0.08195725_8, 0.08118866_8, 0.08036393_8, &
            0.07947895_8, 0.07852934_8, 0.07751036_8, 0.07641695_8, 0.07524368_8, 0.07398470_8, &
            0.07263375_8, 0.07118414_8, 0.06962863_8, 0.06795950_8, 0.06616846_8, 0.06424658_8, &
            0.06218433_8, 0.05997144_8, 0.05759690_8, 0.05504892_8, 0.05231483_8, 0.04938102_8, &
            0.04623292_8, 0.04285487_8, 0.03923006_8, 0.03534049_8, 0.03116681_8, 0.02668825_8, &
            0.02188257_8, 0.01676371_8, 0.01208171_8, 0.007959612_8, 0.004510297_8, 0.001831215_8, &
            0.0_8, 0.0_8 /)
       a_vert=a_vert(50:1:-1)
       b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
            0.006755112_8, 0.01400364_8, 0.02178164_8, 0.03012778_8, 0.03908356_8, 0.04869352_8, &
            0.05900542_8, 0.07007056_8, 0.08194394_8, 0.09468459_8, 0.1083559_8, 0.1230258_8, &
            0.1387673_8, 0.1556586_8, 0.1737837_8, 0.1932327_8, 0.2141024_8, 0.2364965_8, &
            0.2605264_8, 0.2863115_8, 0.3139801_8, 0.3436697_8, 0.3755280_8, 0.4097133_8, &
            0.4463958_8, 0.4857576_8, 0.5279946_8, 0.5733168_8, 0.6219495_8, 0.6741346_8, &
            0.7301315_8, 0.7897776_8, 0.8443334_8, 0.8923650_8, 0.9325572_8, 0.9637744_8, &
            0.9851122_8, 1.0_8/)
       b_vert=b_vert(50:1:-1)
    else
       write(0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
       stop
    end if
  end subroutine initialize_a_b_vert

  subroutine init_sol(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, d, idN, idE, idNE
    real(8) lon, lat, lev_press, lev_spec_volume, mean_geopot, pert_geopot, eta_c, eta_v
    real(8) pert_z, mean_z, total_z, topo_z !all layer thicknesses (not absolute z coordinates)
    !we use eta_c instead of just eta (cause eta is already defined) as the nondimensional vertical coordinate
    !eta_c is one near the bottom of the domain and 0 near the top (counter-intuitively)

    d = dom%id+1
    id = idx(i, j, offs, dims)
    idN = idx(i, j + 1, offs, dims)
    idE = idx(i + 1, j, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    call cart2sph(dom%node%elts(id+1), lon, lat) !JEMF: check that we need not multiply by MATH_PI

    lev_press = (a_vert(zlev+1)+b_vert(zlev+1))*ref_press_t !see eq (117) in NCAR paper
    PRINT *, 'lev_press', lev_press

    !preliminary (unused) calculation of total_z from (92) in NCAR_ASP_2008_idealized_testcases_29May08.pdf
    !this is the z coordinate of the interface starting from the topography
    total_z = log(ref_press_t/lev_press)*T_0_t*R_d_t/grav_accel_t
    PRINT *, 'total_z', total_z

    eta_c = lev_press/ref_press_t
    PRINT *, 'eta_c', eta_c

    eta_v = (eta_c - eta_0) * MATH_PI / 2.0_8
    PRINT *, 'eta_v', eta_v

    !!! calculate masses
    !define mean geopotential so we can determine mass distribution
    if (eta_c .gt. eta_t) then
        mean_geopot = (T_0_t*grav_accel_t/Gamma_t)*(1.0_8-eta_c**(R_d_t*Gamma_t/grav_accel_t))
        PRINT *, 'active'
    else
        mean_geopot = (T_0_t*grav_accel_t/Gamma_t)*(1.0_8-eta_c**(R_d_t*Gamma_t/grav_accel_t)) - &
            R_d_t*DeltaT_t*((log(eta_c/eta_t)+137.0_8/60.0_8)*eta_t**5-5.0_8*eta_c*eta_t**4+ &
            5.0_8*eta_t**3*eta_c**2-(10.0_8/3.0_8)*eta_t**2*eta_c**3+(5.0_8/4.0_8)*eta_t*eta_c**4-0.2_8*eta_c**5)
    end if

    !define pertubation geopotential so we can determine mass distribution
    pert_geopot = u_0_t*cos(eta_v)**(3.0_8/2.0_8)* &
        ((-2.0_8*(sin(lat)**6)*(cos(lat)**2+1.0_8/3.0_8)+10.0_8/63.0_8)*u_0_t*(cos(eta_v)**(3.0_8/2.0_8)) &
        +(8.0_8/5.0_8*(cos(lat)**3)*((sin(lat)**2)+2.0_8/3.0_8)-MATH_PI/4.0_8)*a_t*Omega_t)

    !define surface geopotential
    dom%surf_geopot%elts(id+1) = u_0_t*cos(eta_v)**(3.0_8/2.0_8)* &
        ((-2.0_8*(sin(lat)**6)*(cos(lat)**2+1.0_8/3.0_8)+10.0_8/63.0_8)*u_0_t*(cos(eta_v)**(3.0_8/2.0_8)) &
        +(8.0_8/5.0_8*(cos(lat)**3)*((sin(lat)**2)+2.0_8/3.0_8)-MATH_PI/4.0_8)*a_t*Omega_t)

    !define associated surface topography height
    topo_z = dom%surf_geopot%elts(id+1)/grav_accel_t
    PRINT *, 'topo_z', topo_z

    !now we set the heights (not masses yet)
    if (zlev.eq.1) then !bottom layer
        mean_z = (mean_geopot-dom%surf_geopot%elts(id+1))/grav_accel_t
        pert_z = pert_geopot/grav_accel_t
        PRINT *, 'mean height', mean_z, 'pert height', pert_z, 'sum is', mean_z + pert_z
        dom%spec_vol%elts(id+1)=2.0_8*kappa_t*c_p_t*T_0_t/lev_press !JEMF: probably need to interpolate denominator
        sol(S_MASS,zlev)%data(d)%elts(id+1)=pert_z/dom%spec_vol%elts(id+1)
        mean(S_MASS,zlev)=mean_z/dom%spec_vol%elts(id+1)
    else !other layers need differencing
        mean_z = (mean_geopot-dom%surf_geopot%elts(id+1))/grav_accel_t - dom%adj_mass%elts(id+1)
        pert_z = pert_geopot/grav_accel_t
        PRINT *, 'mean height', mean_z, 'pert height', pert_z, 'sum is', mean_z + pert_z
        dom%spec_vol%elts(id+1)=2.0_8*kappa_t*c_p_t*T_0_t/lev_press !JEMF: probably need to interpolate denominator
        sol(S_MASS,zlev)%data(d)%elts(id+1)=pert_z/dom%spec_vol%elts(id+1)
        mean(S_MASS,zlev)=mean_z/dom%spec_vol%elts(id+1)
    end if

    !!! calculate temperatures
    !set real, non-mass weighted mean temperature from (8) and (9)
    if (eta_c .gt. eta_t) then
        mean(S_TEMP,zlev)=T_0_t*eta_c**(R_d_t*Gamma_t/grav_accel_t)
    else
        mean(S_TEMP,zlev)=T_0_t*eta_c**(R_d_t*Gamma_t/grav_accel_t)+DeltaT_t*((eta_t-eta_c)**5)
    end if

    !set real, non-mass weighted perturbation temperature from (10)
    sol(S_TEMP,zlev)%data(d)%elts(id+1)=0.75_8*(eta_c*MATH_PI*u_0_t/R_d_t)*sin(eta_v)*sqrt(cos(eta_v))* &
        ((-2.0_8*(sin(lat)**6)*(cos(lat)**2+1.0_8/3.0_8)+10.0_8/63.0_8)*2.0_8*u_0_t*(cos(eta_v)**(3.0_8/2.0_8)) &
        +(8.0_8/5.0_8*(cos(lat)**3)*((sin(lat)**2)+2.0_8/3.0_8)-MATH_PI/4.0_8)*a_t*Omega_t)

    !convert non-mass-weighted temperatures to potential temperatures
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_TEMP,zlev)%data(d)%elts(id+1)*(lev_press / ref_press_t)**(-kappa_t)
    mean(S_TEMP,zlev) = mean(S_TEMP,zlev)*(lev_press / ref_press_t)**(-kappa_t)

    !now make these quantities mass-weighted, note that the perturbation comprises three terms
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1)*sol(S_TEMP,zlev)%data(d)%elts(id+1)+ &
        mean(S_MASS,zlev)*sol(S_TEMP,zlev)%data(d)%elts(id+1)+sol(S_MASS,zlev)%data(d)%elts(id+1)*mean(S_TEMP,zlev)
    mean(S_TEMP,zlev) = mean(S_TEMP,zlev)*mean(S_MASS,zlev)

    !!! calculate velocities
    !set initial velocity field
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel_eta(vel_fun, dom%node%elts(id+1), &
         dom%node%elts(idE+1), eta_v)
    sol(S_VELO,zlev)%data(d)%elts(DG+EDGE*id+1) = proj_vel_eta(vel_fun, dom%node%elts(idNE+1), &
         dom%node%elts(id+1), eta_v)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel_eta(vel_fun, dom%node%elts(id+1), &
         dom%node%elts(idN+1), eta_v)

    dom%adj_mass%elts(id+1) = mean_z !adj_mass is used to save adjacent mean_geopot
    !dom%adj_temp%elts(id+1) = lev_press !adj_temp is used to save adjacent lev_press
  end subroutine init_sol

  subroutine vel_fun(lon, lat, u, v, eta_v) !lon plays role of eta_v
    real(8) lon, lat
    real(8) u, v
    real(8) eta_v
    u = u_0_t*(cos(eta_v)**(3.0_8/2.0_8))*(sin(2.0_8*lat)**2)
    v = 0.0_8
  end subroutine vel_fun

  subroutine read_test_case_parameters(filename)
    character(*) filename
    integer :: fid = 500
    character(255) varname
    open(unit=fid, file=filename, action='READ')
    read(fid,*) varname, max_level
    read(fid,*) varname, zlevels
    read(fid,*) varname, threshold 
    read(fid,*) varname, optimize_grid 
    read(fid,*) varname, const_bathymetry 
    read(fid,*) varname, Hmin 
    read(fid,*) varname, eta 
    read(fid,*) varname, alpha 
    read(fid,*) varname, dt_write
    read(fid,*) varname, CP_EVERY
    read(fid,*) varname, time_end
    read(fid,*) varname, resume

    if (rank.eq.0) then
       write(*,'(A,i3)')     "max_level        = ", max_level
       write(*,'(A,i3)')     "zlevels          = ", zlevels
       write(*,'(A,es11.4)') "threshold        = ", threshold
       write(*,'(A,i2)')     "optimize_grid    = ", optimize_grid 
       write(*,'(A,L3)')     "const_bathymetry = ", const_bathymetry
       write(*,'(A,es11.4)') "Hmin             = ", Hmin
       write(*,'(A,es11.4)') "eta              = ", eta
       write(*,'(A,es11.4)') "alpha            = ", alpha
       write(*,'(A,es11.4)') "dt_write         = ", dt_write
       write(*,'(A,i3)')     "CP_EVERY         = ", CP_EVERY
       write(*,'(A,es11.4)') "time_end         = ", time_end 
       write(*,'(A,i6)')     "resume           = ", resume
       write(*,*) ' '
    end if
    dt_write = dt_write * 60_8!/Tdim
    time_end = time_end * 60_8**2!/Tdim

    allocate (mean(S_MASS:S_VELO,1:zlevels))
    call initialize_a_b_vert()
    close(fid)
  end subroutine read_test_case_parameters

  subroutine write_and_export(k)
    integer l, k, zlev
    integer u, i

    call trend_ml(sol, trend)
    call pre_levelout()

    zlev = zlevels ! export only one vertical level

    do l = level_start, level_end
       minv = 1.d63;
       maxv = -1.d63;
       u = 100000+100*k

       call write_level_mpi(write_primal, u+l, l, zlev, .True.)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
       if (rank .eq. 0) write(u,'(A, 4(E15.5E2, 1X), I3)') &
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
       if (rank .eq. 0) write(u,'(A, 4(E15.5E2, 1X), I3)') &
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
       u = 200000+100*k
    end do

    call post_levelout()
    call barrier
    !if (rank .eq. 0) call compress_files(k)
  end subroutine write_and_export

  subroutine cart2sph2(cin, cout)
    type(Coord) cin
    real(8), intent(out) :: cout(2)
    call cart2sph(cin, cout(1), cout(2))
  end subroutine cart2sph2

  subroutine DCMIP2008c1_dump(fid)
    integer fid
    write(fid) VELO_SCALE
    write(fid) iwrite
  end subroutine DCMIP2008c1_dump

  subroutine DCMIP2008c1_load(fid)
    integer fid
    read(fid) VELO_SCALE
    read(fid) iwrite
  end subroutine DCMIP2008c1_load

  subroutine set_thresholds() ! inertia-gravity wave
    tol_mass = VELO_SCALE * c_p/grav_accel * threshold**(1.5_8)
    tol_velo = VELO_SCALE                        * threshold**(1.5_8)
    tol_temp = tol_mass
  end subroutine set_thresholds
end module DCMIP2008c1_mod

program DCMIP2008c1
  use main_mod
  use DCMIP2008c1_mod
  implicit none

  integer, parameter :: len_cmd_files = 12 + 4 + 12 + 4
  integer, parameter :: len_cmd_archive = 11 + 4 + 4
  character(len_cmd_files) cmd_files
  character(len_cmd_archive) cmd_archive
  character(9+len_cmd_archive) command1
  character(6+len_cmd_files) command2

  integer j_gauge, k, l, d
  logical aligned
  character(8+8+29+14) command
  integer ierr, num
  logical write_init

  call init_main_mod()
  call read_test_case_parameters("DCMIP2008c1.in")

  ! for this testcase, set the pressure at infinity (which is usually close to zero)
  press_infty_t = a_vert(zlevels+1)*ref_press_t ! note that b_vert at top level is 0, a_vert is small but non-zero

  ! Shared non-dimensional parameters, these are set AFTER those in shared.f90
  omega = omega_t !* Tdim
  grav_accel = grav_accel_t !/ acceldim
  radius = a_t !/ Ldim
  press_infty = press_infty_t !/ pdim
  R_d = R_d_t !/ R_d_t
  c_p = c_p_t !/ R_d_t
  kappa = kappa_t
  ref_press = ref_press_t !/ pdim
  cst_density = cst_density_t !not non-dimensionalized at present, NOT USED HERE

  PRINT *, 'massdim=', massdim
  PRINT *, 'Tempdim=', Tempdim
  PRINT *, 'Tdim=', Tdim

  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8)) ! Average minimum grid size
  dx_max = 2.0_8*MATH_PI * radius

  kmin = MATH_PI/dx_max ; kmax = 2.0_8*MATH_PI/dx_max

  csq = grav_accel*a_t !JEMF
  c_p = sqrt(csq)
  VELO_SCALE   = grav_accel*dh/sqrt(csq)  ! Characteristic velocity based on initial perturbation !JEMF must set dh

  wind_stress      = .False.
  penalize         = .False.
  bottom_friction  = .False.
  const_bathymetry = .True.
  compressible     = .True.

  if (rank.eq.0) then
     write(*,'(A,L1)') "wind_stress      = ", wind_stress
     write(*,'(A,L1)') "penalize         = ", penalize
     write(*,'(A,L1)') "bottom friction  = ", bottom_friction
     write(*,'(A,L1)') "const_bathymetry = ", const_bathymetry
     write(*,'(A,L1)') "compressible = ", compressible
     if (rank .eq. 0) write (*,*) 'running without bathymetry and continents'
  end if

  viscosity = 0.0_8 !1.0_8/((2.0_8*MATH_PI/dx_min)/64.0_8)**2     ! grid scale viscosity
  friction_coeff = 3e-3_8 ! Bottom friction
  if (rank .eq. 0) write (*,'(A,es11.4)') 'Viscosity = ',  viscosity

  write_init = (resume .eq. NONE)
  iwrite = 0

  call initialize(apply_initial_conditions, 1, set_thresholds, DCMIP2008c1_dump, DCMIP2008c1_load)
  call sum_total_mass(.True.)

  !PRINT *, 'mean_mass', mean_mass
  !PRINT *, 'mean_temp', mean_temp
  !PRINT *, 'nr_nodes', nr_nodes
  !PRINT *, 'real nr_nodes', nr_nodes/18
  !    PRINT *, 'real dimensional mean_mass', mean_mass/(nr_nodes/18)
  !PRINT *, 'real dimensional mean_temp', mean_temp/(nr_nodes/18)
  !stop

  if (rank .eq. 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:',  tol_mass, tol_temp, tol_velo
  call barrier()

  if (rank .eq. 0) write(6,*) 'Write initial values and grid'
  if (write_init) call write_and_export(iwrite)

  total_time = 0_8
  do while (time .lt. time_end)
     ! Set thresholds dynamically
     max_dh = 0
     do l = level_start, level_end
        do k = 1, zlevels
           call apply_onescale(cpt_max_dh, l, k, 0, 1) !compute averages here
        end do
     end do
     max_dh = sync_max_d(max_dh)
     VELO_SCALE = max(VELO_SCALE*0.99, min(VELO_SCALE, grav_accel*max_dh/c_p))
     call set_thresholds()

     call start_timing()
     call time_step(dt_write, aligned)
     call stop_timing()
     timing = get_timing()
     total_time = total_time + timing

     call write_and_print_step()

     if (rank .eq. 0) write(*,'(A,F9.5,A,F9.5,2(A,ES13.5),A,I9,A,ES11.4)') &
          'time [h] =', time/3600.0_8, &
          ', dt [s] =', dt, &
          ', min. depth =', fd, &
          ', VELO_SCALE =', VELO_SCALE, &
          ', d.o.f. =', sum(n_active), &
          ', cpu = ', timing

     call print_load_balance()

     if (aligned) then
        iwrite = iwrite + 1
        call write_and_export(iwrite)
        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle
        ierr = writ_checkpoint(DCMIP2008c1_dump)

        ! let all cpus exit gracefully if NaN has been produced
        ierr = sync_max(ierr)
        if (ierr .eq. 1) then ! NaN
           write(0,*) "NaN when writing checkpoint"
           call finalize()
           stop
        end if

        call restart_full(set_thresholds, DCMIP2008c1_load)
        call barrier()
     end if

     call sum_total_mass(.False.)

     !call remap_coordinates()
  end do

  if (rank .eq. 0) then
     write(6,'(A,ES11.4)') 'Total cpu time = ', total_time
     close(1011)
     close(8450)
  end if

  call finalize()
end program DCMIP2008c1
