!to compile low-res: go to wavetrisk_hydrostatic and execute
!                   make BUILD_DIR=build TEST_CASE=DCMIP2008c5 PARAM=param_J5
!to run: go to wavetrisk_hydrostatic/jobs/DCMIP2008c5 and execute
!                   ../../bin/DCMIP2008c5

!check steady state case
!check one-layer case which should reduce to shallow water calculation
!are pentagon points the culprit?
!euler step for time integration

module DCMIP2008c5_mod
  use main_mod

  implicit none

  ! DCMIP2008 test case 5 parameters, all dimensional, subscript _t denotes test case values in contrast
  ! to those set in shared.f90
  real(8), parameter :: grav_accel_t = 9.80616_8 !gravitational acceleration in meters per second squared
  real(8), parameter :: omega_t      = 7.29211e-5 !Earth’s angular velocity in radians per second
  real(8), parameter :: f0_t         = 2.0_8*omega_t !Coriolis parameter
  real(8), parameter :: u_0_t        = 0.0_8 !20.0_8 !velocity in meters per second
  real(8), parameter :: T_0_t        = 288.0_8 !temperature in Kelvin
  real(8), parameter :: d2_t         = (1500.0e+03)**2 !square of half width of Gaussian mountain profile in meters
  real(8), parameter :: h_0_t        = 0.0_8 !2000.0_8 !mountain height in meters
  real(8), parameter :: lon_c_t      = MATH_PI/2.0_8 !mountain peak longitudinal location in radians
  real(8), parameter :: lat_c_t      = MATH_PI/6.0_8 !mountain peak latitudinal location in radians
  real(8), parameter :: p_sp_t       = 930.0e2 !South Pole surface pressure in Pascals
  real(8), parameter :: a_t          = 6.371229e6 !mean radius of the Earth in meters
  real(8), parameter :: ref_press_t  = 100145.6_8 !reference pressure (mean surface pressure) in Pascals
  real(8), parameter :: R_d_t        = 287.04_8 !ideal gas constant for dry air in joules per kilogram Kelvin
  real(8), parameter :: c_p_t        = 1004.64_8 !specific heat at constant pressure in joules per kilogram Kelvin
  real(8), parameter :: kappa_t      = R_d_t/c_p_t !kappa=R_d/c_p
  real(8), parameter :: N_t          = sqrt(grav_accel_t*grav_accel_t/(c_p_t*T_0_t)) !Brunt-Vaisala buoyancy frequency
  real(8), parameter :: cst_density_t= 100.0_8 !density for the incompressible case in kilogram per metres cubed
  real(8)            :: press_infty_t !pressure at the top of the model in Pascals, is set later using a's and b's

  ! Dimensional scaling
  real(8), parameter :: Ldim = sqrt(d2_t)                                ! horizontal length scale
  real(8), parameter :: Hdim = 100.0_8 !h_0_t                                     ! vertical length scale
  real(8), parameter :: Udim = 20.0_8 !u_0_t                                     ! velocity scale
  real(8), parameter :: Tdim = Ldim/Udim                                 ! time scale
  real(8), parameter :: Tempdim = T_0_t                                  ! temperature (=T, not theta, from DYNAMICO) scale
  real(8), parameter :: pdim = ref_press_t                               ! pressure scale
  real(8), parameter :: R_ddim = R_d_t                                   ! R_d scale
  real(8), parameter :: massdim = pdim*Hdim/(Tempdim*R_d_t)              ! mass (=rho*dz following DYNAMICO) scale
  real(8), parameter :: specvoldim = (R_d_t*Tempdim)/pdim                ! specific volume scale
  real(8), parameter :: geopotdim = (Hdim/(Udim*Udim))*massdim*specvoldim*Hdim ! geopotential scale

  ! Non-dimensional parameters
  real(8) :: VELO_SCALE !probably remove

  character(255) IC_file

  integer :: CP_EVERY, iwrite, j

  real(8) :: Hmin, eta, alpha, dh_min, dh_max, the_area

  real(8), allocatable :: a_vert(:), b_vert(:), the_mass(:), prev_mass(:)

  logical const_bathymetry, wasprinted

contains
  subroutine apply_initial_conditions()
    integer l, k

    ! check if the default values for the physical parameters in shared.f90 were properly reset
    !if ((abs(omega-1.0_8).lt.1.0e-10).or.(abs(grav_accel-1.0_8).lt.1.0e-10).or.(abs(radius-1.0_8).lt.1.0e-10).or. &
    !        (abs(press_infty-1.0_8).lt.1.0e-10).or.(abs(R_d-1.0_8).lt.1.0e-10).or.(abs(c_p-1.0_8).lt.1.0e-10).or.&
    !        (abs(kappa-1.0_8).lt.1.0e-10).or.(abs(ref_press-1.0_8).lt.1.0e-10)) then
    !    write(0,*) "One of the physical parameters' default values of 1.0_8 in shared.f90 were not reset properly"
    !    stop
    !end if

    ! set the initial conditions
    wasprinted=.false.
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale(init_sol, l, k, 0, 1)
          call apply_onescale(nondim_sol, l, k, 0, 1)
          wasprinted=.false.
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine check_initial_conditions()
    integer k, d, p, l
    do k = zlevels, 1, -1
        do d = 1, size(grid)
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch(show_grid_data, grid(d), p - 1, k, 0, 0)
          end do
       end do
    end do
  end subroutine check_initial_conditions

  subroutine write_and_print_step()
    real(4) timing
    timing = get_timing()
    if (rank .eq. 0) write(1011,'(3(ES13.4,1X), I3, 2(1X, I9), 1(1X,ES13.4))') &
         time, dt, timing, level_end, n_active, VELO_SCALE
  end subroutine write_and_print_step

  subroutine initialize_a_b_vert()
    allocate(a_vert(zlevels+1), b_vert(zlevels+1), the_mass(zlevels), prev_mass(zlevels))

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
    real(8) lon, lat, rgrc, lev_z, lev_press, lev_spec_volume

    d = dom%id+1
    id = idx(i, j, offs, dims)
    idN = idx(i, j + 1, offs, dims)
    idE = idx(i + 1, j, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    call cart2sph(dom%node%elts(id+1), lon, lat) !lambda is longitude, phi is latitude (DCMIP2008)

    !set surface geopotential (note that this really only needs to be done once (improve this later))
    rgrc = a_t*acos(sin(lat_c_t)*sin(lat)+cos(lat_c_t)*cos(lat)*cos(lon-lon_c_t))
    dom%surf_geopot%elts(id+1) = grav_accel_t*h_0_t * exp(-rgrc*rgrc/d2_t)

    !set initial surface pressure; it is given and does not need to be computed via a downward integration
    dom%surf_press%elts(id+1) = p_sp_t*exp(-((a_t*N_t*N_t*u_0_t)/(2.0_8*grav_accel_t*grav_accel_t*kappa_t))* &
        (u_0_t/a_t+2.0_8*omega_t)*(sin(lat)*sin(lat)-1.0_8)-((N_t*N_t)/(grav_accel_t*grav_accel_t*kappa_t))* &
        dom%surf_geopot%elts(id+1))

    !now initialize all the layers
    !set initial mass field
    lev_press=a_vert(zlev+1)*ref_press_t+b_vert(zlev+1)*dom%surf_press%elts(id+1)
    lev_z=log(dom%surf_press%elts(id+1)/lev_press)*T_0_t*R_d_t/grav_accel_t !see (92) in NCAR_ASP_2008_idealized_testcases_29May08.pdf

    if (zlev.eq.1) then
        sol(S_MASS,zlev)%data(d)%elts(id+1) = lev_z
        dom%spec_vol%elts(id+1)=2.0_8*kappa_t*c_p_t*T_0_t/(lev_press+dom%surf_press%elts(id+1))
    elseif (zlev.eq.zlevels) then
        sol(S_MASS,zlev)%data(d)%elts(id+1) = lev_z - dom%adj_mass%elts(id+1)
        dom%spec_vol%elts(id+1)=2.0_8*kappa_t*c_p_t*T_0_t/(lev_press+press_infty_t)
    else
        sol(S_MASS,zlev)%data(d)%elts(id+1) = lev_z - dom%adj_mass%elts(id+1)
        dom%spec_vol%elts(id+1)=2.0_8*kappa_t*c_p_t*T_0_t/(lev_press+dom%adj_temp%elts(id+1))
    end if
    !error of integration: deltap^3/(12*g*p^2)

    !print out representative layer thicknesses
    if (.not.wasprinted) then
        if (rank .eq. 0) write(*,'(A,I2,A,F8.2,A,F9.2,A,F9.2)') &
            'zlev=', zlev, &
            ',   thickness (in metres)=', sol(S_MASS,zlev)%data(d)%elts(id+1), &
            ',   z coord (in metres)=', lev_z, &
            ',   layer upper interface pressure (in Pascals)=', lev_press
        wasprinted=.true.
    end if

    !note that sol(S_MASS) is not exactly mass
    !it is pseudo-density in the DYNAMICO paper but rhodz in the DYNAMICO code
    !we choose it to be rho*dz here, as well, as to facilitate vertical integration
    sol(S_MASS,zlev)%data(d)%elts(id+1)=sol(S_MASS,zlev)%data(d)%elts(id+1)/dom%spec_vol%elts(id+1)

    !set initial velocity field
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
        dom%node%elts(idE+1))
    sol(S_VELO,zlev)%data(d)%elts(DG+EDGE*id+1) = proj_vel(vel_fun, dom%node%elts(idNE+1), &
        dom%node%elts(id+1))
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
        dom%node%elts(idN+1))

    !set initial mass-weighted potential temperature; uniform
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = T_0_t*sol(S_MASS,zlev)%data(d)%elts(id+1)* &
            (lev_press / ref_press_t)**(-kappa_t) !mass weighted potential temperature
    !PRINT *,' potential temperature', sol(S_TEMP,zlev)%data(d)%elts(id+1)
    !PRINT *, (lev_press / ref_press)**(-kappa)
    !PRINT *, lev_press / ref_press
    !PRINT *, '--'

    dom%adj_mass%elts(id+1) = lev_z !adj_mass is used to save adjacent lev_z
    dom%adj_temp%elts(id+1) = lev_press !adj_temp is used to save adjacent lev_press
  end subroutine init_sol

  subroutine nondim_sol(dom, i, j, zlev, offs, dims) !nondimensionalize the initial conditions
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, d, idN, idE, idNE

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (zlev.eq.1) then
        dom%surf_geopot%elts(id+1) = dom%surf_geopot%elts(id+1)/geopotdim
    end if

    sol(S_MASS,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1)/massdim

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1)/Udim
    sol(S_VELO,zlev)%data(d)%elts(DG+EDGE*id+1) = sol(S_VELO,zlev)%data(d)%elts(DG+EDGE*id+1)/Udim
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1)/Udim

    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_TEMP,zlev)%data(d)%elts(id+1)/Tempdim
  end subroutine nondim_sol

  subroutine show_grid_data(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, d, idN, idE, idNE

    d = dom%id+1
    id = idx(i, j, offs, dims)

    PRINT *, 'dom%surf_press%elts(id+1)', dom%surf_press%elts(id+1)

  end subroutine show_grid_data

  subroutine read_test_case_parameters(filename)
    !read in test case parameters
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

    dt_write = dt_write * 60_8/Tdim
    time_end = time_end * 60_8**2/Tdim

    call initialize_a_b_vert()
    close(fid)
  end subroutine read_test_case_parameters

  subroutine vel_fun(lon, lat, u, v)
      real(8) lon, lat
      real(8) u, v
      u = u_0_t*cos(lat)
      v = 0.0_8
  end subroutine

  subroutine total_mass(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id = idx(i, j, offs, dims)

    !if (dom%mask_n%elts(id+1) .gt. 0) then !vertical summation in the column
        the_mass(zlev)=the_mass(zlev)+sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)/dom%areas%elts(id+1)%hex_inv
    !end if
  end subroutine total_mass

  subroutine sphere_area(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id

    id = idx(i, j, offs, dims)

    the_area=the_area+1.0_8/dom%areas%elts(id+1)%hex_inv
  end subroutine sphere_area

  subroutine write_and_export(k,zlev)
    integer l, k, zlev
    integer u, i

    !zlev = 1 ! export only one vertical level
    write(*,*) 'writing and exporting primal grid data for zlev=', zlev

    !call ctrend(sol, trend) !JEMF

    call pre_levelout()

    do l = level_start, level_end
       minv = 1.d63;
       maxv = -1.d63;
       u = 10000000+10000*k+100*l+zlev

       call write_level_mpi(write_primal, u, l, zlev, .True.)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
    end do

    call post_levelout()
    call barrier
    if (rank .eq. 0) call compress_files(k) 
  end subroutine write_and_export

  subroutine cart2sph2(cin, cout)
    type(Coord) cin
    real(8), intent(out) :: cout(2)
    call cart2sph(cin, cout(1), cout(2))
  end subroutine cart2sph2

  subroutine DCMIP2008c5_dump(fid)
    integer fid
    write(fid) VELO_SCALE
    write(fid) iwrite
  end subroutine DCMIP2008c5_dump

  subroutine DCMIP2008c5_load(fid)
    integer fid
    read(fid) VELO_SCALE
    read(fid) iwrite
  end subroutine DCMIP2008c5_load

  subroutine set_thresholds() ! inertia-gravity wave
    tol_mass = VELO_SCALE*c_p/grav_accel * threshold**(3.0_8/2.0_8)
    tol_velo = VELO_SCALE                * threshold**(3.0_8/2.0_8)
    tol_temp = VELO_SCALE*c_p/grav_accel * threshold**(3.0_8/2.0_8) !JEMF: this for now
  end subroutine set_thresholds
end module DCMIP2008c5_mod

program DCMIP2008c5
  use main_mod
  use DCMIP2008c5_mod
  implicit none

  integer, parameter :: len_cmd_files = 12 + 4 + 12 + 4
  integer, parameter :: len_cmd_archive = 11 + 4 + 4
  character(len_cmd_files) cmd_files
  character(len_cmd_archive) cmd_archive
  character(9+len_cmd_archive) command1
  character(6+len_cmd_files) command2

  integer k, l, d
  logical aligned
  character(8+8+29+14) command
  integer ierr, num
  logical write_init

  call init_main_mod()
  call read_test_case_parameters("DCMIP2008c5.in")

  ! for this testcase, set the pressure at infinity (which is usually close to zero)
  press_infty_t = a_vert(zlevels+1)*ref_press_t ! note that b_vert at top level is 0

  ! Shared non-dimensional parameters, these are set AFTER those in shared.f90
  omega = omega_t * Tdim
  grav_accel = grav_accel_t * Hdim/(Udim*Udim)
  radius = a_t / Ldim
  press_infty = press_infty_t / pdim
  R_d = R_d_t / R_d_t
  c_p = c_p_t / R_d_t
  kappa = kappa_t
  ref_press = ref_press_t / pdim
  cst_density = cst_density_t !not non-dimensionalized at present, NOT USED HERE

  PRINT *, 'Tdim=', Tdim

  VELO_SCALE = 1.0_8

  wind_stress      = .False.
  penalize         = .False.
  bottom_friction  = .False.
  compressible  = .True.

  if (rank.eq.0) then
     write(*,'(A,L1)') "wind_stress     = ", wind_stress
     write(*,'(A,L1)') "penalize        = ", penalize
     write(*,'(A,L1)') "bottom friction = ", bottom_friction
     write(*,'(A,L1)') "compressibility = ", compressible
  end if

  write_init = .True. !(resume .eq. NONE)
  iwrite = 0

  if (.not. penalize) then
     const_bathymetry = .True.
     if (rank .eq. 0) write (*,*) 'running without bathymetry and continents'
  else
     if (rank .eq. 0) write (*,*) 'Not supposed to run with penalization or bathymetry'
     call finalize()
     stop
  end if

  call initialize(apply_initial_conditions, 1, set_thresholds, DCMIP2008c5_dump, DCMIP2008c5_load)
  write (*,*) 'DCMIP2008c5.f90 has applied initial conditions'
  !call check_initial_conditions()
  call find_pentagons() !find pentagons on the grid !JEMF: probably move somewhere outside of testfile

  the_mass=0.0_8
  do k = 1, zlevels
    call apply_onescale(total_mass, 5, k, 0, 0)
  end do

       if (rank .eq. 0) write(*,'(A,18(F10.3,1X))') 'mass per layer (J=5) =', the_mass!*massdim !JEMF: 18 is hardcoded
       if (rank .eq. 0) write(*,'(A,(F12.3,1X))') 'total mass (J=5) =', sum(the_mass)!*massdim
prev_mass=the_mass

  if (rank .eq. 0) write (*,*) 'thresholds for mass, velo, temp:',  tol_mass, tol_velo, tol_temp
  if (rank .eq. 0) write(*,*) 'Write initial values and grid'
  if (write_init) call write_and_export(iwrite,2) !write first zlevel to file
        call write_and_export(iwrite,18)

  do while (time .lt. time_end)
     do k = 1, zlevels
     call write_step(0, time, k) !write to terminal some output for the kth vertical level
     end do

     call start_timing()

     do k = 1, zlevels
        call update_bdry(sol(S_MASS,k), NONE)
        call update_bdry(sol(S_TEMP,k), NONE) !JEMF: update velocity?
     end do

     VELO_SCALE = 1.0_8

     call set_thresholds()

     call time_step(dt_write, aligned)

     PRINT *, 'DCMIP2008c5 flag1'

     the_mass=0.0_8
     do k = 1, zlevels
         call apply_onescale(total_mass, 5, k, 0, 0)
     end do

    if (rank .eq. 0) write(*,'(A,18(F10.3,1X))') 'mass per layer (J=5) =', the_mass!*massdim !JEMF: 18 is hardcoded
    if (rank .eq. 0) write(*,'(A,(F12.3,1X))') 'total mass (J=5) =', sum(the_mass)!*massdim
    if (rank .eq. 0) write(*,'(A,18(es10.2,1X))') 'relative change per layer (J=5) =', 100.0_8*(the_mass-prev_mass)/prev_mass !JEMF: 18 is hardcoded
    if (rank .eq. 0) write(*,'(A,18(es10.2,1X))') 'relative total change (J=5) =', 100.0_8*(sum(the_mass)- &
        sum(prev_mass))/sum(prev_mass) !JEMF: 18 is hardcoded

     !prev_mass=the_mass

     call stop_timing()

     call write_and_print_step()

     if (rank .eq. 0) write(*,'(A,F9.5,A,F9.5,A,I9)') &
          '---------- time [h] =', time/3600.0_8*Tdim, &
          ', dt [s] =', dt*Tdim, &
          ', d.o.f. =', sum(n_active)

     call print_load_balance()

     !if (aligned) then
        iwrite = iwrite + 1
        !call write_and_export(iwrite,1)
        call write_and_export(iwrite,2)
        call write_and_export(iwrite,18)
        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle
        ierr = writ_checkpoint(DCMIP2008c5_dump)

        ! let all cpus exit gracefully if NaN has been produced
        ierr = sync_max(ierr)
        if (ierr .eq. 1) then ! NaN
           write(0,*) "NaN when writing checkpoint"
           call finalize()
           stop
        end if

        !call restart_full(set_thresholds, DCMIP2008c5_load)
        call barrier()
     !end if
     !stop
  end do
  if (rank .eq. 0) then
     close(1011)
     close(8450)
  end if
  call finalize()
end program DCMIP2008c5
