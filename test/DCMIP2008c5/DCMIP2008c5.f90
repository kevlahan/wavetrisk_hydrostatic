!TO DO: swap out kinetic energy (wind.f90 DYNAMICO)
! add coriolis correction to timestep calculation
! fix ARK method

module DCMIP2008c5_mod ! DCMIP2008 test case 5 parameters
  use main_mod
  implicit none

  ! all dimensional, subscript _t denotes test case values in contrast to those set in shared.f90
  real(8), parameter :: grav_accel_t = 9.80616_8 !gravitational acceleration in meters per second squared
  real(8), parameter :: omega_t      = 7.29211e-5_8 !Earth’s angular velocity in radians per second
  real(8), parameter :: f0_t         = 2.0_8*omega_t !Coriolis parameter
  real(8), parameter :: u_0_t        = 20.0_8 !velocity in meters per second
  real(8), parameter :: T_0_t        = 288.0_8 !temperature in Kelvin
  real(8), parameter :: d2_t         = (1500.0e+03)**2 !square of half width of Gaussian mountain profile in meters
  real(8), parameter :: h_0_t        = 2000.0_8 !mountain height in meters
  real(8), parameter :: lon_c_t      = MATH_PI/2.0_8 !mountain peak longitudinal location in radians
  real(8), parameter :: lat_c_t      = MATH_PI/6.0_8 !mountain peak latitudinal location in radians
  real(8), parameter :: p_sp_t       = 930.0e2 !South Pole surface pressure in Pascals
  real(8), parameter :: a_t          = 6.371229e6 !mean radius of the Earth in meters
  real(8), parameter :: ref_press_t  = 1e5_8!100145.6_8 !reference pressure (mean surface pressure) in Pascals
  real(8), parameter :: R_d_t        = 287.04_8 !ideal gas constant for dry air in joules per kilogram Kelvin
  real(8), parameter :: c_p_t        = 1004.64_8 !specific heat at constant pressure in joules per kilogram Kelvin
  real(8), parameter :: c_v_t        = 717.6_8 ! specfic heat at constant volume c_v = R_d - c_p
  real(8), parameter :: kappa_t      = R_d_t/c_p_t !kappa=R_d/c_p
  real(8), parameter :: N_t          = sqrt(grav_accel_t*grav_accel_t/(c_p_t*T_0_t)) !Brunt-Vaisala buoyancy frequency
  real(8), parameter :: ref_density_t= 100.0_8 !density for the incompressible case in kilogram per metres cubed

  logical :: uniform  ! Uniform or non-uniform grid in pressure

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

  integer         :: CP_EVERY, id, zlev, iwrite, j
  real(8)         :: Hmin, dh_min, dh_max, dx_min, dx_max, kmin
  real(8)         :: initotalmass, totalmass, timing, total_time, dh
  real(8)         :: csq, VELO_SCALE, max_dh
  logical         :: wasprinted
  character (255) :: IC_file
  
contains
  subroutine apply_initial_conditions()
    integer :: l, d, p, k

    wasprinted=.false.
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, 0, 1)
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

  subroutine initialize_a_b_vert()
    integer :: k

    ! Allocate vertical grid parameters
    if (allocated(a_vert)) deallocate(a_vert)
    if (allocated(b_vert)) deallocate(b_vert)
    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    
    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = real(k-1)/real(zlevels) * press_infty/ref_press
          b_vert(k) = 1.0_8 - real(k-1)/real(zlevels)
       end do
    else
       if (zlevels.eq.18) then
          !a_vert=(/0.0_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0.0_8, 0.0_8 /)
          a_vert=a_vert(19:1:-1) ! DCMIP order is opposite ours
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1.0_8 /)
          b_vert=b_vert(19:1:-1)
       elseif (zlevels.eq.26) then
          !a_vert=(/0.0_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0.0_8, 0.0_8 /)
          a_vert=a_vert(27:1:-1)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
                0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
                0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1.0_8 /)
          b_vert=b_vert(27:1:-1)
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
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
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
    end if
  end subroutine initialize_a_b_vert

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    type (Domain)                  :: dom
    integer                        :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    
    integer :: id, d, idN, idE, idNE
    real(8) :: lon, lat, rgrc, lev_press, pot_temp, p_top, p_bot

    d = dom%id+1
    
    id   = idx(i, j, offs, dims)
    idN  = idx(i, j + 1, offs, dims)
    idE  = idx(i + 1, j, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    ! Initialize vertical grid
    call initialize_a_b_vert()
    
    ! Set pressure at infinity
    press_infty = a_vert(zlevels+1)*ref_press ! note that b_vert at top level is 0, a_vert is small but non-zero

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph(dom%node%elts(id+1), lon, lat)

    ! Surface geopotential for Gaussian mountain (note that this really only needs to be done once)
    rgrc = a_t*acos(sin(lat_c_t)*sin(lat)+cos(lat_c_t)*cos(lat)*cos(lon-lon_c_t))
    dom%surf_geopot%elts(id+1) = grav_accel * h_0_t * exp(-rgrc**2/d2_t)

    ! Surface pressure
    dom%surf_press%elts(id+1) = p_sp_t * exp ( &
         - a_t*N_t**2*u_0_t/(2.0_8*grav_accel_t**2*kappa)*(u_0_t/a_t+2.0_8*omega_t)*(sin(lat)**2-1.0_8) &
         - N_t**2/(grav_accel**2*kappa_t)*dom%surf_geopot%elts(id+1) )

    ! Pressure at level zlev
    lev_press = 0.5_8*(a_vert(zlev)+a_vert(zlev+1))*ref_press + 0.5_8*(b_vert(zlev)+b_vert(zlev+1))*dom%surf_press%elts(id+1)

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = &
         ((a_vert(zlev)-a_vert(zlev+1))*ref_press + (b_vert(zlev)-b_vert(zlev+1))*dom%surf_press%elts(id+1))/grav_accel
    
    ! Horizontally uniform potential temperature
    pot_temp =  T_0_t * (lev_press/ref_press)**(-kappa)
    
    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, dom%node%elts(id+1),   dom%node%elts(idE+1))
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel(vel_fun, dom%node%elts(idNE+1), dom%node%elts(id+1))
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, dom%node%elts(id+1),   dom%node%elts(idN+1))

    ! Print out layer thicknesses
    if (rank.eq.0 .and. .not.wasprinted) then
       if(zlev.eq.1) then
          write(6,'(4(A,es11.4))') &
               ' surf geopot =', dom%surf_geopot%elts(id+1), &
               ' surf press =', dom%surf_press%elts(id+1), &
               ' press_infty =', press_infty
       end if

       write(6,'(A,I2,1x,6(A,es11.4))') &
            ' zlev = ', zlev, &
            ' press =', lev_press, &
            ' mu =', sol(S_MASS,zlev)%data(d)%elts(id+1), &
            ' Theta =', sol(S_TEMP,zlev)%data(d)%elts(id+1), &
            ' U = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1), &
            ' V = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1), &
            ' W = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1)
       wasprinted=.true.
    end if
  end subroutine init_sol

  subroutine vel_fun(lon, lat, u, v)
    ! Zonal latitude-dependent wind
    real(8) :: lon, lat
    real(8) :: u, v
    
    u = u_0_t*cos(lat)  ! Zonal velocity component
    v = 0.0_8           ! Meridional velocity component
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
    read(fid,*) varname, dt_write
    read(fid,*) varname, CP_EVERY
    read(fid,*) varname, time_end
    read(fid,*) varname, resume

    if (rank.eq.0) then
       write(*,'(A,i3)')     "max_level        = ", max_level
       write(*,'(A,i3)')     "zlevels          = ", zlevels
       write(*,'(A,es11.4)') "threshold        = ", threshold
       write(*,'(A,i2)')     "optimize_grid    = ", optimize_grid 
       write(*,'(A,es11.4)') "dt_write         = ", dt_write
       write(*,'(A,i6)')     "CP_EVERY         = ", CP_EVERY
       write(*,'(A,es11.4)') "time_end         = ", time_end 
       write(*,'(A,i6)')     "resume           = ", resume
       write(*,*) ' '
    end if
    dt_write = dt_write * 60_8!/Tdim
    time_end = time_end * 60_8**2!/Tdim

    close(fid)
  end subroutine read_test_case_parameters

  subroutine write_and_export(iwrite)
    integer :: iwrite
    integer :: l, k, zlev, d, u, i, p

    call trend_ml(sol, trend)
    call pre_levelout()

    zlev = 6 ! Only export one vertical level

    do l = level_start, level_end
       minv = 1.d63;
       maxv = -1.d63;
       u = 100000+100*iwrite

       ! Calculate pressure and geopotential at vertical level zlev and scale l
       do k = 1, zlev
          do d = 1, size(grid)
             mass => sol(S_MASS,k)%data(d)%elts
             temp => sol(S_TEMP,k)%data(d)%elts
             exner     => fun(F_EXNER,k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
          end do
       end do

       call write_level_mpi(write_primal, u+l, l, zlev, .True.)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
       if (rank .eq. 0) write(u,'(A, 5(E15.5E2, 1X), I3)') &
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
       if (rank .eq. 0) write(u,'(A, 5(E15.5E2, 1X), I3)') &
            "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
       u = 200000+100*iwrite
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
    tol_mass = VELO_SCALE * sqrt(csq)/grav_accel * threshold**(1.5_8)
    tol_velo = VELO_SCALE                        * threshold**(1.5_8)
    tol_temp = tol_mass
  end subroutine set_thresholds
end module DCMIP2008c5_mod

program DCMIP2008c5
  use main_mod
  use DCMIP2008c5_mod
  implicit none

  integer                      :: ierr, k, l, d
  integer, parameter           :: len_cmd_files = 12 + 4 + 12 + 4
  integer, parameter           :: len_cmd_archive = 11 + 4 + 4
  character(len_cmd_files)     :: cmd_files
  character(len_cmd_archive)   :: cmd_archive
  character(8+8+29+14)         :: command
  character(9+len_cmd_archive) :: command1
  character(6+len_cmd_files)   :: command2
  logical                      :: aligned, write_init

  call init_main_mod()
  call read_test_case_parameters("DCMIP2008c5.in")
  
  compressible = .True.
  if (rank.eq.0) write(6,'(A,L1)') "compressible = ", compressible

  ! Shared non-dimensional parameters, these are set AFTER those in shared.f90
  omega = omega_t !* Tdim
  grav_accel = grav_accel_t !/ acceldim
  radius = a_t !/ Ldim
 
  R_d = R_d_t !/ R_d_t
  c_p = c_p_t !/ R_d_t
  kappa = kappa_t
  ref_press = ref_press_t !/ pdim
  ref_density = ref_density_t !not non-dimensionalized at present, NOT USED HERE
      
  ! Set (non-dimensional) mean values of variables
  allocate (mean(S_MASS:S_VELO,1:zlevels))
  allocate (mean_press(1:zlevels), mean_spec_vol(1:zlevels), mean_exner(1:zlevels))
  allocate (mean_density(1:zlevels), mean_temperature(1:zlevels))

  ! Put both mean and perturbations into perturbation term
  mean = 0.0_8
  mean_press = 0.0_8; mean_spec_vol = 0.0_8; mean_exner = 0.0_8; mean_density = 0.0_8; mean_temperature = 0.0_8

  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8)) ! Average minimum grid size
  dx_max = 2.0_8*MATH_PI * radius

  kmin = MATH_PI/dx_max ; kmax = MATH_PI/dx_min

  csq        = grav_accel*a_t 
  VELO_SCALE = grav_accel*dh/sqrt(csq)  ! Characteristic velocity based on initial perturbation !JEMF must set dh

  !viscosity = 1.0_8/(kmax*40.0_8)**2     ! grid scale viscosity
  viscosity = 0.0_8
  if (rank .eq. 0) write (6,'(A,es11.4)') 'Viscosity  = ',  viscosity
  if (rank .eq. 0) write (6,*) ' '

  write_init = (resume .eq. NONE)
  iwrite = 0
  uniform = .false.

  call initialize(apply_initial_conditions, 1, set_thresholds, DCMIP2008c5_dump, DCMIP2008c5_load)
  
  call sum_total_mass(.True.)

  if (rank .eq. 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:',  tol_mass, tol_temp, tol_velo
  call barrier()

  if (rank .eq. 0) write(6,*) 'Write initial values and grid'
  if (write_init) call write_and_export(iwrite)

  total_time = 0_8
  do while (time .lt. time_end)
     call set_thresholds()

     call start_timing()
     call time_step(dt_write, aligned, set_thresholds)
     call remap_vertical_coordinates()
     call stop_timing()
     timing = get_timing()
     total_time = total_time + timing

     call write_and_print_step()

     if (rank .eq. 0) write(*,'(A,es11.4,A,es11.4,2(A,ES11.4),A,I9,A,ES9.2)') &
          'time [h] = ', time/3600.0_8, &
          ', dt [s] = ', dt, &
          ', min. depth = ', fd, &
          ', VELO_SCALE = ', VELO_SCALE, &
          ', d.o.f. = ', sum(n_active), &
          ', cpu = ', timing

     call print_load_balance()
     
     if (aligned) then
        iwrite = iwrite + 1
!        call remap_vertical_coordinates()
        call write_and_export(iwrite)

        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle
        ierr = writ_checkpoint(DCMIP2008c5_dump)

        ! let all cpus exit gracefully if NaN has been produced
        ierr = sync_max(ierr)
        if (ierr .eq. 1) then ! NaN
           write(0,*) "NaN when writing checkpoint"
           call finalize()
           stop
        end if

        call restart_full(set_thresholds, DCMIP2008c5_load)
        call barrier()
     end if

     call sum_total_mass(.False.)
  end do

  if (rank .eq. 0) then
     write(6,'(A,ES11.4)') 'Total cpu time = ', total_time
     close(1011)
     close(8450)
  end if

  call finalize()
end program DCMIP2008c5
