module tenlayergauss_mod
  use main_mod
  implicit none

  ! Values for troposphere (roughly)
  real(8), parameter :: grav_accel_t      = 9.80616_8 !gravitational acceleration in meters per second squared
  real(8), parameter :: Height_t             = 11e3_8 ! Assumed atmosphere height
  real(8), parameter :: ref_density_t     = 1.225_8 ! Density of air in kg/m^3
  real(8), parameter :: a_t               = 6.371229e6 !mean radius of the Earth in meters
  real(8), parameter :: omega_t           = 0.0_8 !Earth’s angular velocity in radians per second
  real(8), parameter :: f0_t              = 2.0_8*omega_t !Coriolis parameter
  real(8), parameter :: ref_temperature_t = 288.15_8 !temperature in Kelvin
  real(8), parameter :: R_d_t             = 287.04_8 !ideal gas constant for dry air in joules per kilogram Kelvin
  real(8), parameter :: c_p_t             = 1004.64_8 !specific heat at constant pressure in joules per kilogram Kelvin
  real(8), parameter :: kappa_t           = R_d_t/c_p_t !kappa=R_d/c_p
  
  !real(8), parameter :: ref_press_t  = 100145.6_8 ! reference pressure (mean atmospheric surface pressure)
  
  ! Dimensional scaling
  real(8), parameter :: Ldim = 1.0!a_t  ! Horizontal length scale
  real(8), parameter :: Hdim = 1.0!Height_t  ! Vertical length scale
  real(8), parameter :: Udim = 1.0!sqrt(h_0_t*grav_accel_t); ! Velocity scale is unperturbed wave speed
  real(8), parameter :: Tdim = Ldim/Udim            ! Time scale
  
  real(8), parameter :: acceldim = Udim*Udim/Hdim  ! acceleration scale
  real(8), parameter :: temperature_dim = 1.0!T_0_t   ! temperature scale (both theta and T from DYNAMICO)
  real(8), parameter :: pdim = 1.0!ref_press_t    ! pressure scale
  real(8), parameter :: R_dim = 1.0!R_d_t  ! R_d scale
  real(8), parameter :: density_dim = 1.0!ref_density_t

  ! Non-dimensional parameters
  real(8), parameter :: H = Height_t/Hdim
  real(8) :: ref_temperature, press_infty_t

  real(8) :: csq

  real(8), parameter :: LAND = 1
  real(8), parameter :: SEA  = 0
  character(255) IC_file

  integer :: CP_EVERY 

  real(8) :: Hmin, eta, alpha, dh_min, dh_max, dx_min, dx_max, kmin, k_tsu
  real(8) :: initotalmass, totalmass, timing, total_time
  real(8) :: max_mass, max_temp, max_velo, max_dmass, max_dtemp, mass_scale, temp_scale, velo_scale

  logical const_bathymetry

  integer iwrite, j

contains
  subroutine apply_initial_conditions()
    integer l, d, p, k

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale(init_sol, l, k, 0, 1)
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
       if (rank.eq.0) write(*,'(A,ES23.14)') 'integr_hex relative change in mass', abs(totalmass-initotalmass)/initotalmass
    end if
  end subroutine sum_total_mass

  subroutine write_and_print_step()
    real(4) timing
    timing = get_timing()
    if (rank .eq. 0) write(1011,'(3(ES13.4,1X), I3, 2(1X, I9), 1(1X,ES13.4))') &
         time, dt, timing, level_end, n_active, VELO_SCALE
  end subroutine write_and_print_step

  subroutine cpt_max_vars(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, e

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       max_mass = max(max_mass, abs(sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)))
       max_temp = max(max_temp, abs(sol(S_TEMP,zlev)%data(dom%id+1)%elts(id+1)))
       do e = 1, EDGE
          max_velo  = max(max_velo, abs(sol(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+e)))
       end do
    end if
  end subroutine cpt_max_vars

  subroutine cpt_max_trend(dom, i, j, zlev, offs, dims)
    type(Domain) :: dom
    integer :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id, e

    id = idx(i, j, offs, dims)

    ! Maximum trends
    if (dom%mask_n%elts(id+1) .gt. 0) then
       if (dom%level%elts(id+1) .eq. level_end .or. dom%mask_n%elts(id+1) .eq. ADJZONE) then
          max_mass = max(max_mass, abs(trend(S_MASS,zlev)%data(dom%id+1)%elts(id+1)))
          max_temp = max(max_temp, abs(trend(S_TEMP,zlev)%data(dom%id+1)%elts(id+1)))
          do e = 1, EDGE
             max_velo  = max(max_velo, abs(trend(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+e)))
          end do
       end if
    endif
  end subroutine cpt_max_trend

  subroutine init_sol(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i, j, k, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id, d
    real(8) lon, lat
    real(8) s, t, rgrc
    real(8), parameter :: lon_c_t = MATH_PI/2.0_8
    real(8), parameter :: lat_c_t = MATH_PI/6.0_8
    real(8) :: pot_temp

    real(8), parameter :: width = 0.1_8

    d = dom%id+1
    id = idx(i, j, offs, dims)

    call cart2sph(dom%node%elts(id+1), lon, lat)

    rgrc = acos(sin(lat_c_t)*sin(lat)+cos(lat_c_t)*cos(lat)*cos(lon-lon_c_t))

    dom%surf_press%elts(id+1) = 0.0_8

    dom%surf_geopot%elts(id+1) = 0.0_8

    ! Perturbation to mass and potential temperature
    if (zlev.eq.zlevels) then
       sol(S_MASS,zlev)%data(d)%elts(id+1) = mean(S_MASS,zlev) * max_dmass*exp(-(rgrc/width)**2)
       pot_temp = 0.0_8!mean(S_TEMP,zlev)/mean(S_MASS,zlev) * max_dtemp*exp(-(rgrc/width)**2)
    else
       sol(S_MASS,zlev)%data(d)%elts(id+1) = 0.0_8
       pot_temp = 0.0_8
    end if
    
    !--Peturbation to mass weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1)*pot_temp + &
         sol(S_MASS,zlev)%data(d)%elts(id+1)*mean(S_TEMP,zlev)/mean(S_MASS,zlev) + &
         mean(S_MASS,zlev)*pot_temp
    
    !--Perturbation to velocity
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
  end subroutine init_sol

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
    dt_write = dt_write * 60_8/Tdim
    time_end = time_end * 60_8**2/Tdim
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

  subroutine tenlayergauss_dump(fid)
    integer fid
    write(fid) VELO_SCALE
    write(fid) iwrite
  end subroutine tenlayergauss_dump

  subroutine tenlayergauss_load(fid)
    integer fid
    read(fid) velo_scale
    read(fid) iwrite
  end subroutine tenlayergauss_load

  subroutine set_thresholds(itype)
    ! Scaling for inertia-gravity wave
    integer, optional :: itype
    integer :: l, k

    ! tol_mass = max_mass * sqrt(csq)/radius * threshold !  dmass/T = max_mass c/L
    ! tol_velo = 2.0_8*max_mass*sqrt(grav_accel/H) * c_p/radius * threshold ! U/T = U c/L
    ! tol_temp = tol_mass

    ! Set thresholds dynamically
    if (istep.eq.0) then
       dt = cpt_dt_mpi()
       tol_mass = mass_scale * threshold*dt
       tol_temp = temp_scale * threshold*dt
       tol_velo = velo_scale * threshold*dt
    else
       max_mass = 0.0_8
       max_temp = 0.0_8
       max_velo = 0.0_8
       do l = level_start, level_end
          do k = 1, zlevels
             if (adapt_trend) then
                call apply_onescale(cpt_max_trend, l, k, 0, 1)
             else
                call apply_onescale(cpt_max_vars, l, k, 0, 1)
             end if
          end do
       end do
       mass_scale = sync_max_d(max_mass)
       temp_scale = sync_max_d(max_temp)
       velo_scale = sync_max_d(max_velo)

       tol_mass = mass_scale * threshold
       tol_temp = temp_scale * threshold
       tol_velo = velo_scale * threshold
    end if
  end subroutine set_thresholds
end module tenlayergauss_mod

program tenlayergauss
  use main_mod
  use tenlayergauss_mod
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
  call read_test_case_parameters("tenlayergauss.in")

  wind_stress      = .False.
  penalize         = .False.
  bottom_friction  = .False.
  const_bathymetry = .True.
  compressible     = .True.

  press_infty = 14101_8/pdim
  
  omega           = omega_t * Tdim
  grav_accel      = grav_accel_t / acceldim
  radius          = a_t / Ldim
  press_infty     = press_infty_t / pdim
  ref_density     = ref_density_t/density_dim
  kappa           = kappa_t
  ref_temperature = ref_temperature_t/temperature_dim

  R_d = R_d_t
  c_p = c_p_t

  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8)) ! Average minimum grid size
  dx_max = 2.0_8*MATH_PI * radius

  kmin = MATH_PI/dx_max ; kmax = 2.0_8*MATH_PI/dx_max

  ! Set (non-dimensional) mean values of variables
  allocate (mean(S_MASS:S_VELO,1:zlevels))
  allocate (mean_press(1:zlevels), mean_spec_vol(1:zlevels), mean_exner(1:zlevels))
  allocate (mean_density(1:zlevels), mean_temperature(1:zlevels))

  mean(S_VELO,1:zlevels) = 0.0_8
  
  do k = 1, zlevels
     mean_temperature(k) = ref_temperature - 6.5_8 * real(k-1+0.5) * H/real(zlevels)/1e3_8 ! Lapse rate for troposphere per km
     mean_density(k)     = ref_density     - (ref_density-4.66348e-01_8)/real(zlevels) * real(k-1+0.5) ! Density profile for troposphere
     
     mean(S_MASS,k) = mean_density(k) * H/real(zlevels)
     
     if (.not. compressible) mean(S_TEMP,k) = mean(S_MASS,k)*mean_density(k)/ref_density ! constant density ref_density
  end do
  if (rank.eq.0) then
     write(6,'("Mean masses       = ", 100(es11.4,1x))') mean(S_MASS,:)
     write(6,'("Mean temperatures = ", 100(es11.4,1x))') mean_temperature
     write(6,'("Mean densities    = ", 100(es11.4,1x))') mean_density
     write(6,*) ' '
  end if

  ! Set mean pressure at each vertical level starting at top level
  do k = zlevels, 1, -1
     if (compressible) then ! Compressible case
        if (k .eq. zlevels) then 
           mean_press(k) = press_infty + 0.5_8*grav_accel*mean(S_MASS,k)
        else
           mean_press(k) = mean_press(k+1) + 0.5_8*grav_accel*(mean(S_MASS,k) + mean(S_MASS,k+1))
        end if

        ! Mean surface pressure is reference pressure
        if (k .eq. 1) ref_press = mean_press(k) + 0.5_8*grav_accel*mean(S_MASS,k)
     else ! Incompressible case
        if (k .eq. zlevels) then
           mean_press(k) = press_infty + 0.5_8*grav_accel*mean(S_TEMP,k)
        else 
           mean_press(k) = mean_press(k+1) + 0.5_8*grav_accel*(mean(S_TEMP,k)+mean(S_TEMP,k+1))
        end if

        ! Mean surface pressure is reference pressure
        if (k .eq. 1) ref_press = mean_press(k) + 0.5_8*grav_accel*mean(S_TEMP,k)
     end if
  end do

  if (rank.eq.0) then
     write(6,'(A,100(es11.4,1x))'), 'Mean pressures    = ', mean_press
     write(6,'(A,es11.4,/)'),       'Surface pressure  = ', ref_press
  end if
  
  ! Calculate mean mass-weighted potential temperature (assume constant temperature)
  if (compressible) then
     do k = 1, zlevels
        mean(S_TEMP,k) = mean_temperature(k) * (ref_press/mean_press(k))**kappa * mean(S_MASS,k)
     end do
  end if

  ! Set mean exner and mean specific volume once
  do k = 1, zlevels
     mean_exner(k) = c_p*(mean_press(k)/ref_press)**kappa
     mean_spec_vol(k) = kappa*(mean(S_TEMP,k)/mean(S_MASS,k))*mean_exner(k)/mean_press(k)
  end do
  
  if (compressible) then ! Only used in compressible case
     if (rank.eq.0) then
        write(6,'(A,100(es11.4,1x))')   'Mean exner            = ', mean_exner
        write(6,'(A,100(es11.4,1x),/)') 'Mean specific volumes = ', mean_spec_vol
     end if
  end if

  ! Relative perturbation magnitudes
  max_dmass = 5e-2_8
  max_dtemp = 5e-2_8

  csq = grav_accel*H
  
  k_tsu = 2.0_8*MATH_PI/(1e6_8/Ldim) ! Approximate wavelength of tenlayergauss: 100km

  mass_scale = max_dmass * mean(S_MASS,zlevels)
  temp_scale = mass_scale
  velo_scale = 2.0_8*max_dmass*sqrt(grav_accel/H)/sqrt(csq)  ! Characteristic velocity based on initial perturbation
  
  if (rank.eq.0) then
     write(6,'(A,L1)') "wind_stress      = ", wind_stress
     write(6,'(A,L1)') "penalize         = ", penalize
     write(6,'(A,L1)') "bottom friction  = ", bottom_friction
     write(6,'(A,L1)') "const_bathymetry = ", const_bathymetry
     write(*,'(A,L1)') "compressible     = ", compressible
     write(6,*) ' '
  end if

  viscosity = 0.0_8
  !viscosity = 1.0_8/((2.0_8*MATH_PI/dx_min)/64.0_8)**2     ! grid scale viscosity
  !viscosity = 1e2_8/(2.0_8*MATH_PI/dx_min)**2     ! grid scale viscosity

  friction_coeff = 3e-3_8 ! Bottom friction
  if (rank .eq. 0) write (*,'(A,es11.4)') 'Viscosity = ',  viscosity

  write_init = (resume .eq. NONE)
  iwrite = 0

  call initialize(apply_initial_conditions, 1, set_thresholds, tenlayergauss_dump, tenlayergauss_load)
  call sum_total_mass(.True.)

  if (rank .eq. 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:',  tol_mass, tol_temp, tol_velo
  call barrier()

  if (rank .eq. 0) write(6,*) 'Write initial values and grid'
  if (write_init) call write_and_export(iwrite)

  total_time = 0_8

  do while (time .lt. time_end)
     call start_timing()
     call time_step(dt_write, aligned, set_thresholds)
     call stop_timing()
     timing = get_timing()
     total_time = total_time + timing

     call write_and_print_step()

     if (rank .eq. 0) write(*,'(A,es11.4,A,es9.2,3(A,ES9.2),A,I9,A,ES9.2,1x)') &
          'time [h] = ', time/3600.0_8*Tdim, &
          ' dt [s] = ', dt*Tdim, &
          ' min. depth = ', fd, &
          ' mass scale = ', mass_scale, &
          ' velo scale = ', velo_scale, &
          ' d.o.f. = ', sum(n_active), &
          ' cpu = ', timing

     call print_load_balance()

     if (aligned) then
        iwrite = iwrite + 1
        call write_and_export(iwrite)
        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle
        ierr = writ_checkpoint(tenlayergauss_dump)

        ! let all cpus exit gracefully if NaN has been produced
        ierr = sync_max(ierr)
        if (ierr .eq. 1) then ! NaN
           write(0,*) "NaN when writing checkpoint"
           call finalize()
           stop
        end if

        call restart_full(set_thresholds, tenlayergauss_load)
        call barrier()
     end if
  end do

  if (rank .eq. 0) then
     write(6,'(A,ES11.4)') 'Total cpu time = ', total_time
     close(1011)
     close(8450)
  end if

  call finalize()
end program tenlayergauss
