module tenlayergauss_mod
  use main_mod
  implicit none

  ! Dimensional physical parameters
  real(8), parameter :: g_star     = 9.80665
  real(8), parameter :: rho        = 1.027e3
  real(8), parameter :: H_star     = 3.344175893265152e+03 ! mean ocean depth
  real(8), parameter :: R_star     = 6371e3
  real(8), parameter :: L_star     = R_star
  real(8), parameter :: Omega_star = 7.29e-5
  real(8), parameter :: tau_star   = 0.0_8
  real(8), parameter :: f0_star    = 2.0*Omega_star

  ! Dimensional scaling
  real(8), parameter :: Ldim = L_star  ! Horizontal length scale 
  real(8), parameter :: Hdim = H_star  ! Vertical length scale
  real(8), parameter :: Udim = sqrt(H_star*g_star); ! Velocity scale is unperturbed wave speed
  real(8), parameter :: Tdim = Ldim/Udim            ! Time scale

  ! Non-dimensional parameters
  real(8), parameter :: f0   = f0_star * Ldim/Udim
  real(8), parameter :: H    = H_star/Hdim
  real(8), parameter :: dh   = 1e-3_8

  real(8) :: csq

  real(8) :: VELO_SCALE

  real(8), parameter :: LAND = 1
  real(8), parameter :: SEA  = 0
  character(255) IC_file

  integer :: CP_EVERY 

  real(8) :: Hmin, eta, alpha, dh_min, dh_max, dx_min, dx_max, kmin, k_tsu
  real(8) :: initotalmass, totalmass, timing, total_time

  logical const_bathymetry

  real(8) max_dh

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
  end subroutine

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

    d = dom%id+1
    id = idx(i, j, offs, dims)

    call cart2sph(dom%node%elts(id+1), lon, lat)

    rgrc = acos(sin(lat_c_t)*sin(lat)+cos(lat_c_t)*cos(lat)*cos(lon-lon_c_t))

    dom%surf_press%elts(id+1) = 0.0_8

    dom%surf_geopot%elts(id+1) = 0.0_8

    sol(S_MASS,zlev)%data(d)%elts(id+1) = 0.0_8

    if (zlev.eq.zlevels) then
       sol(S_MASS,zlev)%data(d)%elts(id+1) = dh*exp(-1e3_8*rgrc*rgrc)
    else
       sol(S_MASS,zlev)%data(d)%elts(id+1) = 0.0_8
    end if
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1)
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
    read(fid) VELO_SCALE
    read(fid) iwrite
  end subroutine tenlayergauss_load

  subroutine set_thresholds() ! inertia-gravity wave
    tol_mass = VELO_SCALE * c_p/grav_accel * threshold**(1.5_8)
    tol_velo = VELO_SCALE                        * threshold**(1.5_8)
    tol_temp = tol_mass
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

  ! Shared non-dimensional parameters
  radius     = R_star / Ldim
  grav_accel = g_star * Hdim/Udim**2
  omega      = Omega_star * Ldim/Udim

  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8)) ! Average minimum grid size
  dx_max = 2.0_8*MATH_PI * radius

  kmin = MATH_PI/dx_max ; kmax = 2.0_8*MATH_PI/dx_max

  csq = grav_accel*H
  c_p = sqrt(csq)
  k_tsu = 2.0_8*MATH_PI/(1e6_8/Ldim) ! Approximate wavelength of tenlayergauss: 100km

  VELO_SCALE   = grav_accel*dh/sqrt(csq)  ! Characteristic velocity based on initial perturbation

  ! Set (non-dimensional) mean values of variables
  allocate (mean(S_MASS:S_VELO,1:zlevels))
  do k = 1, zlevels
     mean(S_MASS,k) = 1.0_8/real(zlevels)
     mean(S_TEMP,k) = mean(S_MASS,k)
     mean(S_VELO,k) = 0.0_8
  end do

  wind_stress      = .False.
  penalize         = .False.
  bottom_friction  = .False.
  const_bathymetry = .True.
  compressible     = .False.

  if (rank.eq.0) then
     write(*,'(A,L1)') "wind_stress      = ", wind_stress
     write(*,'(A,L1)') "penalize         = ", penalize
     write(*,'(A,L1)') "bottom friction  = ", bottom_friction
     write(*,'(A,L1)') "const_bathymetry = ", const_bathymetry
     if (rank .eq. 0) write (*,*) 'running without bathymetry and continents'
  end if

  viscosity = 0.0_8 !1.0_8/((2.0_8*MATH_PI/dx_min)/64.0_8)**2     ! grid scale viscosity
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
     ! Set thresholds dynamically
     max_dh = 0
     do l = level_start, level_end
        do k = 1, zlevels
           call apply_onescale(cpt_max_dh, l, k, 0, 1)
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
          'time [h] =', time/3600.0_8*Tdim, &
          ', dt [s] =', dt*Tdim, &
          ', min. depth =', fd, &
          ', VELO_SCALE =', VELO_SCALE, &
          ', d.o.f. =', sum(n_active), &
          ', cpu = ', timing

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

     call sum_total_mass(.False.)
  end do

  if (rank .eq. 0) then
     write(6,'(A,ES11.4)') 'Total cpu time = ', total_time
     close(1011)
     close(8450)
  end if

  call finalize()
end program tenlayergauss
