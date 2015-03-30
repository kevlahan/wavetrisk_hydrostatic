include 'gauss_quad.f90'

module turb_mod
  use main_mod
  use gauss_quad

  implicit none

  real(8), parameter :: VELO_SCALE = 80.0_8
  real(8) :: length_scale
  real(8) :: corolis_param
  real(8) :: Ro_number

  real(8), parameter :: GAL_H0   = 10158.2_8/8.0_8
  real(8), parameter :: GAL_ALPHA = 1.0_8/3.0_8
  real(8), parameter :: GAL_BETA  = 1.0_8/15.0_8
  real(8), parameter :: GAL_H_HAT = 120.0_8
  real(8) :: GAL_E_N  
  real(8) :: GAL_PHI0
  real(8) :: GAL_PHI1
  real(8) :: GAL_PHI2
  real(8) :: GAL_LON_PERP = MATH_PI/8.0_8

  integer iwrite

contains
  subroutine apply_initial_conditions()
      integer l
      do l = level_start, level_end
          call init_ring(MATH_PI*0.0_8, MATH_PI*0.16_8)      
          call apply_onescale(init_sol, l, 0, 1)
          call init_ring(MATH_PI*0.1_8, MATH_PI*0.26_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
          call init_ring(MATH_PI*0.2_8, MATH_PI*0.36_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
          call init_ring(MATH_PI*0.3_8, MATH_PI*0.46_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
      
          call init_ring(-MATH_PI*0.16_8, -MATH_PI*0.0_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
      
          call init_ring(-MATH_PI*0.26_8, -MATH_PI*0.1_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
          call init_ring(-MATH_PI*0.36_8, -MATH_PI*0.2_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
          call init_ring(-MATH_PI*0.46_8, -MATH_PI*0.3_8)
          call apply_onescale2(add_sol, l, 0, 0)
          call apply_to_pole(add_sol, l, 0, .True.)
      end do
  end subroutine
  subroutine read_test_case_parameters(filename)
      character(*) filename
      integer :: fid = 500
      character(255) varname
      open(unit=fid, file=filename, action='READ')
      read(fid,*) varname, max_level
      read(fid,*) varname, viscosity
      read(fid,*) varname, threshold
      read(fid,*) varname, time_end
      read(fid,*) varname, dt_write
      read(fid,*) varname, resume
      read(fid,*) varname, optimize_grid
      close(fid)
  end subroutine
  subroutine write_and_print_step()
      real(8) timing, energ1, energ_pot, enstr
!     enstr = integrate_tri(pot_enstr)
      energ1 = integrate_hex(energy,level_start)
      energ_pot = integrate_hex(pot_energy,level_start)
      timing = get_timing()
      if (rank .eq. 0) write(*,'(A,F10.5,A,F10.5,A,I8,I8)') &
              'time [h] = ', time/3600.0_8, &
              ',  dt [s] = ', dt, &
              ',  active :', sum(n_active)
      if (rank .eq. 0) write(1011,'(E16.9, I3, 2(1X, I9), 2(1X, E16.8), 1X, F16.7)') &
              time, level_end, n_active, energ1, energ_pot, timing
!     if (rank .eq. 0) write(1011,'(E16.9, I3, 2(1X, I9), 3(1X, E16.8), 1X, F16.7)') &
!             time, level_end, n_active, enstr, energ1, energ_pot, timing
  end subroutine
  subroutine write_grid_and_vort(k)
      integer l, k
      call trend_ml(sol, trend)
      do l = level_start, level_end
!         call write_level_mpi(write_hex, 100000+100*k+l, l, .True.)
          call write_level_mpi(write_vort, 200000+100*k+l, l, .False.)
      end do
  end subroutine
  subroutine fill_up_grid_and_IWT()
      integer old_level_start
      old_level_start = level_start
      do while (level_start .lt. level_end)
          if (rank .eq. 0) write(*,*) 'Filling up level', level_start+1
          call fill_up_level()
      end do
      call invers_wavelet_transform(sol, old_level_start)
  end subroutine
  subroutine init_ring(phi0,phi1)
      real(8) phi0, phi1
      GAL_PHI0 = phi0
      GAL_PHI1 = phi1
      GAL_PHI2 = 0.5_8*(phi1+phi0)
      GAL_E_N  = exp(-4.0_8/(GAL_PHI1-GAL_PHI0)**2)
  end subroutine
  subroutine vel_fun(lon, lat, u, v)
    real(8), intent(out) ::  u, v
    real(8), intent(in) ::  lon, lat
    u = zonal_vel(lat)
    v = 0.0_8
  end subroutine
  real(8) function zonal_vel(lat)
    real(8) lat
    if (GAL_PHI0 .lt. lat .and. lat .lt. GAL_PHI1) then
        zonal_vel = VELO_SCALE/GAL_E_N * exp__flush(1.0_8/((lat-GAL_PHI0)*(lat-GAL_PHI1)))
    else 
        zonal_vel = 0.0_8
    end if
  end function
  function height_fun(lon, lat)
    real(8) ::  height_fun
    real(8) lon, lat, integral
    if (lat .le. GAL_PHI0) then
        integral = 0.0_8
    else if (GAL_PHI0 .lt. lat .and. lat .lt. GAL_PHI1) then
        integral = eval_gauss(GAL_PHI0, lat);
    else  !(lat .gt. GAL_PHI1) 
        integral = eval_gauss(GAL_PHI0, GAL_PHI1);
    end if
    height_fun = GAL_H0 - integral/GRAV_ACCEL &
               + GAL_H_HAT*cos(lat) &
                 *exp__flush(-((lon-GAL_LON_PERP)/GAL_ALPHA)**2) &
                 *exp__flush(-((GAL_PHI2-lat)/GAL_BETA)**2)
  end function
  real(8) function eval_gauss(a, b)
    real(8) a, b
    integer i
    eval_gauss = 0.5_8*(b-a) * sum(GAUSS_WEIGHTS * (/ &
                 (vel_integrant(0.5_8*(b-a)*GAUSS_ABSC(i) + 0.5_8*(b+a)), i = 1,size(GAUSS_ABSC)) &
                 /))
  end function
  real(8) function vel_integrant(lat)
    real(8) lat
    real(8) u
    u = zonal_vel(lat)
    vel_integrant = u*(RADIUS*2.0_8*OMEGA*sin(lat) + tan(lat)*u)
  end function
  subroutine init_sol(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      integer idN
      integer idE
      integer idNE
      real(8) lon
      real(8) lat
      d = dom%id+1
      id = idx(i, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      sol(S_HEIGHT)%data(d)%elts(id+1) = height_fun(lon, lat)
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idE+1))
      sol(S_VELO)%data(d)%elts(DG+EDGE*id+1) = proj_vel(vel_fun, dom%node%elts(idNE+1), &
              dom%node%elts(id+1))
      sol(S_VELO)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idN+1))
  end subroutine
  subroutine add_sol(dom, p, i, j, offs, dims)
      type(Domain) dom
      integer p
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      integer idN
      integer idE
      integer idNE
      real(8) lon
      real(8) lat
      d = dom%id+1
      id = idx(i, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      sol(S_HEIGHT)%data(d)%elts(id+1) = sol(S_HEIGHT)%data(d)%elts(id+1) + height_fun(lon, lat)
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1) = sol(S_VELO)%data(d)%elts(EDGE*id+RT+1) + &
              proj_vel(vel_fun, dom%node%elts(id+1), dom%node%elts(idE+1))
      sol(S_VELO)%data(d)%elts(DG+EDGE*id+1) = sol(S_VELO)%data(d)%elts(DG+EDGE*id+1) + &
              proj_vel(vel_fun, dom%node%elts(idNE+1), dom%node%elts(id+1))
      sol(S_VELO)%data(d)%elts(EDGE*id+UP+1) = sol(S_VELO)%data(d)%elts(EDGE*id+UP+1) + &
              proj_vel(vel_fun, dom%node%elts(id+1), dom%node%elts(idN+1))
  end subroutine
  subroutine turb_dump(fid)
      integer fid
      write(fid) iwrite
  end subroutine
  subroutine turb_load(fid)
      integer fid
      read(fid) iwrite
  end subroutine
  subroutine set_thresholds()
      toll_velo  = Ro_number*velo_scale*threshold**(3.0_8/2.0_8)
      toll_height = Ro_number*velo_scale*length_scale*corolis_param*threshold**(3.0_8/2.0_8)
  end subroutine

end module
program turb
  use main_mod
  use turb_mod

  implicit none

  logical aligned
  real(8) last_checkpoint

  integer k_p, k_v, recommanded_level_start, ii, ierr
  character(8+8+29) command

  call init_main_mod()

  call read_test_case_parameters("turb.in")
  length_scale = MATH_PI*radius*0.16_8
  corolis_param = 2.0_8*omega*sqrt(2.0_8)/2.0_8
  Ro_number = VELO_SCALE/(length_scale*corolis_param)

  iwrite = 0
  call set_thresholds()
  if (rank .eq. 0) write (*,*) 'thresholds p, u:',  toll_height, toll_velo
  call initialize(apply_initial_conditions, 1, set_thresholds, turb_dump, turb_load)

  recommanded_level_start = write_active_per_level()
! do while (level_start .lt. recommanded_level_start)
!     call fill_up_level()
! end do
! call invers_wavelet_transform(sol, min_level)
! call apply_onescale__int(set_masks, level_start, -2, 2, ADJZONE)
! dt = cpt_dt_mpi() ! to count active nodes
! recommanded_level_start = write_active_per_level()

  if (rank .eq. 0) write(*,*) 'Write initial values and grid'
  ii = 0
! call write_grid_and_vort(iwrite)

  if (istep .eq. 0) then
      if (rank .eq. 0) open(unit=1011,file='verlauf.out')
      call write_and_print_step()
  else
      if (rank .eq. 0) then
          ! append to verlauf.out after restart
          write(command, '(A,I8,A)') "head -n ", istep, " verlauf.out > tmp && mv -f tmp verlauf.out"
          call system(command) 
          open(unit=1011,file='verlauf.out',access='APPEND')
      end if
  end if
  last_checkpoint = time
  do while (time .lt. time_end)
      call start_timing()
      call time_step(dt_write, aligned)
      call stop_timing()
      call write_and_print_step()
      if (istep .eq. 20) then
          close(1011) ! close and reopen to force write
          open(unit=1011,file='verlauf.out',access='APPEND')
      end if
      if (time-last_checkpoint .gt. 3600.0_8) then
          if (rank .eq. 0) write(*,'(A, I4, A, L1, A)') &
                  'clean and checkpoint', cp_idx+1, ' (exceed = ', max_level_exceeded, ')'      
          close(1011) ! close and reopen to force write
!         ierr = writ_checkpoint(turb_dump)
!         last_checkpoint = time
!         call restart_full(set_thresholds, turb_load)
          open(unit=1011,file='verlauf.out',access='APPEND')
      end if
      if (aligned) then
          iwrite = iwrite + 1
          if (rank .eq. 0) write(*,*) 'write grid and vorticity'
          call write_grid_and_vort(iwrite)
      end if
  end do
  close(1011)
  call finalize()
end program

