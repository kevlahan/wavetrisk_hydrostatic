module will2_mod
  use main_mod

  implicit none

  real(8) :: velo_scale
  real(8) :: length_scale
  real(8) :: corolis_param
  real(8) :: Ro_number

  real(8) :: alpha

  real(8) max_height_error, max_height
  real(8) L2_height_error, L2_height
  logical aligned

contains
  real(8) function press_fun(lon, lat)
      real(8) lon, lat
      real(8) gh0
      gh0 = 2.94e4_8
      press_fun = gh0 - (RADIUS*OMEGA*VELO_SCALE + VELO_SCALE**2/2.0_8)* &
                  (-cos(lon)*cos(lat)*sin(alpha) + sin(lat)*cos(alpha))**2
  end function

  subroutine vel_fun(lon, lat, u, v)
      real(8) lon
      real(8) lat
      real(8) u
      real(8) v
      u = VELO_SCALE*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha))
      v = -VELO_SCALE*sin(lon)*sin(alpha)
  end subroutine

  subroutine init_sol(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      integer idN, idE, idNE
      real(8) lon, lat
      d = dom%id+1
      id = idx(i, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      sol(S_HEIGHT)%data(d)%elts(id+1) = press_fun(lon, lat)/grav_accel
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idE+1))
      sol(S_VELO)%data(d)%elts(DG+EDGE*id+1) = proj_vel(vel_fun, dom%node%elts(idNE+1), &
              dom%node%elts(id+1))
      sol(S_VELO)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idN+1))
  end subroutine

  subroutine apply_initial_conditions()
      integer l
      do l = level_start, level_end
           call apply_onescale(init_sol, l, 0, 1)
      end do
  end subroutine

  subroutine read_test_case_parameters(filename)
      character(*) filename
      character(255) varname
      integer :: fid
      open(unit=fid, file=filename, action='READ')
      read(fid,*) varname, max_level
      read(fid,*) varname, threshold
      read(fid,*) varname, alpha
      read(fid,*) varname, time_end
      close(fid)
  end subroutine

  subroutine cpt_height_err(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      real(8) height, height_error, lon, lat
      integer id
      id = idx(i, j, offs, dims)
      if (dom%mask_p%elts(id+1) .lt. ADJZONE) return
      if (dom%mask_p%elts(id+1) .eq. TOLLRNZ .and. dom%level%elts(id+1) .ne. level_end) return
      call cart2sph(dom%node%elts(id+1), lon, lat)
      height = sol(S_HEIGHT)%data(dom%id+1)%elts(id+1) 
      height_error = height - press_fun(lon, lat)/grav_accel
      max_height_error = max(max_height_error, abs(height_error))
      max_height       = max(max_height,       abs(height))
      L2_height_error  = L2_height_error + height_error**2
      L2_height      = L2_height     + height**2
  end subroutine

  subroutine set_thresholds()
      toll_velo  = Ro_number*velo_scale*threshold**(3.0_8/2.0_8)
      toll_height = Ro_number*velo_scale*length_scale*corolis_param*threshold**(3.0_8/2.0_8)
  end subroutine
end module

program will2
  use main_mod
  use will2_mod

  implicit none

  integer k, l, avg_active

  call init_main_mod()

  call read_test_case_parameters("will2.in")
  optimize_grid = HR_GRID
  corolis_param = 1.0_8*omega
  velo_scale = (2.0_8*MATH_PI*radius)/(12.0_8*DAY)
  length_scale = MATH_PI/4.0_8*radius
  Ro_number = VELO_SCALE/(LENGTH_SCALE*corolis_param)

  call initialize(apply_initial_conditions, 1, set_thresholds, default_dump, default_load)
  k = 0
  avg_active = sum(n_active)

  if (rank .eq. 0) open(unit=1011, file="history.out")
  do while (time .lt. time_end)
      call time_step(time_end, aligned)
      max_height_error = 0.0_8
      max_height = 0.0_8
      do l = level_start, level_end
          call apply_onescale(cpt_height_err, l, 0, 0)
      end do
      max_height_error = sync_max_d(max_height_error)/sync_max_d(max_height)
      L2_height_error  = sqrt(sum_real(L2_height_error)/sum_real(L2_height))
      if (rank .eq. 0) write(1011,'(1(E12.5,1X),I9,2(1X,E12.5))') threshold, sum(n_active), max_height_error, L2_height_error
      if (rank .eq. 0) write(*,*) time/3600.0_8, sum(n_active), max_height_error, L2_height_error
      avg_active = avg_active + sum(n_active)
      k = k + 1
  end do
  avg_active = avg_active / (k+1)
  if (rank .eq. 0) close(1011)

  call finalize()
end program
