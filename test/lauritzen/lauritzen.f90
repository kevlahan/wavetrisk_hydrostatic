module lauritzen_mod
  use main_mod

  implicit none

  real(8) :: VELO_SCALE
  real(8) :: LENGTH_SCALE
  real(8) :: HEIGHT_SCALE

  real(8) :: lat_1, lon_1
  real(8) :: lat_2, lon_2

contains
  subroutine apply_initial_conditions()
      integer l
      call apply(init_sol)
  end subroutine

  subroutine read_test_case_parameters(filename)
      character(*) filename
      integer :: fid = 500
      character(255) varname
      open(unit=fid, file=filename, action='READ')
      read(fid,*) varname, max_level
      read(fid,*) varname, threshold
      read(fid,*) varname, time_end
      read(fid,*) varname, dt_write
      close(fid)
  end subroutine

  subroutine init_sol(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      real(8) lon, lat
      d = dom%id+1
      id = idx(i, j, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      sol(S_HEIGHT)%data(d)%elts(id+1) = height_fun(lon, lat)
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0
  end subroutine

  subroutine substract_init_cond(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      real(8) lon, lat
      d = dom%id+1
      id = idx(i, j, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      sol(S_HEIGHT)%data(d)%elts(id+1) = sol(S_HEIGHT)%data(d)%elts(id+1) - height_fun(lon, lat)
  end subroutine

  subroutine lauritzen_dump(fid)
      integer fid
  end subroutine

  subroutine lauritzen_load(fid)
      integer fid
  end subroutine

  subroutine set_threshold() 
      toll_height = threshold
      toll_velo = 1.0e30_8 ! do not refine from velocity
  end subroutine

  subroutine write_and_print_step()
      if (rank .eq. 0) write(*,'(A,F10.5,A,F10.5,A,I8)') &
              'time [h] = ', time, &
              ',  dt [s] = ', dt, &
              ',  active', n_active(S_HEIGHT)
  end subroutine

  subroutine write_and_export(k)
      integer k, l
      do l = level_start, level_end
          call write_level_mpi(write_primal, 100000+100*k+l, l, .True.)
      end do
  end subroutine

  real(8) function height_fun(lon, lat)
      real(8) lon, lat
      real(8) r1, r2
      r1 = RADIUS*acos(sin(lat_1)*sin(lat) + cos(lat_1)*cos(lat)*cos(lon - lon_1))
      r2 = RADIUS*acos(sin(lat_2)*sin(lat) + cos(lat_2)*cos(lat)*cos(lon - lon_2))
      if (r1 .lt. LENGTH_SCALE) then
          height_fun = HEIGHT_SCALE/2.0_8*(1.0_8 + cos((MATH_PI*r1)/LENGTH_SCALE))
      else if (r2 .lt. LENGTH_SCALE) then
          height_fun = HEIGHT_SCALE/2.0_8*(1.0_8 + cos((MATH_PI*r2)/LENGTH_SCALE))
      else
          height_fun = 1.0e-30_8 ! fewer zero divisions for performance
      end if
  end function

  subroutine vel_fun(lon, lat, u, v)
      real(8) lon, lat
      real(8) u, v
      u = VELO_SCALE * sin(lon/2.0_8)**2 * sin(2.0_8*lat) * cos(MATH_PI*time/time_end)
      v = VELO_SCALE/2.0_8 * sin(lon) * cos(lat) * cos(MATH_PI*time/time_end)
  end subroutine

  subroutine adv_velocity(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      integer idN
      integer idE
      integer idNE
      d = dom%id+1
      id = idx(i, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idE+1))
      sol(S_VELO)%data(d)%elts(DG+EDGE*id+1) = proj_vel(vel_fun, dom%node%elts(idNE+1), &
              dom%node%elts(id+1))
      sol(S_VELO)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, dom%node%elts(id+1), &
              dom%node%elts(idN+1))
  end subroutine

  subroutine set_velo_at_all_patches()
      integer d, p
      do d = 1, size(grid)
          do p = 1+1, grid(d)%patch%length
              call apply_onescale_to_patch(adv_velocity, grid(d), p-1, -1, 1)
          end do
      end do
  end subroutine
end module

program lauritzen
  use main_mod
  use lauritzen_mod
  implicit none

  integer avg_active

  logical :: aligned
  real(8) ierr
  integer d

  call init_main_mod()
  call read_test_case_parameters("lauritzen.in")
  VELO_SCALE = 2.4_8
  LENGTH_SCALE = 0.5_8
  HEIGHT_SCALE = 1.0_8
  radius = 1
  optimize_grid = HR_GRID
  advect_only = .True.
  lon_1 = MATH_PI
  lat_1 = MATH_PI/3.0_8
  lon_2 = MATH_PI
  lat_2 = -MATH_PI/3.0_8
  call set_threshold()
  if (rank .eq. 0) write (*,*) 'threshold',  toll_height
  call initialize(apply_initial_conditions, 1, set_threshold, lauritzen_dump, lauritzen_load)

  call set_velo_at_all_patches()

  avg_active = n_active(S_HEIGHT)
  call write_and_export(cp_idx)
  do while (time .lt. time_end)
      call time_step(dt_write, aligned)
      call write_and_print_step()
      if (aligned) then
          avg_active = avg_active + n_active(S_HEIGHT)
          call write_and_export(cp_idx+1)
          ierr = writ_checkpoint(lauritzen_dump)
          call restart_full(set_threshold, lauritzen_load)
      end if
      call set_velo_at_all_patches()
  end do
  avg_active = avg_active / (cp_idx)

  if (rank .eq. 0) write(*,*) 'Writing out error'
  call apply(substract_init_cond)
  call write_and_export(cp_idx+10)

  if (rank .eq. 0) write(*,*) 'Average number of active nodes:', avg_active
  call finalize()
end program
