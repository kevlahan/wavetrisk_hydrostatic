module will1_mod
  use main_mod

  implicit none

  real(8) :: VELO_SCALE
  real(8) :: LENGTH_SCALE
  real(8) :: HEIGHT_SCALE

  real(8) :: lat_c, lon_c
  integer, allocatable :: n_patch(:)

contains
  subroutine apply_initial_conditions()
      integer l
      do l = level_start, level_end
          call apply_onescale(init_sol, l, 0, 1)
      end do
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
      real(8) lon, lat, r
      d = dom%id+1
      id = idx(i, j, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      r = RADIUS*acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon - lon_c))
      if (r .lt. LENGTH_SCALE) then
          sol(S_HEIGHT)%data(d)%elts(id+1) = HEIGHT_SCALE/2.0_8*(1.0_8 + cos((MATH_PI*r)/LENGTH_SCALE))
      else
          sol(S_HEIGHT)%data(d)%elts(id+1) = 1.0e-30_8 ! fewer zero divisions for performance
      end if
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0
  end subroutine

  subroutine will1_dump(fid)
      integer fid
  end subroutine

  subroutine will1_load(fid)
      integer fid
  end subroutine

  subroutine set_threshold() 
      real(8) corolis_param
      real(8) Ro_number
      corolis_param = 2.0_8*omega
      Ro_number = VELO_SCALE/(LENGTH_SCALE*corolis_param)
      toll_height = Ro_number*VELO_SCALE*LENGTH_SCALE*corolis_param*threshold**(3.0_8/2.0_8)
      toll_velo = 1.0e30_8 ! do not refine from velocity
  end subroutine

  subroutine write_and_print_step()
      if (rank .eq. 0) write(*,'(A,F10.5,A,F10.5,A,I8)') &
              'time [h] = ', time/3600.0_8, &
              ',  dt [s] = ', dt, &
              ',  active', n_active(S_HEIGHT)
  end subroutine

  subroutine write_and_export(k)
      integer k, l
      do l = level_start, level_end
          call write_level_mpi(write_primal, 100000+100*k+l, l, .True.)
      end do
  end subroutine

  subroutine vel_fun(lon, lat, u, v)
      real(8) lon, lat
      real(8) u, v
      real(8) :: alpha = 0
      u = VELO_SCALE*(cos(lat)*cos(alpha) + sin(lat)*cos(lon)*sin(alpha))
      v = -VELO_SCALE*sin(lon)*sin(alpha)
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

  subroutine set_velo_at_new_patches()
      integer d, p
      do d = 1, size(grid)
          do p = n_patch(d)+1, grid(d)%patch%length
              call apply_onescale_to_patch(adv_velocity, grid(d), p-1, -1, 1)
          end do
      end do
  end subroutine
end module

program will1
  use main_mod
  use will1_mod
  implicit none

  integer avg_active

  logical :: aligned
  real(8) ierr
  integer d

  call init_main_mod()
  call read_test_case_parameters("will1.in")
  VELO_SCALE = (2.0_8*MATH_PI*radius)/(12.0_8*DAY)
  LENGTH_SCALE = radius/3.0_8
  HEIGHT_SCALE = 1000.0_8
  optimize_grid = HR_GRID
  advect_only = .True.
  lon_c = (3.0_8*MATH_PI)/2.0_8
  lat_c = 0.0_8
  call set_threshold()
  if (rank .eq. 0) write (*,*) 'threshold',  toll_height
  call initialize(apply_initial_conditions, 1, set_threshold, will1_dump, will1_load)

  allocate(n_patch(size(grid)))
  n_patch = 1
  call set_velo_at_new_patches()

  avg_active = n_active(S_HEIGHT)
  call write_and_export(cp_idx)
  do while (time .lt. time_end)
      n_patch = grid(:)%patch%length
      call time_step(dt_write, aligned)
      call write_and_print_step()
      if (aligned) then
          avg_active = avg_active + n_active(S_HEIGHT)
          call write_and_export(cp_idx+1)
          ierr = writ_checkpoint(will1_dump)
          deallocate(n_patch)
          call restart_full(set_threshold, will1_load)
          allocate(n_patch(size(grid)))
          n_patch = 1
      end if
      call set_velo_at_new_patches()
  end do
  avg_active = avg_active / (cp_idx)

  if (rank .eq. 0) write(*,*) 'Average number of active nodes:', avg_active
  call finalize()
end program
