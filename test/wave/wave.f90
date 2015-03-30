module wave_mod
  use main_mod

  implicit none
  real(8) :: WAVE_SPEED
  real(8) :: VELO_SCALE
  real(8) :: HEIGHT_SCALE

  real(8) max_vort

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
      read(fid,*) varname, optimize_grid
      close(fid)
  end subroutine

  subroutine max_vort_fun(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id
      id = idx(i, j, offs, dims)
      max_vort = max(max_vort, abs(dom%vort%elts(id+1)))
  end subroutine

  subroutine write_and_print_step()
      real(4) timing
      real(8) tot_h
      integer l
      timing = get_timing()
      tot_h = integrate_hex(height_pert, level_start)
      max_vort = 0.0_8
      do l = level_start, level_end
        call apply_onescale(max_vort_fun, l, 0, 0)
      end do
      if (rank .eq. 0) write(*,'(A,F10.5,A,F10.5,A,I8,I8)') &
              'time [h] = ', time/3600.0_8, &
              ',  dt [s] = ', dt, &
              ',  active :', sum(n_active)
      if (rank .eq. 0) write(1011,'(E16.9, I3, 2(1X, I9), 2(1X, E24.15), 1X, F16.7)') &
              time, level_end, n_active, tot_h, max_vort, timing
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
  subroutine init_sol(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      real(8) r2, alph2
      d = dom%id+1
      id = idx(i, j, offs, dims)
      alph2 = (1.0_8 / 200000.0)**2
      if (dom%node%elts(id+1)%y .gt. 0) then
          r2 = dom%node%elts(id+1)%x**2 + dom%node%elts(id+1)%z**2 
          sol(S_HEIGHT)%data(d)%elts(id+1) = HEIGHT_SCALE + exp(-alph2*r2)
      else
          sol(S_HEIGHT)%data(d)%elts(id+1) = HEIGHT_SCALE
      end if
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0
  end subroutine

  subroutine cart2sph2(cin, cout)
      type(Coord) cin
      real(8), intent(out) :: cout(2)
      call cart2sph(cin, cout(1), cout(2))
  end subroutine

  subroutine wave_dump(fid)
      integer fid
  end subroutine

  subroutine wave_load(fid)
      integer fid
  end subroutine

  subroutine set_thresholds() ! inertia-gravity wave
      toll_height = VELO_SCALE*WAVE_SPEED * threshold**(3.0_8/2.0_8)
      toll_velo   = VELO_SCALE            * threshold**(3.0_8/2.0_8)
  end subroutine

end module

program wave
  use main_mod
  use wave_mod
  implicit none

  logical :: aligned
  integer d

  call init_main_mod()
  WAVE_SPEED = 1.0
  VELO_SCALE = GRAV_ACCEL/WAVE_SPEED
  HEIGHT_SCALE = 1000.0

  call read_test_case_parameters("wave.in")
  call set_thresholds()
  call initialize(apply_initial_conditions, 1, set_thresholds, wave_dump, wave_load)
  if (rank .eq. 0) write (*,*) 'thresholds p, u:',  toll_height, toll_velo
  if (rank .eq. 0) write (*,*) 'final time:', time_end

  if (rank .eq. 0) open(unit=1011,file='verlauf.out')

  if (threshold .eq. 0) level_start = level_end
  do while (time .lt. time_end - time_end*1.0e-13_8)
      call start_timing()
      call time_step(time_end, aligned)
      call stop_timing()
      call write_and_print_step()
      call print_load_balance()
  end do
  close(1011)
  call fill_up_grid_and_IWT()
  call trend_ml(sol, trend) ! to obtain vorticity on finest grid
  do d = 1, size(sol(S_HEIGHT)%data)
      sol(S_HEIGHT)%data(d)%elts = sol(S_HEIGHT)%data(d)%elts - HEIGHT_SCALE
  end do
  call export_2d(cart2sph2, (/sol(S_HEIGHT)/), 1, 100000+min_level*1000+level_end*10, level_end, &
                         (/-768, 768/), (/-384, 384/), (/2.0_8*MATH_PI, MATH_PI/), (/0.0_8, 0.0_8/))

  call finalize()
end program
