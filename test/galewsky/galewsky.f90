include 'gauss_quad.f90'

module galewsky_mod
  use main_mod
  use gauss_quad

  implicit none

  real(8), parameter :: VELO_SCALE = 80.0_8
  real(8) :: length_scale
  real(8) :: corolis_param
  real(8) :: Ro_number

  real(8), parameter :: GAL_H0   = 10157.9_8
  real(8), parameter :: GAL_ALPHA = 1.0_8/3.0_8
  real(8), parameter :: GAL_BETA  = 1.0_8/15.0_8
  real(8), parameter :: GAL_H_HAT = 120.0_8
  real(8) :: GAL_E_N  
  real(8) :: GAL_PHI0
  real(8) :: GAL_PHI1
  real(8) :: GAL_PHI2
  real(8) :: GAL_LON_PERP = 0.0_8

  type(Float_Field) pv_smooth

  logical perturbed
  integer iwrite

contains
  subroutine apply_initial_conditions()
      integer l
      do l = level_start, level_end
          call init_ring(MATH_PI*2.0_8/14.0_8, MATH_PI*5.0_8/14.0_8)      
          call apply_onescale(init_sol, l, 0, 1)
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
      read(fid,*) varname, perturbed
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

  subroutine postproc_pv(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, idW, idS, idSW
      id = idx(i, j, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      pv_smooth%data(dom%id+1)%elts(id+1) = 0.5 * ( &
              (dom%vort%elts(idW*TRIAG+LORT+1)+dom%vort%elts(id*TRIAG+UPLT+1)) &
              *dom%len%elts(id*EDGE+UP+1)*dom%pedlen%elts(id*EDGE+UP+1)+ &
              (dom%vort%elts(id*TRIAG+LORT+1)+dom%vort%elts(id*TRIAG+UPLT+1)) &
              *dom%len%elts(id*EDGE+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
              (dom%vort%elts(idS*TRIAG+UPLT+1)+dom%vort%elts(id*TRIAG+LORT+1)) &
              *dom%len%elts(id*EDGE+RT+1)*dom%pedlen%elts(id*EDGE+RT+1)+ &
              (dom%vort%elts(idS*TRIAG+UPLT+1)+dom%vort%elts(idSW*TRIAG+LORT+1)) &
              *dom%len%elts(idS*EDGE+UP+1)*dom%pedlen%elts(idS*EDGE+UP+1)+ &
              (dom%vort%elts(idSW*TRIAG+UPLT+1)+dom%vort%elts(idSW*TRIAG+LORT+1)) &
              *dom%len%elts(idSW*EDGE+DG+1)*dom%pedlen%elts(idSW*EDGE+DG+1)+ &
              (dom%vort%elts(idSW*TRIAG+UPLT+1)+dom%vort%elts(idW*TRIAG+LORT+1)) &
              *dom%len%elts(idW*EDGE+RT+1)*dom%pedlen%elts(idW*EDGE+RT+1))/ &
              (dom%len%elts(id*EDGE+UP+1)*dom%pedlen%elts(id*EDGE+UP+1)+ &
              dom%len%elts(id*EDGE+DG+1)*dom%pedlen%elts(id*EDGE+DG+1)+ &
              dom%len%elts(id*EDGE+RT+1)*dom%pedlen%elts(id*EDGE+RT+1)+ &
              dom%len%elts(idS*EDGE+UP+1)*dom%pedlen%elts(idS*EDGE+UP+1)+ &
              dom%len%elts(idSW*EDGE+DG+1)*dom%pedlen%elts(idSW*EDGE+DG+1)+ &
              dom%len%elts(idW*EDGE+RT+1)*dom%pedlen%elts(idW*EDGE+RT+1))
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
    height_fun = GAL_H0 - integral/GRAV_ACCEL
    if (perturbed) height_fun = height_fun &
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
  subroutine cart2sph2(cin, cout)
      type(Coord) cin
      real(8), intent(out) :: cout(2)
      call cart2sph(cin, cout(1), cout(2))
  end subroutine
  subroutine galewsky_dump(fid)
      integer fid
      write(fid) iwrite
  end subroutine
  subroutine galewsky_load(fid)
      integer fid
      read(fid) iwrite
  end subroutine
  subroutine set_thresholds()
      toll_velo  = Ro_number*velo_scale*threshold**(3.0_8/2.0_8)
      toll_height = Ro_number*velo_scale*length_scale*corolis_param*threshold**(3.0_8/2.0_8)
  end subroutine

end module
program galewsky
  use main_mod
  use galewsky_mod

  implicit none

  logical aligned
  real(8) last_checkpoint

  integer recommanded_level_start, d
  character(8+8+29) command

  call init_main_mod()

  call read_test_case_parameters("galewsky.in")
  length_scale = MATH_PI*radius*3.0_8/14.0_8
  corolis_param = 2.0_8*omega*sqrt(2.0_8)/2.0_8
  Ro_number = VELO_SCALE/(length_scale*corolis_param)

  iwrite = 0
  call set_thresholds()
  if (rank .eq. 0) write (*,*) 'thresholds p, u:',  toll_height, toll_velo
  call initialize(apply_initial_conditions, 1, set_thresholds, galewsky_dump, galewsky_load)

  recommanded_level_start = write_active_per_level()

  if (rank .eq. 0) open(unit=1011,file='verlauf.out')
  call write_and_print_step()
  do while (time .lt. time_end)
      call start_timing()
      call time_step(dt_write, aligned)
      call stop_timing()
      call write_and_print_step()
      call print_load_balance()
!     if (istep .eq. 20) then
!         close(1011) ! close and reopen to force write
!         open(unit=1011,file='verlauf.out',access='APPEND')
!     end if
!     if (time-last_checkpoint .gt. 3600.0_8) then
!         if (rank .eq. 0) write(*,'(A, I4, A, L1, A)') &
!                 'clean and checkpoint', cp_idx+1, ' (exceed = ', max_level_exceeded, ')'      
!         close(1011) ! close and reopen to force write
!         call checkpoint(galewsky_dump, galewsky_load)
!         last_checkpoint = time
!         open(unit=1011,file='verlauf.out',access='APPEND')
!     end if
  end do
  close(1011)
  call fill_up_grid_and_IWT()
  call trend_ml(sol, trend) ! to obtain vorticity on finest grid

  call init_Float_Field(pv_smooth, S_HEIGHT)
  do d = 1, n_domain(rank+1)
      call init(pv_smooth%data(d), grid(d)%node%length)
  end do
  call apply_onescale(postproc_pv, level_end, 0, 0)
  call update_bdry(pv_smooth, level_end)

  call export_2d(cart2sph2, (/sol(S_HEIGHT), pv_smooth/), 2, 100000+min_level*1000+level_end*10, level_end, &
                         (/-768, 768/), (/-384, 384/), (/2.0_8*MATH_PI, MATH_PI/), (/0.0_8, 0.0_8/))

  call finalize()
end program

