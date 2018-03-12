module DCMIP2012c4_mod
  ! DCMIP2012 test case 4: baroclinic instability
  use main_mod
  use remap_mod
  implicit none

  integer                            :: CP_EVERY, iwrite, N_node
  integer, dimension(:), allocatable :: n_patch_old, n_node_old
  real(8)                            :: initotalmass, totalmass, timing, total_cpu_time
  logical                            :: wasprinted, uniform
  character (255)                    :: IC_file

  real(8)                            :: c_v, d2, h_0, lat_c, lon_c, N_freq, T_0
  real(8)                            :: acceldim, f0, geopotdim, Ldim, Hdim, massdim, Tdim, Tempdim, Udim, pdim, R_ddim, specvoldim
  real(8)                            :: norm_mass, norm_temp, norm_velo, norm_mass_trend, norm_temp_trend, norm_velo_trend
  real(8)                            :: mass_scale, temp_scale, velo_scale, mass_scale_trend, temp_scale_trend, velo_scale_trend
  real(8)                            :: l2_mass, l2_temp, l2_velo, mass_error
  real(8)                            :: visc, viscosity_divu, viscosity_rotu, viscosity_mass, viscosity_temp, ray_friction
  real(8)                            :: delta_T, eta, eta_t, eta_v, eta_0, gamma_T, R_pert, u_p, u_0

  type(Float_Field)                  :: rel_vort 
contains
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    type (Domain)                  :: dom
    integer                        :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: rgrc, lev_press, pot_temp, p_top, p_bot

    d = dom%id+1

    id   = idx(i, j, offs, dims)
    idN  = idx(i, j + 1, offs, dims)
    idE  = idx(i + 1, j, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    ! Surface pressure
    dom%surf_press%elts(id+1) = surf_pressure_fun (x_i)

    ! Pressure at level zlev
    lev_press = 0.5_8*(a_vert(zlev)+a_vert(zlev+1))*ref_press + 0.5_8*(b_vert(zlev)+b_vert(zlev+1))*dom%surf_press%elts(id+1)

    ! Normalized pressure
    eta = lev_press/dom%surf_press%elts(id+1)
    eta_v = (eta - eta_0) * MATH_PI/2.0_8

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = &
         ((a_vert(zlev)-a_vert(zlev+1))*ref_press + (b_vert(zlev)-b_vert(zlev+1))*dom%surf_press%elts(id+1))/grav_accel

    ! Horizontally uniform potential temperature
    pot_temp =  set_temp(x_i) * (lev_press/ref_press)**(-kappa)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  function set_temp(x_i)
    real(8) :: set_temp
    type(Coord) :: x_i
    real(8) :: lon, lat, Tmean

    call cart2sph (x_i, lon, lat)
    
    if (eta.ge.eta_t) then
       Tmean = T_0*eta**(R_d*Gamma_T/grav_accel)
    else
       Tmean = T_0*eta**(R_d*Gamma_T/grav_accel) + delta_T * (eta_t - eta)**5
    end if

    set_temp = Tmean + 0.75_8 * eta*MATH_PI*u_0/R_d * sin(eta_v) * sqrt(cos(eta_v)) * &
         (2.0_8*cos(eta_v)**1.5*(-2.0_8*sin(lat)**6*(cos(lat)**2+1.0_8/3.0_8) + 10.0_8/63.0_8) + &
         radius*omega*(8.0_8/5.0_8*cos(lat)**3*(sin(lat)**2+2.0_8/3.0_8) - MATH_PI/4.0_8))
  end function set_temp
 
  subroutine set_surfgeopot (dom, i, j, zlev, offs, dims)
    ! Initialize surface geopotential after restart
    type (Domain)                   :: dom
    integer                         :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type (Coord) :: x_i
    integer      :: id

    id   = idx(i, j, offs, dims)
    x_i  = dom%node%elts(id+1)

    ! Surfaced geopotential
    dom%surf_geopot%elts(id+1) = surf_geopot_fun(x_i)
  end subroutine set_surfgeopot
 
  function geopot_fun (x_i)
    ! Geopotential 
    type(Coord) :: x_i
    real(8)     :: geopot_fun
    
    real(8) :: lon, lat, phi_mean, rgrc

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    if (eta.ge.eta_t) then
       phi_mean = T_0*grav_accel/gamma_T * (1.0_8 - eta**(R_d*gamma_T/grav_accel))
    else
       phi_mean = T_0*grav_accel/gamma_T * (1.0_8 - eta**(R_d*gamma_T/grav_accel)) - delta_phi(eta)
    end if

    geopot_fun = phi_mean + u_0*cos(eta_v)**1.5*(u_0*cos(eta_v)**1.5* &
         (-2.0_8*sin(lat)**6*(cos(lat)**2 + 1.0_8/3.0_8) + 10.0_8/63.0_8) + &
         radius*omega*(8.0_8/5.0_8*cos(lat)**3*(sin(lat)**2 + 2.0_8/3.0_8) - MATH_PI/4.0_8))
  end function geopot_fun

  function delta_phi (eta)
    real(8) :: delta_phi
    real(8) :: eta

    delta_phi = R_d * delta_T*((log(eta/eta_t) + 137.0_8/60.0_8)*eta_t**5 - 5.0_8*eta_t**4*eta + 5.0_8*eta_t**3*eta**2 &
         - 10.0_8/3.0_8 * eta_t**2*eta**3 + 5.0_8/4.0_8*eta_t*eta**4 - eta**5/5.0_8)
  end function delta_phi
  
  function surf_geopot_fun (x_i)
    ! Surface geopotential
    Type(Coord) :: x_i
    real(8)     :: surf_geopot_fun
    real(8)     :: lon, lat

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    surf_geopot_fun = u_0*cos((1.0_8-eta_0)*MATH_PI/2.0_8)**1.5 * &
         (u_0*cos((1.0_8-eta_0)*MATH_PI/2.0_8)**1.5 * &
         (-2.0_8*sin(lat)**6*(cos(lat)**2 + 1.0_8/3.0_8) + 10.0_8/63.0_8)  + &
         radius*omega*(8.0_8/5.0_8*cos(lat)**3*(sin(lat)**2 + 2.0_8/3.0_8) - MATH_PI/4.0_8))
  end function surf_geopot_fun

   function surf_pressure_fun (x_i)
    ! Surface pressure
    type(Coord) :: x_i
    real(8)     :: surf_pressure_fun

    surf_pressure_fun = ref_press
  end function surf_pressure_fun

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    real(8) :: lon, lat, u, v
    real(8) :: rgrc

    ! Great circle distance
    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))
    
    u = u_0*cos(eta_v)**1.5*sin(2.0_8*lat)**2 + u_p*exp__flush(-(rgrc/R_pert)**2)  ! Zonal velocity component
    v = 0.0_8         ! Meridional velocity component
  end subroutine vel_fun

  subroutine initialize_a_b_vert
    integer :: k

    ! Allocate vertical grid parameters
    if (allocated(a_vert)) deallocate(a_vert)
    if (allocated(b_vert)) deallocate(b_vert)
    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))

    if (uniform) then
       do k = 1, zlevels+1
          press_infty = 0.0_8
          a_vert(k) = real(k-1)/real(zlevels) * press_infty/ref_press
          b_vert(k) = 1.0_8 - real(k-1)/real(zlevels)
       end do
    else
       if (zlevels.eq.30) then
          a_vert = (/ 0.00225523952394724, 0.00503169186413288, 0.0101579474285245, 0.0185553170740604, 0.0306691229343414, &
               0.0458674766123295, 0.0633234828710556, 0.0807014182209969, 0.0949410423636436, 0.11169321089983, & 
               0.131401270627975, 0.154586806893349, 0.181863352656364, 0.17459799349308, 0.166050657629967, &
               0.155995160341263, 0.14416541159153, 0.130248308181763, 0.113875567913055, 0.0946138575673103, &
               0.0753444507718086, 0.0576589405536652, 0.0427346378564835, 0.0316426791250706, 0.0252212174236774, &
               0.0191967375576496, 0.0136180268600583, 0.00853108894079924, 0.00397881818935275, 0.0, 0.0 /)
          b_vert = (/ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0393548272550106, &
               0.0856537595391273, 0.140122056007385, 0.204201176762581, 0.279586911201477, 0.368274360895157,  &
               0.47261056303978, 0.576988518238068, 0.672786951065063, 0.753628432750702, 0.813710987567902, &
               0.848494648933411, 0.881127893924713, 0.911346435546875, 0.938901245594025, 0.963559806346893, &
               0.985112190246582, 1.0 /)
       else
          write(0,*) "For this number of zlevels, no rule has been defined for a_vert and b_vert"
          stop
       end if

       ! DCMIP order is opposite to ours
       if (.not. uniform) then
          a_vert = a_vert(zlevels+1:1:-1)
          b_vert = b_vert(zlevels+1:1:-1)
       end if

       ! Set pressure at infinity
       press_infty = a_vert(zlevels+1)*ref_press ! note that b_vert at top level is 0, a_vert is small but non-zero
    end if
  end subroutine initialize_a_b_vert

  subroutine read_test_case_parameters (filename)
    character(*)   :: filename
    integer        :: fid = 500
    character(255) :: varname
    
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
       write(*,'(A,es10.4)') "threshold        = ", threshold
       write(*,'(A,i2)')     "optimize_grid    = ", optimize_grid 
       write(*,'(A,es10.4)') "dt_write         = ", dt_write
       write(*,'(A,i6)')     "CP_EVERY         = ", CP_EVERY
       write(*,'(A,es10.4)') "time_end         = ", time_end 
       write(*,'(A,i6)')     "resume           = ", resume
       write(*,*) ' '
    end if
    dt_write = dt_write * 60.0_8
    time_end = time_end * 60.0_8**2

    close(fid)
  end subroutine read_test_case_parameters

  subroutine write_and_export (iwrite, zlev)
    integer :: iwrite, zlev

    integer :: d, i, j, k, l, p, u

    if (rank.eq.0) write(6,*) 'Saving fields'

    call update_array_bdry (sol, NONE)

    call pre_levelout

    ! First integrate pressure down across all grid points in order to compute surface pressure
    do k = zlevels, 1, -1
       do d = 1, size(grid)
          mass => sol(S_MASS,k)%data(d)%elts
          temp => sol(S_TEMP,k)%data(d)%elts

          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (integrate_pressure_down, grid(d), p-1, k, 0, 1)
          end do

          nullify (mass, temp)
       end do
    end do

    do l = level_start, level_end
       minv = 1.0d63; maxv = -1.0d63
       u = 100000+100*iwrite

       ! Calculate pressure, exner and geopotential at vertical level zlev and scale l
       do k = 1, zlev
          do d = 1, size(grid)
             mass  => sol(S_MASS,k)%data(d)%elts
             temp  => sol(S_TEMP,k)%data(d)%elts
             exner => exner_fun(k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
             nullify (mass, temp, exner)
          end do
       end do

       ! Calculate zonal and meridional velocities and vorticity for vertical level zlev
       do d = 1, size(grid)
          velo => sol(S_VELO,zlev)%data(d)%elts
          vort => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), zlev,    0, 0)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, zlev)
          nullify (velo, vort)
       end do

       ! Calculate vorticity at hexagon points (stored in adj_mass)
       call apply_onescale (vort_triag_to_hex, l, z_null, 0, 1)

       Call write_level_mpi (write_primal, u+l, l, zlev, .True.)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
       if (rank .eq. 0) write(u,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
       if (rank .eq. 0) write(u,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
       u = 200000+100*iwrite
       call write_level_mpi (write_dual, u+l, l, zlev, .False.)
    end do

    call post_levelout
    call barrier
    if (rank .eq. 0) call compress_files (iwrite)

    ! Save 2D projection
    call export_2d (cart2sph2, 300000+100*iwrite, (/-96, 96/), (/-48, 48/), (/2.0_8*MATH_PI, MATH_PI/), set_thresholds)
    !call export_2d (cart2sph2, 300000+100*iwrite, (/-768, 768/), (/-384, 384/), (/2.0_8*MATH_PI, MATH_PI/), set_thresholds)
  end subroutine write_and_export

  subroutine DCMIP2012c4_dump (fid)
    integer :: fid
    
    write(fid) itime
    write(fid) iwrite
    write(fid) tol_mass, tol_temp, tol_velo
  end subroutine DCMIP2012c4_dump

  subroutine DCMIP2012c4_load (fid)
    integer :: fid
    
    read(fid) itime
    read(fid) iwrite
    read(fid) tol_mass, tol_temp, tol_velo
  end subroutine DCMIP2012c4_load

  subroutine set_thresholds (itype)
    integer, optional :: itype

    integer :: l, k

    ! Set thresholds dynamically (trend or sol must be known)
    if (itype.eq.0) then ! Adapt on trend
       norm_mass_trend = 0.0_8
       norm_temp_trend = 0.0_8
       norm_velo_trend = 0.0_8
       do l = level_start, level_end
          call apply_onescale (linf_trend, l, z_null, 0, 0)
       end do

       mass_scale = sync_max_d (norm_mass_trend)
       temp_scale = sync_max_d (norm_temp_trend)
       velo_scale = sync_max_d (norm_velo_trend)
    else ! Adapt on variables
       norm_mass = 0.0_8
       norm_temp = 0.0_8
       norm_velo = 0.0_8
       do l = level_start, level_end
          call apply_onescale (linf_vars, l, z_null, 0, 0)
       end do

       mass_scale = sync_max_d (norm_mass)
       temp_scale = sync_max_d (norm_temp)
       velo_scale = sync_max_d (norm_velo)
    end if

    if (istep.ne.0) then
       tol_mass = 0.99_8*tol_mass + 0.01_8*threshold * mass_scale
       tol_temp = 0.99_8*tol_temp + 0.01_8*threshold * temp_scale
       tol_velo = 0.99_8*tol_velo + 0.01_8*threshold * velo_scale
    elseif (istep.eq.0) then
       tol_mass = threshold * mass_scale
       tol_temp = threshold * temp_scale
       tol_velo = threshold * velo_scale
       if (adapt_trend .and. itype.eq.1) then ! Re-scale trend threshold for variables
          tol_mass = threshold**1.5_8 * mass_scale/5.0d1
          tol_temp = threshold**1.5_8 * temp_scale/5.0d1
          tol_velo = threshold**1.5_8 * velo_scale/5.0d1
       end if
    end if
  end subroutine set_thresholds

  subroutine linf_trend (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, e, k

    id = idx(i, j, offs, dims)

    ! Maximum trends
    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       do k = 1, zlevels
          norm_mass_trend = max(norm_mass_trend, abs(trend(S_MASS,k)%data(dom%id+1)%elts(id+1)))
          norm_temp_trend = max(norm_temp_trend, abs(trend(S_TEMP,k)%data(dom%id+1)%elts(id+1)))
          do e = 1, EDGE
             norm_velo_trend  = max(norm_velo_trend, abs(trend(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e)))
          end do
       end do
    end if
  end subroutine linf_trend

  subroutine linf_vars (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, e, k

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       do k = 1, zlevels
          norm_mass = max(norm_mass, abs(sol(S_MASS,k)%data(d)%elts(id+1)))
          norm_temp = max(norm_temp, abs(sol(S_TEMP,k)%data(d)%elts(id+1)))
          do e = 1, EDGE
             norm_velo  = max(norm_velo, abs(sol(S_VELO,k)%data(d)%elts(EDGE*id+e)))
          end do
       end do
    end if
  end subroutine linf_vars

  subroutine l2_trend (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, k

    id = idx(i, j, offs, dims)
    d = dom%id+1

    ! L2 norms of trends
    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       N_node = N_node + 1
       do k = 1, zlevels
          norm_mass_trend = norm_mass_trend + trend(S_MASS,k)%data(d)%elts(id+1)**2
          norm_temp_trend = norm_temp_trend + trend(S_TEMP,k)%data(d)%elts(id+1)**2
          do e = 1, EDGE
             norm_velo_trend  = norm_velo_trend + trend(S_VELO,k)%data(d)%elts(EDGE*id+e)**2
          end do
       end do
    endif
  end subroutine l2_trend

  subroutine l2_vars (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, k

    d = dom%id+1
    id = idx(i, j, offs, dims)

    ! L2 norms of trends
    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       N_node = N_node + 1
       do k = 1, zlevels
          norm_mass = norm_mass + sol(S_MASS,zlev)%data(d)%elts(id+1)**2
          norm_temp = norm_temp + sol(S_TEMP,zlev)%data(d)%elts(id+1)**2
          do e = 1, EDGE
             norm_velo  = norm_velo + sol(S_VELO,zlev)%data(d)%elts(EDGE*id+e)**2
          end do
       end do
    endif
  end subroutine l2_vars
  subroutine apply_initial_conditions
    integer :: k, l

    wasprinted=.false.
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -1, 1)
          wasprinted=.false.
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine set_surf_geopot
    integer ::  d, p

    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (set_surfgeopot, grid(d), p-1, z_null, 0, 1)
       end do
    end do
  end subroutine set_surf_geopot

  subroutine sum_total_mass (initialgo)
    ! Total mass over all vertical layers
    logical :: initialgo

    integer :: k

    if (initialgo) then
       initotalmass = 0.0_8
       do k = 1, zlevels
          initotalmass = initotalmass + integrate_hex (mu, level_start, k)
       end do
    else
       totalmass = 0.0_8
       do k = 1, zlevels
          totalmass = totalmass + integrate_hex (mu, level_start, k)
       end do
       mass_error = abs(totalmass-initotalmass)/initotalmass
    end if
  end subroutine sum_total_mass
end module DCMIP2012c4_mod

program DCMIP2012c4
  use main_mod
  use DCMIP2012c4_mod
  implicit none

  integer                      :: d, ierr, k, l, v, zlev
  integer, parameter           :: len_cmd_files = 12 + 4 + 12 + 4
  integer, parameter           :: len_cmd_archive = 11 + 4 + 4
  character(len_cmd_files)     :: cmd_files
  character(len_cmd_archive)   :: cmd_archive
  character(8+8+29+14)         :: command
  character(9+len_cmd_archive) :: command1
  character(6+len_cmd_files)   :: command2
  logical                      :: aligned, remap, write_init

  ! Initialize grid etc
  call init_main_mod 

  ! Nullify all pointers initially
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, qe, vort)

  ! Read test case parameters
  call read_test_case_parameters ("DCMIP2012c4.in")

  ! Average minimum grid size and maximum wavenumber
  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8))

  ! Parameters for the simulation
  radius         = 6.371229d6    ! mean radius of the Earth in meters
  grav_accel     = 9.80616_8     ! gravitational acceleration in meters per second squared
  omega          = 7.29212d-5    ! Earth’s angular velocity in radians per second
  f0             = 2.0_8*omega   ! Coriolis parameter
  u_0            = 35.0_8        ! maximum velocity of zonal wind
  u_p            = 1.0_8         ! maximum perturbation to zonal wind
  R_pert         = radius/10.0_8 ! 
  T_0            = 288.0_8       ! temperature in Kelvin
  gamma_T        = 0.005         ! temperature lapse rate
  delta_T        = 4.8d5         ! empirical temperature difference
  eta_0          = 0.252_8       ! value of eta at reference level (level of the jet)
  eta_t          = 0.2_8         ! value of eta at the tropopause
  lon_c          = MATH_PI/9.0_8 ! longitude location of perturbation to zonal wind
  lat_c          = 2.0_8*MATH_PI/9.0_8 ! latitude location of perturbation to zonal wind
  ref_press      = 1.0d5        ! reference pressure (mean surface pressure) in Pascals
  ref_surf_press = ref_press    ! reference surface pressure
  R_d            = 287.0_8      ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_p            = 1004.64_8     ! specific heat at constant pressure in joules per kilogram Kelvin
  c_v            = 717.6_8       ! specfic heat at constant volume c_v = R_d - c_p
  gamma          = c_p/c_v       ! heat capacity ratio
  kappa          = 2.0_8/7.0_8   ! kappa=R_d/c_p
  N_freq         = sqrt(grav_accel**2/(c_p*T_0)) ! Brunt-Vaisala buoyancy frequency

  ! Dimensional scaling
  Ldim           = sqrt(d2)                         ! horizontal length scale
  Hdim           = h_0                              ! vertical length scale
  Udim           = u_0                              ! velocity scale
  Tdim           = Hdim/Udim                        ! time scale
  Tempdim        = T_0                              ! temperature scale (both theta and T from DYNAMICO)

  acceldim       = Udim**2/Hdim                     ! acceleration scale
  pdim           = ref_press                        ! pressure scale
  R_ddim         = R_d                              ! R_d scale
  massdim        = pdim*Hdim/(Tempdim*R_d)          ! mass (=rho*dz following DYNAMICO) scale
  specvoldim     = (R_d*Tempdim)/pdim               ! specific volume scale
  geopotdim      = acceldim*massdim*specvoldim/Hdim ! geopotential scale
  wave_speed     = sqrt(gamma*pdim*specvoldim)      ! acoustic wave speed
  cfl_num        = 0.8_8                            ! cfl number
  n_remap        = 10                               ! Vertical remap interval

  ray_friction   = 0.0_8                            ! Rayleigh friction

  zlev           = 8 ! about 820 hPa
  save_levels    = 1; allocate(pressure_save(1:save_levels))  ! number of vertical levels to save
  level_save     = level_end                                  ! resolution level at which to save lat-lon data
  pressure_save  = (/850.0d2/)                                ! interpolate values to this pressure level when interpolating to lat-lon grid

  ! Set logical switches
  adapt_trend  = .false. ! Adapt on trend or on variables
  adapt_dt     = .true.  ! Adapt time step
  compressible = .true.  ! Compressible equations
  remap        = .false. ! Remap vertical coordinates (always remap when saving results)
  uniform      = .false. ! Type of vertical grid

  ! Set viscosity
  visc = 2.0d-4 ! Constant for viscosity

  viscosity_mass = visc * dx_min**2 ! viscosity for mass equation
  viscosity_temp = visc * dx_min**2 ! viscosity for mass-weighted potential temperature equation
  viscosity_divu = visc * dx_min**2 ! viscosity for divergent part of momentum equation
  viscosity_rotu = visc/1.0d2 * dx_min**2 ! viscosity for rotational part of momentum equation
  viscosity = max (viscosity_mass, viscosity_temp, viscosity_divu, viscosity_rotu)

  ! Time step based on acoustic wave speed and hexagon edge length (not used if adaptive dt)  
  dt_init = min(cfl_num*dx_min/wave_speed, 0.25_8*dx_min**2/viscosity)  
  if (rank.eq.0) write(6,'(2(A,es10.4,1x))') "dt_cfl = ", cfl_num*dx_min/(wave_speed+u_0), " dt_visc = ", 0.25_8*dx_min**2/viscosity

  if (rank .eq. 0) then
     write(6,'(A,es10.4)') 'Viscosity_mass   = ', viscosity_mass
     write(6,'(A,es10.4)') 'Viscosity_temp   = ', viscosity_temp
     write(6,'(A,es10.4)') 'Viscosity_divu   = ', viscosity_divu
     write(6,'(A,es10.4)') 'Viscosity_rotu   = ', viscosity_rotu
     write(6,'(A,es10.4)') ' '
  end if

  ! Initialize vertical grid
  call initialize_a_b_vert

  ! Initialize variables
  call initialize (apply_initial_conditions, 1, set_thresholds, DCMIP2012c4_dump, DCMIP2012c4_load)

  allocate (n_patch_old(size(grid)), n_node_old(size(grid)))
  n_patch_old = 2;  call set_surf_geopot 

  call sum_total_mass (.True.)

  if (rank .eq. 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:', tol_mass, tol_temp, tol_velo
  call barrier

  if (rank .eq. 0) write(6,*) 'Write initial values and grid'
  call write_and_export (iwrite, zlev)

  if (resume.le.0) iwrite = 0
  total_cpu_time = 0.0_8

  open(unit=12, file='DCMIP2012c4_log', action='WRITE', form='FORMATTED')
  if (rank .eq. 0) then
     write (6,'(A,ES12.6,3(A,ES10.4),A,I2,A,I9)') &
          ' time [h] = ', time/3600.0_8, &
          '  mass tol = ', tol_mass, &
          ' temp tol = ', tol_temp, &
          ' velo tol = ', tol_velo, &
          ' Jmax =', level_end, &
          '  dof = ', sum(n_active)
  end if

  do while (time .lt. time_end)
     call update_array_bdry (sol, NONE)
     n_patch_old = grid(:)%patch%length
     n_node_old = grid(:)%node%length

     if (remap .and. mod(istep, n_remap).eq.0 .and. istep.gt.1) call remap_vertical_coordinates (set_thresholds)

     call start_timing
     call time_step (dt_write, aligned, set_thresholds)
     call stop_timing

     call set_surf_geopot
     timing = get_timing()
     total_cpu_time = total_cpu_time + timing

     if (rank .eq. 0) then
        write (6,'(A,ES12.6,4(A,ES10.4),A,I2,A,I9,A,ES8.2,1x,A,ES8.2)') &
             ' time [h] = ', time/60.0_8**2, &
             ' dt [s] = ', dt, &
             '  mass tol = ', tol_mass, &
             ' temp tol = ', tol_temp, &
             ' velo tol = ', tol_velo, &
             ' Jmax = ', level_end, &
             '  dof = ', sum(n_active), &
             ' mass error = ', mass_error, &
             ' cpu = ', timing

        write (12,'(5(ES15.9,1x),I2,1X,I9,1X,2(ES15.9,1x))')  &
             time/3600.0_8, dt, tol_mass, tol_temp, tol_velo, level_end, sum(n_active), mass_error,  timing
     end if

     if (aligned) then
        iwrite = iwrite + 1

        ! Save fields
        if (remap) call remap_vertical_coordinates (set_thresholds)
        call write_and_export (iwrite, zlev)

        call sum_total_mass (.False.)

        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle ! Do not write checkpoint

        ierr = write_checkpoint (DCMIP2012c4_dump)

        ! Let all cpus exit gracefully if NaN has been produced
        ierr = sync_max (ierr)
        if (ierr .eq. 1) then ! NaN
           write (0,*) "NaN when writing checkpoint"
           call finalize
           stop
        end if

        ! Restart after checkpoint and load balance
        call restart_full (set_thresholds, DCMIP2012c4_load)
        call print_load_balance

        call barrier
     end if
     call sum_total_mass (.False.)
  end do

  if (rank .eq. 0) then
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
     close (12)
     close (1011)
     close (8450)
     command = '\rm tmp tmp1 tmp2'; call system (command)
  end if

  call finalize
end program DCMIP2012c4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (dom, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  use DCMIP2012c4_mod
  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: grad
  logical :: local_type

  if (present(type)) then
     local_type = type
  else
     local_type = .false.
  end if

  if (max(viscosity_mass, viscosity_temp).eq.0.0_8) then
     physics_scalar_flux = 0.0_8
  else
     ! Calculate gradients
     if (.not.local_type) then ! Usual gradient at edges of hexagon E, NE, N
        grad(S_MASS,RT+1) = (mass(idE+1) - mass(id+1))  /dom%len%elts(EDGE*id+RT+1) 
        grad(S_MASS,DG+1) = (mass(id+1)  - mass(idNE+1))/dom%len%elts(EDGE*id+DG+1) 
        grad(S_MASS,UP+1) = (mass(idN+1) - mass(id+1))  /dom%len%elts(EDGE*id+UP+1) 

        grad(S_TEMP,RT+1) = (temp(idE+1) - temp(id+1))  /dom%len%elts(EDGE*id+RT+1) 
        grad(S_TEMP,DG+1) = (temp(id+1)  - temp(idNE+1))/dom%len%elts(EDGE*id+DG+1) 
        grad(S_TEMP,UP+1) = (temp(idN+1) - temp(id+1))  /dom%len%elts(EDGE*id+UP+1) 
     else ! Gradient for southwest edges of hexagon W, SW, S
        grad(S_MASS,RT+1) = -(mass(idE+1) - mass(id+1))  /dom%len%elts(EDGE*idE+RT+1) 
        grad(S_MASS,DG+1) = -(mass(id+1)  - mass(idNE+1))/dom%len%elts(EDGE*idNE+DG+1)
        grad(S_MASS,UP+1) = -(mass(idN+1) - mass(id+1))  /dom%len%elts(EDGE*idN+UP+1) 

        grad(S_TEMP,RT+1) = -(temp(idE+1) - temp(id+1))  /dom%len%elts(EDGE*idE+RT+1) 
        grad(S_TEMP,DG+1) = -(temp(id+1)  - temp(idNE+1))/dom%len%elts(EDGE*idNE+DG+1)
        grad(S_TEMP,UP+1) = -(temp(idN+1) - temp(id+1))  /dom%len%elts(EDGE*idN+UP+1) 
     end if

     ! Fluxes of physics
     if (.not.local_type) then ! Usual flux at edges E, NE, N
        physics_scalar_flux(S_MASS,RT+1) = -viscosity_mass * grad(S_MASS,RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
        physics_scalar_flux(S_MASS,DG+1) = -viscosity_mass * grad(S_MASS,DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
        physics_scalar_flux(S_MASS,UP+1) = -viscosity_mass * grad(S_MASS,UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

        physics_scalar_flux(S_TEMP,RT+1) = -viscosity_temp * grad(S_TEMP,RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
        physics_scalar_flux(S_TEMP,DG+1) = -viscosity_temp * grad(S_TEMP,DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
        physics_scalar_flux(S_TEMP,UP+1) = -viscosity_temp * grad(S_TEMP,UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
     else ! Flux at edges W, SW, S
        physics_scalar_flux(S_MASS,RT+1) = -viscosity_mass * grad(S_MASS,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(S_MASS,DG+1) = -viscosity_mass * grad(S_MASS,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(S_MASS,UP+1) = -viscosity_mass * grad(S_MASS,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)

        physics_scalar_flux(S_TEMP,RT+1) = -viscosity_temp * grad(S_TEMP,RT+1) * dom%pedlen%elts(EDGE*idE+RT+1)
        physics_scalar_flux(S_TEMP,DG+1) = -viscosity_temp * grad(S_TEMP,DG+1) * dom%pedlen%elts(EDGE*idNE+DG+1)
        physics_scalar_flux(S_TEMP,UP+1) = -viscosity_temp * grad(S_TEMP,UP+1) * dom%pedlen%elts(EDGE*idN+UP+1)
     end if
  end if
end function physics_scalar_flux

function physics_scalar_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the scalar trend
  ! In this test case there is no scalar source term
  use domain_mod
  use DCMIP2012c4_mod
  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source = 0.0_8
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  !
  ! In this test case we add Rayleigh friction and Laplacian diffusion
  use domain_mod
  use ops_mod
  use DCMIP2012c4_mod
  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(Domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                    :: e, id
  real(8), dimension(1:EDGE) :: diffusion, friction, curl_rotu, grad_divu

  interface
     function velo_diffusion (dom, i, j, zlev, offs, dims)
       use domain_mod
       real(8), dimension(1:EDGE)     :: velo_diffusion
       type(Domain)                   :: dom
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function velo_diffusion
  end interface

  id = idx(i, j, offs, dims)
  
  if (max(viscosity_divu, viscosity_rotu).eq.0.0_8) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     grad_divu = gradi_e (divu, dom, i, j, offs, dims)
     curl_rotu = curlv_e (vort, dom, i, j, offs, dims)
     do e = 1, EDGE 
        diffusion(e) = viscosity_divu * grad_divu(e) - viscosity_rotu * curl_rotu(e)
     end do
  end if

  ! Calculate Rayleigh friction
  do e = 1, EDGE 
     friction(e) = -ray_friction*velo(EDGE*id+e)
  end do

  ! Total physics for source term of velocity trend
  do e = 1, EDGE
     physics_velo_source(e) = friction(e) + diffusion(e)
  end do
end function physics_velo_source


