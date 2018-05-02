module DCMIP2008c5_mod
  ! DCMIP2008c5 test case 5 parameters
  use main_mod
  use remap_mod
  implicit none

  character(*), parameter            :: test_case = "DCMIP2008c5"               
  integer                            :: CP_EVERY, iwrite, N_node
  integer, dimension(:), allocatable :: n_patch_old, n_node_old
  real(8)                            :: initotalmass, totalmass, timing, total_cpu_time
  logical                            :: wasprinted, uniform
  character (255)                    :: IC_file

  real(8)                            :: c_v, d2, h_0, lat_c, lon_c, N_freq, T_0, u_0
  real(8)                            :: acceldim, f0, geopotdim, Ldim, Hdim, massdim, Tdim, Tempdim, Udim, pdim, R_ddim, specvoldim
  real(8)                            :: norm_mass, norm_temp, norm_velo, norm_mass_trend, norm_temp_trend, norm_velo_trend
  real(8)                            :: mass_scale, temp_scale, velo_scale, mass_scale_trend, temp_scale_trend, velo_scale_trend
  real(8)                            :: l2_mass, l2_temp, l2_velo, mass_error
  real(8)                            :: visc, ray_friction

  type(Float_Field)                  :: rel_vort 
contains
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: column_mass, lev_press, pot_temp, p_top, p_bot

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
    column_mass = dom%surf_press%elts(id+1)/grav_accel

    ! Pressure at level zlev
    lev_press = 0.5_8*(a_vert(zlev)+a_vert(zlev+1))*ref_press + 0.5_8*(b_vert(zlev)+b_vert(zlev+1))*dom%surf_press%elts(id+1)

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*column_mass

    ! Horizontally uniform potential temperature
    pot_temp =  T_0 * (lev_press/ref_press)**(-kappa)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  subroutine set_surfgeopot (dom, i, j, zlev, offs, dims)
    ! Initialize surface geopotential after restart
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type (Coord) :: x_i
    integer      :: id
    real (8)     :: lat, lon

    id   = idx(i, j, offs, dims)
    x_i  = dom%node%elts(id+1)

    ! Surfaced geopotential
    dom%surf_geopot%elts(id+1) = surf_geopot_fun(x_i)
  end subroutine set_surfgeopot

  function surf_pressure_fun (x_i)
    ! Surface pressure
    implicit none
    type(Coord) :: x_i
    real(8)     :: surf_pressure_fun

    real(8) :: lon, lat, surf_geopot

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    surf_pressure_fun = ref_surf_press * exp__flush ( &
         - radius*N_freq**2*u_0/(2.0_8*grav_accel**2*kappa)*(u_0/radius+f0)*(sin(lat)**2-1.0_8) &
         - N_freq**2/(grav_accel**2*kappa)*surf_geopot_fun (x_i) )
  end function surf_pressure_fun

  function surf_geopot_fun (x_i)
    ! Surface geopotential for Gaussian mountain (note that this really only needs to be done once)
    implicit none
    type(Coord) :: x_i
    real(8)     :: surf_geopot_fun

    real(8) :: lon, lat, rgrc

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))

    surf_geopot_fun = grav_accel * h_0 * exp__flush(-rgrc**2/d2)
  end function surf_geopot_fun

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    u = u_0*cos(lat)  ! Zonal velocity component
    v = 0.0_8         ! Meridional velocity component
  end subroutine vel_fun

  subroutine initialize_a_b_vert
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    if (allocated(a_vert)) deallocate(a_vert)
    if (allocated(b_vert)) deallocate(b_vert)
    if (allocated(a_vert_mass)) deallocate(a_vert_mass)
    if (allocated(b_vert_mass)) deallocate(b_vert_mass)
    allocate (a_vert(1:zlevels), b_vert(1:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          press_infty = 0.0_8
          a_vert(k) = real(k-1)/real(zlevels) * press_infty/ref_press
          b_vert(k) = 1.0_8 - real(k-1)/real(zlevels)
       end do
    else
       if (zlevels==18) then
          !a_vert=(/0.0_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
          a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1.0_8 /)
       elseif (zlevels==26) then
          !a_vert=(/0.0_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
               0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
               0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1.0_8 /)
       elseif (zlevels==49) then
          a_vert=(/0.002251865_8, 0.003983890_8, 0.006704364_8, 0.01073231_8, 0.01634233_8, 0.02367119_8, &
               0.03261456_8, 0.04274527_8, 0.05382610_8, 0.06512175_8, 0.07569850_8, 0.08454283_8, &
               0.08396310_8, 0.08334103_8, 0.08267352_8, 0.08195725_8, 0.08118866_8, 0.08036393_8, &
               0.07947895_8, 0.07852934_8, 0.07751036_8, 0.07641695_8, 0.07524368_8, 0.07398470_8, &
               0.07263375_8, 0.07118414_8, 0.06962863_8, 0.06795950_8, 0.06616846_8, 0.06424658_8, &
               0.06218433_8, 0.05997144_8, 0.05759690_8, 0.05504892_8, 0.05231483_8, 0.04938102_8, &
               0.04623292_8, 0.04285487_8, 0.03923006_8, 0.03534049_8, 0.03116681_8, 0.02668825_8, &
               0.02188257_8, 0.01676371_8, 0.01208171_8, 0.007959612_8, 0.004510297_8, 0.001831215_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
               0.006755112_8, 0.01400364_8, 0.02178164_8, 0.03012778_8, 0.03908356_8, 0.04869352_8, &
               0.05900542_8, 0.07007056_8, 0.08194394_8, 0.09468459_8, 0.1083559_8, 0.1230258_8, &
               0.1387673_8, 0.1556586_8, 0.1737837_8, 0.1932327_8, 0.2141024_8, 0.2364965_8, &
               0.2605264_8, 0.2863115_8, 0.3139801_8, 0.3436697_8, 0.3755280_8, 0.4097133_8, &
               0.4463958_8, 0.4857576_8, 0.5279946_8, 0.5733168_8, 0.6219495_8, 0.6741346_8, &
               0.7301315_8, 0.7897776_8, 0.8443334_8, 0.8923650_8, 0.9325572_8, 0.9637744_8, &
               0.9851122_8, 1.0_8/)
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

       ! Set mass coefficients
       b_vert_mass = b_vert(1:zlevels)-b_vert(2:zlevels+1)
       a_vert_mass = ((a_vert(1:zlevels)-a_vert(2:zlevels+1))*ref_press + b_vert_mass*press_infty)/grav_accel
    end if
  end subroutine initialize_a_b_vert

  subroutine read_test_case_parameters (filename)
    implicit none
    character(*)       :: filename
    integer, parameter :: fid = 500
    character(255)     :: varname
    
    open(unit=fid, file=filename, action='READ')
    read(fid,*) varname, max_level
    read(fid,*) varname, zlevels
    read(fid,*) varname, threshold 
    read(fid,*) varname, optimize_grid 
    read(fid,*) varname, dt_write
    read(fid,*) varname, CP_EVERY
    read(fid,*) varname, time_end
    read(fid,*) varname, resume

    if (rank==0) then
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
    implicit none
    integer :: iwrite, zlev

    integer      :: d, i, j, k, l, p, u
    character(7) :: var_file

    if (rank == 0) write(6,*) 'Saving fields'

    call update_array_bdry (sol, NONE)

    call pre_levelout

    ! Compute surface pressure
    call cal_surf_press (sol)

    do l = level_start, level_end
       minv = 1.0d63; maxv = -1.0d63
       u = 1000000+100*iwrite

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
          grid(d)%adj_mass%elts = 100.0_8
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), zlev,    0, 0)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, zlev)
          nullify (velo, vort)
       end do

       ! Calculate vorticity at hexagon points (stored in adj_mass)
       call apply_onescale (vort_triag_to_hex, l, z_null, 0, 1)

       Call write_level_mpi (write_primal, u+l, l, zlev, .True., test_case)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
       if (rank == 0) then
          write (var_file, '(i7)') u
          open(unit=50, file=trim(test_case)//'.'//var_file)
          write(50,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
          write(50,'(A, 7(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
          close(50)
       end if
       u = 2000000+100*iwrite
       call write_level_mpi (write_dual, u+l, l, zlev, .False., test_case)
    end do

    call post_levelout
    call barrier
    if (rank == 0) call compress_files (iwrite, test_case)

    ! Save 2D projection
    call export_2d (cart2sph2, 3000000+100*iwrite, (/-96, 96/), (/-48, 48/), (/2.0_8*MATH_PI, MATH_PI/), set_thresholds, test_case)
    !call export_2d (cart2sph2, 3000000+100*iwrite, (/-768, 768/), (/-384, 384/), (/2.0_8*MATH_PI, MATH_PI/), set_thresholds)
  end subroutine write_and_export

  subroutine dump (fid)
    implicit none
    integer :: fid
    
    write(fid) itime
    write(fid) iwrite
    write(fid) tol_mass, tol_temp, tol_velo
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid
    
    read(fid) itime
    read(fid) iwrite
    read(fid) tol_mass, tol_temp, tol_velo
  end subroutine load

  subroutine set_thresholds (itype)
    implicit none
    integer, optional :: itype

    integer :: l, k

    ! Set thresholds dynamically (trend or sol must be known)
    if (itype==0) then ! Adapt on trend
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

    if (istep/=0) then
       tol_mass = 0.99_8*tol_mass + 0.01_8*threshold * mass_scale
       tol_temp = 0.99_8*tol_temp + 0.01_8*threshold * temp_scale
       tol_velo = 0.99_8*tol_velo + 0.01_8*threshold * velo_scale
    elseif (istep==0) then
       tol_mass = threshold * mass_scale
       tol_temp = threshold * temp_scale
       tol_velo = threshold * velo_scale
       if (adapt_trend .and. itype==1) then ! Re-scale trend threshold for variables
          tol_mass = threshold**1.5_8 * mass_scale/5.0d1
          tol_temp = threshold**1.5_8 * temp_scale/5.0d1
          tol_velo = threshold**1.5_8 * velo_scale/5.0d1
       end if
    end if
  end subroutine set_thresholds

  subroutine linf_trend (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, e, k

    id = idx(i, j, offs, dims)

    ! Maximum trends
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
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
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, e, k

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
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
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, k

    id = idx(i, j, offs, dims)
    d = dom%id+1

    ! L2 norms of trends
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
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
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, k

    d = dom%id+1
    id = idx(i, j, offs, dims)

    ! L2 norms of trends
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
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
    implicit none
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
    implicit none
    integer ::  d, p

    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (set_surfgeopot, grid(d), p-1, z_null, 0, 1)
       end do
    end do
  end subroutine set_surf_geopot

  subroutine sum_total_mass (initialgo)
    ! Total mass over all vertical layers
    implicit none
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
end module DCMIP2008c5_mod

program DCMIP2008c5
  use main_mod
  use DCMIP2008c5_mod
  implicit none

  integer        :: d, ierr, k, l, v, zlev
  real(8)        :: dt_cfl, dt_visc
  character(255) :: command
  logical        :: aligned, remap, write_init

  ! Initialize grid etc
  call init_main_mod 

  ! Nullify all pointers initially
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, qe, vort)

  ! Read test case parameters
  call read_test_case_parameters (trim(test_case)//".in")

  ! Average minimum grid size and maximum wavenumber
  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8))

  ! Parameters for the simulation
  grav_accel     = 9.80616_8     ! gravitational acceleration in meters per second squared
  omega          = 7.29211d-5    ! Earth’s angular velocity in radians per second
  f0             = 2.0_8*omega   ! Coriolis parameter
  u_0            = 20.0_8        ! velocity in meters per second
  T_0            = 288.0_8       ! temperature in Kelvin
  d2             = 1.5d6**2      ! square of half width of Gaussian mountain profile in meters
  h_0            = 2.0d3         ! mountain height in meters
  lon_c          = MATH_PI/2.0_8 ! mountain peak longitudinal location in radians
  lat_c          = MATH_PI/6.0_8 ! mountain peak latitudinal location in radians
  radius         = 6.371229d6    ! mean radius of the Earth in meters
  ref_press      = 100145.6_8    ! reference pressure (mean surface pressure) in Pascals
  ref_surf_press = 930.0d2       ! South Pole surface pressure in Pascals
  R_d            = 287.04_8      ! ideal gas constant for dry air in joules per kilogram Kelvin
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

  zlev           = 6
  save_levels    = 1; allocate(pressure_save(1:save_levels))  ! number of vertical levels to save
  level_save     = level_end                                  ! resolution level at which to save lat-lon data
  pressure_save  = (/700.0d2/)                                ! interpolate values to this pressure level when interpolating to lat-lon grid

  ! Set logical switches
  adapt_trend  = .false. ! Adapt on trend or on variables
  adapt_dt     = .false.  ! Adapt time step
  compressible = .true.  ! Compressible equations
  remap        = .true. ! Remap vertical coordinates (always remap when saving results)
  uniform      = .false. ! Type of vertical grid
  
  ! Set viscosity
  Laplace_order = 1 ! Usual Laplacian diffusion
  !Laplace_order = 2 ! Iterated Laplacian diffusion
  
  if (Laplace_order==1) then ! Usual Laplacian diffusion
     viscosity_mass = 0.0_8
     viscosity_temp = viscosity_mass
     viscosity_divu = 0.0_8
     viscosity_rotu = viscosity_divu/1.0d2!visc * dx_min**2 ! viscosity for rotational part of momentum equation
  elseif (Laplace_order==2) then ! Second-order iterated Laplacian for diffusion
     viscosity_mass = 1.0d14!visc * dx_min**4/1.0e3 ! viscosity for mass equation
     viscosity_temp = viscosity_mass
     viscosity_divu = 1.0d14!visc * dx_min**4/1.0e3 ! viscosity for mass equation
     viscosity_rotu = viscosity_divu
  else
     write(6,*) 'Unsupported iterated Laplacian (only 1 or 2 supported)'
     stop
  end if
  viscosity = max (viscosity_mass, viscosity_temp, viscosity_divu, viscosity_rotu)

  ! Time step based on acoustic wave speed and hexagon edge length (not used if adaptive dt)
  dt_cfl = cfl_num*dx_min/(wave_speed+u_0)
  if (viscosity/=0.0_8) then
     dt_visc = 0.25_8*dx_min**2/viscosity
     dt_init = min(dt_cfl, dt_visc)
  else
     dt_init = dt_cfl
  end if
  if (rank==0)                      write(6,'(1(A,es10.4,1x))') "dt_cfl = ", dt_cfl
  if (rank==0.and.viscosity/=0.0_8) write(6,'(1(A,es10.4,1x))')" dt_visc = ", dt_visc

  if (rank == 0) then
     write(6,'(A,es10.4)') 'Viscosity_mass   = ', viscosity_mass
     write(6,'(A,es10.4)') 'Viscosity_temp   = ', viscosity_temp
     write(6,'(A,es10.4)') 'Viscosity_divu   = ', viscosity_divu
     write(6,'(A,es10.4)') 'Viscosity_rotu   = ', viscosity_rotu
     write(6,'(A,es10.4)') ' '
  end if

  ! Initialize vertical grid
  call initialize_a_b_vert

  ! Initialize variables
  call initialize (apply_initial_conditions, 1, set_thresholds, dump, load, test_case)

  allocate (n_patch_old(size(grid)), n_node_old(size(grid)))
  n_patch_old = 2;  call set_surf_geopot 

  call sum_total_mass (.True.)

  if (rank == 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:', tol_mass, tol_temp, tol_velo
  call barrier

  if (rank == 0) write(6,*) 'Write initial values and grid'
  call write_and_export (iwrite, zlev)

  if (resume<=0) iwrite = 0
  total_cpu_time = 0.0_8

  open(unit=12, file=trim(test_case)//'_log', action='WRITE', form='FORMATTED')
  if (rank == 0) then
     write (6,'(A,ES12.6,3(A,ES10.4),A,I2,A,I9)') &
          ' time [h] = ', time/3600.0_8, &
          '  mass tol = ', tol_mass, &
          ' temp tol = ', tol_temp, &
          ' velo tol = ', tol_velo, &
          ' Jmax =', level_end, &
          '  dof = ', sum(n_active)
  end if

  do while (time < time_end)
     call update_array_bdry (sol, NONE)
     n_patch_old = grid(:)%patch%length
     n_node_old = grid(:)%node%length

     if (remap .and. mod(istep, n_remap)==0 .and. istep>1) call remap_vertical_coordinates (set_thresholds)

     call start_timing
     call time_step (dt_write, aligned, set_thresholds)
     call stop_timing

     call set_surf_geopot
     timing = get_timing()
     total_cpu_time = total_cpu_time + timing

     if (rank == 0) then
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
             time/3600.0_8, dt, tol_mass, tol_temp, tol_velo, level_end, sum(n_active), mass_error, timing
     end if

     if (aligned) then
        iwrite = iwrite + 1

        ! Save fields
        if (remap) call remap_vertical_coordinates (set_thresholds)
        call write_and_export (iwrite, zlev)

        call sum_total_mass (.False.)

        if (modulo(iwrite,CP_EVERY) /= 0) cycle ! Do not write checkpoint

        ierr = write_checkpoint (dump)

        ! Let all cpus exit gracefully if NaN has been produced
        ierr = sync_max (ierr)
        if (ierr == 1) then ! NaN
           write (0,*) "NaN when writing checkpoint"
           call finalize
           stop
        end if

        ! Restart after checkpoint and load balance
        call restart_full (set_thresholds, load, test_case)
        call print_load_balance

        call barrier
     end if
     call sum_total_mass (.False.)
  end do

  if (rank == 0) then
     close (12)
     write (6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
     command = '\rm tmp tmp1 tmp2'; call system (command)
  end if

  call finalize
end program DCMIP2008c5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  ! Additional physics for the flux term of the scalar trend
  ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
  !
  ! NOTE: call with arguments (dom, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
  use domain_mod
  use DCMIP2008c5_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  integer                                  :: d, v
  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: grad
  real(8), dimension(:), pointer           :: flx
  logical                                  :: local_type
  
  interface
     function grad_physics (scalar, dom, id, idE, idNE, idN, type)
       use domain_mod
       use DCMIP2008c5_mod
       implicit none
       real(8), dimension(1:EDGE)               :: grad_physics
       real(8), dimension(:), pointer           :: scalar
       type(Domain)                             :: dom
       integer                                  :: id, idE, idNE, idN
       logical                                  :: type
     end function grad_physics
  end interface
  
  if (present(type)) then
     local_type = type
  else
     local_type = .false.
  end if
  
  if (max(viscosity_mass, viscosity_temp)==0.0_8) then
     physics_scalar_flux = 0.0_8
  else
     ! Calculate gradients
     if (Laplace_order==1) then
        grad(S_MASS,:) = grad_physics (mass, dom, id, idE, idNE, idN, local_type)
        grad(S_TEMP,:) = grad_physics (temp, dom, id, idE, idNE, idN, local_type)
     elseif (Laplace_order==2) then
        d = dom%id+1
        do v = S_MASS, S_TEMP
           flx => Laplacian_scalar(v)%data(d)%elts
           grad(v,:) = grad_physics (flx, dom, id, idE, idNE, idN, local_type)
           nullify (flx)
        end do
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

function grad_physics (scalar, dom, id, idE, idNE, idN, local_type)
  use domain_mod
  use DCMIP2008c5_mod
  implicit none

  real(8), dimension(1:EDGE)               :: grad_physics
  real(8), dimension(:), pointer           :: scalar
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical                                  :: local_type
  
  if (.not.local_type) then ! Usual gradient at edges of hexagon E, NE, N
     grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1) 
     grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1) 
     grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1) 
  else ! Gradient for southwest edges of hexagon W, SW, S
     grad_physics(RT+1) = -(scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*idE+RT+1) 
     grad_physics(DG+1) = -(scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*idNE+DG+1)
     grad_physics(UP+1) = -(scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*idN+UP+1) 
  end if
end function grad_physics

function physics_scalar_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the scalar trend
  ! Newton cooling to equilibrium potential temperature theta_equil
  use domain_mod
  use DCMIP2008c5_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source(S_MASS) = 0.0_8
  physics_scalar_source(S_TEMP) = 0.0_8
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the velocity trend
  !
  ! In this test case we add Rayleigh friction and Laplacian diffusion
  use domain_mod
  use ops_mod
  use DCMIP2008c5_mod
  implicit none
  
  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(Domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  integer                      :: e
    real(8), dimension(1:EDGE) :: diffusion,  curl_rotu, grad_divu
  
  if (max(viscosity_divu, viscosity_rotu)==0.0_8) then
     diffusion = 0.0_8
  else
     ! Calculate Laplacian of velocity
     grad_divu = gradi_e (divu, dom, i, j, offs, dims)
     curl_rotu = curlv_e (vort, dom, i, j, offs, dims)
     do e = 1, EDGE 
        diffusion(e) = viscosity_divu * grad_divu(e) - viscosity_rotu * curl_rotu(e)
     end do
  end if

  ! Total physics for source term of velocity trend
  do e = 1, EDGE
     physics_velo_source(e) =  diffusion(e)
  end do
end function physics_velo_source

