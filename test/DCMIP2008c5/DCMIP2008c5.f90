module DCMIP2008c5_mod
  ! DCMIP2008 test case 5 parameters
  use main_mod
  use remap_mod
  implicit none

  integer                            :: CP_EVERY, id, zlev, iwrite, j, N_node
  integer, dimension(:), allocatable :: n_patch_old, n_node_old
  real(8)                            :: Hmin, dh_min, dh_max, kmin
  real(8)                            :: initotalmass, totalmass, timing, total_cpu_time, dh
  logical                            :: wasprinted, uniform
  character (255)                    :: IC_file

  real(8)                            :: c_v, d2, h_0, lat_c, lon_c, N_freq, p_sp, T_0, u_0
  real(8)                            :: acceldim, f0, geopotdim, Ldim, Hdim, massdim, Tdim, Tempdim, Udim, pdim, R_ddim, specvoldim
  real(8)                            :: norm_mass, norm_temp, norm_velo, norm_mass_trend, norm_temp_trend, norm_velo_trend
  real(8)                            :: mass_scale, temp_scale, velo_scale, mass_scale_trend, temp_scale_trend, velo_scale_trend
  real(8)                            :: l2_mass, l2_temp, l2_velo
contains
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
    
    ! do d = 1, size(grid)
    !    do p = n_patch_old(d)+1, grid(d)%patch%length
    !       call apply_onescale_to_patch(set_surfgeopot, grid(d), p-1, z_null, 0, 1)
    !    end do
    ! end do

     do d = 1, size(grid)
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (set_surfgeopot, grid(d), p-1, z_null, 0, 1)
          end do
          nullify (mass, temp, dvelo, bernoulli, exner)
       end do
  end subroutine set_surf_geopot

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

  subroutine write_and_print_step
    real(4) timing
    timing = get_timing()
    if (rank .eq. 0) write(1011,'(3(ES13.4,1X), I3, 2(1X, I9), 1(1X,ES13.4))') &
         time, dt, timing, level_end, n_active
  end subroutine write_and_print_step

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
       if (zlevels.eq.18) then
          !a_vert=(/0.0_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               a_vert=(/0.00251499_8, 0.00710361_8, 0.01904260_8, 0.04607560_8, 0.08181860_8, &
               0.07869805_8, 0.07463175_8, 0.06955308_8, 0.06339061_8, 0.05621774_8, 0.04815296_8, &
               0.03949230_8, 0.03058456_8, 0.02193336_8, 0.01403670_8, 0.007458598_8, 0.002646866_8, &
               0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.03756984_8, 0.08652625_8, 0.1476709_8, 0.221864_8, &
               0.308222_8, 0.4053179_8, 0.509588_8, 0.6168328_8, 0.7209891_8, 0.816061_8, 0.8952581_8, &
               0.953189_8, 0.985056_8, 1.0_8 /)
       elseif (zlevels.eq.26) then
          !a_vert=(/0.0_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
          a_vert=(/0.002194067_8, 0.004895209_8, 0.009882418_8, 0.01805201_8, 0.02983724_8, 0.04462334_8, 0.06160587_8, &
               0.07851243_8, 0.07731271_8, 0.07590131_8, 0.07424086_8, 0.07228744_8, 0.06998933_8, 0.06728574_8, 0.06410509_8, &
               0.06036322_8, 0.05596111_8, 0.05078225_8, 0.04468960_8, 0.03752191_8, 0.02908949_8, 0.02084739_8, 0.01334443_8, &
               0.00708499_8, 0.00252136_8, 0.0_8, 0.0_8 /)
          b_vert=(/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.01505309_8, 0.03276228_8, 0.05359622_8, &
                0.07810627_8, 0.1069411_8, 0.1408637_8, 0.1807720_8, 0.2277220_8, 0.2829562_8, 0.3479364_8, 0.4243822_8, &
                0.5143168_8, 0.6201202_8, 0.7235355_8, 0.8176768_8, 0.8962153_8, 0.9534761_8, 0.9851122_8, 1.0_8 /)
       elseif (zlevels.eq.49) then
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
    end if
  end subroutine initialize_a_b_vert

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

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = &
         ((a_vert(zlev)-a_vert(zlev+1))*ref_press + (b_vert(zlev)-b_vert(zlev+1))*dom%surf_press%elts(id+1))/grav_accel
    
    ! Horizontally uniform potential temperature
    pot_temp =  T_0 * (lev_press/ref_press)**(-kappa)
    
    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)

    ! Print out layer thicknesses
    ! if (rank.eq.0 .and. .not.wasprinted) then
    !    if(zlev.eq.1) then
    !       write(6,'(4(A,es11.4))') &
    !            ' surf geopot =', dom%surf_geopot%elts(id+1), &
    !            ' surf press =', dom%surf_press%elts(id+1), &
    !            ' press_infty =', press_infty
    !    end if

    !    write(6,'(A,I2,1x,6(A,es11.4))') &
    !         ' zlev = ', zlev, &
    !         ' press =', lev_press, &
    !         ' mu =', sol(S_MASS,zlev)%data(d)%elts(id+1), &
    !         ' Theta =', sol(S_TEMP,zlev)%data(d)%elts(id+1), &
    !         ' U = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1), &
    !         ' V = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1), &
    !         ' W = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1)
    !    wasprinted=.true.
    ! end if
  end subroutine init_sol

  subroutine set_surfgeopot (dom, i, j, zlev, offs, dims)
    ! Initialize surface geopotential after restart
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

  function surf_geopot_fun (x_i)
    ! Surface geopotential for Gaussian mountain (note that this really only needs to be done once)
    type(Coord) :: x_i
    real(8)     :: surf_geopot_fun
    
    real(8) :: lon, lat, rgrc
    
    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph(x_i, lon, lat)

    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))
    
    surf_geopot_fun = grav_accel * h_0 * exp__flush(-rgrc**2/d2)
  end function surf_geopot_fun

  function surf_pressure_fun (x_i)
    ! Surface pressure
    type(Coord) :: x_i
    real(8)     :: surf_pressure_fun
    
    real(8) :: lon, lat, surf_geopot

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph(x_i, lon, lat)
    
    surf_pressure_fun = p_sp * exp__flush ( &
         - radius*N_freq**2*u_0/(2.0_8*grav_accel**2*kappa)*(u_0/radius+f0)*(sin(lat)**2-1.0_8) &
         - N_freq**2/(grav_accel**2*kappa)*surf_geopot_fun (x_i) )
  end function surf_pressure_fun
  
  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    real(8) :: lon, lat, u, v
    
    u = u_0*cos(lat)  ! Zonal velocity component
    v = 0.0_8         ! Meridional velocity component
  end subroutine vel_fun

  subroutine read_test_case_parameters(filename)
    character(*) filename
    integer :: fid = 500
    character(255) varname
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
    dt_write = dt_write * 60_8!/Tdim
    time_end = time_end * 60_8**2!/Tdim

    close(fid)
  end subroutine read_test_case_parameters

  subroutine write_and_export(iwrite)
    integer :: iwrite
    
    integer :: l, k, zlev, d, u, i, p

    call update_array_bdry (sol, NONE)

    call pre_levelout

    zlev = 6 ! Only export one vertical level

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
             exner => grid(d)%exner%elts
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
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, zlev)
          nullify (velo, vort)
       end do

       Call write_level_mpi (write_primal, u+l, l, zlev, .True.)

       do i = 1, N_VAR_OUT
          minv(i) = -sync_max_d(-minv(i))
          maxv(i) =  sync_max_d( maxv(i))
       end do
       if (rank .eq. 0) write(u,'(A, 5(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
       if (rank .eq. 0) write(u,'(A, 5(E15.5E2, 1X), I3)') "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
       u = 200000+100*iwrite
       call write_level_mpi (write_dual, u+l, l, zlev, .False.)
    end do

    call post_levelout
    call barrier
    !if (rank .eq. 0) call compress_files(k)
  end subroutine write_and_export

  subroutine DCMIP2008c5_dump(fid)
    integer :: fid
    write(fid) itime
    write(fid) iwrite
    write(fid) tol_mass, tol_temp, tol_velo
  end subroutine DCMIP2008c5_dump

  subroutine DCMIP2008c5_load(fid)
    integer :: fid
    read(fid) itime
    read(fid) iwrite
    read(fid) tol_mass, tol_temp, tol_velo
  end subroutine DCMIP2008c5_load
  
  subroutine set_thresholds (itype)
    integer, optional :: itype
    
    integer :: l, k

    ! Set thresholds dynamically (trend or sol must be known)
    if (adapt_trend) then
       norm_mass_trend = 0.0_8
       norm_temp_trend = 0.0_8
       norm_velo_trend = 0.0_8
       do l = level_start, level_end
          call apply_onescale (linf_trend, l, z_null, 0, 0)
       end do
       
       mass_scale = sync_max_d(norm_mass_trend)
       temp_scale = sync_max_d(norm_temp_trend)
       velo_scale = sync_max_d(norm_velo_trend)
    else
       norm_mass = 0.0_8
       norm_temp = 0.0_8
       norm_velo = 0.0_8
       do l = level_start, level_end
          call apply_onescale (linf_vars, l, z_null, 0, 0)
       end do
       
       mass_scale = sync_max_d(norm_mass)
       temp_scale = sync_max_d(norm_temp)
       velo_scale = sync_max_d(norm_velo)
    end if
   
    if (istep.ne.0) then
       tol_mass = 0.9_8*tol_mass + 0.1_8*threshold * mass_scale
       tol_temp = 0.9_8*tol_temp + 0.1_8*threshold * temp_scale
       tol_velo = 0.9_8*tol_velo + 0.1_8*threshold * velo_scale
    elseif (istep.eq.0) then
       if (adapt_trend .and. itype.eq.1) then
          tol_mass = threshold * dt_init*mass_scale * 8_8
          tol_temp = threshold * dt_init*temp_scale * 8_8
          tol_velo = threshold * dt_init*velo_scale * 8_8
       else
          tol_mass = threshold * mass_scale
          tol_temp = threshold * temp_scale
          tol_velo = threshold * velo_scale
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

    d = dom%id+1
    id = idx(i, j, offs, dims)

    ! L2 norms of trends
    if (dom%mask_n%elts(id+1) .ge. ADJZONE) then
       N_node = N_node + 1
       do k = 1, zlevels
          norm_mass_trend = norm_mass_trend + trend(S_MASS,zlev)%data(d)%elts(id+1)**2
          norm_temp_trend = norm_temp_trend + trend(S_TEMP,zlev)%data(d)%elts(id+1)**2
          do e = 1, EDGE
             norm_velo_trend  = norm_velo_trend + trend(S_VELO,zlev)%data(d)%elts(EDGE*id+e)**2
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

end module DCMIP2008c5_mod

program DCMIP2008c5
  use main_mod
  use DCMIP2008c5_mod
  implicit none

  integer                      :: d, ierr, k, l, v
  integer, parameter           :: len_cmd_files = 12 + 4 + 12 + 4
  integer, parameter           :: len_cmd_archive = 11 + 4 + 4
  real(8)                      :: Area_min, total_mass, visc, wave_speed
  character(len_cmd_files)     :: cmd_files
  character(len_cmd_archive)   :: cmd_archive
  character(8+8+29+14)         :: command
  character(9+len_cmd_archive) :: command1
  character(6+len_cmd_files)   :: command2
  logical                      :: aligned, remap, write_init

  ! Initialize grid etc
  call init_main_mod 

  ! Nullify all pointers initially
  nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t)
  nullify (bernoulli, exner, qe)
  nullify (adj_temp_up, adj_mass_up, adj_velo_up, v_mflux)

  ! Read test case parameters
  call read_test_case_parameters("DCMIP2008c5.in")

  call initialize_a_b_vert

  ! Average minimum grid size and maximum wavenumber
  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8))
  Area_min = sqrt(3.0_8)/2.0_8*dx_min**2
  dx_max = 2.0_8*MATH_PI * radius
  kmin = MATH_PI/dx_max ; kmax = MATH_PI/dx_min

  ! Parameters for the simulation
  grav_accel  = 9.80616_8     ! gravitational acceleration in meters per second squared
  omega       = 7.29211d-5    ! Earth’s angular velocity in radians per second
  f0          = 2.0_8*omega   ! Coriolis parameter
  u_0         = 20.0_8        ! velocity in meters per second
  T_0         = 288.0_8       ! temperature in Kelvin
  d2          = 1500d3**2     ! square of half width of Gaussian mountain profile in meters
  h_0         = 2000.0_8      ! mountain height in meters
  lon_c       = MATH_PI/2.0_8 ! mountain peak longitudinal location in radians
  lat_c       = MATH_PI/6.0_8 ! mountain peak latitudinal location in radians
  p_sp        = 930.0d2       ! South Pole surface pressure in Pascals
  radius      = 6.371229d6    ! mean radius of the Earth in meters
  ref_press   = 100145.6_8    ! reference pressure (mean surface pressure) in Pascals
  R_d         = 287.04_8      ! ideal gas constant for dry air in joules per kilogram Kelvin
  c_p         = 1004.64_8     ! specific heat at constant pressure in joules per kilogram Kelvin
  c_v         = 717.6_8       ! specfic heat at constant volume c_v = R_d - c_p
  kappa       = 2.0_8/7.0_8   ! kappa=R_d/c_p
  N_freq      = sqrt(grav_accel**2/(c_p*T_0)) ! Brunt-Vaisala buoyancy frequency

  ! Dimensional scaling
  Ldim        = sqrt(d2)                         ! horizontal length scale
  Hdim        = h_0                              ! vertical length scale
  Udim        = u_0                              ! velocity scale
  Tdim        = Hdim/Udim                        ! time scale
  Tempdim     = T_0                              ! temperature scale (both theta and T from DYNAMICO)

  acceldim    = Udim*Udim/Hdim                   ! acceleration scale
  pdim        = ref_press                        ! pressure scale
  R_ddim      = R_d                              ! R_d scale
  massdim     = pdim*Hdim/(Tempdim*R_d)          ! mass (=rho*dz following DYNAMICO) scale
  specvoldim  = (R_d*Tempdim)/pdim               ! specific volume scale
  geopotdim   = acceldim*massdim*specvoldim/Hdim ! geopotential scale

  cfl_num     = 0.8d0                            ! cfl number
  n_diffuse   = 1                                ! Diffusion step interval
  n_remap     = 1                                ! Vertical remap interval
  
  dt_init     = 500.0_8                          ! Time step (not used if adapt_dt is true)

  ray_friction = 0.0_8!1_8/25_8                        ! Rayleigh friction

  viscosity_mass = 1.0d-3/kmax**2                ! viscosity for mass equation
  viscosity_temp = 1.0d-3/kmax**2                ! viscosity for mass-weighted potential temperature equation
  viscosity_divu = 1.0d-3/kmax**2                ! viscosity for divergent part of momentum equation
  viscosity_rotu = 1.0d-3/kmax**2                ! viscosity for divergent part of momentum equation

  if (rank .eq. 0) then
     write(6,'(A,es10.4)') 'Viscosity_mass   = ', viscosity_mass
     write(6,'(A,es10.4)') 'Viscosity_temp   = ', viscosity_temp
     write(6,'(A,es10.4)') 'Viscosity_divu   = ', viscosity_divu
     write(6,'(A,es10.4)') 'Viscosity_rotu   = ', viscosity_rotu
     write(6,'(A,es10.4)') ' '
  end if

  ! Set logical switches
  adapt_trend      = .false. ! Adapt on trend or on variables
  adapt_dt         = .true.  ! Adapt time step
  diffuse          = .true.  ! Diffuse scalars
  compressible     = .true.  ! Compressible equations
  remap            = .false.  ! Remap vertical coordinates
  uniform          = .false. ! Type of vertical grid

  ! Initialize variables
  call initialize (apply_initial_conditions, 1, set_thresholds, DCMIP2008c5_dump, DCMIP2008c5_load)

  allocate (n_patch_old(size(grid)), n_node_old(size(grid)))
  n_patch_old = 2;  call set_surf_geopot 

  call sum_total_mass (.True.)

  if (rank .eq. 0) write (6,'(A,3(ES12.4,1x))') 'Thresholds for mass, temperature, velocity:',  tol_mass, tol_temp, tol_velo
  call barrier

  if (rank .eq. 0) write(6,*) 'Write initial values and grid'
  call write_and_export (iwrite)

  if(resume.le.0) iwrite = 0
  total_cpu_time = 0.0_8

  open(unit=12, file='DCMIP2008c5_log', action='WRITE', form='FORMATTED')
  if (rank .eq. 0) then
     write (6,'(A,ES12.6,4(A,ES10.4),A,I2,A,I9,A,ES9.2)') &
          ' time [h] = ', time/3600.0_8, &
          ' dt [s] = ', dt_init, &
          '  mass tol = ', tol_mass, &
          ' temp tol = ', tol_temp, &
          ' velo tol = ', tol_velo, &
          ' Jmax =', level_end, &
          '  dof = ', sum(n_active), &
          ' cpu = ', timing

     write (12,'(5(ES15.9,1x),I2,1X,I9,1X,ES14.8)')  &
          time/3600.0_8, dt_init, tol_mass, tol_temp, tol_velo, level_end, sum(n_active), timing
  end if
  do while (time .lt. time_end)
     call start_timing
     call update_array_bdry (sol, NONE)
     n_patch_old = grid(:)%patch%length
     n_node_old = grid(:)%node%length

     call time_step (dt_write, aligned, set_thresholds)

     if (remap .and. mod(istep, n_remap).eq.0) then
        if (rank.eq.0) write(6,*) 'Remapping vertical coordinates'
       call remap_vertical_coordinates
     end if

     call set_surf_geopot
     call stop_timing
     timing = get_timing()
     total_cpu_time = total_cpu_time + timing
     
     call write_and_print_step
     
     if (rank .eq. 0) then
        write (6,'(A,ES12.6,4(A,ES10.4),A,I2,A,I9,A,ES9.2)') &
             ' time [h] = ', time/3600.0_8, &
             ' dt [s] = ', dt, &
             '  mass tol = ', tol_mass, &
             ' temp tol = ', tol_temp, &
             ' velo tol = ', tol_velo, &
             ' Jmax =', level_end, &
             '  dof = ', sum(n_active), &
             ' cpu = ', timing

        write (12,'(5(ES15.9,1x),I2,1X,I9,1X,ES14.8)')  &
             time/3600.0_8, dt, tol_mass, tol_temp, tol_velo, level_end, sum(n_active), timing
     end if

     call print_load_balance
     
     if (aligned) then
        iwrite = iwrite + 1
        ! Remap to original vertical coordinates before saving data or checkpoint
        !call remap_vertical_coordinates
        if (rank.eq.0) write(6,*) 'Saving fields'
        call write_and_export (iwrite)

        if (modulo(iwrite,CP_EVERY) .ne. 0) cycle 
        ierr = write_checkpoint (DCMIP2008c5_dump)

        ! let all cpus exit gracefully if NaN has been produced
        ierr = sync_max(ierr)
        if (ierr .eq. 1) then ! NaN
           write (0,*) "NaN when writing checkpoint"
           call finalize
           stop
        end if

        call restart_full (set_thresholds, DCMIP2008c5_load)

        ! deallocate(n_patch_old); allocate(n_patch_old(size(grid)))
        ! deallocate(n_node_old);  allocate(n_node_old(size(grid)))
        ! n_patch_old = 2; call set_surf_geopot

        call barrier
     end if

     call sum_total_mass(.False.)
  end do

  if (rank .eq. 0) then
     write(6,'(A,ES11.4)') 'Total cpu time = ', total_cpu_time
     close(12)
     close(1011)
     close(8450)
  end if

  call finalize
end program DCMIP2008c5
