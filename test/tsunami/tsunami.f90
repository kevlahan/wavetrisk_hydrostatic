module tsunami_mod
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

  real(8) :: csq, c_p

  real(8) :: U 
  real(8) :: Fr

  real(8) :: VELO_SCALE

  ! Coordinates of initial condition
  real(8) :: LON_MIN, LON_MAX, LAT_MIN, LAT_MAX

  real(8), parameter :: LAND = 1
  real(8), parameter :: SEA  = 0
  character(255) IC_file

  integer :: CP_EVERY 

  real(8) :: Hmin, eta, alpha, dh_min, dh_max, dx_min, dx_max, kmin, k_tsu

  integer, dimension(2) :: OKADA_DIM
  integer :: BATHY_PER_DEG, npts_chi, npts_topo
  real(8), allocatable :: okada_data(:,:)
  real(4), allocatable :: bathy_data(:,:)

  type(Float_Field) arrival, wave_h

  integer, allocatable :: n_patch_old(:), n_node_old(:)

  logical const_bathymetry

  logical :: calc_tide_gauges
  integer iwrite, j
  integer, parameter :: n_gauge = 22
  real(8) max_height
  real(8) cur_gauge, last_gauge_dist, glo_gauge_dist
  real(8), dimension (1:n_gauge) :: tide_record
  real(8), dimension (1:n_gauge,1:2) :: station_coord
  type(Coord) :: gauge_coord

contains
  subroutine apply_initial_conditions()
      integer l, d, p
      do l = level_start, level_end
          call apply_onescale(init_sol, l, 0, 1)
      end do
      do d = 1, size(grid)
          do p = 3, grid(d)%patch%length
              call apply_onescale_to_patch(cpt_topo_penal, grid(d), p-1, -2, 3)
          end do
      end do
      do l = level_start, level_end
          ! FIXME: works only for zero initial heigth perturbation at poles
          if (penalize) call apply_onescale(penalize_ic, l, 0, 0) 
      end do
  end subroutine

  subroutine write_and_print_step()
      real(4) timing
      timing = get_timing()
      if (rank .eq. 0) write(1011,'(3(ES13.4,1X), I3, 2(1X, I9), 2(1X,ES13.4))') &
              time, dt, timing, level_end, n_active, max_height, VELO_SCALE
  end subroutine

  subroutine cpt_max_dh(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      id = idx(i, j, offs, dims)
      if (dom%mask_p%elts(id+1) .gt. 0) then
          if (dom%level%elts(id+1) .eq. level_end .or. dom%mask_p%elts(id+1) .eq. ADJZONE) &
              max_height = max(max_height, abs(sol(S_HEIGHT)%data(dom%id+1)%elts(id+1)))
      end if
  end subroutine

  subroutine penalize_ic(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      d = dom%id+1
      id = idx(i, j, offs, dims)
      sol(S_HEIGHT)%data(d)%elts(id+1) = sol(S_HEIGHT)%data(d)%elts(id+1)* &
                  (1+alpha_m1*penal%data(d)%elts(id+1))
  end subroutine

  subroutine init_sol(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id, d
      real(8) lon, lat
      real(8) s, t
      d = dom%id+1
      id = idx(i, j, offs, dims)
      call cart2sph(dom%node%elts(id+1), lon, lat)
      lon = lon*180.0_8/MATH_PI
      lat = lat*180.0_8/MATH_PI
      if (LON_MIN .lt. lon .and. lon .lt. LON_MAX .and. LAT_MIN .lt. lat .and. lat .lt. LAT_MAX) then
          s = (lon-LON_MIN)/(LON_MAX-LON_MIN)*OKADA_DIM(1)
          t = (lat-LAT_MIN)/(LAT_MAX-LAT_MIN)*OKADA_DIM(2)
          s = ceiling(s)
          t = ceiling(t)
          if (nint(s) .lt. lbound(okada_data,1) .or. nint(s) .gt. ubound(okada_data,1) .or. &
              nint(t) .lt. lbound(okada_data,2) .or. nint(t) .gt. ubound(okada_data,2)) stop "okada out of bound"
          sol(S_HEIGHT)%data(d)%elts(id+1) = okada_data(nint(s),nint(t))
      else
          sol(S_HEIGHT)%data(d)%elts(id+1) = 0
      end if
      ! to be double save (u is initialized to zero by default):
      sol(S_VELO)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0
  end subroutine
    
  subroutine read_test_case_parameters(filename)
      character(*) filename
      integer :: fid = 500
      character(255) varname
      open(unit=fid, file=filename, action='READ')
      read(fid,*) varname, max_level 
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
      read(fid,*) varname, IC_file
      read(fid,*) varname, LON_MIN
      read(fid,*) varname, LON_MAX
      read(fid,*) varname, LAT_MIN
      read(fid,*) varname, LAT_MAX
      read(fid,*) varname, BATHY_PER_DEG
      read(fid,*) varname, OKADA_DIM(1)
      read(fid,*) varname, OKADA_DIM(2)
      read(fid,*) varname, npts_chi
      read(fid,*) varname, npts_topo

      if (rank.eq.0) then
         write(*,'(A,i3)')     "max_level        = ", max_level
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
         write(*,'(A,A)')      "IC_file          = ", IC_file
         write(*,'(A,es11.4)') "LON_MIN          = ", LON_MIN
         write(*,'(A,es11.4)') "LON_MAX          = ", LON_MAX
         write(*,'(A,es11.4)') "LAT_MIN          = ", LAT_MIN
         write(*,'(A,es11.4)') "LAT_MAX          = ", LAT_MAX
         write(*,'(A,i6)')     "BATHY_PER_DEG    = ", BATHY_PER_DEG
         write(*,'(A,i6)')     "OKADA_DIM(1)     = ", OKADA_DIM(1)
         write(*,'(A,i6)')     "OKADA_DIM(2)     = ", OKADA_DIM(2)
         write(*,'(A,i4)')     "npts_chi         = ", npts_chi
         write(*,'(A,i4)')     "npts_topo        = ", npts_topo
         write(*,*) ' '
      end if
      dt_write = dt_write * 60_8/Tdim
      time_end = time_end * 60_8**2/Tdim
      close(fid)
  end subroutine

  subroutine cpt_topo_penal(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer e, grid_level
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      real(8) lon, lat
      real(8) s, t
      real(8) geo, M_chi, b, h, dx_local
      real(8), dimension(3) :: dx_primal, dx_dual

      id = idx(i, j, offs, dims)
      if (.not. penalize) then
          dom%topo%elts(id+1) = 1
          return
      end if
      
      do e = 1, 3
         dx_primal(e) = dom%len%elts(id*EDGE+e)
         dx_dual(e) = dom%pedlen%elts(id*EDGE+e)
      end do
      dx_local = max(maxval(dx_primal), maxval(dx_dual))
      if (dx_local.eq.0.0_8) dx_local = dx_min

      call cart2sph(dom%node%elts(id+1), lon, lat)
      s = lon/MATH_PI*dble(180*BATHY_PER_DEG)
      t = lat/MATH_PI*dble(180*BATHY_PER_DEG)
      call smoothed_penal_topo(s, t, b, penal%data(dom%id+1)%elts(id+1), dx_local)
      if (const_bathymetry) b = 0
      dom%topo%elts(id+1) = 1.0_8 + b
  end subroutine

  subroutine smoothed_penal_topo(s0, t0, stopo, spenal, dx_smooth)
      real(8) s0, t0, stopo, spenal, dx_smooth
      integer s, t, is0, it0
      integer i, j
      real(8) chi_sum, topo_sum, sw_chi, sw_topo, M_chi, M_topo, r, wgt_chi, wgt_topo
      type(Coord) :: p, q

      p = proj_lon_lat(s0,t0)
      is0=nint(s0); it0=nint(t0)

      ! Smooth penalization mask
      if (npts_chi.eq.0) then ! no smooth
         if (bathy_data(is0,it0) > 0.0_8) then
            spenal = LAND
         else
            spenal = SEA
         end if
      else ! smooth
         sw_chi   = 0.0_8
         chi_sum  = 0.0_8
         do i = -npts_chi, npts_chi
            do j = -npts_chi, npts_chi
               s = is0+i ; t = it0+j
               call wrap_lonlat(s, t)
               q = proj_lon_lat(dble(s), dble(t))
               
               r = norm(vector(p,q))
               wgt_chi  = radial_basis_fun(r, npts_chi,  dx_smooth)
               
               if (bathy_data(s, t) > 0.0_8) then ! land
                  M_chi = LAND ! = 1
               else ! sea: geo < 0
                  M_chi = SEA ! = 0
               end if
               chi_sum  = chi_sum  + M_chi *wgt_chi
               sw_chi  = sw_chi  + wgt_chi
            end do
         end do
         spenal = chi_sum/sw_chi
      end if

      if (npts_topo.eq.0) then 
         if (bathy_data(is0,it0) > 0.0_8) then
            stopo = 0.0_8
         else
            stopo = (-min(bathy_data(s, t), -Hmin) -H_star)/Hdim
         end if
      else
         sw_topo  = 0.0_8
         topo_sum = 0.0_8
         do i = -npts_topo, npts_topo
            do j = -npts_topo, npts_topo
               s = is0+i ; t = it0+j
               call wrap_lonlat(s, t)
               q = proj_lon_lat(dble(s), dble(t))
               
               r = norm(vector(p,q))
               wgt_topo = radial_basis_fun(r, npts_topo, dx_smooth)
               
               if (bathy_data(s, t) > 0.0_8) then ! land
                  M_topo = 0.0_8
               else ! sea: geo < 0
                  M_topo = (-min(bathy_data(s, t), -Hmin) -H_star)/Hdim
               end if
               topo_sum = topo_sum + M_topo*wgt_topo
               sw_topo = sw_topo + wgt_topo
            end do
         end do
         stopo  = topo_sum/sw_topo
      end if
  end subroutine

  type(Coord) function proj_lon_lat(s,t)
    real(8) :: s, t
    real(8) :: lon, lat
    
    lon = s*MATH_PI/dble(180*BATHY_PER_DEG)
    lat = t*MATH_PI/dble(180*BATHY_PER_DEG)
    proj_lon_lat = project_on_sphere(sph2cart(lon, lat))
  end function proj_lon_lat 

  real(8) function radial_basis_fun(r, npts, dx_local)
      real(8) r, alph, dx_local
      integer :: npts
      alph = 1.0_8 / (npts/2 * dx_local)
            
      radial_basis_fun = exp(-(alph*r)**2)
  end function

  subroutine wrap_lonlat(s, t)
  ! longitude: wraparound allows for values outside [-180,180]
  ! latitude: works only if there is no cost at the pole
      integer s, t
      if (t .lt. lbound(bathy_data,2)) t = lbound(bathy_data,2) ! pole
      if (t .gt. ubound(bathy_data,2)) t = ubound(bathy_data,2) ! pole
      if (s .lt. lbound(bathy_data,1)) s = s + 360*BATHY_PER_DEG
      if (s .gt. ubound(bathy_data,1)) s = s - 360*BATHY_PER_DEG
  end subroutine

  subroutine finish_new_patches()
      integer d, p
      do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
              call apply_onescale_to_patch(cpt_topo_penal, grid(d), p-1, -2, 3)
          end do
      end do
  end subroutine

  subroutine first_arrival(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      id = idx(i, j, offs, dims)
      if (sol(S_HEIGHT)%data(dom%id+1)%elts(id+1)*Hdim .gt. 0.05 &
              .and. dom%mask_p%elts(id+1) .ge. ADJZONE) &
          arrival%data(dom%id+1)%elts(id+1) = min(time*Tdim, arrival%data(dom%id+1)%elts(id+1))
  end subroutine

  subroutine wave_height(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      id = idx(i, j, offs, dims)
      if (dom%mask_p%elts(id+1) .gt. 0) then
          wave_h%data(dom%id+1)%elts(id+1) = &
              max(sol(S_HEIGHT)%data(dom%id+1)%elts(id+1)*Hdim, wave_h%data(dom%id+1)%elts(id+1))
      end if
  end subroutine

  subroutine tide_gauge(dom, i, j, offs, dims)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      real(8) :: cur_gauge_dist

      id = idx(i, j, offs, dims)

      if (dom%mask_p%elts(id+1) .gt. 0) then
         cur_gauge_dist = norm(vector(dom%node%elts(id+1), gauge_coord))
         if (cur_gauge_dist .lt. last_gauge_dist) then
            cur_gauge = sol(S_HEIGHT)%data(dom%id+1)%elts(id+1)*Hdim
            last_gauge_dist = cur_gauge_dist
         end if
      end if
  end subroutine

  subroutine write_and_export(k)
      integer l, k
      integer u, i
      call trend_ml(sol, trend)
      call pre_levelout()
      do l = level_start, level_end
          minv = 1.d63;
          maxv = -1.d63;
          u = 100000+100*k
          call write_level_mpi(write_primal, u+l, l, .True.)
          do i = 1, N_VAR_OUT
              minv(i) = -sync_max_d(-minv(i))
              maxv(i) = sync_max_d(maxv(i))
          end do
          if (rank .eq. 0) write(u,'(A, 4(E15.5E2, 1X), I3)') &
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", minv, l
          if (rank .eq. 0) write(u,'(A, 4(E15.5E2, 1X), I3)') &
                  "0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 ", maxv, l
          u = 200000+100*k
          ! call write_level_mpi(write_dual, u+l, l, .False.)
      end do
      call post_levelout()
      call barrier
      if (rank .eq. 0) call compress_files(k) 
  end subroutine

  subroutine cart2sph2(cin, cout)
      type(Coord) cin
      real(8), intent(out) :: cout(2)
      call cart2sph(cin, cout(1), cout(2))
  end subroutine

  subroutine tsunami_dump(fid)
      integer fid
      write(fid) VELO_SCALE
      write(fid) iwrite
  end subroutine

  subroutine tsunami_load(fid)
      integer fid
      read(fid) VELO_SCALE
      read(fid) iwrite
  end subroutine

  subroutine set_thresholds() ! inertia-gravity wave
      toll_height = VELO_SCALE*c_p/grav_accel * threshold**(3.0_8/2.0_8)
      toll_velo   = VELO_SCALE     * threshold**(3.0_8/2.0_8)
  end subroutine
end module

program tsunami
  use main_mod
  use tsunami_mod
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
  call read_test_case_parameters("tsunami.in")

  if (resume .eq. NONE) then
      if (rank .eq. 0) write(*,*) 'Reading okada initial height data'
      allocate(okada_data(OKADA_DIM(1),OKADA_DIM(2)))
      open(unit=1085,file=IC_file)
      do k = OKADA_DIM(2), 1, -1 ! north to south (as read from file)
          read(1085,*) okada_data(:,k)
      end do
      close(1085)
      dh_max = maxval(okada_data)
      dh_min = minval(okada_data)
  end if

  ! Shared non-dimensional parameters
  radius     = R_star / Ldim
  grav_accel = g_star * Hdim/Udim**2
  omega      = Omega_star * Ldim/Udim
  
  dx_min = sqrt(4.0_8*MATH_PI*radius**2/(10.0_8*4**max_level+2.0_8)) ! Average minimum grid size
  dx_max = 2.0_8*MATH_PI * radius

  kmin = MATH_PI/dx_max ; kmax = 2.0_8*MATH_PI/dx_max
  
  csq = grav_accel*H
  k_tsu = 2.0_8*MATH_PI/(1e6_8/Ldim) ! Approximate wavelength of tsunami: 100km
  c_p = sqrt(f0**2/k_tsu**2 + csq) ! Maximum phase wave speed
   
  U = grav_accel*0.5_8*(dh_max-dh_min)/c_p ! Characteristic velocity based on initial perturbation
  Fr = U/c_p ! Froude number 
  VELO_SCALE   = U
 
  wind_stress      = .False.
  penalize         = .True.
  bottom_friction  = .False.
  calc_tide_gauges = .False.

  if (rank.eq.0) then
     write(*,'(A,L1)') "wind_stress     = ", wind_stress
     write(*,'(A,L1)') "penalize        = ", penalize
     write(*,'(A,L1)') "bottom friction = ", bottom_friction
     write(*,'(A,L1)') "tide gauges     = ", calc_tide_gauges
  end if

  viscosity = 1.0_8/((2.0_8*MATH_PI/dx_min)/64.0_8)**2     ! grid scale viscosity
  friction_coeff = 3e-3_8 ! Bottom friction
  if (rank .eq. 0) write (*,'(A,es11.4)') 'Viscosity = ',  viscosity

  write_init = (resume .eq. NONE)
  iwrite = 0

  ! Tide station coordinates (shifted by a few km off-shore) from Merrifield et al (2005)
  station_coord(1,:)  = (/   3.818303_8,  98.744346_8 /) !  Belawan, Indonesia
  station_coord(2,:)  = (/  -8.671252_8, 115.991079_8/) !  Lembar, Indonesia
  station_coord(3,:)  = (/  -8.475624_8, 111.765139_8 /) !  Prigi, Indonesia
  station_coord(4,:)  = (/   1.654820_8,  98.728099_8 /) !  Sibolga, Indonesia
  station_coord(5,:)  = (/  -12.113311_8,  96.953319_8 /) ! Cocos Islands
  station_coord(6,:)  = (/  -33.940205_8, 121.957036_8 /) ! Esperance, Australia
  station_coord(7,:)  = (/  -31.820068_8, 115.648509_8 /) ! Hillarys, Australia
  station_coord(8,:)  = (/ -38.434221_8, 141.687577_8 /) ! Portland, Australia
  station_coord(9,:)  = (/ -66.638126_8, 140.051264_8 /) ! Dumont d’Urville, France
  station_coord(10,:) = (/   6.924087_8,  79.787262_8 /) ! Colombo, Sri Lanka
  station_coord(11,:) = (/  -0.721699_8,  73.175345_8 /) ! Gan, Maldives
  station_coord(12,:) = (/   6.734497_8,  73.233375_8 /) ! Hanimaadhoo, Maldives
  station_coord(13,:) = (/   4.167876_8,  73.576413_8 /) ! Hulule, Male, Maldives
  station_coord(14,:) = (/  -7.359737_8,  72.537204_8 /) ! Diego Garcia, UK
  station_coord(15,:) = (/ -20.109957_8,  57.423372_8 /) ! Port Louis, Mauritius
  station_coord(16,:) = (/  -4.649683_8,  55.585582_8 /) ! Point La Rue, Seychelles 
  station_coord(17,:) = (/  16.954427_8,  54.153624_8 /) ! Salalah, Oman
  station_coord(18,:) = (/  -2.357886_8,  40.960717_8 /) ! Lamu, Kenya
  station_coord(19,:) = (/  -6.081941_8,  39.659947_8 /) ! Zanzibar, Tanzania
  station_coord(20,:) = (/ -28.865364_8,  32.137647_8 /) ! Richard's Bay, South Africa
  station_coord(21,:) = (/ -33.934105_8,  25.710086_8 /) ! Port Elizabeth, South Africa
  station_coord(22,:) = (/   3.960671_8,  93.461125_8 /) ! East of Banda Aceh

  ! Convert to radians for sph2cart function
  station_coord = station_coord * MATH_PI/180_8

  if (.not. penalize) then
      const_bathymetry = .True.
      if (rank .eq. 0) write (*,*) 'running without bathymetry and continents'
  else
      ieta = 1.0_8/eta
      alpha_m1 = alpha - 1.0_8
  end if

  if (penalize) then
      if (rank .eq. 0) write(*,*) 'Reading bathymetry data'
      allocate(bathy_data(-180*BATHY_PER_DEG:180*BATHY_PER_DEG,-90*BATHY_PER_DEG:90*BATHY_PER_DEG))
      open(unit=1086,file='bathymetry')
      do k = ubound(bathy_data,2), lbound(bathy_data,2), -1 ! north to south (as read from file)
          read(1086,*) bathy_data(:,k)
      end do
      close(1086)
  end if
  call initialize(apply_initial_conditions, 1, set_thresholds, tsunami_dump, tsunami_load)
  if (allocated(okada_data)) deallocate(okada_data)

  call init_Float_Field(wave_h, S_HEIGHT)
  call init_Float_Field(arrival, S_HEIGHT)
  do d = 1, n_domain(rank+1)
      call init(wave_h%data(d), grid(d)%node%length)
      call init(arrival%data(d), grid(d)%node%length)
      wave_h%data(d)%elts = 0.0_8
      arrival%data(d)%elts = 1.0e8_8
  end do

  if (rank .eq. 0) write (*,*) 'thresholds p, u:',  toll_height, toll_velo
  allocate(n_patch_old(size(grid)), n_node_old(size(grid)))
  n_patch_old = 2; call finish_new_patches(); call barrier()

  if (rank .eq. 0) write(*,*) 'Write initial values and grid'
  if (write_init) call write_and_export(iwrite)

  if (istep .eq. 0) then
      if (rank .eq. 0) open(unit=1011,file='verlauf.out',form='formatted',STATUS='replace')
      if (rank .eq. 0) open(unit=8450,file='tide_gauges',form='formatted',STATUS='replace')
  else
      if (rank .eq. 0) then
          open(unit=1011,file='verlauf.out',form='formatted',STATUS='old', POSITION='append')
          open(unit=8450,file='tide_gauges',form='formatted',STATUS='old', POSITION='append')
      end if
  end if
  do while (time .lt. time_end)
      call start_timing()
      call update_bdry(sol(S_HEIGHT), NONE)
      max_height = 0
      do l = level_start, level_end
          call apply_onescale(cpt_max_dh, l, 0, 1)
      end do
      max_height = sync_max_d(max_height)
      VELO_SCALE = max(VELO_SCALE*0.99, min(VELO_SCALE, grav_accel * max_height / c_p))
      call set_thresholds()
      n_patch_old = grid(:)%patch%length
      n_node_old = grid(:)%node%length
      call time_step(dt_write, aligned)
      call finish_new_patches()
      do d = 1, size(grid)
          num = grid(d)%node%length - n_node_old(d)
          if (num .gt. 0) then
              call extend(wave_h%data(d), num, 0.0_8)
              call extend(arrival%data(d), num, 1.0e8_8)
          end if
      end do
      call apply_onescale(wave_height, 9, 0, 0)
      call apply_onescale(first_arrival, 9, 0, 0)

      if (calc_tide_gauges) then
        tide_record(1:n_gauge) = -1E16_8
        do j_gauge = 1, n_gauge
           last_gauge_dist = 1E16_8
           cur_gauge = 1E3_8
           gauge_coord = project_on_sphere(sph2cart(station_coord(j_gauge,2), station_coord(j_gauge,1)))
           call apply_onescale(tide_gauge, 9, 0, 0)
           glo_gauge_dist = - sync_max_d( -last_gauge_dist ) 
           if (glo_gauge_dist .eq. last_gauge_dist) tide_record(j_gauge) = cur_gauge 
           tide_record(j_gauge) = sync_max_d(tide_record(j_gauge))
        end do
        if (rank .eq. 0) write(8450, '(50(ES13.5,1x))') time/3600.0_8*Tdim, tide_record(1:n_gauge)
     end if

      call stop_timing()
      call write_and_print_step()
      if (rank .eq. 0) write(*,'(A,F9.5,A,F9.5,2(A,E13.5),A,I9)') &
              'time [h] =', time/3600.0_8*Tdim, &
              ', dt [s] =', dt*Tdim, &
              ', min. depth =', fd, &
              ', U =', VELO_SCALE, &
              ', d.o.f. =', sum(n_active)
      call print_load_balance()
      if (aligned) then
          iwrite = iwrite + 1
          call write_and_export(iwrite)
          if (modulo(iwrite,CP_EVERY) .ne. 0) cycle
          call update_bdry(arrival, 9)
          call update_bdry(wave_h, 9)
          call export_2d(cart2sph2, (/arrival, wave_h/), 2, 10000+10*iwrite/CP_EVERY, 9, & 
                  (/-768, 768/), (/-384, 384/), (/2.0_8*MATH_PI, MATH_PI/), (/1.0e8_8, 0.0_8/))
          ierr = writ_checkpoint(tsunami_dump)

          ! let all cpus exit gracefully if NaN has been produced
          ierr = sync_max(ierr)
          if (ierr .eq. 1) then ! NaN
              write(0,*) "NaN when writing checkpoint"
              call finalize()
              stop
          end if

          do d = 1, n_domain(rank+1)
              deallocate(wave_h%data(d)%elts)
              deallocate(arrival%data(d)%elts)
          end do
          deallocate(wave_h%data)
          deallocate(arrival%data)
          call restart_full(set_thresholds, tsunami_load)
          deallocate(n_patch_old); allocate(n_patch_old(size(grid)))
          deallocate(n_node_old);  allocate(n_node_old(size(grid)))
          call init_Float_Field(wave_h, S_HEIGHT)
          call init_Float_Field(arrival, S_HEIGHT)
          do d = 1, n_domain(rank+1)
              call init(wave_h%data(d), grid(d)%node%length)
              call init(arrival%data(d), grid(d)%node%length)
              wave_h%data(d)%elts = 0.0_8
              arrival%data(d)%elts = 1.0e8_8
          end do
          n_patch_old = 2; call finish_new_patches()
          ! finish_new_patches takes long time (the smoothing of penalization)
          ! barrier here so that this does not affect following timing
          call barrier()
      end if
  end do
  if (rank .eq. 0) then
     close(1011)
     close(8450)
  end if
  call finalize()
end program
