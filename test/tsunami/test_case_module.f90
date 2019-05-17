Module test_case_mod
  ! Module file for DCMIP2008c5
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  ! Standard variables
  integer                              :: CP_EVERY, resume_init, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Local variables
  integer                              :: bathy_per_deg, npts_chi, npts_topo
  integer, dimension(:), allocatable   :: n_patch_old
  real(8)                              :: dH, lon_c, lat_c, min_depth, pert_radius
  real(4), allocatable, dimension(:,:) :: bathy_data
  logical                              :: const_bathymetry
contains
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id_i
    real (8)    :: dz, porosity, porous_density, z
    type(Coord) :: x_i
    
    d    = dom%id+1
    id_i = idx (i, j, offs, dims) + 1
    x_i  = dom%node%elts(id_i)

    ! Initial vertical grid size (uniform)
    dz = (free_surface (x_i) - dom%topo%elts(id_i)) / zlevels
 
    ! Local z-coordinate (mid layer)
    z = free_surface (x_i) - (zlevels - zlev + 0.5_8) * dz 

    ! Penalization
    porosity = 1.0_8 - (1.0_8-alpha) * penal(zlev)%data(d)%elts(id_i)
    porous_density = ref_density * porosity
    
    ! Equally spaced sigma coordinates in z: sol(S_MASS) = ref_density * dz
    sol(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * dz

    ! Density * dz
    sol(S_TEMP,zlev)%data(d)%elts(id_i) = sol(S_MASS,zlev)%data(d)%elts(id_i) * density (x_i, z)/porous_density

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  real(8) function surf_geopot (x_i)
    ! Surface geopotential: postive if greater than mean seafloor
    implicit none
    type(Coord) :: x_i

    surf_geopot = grav_accel * 0.0_8
  end function surf_geopot

  real(8) function free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    real(8) :: lon, lat, rgrc
    
    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    rgrc = radius * acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

    free_surface = dH * exp__flush (-(rgrc/pert_radius)**2)
  end function free_surface

  real(8) function density (x_i, z)
    ! Density profile
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    real(8), parameter :: drho = 0.05_8

    ! density = ref_density
    
    if (z > 0.8 * mean_depth) then ! less dense fluid in upper layer
       density = ref_density * (1.0_8 - drho)
    else
       density = ref_density
    end if
  end function density

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    u = 0.0_8 ! Zonal velocity component
    v = 0.0_8 ! Meridional velocity component
  end subroutine vel_fun
  
  subroutine set_thresholds
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                     :: k
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: threshold_new
    character(3), parameter                     :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       if (adapt_trend) then
          call cal_lnorm_trend (trend, order)
       else
          call cal_lnorm_sol (sol, order)
       end if
       threshold_new = tol*lnorm
       ! Correct for zero velocity case
       do k = 1, zlevels
          if (threshold_new(S_VELO,k) == 0.0_8) threshold_new(S_VELO,k) = 1d16 
       end do
    end if

    if (istep >= 10) then
       threshold = 0.01*threshold_new + 0.99*threshold
    else
       threshold = threshold_new
    end if
  end subroutine set_thresholds
  
  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    real(8) :: dz

    allocate (threshold(S_MASS:S_VELO,1:zlevels));     threshold     = 0.0_8
    allocate (threshold_def(S_MASS:S_VELO,1:zlevels)); threshold_def = 0.0_8

    dz = Hdim/zlevels
    
    lnorm(S_MASS,:) = ref_density * dz
    lnorm(S_TEMP,:) = ref_density * dz
    lnorm(S_VELO,:) = Udim
    if (adapt_trend) lnorm = lnorm/Tdim
    threshold_def = tol * lnorm
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, tau_divu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle
      
    ! Diffusion constants
    C_sclr = 5d-3    ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
    C_divu = 5d-3    ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
    C_rotu = C_sclr / 4**Laplace_order_init ! <= 1.09e-3 for hyperdiffusion (lower than exact limit 1/24^2 = 1.7e-3 due to non-uniform grid)
    
    ! CFL limit for time step
    dt_cfl = cfl_num*dx_min/(wave_speed+Udim) * 0.85 ! corrected for dynamic value
    dt_init = dt_cfl

    tau_sclr = dt_cfl / C_sclr
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,3(a,es8.2),a,/)') "dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2))') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
          " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
          " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity

  subroutine set_depth_penal (dom, i, j, zlev, offs, dims)
    ! Sets depth and penalization mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i, k
    logical, parameter :: etopo = .false.

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    ! Set topography
    if (etopo) then ! set bathymetry coordinates using etopo data
       call etopo_topography (id, dom, dom%topo%elts(id_i), 'bathymetry')
    else
       call analytic_topography (dom%node%elts(id_i), dom%topo%elts(id_i), 'bathymetry')
    end if

    ! Set penalization
    if (penalize) then
       if (etopo) then
          call etopo_topography (id, dom, dom%topo%elts(id_i), 'penalize')
       else
          call analytic_topography (dom%node%elts(id_i), penal(1)%data(d)%elts(id_i), 'penalize')
       end if
    else
       penal(1)%data(d)%elts(id_i) = 0.0_8       
    end if

    ! For now, penalization is the same at all vertical levels (only used to define land areas)
    do k = 2, zlevels
       penal(k)%data(dom%id+1)%elts(id_i) = penal(1)%data(dom%id+1)%elts(id_i)
    end do
  end subroutine set_depth_penal

  subroutine analytic_topography (x_i, mask, itype)
    implicit none
    real(8)      :: mask
    character(*) :: itype
    type(Coord)  :: x_i

    real(8)            :: lon, lat, rgrc
    real(8), parameter :: land_radius = 4000*KM

    if (itype == 'bathymetry') then
       mask = mean_depth + surf_geopot (x_i) / grav_accel
    else ! penalization
       call cart2sph (x_i, lon, lat)
       rgrc = radius * acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))
       if (dist (x_i, sph2cart (lon, lat)) < land_radius) then
          mask = 1.0_8
       else
          mask = 0.0_8
       end if
    end if
  end subroutine analytic_topography

  subroutine etopo_topography (id, dom, topo, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo using etopo data 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    integer      :: id
    real(8)      :: penal, topo
    character(*) :: itype
    type(Domain) :: dom

    integer               :: e, i, id_e, id_i, j,  is0, it0, k, s, t
    real(8)               :: dx_local, dx_smooth, lat, lon, M_topo, r, s0, t0, sw_topo, topo_sum, wgt
    real(8), dimension(3) :: dx_primal, dx_dual
    type(Coord)           :: p, q

    id_i = id+1

    ! Determine effective grid size
    do e = 1, EDGE
       id_e = id*EDGE + e
       dx_primal(e) = dom%len%elts(id_e)
       dx_dual(e)   = dom%pedlen%elts(id_e)
    end do
    dx_min = sqrt (4/sqrt(3.0_8) * 4*MATH_PI*radius**2/(20*4**max_level))
    dx_local = max (maxval(dx_primal), maxval(dx_dual))
    if (dx_local == 0.0_8) dx_local = dx_min

    ! Find longitude and latitude coordinates
    call cart2sph (dom%node%elts(id_i), lon, lat)
    s0 = lon/MATH_PI * dble (180*BATHY_PER_DEG)
    t0 = lat/MATH_PI * dble (180*BATHY_PER_DEG)
    p = proj_lon_lat (s0, t0)
    is0 = nint (s0); it0 = nint (t0)

    if (npts_topo == 0) then ! no smoothing
       topo = topo_value ()
    else ! smoothing
       sw_topo  = 0.0_8
       topo_sum = 0.0_8
       do i = -npts_topo, npts_topo
          do j = -npts_topo, npts_topo
             s = is0+i ; t = it0+j

             call wrap_lonlat (s, t)
             q = proj_lon_lat (dble(s), dble(t))
             r = norm (vector(p, q))
             wgt = radial_basis_fun (r, npts_topo, dx_smooth)

             M_topo = topo_value ()
             topo_sum = topo_sum + wgt * M_topo
             sw_topo  = sw_topo  + wgt
          end do
       end do
       topo = topo_sum / sw_topo
    end if
  contains
    real(8) function topo_value ()
      implicit none

      if (bathy_data(s, t) > 0.0_8) then ! land
         if (itype == 'bathymetry') then ! topography
            topo_value = 0.0_8
         else ! penalization
            topo_value = LAND
         end if
      else ! sea: bathymetry is less than zero
         if (itype == 'topography') then ! topography
            topo_value = min (bathy_data(s, t), min_depth)
         else ! penalization
            topo_value = SEA
         end if
      end if
    end function topo_value
  end subroutine etopo_topography

  real(8) function radial_basis_fun (r, npts, dx_local)
    ! Radial basis function for smoothing topography
    implicit none
    integer :: npts
    real(8) :: r, dx_local

    real(8) :: alph

    alph = 1.0_8 / (npts/2 * dx_local)

    radial_basis_fun = exp (-(alph*r)**2)
  end function radial_basis_fun

  subroutine wrap_lonlat (s, t)
    ! Longitude: wraparound allows for values outside [-180,180]
    ! Latitude: works only if there is no cost at the pole
    implicit none
    integer :: s, t
    
    if (t < lbound (bathy_data,2)) t = lbound (bathy_data,2) ! pole
    if (t > ubound (bathy_data,2)) t = ubound (bathy_data,2) ! pole
    if (s < lbound (bathy_data,1)) s = s + 360*BATHY_PER_DEG
    if (s > ubound (bathy_data,1)) s = s - 360*BATHY_PER_DEG
  end subroutine wrap_lonlat

  type(Coord) function proj_lon_lat (s, t)
    implicit none
    real(8) :: s, t
    real(8) :: lon, lat
    
    lon = s * MATH_PI / dble (180*BATHY_PER_DEG)
    lat = t * MATH_PI / dble (180*BATHY_PER_DEG)
    proj_lon_lat = project_on_sphere (sph2cart(lon, lat))
  end function proj_lon_lat 

  subroutine read_test_case_parameters
    implicit none
    integer            :: k
    integer, parameter :: fid = 500
    real(8)            :: press_save
    character(255)     :: filename, varname

    ! Find input parameters file name
    if (command_argument_count() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'test_case.in'
    end if
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

    open(unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, compressible
    read (fid,*) varname, penalize
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, remap
    read (fid,*) varname, remapscalar_type
    read (fid,*) varname, remapvelo_type
    read (fid,*) varname, iremap
    read (fid,*) varname, adapt_trend
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, perfect
    read (fid,*) varname, tol
    read (fid,*) varname, optimize_grid
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, cfl_num
    read (fid,*) varname, timeint_type
    read (fid,*) varname, press_save
    read (fid,*) varname, Laplace_order_init
    read (fid,*) varname, n_diffuse
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, rebalance
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    read (fid,*) varname, eta
    read (fid,*) varname, alpha
    read (fid,*) varname, min_depth
    read (fid,*) varname, bathy_per_deg
    read (fid,*) varname, npts_chi
    read (fid,*) varname, npts_topo
    close(fid)
    
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    call set_save_level
    dt_write = dt_write * MINUTE
    time_end = time_end * HOUR
    resume   = resume_init
    Laplace_order = Laplace_order_init

    if (penalize) then
       const_bathymetry = .false.
       if (rank == 0) write(6,'(a)') 'Reading bathymetry data'
       allocate (bathy_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))
       open (unit=1086,file='bathymetry')
       do k = ubound (bathy_data,2), lbound (bathy_data,2), -1 ! north to south (as read from file)
          read (1086,*) bathy_data(:,k)
       end do
       close (1086)
    else
       const_bathymetry = .true.
    end if
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none
    
     if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case            = ", trim (test_case)
       write (6,'(A,A)')      "run_id               = ", trim (run_id)
       write (6,'(A,L1)')     "compressible         = ", compressible
       write (6,'(A,L1)')     "penalize             = ", penalize
       write (6,'(A,i3)')     "min_level            = ", min_level
       write (6,'(A,i3)')     "max_level            = ", max_level
       write (6,'(A,i5)')     "number of domains    = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL         = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL          = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels              = ", zlevels
       write (6,'(A,L1)')     "remap                = ", remap
       write (6,'(a,a)')      "remapscalar_type     = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type       = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap               = ", iremap
       write (6,'(A,L1)')     "adapt_trend          = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds   = ", default_thresholds
       write (6,'(A,L1)')     "perfect              = ", perfect
       write (6,'(A,es10.4)') "tolerance            = ", tol
       write (6,'(A,i1)')     "optimize_grid        = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt             = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num              = ", cfl_num
       write (6,'(a,a)')      "timeint_type         = ", trim (timeint_type)
       write (6,'(A,es10.4)') "pressure_save [hPa]  = ", pressure_save(1)/100
       write (6,'(A,i1)')     "Laplace_order        = ", Laplace_order_init
       write (6,'(A,i4)')     "n_diffuse            = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write [min]       = ", dt_write/MINUTE
       write (6,'(A,i6)')     "CP_EVERY             = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance            = ", rebalance
       write (6,'(A,es10.4)') "time_end [h]         = ", time_end/HOUR
       write (6,'(A,i6)')     "resume               = ", resume_init
        
       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius               = ", radius
       write (6,'(A,es10.4)') "omega                = ", omega
       write (6,'(A,es10.4)') "p_top (hPa)          = ", p_top/100

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es10.4)') "pert_radius          = ", pert_radius
       write (6,'(A,es10.4)') "dH                   = ", dH
       write (6,'(A,es11.4)') "mean_depth           = ", mean_depth
       write (6,'(A,es10.4)') "lon_c                = ", lon_c
       write (6,'(A,es10.4)') "lat_c                = ", lat_c
       write (6,'(A,es10.4)') "eta                  = ", eta
       write (6,'(A,es10.4)') "alpha                = ", alpha
       write (6,'(A,es10.4)') "min_depth            = ", min_depth
       write (6,'(A,i5)')     "bathy_per_deg        = ", bathy_per_deg
       write (6,'(A,i5)')     "npts_chi             = ", npts_chi
       write (6,'(A,i5)')     "npts_topo            = ", npts_topo

       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    
    if (rank == 0) then
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i9,4(a,es9.2,1x))') &
            'time [h] = ', time/HOUR, &
            ' dt [s] = ', dt, &
            '  mass tol = ', sum (threshold(S_MASS,:))/zlevels, &
            ' temp tol = ', sum (threshold(S_TEMP,:))/zlevels, &
            ' velo tol = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(5(es15.9,1x),i2,1x,i9,1x,4(es15.9,1x))')  &
            time/HOUR, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine apply_initial_conditions
    implicit none
    integer :: d, k, l, p

    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (set_depth_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine update_depth_penalization
    implicit none
    integer :: d, k, p

    do d = 1, size(grid)
       do k = 1, zlevels
          do p = n_patch_old(d)+1, grid(d)%patch%length
             call apply_onescale_to_patch (set_depth_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do
  end subroutine update_depth_penalization
  
  subroutine vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
    ! Sets the velocities on the computational grid given a function vel_fun that provides zonal and meridional velocities
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    external                        :: vel_fun

    integer      :: d, id, idE, idN, idNE
    type (Coord) :: x_i, x_E, x_N, x_NE

    d = dom%id+1

    id   = idx (i,   j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel (vel_fun, x_i,  x_E)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel (vel_fun, x_NE, x_i)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel (vel_fun, x_i,  x_N)
  end subroutine vel2uvw

  subroutine set_save_level
    ! Determines closest vertical level to desired height z
    implicit none
    integer :: k
    real(8) :: dz, p, save_press, z

    dz = abs (mean_depth) / zlevels; save_zlev = 0
    do k = 1, zlevels
       z = (zlevels - k) * dz 
       if (abs (z - pressure_save(1)) < dz) then
          save_zlev = k
          save_press = z
          exit
       end if
    end do
    if (rank==0) write (6,'(/,A,i2,A,es10.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_press, " [m])"
  end subroutine set_save_level
  
  subroutine initialize_a_b_vert
    implicit none

    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
  end subroutine initialize_a_b_vert

  subroutine dump (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load
end module test_case_mod
