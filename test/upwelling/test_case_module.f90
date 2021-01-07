Module test_case_mod
  ! Module file for Drake passage test case
  use shared_mod
  use comm_mpi_mod
  use utils_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, resume_init, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Local variables
  real(8)                              :: beta, beta0, drho, f0, Rd, ref_temp
  real(8)                              :: bottom_friction, max_depth, r_max, r_max_loc
  real(8)                              :: npts_penal, u_wbc
  real(8)                              :: lat_c, lat_width, tau_0, width
  real(4), allocatable, dimension(:,:) :: topo_data
  character(255)                       :: coords
contains
  subroutine read_test_case_parameters
    implicit none
    integer            :: ilat, ilon, k
    integer, parameter :: fid = 500
    real(8)            :: lat, lon, press_save
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
    read (fid,*) varname, max_level
    read (fid,*) varname, remap
    read (fid,*) varname, iremap
    read (fid,*) varname, log_iter
    read (fid,*) varname, tol
    read (fid,*) varname, cfl_num
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)

    press_save = 0.0_8
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    call set_save_level

    call cal_r_max

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case                      = ", trim (test_case)
       write (6,'(A,A)')      "run_id                         = ", trim (run_id)
       write (6,'(A,L1)')     "compressible                   = ", compressible
       write (6,'(A,L1)')     "mode_split                     = ", mode_split
       write (6,'(A,L1)')     "penalize                       = ", penalize
       write (6,'(A,es10.4)') "npts_penal                     = ", npts_penal
       write (6,'(A,i3)')     "min_level                      = ", min_level
       write (6,'(A,i3)')     "max_level                      = ", max_level
       write (6,'(A,i3)')     "level_fill                     = ", level_fill
       write (6,'(A,i5)')     "number of domains              = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors           = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL                   = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL                    = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels                        = ", zlevels
       write (6,'(A,L1)')     "remap                          = ", remap
       write (6,'(a,a)')      "remapscalar_type               = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type                 = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap                         = ", iremap
       write (6,'(A,es10.4)') "tol_elliptic                   = ", tol_elliptic
       write (6,'(a,i3)')     "coarse_iter                    = ", coarse_iter
       write (6,'(a,i3)')     "fine_iter                      = ", fine_iter
       write (6,'(a,l1)')     "log_iter                       = ", log_iter
       write (6,'(A,L1)')     "adapt_trend                    = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds             = ", default_thresholds
       write (6,'(A,L1)')     "perfect                        = ", perfect
       write (6,'(A,es10.4)') "tolerance                      = ", tol
       write (6,'(A,i1)')     "optimize_grid                  = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt                       = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num                        = ", cfl_num
       write (6,'(a,a)')      "timeint_type                   = ", trim (timeint_type)
       write (6,'(A,i1)')     "Laplace_order                  = ", Laplace_order_init
       write (6,'(A,i1)')     "n_diffuse                      = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write [d]                   = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY                       = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance                      = ", rebalance
       write (6,'(A,es10.4)') "time_end [d]                   = ", time_end/DAY
       write (6,'(A,i6)')     "resume                         = ", resume_init

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius                   [km]  = ", radius / KM
       write (6,'(A,es10.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", abs (max_depth)
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/bottom_friction / DAY
       write (6,'(A,es11.4)') "f0 at 45 deg          [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 45 deg       [rad/ms]  = ", beta
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(A,es11.4)') "barotropic Rossby radius [km]  = ", Rd / KM
       write (6,'(a,es11.4)') "r_max                          = ", r_max
       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'

       call print_density
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
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i12,4(a,es9.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', threshold(S_MASS,zlevels), &
            ' temp tol = ', threshold(S_TEMP,zlevels), &
            ' velo tol = ', threshold(S_VELO,zlevels), &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(5(es15.9,1x),i2,1x,i12,1x,4(es15.9,1x))')  time/DAY, dt, &
            threshold(S_MASS,zlevels), threshold(S_TEMP,zlevels), threshold(S_VELO,zlevels), &
            level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine apply_initial_conditions
    implicit none
    integer :: d, k, l

    do l = level_start, level_end
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_start, level_end
       do k = 1, zmax
          call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          call apply_onescale (init_sol,  l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initial perturbation to mean 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id, id_i
    real (8)    :: dz, eta, phi, z_s
    type(Coord) :: p

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    p = dom%node%elts(id_i)
    eta = init_free_surface (p)
    z_s = dom%topo%elts(id_i)

    if (zlev == zlevels+1) then ! 2D barotropic mode
       phi = 1.0_8 + (alpha - 1.0_8) * penal_node(zlevels)%data(d)%elts(id_i)
       sol(S_MASS,zlev)%data(d)%elts(id_i) = phi * eta ! free surface perturbation
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
    else ! 3D layers
       dz = a_vert_mass(zlev) * eta + b_vert_mass(zlev) * z_s
       if (zlev == zlevels) then
          sol(S_MASS,zlev)%data(d)%elts(id_i) = porous_density (dom, i, j, zlev, offs, dims) * eta
       else
          sol(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
       end if
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
    end if
    ! Set initial velocity field to zero
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id:EDGE*id_i) = 0.0_8
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id, id_i
    real (8)    :: dz, eta, z_s
    type(Coord) :: p

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    p = dom%node%elts(id_i)

    eta = init_free_surface (p)
    z_s = dom%topo%elts(id_i)

    if (zlev == zlevels+1) then
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
    else
       dz = a_vert_mass(zlev) * eta + b_vert_mass(zlev) * z_s
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = porous_density (dom, i, j, zlev, offs, dims) * dz
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy (eta, z_s, zlev)
    end if
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
  end subroutine init_mean

  subroutine update
    ! Update means, bathymetry and penalization mask
    implicit none
    integer :: d, k, p

    do d = 1, size(grid)
       do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches
          call apply_onescale_to_patch (set_bathymetry, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          do k = 1, zmax
             call apply_onescale_to_patch (set_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do

    do k = 1, zmax
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches
             call apply_onescale_to_patch (init_mean, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do
  end subroutine update

  real(8) function surf_geopot (p)
    ! Surface geopotential: postive if greater than mean seafloor
    ! Gives minimum depth of approximately 20 m.
    implicit none
    type(Coord) :: p

    real(8)            :: lat, lon, y
    real(8), parameter :: d_min = -14.3, b1 = 0.38, b2 = 0.7

    call cart2sph (p, lon, lat)

    lat = lat / DEG
    if (abs(lat-lat_c) <= lat_width/2) then
       y = (lat - (lat_c - lat_width/2))/180 * MATH_PI*radius
       surf_geopot = 65.5_8 - 66.526 * tanh (1.5d-4 * (f(y) - width/8))
    else
       surf_geopot = 65.5_8 - 66.526 * tanh (-1.875d-5*width)
    end if
    surf_geopot = grav_accel * surf_geopot
  end function surf_geopot

  real(8) function f (y)
    implicit none
    real(8) :: y

    if (y <= width/2) then
       f = y
    else
       f = width - y
    end if
  end function f

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0.0_8
  end function init_free_surface

   real(8) function buoyancy (eta, z_s, zlev)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    integer :: zlev
    real(8) :: eta, z_s

    real(8) :: rho, z1, z2

    z1 = a_vert(zlev-1) * eta + b_vert(zlev-1) * z_s
    z2 = a_vert(zlev)   * eta + b_vert(zlev)   * z_s

    rho = 0.5 * (density (z1) + density (z2))
    
    buoyancy = (ref_density - rho) / ref_density 
  end function buoyancy

  real(8) function density (z)
    implicit none
    real(8) :: z

    density = eqn_of_state (temp_profile(z))
  end function density

  real(8) function temp_profile (z)
    implicit none
    real(8) :: z
    
    real(8), parameter :: h_z = 6.5_8, T0 = 14_8, strat = 150_8, z0 = -35_8, z1 = -75_8

    temp_profile = T0 + 4*tanh ((z - z0) / h_z) + (z - z1) / strat
  end function temp_profile

  real(8) function eqn_of_state (temperature)
    implicit none
    real(8) :: temperature

    real(8), parameter :: beta = 0.28, T0 = 14_8

    eqn_of_state = ref_density - beta * (temperature - T0)
  end function eqn_of_state

  subroutine print_density
    implicit none
    integer     :: k
    real(8)     :: bv, c1, drho, dz, eta, rho, rho_above, z, z_s, z_above
    type(Coord) :: p

    p = Coord (radius, 0.0_8, 0.0_8)

    eta = 0.0_8
    z_s = max_depth

    write (6,'(a)') " Layer    z        dz            rho "
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta + b_vert_mass(k) * z_s
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta + (b_vert(k)+b_vert(k-1)) * z_s)
       write (6, '(2x, i2, 1x, 3(es11.4, 1x))') k, z, dz, ref_density * (1.0_8 - buoyancy (eta, z_s, k))
    end do
    
    write (6,'(/,a)') " Interface     a_vert      b_vert        z"
    do k = 0, zlevels
       write (6, '(3x, i3, 5x, 3(1x, es11.4))') k, a_vert(k), b_vert(k), a_vert(k) * eta + b_vert(k) * z_s
    end do

    write (6,'(/,a)') " Interface    V        c1     CFL_c1"
    do k = 1, zlevels-1
       z_above = 0.5 * ((a_vert(k+1)+a_vert(k)) * eta + (b_vert(k+1)+b_vert(k)) * z_s)
       z       = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta + (b_vert(k)+b_vert(k-1)) * z_s)
       dz = z_above - z
       rho_above = ref_density * (1.0_8 - buoyancy (eta, z_s, k+1))
       rho  = ref_density * (1.0_8 - buoyancy (eta, z_s, k))
       drho = rho_above - rho
       bv = sqrt(- grav_accel * drho/dz/rho)
       c1 = bv * abs(max_depth) / MATH_PI
       write (6, '(3x, i3, 5x,3(es8.2,1x))') k, bv, c1, c1*dt_init/dx_min
    end do
    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density

  subroutine set_thresholds
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                 :: k, v
    real(8), dimension(1:N_VARIABLE,1:zmax) :: threshold_new
    character(3), parameter                 :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       if (adapt_trend) then
          call cal_lnorm_trend (trend, order)
       else
          call cal_lnorm_sol (sol, order)
       end if
       threshold_new = tol * lnorm
       ! Correct very small values
       do k = 1, zmax
          if (threshold_new(S_MASS,k) < threshold_def(S_MASS,k)/10) threshold_new(S_MASS,k) = threshold_def(S_MASS,k)
          if (threshold_new(S_TEMP,k) < threshold_def(S_TEMP,k)/10) threshold_new(S_TEMP,k) = threshold_def(S_TEMP,k)
          if (threshold_new(S_VELO,k) < threshold_def(S_VELO,k)/10) threshold_new(S_VELO,k) = threshold_def(S_VELO,k)
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
    integer     :: k
    real(8)     :: dz, eta, z_s
    type(Coord) :: x_i

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0.0_8
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0.0_8

    x_i = Coord (radius, 0.0_8, 0.0_8)
    eta = 0.0_8
    z_s = max_depth
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta + b_vert_mass(k) * z_s
       lnorm(S_MASS,k) = ref_density * dz

!!$       lnorm(S_TEMP,k) = ref_density*dz * buoyancy (max_depth, x_i, k)
       lnorm(S_TEMP,k) = 1d16

       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8)            :: area, C_divu, C_sclr, C_rotu, C_visc, tau_divu, tau_rotu, tau_sclr
    logical, parameter :: munk = .true.

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

    area = 4*MATH_PI*radius**2/(20*4**min_level)
    dx_max = sqrt (4/sqrt(3.0_8) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, 1.4*dx_min/u_wbc, dx_min/c1)
    dt_init = dt_cfl

    C_rotu = 1d-3
    C_divu = 1d-3
    C_sclr = 1d-3

    ! Diffusion time scales
    tau_sclr = dt_cfl / C_sclr
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr(S_MASS) = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_sclr(S_TEMP) = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,4(a,es8.2),a,/)') &
            "dx_max  = ", dx_max/KM, " dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topography (dom, i, j, zlev, offs, dims, 'bathymetry')
  end subroutine set_bathymetry

  subroutine set_penal (dom, i, j, zlev, offs, dims)
    ! Set penalization mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (penalize) then
       call topography (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
    end if
  end subroutine set_penal

  subroutine topography (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer     :: d, e, id, id_e, id_i, idE, idN, idNE
    real(8)     :: lat, lon, mask, n_smth_N, n_smth_S, width_N, width_S
    type(Coord) :: p

    id = idx (i, j, offs, dims)
    id_i = id + 1

    p = dom%node%elts(id_i)

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth + surf_geopot (p) / grav_accel
    case ("penalize") ! analytic land mass with smoothing
       d = dom%id + 1
       call cart2sph (dom%node%elts(id_i), lon, lat)

       width_S = lat_c - lat_width/2 + 90_8
       width_N = lat_c - lat_width/2

       ! Smoothing exponent for land mass
       n_smth_S = 4*radius * width_S*DEG / (dx_max * npts_penal)
       n_smth_N = 4*radius * width_N*DEG / (dx_max * npts_penal)
       
       mask = exp__flush (- abs((lat/DEG+90_8)/width_S)**n_smth_S) + exp__flush (- abs((lat/DEG-90_8)/width_N)**n_smth_N)
       penal_node(zlev)%data(d)%elts(id_i) = mask
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = max (penal_edge(zlev)%data(d)%elts(id_e), mask)
       end do
    end select
  end subroutine topography

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid
    
    if (time/DAY <= 2.0_8) then
       tau_zonal = tau_0 * sin (MATH_PI * time/DAY)
    else
       tau_zonal = tau_0
    end if
    tau_merid = 0.0_8
  end subroutine wind_stress

  subroutine set_save_level
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level

  subroutine initialize_a_b_vert
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    b_vert(0) = 1.0_8 ; b_vert(zlevels) = 0.0_8
    do k = 1, zlevels-1
       if (trim (coords) == "uniform") then 
          b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
       elseif (trim (coords) == "chebyshev") then
          b_vert(k) = (1.0_8 + cos (dble(2*k-1)/dble(2*(zlevels-1)) * MATH_PI)) / 2
       elseif (trim(coords) == "roms") then
          b_vert(0) = -150
          b_vert(1) = -103.935
          b_vert(2) =  -73.66
          b_vert(3) =  -53.57
          b_vert(4) =  -40.06
          b_vert(5) =  -30.80
          b_vert(6) =  -24.28
          b_vert(7) =  -19.54
          b_vert(8) =  -15.94
          b_vert(9) =  -13.07
          b_vert(10) = -10.68 
          b_vert(11) =  -8.60
          b_vert(12) =  -6.71
          b_vert(13) =  -4.95
          b_vert(14) =  -3.26
          b_vert(15) =  -1.62
          b_vert(16) =   0
          b_vert = b_vert/max_depth
       end if
    end do
    a_vert = 1.0_8 - b_vert
       
    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert

  subroutine update_diagnostics
    ! Update diagnostics
    implicit none
    
  end subroutine update_diagnostics

  subroutine init_diagnostics
    ! Initialize diagnostics
    implicit none 
   
  end subroutine init_diagnostics

  subroutine deallocate_diagnostics
    implicit none
    
  end subroutine deallocate_diagnostics
  
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

  subroutine cal_r_max
    ! Calculates minimum relative mass and checks diffusion stability limits
    use mpi
    implicit none
    integer :: ierror, k, l

    r_max_loc = 1d-16
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (cal_rmax_loc, l, k, 0, 0)
       end do
    end do

    call MPI_Allreduce (r_max_loc, r_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine cal_r_max

  subroutine cal_rmax_loc (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idE, idN, idNE, idS, idSW, idW, k
    real(8) :: dz0, dz_e, r_loc

    id   = idx (i,   j,   offs, dims)
    
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    
    d    = dom%id + 1
    
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       dz0  = sol_mean(S_MASS,zlev)%data(d)%elts(id+1) / porous_density (dom, i, j, zlev, offs, dims)
       
       dz_e = sol_mean(S_MASS,zlev)%data(d)%elts(idE+1) / porous_density (dom, i+1, j, zlev, offs, dims)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1) / porous_density (dom, i+1, j+1, zlev, offs, dims)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = sol_mean(S_MASS,zlev)%data(d)%elts(idN+1) / porous_density(dom, i, j+1, zlev, offs, dims)
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc

  real(8) function eddy_diffusivity (eta, ri, z)
    ! Eddy diffusivity at nodes
    implicit none
    real(8) :: eta, ri, z

    real(8), parameter :: K_t = 1d-6
    real(8), parameter :: a = 5, A_ric = 1d-4, At_b = 1.2d-5

    eddy_diffusivity = K_t
!!$    eddy_diffusivity = A_ric/(1.0_8 + a*ri**2) + At_b
  end function eddy_diffusivity

  function eddy_viscosity (eta, ri, z)
    ! Eddy viscosity at at edges
    implicit none
    real(8), dimension(1:EDGE) :: eddy_viscosity
    real(8)                    :: ri
    real(8), dimension(1:EDGE) :: eta, z

    real(8), parameter :: a = 5, Av_b = 1.2d-4, K_m = 2d-3
    
    eddy_viscosity = K_m * (1.0_8 + 4 * exp ( (z - eta) / abs(max_depth) ))
!!$    eddy_viscosity = eddy_diffusivity (eta(1), ri, z(1)) / (1.0_8 + a*ri) + Av_b
  end function eddy_viscosity

  real(8) function bottom_temp_source (dom, i, j, zlev, offs, dims)
    ! Top boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    bottom_temp_source = 0.0_8
  end function bottom_temp_source

  real(8) function top_temp_source (dom, i, j, zlev, offs, dims)
    ! Bottom boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    top_temp_source = 0.0_8
  end function top_temp_source
  
  function bottom_velo_source (dom, i, j, zlev, offs, dims)
    ! Linear bottom friction source term (bottom boundary condition for vertical diffusion of velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: bottom_velo_source

    integer                     :: id
    real(8), dimension(1:EDGE)  :: dz


    id = idx (i, j, offs, dims)

    dz = dz_e (dom, i, j, zlev, offs, dims)

    bottom_velo_source = - bottom_friction * velo(EDGE*id+RT+1:EDGE*id+UP+1) / dz
  end function bottom_velo_source

  function top_velo_source (dom, i, j, zlev, offs, dims)
    ! Wind stress velocity source term evaluated at edges (top boundary condition for vertical diffusion of velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: top_velo_source

    integer                         :: id, idE, idN, idNE               
    real(8), dimension(1:EDGE)      :: mass_e, tau_wind
    real(8), dimension(0:NORTHEAST) :: full_mass

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    full_mass(0:NORTHEAST) = mean_m((/id,idN,idE,id,id,idNE/)+1) + mass((/id,idN,idE,id,id,idNE/)+1)

    mass_e(RT+1) = 0.5 * (full_mass(0) + full_mass(EAST))
    mass_e(DG+1) = 0.5 * (full_mass(0) + full_mass(NORTHEAST))
    mass_e(UP+1) = 0.5 * (full_mass(0) + full_mass(NORTH))

    tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idE+1))
    tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE+1), dom%node%elts(id+1))
    tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idN+1))

    top_velo_source = tau_wind / mass_e
  end function top_velo_source
end module test_case_mod
