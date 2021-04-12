Module test_case_mod
  ! Module file for upwelling test case
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
  real(8)                              :: friction_upwelling, r_max, r_max_loc
  real(8)                              :: n_smth_N, n_smth_S, width_N, width_S
  real(8)                              :: npts_penal, u_wbc
  real(8)                              :: lat_c, lat_width, tau_0, Tcline, width
  real(4), allocatable, dimension(:,:) :: topo_data
  character(255)                       :: coords

  real(8), parameter :: slope = 1.75d-4 ! slope parameter (larger value -> steeper slope)
  real(8), parameter :: shift = 8
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
    read (fid,*) varname, zlevels
    read (fid,*) varname, no_slip
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
       write (6,'(A,L1)')     "no_slip                        = ", no_slip
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
       write (6,'(A,L1)')     "sigma_z                        = ", sigma_z
       write (6,'(a,a)')      "remapscalar_type               = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type                 = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap                         = ", iremap
       write (6,'(A,es10.4)') "tol_elliptic                   = ", tol_elliptic
       write (6,'(a,i3)')     "coarse_iter                    = ", coarse_iter
       write (6,'(a,i3)')     "fine_iter                      = ", fine_iter
       write (6,'(a,l1)')     "log_iter                       = ", log_iter
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
       write (6,'(A,es11.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", max_depth
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", min_depth
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", Tcline
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", friction_upwelling
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/friction_upwelling / DAY
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
            '  mass threshold = ', sum (threshold(S_MASS,:))/zlevels, &
            ' temp threshold = ', sum (threshold(S_TEMP,:))/zlevels, &
            ' velo threshold = ', sum (threshold(S_VELO,:))/zlevels, &
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
  
  real(8) function buoyancy (z, k)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    integer                       :: k
    real(8), dimension(0:zlevels) :: z

    real(8) :: rho

    rho = interp (density (z(k-1)), density (z(k)))
    buoyancy = (ref_density - rho) / ref_density
  end function buoyancy
  
  subroutine apply_initial_conditions
    implicit none
    integer :: d, k, l

    do l = level_end, level_start, -1
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_end, level_start, -1
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_sol,  l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply_initial_conditions

  subroutine update
    ! Update means, bathymetry and penalization mask
    implicit none
    integer :: d, k, l, p

    do d = 1, size(grid)
       do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches
          call apply_onescale_to_patch (set_bathymetry, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          do k = 1, zmax
             call apply_onescale_to_patch (set_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do

    do d = 1, size(grid)
       do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches
          call apply_onescale_to_patch (init_mean, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine update

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initial perturbation to mean for an entire vertical column 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    eta = init_free_surface (dom%node%elts(id_i))
    z_s = dom%topo%elts(id_i)
    
    if (sigma_z) then
       z = z_coords (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (dom, i, j, k, offs, dims)

       if (k == zlevels) then
          sol(S_MASS,zlevels)%data(d)%elts(id_i) = rho * eta
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
       end if
       sol(S_TEMP,k)%data(d)%elts(id_i)                      = rho * dz(k) * buoyancy (z, k)
       sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end do

    if (mode_split) then
       phi = phi_node (d, id_i, zlevels)
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i) = phi * eta ! free surface perturbation
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0.0_8
    end if
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values for an entire vertical column
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    eta = 0.0_8
    z_s = dom%topo%elts(id_i)
    
    if (sigma_z) then
       z = z_coords (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (dom, i, j, k, offs, dims)

       sol_mean(S_MASS,k)%data(d)%elts(id_i)                      = rho * dz(k)
       sol_mean(S_TEMP,k)%data(d)%elts(id_i)                      = 0.0_8
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end do

    if (mode_split) then
       sol_mean(S_MASS,zlevels+1)%data(d)%elts(id_i) = 0.0_8
       sol_mean(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0.0_8
    end if
  end subroutine init_mean
  
  subroutine initialize_a_b_vert
    ! Initialize hybrid sigma-coordinate vertical grid 
    ! (a_vert, b_vert not used if sigma_z = .true.)
    implicit none
    integer               :: k
    real(8)               :: z
    real(8), dimension(6) :: p
    
    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    b_vert(0) = 1.0_8 ; b_vert(zlevels) = 0.0_8
    if (trim (coords) == "uniform") then
       do k = 2, zlevels-1
          b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
       end do
    elseif (trim(coords) == "croco") then
       p = (/ -5.7831,  18.9754, -24.6521,  16.1698, -5.7092, 0.9972 /)
       do k = 1, zlevels-1
          z = dble(k)/dble(zlevels)
          b_vert(k) = p(1)*z**5 + p(2)*z**4 + p(3)*z**3 + p(4)*z**2 + p(5)*z + p(6)
       end do
    end if
    a_vert = 1.0_8 - b_vert

    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert

  subroutine init_tke (dom, i, j, zlev, offs, dims)
    ! Initialize TKE
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer  :: d, id

    d    = dom%id+1
    id   = idx (i, j, offs, dims) + 1

    tke(zlev)%data(d)%elts(id) = 1d-6
  end subroutine init_tke

  real(8) function surf_geopot (p)
    ! Surface geopotential: postive if greater than mean seafloor                                                                                        
    ! Gives minimum depth of approximately 24 m                                                                                                          
    implicit none
    type(Coord) :: p

    real(8) :: amp, b_max, lat, lon, y

    call cart2sph (p, lon, lat)

    b_max = abs (max_depth - min_depth)
    lat = lat / DEG
    
    if (abs(lat-lat_c) <= lat_width/2) then
        amp = b_max / (1.0_8 - tanh (-slope*width/shift))
        y = (lat - (lat_c - lat_width/2))/180 * MATH_PI*radius ! y = 0 at low latitude boundary of channel                                               
	surf_geopot = amp * (1.0_8 - tanh (slope * (f(y) - width/shift)))
     else
      surf_geopot = b_max
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

  real(8) function density (z)
    implicit none
    real(8) :: z

    real(8)            :: delta
    real(8), parameter :: drho = -2.5d0
    character(255)     :: stratification

    delta = abs(max_depth)/3

    stratification = "temp"

    if (trim(stratification) == "temp") then
       density = eqn_of_state (temp_profile(z))
    elseif (trim(stratification) == "linear") then 
       density = ref_density + drho * (max_depth-z)/max_depth
    elseif (trim(stratification) == "exponential") then
       density = ref_density + drho * exp__flush (z/delta)
    elseif (trim(stratification) == "none") then
       density = ref_density 
    end if
  end function density

  real(8) function temp_profile (z)
    implicit none
    real(8) :: z

    real(8)            :: hz, strat, z0, z1
    real(8), parameter :: h0 = 150_8, hz_0 = 6.5_8, T0 = 14_8, z0_0 = -35_8, z1_0 = -75_8

    strat = abs(max_depth)
    hz = hz_0 * abs(max_depth/h0)
    z0 = z0_0 * abs(max_depth/h0)
    z1 = z1_0 * abs(max_depth/h0)

    temp_profile = T0 + 4*tanh ((z - z0) / hz) + (z - z1) / strat
  end function temp_profile

  real(8) function eqn_of_state (temperature)
    implicit none
    real(8) :: temperature

    real(8), parameter :: beta = 0.28, T0 = 14_8

    eqn_of_state = ref_density - beta * (temperature - T0)
  end function eqn_of_state

  subroutine print_density
    implicit none
    integer                       :: k
    real(8)                       :: bv, c_k, c1, drho, dz_l, eta, rho, rho_above, z_k, z_s, z_above
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    eta = 0.0_8
    z_s = max_depth
    
    if (sigma_z) then
       z = z_coords (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    write (6,'(a)') " Layer      z         dz         rho "
    do k = 1, zlevels
       z_k = interp (z(k-1), z(k))
       write (6, '(2x, i2, 4x, 2(es9.2, 1x), es11.5)') &
            k, z_k, dz(k), ref_density * (1.0_8 - buoyancy (z, k))
    end do
    
    write (6,'(/,a)') " Interface     z"
    do k = 0, zlevels
       write (6, '(3x, i3, 5x, es9.2)') k, z(k)
    end do

    write (6,'(/,a)') " Interface      V        c1      CFL_c1"
    c1 = 0.0_8
    do k = 1, zlevels-1
       z_above = interp (z(k),   z(k+1))
       z_k     = interp (z(k-1), z(k))
       dz_l    = z_above - z_k
       
       rho_above = ref_density * (1.0_8 - buoyancy (z, k+1))
       rho  = ref_density * (1.0_8 - buoyancy (z, k))
       drho = rho_above - rho
       
       bv = sqrt(- grav_accel * drho/dz_l/rho)
       c_k = bv * abs(max_depth) / MATH_PI
       c1 = max (c1, c_k)
       
       write (6, '(3x, i3, 5x,3(es9.2,1x))') k, bv, c_k, c1*dt_init/dx_min
    end do
    write (6,'(/,a,es10.4)') "Maximum internal wave speed = ", c1
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

    if (default_thresholds) then ! initialize once
       threshold = threshold_def
    else
       call cal_lnorm_sol (sol, order)
       threshold_new = tol * lnorm
       
       ! Correct very small values
       do k = 1, zmax
          if (threshold_new(S_MASS,k) < threshold_def(S_MASS,k)/10) threshold_new(S_MASS,k) = threshold_def(S_MASS,k)
          if (threshold_new(S_TEMP,k) < threshold_def(S_TEMP,k)/10) threshold_new(S_TEMP,k) = threshold_def(S_TEMP,k)
          if (threshold_new(S_VELO,k) < threshold_def(S_VELO,k)/10) threshold_new(S_VELO,k) = threshold_def(S_VELO,k)
       end do
       
       if (istep >= 10) then
          threshold = 0.01*threshold_new + 0.99*threshold
       else
          threshold = threshold_new
       end if
    end if
  end subroutine set_thresholds

  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer                       :: k
    real(8)                       :: eta, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord)                   :: x_i

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0.0_8
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0.0_8

    eta = 0.0_8
    z_s = max_depth

    if (sigma_z) then
       z = z_coords (eta, z_s)
    else
       z = eta * a_vert + z_s * b_vert
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)
    
    do k = 1, zlevels
       lnorm(S_MASS,k) = ref_density * dz(k)
       lnorm(S_TEMP,k) = drho * dz(k)
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C, C_b, C_divu, C_mu, C_rotu, C_visc, dlat, tau_b, tau_divu, tau_mu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = 0.891 * sqrt (4/sqrt(3.0_8) * area) ! edge length of average triangle

    area = 4*MATH_PI*radius**2/(20*4**min_level)
    dx_max = sqrt (4/sqrt(3.0_8) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, 1.4*dx_min/u_wbc, 1.2*dx_min/c1)
    dt_init = dt_cfl

    C = 0.0_8 ! <= 1/2 if explicit
    C_rotu = C / 4**Laplace_order_init
    C_divu = C
    C_mu   = C
    C_b    = C
    
    ! Diffusion time scales
    tau_mu   = dt_cfl / C_mu
    tau_b    = dt_cfl / C_b
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr(S_MASS) = dx_min**(2*Laplace_order_init) / tau_mu
       visc_sclr(S_TEMP) = dx_min**(2*Laplace_order_init) / tau_b
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,4(a,es8.2),a,/)') &
            "dx_max  = ", dx_max/KM, " dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_mu = ", tau_mu/HOUR, " [h]"
       write (6,'(4(a,es8.2),/)') "C_mu = ", C_mu,  " C_b = ", C_mu, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
    
    ! Penalization parameterss
    dlat = 0.5*npts_penal * (dx_max/radius) / DEG ! widen channel to account for boundary smoothing

    width_S = lat_c - (lat_width/2 + dlat) + 90_8
    width_N = lat_c - (lat_width/2 + dlat)       

    ! Smoothing exponent for land mass
    n_smth_S = 4*radius * width_S*DEG / (dx_max * npts_penal)
    n_smth_N = 4*radius * width_N*DEG / (dx_max * npts_penal)
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

    integer                        :: d, e, id, id_e, id_i, l, nsmth
    real(8)                        :: dx
    type(Coord)                    :: p
    type(Coord), dimension(1:EDGE) :: q

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1
    
    p = dom%node%elts(id_i)
    l = dom%level%elts(id_i)

!!$    dx = max (dx_min, maxval (dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1))) ! local grid size
    dx = dx_max
!!$    dx = dx_min
    
    select case (itype)
    case ("bathymetry")
!!$       nsmth = 2 * (l - min_level)
       nsmth = 1
       dom%topo%elts(id_i) = max_depth + smooth (surf_geopot, p, dx, nsmth) / grav_accel
    case ("penalize") ! analytic land mass with smoothing
       nsmth = 0
       
       penal_node(zlev)%data(d)%elts(id_i) = smooth (mask, p, dx, nsmth)
      
       q(RT+1) = dom%node%elts(idx(i+1, j,   offs, dims)+1)
       q(DG+1) = dom%node%elts(idx(i+1, j+1, offs, dims)+1)
       q(UP+1) = dom%node%elts(idx(i,   j+1, offs, dims)+1)
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = interp (smooth (mask, p, dx, nsmth), smooth (mask, q(e), dx, nsmth))
       end do
    end select
  end subroutine topography

  real(8) function smooth (fun, p, dx, npts)
    ! Smooth a function using radial basis functions
    implicit none
    integer     :: npts
    real(8)     :: dx
    type(Coord) :: p

    integer     :: ii, jj
    real(8)     :: dtheta, lat, lat0, lon, lon0, nrm, r, rbf, wgt
    type(Coord) :: q

    interface
       real(8) function fun (q)
         use geom_mod
         type(Coord) :: q
       end function fun
    end interface

    if (npts == 0) then
       smooth = fun (p)
    else
       dtheta = dx/radius
       call cart2sph (p, lon0, lat0)
       nrm = 0.0_8
       rbf = 0.0_8
       do ii = -npts, npts
          lat = lat0 + dtheta * ii
          do jj = -npts, npts
             lon = lon0 + dtheta * jj

             q = project_on_sphere (sph2cart(lon, lat))
             r = norm (vector(p, q))
             wgt = radial_basis_fun ()

             nrm = nrm + wgt
             rbf = rbf + wgt * fun (q)
          end do
       end do
       smooth = rbf /nrm
    end if
  contains
    real(8) function radial_basis_fun ()
      ! Radial basis function for smoothing topography
      implicit none

      radial_basis_fun = exp (-(r/(npts*dx/2))**2)
    end function radial_basis_fun
  end function smooth

  real(8) function mask (p)
    implicit none
    type(Coord) :: p

    real(8) :: lat, lon

    call cart2sph (p, lon, lat)

    mask = exp__flush (- abs((lat/DEG+90_8)/width_S)**n_smth_S) + exp__flush (- abs((lat/DEG-90_8)/width_N)**n_smth_N)
  end function mask

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid

    if (time/DAY <= 2.0_8) then
       tau_zonal = tau_0 * sin (MATH_PI/4 * time/DAY)
    else
       tau_zonal = tau_0
    end if
    tau_merid = 0.0_8
  end subroutine wind_stress

  real(8) function tau (p)
    ! Magnitude of wind stress at node p
    implicit none
    type(Coord) :: p

    real(8) :: lat, lon, tau_zonal, tau_merid

    call cart2sph (p, lon, lat)

    call wind_stress (lon, lat, tau_zonal, tau_merid)

    tau = sqrt (tau_zonal**2 + tau_merid**2)
  end function tau

  subroutine set_save_level
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level

  function z_coords (eta_surf, z_s)
    ! Hybrid sigma-z vertical coordinates to minimize inclination of layers to geopotential
    ! near the free surface over strong bathymetry gradients.
    ! Reference: similar to Shchepetkin and McWilliams (JCP vol 228, 8985-9000, 2009)
    !
    ! Sets the a_vert parameter that depends on eta_surf (but not b_vert).
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords

    integer                       :: k
    real(8)                       :: cff, cff1, cff2, hc, z_0
    real(8), parameter            :: theta_b = 0d0, theta_s = 7d0
    real(8), dimension(0:zlevels) :: Cs, sc
    
    hc = min (abs(min_depth), abs(Tcline))
    
    cff1 = 1.0_8 / sinh (theta_s)
    cff2 = 0.5d0 / tanh (0.50 * theta_s)
    
    sc(0) = -1.0_8
    Cs(0) = -1.0_8
    cff = 1d0 / dble(zlevels)
    do k = 1, zlevels
       sc(k) = cff * dble (k - zlevels)
       Cs(k) = (1.0_8 - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5d0)) - 0.5d0)
    end do

    z_coords(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1.0_8 - z_0 / z_s
       z_coords(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords

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
       dz0  = (sol(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(id+1)) &
            / porous_density (dom, i, j, zlev, offs, dims)
       
       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)) &
            / porous_density (dom, i+1, j, zlev, offs, dims)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idNE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1)) &
            / porous_density (dom, i+1, j+1, zlev, offs, dims)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idN+1)  + sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)) &
            / porous_density(dom, i, j+1, zlev, offs, dims)
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc
  
  real(8) function upwelling_bottom (dom, i, j, z_null, offs, dims)
    ! Bottom boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    upwelling_bottom = 0.0_8
  end function upwelling_bottom

  real(8) function upwelling_top (dom, i, j, z_null, offs, dims)
    ! Top boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    upwelling_top = 0.0_8
  end function upwelling_top

  function upwelling_drag (dom, i, j, zlev, offs, dims)
    ! Wind stress velocity source term evaluated at edges (top boundary condition for vertical diffusion of velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: upwelling_drag

    integer                         :: d, id, idE, idN, idNE               
    real(8), dimension(1:EDGE)      :: mass_e, tau_wind
    real(8), dimension(0:NORTHEAST) :: full_mass


    id = idx (i, j, offs, dims)

    if (maxval (dom%mask_e%elts(EDGE*id+RT+1:EDGE*id+UP+1)) >= ADJZONE) then
       d = dom%id + 1
       idE  = idx (i+1, j,   offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idNE = idx (i+1, j+1, offs, dims)

       full_mass(0:NORTHEAST) = sol_mean(S_MASS,zlevels)%data(d)%elts((/id,idN,idE,id,id,idNE/)+1) &
            + sol(S_MASS,zlevels)%data(d)%elts((/id,idN,idE,id,id,idNE/)+1)

       mass_e(RT+1) = 0.5 * (full_mass(0) + full_mass(EAST))
       mass_e(DG+1) = 0.5 * (full_mass(0) + full_mass(NORTHEAST))
       mass_e(UP+1) = 0.5 * (full_mass(0) + full_mass(NORTH))

       tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idE+1))
       tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE+1), dom%node%elts(id+1))
       tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idN+1))

       upwelling_drag = tau_wind / mass_e * (1 - penal_edge(zlevels)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
    else
       upwelling_drag = 0.0_8
    end if
  end function upwelling_drag
end module test_case_mod

!  LocalWords:  zlevels
