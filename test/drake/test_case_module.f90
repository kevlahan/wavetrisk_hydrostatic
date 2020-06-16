Module test_case_mod
  ! Module file for Drake passage test case
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, etopo_res, npts_penal, resume_init, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Local variables
  real(8)                              :: beta, bv, delta_I, delta_M, delta_S, delta_sm, drho, drho_dz, f0, L_R, Rey, Ro
  real(8)                              :: bottom_friction, max_depth, min_depth, mixed_layer, scale, halocline, u_wbc
  real(8)                              :: resolution, tau_0, wave_friction
  real(4), allocatable, dimension(:,:) :: topo_data
  logical                              :: drag, etopo_coast, mean_split
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
    read (fid,*) varname, level_fill
    read (fid,*) varname, zlevels
    read (fid,*) varname, remap
    read (fid,*) varname, iremap
    read (fid,*) varname, coarse_iter
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, tol
    read (fid,*) varname, cfl_num
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    read (fid,*) varname, alpha
    read (fid,*) varname, drag
    read (fid,*) varname, npts_penal
    read (fid,*) varname, resolution

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

    delta_M = (visc_rotu/beta)**(1.0_8/(2*Laplace_order_init+1)) * METRE ! Munk layer scale

    ! Bottom drag
    if (drag) then
       bottom_friction = beta * delta_M/4
!!$      bottom_friction = 1.0_8/(110*DAY) * METRE/SECOND ! linear bottom friction coefficient (nemo is 4e-4/H for H = 4000m)
    else
       bottom_friction = 0.0_8
    end if

    ! Internal wave drag to reduce oscillation of internal wave (ensure stability)
    if (drho == 0.0_8) then
       wave_friction = 0.0_8
    else
       wave_friction = min (1/(1000/bv), 1/dt_init) 
    end if

    delta_S = bottom_friction / beta      ! Stommel layer (want delta_S = delta_M/4)
    Rey     = u_wbc * delta_I / visc_rotu ! Reynolds number of western boundary current
    Ro      = u_wbc / (delta_M*f0)        ! Rossby number (based on boundary current)

    call set_save_level

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case            = ", trim (test_case)
       write (6,'(A,A)')      "run_id               = ", trim (run_id)
       write (6,'(A,L1)')     "compressible         = ", compressible
       write (6,'(A,L1)')     "mean_split           = ", mean_split
       write (6,'(A,L1)')     "mode_split           = ", mode_split
       write (6,'(A,L1)')     "penalize             = ", penalize
       write (6,'(A,i3)')     "npts_penal           = ", npts_penal
       write (6,'(A,i3)')     "min_level            = ", min_level
       write (6,'(A,i3)')     "max_level            = ", max_level
       write (6,'(A,i3)')     "level_fill           = ", level_fill
       write (6,'(A,i5)')     "number of domains    = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL         = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL          = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels              = ", zlevels
       write (6,'(A,L1)')     "remap                = ", remap
       write (6,'(a,a)')      "remapscalar_type     = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type       = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap               = ", iremap
       write (6,'(a,i3)')     "coarse_iter          = ", coarse_iter
       write (6,'(A,L1)')     "adapt_trend          = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds   = ", default_thresholds
       write (6,'(A,L1)')     "perfect              = ", perfect
       write (6,'(A,es10.4)') "tolerance            = ", tol
       write (6,'(A,i1)')     "optimize_grid        = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt             = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num              = ", cfl_num
       write (6,'(a,a)')      "timeint_type         = ", trim (timeint_type)
       write (6,'(A,i1)')     "Laplace_order        = ", Laplace_order_init
       write (6,'(A,i1)')     "n_diffuse            = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write [d]         = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY             = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance            = ", rebalance
       write (6,'(A,es10.4)') "time_end [d]         = ", time_end/DAY
       write (6,'(A,i6)')     "resume               = ", resume_init
       write (6,'(A,L1)')     "bottom drag          = ", drag
       write (6,'(A,i3)')     "npts_penal           = ", npts_penal
       write (6,'(A,i3)')     "etopo_res            = ", etopo_res
       write (6,'(A,es10.4)') "resolution           = ", resolution

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius             [km]     = ", radius / KM
       write (6,'(A,es10.4)') "omega              [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density        [kg/m^3] = ", ref_density
       write (6,'(A,es10.4)') "grav accel         [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "min_depth          [m]      = ", abs (min_depth)
       write (6,'(A,es11.4)') "max_depth          [m]      = ", abs (max_depth)
       write (6,'(A,es11.4)') "halocline          [m]      = ", abs (halocline)
       write (6,'(A,es11.4)') "mixed layer        [m]      = ", abs (mixed_layer)
       write (6,'(A,es11.4)') "density difference [kg/m^3] = ", drho
       write (6,'(A,es11.4)') "Brunt-Vaisala freq [1/s]    = ", bv
       write (6,'(A,es11.4)') "c0 wave speed      [m/s]    = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed      [m/s]    = ", c1
       write (6,'(A,es11.4)') "max wind stress    [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "eta (permeability) [s]      = ", eta
       write (6,'(A,es11.4)') "alpha (porosity)            = ", alpha
       write (6,'(A,es11.4)') "bottom friction    [m/s]    = ", bottom_friction
       write (6,'(A,es11.4)') "bottom drag decay  [d]      = ", 1/bottom_friction / DAY
       write (6,'(A,es11.4)') "wave drag decay    [h]      = ", 1/wave_friction / HOUR
       write (6,'(A,es11.4)') "f0 at 30 deg       [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 30 deg     [rad/ms] = ", beta
       write (6,'(A,es11.4)') "dx_max             [km]     = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min             [km]     = ", dx_min   / KM
       write (6,'(A,es11.4)') "L_R at 30 deg      [km]     = ", L_R      / KM
       write (6,'(A,es11.4)') "Inertial layer     [km]     = ", delta_I  / KM
       write (6,'(A,es11.4)') "Munk layer         [km]     = ", delta_M  / KM
       write (6,'(A,es11.4)') "Stommel layer      [km]     = ", delta_S  / KM
       write (6,'(A,es11.4,/)') "submesoscale       [km]     = ", delta_sm / KM
       write (6,'(A,es11.4)') "Rossby number               = ", Ro
       write (6,'(A,es11.4)') "Resolution of Munk layer    = ", resolution
       write (6,'(A,es11.4)') "Re (delta_I u_wbc / nu)     = ", Rey 
       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'

       call print_density_pert
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
    use wavelet_mod
    implicit none
    integer :: d, k, l

    call topography_data

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
    real (8)    :: dz, eta_surf, phi, porous_density, z
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)
    eta_surf = init_free_surface (x_i)

    if (zlev == zlevels+1) then ! 2D barotropic mode
       phi = 1.0_8 + (alpha - 1.0_8) * penal_node(zlevels)%data(d)%elts(id_i)

       sol(S_MASS,zlev)%data(d)%elts(id_i) = phi * eta_surf ! free surface perturbation
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
    else ! 3D layers
       dz = a_vert_mass(zlev) * eta_surf + b_vert_mass(zlev) * dom%topo%elts(id_i)
       z = 0.5 * ((a_vert(zlev)+a_vert(zlev-1)) * eta_surf + (b_vert(zlev)+b_vert(zlev-1)) * dom%topo%elts(id_i))

       porous_density = ref_density * (1.0_8 + (alpha - 1.0_8) * penal_node(zlev)%data(d)%elts(id_i))

       if (mean_split) then
          if (zlev == zlevels) then
             sol(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * eta_surf
          else
             sol(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
          end if
          sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
       else
          sol(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * dz
          sol(S_TEMP,zlev)%data(d)%elts(id_i) = sol(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy (x_i, z)
       end if
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
    real (8)    :: dz, eta_surf, porous_density, z
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)
    eta_surf  = init_free_surface (x_i)

    if (zlev == zlevels+1) then
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
    else
       dz = a_vert_mass(zlev) * eta_surf + b_vert_mass(zlev) * dom%topo%elts(id_i)
       z = 0.5 * ((a_vert(zlev)+a_vert(zlev-1)) * eta_surf + (b_vert(zlev)+b_vert(zlev-1)) * dom%topo%elts(id_i))

       porous_density = ref_density * (1.0_8 + (alpha - 1.0_8) * penal_node(zlev)%data(d)%elts(id_i))

       if (mean_split) then
          sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * dz
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy (x_i, z)
       else
          sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
       end if
    end if
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
  end subroutine init_mean

  subroutine update
    ! Update means, bathymetry and penalization mask
    implicit none
    integer :: d, k, p

    if (resume /= NONE) call topography_data

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

  real(8) function surf_geopot (x_i)
    ! Surface geopotential: postive if greater than mean seafloor
    ! MUST BE SET EQUAL TO ZERO FOR THIS test case
    implicit none
    type(Coord) :: x_i

    surf_geopot = grav_accel * 0.0_8
  end function surf_geopot

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0.0_8
  end function init_free_surface

  real(8) function buoyancy (x_i, z)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    if (zlevels /= 1 .and. z >= halocline) then
       buoyancy = - (1.0_8 - z/halocline) * drho/ref_density
    else
       buoyancy = 0.0_8
    end if
  end function buoyancy

  subroutine print_density_pert
    implicit none
    integer :: k
    real(8) :: eta_surf, z
    type(Coord), parameter :: x_i = Coord (0.0_8, 0.0_8, 0.0_8)

    eta_surf = 0.0_8

    write (6,'(a)') " Layer    z       drho"      
    do k = 1, zlevels
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta_surf + (b_vert(k)+b_vert(k-1)) * max_depth)
       write (6, '(2x,i2, 1x, 2(es9.2,1x))') k, z, - buoyancy (x_i, z)*ref_density
    end do
    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density_pert

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
    integer :: k
    real(8) :: dz, eta_surf, z

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0.0_8
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0.0_8

    eta_surf = 0.0_8
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta_surf + b_vert_mass(k) * max_depth
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta_surf + (b_vert(k)+b_vert(k-1)) * max_depth)
       lnorm(S_MASS,k) = ref_density*dz
       if (drho /= 0.0_8) then
          lnorm(S_TEMP,k) = abs(drho)*dz
       else
          lnorm(S_TEMP,k) = 1d16
       end if
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, C_visc, tau_divu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

    area = 4*MATH_PI*radius**2/(20*4**min_level)
    dx_max = sqrt (4/sqrt(3.0_8) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, dx_min/c1)
    dt_init = dt_cfl

    ! Permeability penalization parameter
    eta = dt_cfl

    ! Diffusion constants
    C_visc = min (1/31d0, max (dt_cfl * beta * dx_min * resolution**3, 1d-4))  ! ensure stability and that Munk layer is resolved with resolution grid points
    resolution = (C_visc/(dt_cfl * beta * dx_min))**(1/3d0)
    C_rotu = C_visc
    C_divu = C_visc * 4
    C_sclr = C_visc * 4

    ! Diffusion time scales
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
       write (6,'(/,4(a,es8.2),a,/)') &
            "dx_max  = ", dx_max/KM, " dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
       write (6,'(a,es10.4,a)') "eta = ", eta," [s]"
    end if
  end subroutine initialize_dt_viscosity

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topography (dom, i, j, zlev, offs, dims, 'bathymetry', npts_penal)
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
       call topography (dom, i, j, zlev, offs, dims, "penalize", npts_penal)
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
    end if
  end subroutine set_penal

  subroutine topography (dom, i, j, zlev, offs, dims, itype, npts)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: npts
    character(*)                   :: itype

    integer     :: d, e, id, id_e, id_i, ii, is0, it0, jj, s, t
    real(8)     :: dx, lat, lon, mask, M_topo, r, s0, t0, sw_topo, topo_sum, wgt
    type(Coord) :: p, q

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    call cart2sph (dom%node%elts(id_i), lon, lat)
    s0 = lon/DEG * BATHY_PER_DEG
    t0 = lat/DEG * BATHY_PER_DEG
    is0 = nint (s0); it0 = nint (t0)
    p = proj_lon_lat (s0, t0)

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth + surf_geopot (p) / grav_accel
    case ("penalize")
       if (npts == 0) then ! no smoothing
          mask = topo_data (is0, it0)
       else ! smoothing
          dx = max (dx_min, maxval (dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1)))
          sw_topo  = 0.0_8
          topo_sum = 0.0_8
          do ii = -npts, npts
             s = is0+ii
             do jj = -npts, npts
                t = it0+jj
                call wrap_lonlat (s, t)
                q = proj_lon_lat (dble(s), dble(t))
                r = norm (vector(p, q))
                wgt = radial_basis_fun (r, npts, dx)
                M_topo = topo_data (s, t)
                topo_sum = topo_sum + wgt * M_topo
                sw_topo  = sw_topo  + wgt
             end do
          end do
          mask = topo_sum / sw_topo
       end if
       penal_node(zlev)%data(d)%elts(id_i) = mask
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = max (penal_edge(zlev)%data(d)%elts(id_e), mask)
       end do
    end select
  end subroutine topography

  real(8) function radial_basis_fun (r, npts, dx)
    ! Radial basis function for smoothing topography
    implicit none
    integer :: npts
    real(8) :: r, dx

    radial_basis_fun = exp (-(r/(npts*dx/2))**2)
  end function radial_basis_fun

  subroutine wrap_lonlat (s, t)
    ! Longitude: wraparound allows for values outside [-180,180]
    ! Latitude: works only if there is no coast at the pole
    implicit none
    integer :: s, t

    if (t < lbound (topo_data,2)) t = lbound (topo_data,2) ! pole
    if (t > ubound (topo_data,2)) t = ubound (topo_data,2) ! pole
    if (s < lbound (topo_data,1)) s = s + 360*BATHY_PER_DEG
    if (s > ubound (topo_data,1)) s = s - 360*BATHY_PER_DEG
  end subroutine wrap_lonlat

  type(Coord) function proj_lon_lat (s, t)
    implicit none
    real(8) :: s, t
    real(8) :: lon, lat
    
    lon = s * DEG / BATHY_PER_DEG
    lat = t * DEG / BATHY_PER_DEG
    proj_lon_lat = project_on_sphere (sph2cart(lon, lat))
  end function proj_lon_lat

  subroutine topography_data
    ! Defines analytic latitude-longitude topography data for Drake passage case
    implicit none
    integer                              :: ii, ilat, ilon, lat_avg, lon_avg, jj, kk
    integer, parameter                   :: npts = 2, n_smooth = 20
    real(8)                              :: avg, lat, lon
    real(8), parameter                   :: width = 30, lat_max = 70, lat_min = -35
    real(8), dimension(2)                :: sz
    real(8), dimension(:,:), allocatable :: temp_data

    if (etopo_coast) then
       if (rank == 0) write(6,'(a)') 'Reading bathymetry data'
       bathy_per_deg = 60/etopo_res
       
       allocate (topo_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))
       open (unit=1086,file='bathymetry') ! "bathymetry" is symbolic link to appropriate etopo bathymetry data
       do kk = ubound (topo_data,2), lbound (topo_data,2), -1 ! north to south (as read from file)
          read (1086,*) topo_data(:,kk)
       end do
       close (1086)
    else
       ! Resolution
       bathy_per_deg = max (1, nint (2*MATH_PI*radius/dx_min/360))
       etopo_res = nint (60/(2*MATH_PI*radius/dx_min/360)) ! effective resolution of topography in arcminutes

       if (.not. allocated (topo_data)) &
            allocate (topo_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))

       if (.not. allocated (temp_data)) &
            allocate (temp_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))

       sz(1) = ubound (topo_data,1) - lbound (topo_data,1)
       sz(2) = ubound (topo_data,2) - lbound (topo_data,2)

       ! Heaviside mask
       temp_data = 0.0_8
       do ilon = lbound (topo_data,1), ubound (topo_data,1)
          lon = ilon / dble(bathy_per_deg)
          do ilat = lbound (topo_data,2), ubound (topo_data,2)
             lat = ilat / dble(bathy_per_deg)
             if (abs(lon) < width/2/cos(lat*DEG) &
                  .and. lat < lat_max .and. lat > lat_min &
                  .and. lat < acos(cos(lat_max*DEG)/cos(lon*DEG)) * 180/MATH_PI &
                  .and. lat > (acos(cos(lat_min*DEG)/cos(lon*DEG))-MATH_PI/2) * 180/MATH_PI) &
                  temp_data(ilon,ilat) = land
          end do
       end do
       topo_data = temp_data

       ! Smoothing
       do kk = 1, n_smooth
          do ilon = lbound (topo_data,1), ubound (topo_data,1)
             do ilat = lbound (topo_data,2), ubound (topo_data,2)
                avg = 0.0_8
                do ii = -npts, npts
                   lon_avg = ilon + ii
                   if (lon_avg < lbound (topo_data,1)) lon_avg = lon_avg + sz(1)
                   if (lon_avg > ubound (topo_data,1)) lon_avg = lon_avg - sz(1)
                   do jj = -npts, npts
                      lat_avg = ilat + jj
                      if (lat_avg < lbound (topo_data,2)) lat_avg = lat_avg + sz(2)
                      if (lat_avg > ubound (topo_data,2)) lat_avg = lat_avg - sz(2)
                      avg = avg + temp_data(lon_avg,lat_avg)
                   end do
                end do
                topo_data(ilon,ilat) = avg/(2*npts+1)**2
             end do
          end do
          temp_data = topo_data
       end do
    end if
  end subroutine topography_data

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    ! Idealized zonally and temporally averaged zonal and meridional wind stresses
    ! (based on Figure 4 from Ferreira et al J Climate 24, 992-1012 (2011) and Figure 4 from Gille J Atmos Ocean Tech 22, 1353-1372 respectively)

    implicit none
    real(8) :: lat, lon, peak, tau_zonal, tau_merid

    peak = (abs(lat)*180/MATH_PI - 35.0_8) / 20
    tau_zonal = -tau_0 * 1.2 * exp (-peak**2) * sin (abs(lat)*6) - 5d-3*exp(-(lat*180/MATH_PI/10)**2)

    !peak = lat*180/MATH_PI / 15
    !tau_merid = -tau_0 * 2.5 * exp (-peak**2) * sin (2*lat) * peak**2
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

    if (zlevels == 2) then 
       a_vert(0) = 0.0_8; a_vert(1) = 0.0_8;               a_vert(2) = 1.0_8
       b_vert(0) = 1.0_8; b_vert(1) = halocline/max_depth; b_vert(2) = 0.0_8
    elseif (zlevels == 3) then
       a_vert(0) = 0.0_8; a_vert(1) = 0.0_8;               a_vert(2) = 0.0_8;                 a_vert(3) = 1.0_8 
       b_vert(0) = 1.0_8; b_vert(1) = halocline/max_depth; b_vert(2) = mixed_layer/max_depth; b_vert(3) = 0.0_8
    else ! uniform sigma grid: z = a_vert*eta + b_vert*z_s
       do k = 0, zlevels
          a_vert(k) = dble(k)/dble(zlevels)
          b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
       end do
    end if
    
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
end module test_case_mod
