Module test_case_mod
  ! Module file for Drake passage test case
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, resume_init, save_zlev
  real(8)                              :: dt_cfl, initotalmass, k_T, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Local variables
  real(8)                              :: beta, bv, drho, drho_dz, f0, Rb, Rd, Rey, Ro
  real(8)                              :: bottom_friction, delta, h0, lat_c, lon_c, max_depth, min_depth, mixed_layer, width
  real(8)                              :: radius_earth, omega_earth, scale, halocline, visc
  real(8)                              :: tau_0, wave_friction
  real(4), allocatable, dimension(:,:) :: topo_data
  logical                              :: drag, mean_split
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
    read (fid,*) varname, tol_elliptic
    read (fid,*) varname, coarse_iter
    read (fid,*) varname, fine_iter
    read (fid,*) varname, log_iter
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, tol
    read (fid,*) varname, cfl_num
    read (fid,*) varname, adapt_dt
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

    Rey  = Udim * delta / visc_rotu ! Reynolds number 

    call set_save_level

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case                      = ", trim (test_case)
       write (6,'(A,A)')      "run_id                         = ", trim (run_id)
       write (6,'(A,es10.4)') "scale                          = ", scale
       write (6,'(A,L1)')     "compressible                   = ", compressible
       write (6,'(A,L1)')     "mean_split                     = ", mean_split
       write (6,'(A,L1)')     "mode_split                     = ", mode_split
       write (6,'(A,L1)')     "penalize                       = ", penalize
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
       write (6,'(A,L1)')     "bottom drag                    = ", drag

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius                   [km]  = ", radius / KM
       write (6,'(A,es10.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", abs (min_depth)
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", abs (max_depth)
       write (6,'(A,es11.4)') "halocline                 [m]  = ", abs (halocline)
       write (6,'(A,es11.4)') "mixed layer               [m]  = ", abs (mixed_layer)
       write (6,'(A,es11.4)') "density difference   [kg/m^3]  = ", drho
       write (6,'(A,es11.4)') "Brunt-Vaisala freq      [1/s]  = ", bv
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed           [m/s]  = ", c1
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/bottom_friction / DAY
       write (6,'(A,es11.4)') "wave drag decay           [h]  = ", 1/wave_friction / HOUR
       write (6,'(A,es11.4)') "buoyancy relaxation       [d]  = ", 1/k_T / DAY
       write (6,'(A,es11.4)') "f0 at 45 deg          [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 45 deg       [rad/ms]  = ", beta
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(A,es11.4)') "barotropic Rossby radius [km]  = ", Rd / KM
       write (6,'(A,es11.4,/)') "baroclinic Rossby radius [km]  = ", Rb / KM
       write (6,'(A,es11.4)') "Rossby number                  = ", Ro
       write (6,'(A,es11.4)') "Re (delta_I u_wbc / nu)        = ", Rey 
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
    eta_surf = init_free_surface (x_i)

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
    implicit none
    type(Coord) :: x_i

    real(8) :: lon, lat, rgrc

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))

    surf_geopot = grav_accel*h0 * exp__flush (-(rgrc/width)**2)
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

    real(8) :: density
    
    density = ref_density + drho * exp__flush (z/delta)
    
    buoyancy = (density - ref_density) / ref_density 
  end function buoyancy

  subroutine print_density_pert
    implicit none
    integer     :: k
    real(8)     :: depth, eta_surf, lat, lon, z
    type(Coord) :: p 

    eta_surf = 0.0_8

    p = sph2cart (lon_c, lat_c) ! centre of seamount
!!$    p = sph2cart (lon_c + width/radius / DEG, lat_c) ! edge of seamount
    
    p%x = radius * p%x ; p%y = radius * p%y ; p%z = radius * p%z  

    depth = max_depth + surf_geopot (p) / grav_accel
    write (6,'(a)') " Layer    z       drho"      
    do k = 1, zlevels
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta_surf + (b_vert(k)+b_vert(k-1)) * depth)
       write (6, '(2x,i2, 1x, 2(es9.2,1x))') k, z, buoyancy (p, z)*ref_density
    end do
    write (6,'(/)')
    write (6,'(a)') " Layer    z"    
    do k = 0, zlevels
       z = a_vert(k) * eta_surf + b_vert(k) * depth
       write (6, '(2x,i2, 1x, 2(es9.2,1x))') k, z
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
    integer     :: k
    real(8)     :: dz, eta_surf, z
    type(Coord) :: x_i

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0.0_8
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0.0_8

    x_i = Coord (radius, 0.0_8, 0.0_8)
    eta_surf = 0.0_8
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta_surf + b_vert_mass(k) * max_depth
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta_surf + (b_vert(k)+b_vert(k-1)) * max_depth)

       lnorm(S_MASS,k) = ref_density*dz

       lnorm(S_TEMP,k) = ref_density*dz * abs(buoyancy (x_i, z))
       if (lnorm(S_TEMP,k) == 0.0_8) lnorm(S_TEMP,k) = 1d16

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

    C_sclr = 1d-3
    C_divu = 1d-3
    C_rotu = 1d-3
    
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
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
!!$       visc_rotu = visc
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

    penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
    penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
  end subroutine set_penal

  subroutine topography (dom, i, j, zlev, offs, dims, itype)
    ! Returns bathymetrys
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer            :: id, id_i
    type(Coord)        :: p

    id_i = idx (i, j, offs, dims)

    dom%topo%elts(id_i) = max_depth + surf_geopot (p) / grav_accel
  end subroutine topography

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    ! Idealized zonally and temporally averaged zonal and meridional wind stresses
    ! (based on Figure 4 from Ferreira et al J Climate 24, 992-1012 (2011) and Figure 4 from Gille J Atmos Ocean Tech 22, 1353-1372 respectively)

    implicit none
    real(8) :: lat, lon, peak, tau_zonal, tau_merid
    
    logical, parameter :: merid_stress = .false.

    peak = (abs(lat)*180/MATH_PI - 35.0_8) / 20
    tau_zonal = -tau_0 * 1.2 * exp (-peak**2) * sin (abs(lat)*6) - 5d-3*exp(-(lat*180/MATH_PI/10)**2)

    if (merid_stress) then
       peak = lat*180/MATH_PI / 15
       tau_merid = -tau_0 * 2.5 * exp (-peak**2) * sin (2*lat) * peak**2
    else
       tau_merid = 0.0_8
    end if
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

  subroutine trend_relax (q, dq)
    ! Trend relaxation to mean buoyancy
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (q, NONE, 27)
    
    do k = 1, zlevels
       ! Scalars
       do d = 1, size(grid)
          temp  =>  q(S_TEMP,k)%data(d)%elts
          dmass => dq(S_MASS,k)%data(d)%elts
          dtemp => dq(S_TEMP,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
          end do
          nullify (dmass, dscalar, temp)
       end do

       ! Velocity and mass 
       do d = 1, size(grid)
          dvelo => dq(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (trend_velo, grid(d), p-1, k, 0, 0)
          end do
          nullify (dvelo)
       end do
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_relax

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Relax buoyancy to mean
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    
    dmass(id_i) = 0.0_8
    dtemp(id_i) = - temp(id_i) * k_T
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims)

    dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
  end subroutine trend_velo
end module test_case_mod
