Module test_case_mod
  ! Module file for tsunami test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use vert_diffusion_mod
  implicit none

  ! Standard variables
  integer :: CP_EVERY, resume_init
  real(8) :: dt_cfl, total_cpu_time, Hdim, Ldim, Tdim, Udim
  real(8) :: g_earth, H_earth, H_norm, L_norm, U_norm, T_norm

  ! Local variables
  integer                              :: bathy_per_deg, etopo_res
  real(4), allocatable, dimension(:,:) :: etopo_data
  real(8)                              :: beta, bv, delta_I, delta_M, delta_S, delta_sm
  real(8)                              :: drho, drho_dz, f0, Fr, Ku, k_T, lambda0, lambda1, Rb, Rd, Rey, Ro, radius_earth
  real(8)                              :: omega_earth, scale, scale_omega, mixed_layer, tau_0, thermocline, u_wbc 
  real(8),                      target :: bottom_friction_case
  logical                              :: etopo_bathy, etopo_coast
  logical                              :: aligned
  logical                              :: normalized
  character(255)                       :: coords

  real(8)                              :: dH, lon_c, lat_c, pert_radius
  real(8), allocatable, dimension(:,:) :: analytic_data

  ! Drake land boundaries
  real(8),                   parameter :: lat_max = 60d0*DEG, lat_min = -35d0*DEG, lon_min = -15*DEG, lon_max = 15d0*DEG
contains
  subroutine assign_functions
    ! Assigns generic pointer functions to functions defined in test cases
    implicit none

    ! Standard functions
    apply_initial_conditions => apply_initial_conditions_case
    dump                     => dump_case
    load                     => load_case
    initialize_a_b_vert      => initialize_a_b_vert_case
    initialize_dt_viscosity  => initialize_dt_viscosity_case
    initialize_thresholds    => initialize_thresholds_case
    physics_scalar_flux      => physics_scalar_flux_case
    physics_velo_source      => physics_velo_source_case
    set_save_level           => set_save_level_case
    set_thresholds           => set_thresholds_case
    surf_geopot              => surf_geopot_case
    update                   => update_case
    z_coords                 => z_coords_case

    ! Needed for vertical diffusion
    bottom_friction  => bottom_friction_case
    bottom_buoy_flux => bottom_buoy_flux_case
    top_buoy_flux    => top_buoy_flux_case
    wind_flux        => wind_flux_case
    tau_mag          => tau_mag_case
  end subroutine assign_functions

  subroutine read_test_case_parameters
#ifdef MPI
    use mpi
#endif
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
    if (rank == 0) then
       write (6,'(A,A)') "Input file = ", trim (filename)
       open(unit=fid, file=filename, action='READ')
       
       read (fid,*) varname, test_case
       read (fid,*) varname, run_id
       read (fid,*) varname, scale
       read (fid,*) varname, scale_omega
       read (fid,*) varname, max_level
       read (fid,*) varname, zlevels
       read (fid,*) varname, tol
       read (fid,*) varname, dt_write
       read (fid,*) varname, CP_EVERY
       read (fid,*) varname, time_end
       read (fid,*) varname, resume_init
       close(fid)
    end if

#ifdef MPI
    call MPI_Bcast (test_case,        255, MPI_BYTE,             0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (run_id,           255, MPI_BYTE,             0, MPI_COMM_WORLD, ierror)    
    call MPI_Bcast (scale,              1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (scale_omega,        1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (max_level,          1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (zlevels,            1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (tol,                1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (dt_write,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (CP_EVERY,           1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (time_end,           1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (resume_init,        1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror) 
#endif
    
    press_save = 0d0
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    call set_save_level_case

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case                      = ", trim (test_case)
       write (6,'(A,A)')      "run_id                         = ", trim (run_id)
       write (6,'(A,es10.4)') "scale                          = ", scale
       write (6,'(A,es10.4)') "scale_omega                    = ", scale_omega
       write (6,'(A,L1)')     "compressible                   = ", compressible
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
       write (6,'(A,i3)')     "save_zlev                      = ", save_zlev
       write (6,'(A,L1)')     "remap                          = ", remap
       write (6,'(a,a)')      "remapscalar_type               = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type                 = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap                         = ", iremap
       write (6,'(A,L1)')     "sigma_z                        = ", sigma_z
       write (6,'(A,es10.4)') "tol_elliptic                   = ", tol_elliptic
       write (6,'(a,i3)')     "coarse_iter                    = ", coarse_iter
       write (6,'(a,i3)')     "fine_iter                      = ", fine_iter
       write (6,'(a,l1)')     "log_iter                       = ", log_iter
       write (6,'(A,L1)')     "default_thresholds             = ", default_thresholds
       write (6,'(A,es10.4)') "tolerance                      = ", tol
       write (6,'(A,i1)')     "optimize_grid                  = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt                       = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num                        = ", cfl_num
       write (6,'(a,a)')      "timeint_type                   = ", trim (timeint_type)
       write (6,'(A,i1)')     "Laplace_order                  = ", Laplace_order_init
       write (6,'(/,a,/)') "Scale-aware horizontal diffusion"
       write (6,'(3(a,es8.2/))') "C_visc(S_MASS) = ", C_visc(S_MASS), "C_visc(S_TEMP) = ", C_visc(S_TEMP), &
            "C_visc(S_VELO) = ", C_visc(S_VELO)
       write (6,'(a,/,a,/,/,a,/,a,/)') "Stability limits:", &
            "[Klemp 2017 Damping Characteristics of Horizontal Laplacian Diffusion Filters Mon Weather Rev 145, 4365-4379.]", &
            "C_visc(S_MASS) and C_visc(S_TEMP) <  (1/6)**Laplace_order", &
            "                   C_visc(S_VELO) < (1/24)**Laplace_order"
       write (6,'(A,i1)')     "n_diffuse                      = ", n_diffuse
       write (6,'(A,L1)')     "vert_diffuse                   = ", vert_diffuse
       write (6,'(A,L1)')     "tke_closure                    = ", tke_closure
       write (6,'(A,es10.4)') "dt_write [d]                   = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY                       = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance                      = ", rebalance
       write (6,'(A,es10.4)') "time_end [d]                   = ", time_end/DAY
       write (6,'(A,i6)')     "resume                         = ", resume_init
       write (6,'(A,i3)')     "etopo_res                      = ", etopo_res

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius                   [km]  = ", radius / KM
       write (6,'(A,es10.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", abs (max_depth)
       write (6,'(A,es11.4)') "mixed_layer               [m]  = ", abs (mixed_layer)
       write (6,'(A,es11.4)') "thermocline               [m]  = ", abs (thermocline)
       write (6,'(a,a)')      "vertical coordinates           = ", trim (coords)
       write (6,'(A,es11.4)') "density perturbation [kg/m^3]  = ", drho
       write (6,'(A,es11.4)') "Brunt-Vaisala freq      [1/s]  = ", bv
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed           [m/s]  = ", c1
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay time    [d]  = ", abs(max_depth)/bottom_friction_case / DAY
       if (vert_diffuse)  then
          if (tke_closure) then
             write (6,'(A,es11.4)') "Kv_0                  [m^2/s]  = ", Kv_0
             write (6,'(A,es11.4)') "Kt_0                  [m^2/s]  = ", Kt_0
             write (6,'(A,es11.4)') "Kt_max                [m^2/s]  = ", Kt_max
             write (6,'(a,l1)')     "tke_closure                    = ", tke_closure
             write (6,'(a,l1)')     "patankar                       = ", patankar
             write (6,'(a,l1)')     "enhance_diff                   = ", enhance_diff
          else
             write (6,'(A,es11.4)') "Kv_bottom             [m^2/s]  = ", Kv_bottom
             write (6,'(A,es11.4)') "Kt_const              [m^2/s]  = ", Kt_const
          end if
       elseif (zlevels == 2) then
          write (6,'(A,es11.4)') "Ku                    [m^2/s]  = ", Ku
       end if
       write (6,'(A,es11.4)') "buoyancy relaxation       [d]  = ", 1d0/k_T / DAY
       write (6,'(A,es11.4)') "f0 at 45 deg          [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 45 deg       [rad/ms]  = ", beta
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(A,es11.4)') "dt_cfl                   [s]   = ", dt_cfl
       write (6,'(A,es11.4)') "External scale l0        [km]  = ", lambda0  / KM
       write (6,'(A,es11.4)') "Internal scale l1        [km]  = ", lambda1  / KM
       write (6,'(A,es11.4)') "Inertial layer           [km]  = ", delta_I  / KM
       write (6,'(A,es11.4)') "Munk layer               [km]  = ", delta_M  / KM
       write (6,'(A,es11.4)') "Stommel layer            [km]  = ", delta_S  / KM
       write (6,'(A,es11.4)') "submesoscale             [km]  = ", delta_sm / KM
       write (6,'(A,es11.4)') "barotropic Rossby radius [km]  = ", Rd / KM
       write (6,'(A,es11.4,/)') "baroclinic Rossby radius [km]  = ", Rb / KM
       write (6,'(A,es11.4)') "Rossby number                  = ", Ro
       write (6,'(A,es11.4)') "Froude number                  = ", Fr
       write (6,'(A,es11.4)') "Re (delta_I u_wbc / nu)        = ", Rey 
       write (6,'(A,es11.4)') "Resolution of Munk layer       = ", delta_M / dx_min
       write (6,'(A,es11.4)') "Resolution of Taylor scale     = ", delta_I / sqrt(Rey) / dx_min
       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'

       call print_density
    end if
  end subroutine print_test_case_parameters

  subroutine print_density
    implicit none
    integer                       :: k
    real(8)                       :: bv, c_k, c1, drho, dz_l, eta, lat, rho, rho_above, z_k, z_s, z_above
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord) :: x_i

    lat = 0d0 ! latitude to evaluate buoyancy
    
    eta = 0d0
    z_s = max_depth
    !x_i = Coord (radius, 0d0, 0d0)
    x_i = sph2cart (0d0, 0d0)
    
    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    write (6,'(a)') " Layer      z         dz         drho "
    do k = 1, zlevels
       z_k = interp (z(k-1), z(k))
       write (6, '(2x, i2, 4x, 2(es9.2, 1x), es12.5)') &
            k, z_k, dz(k), -ref_density * buoyancy_init (x_i, z_k)
    end do
    
    write (6,'(/,a)') " Interface     z"
    do k = 0, zlevels
       write (6, '(3x, i3, 5x, es9.2)') k, z(k)
    end do

    write (6,'(/,a)') " Interface      V        c1      CFL_c1"
    c1 = 0d0
    do k = 1, zlevels-1
       z_above = interp (z(k),   z(k+1))
       z_k     = interp (z(k-1), z(k))
       dz_l    = z_above - z_k
       
       rho_above = ref_density * (1d0 - buoyancy_init (x_i, z_above))
       rho  = ref_density * (1d0 - buoyancy_init (x_i, z_k))
       drho = rho_above - rho
       
       bv = sqrt(grav_accel * abs(drho)/dz_l/rho)
       c_k = bv * abs(max_depth) / MATH_PI
       c1 = max (c1, c_k)
       
       write (6, '(3x, i3, 5x,3(es9.2,1x))') k, bv, c_k, c_k*dt_init/dx_min
    end do
    write (6,'(/,a,es8.2)') "Maximum internal wave speed [m/s] = ", c1
    write (6,'(a,es8.2)')   "Maximum baroclinic CFL number     = ", c1*dt_init/dx_min
    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

    if (rank == 0) then
       open (unit=12, file=trim (run_id)//'_log', action='WRITE', form='FORMATTED', position='APPEND')
       if (log_mass) then
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
       else
          write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i12,2(a,es9.2,1x))') &
               'time [d] = ', time/DAY, &
               ' dt [s] = ', dt, &
               '  mass tol = ', threshold(S_MASS,zlevels), &
               ' temp tol = ', threshold(S_TEMP,zlevels), &
               ' velo tol = ', threshold(S_VELO,zlevels), &
               ' Jmax = ', level_end, &
               ' dof = ', sum (n_active), &
               ' balance = ', rel_imbalance, &
               ' cpu = ', timing

          write (12,'(5(es15.9,1x),i2,1x,i12,1x,2(es15.9,1x))')  time/DAY, dt, &
               threshold(S_MASS,zlevels), threshold(S_TEMP,zlevels), threshold(S_VELO,zlevels), &
               level_end, sum (n_active), rel_imbalance, timing
       end if
       close (12)
    end if
  end subroutine print_log

  subroutine read_etopo_data
    ! Read in etopo bathymetry data
    implicit none
    integer :: k

    bathy_per_deg  = 60/etopo_res                     ! bathymetry data poin
    
    if (rank == 0) write(6,'(a)') 'Reading bathymetry data ...'
    
    allocate (etopo_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))
    
    open (unit=1086,file='bathymetry')
    
    do k = ubound (etopo_data,2), lbound (etopo_data,2), -1 ! north to south (as read from file)
       read (1086,*) etopo_data(:,k)
    end do
    close (1086)
  end subroutine read_etopo_data

  subroutine analytic_topo_data
    implicit none
    integer            :: ii, jj, nlat, nlon
    integer, parameter :: npts = 1, nsmth = 22 ! yes, but internal grid refinement
    !integer, parameter :: npts = 2, nsmth = 12 ! yes, but internal grid refinement 
    !integer, parameter :: npts = 3, nsmth = 8 ! yes
    !integer, parameter :: nsmth = 20 ! shapiro
    real(8)            :: lat, lon

    bathy_per_deg  = 60/etopo_res ! topography data per degree

    ! Set up land data
    nlat = 180 * BATHY_PER_DEG
    nlon = 360 * BATHY_PER_DEG

    allocate (analytic_data(nlon,nlat))

    do jj = 1, nlat
       lat = -90d0 + dble (jj-1) / dble (BATHY_PER_DEG)
       do ii = 1, nlon
          lon = -180d0 + dble (ii-1) / dble (BATHY_PER_DEG)

          analytic_data(ii,jj) = drake_land (lon*DEG, lat*DEG)
       end do
    end do

    call smoothing_rbf (dx_max, npts, nsmth, analytic_data)
  end subroutine analytic_topo_data

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initialize mean values for an entire vertical column
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)

    eta = init_free_surface (x_i)
    z_s = dom%topo%elts(id_i)
    
    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (d, id_i, k)
       z_k = interp (z(k-1), z(k))

       if (k == zlevels) then
          sol(S_MASS,k)%data(d)%elts(id_i) = rho * eta
          !sol(S_MASS,k)%data(d)%elts(id_i) = rho * (eta + dz(k))
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = 0d0
          !sol(S_MASS,k)%data(d)%elts(id_i) = rho * dz(k)
       end if
       sol(S_TEMP,k)%data(d)%elts(id_i)                      = rho * dz(k) * buoyancy_init (x_i, z_k)
       sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do

    if (mode_split) then
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i)                      = eta
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i)                      = 0d0
       sol(S_VELO,zlevels+1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
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
    real(8)                       :: eta, phi, rho, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)

    eta = init_free_surface (x_i)
    z_s = dom%topo%elts(id_i)
    
    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (d, id_i, k)
       z_k = interp (z(k-1), z(k))
       
       sol_mean(S_MASS,k)%data(d)%elts(id_i)                      = rho * dz(k)
       !sol_mean(S_MASS,k)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_TEMP,k)%data(d)%elts(id_i)                      = rho * dz(k) * buoyancy_init (x_i, z_k)
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do

    if (mode_split) then
       sol_mean(S_MASS,zlevels+1)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_TEMP,zlevels+1)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_VELO,zlevels+1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end if
  end subroutine init_mean

  real(8) function surf_geopot_case (x_i)
    ! Surface geopotential: postive if greater than mean seafloor
    implicit none
    type(Coord) :: x_i

    surf_geopot_case = grav_accel * 0d0
  end function surf_geopot_case

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    real(8)            :: lon, lat, rgrc
    real(8), parameter :: n_lat = 4d0, n_lon = 4d0

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    rgrc = radius * acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))
    
    init_free_surface = dH * exp__flush (- abs((lat-lat_c)/(pert_radius/radius))**n_lat &
         - abs((lon-lon_c)/(pert_radius/radius/10d0))**n_lon)

    !init_free_surface = dH * exp__flush (-(rgrc/pert_radius)**6)
  end function init_free_surface

  subroutine init_tke (dom, i, j, zlev, offs, dims)
    ! Initialize TKE
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer  :: d, id

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    tke(zlev)%data(d)%elts(id) = e_min
  end subroutine init_tke

  real(8) function buoyancy_init (x_i, z)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    real(8) :: eps_rho, lat, lon

    call cart2sph (x_i, lon, lat)
    
    if (zlevels == 1) then
       buoyancy_init = 0d0
    elseif (zlevels == 2) then
       if (z >= mixed_layer) then
          buoyancy_init = - drho / ref_density
       else
          buoyancy_init = 0d0
       end if
    elseif (zlevels >= 3) then
       if (z >= mixed_layer) then                           ! constant density perturbation near surface
          buoyancy_init = - drho / ref_density
       elseif (z <= mixed_layer .and. z > thermocline) then ! linear stratification
          buoyancy_init = - (z - thermocline)/(mixed_layer - thermocline) * drho / ref_density
       elseif (z <= thermocline) then                       ! zero density perturbation at depth
          buoyancy_init = 0d0
       end if
    end if
  end function buoyancy_init
  
  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                 :: k
    real(8), dimension(1:N_VARIABLE,1:zmax) :: threshold_new
    character(3), parameter                 :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       call cal_lnorm_sol (sol, order)
       threshold_new = tol*lnorm
       ! Correct for zero velocity case
       do k = 1, zlevels
          if (threshold_new(S_MASS,k) == 0d0) threshold_new(S_MASS,k) = 1d16
          if (threshold_new(S_TEMP,k) == 0d0) threshold_new(S_TEMP,k) = 1d16
          if (threshold_new(S_VELO,k) == 0d0) threshold_new(S_VELO,k) = 1d16
       end do
    end if

    if (istep >= 10) then
       threshold = 0.01d0*threshold_new + 0.99d0*threshold
    else
       threshold = threshold_new
    end if
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer     :: k
    real(8)     :: dz, z
    type(Coord) :: x_i

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0d0

    x_i = Coord (radius, 0d0, 0d0)
    do k = 1, zlevels
       dz = b_vert_mass(k) * max_depth
       z = 0.5d0 * (b_vert(k)+b_vert(k-1)) * max_depth

       ! if (dH == 0d0) then
       !    lnorm(S_MASS,k) = ref_density*dz
       ! else
       !    lnorm(S_MASS,k) = ref_density*dH
       ! end if
       
       lnorm(S_TEMP,k) = ref_density*dz * abs(buoyancy_init (x_i, z))
       if (lnorm(S_TEMP,k) == 0d0) lnorm(S_TEMP,k) = 1d20

       lnorm(S_VELO,k) = Udim
    end do
    lnorm(S_MASS,:) = 1d20

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case
  
  subroutine initialize_dt_viscosity_case 
    ! Not used
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set depth 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i

    id = idx (i, j, offs, dims)
    id_i = id + 1

    ! Set bathymetry
    if (etopo_bathy) then ! set bathymetry coordinates using etopo data
       call etopo_topography (dom, i, j, zlev, offs, dims, 'bathymetry')
    else
       call analytic_topography (dom, i, j, zlev, offs, dims, "bathymetry")
    end if
  end subroutine set_bathymetry

  subroutine set_penal (dom, i, j, zlev, offs, dims)
    ! Set penalization mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i, k

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (penalize) then
       if (etopo_coast) then
          call etopo_topography (dom, i, j, zlev, offs, dims, "penalize")
       else
          call analytic_topography (dom, i, j, zlev, offs, dims, "penalize")
       end if
    else
       penal_node(zlev)%data(d)%elts(id_i)                = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+1:EDGE*id_i) = 0d0   
    end if
  end subroutine set_penal

  subroutine analytic_topography (dom, i, j, zlev, offs, dims, itype)
    ! Only sets penal_node
    use utils_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer     :: d, e, id, id_i, id_e
    real(8)     :: dx, land_radius, lat, lon, mask, r, rgrc, shift, width
    type(Coord) :: p

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    p = dom%node%elts(id_i)
    
    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth + surf_geopot_case (p) / grav_accel
       topography%data(d)%elts(id_i) = max_depth + surf_geopot_case (p) / grav_accel
    case ("penalize")
       call cart2sph (p, lon, lat)

       ! Non-smoothed
       penal_node(zlev)%data(d)%elts(id_i) = drake_land (lon, lat)

       ! Smoothed
       !width = dx_max
       !shift = 4d0*width/radius
       !penal_node(zlev)%data(d)%elts(id_i) = profile2d ()

       ! i = min (int ((lon/DEG + 180d0) * BATHY_PER_DEG) + 1, size(analytic_data,1))
       ! j = min (int ((lat/DEG + 90D0 ) * BATHY_PER_DEG) + 1, size(analytic_data,2))
       ! penal_node(zlev)%data(d)%elts(id_i) = analytic_data(i,j)
    end select
  contains
    real(8) function profile2d ()
      implicit none

      profile2d = profile1d (lat, lat_min+shift, lat_max-shift) * profile1d (lon, lon_min+shift, lon_max-shift)
    end function profile2d

    real(8) function profile1d (x, xmin, xmax)
      implicit none
      real(8) :: x, xmin, xmax

      profile1d = prof (x, xmax) - prof (x, xmin)
    end function profile1d

    real(8) function prof (x, x0)
      implicit none
      real(8) :: x, x0

      prof = 0.5d0 * (1d0 - tanh ((x - x0)/(width/radius)))
    end function prof
  end subroutine analytic_topography

  real(8) function drake_land (lon, lat)
    implicit none
    real(8)           :: lat, lon

    if ((lat<lat_max .and. lat>lat_min) .and. (lon<lon_max .and. lon>lon_min)) then
       drake_land = 1d0
    else
       drake_land = 0d0
    end if
  end function drake_land

  subroutine etopo_topography (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo using etopo data 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer :: d, id, id_i, is0, it0, s, t
    real(8) :: lat, lon, mask, s0, t0

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    ! Find longitude and latitude coordinates
    call cart2sph (dom%node%elts(id_i), lon, lat)
    s0 = lon/DEG * BATHY_PER_DEG
    t0 = lat/DEG * BATHY_PER_DEG
    is0 = nint (s0); it0 = nint (t0)

    mask = topo_value (is0, it0, itype)

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = mask
       !       topography%data(d)%elts(id_i) = mask
    case ("penalize")
       penal_node(zlev)%data(d)%elts(id_i) = mask
    end select
  end subroutine etopo_topography

  real(8) function topo_value (s, t, itype)
    implicit none
    integer      :: s, t
    character(*) :: itype

    select case (itype)
    case ("bathymetry")
       if (etopo_data(s, t) > 0d0) then ! land
          topo_value = min_depth
       else ! sea: topography is less than zero
          topo_value = max (min (etopo_data(s, t), min_depth), max_depth)
       end if
    case ("penalize")
       if (etopo_data(s, t) > 0d0) then ! land
          topo_value = 1d0
       else ! sea: topography is less than zero
          topo_value = 0d0
       end if
    end select
  end function topo_value

  type(Coord) function proj_lon_lat (i, j)
    integer :: i, j
    
    real(8) :: lon, lat

    lon = dble(i) * (BATHY_PER_DEG * DEG)
    lat = dble(j) * (BATHY_PER_DEG * DEG)

    proj_lon_lat = project_on_sphere (sph2cart (lon, lat))
  end function proj_lon_lat 

  subroutine apply_initial_conditions_case
    implicit none
    integer :: d, k, l

    do l = level_start, level_end
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_start, level_end
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_sol,  l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       if (vert_diffuse) then
          do k = 1, zlevels
             call apply_onescale (init_tke,  l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end if
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update means, bathymetry and penalization mask
    use wavelet_mod
    implicit none
    integer :: d, k, l, p

    if (istep /= 0) then
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             call apply_onescale_to_patch (set_bathymetry, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
             do k = 1, zlevels
                call apply_onescale_to_patch (set_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
             end do
          end do
       end do
       call barrier

       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             call apply_onescale_to_patch (init_mean, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    else ! need to set values over entire grid on restart
       do l = level_start, level_end
          call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          do k = 1, zmax
             call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do

       do l = level_start, level_end
          call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
   end if
  end subroutine update_case

  subroutine update_diagnostics
    ! Update diagnostics
    implicit none

  end subroutine update_diagnostics

  subroutine init_diagnostics
    ! Initialize diagnostics
    implicit none 

  end subroutine init_diagnostics

  subroutine set_save_level_case
    implicit none

    save_zlev = zlevels
  end subroutine set_save_level_case

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (zlevels == 1) then
       a_vert(0) = 0d0; a_vert(1) = 1d0
       b_vert(0) = 1d0; b_vert(1) = 0d0
    elseif (zlevels == 2) then 
       a_vert(0) = 0d0; a_vert(1) = 0d0;                 a_vert(2) = 1d0
       b_vert(0) = 1d0; b_vert(1) = mixed_layer/max_depth; b_vert(2) = 0d0
    elseif (zlevels >= 3) then
       if (trim (coords) == "chebyshev") then
          do k = 0, zlevels
             b_vert(k) = (1d0 + cos (dble(k)/dble(zlevels) * MATH_PI)) / 2d0
          end do
       elseif (trim (coords) == "chebyshev_half") then
          do k = 0, zlevels
             b_vert(k) = 1d0 - sin (dble(k)/dble(zlevels) * MATH_PI/2d0)
          end do
       else ! default coordinates are uniform (not used if sigma_z = .true.)
          coords = "uniform"
          do k = 0, zlevels
             b_vert(k) = 1d0 - dble(k)/dble(zlevels)
          end do
       end if
       a_vert = 1d0 - b_vert
    end if

    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

  subroutine dump_case (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump_case

  subroutine load_case (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load_case

  subroutine trend_relax (q, dq)
    ! Trend relaxation to mean buoyancy
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (q, NONE, 27)

    do k = 1, zlevels
       ! Scalars
       do d = 1, size(grid)
          mass   =>  q(S_MASS,k)%data(d)%elts
          temp   =>  q(S_TEMP,k)%data(d)%elts
          mean_m =>  sol_mean(S_MASS,k)%data(d)%elts
          mean_t =>  sol_mean(S_TEMP,k)%data(d)%elts
          dmass  => dq(S_MASS,k)%data(d)%elts
          dtemp  => dq(S_TEMP,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
          end do
          nullify (dmass, dscalar, mass, temp, mean_m, mean_t)
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

    dmass(id_i) = 0d0
    dtemp(id_i) = - k_T * (temp(id_i) - mass(id_i) * mean_t(id_i)/mean_m(id_i))
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims)

    dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
  end subroutine trend_velo

  function z_coords_case (eta_surf, z_s)
    ! Hybrid sigma-z vertical coordinates to minimize inclination of layers to geopotential
    ! near the free surface over strong bathymetry gradients.
    ! Reference: similar to Shchepetkin and McWilliams (JCP vol 228, 8985-9000, 2009)
    !
    ! Sets the a_vert parameter that depends on eta_surf (but not b_vert).
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    integer                       :: k
    real(8)                       :: cff, cff1, cff2, hc, z_0
    real(8), dimension(0:zlevels) :: Cs, sc

    real(8), parameter            :: theta_b = 0d0, theta_s = 7d0
    real(8), parameter            :: hc_min = -200d0 * METRE ! minimum depth of uniform layer region
    
    hc = abs (min (mixed_layer, hc_min)) 
    
    cff1 = 1d0 / sinh (theta_s)
    cff2 = 0.5d0 / tanh (0.5d0 * theta_s)
    
    sc(0) = -1d0
    Cs(0) = -1d0
    cff = 1d0 / dble(zlevels)
    do k = 1, zlevels
       sc(k) = cff * dble (k - zlevels)
       Cs(k) = (1d0 - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5d0)) - 0.5d0)
    end do

    z_coords_case(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1d0 - z_0 / z_s
       z_coords_case(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords_case

  function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type)
    ! Additional physics for the flux term of the scalar trend
    ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
    !
    ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
    type(domain)                                         :: dom
    integer                                              :: d, id, idE, idNE, idN, v, zlev
    logical, optional                                    :: type

    integer                    :: id_i
    real(8)                    :: visc
    real(8), dimension(1:EDGE) :: d_e, grad, l_e
    logical                    :: local_type

    if (present(type)) then
       local_type = type
    else
       local_type = .false.
    end if

    id_i = id + 1
    d = dom%id + 1

    if (Laplace_order == 0 .or. maxval (C_visc(scalars(1):scalars(2))) == 0d0) then
       physics_scalar_flux_case = 0d0
    else
       if (.not.local_type) then ! usual flux at edges E, NE, N
          l_e =  dom%pedlen%elts(EDGE*id+1:EDGE*id_i)
          d_e =  dom%len%elts(EDGE*id+1:EDGE*id_i)
       else ! flux at SW corner
          l_e(RT+1) = dom%pedlen%elts(EDGE*idE+RT+1)
          l_e(DG+1) = dom%pedlen%elts(EDGE*idNE+DG+1)
          l_e(UP+1) = dom%pedlen%elts(EDGE*idN+UP+1)
          d_e(RT+1) = -dom%len%elts(EDGE*idE+RT+1)
          d_e(DG+1) = -dom%len%elts(EDGE*idNE+DG+1)
          d_e(UP+1) = -dom%len%elts(EDGE*idN+UP+1)
       end if

       ! Calculate gradients
       if (Laplace_order == 1) then
          grad = grad_physics (q(v,zlev)%data(d)%elts)
       elseif (Laplace_order == 2) then
          grad = grad_physics (Laplacian_scalar(v)%data(d)%elts)
       end if

       ! Scale aware viscosity
       visc =  C_visc(v) * dom%len%elts(EDGE*id+RT+1)**(2d0*Laplace_order_init)/dt

       ! Flux on lhs of scalar equation (hence negative)
       physics_scalar_flux_case = visc * (-1d0)**Laplace_order * grad * l_e
    end if
  contains
    function grad_physics (scalar)
      implicit none
      real(8), dimension(1:EDGE) :: grad_physics
      real(8), dimension(:)      :: scalar

      grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))   / d_e(RT+1)
      grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1)) / d_e(DG+1)
      grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))   / d_e(UP+1)
    end function grad_physics
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    ! Additional physics for the source term of the velocity trend
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id, id_i, idE, idN, idNE
    real(8)                    :: lat, lon, visc
    real(8), dimension(1:EDGE) :: horiz_diffusion, h1, h2, penal, u1, u2, vert_diffusion

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    
    idE  = idx (i+1, j,   offs, dims) 
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    ! Scale aware viscosity
    visc = C_visc(S_VELO) * dom%len%elts(EDGE*id+RT+1)**(2d0*Laplace_order)/dt

    ! Only diffuse rotu
    horiz_diffusion = visc * (-1d0)**Laplace_order * (curl_rot () - grad_div ())

    ! Vertical diffusion
    if (vert_diffuse) then ! using vertical diffusion module
       vert_diffusion = 0d0
    else
       ! Layer thicknesses and velocities
       if (zlevels == 2) then
          h1 = dz_e (dom, i, j, 1, offs, dims, sol)
          h2 = dz_e (dom, i, j, 2, offs, dims, sol)

          u1 = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
          u2 = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
       elseif (zlevels == 1) then
          h1 = dz_e (dom, i, j, 1, offs, dims, sol);                  h2 = h1
          u1 = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1); u2 = u1
       end if

       ! Vertical diffusion
       if (zlevels == 2) then
          if (zlev == 1) then
             vert_diffusion = - Ku / (h1 * (h1 + h2)/2d0) * (u1 - u2) + bottom_drag ()
          elseif (zlev == 2) then
             vert_diffusion = - Ku / (h2 * (h1 + h2)/2d0) * (u2 - u1) + wind_drag ()
          end if
       elseif (zlevels == 1) then
          vert_diffusion = bottom_drag () + wind_drag ()
       end if
    end if

    ! Penalization (non-smoothed)
    call cart2sph (dom%node%elts(id_i), lon, lat)
    penal = - drake_land (lon, lat)/dt * velo(EDGE*id+RT+1:EDGE*id+UP+1)

    physics_velo_source_case = horiz_diffusion + vert_diffusion + penal
  contains
    function bottom_drag ()
      implicit none
      real(8), dimension(3) :: bottom_drag

      bottom_drag = - bottom_friction * u1 / h1
    end function bottom_drag
    
    function wind_drag ()
      implicit none
      real(8), dimension(3) :: wind_drag
      
      real(8), dimension(3) :: tau_wind

      tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id_i),   dom%node%elts(idE+1))
      tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE+1), dom%node%elts(id_i))
      tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id_i),   dom%node%elts(idN+1))

      wind_drag = tau_wind / (ref_density * h2)
    end function wind_drag
    
    function grad_div ()
      implicit none
      real(8), dimension(3) :: grad_div

      grad_div(RT+1) = (divu(idE+1) - divu(id+1))   / dom%len%elts(EDGE*id+RT+1)
      grad_div(DG+1) = (divu(id+1)  - divu(idNE+1)) / dom%len%elts(EDGE*id+DG+1)
      grad_div(UP+1) = (divu(idN+1) - divu(id+1))   / dom%len%elts(EDGE*id+UP+1)
    end function grad_div

    function curl_rot ()
      implicit none
      real(8), dimension(3) :: curl_rot

      integer :: idS, idW

      idS = idx (i,   j-1, offs, dims)
      idW = idx (i-1, j,   offs, dims)

      curl_rot(RT+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)
      curl_rot(DG+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
      curl_rot(UP+1) = (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
    end function curl_rot
  end function physics_velo_source_case

  real(8) function bottom_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Bottom boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    bottom_buoy_flux_case = 0d0
  end function bottom_buoy_flux_case

  real(8) function top_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Top boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    top_buoy_flux_case = 0d0
  end function top_buoy_flux_case

  function wind_flux_case (dom, i, j, zlev, offs, dims)
    ! Wind stress velocity source term evaluated at edges (top boundary condition for vertical diffusion of velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: wind_flux_case

    integer                    :: d, id, idE, idN, idNE
    real(8)                    :: rho
    real(8), dimension(1:EDGE) :: tau_wind

    id = idx (i, j, offs, dims)

    if (maxval (dom%mask_e%elts(EDGE*id+RT+1:EDGE*id+UP+1)) >= ADJZONE) then
       d = dom%id + 1
       idE  = idx (i+1, j,   offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idNE = idx (i+1, j+1, offs, dims)

       tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idE+1))
       tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE+1), dom%node%elts(id+1))
       tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idN+1))

       rho = porous_density (d, id+1, zlevels)

       wind_flux_case = tau_wind / rho  * (1d0 - penal_edge(zlevels)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
    else
       wind_flux_case = 0d0
    end if
  end function wind_flux_case

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    ! Idealized zonally and temporally averaged zonal and meridional wind stresses
    ! (based on Figure 4 from Ferreira et al J Climate 24, 992-1012 (2011) and Figure 4 from Gille J Atmos Ocean Tech 22, 1353-1372 respectively)
    implicit none
    real(8) :: lat, lon, peak, tau_zonal, tau_merid

    logical, parameter :: merid_stress = .false.

    peak = (abs(lat)*180d0/MATH_PI - 35d0) / 20d0
    tau_zonal = -tau_0 * 1.2d0 * exp (-peak**2) * sin (abs(lat)*6d0) - 5d-3*exp(-(lat*180d0/MATH_PI/10d0)**2)

    if (merid_stress) then
       peak = lat*180d0/MATH_PI / 15d0
       tau_merid = -tau_0 * 2.5d0 * exp (-peak**2) * sin (2d0*lat) * peak**2
    else
       tau_merid = 0d0
    end if
  end subroutine wind_stress

  real(8) function tau_mag_case (p)
    ! Magnitude of wind stress at node p
    implicit none
    type(Coord) :: p

    real(8) :: lat, lon, tau_zonal, tau_merid

    call cart2sph (p, lon, lat)

    call wind_stress (lon, lat, tau_zonal, tau_merid)

    tau_mag_case = sqrt (tau_zonal**2 + tau_merid**2)
  end function tau_mag_case
end module test_case_mod
