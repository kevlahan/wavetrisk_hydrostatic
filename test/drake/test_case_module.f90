Module test_case_mod
  ! Module file for Drake passage test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, etopo_res, resume_init
  real(8)                              :: dt_cfl, k_T, total_cpu_time
  real(8)                              :: Hdim, Ldim, Tdim, Udim

  ! Local variables
  real(8)                              :: beta, bv
  real(8)                              :: delta_I, delta_M, delta_S, delta_sm, drho, drho_dz, f0, Ku, mixed_layer, Rb, Rd
  real(8)                              :: Rey, Ro, radius_earth, omega_earth, scale, scale_omega, halocline, npts_penal, u_wbc 
  real(8)                              :: resolution, tau_0
  real(8),                      target :: bottom_friction_case  
  real(4), allocatable, dimension(:,:) :: topo_data
  logical                              :: aligned, etopo_coast
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

    bottom_friction  => bottom_friction_case
  end subroutine assign_functions

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
    real(8), dimension(1:EDGE) :: d_e, grad, l_e
    logical                    :: local_type

    if (present(type)) then
       local_type = type
    else
       local_type = .false.
    end if

    id_i = id + 1
    d = dom%id + 1

    if (Laplace_order == 0 .or. maxval (visc_sclr) == 0d0) then
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

       ! Complete scalar diffusion
       physics_scalar_flux_case = (-1)**Laplace_order * visc_sclr(v) * grad * l_e
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
    real(8)                    :: dx, visc
    real(8), dimension(1:EDGE) :: horiz_diffusion, h1, h2, u1, u2, vert_diffusion

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    
    idE  = idx (i+1, j,   offs, dims) 
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    ! Layer thicknesses and velocities
    if (zlevels == 2) then
       h1 = dz_e (dom, i, j, 1, offs, dims)
       h2 = dz_e (dom, i, j, 2, offs, dims)
       
       u1 = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
       u2 = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    elseif (zlevels == 1) then
        h1 = dz_e (dom, i, j, 1, offs, dims);                       h2 = h1
        u1 = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1); u2 = u1
    end if

    ! Horizontal diffusion
    if ((dom%node%elts(id_i)%x**2 + dom%node%elts(id_i)%y**2)/(4d0*dx_max)**2 < 1d0) then ! increase diffusion near poles to remove noise at these lower accuracy points
       dx = sqrt (4d0/sqrt(3d0) * 4d0*MATH_PI*radius**2/(20d0*4**level_end)) 
       visc = dx**2/dt/32d0
       horiz_diffusion = (-1d0)**(Laplace_order-1) * visc * (grad_divu() - curl_rotu())
    else
       horiz_diffusion = (-1d0)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
    end if
    
    ! Vertical diffusion
    if (zlevels == 2) then
       if (zlev == 1) then
          vert_diffusion = - Ku / (h1 * (h1 + h2)/2d0) * (u1 - u2) + bottom_drag()
       elseif (zlev == 2) then
          vert_diffusion = - Ku / (h2 * (h1 + h2)/2d0) * (u2 - u1) + wind_drag()
       end if
    elseif (zlevels == 1) then
       vert_diffusion = bottom_drag() + wind_drag()
    end if

    ! Complete source term for velocity trend (do not include drag and wind stress in solid regions)
    if (penal_node(zlevels)%data(d)%elts(id_i) < 1d-3) then
       physics_velo_source_case = horiz_diffusion + vert_diffusion
    else
       physics_velo_source_case = horiz_diffusion
    end if
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
    
    function grad_divu ()
      implicit none
      real(8), dimension(3) :: grad_divu

      grad_divu(RT+1) = (divu(idE+1) - divu(id+1))   / dom%len%elts(EDGE*id+RT+1)
      grad_divu(DG+1) = (divu(id+1)  - divu(idNE+1)) / dom%len%elts(EDGE*id+DG+1)
      grad_divu(UP+1) = (divu(idN+1) - divu(id+1))   / dom%len%elts(EDGE*id+UP+1)
    end function grad_divu

    function curl_rotu()
      implicit none
      real(8), dimension(3) :: curl_rotu

      integer :: idS, idW

      idS = idx (i,   j-1, offs, dims)
      idW = idx (i-1, j,   offs, dims)

      curl_rotu(RT+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)
      curl_rotu(DG+1) = (vort(TRIAG*id +LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
      curl_rotu(UP+1) = (vort(TRIAG*idW+LORT+1) - vort(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
    end function curl_rotu
  end function physics_velo_source_case

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
       write (6,'(A,L1)')     "default_thresholds             = ", default_thresholds
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
       write (6,'(A,i3)')     "etopo_res                      = ", etopo_res
       write (6,'(A,es10.4)') "resolution                     = ", resolution

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
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay time    [d]  = ", abs(max_depth)/bottom_friction_case / DAY
       if (zlevels == 2)  write (6,'(A,es11.4)') "Ku                    [m^2/s]  = ", Ku
       write (6,'(A,es11.4)') "buoyancy relaxation       [d]  = ", 1d0/k_T / DAY
       write (6,'(A,es11.4)') "f0 at 45 deg          [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 45 deg       [rad/ms]  = ", beta
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(A,es11.4)') "Inertial layer           [km]  = ", delta_I  / KM
       write (6,'(A,es11.4)') "Munk layer               [km]  = ", delta_M  / KM
       write (6,'(A,es11.4)') "Stommel layer            [km]  = ", delta_S  / KM
       write (6,'(A,es11.4)') "submesoscale             [km]  = ", delta_sm / KM
       write (6,'(A,es11.4)') "barotropic Rossby radius [km]  = ", Rd / KM
       write (6,'(A,es11.4,/)') "baroclinic Rossby radius [km]  = ", Rb / KM
       write (6,'(A,es11.4)') "Rossby number                  = ", Ro
       write (6,'(A,es11.4)') "Re (delta_I u_wbc / nu)        = ", Rey 
       write (6,'(A,es11.4)') "Resolution of Munk layer       = ", delta_M / dx_min
       write (6,'(A,es11.4)') "Resolution of Taylor scale     = ", delta_I / sqrt(Rey) / dx_min
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

  subroutine apply_initial_conditions_case
    implicit none
    integer :: d, k, l

    call topo_drake_data

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
  end subroutine apply_initial_conditions_case

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initial perturbation to mean 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: d, id, id_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    if (zlev == zlevels+1) then ! 2D barotropic mode
       sol(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
    else ! 3D layers
       sol(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
    end if
    ! Set initial velocity field to zero
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id:EDGE*id_i) = 0d0
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id, id_i
    real (8)    :: dz, rho, z
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)

    if (zlev == zlevels+1) then
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
    else
       dz =  b_vert_mass(zlev) * dom%topo%elts(id_i)
       z = 0.5d0 * (b_vert(zlev) + b_vert(zlev-1)) * dom%topo%elts(id_i)

       rho = porous_density (d, id_i, zlev)
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = rho * dz
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy_init (x_i, z)
    end if
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
  end subroutine init_mean

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

       do k = 1, zlevels
          do d = 1, size(grid)
             do p = n_patch_old(d)+1, grid(d)%patch%length
                call apply_onescale_to_patch (init_mean, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
             end do
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
          do k = 1, zmax
             call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end if
  end subroutine update_case

  real(8) function surf_geopot_case (x_i)
    ! Surface geopotential: postive if greater than mean seafloor
    ! MUST BE SET EQUAL TO ZERO FOR THIS test case
    implicit none
    type(Coord) :: x_i

    surf_geopot_case = grav_accel * 0d0
  end function surf_geopot_case

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0d0
  end function init_free_surface

  real(8) function buoyancy_init (x_i, z)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    if (zlevels /= 1 .and. z >= halocline) then
       buoyancy_init = - (1d0 - z/halocline) * drho/ref_density
    else
       buoyancy_init = 0d0
    end if
  end function buoyancy_init

  subroutine print_density_pert
    implicit none
    integer     :: k
    real(8)     :: z
    type(Coord) :: x_i 

    x_i = Coord (radius, 0d0, 0d0)

    write (6,'(a)') " Layer    z       drho"      
    do k = 1, zlevels
       z = 0.5d0 * (b_vert(k)+b_vert(k-1)) * max_depth
       write (6, '(2x,i2, 1x, 2(es9.2,1x))') k, z, - buoyancy_init (x_i, z)*ref_density
    end do
    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density_pert

  subroutine set_thresholds_case
    ! Set thresholds dynamically
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                 :: k, v
    real(8), dimension(1:N_VARIABLE,1:zmax) :: threshold_new
    character(3), parameter                 :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       call cal_lnorm_sol (sol, order)
       threshold_new = tol * lnorm
       ! Correct very small values
       do k = 1, zmax
          if (threshold_new(S_MASS,k) < threshold_def(S_MASS,k)/1d1) threshold_new(S_MASS,k) = threshold_def(S_MASS,k)
          if (threshold_new(S_TEMP,k) < threshold_def(S_TEMP,k)/1d1) threshold_new(S_TEMP,k) = threshold_def(S_TEMP,k)
          if (threshold_new(S_VELO,k) < threshold_def(S_VELO,k)/1d1) threshold_new(S_VELO,k) = threshold_def(S_VELO,k)
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

       lnorm(S_MASS,k) = ref_density*dz

       lnorm(S_TEMP,k) = ref_density*dz * abs(buoyancy_init (x_i, z))
       if (lnorm(S_TEMP,k) == 0d0) lnorm(S_TEMP,k) = 1d16

       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8)            :: area, C_divu, C_sclr, C_rotu, C_visc, tau_divu, tau_rotu, tau_sclr
    logical, parameter :: munk = .true.

    area = 4d0*MATH_PI*radius**2/(20d0*4d0**max_level) ! average area of a triangle
    dx_min = sqrt (4d0/sqrt(3d0) * area)               ! edge length of average triangle

    area = 4d0*MATH_PI*radius**2/(20d0*4d0**min_level)
    dx_max = sqrt (4d0/sqrt(3d0) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, 1.4d0*dx_min/u_wbc, dx_min/c1)
    dt_init = dt_cfl

    ! Diffusion constants
    if (munk) then
       C_visc = dt_cfl * beta * dx_min * resolution**(2d0*Laplace_order_init+1)                ! resolve Munk layer with "resolution" points
    else
       C_visc = resolution**2 * dx_min**(2*(1-Laplace_order_init)) * dt_cfl * u_wbc / delta_I  ! resolve Taylor scale with "resolution" points
    end if

    ! Ensure stability
    C_visc = min ((1d0/32d0)**Laplace_order_init, C_visc)

    C_visc = 5d-3
    
    C_rotu = C_visc
    C_divu = C_visc
    C_sclr = C_visc

    ! Diffusion time scales
    tau_sclr = dt_cfl / C_sclr
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0d0
       visc_divu = 0d0
       visc_rotu = 0d0
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr(S_MASS) = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_sclr(S_TEMP) = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_rotu = dx_min**(2d0*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2d0*Laplace_order_init) / tau_divu
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

    Laplace_order = Laplace_order_init
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_drake (dom, i, j, zlev, offs, dims, 'bathymetry')
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
       call topo_drake (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0       
    end if
  end subroutine set_penal

  subroutine topo_drake (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer            :: d, e, id, id_e, id_i
    real(8)            :: dx, lat, lat0, lat_width, lon, mask, n_lat, n_lon
    real(8), parameter :: lat_max = 60, lat_min = -35, lon_width = 15
    type(Coord)        :: p

    id = idx (i, j, offs, dims)
    id_i = id + 1

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth + surf_geopot_case (dom%node%elts(id_i)) / grav_accel
    case ("penalize")
       call cart2sph (dom%node%elts(id_i), lon, lat)
       dx = dx_max

       ! Analytic land mass with smoothing
       lat_width = (lat_max - lat_min) / 2d0
       lat0 = lat_max - lat_width

       n_lat = 4d0*radius * lat_width*DEG / (dx * npts_penal)
       n_lon = 4d0*radius * lon_width*DEG / (dx * npts_penal)

       mask = exp__flush (- abs((lat/DEG-lat0)/lat_width)**n_lat - abs(lon/DEG/(lon_width))**n_lon) ! constant longitude width

       d = dom%id + 1
       penal_node(zlev)%data(d)%elts(id_i) = mask
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = max (penal_edge(zlev)%data(d)%elts(id_e), mask)
       end do
    end select
  end subroutine topo_drake

  subroutine topo_drake_data
    ! Defines analytic latitude-longitude topography data for Drake passage case
    implicit none
    integer                              :: ii, ilat, ilon, lat_avg, lon_avg, jj, kk
    integer, parameter                   :: npts = 2, n_smooth = 20
    real(8)                              :: avg, lat, lon
    real(8), parameter                   :: width = 30d0, lat_max = 70d0, lat_min = -35d0
    real(8), dimension(2)                :: sz
    real(8), dimension(:,:), allocatable :: temp_data

    if (etopo_coast) then
       if (rank == 0) write(6,'(a)') 'Reading bathymetry data'
       bathy_per_deg = 60d0/etopo_res

       allocate (topo_data(-180*bathy_per_deg:180*bathy_per_deg, -90*bathy_per_deg:90*bathy_per_deg))
       open (unit=1086,file='bathymetry') ! "bathymetry" is symbolic link to appropriate etopo bathymetry data
       do kk = ubound (topo_data,2), lbound (topo_data,2), -1 ! north to south (as read from file)
          read (1086,*) topo_data(:,kk)
       end do
       close (1086)
    end if
  end subroutine topo_drake_data

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

  subroutine set_save_level_case
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5d0 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level_case

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (zlevels == 2) then 
       a_vert(0) = 0d0; a_vert(1) = 0d0;                 a_vert(2) = 1d0
       b_vert(0) = 1d0; b_vert(1) = halocline/max_depth; b_vert(2) = 0d0
    elseif (zlevels == 3) then
       a_vert(0) = 0d0; a_vert(1) = 0d0;                 a_vert(2) = 0d0;                   a_vert(3) = 1d0 
       b_vert(0) = 1d0; b_vert(1) = halocline/max_depth; b_vert(2) = mixed_layer/max_depth; b_vert(3) = 0d0
    else ! uniform sigma grid: z = a_vert*eta + b_vert*z_s
       do k = 0, zlevels
          a_vert(k) = dble(k)/dble(zlevels)
          b_vert(k) = 1d0 - dble(k)/dble(zlevels)
       end do
    end if

    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

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

    dmass(id_i) = 0d0
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

    dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
  end subroutine trend_velo

  function z_coords_case (eta_surf, z_s)
    ! Dummy routine
    ! (see upwelling test case for example)
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    z_coords_case = 0d0
  end function z_coords_case
end module test_case_mod
