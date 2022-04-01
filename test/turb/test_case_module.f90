Module test_case_mod
  ! Module file for homogeneous test case
  use shared_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use equation_of_state_mod
  use projection_mod
  use spline_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, resume_init
  real(8)                              :: dt_cfl, tau_diffusion, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim

  ! Local variables
  real(8)                              :: bv, drho, drho_dz, grav_reduced, tau_nudge, Tcline
  real(8)                              :: r_max, r_max_loc, tau_0
  real(8), target                      :: bottom_friction_case
  character(255)                       :: coords

  ! 2D projection
  integer                                :: Nproj
  real(8), dimension(:,:,:), allocatable :: y2, y2_0, zonal, zonal_0
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

    if (implicit_diff_sclr .or. Laplace_order == 0 .or. maxval (visc_sclr) == 0d0) then
       physics_scalar_flux_case = 0.0_8
    else
       id_i = id + 1
       d = dom%id + 1

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

  function physics_scalar_source (q, id, zlev)
    ! Additional physics for the source term of the scalar trend
    use domain_mod
    implicit none
    real(8), dimension(scalars(1):scalars(2))            :: physics_scalar_source
    integer                                              :: id, zlev
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q

    physics_scalar_source = 0d0
  end function physics_scalar_source

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    ! Additional physics for the source term of the velocity trend
    ! wind stress and bottom friction are included as surface fluxes in the split eddy viscosity split step
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id
    real(8), dimension(1:EDGE) :: diffusion

    id = idx (i, j, offs, dims)

    ! Horizontal diffusion
    if (implicit_diff_divu) then
       diffusion = - visc_rotu * curl_rotu()
    else
       diffusion = (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
    end if

    physics_velo_source_case = diffusion 
  contains
    function grad_divu()
      implicit none
      real(8), dimension(3) :: grad_divu

      integer :: idE, idN, idNE

      idE  = idx (i+1, j,   offs, dims)
      idN  = idx (i,   j+1, offs, dims)
      idNE = idx (i+1, j+1, offs, dims)

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
    end if
#ifdef MPI
    call MPI_Bcast (test_case, 255, MPI_BYTE,             0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (run_id,    255, MPI_BYTE,             0, MPI_COMM_WORLD, ierror)    
    call MPI_Bcast (max_level,   1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (zlevels,     1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (no_slip,     1, MPI_LOGICAL,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (remap,       1, MPI_LOGICAL,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (iremap,      1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (log_iter,    1, MPI_LOGICAL,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (tol,         1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (cfl_num,     1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (dt_write,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (CP_EVERY,    1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (time_end,    1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast (resume_init, 1, MPI_INTEGER,          0, MPI_COMM_WORLD, ierror)
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
       write (6,'(A,L1)')     "penalize                       = ", penalize
       write (6,'(A,L1)')     "no_slip                        = ", no_slip
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
       write (6,'(a,i5)')     "iadapt                         = ", iadapt
       write (6,'(A,i1)')     "optimize_grid                  = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt                       = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num                        = ", cfl_num
       write (6,'(a,a)')      "timeint_type                   = ", trim (timeint_type)
       write (6,'(A,i1)')     "Laplace_order                  = ", Laplace_order_init
       write (6,'(A,i1)')     "n_diffuse                      = ", n_diffuse
       write (6,'(A,L1)')     "vert_diffuse                   = ", vert_diffuse
       write (6,'(A,L1)')     "tke_closure                    = ", tke_closure
       write (6,'(A,es10.4)') "dt_write [d]                   = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY                       = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance                      = ", rebalance
       write (6,'(A,es10.4)') "time_end [d]                   = ", time_end/DAY
       write (6,'(A,i6)')     "resume                         = ", resume_init

       write (6,'(/,A)')      "TIME INTEGRATION PARAMETERS"
       write (6,'(A,L1)')   "mode_split                      = ", mode_split
       write (6,'(A,F4.2)') "theta1                          = ", theta1
       write (6,'(A,F4.2)') "theta2                          = ", theta2

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius                   [km]  = ", radius / KM
       write (6,'(A,es11.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", max_depth
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", min_depth
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed           [m/s]  = ", c1
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1d0/bottom_friction_case / DAY
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(a,es11.4)') "r_max                          = ", r_max
       write (6,'(a,es11.4)') "tau_nudge                [d]   = ", tau_nudge / DAY
       write (6,'(A,i4)')     "Nproj                          = ", Nproj
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

  subroutine apply_initial_conditions_case
    use ops_mod
    implicit none
    integer :: d, k, l, p

    do l = level_end, level_start, -1
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_end, level_start, -1
       call apply_onescale (init_mean,    l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_scalars, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do

    ! Initial velocity is given by thermal wind geostrophic balance with density
    do k = 1, zlevels
       call thermal_wind (k)
       do l = level_end, level_start, -1
          call apply_onescale (init_velo, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do
  end subroutine apply_initial_conditions_case

  subroutine thermal_wind (k)
    ! Computes thermal wind geostrophic balance based on initial density at vertical level k
    ! by upwards integration assuming zero velocity at bathymetry
    ! results are stored in u_zonal, v_merid
    implicit none
    integer :: k

    integer :: d, j, l

    do l = level_end, level_start, -1
       do d = 1, size (grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (thermal_wind_integration, grid(d), grid(d)%lev(l)%elts(j), k, 0, 0)
          end do
       end do
       call update_bdry (sol_mean(S_VELO,k), l, 10)

       do d = 1, size (grid)
          velo  => sol_mean(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_edge_node,    grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             call apply_onescale_to_patch (geostrophic_balance, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
    end do
  end subroutine thermal_wind

  subroutine thermal_wind_integration (dom, i, j, zlev, offs, dims)
    ! Integrates thermal wind geostrophic equations upwards with zero velocity boundary condition
    use ops_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id, idE, idN, idNE
    real(8), dimension(1:EDGE) :: f

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    f = f_coriolis_edge (dom, i, j, zlev, offs, dims)
    if (zlev > 1) then
       sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = &
            sol_mean(S_VELO,zlev-1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) + 0.5d0 * (increment(zlev) + increment(zlev-1)) 
    else
       sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.5d0 * increment(zlev)
    end if
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
         * (1d0 - penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
  contains
    function increment (k)
      implicit none
      integer                    :: k
      real(8), dimension(1:EDGE) :: increment

      real(8), dimension(1:EDGE) :: drho, dz, rho_0
      real(8), dimension(0:EDGE) :: rho

      dz = dz_e (dom, i, j, k, offs, dims)

      rho_0 = porous_density_edge (d, id+1, k)

      rho(0)    = rho_i (i,   j,   k)
      rho(RT+1) = rho_i (i+1, j,   k)
      rho(DG+1) = rho_i (i+1, j+1, k)
      rho(UP+1) = rho_i (i,   j+1, k)

      drho(RT+1) = (rho(1) - rho(0)) / dom%len%elts(EDGE*id+RT+1)
      drho(DG+1) = (rho(0) - rho(2)) / dom%len%elts(EDGE*id+DG+1)
      drho(UP+1) = (rho(3) - rho(0)) / dom%len%elts(EDGE*id+UP+1)

      increment = drho * grav_accel / (rho_0 * f) * dz
    end function increment

    real(8) function rho_i (i, j, k)
      ! Initial density at node corresponding to coordinates (i,j)
      implicit none
      integer :: i, j, k

      integer :: id
      real(8) :: lat, lon, z

      id = idx (i, j, offs, dims)

      call cart2sph (dom%node%elts(id+1), lon, lat)
      z = z_i (dom, i, j, k, offs, dims)

      rho_i = density_init (lat/DEG, z)
    end function rho_i
  end subroutine thermal_wind_integration

  subroutine geostrophic_balance (dom, i, j, zlev, offs, dims)
    ! Geostrophic balance from thermal wind relation
    ! d_z u = g/(rho_0 f) (d_y rho, -d_x rho)
    ! (rhs provided in dom%u_zonal, dom%v_merid)

    use ops_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    real(8) :: u_z, v_m

    id = idx (i, j, offs, dims)

    u_z = dom%u_zonal%elts(id+1)
    v_m = dom%v_merid%elts(id+1)

    dom%u_zonal%elts(id+1) =   v_m
    dom%v_merid%elts(id+1) = - u_z
  end subroutine geostrophic_balance

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

  subroutine init_scalars (dom, i, j, zlev, offs, dims)
    ! Initial scalars fors entire vertical column 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: density, eta, lat, lon, phi, rho, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord)                   :: p 

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    p = dom%node%elts(id_i)
    call cart2sph (p, lon, lat)

    eta = init_free_surface (dom%node%elts(id_i))
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
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = 0d0
       end if
       sol(S_TEMP,k)%data(d)%elts(id_i) = rho * dz(k) * buoyancy_init (lat/DEG, z_k)
    end do

    if (mode_split) then
       phi = phi_node (d, id_i, zlevels)
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i) = phi * eta ! free surface perturbation
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0d0
    end if
  end subroutine init_scalars

  subroutine init_velo (dom, i, j, zlev, offs, dims, vel_fun)
    ! Sets the velocities on the computational grid from zonal and meridional velocities dom%u_zonal and dom%v_merid
    ! (also sets sol_mean to be used in nudging)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    external                        :: vel_fun

    integer :: d, id

    d  = dom%id+1
    id = idx (i, j,offs, dims)

    call interp_node_edge (dom, i, j, z_null, offs, dims, sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
  end subroutine init_velo

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

    eta = 0d0
    z_s = dom%topo%elts(id_i)

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (d, id_i, k)

       sol_mean(S_MASS,k)%data(d)%elts(id_i)                      = rho * dz(k)
       sol_mean(S_TEMP,k)%data(d)%elts(id_i)                      = 0d0
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do

    if (mode_split) then
       sol_mean(S_MASS,zlevels+1)%data(d)%elts(id_i) = 0d0
       sol_mean(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0d0
       sol_mean(S_VELO,zlevels+1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end if
  end subroutine init_mean

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid 
    ! (a_vert, b_vert not used if sigma_z = .true.)
    implicit none
    integer               :: k
    real(8)               :: z
    real(8), dimension(6) :: p

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    b_vert(0) = 1.0_8 ; b_vert(zlevels) = 0d0
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
  end subroutine initialize_a_b_vert_case

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

  real(8) function surf_geopot_case (p)
    ! Surface geopotential: postive if greater than mean seafloor                                                                                        
    implicit none
    type(Coord) :: p

    real(8) :: lat, lon

    call cart2sph (p, lon, lat)

    surf_geopot_case = grav_accel * 0
  end function surf_geopot_case

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0d0
  end function init_free_surface

  real(8) function buoyancy_init (lat, z)
    ! Initial buoyancy at depth z at latitude lat (in degrees)
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8) :: lat, z

    buoyancy_init = (ref_density - density_init (lat, z)) / ref_density
  end function buoyancy_init

  real(8) function density_init (lat, z)
    ! Initial density given latitude (in degrees) and depth (in metres)
    ! Initials conditions are a set of 17 baroclinic jets in geostrophic (thermal wind) balance.
    implicit none
    real(8) :: lat, z

    integer                          :: i
    real(8)                          :: amp, amp0, density_pert, jet_depth
    real(8),               parameter :: jet_width = 10d0                                              ! width of density jet transition regions (in degrees)
    real(8), dimension(8), parameter :: lat0 = (/ 10d0, 20d0, 30d0, 40d0, 50d0, 60d0,  70d0,  80d0 /) ! latitudes of density jets

    ! Amplitudes of density jet density perturbations to give a maximum velocity Udim in thermal wind balance (assume drho = -2 kg/m^3)
    amp0 = 2.30d-2 * Udim
    amp  = 1.25d+0 * Udim

    ! Penetration depth of density perturbations
    jet_depth = max_depth / 4d0  

    density_pert = density_jet (0d0) ! equatorial jet
    do i = 1, size (lat0)
       density_pert = density_pert + density_jet (lat0(i)) + density_jet (-lat0(i))
    end do
    density_init = ref_density + density_pert
  contains
    real(8) function density_jet (lat0)
      ! Gaussian density profile in latitude and depth
      implicit none
      real(8) :: lat0

      density_jet = (amp0 + amp * abs(sin(lat0*DEG))) * drho &
           * exp (- ((lat - lat0) / (0.25d0*jet_width))**2 - (z/(0.6d0*jet_depth))**2)
    end function density_jet
  end function density_init

  subroutine print_density
    implicit none
    integer                       :: k
    real(8)                       :: bv, c_k, c1, drho, dz_l, eta, lat, rho, rho_above, z_k, z_s, z_above
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    lat = 10d0 ! latitude to evaluate buoyancy

    eta = 0d0
    z_s = max_depth

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    write (6,'(a)') " Layer      z         dz         rho "
    do k = 1, zlevels
       z_k = interp (z(k-1), z(k))
       write (6, '(2x, i2, 4x, 2(es9.2, 1x), es11.5)') &
            k, z_k, dz(k), ref_density * (1d0 - buoyancy_init (lat, z_k))
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

       rho_above = ref_density * (1d0 - buoyancy_init (lat, z_above))
       rho  = ref_density * (1d0 - buoyancy_init (lat, z_k))
       drho = rho_above - rho

       bv = sqrt(- grav_accel * drho/dz_l/rho)
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

  subroutine set_thresholds_case
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
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer                       :: k
    real(8)                       :: eta, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord)                   :: x_i

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0d0

    eta = 0d0
    z_s = max_depth

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = eta * a_vert + z_s * b_vert
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       lnorm(S_MASS,k) = ref_density * dz(k)
       lnorm(S_TEMP,k) = abs (drho) * dz(k)
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C, C_b, C_divu, C_mu, C_rotu, C_visc, dlat, tau_b, tau_divu, tau_mu, tau_rotu, tau_sclr

    area = 4d0*MATH_PI*radius**2/(20d0*4**max_level) ! average area of a triangle
    dx_min = 0.891d0 * sqrt (4d0/sqrt(3d0) * area)   ! edge length of average triangle

    area = 4d0*MATH_PI*radius**2/(20d0*4**min_level)
    dx_max = sqrt (4d0/sqrt(3d0) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, 1.4d0*dx_min/Udim)
    dt_init = dt_cfl

    C = 5d-3 ! <= 0.5 if explicit for Laplacian diffusion, <= 0.1 for hyperdiffusion
    C_rotu = C / 4**Laplace_order_init  ! <= 1.09e-3 for hyperdiffusion (lower than exact limit 1/24^2 = 1.7e-3 due to non-uniform grid)
    C_divu = C
    C_mu   = 0d0
    C_b    = 0d0

    ! Diffusion time scales
    tau_mu   = dt_cfl / C_mu
    tau_b    = dt_cfl / C_b
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0d0
       visc_divu = 0d0
       visc_rotu = 0d0
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
       write (6,'(4(a,es8.2),/)') "C_mu = ", C_mu,  " C_b = ", C_b, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_turb (dom, i, j, zlev, offs, dims, 'bathymetry')
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
       call topo_turb (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0       
    end if
  end subroutine set_penal

  subroutine topo_turb (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer :: d, id, id_i

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth
    case ("penalize") ! analytic land mass with smoothing
       penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0   
    end select
  end subroutine topo_turb

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid

    if (time/DAY <= 2.0_8) then
       tau_zonal = tau_0 * sin (MATH_PI/4 * time/DAY)
    else
       tau_zonal = tau_0
    end if
    tau_merid = 0d0
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

  subroutine set_save_level_case
    ! Compute depth of saved layer
    implicit none
    real(8) :: eta, z_s
    real(8), dimension(0:zlevels) :: z

    eta = 0d0; z_s = max_depth
    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", z(save_zlev), " [m])"
  end subroutine set_save_level_case

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

    z_coords_case(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1.0_8 - z_0 / z_s
       z_coords_case(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords_case

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

  subroutine cal_r_max
    ! Calculates minimum relative mass and checks diffusion stability limits
#ifdef MPI
    use mpi
#endif
    implicit none
    integer :: ierror, k, l

    r_max_loc = 1d-16
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (cal_rmax_loc, l, k, 0, 0)
       end do
    end do

#ifdef MPI
    call MPI_Allreduce (r_max_loc, r_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
#else
    r_max = r_max_loc
#endif
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
            / porous_density (d, id+1, zlev)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)) &
            / porous_density (d, idE+1, zlev)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idNE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1)) &
            / porous_density (d, idNE+1, zlev)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idN+1)  + sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)) &
            / porous_density(d, idN+1, zlev)
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc

  real(8) function bottom_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Bottom boundary flux boundary condition for vertical diffusion of buoyancy (e.g. heat source)
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

  subroutine zonal_mean (avg, y2_avg)
    ! Computes zonal means of all prognostic variables and cubic spline interpolant
    ! (projects coarsest resolution onto lat-lon plane)
    implicit none
    real(8), dimension(Ny(1):Ny(2),1:zlevels,1:4) :: avg, y2_avg

    integer :: d, j, k, l

    l = min_level ! coarsest level

    do k = 1, zlevels
       ! Mass perturbation
       call project_field_onto_plane (sol(S_MASS,k), l, 0d0) 
       avg(:,k,1) = sum (field2d, 1) / size (field2d, 1)
       call spline (lat, avg(:,k,1), Nproj/2, 1d35, 1d35, y2_avg(:,k,1))

       ! Buoyancy
       do d = 1, size(grid)
          scalar =>    trend(S_TEMP,k)%data(d)%elts
          mass   =>      sol(S_MASS,k)%data(d)%elts
          temp   =>      sol(S_TEMP,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_buoyancy, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
          end do
          nullify (mass, mean_m, temp, mean_t, scalar)
       end do
       call project_field_onto_plane (trend(S_TEMP,k), l, 0d0)
       avg(:,k,2) = sum (field2d, 1) / size (field2d, 1)
       call spline (lat, avg(:,k,2), Nproj/2, 1d35, 1d35, y2_avg(:,k,2))

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo   => sol(S_VELO,k)%data(d)%elts
          velo1  => grid(d)%u_zonal%elts
          velo2  => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_edge_node, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       call project_array_onto_plane ("u_zonal", l, 0d0)
       avg(:,k,3) = sum (field2d, 1) / size (field2d, 1)
       call spline (lat, avg(:,k,3), Nproj/2, 1d35, 1d35, y2_avg(:,k,3))

       call project_array_onto_plane ("v_merid", l, 0d0)
       avg(:,k,4) = sum (field2d, 1) / size (field2d, 1)
       call spline (lat, avg(:,k,4), Nproj/2, 1d35, 1d35, y2_avg(:,k,4))
    end do
  end subroutine zonal_mean

  subroutine write_zonal_avg (var, y2, filename)
    ! Tests cubic spline interpolation of zonally averaged values
    implicit none
    real(8), dimension(Ny(1):Ny(2)) :: var, y2
    character(*)                    :: filename

    integer :: i, N_interp 
    real(8) :: lat_interp, var_interp

    open (unit=10, file='zonal_'//trim(filename), status='REPLACE')
    do i = Ny(1), Ny(2)
       write (10,'(2(es11.4,1x))') lat(i), var(i)
    end do
    close (10)

    open (unit=10, file='interp_'//trim(filename), status='REPLACE')
    N_interp = Nproj*5
    do i = -N_interp/2, N_interp/2
       lat_interp = -90d0 + lon_lat_range(2) / (N_interp + 1) * (i + N_interp/2)/MATH_PI * 180d0
       call splint (lat, var, y2, Nproj/2, lat_interp, var_interp)
       write (10,'(2(es11.4,1x))') lat_interp, var_interp
    end do
    close (10)
  end subroutine write_zonal_avg

  subroutine trend_nudge (q, dq)
    ! Trend relaxation to mean buoyancy
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq

    integer :: d, k, p

    do k = 1, zlevels
       do d = 1, size(grid)
          mass => q(S_MASS,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          dmass => dq(S_MASS,k)%data(d)%elts
          dtemp => dq(S_TEMP,k)%data(d)%elts
          dvelo => dq(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (nudging_scalars, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (nudging_velo,    grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dscalar, dvelo, mass, mean_m)
       end do
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_nudge

  subroutine nudging_scalars (dom, i, j, zlev, offs, dims)
    ! Relax buoyancy to mean
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: lat_i, lon_i, b_zonal, b_zonal_0

    id_i = idx (i, j, offs, dims) + 1

    call cart2sph (dom%node%elts(id_i), lon_i, lat_i)

    call splint (lat, zonal_0(:,zlev,2), y2_0(:,zlev,2), Nproj/2, lat_i/DEG, b_zonal_0)
    call splint (lat, zonal  (:,zlev,2), y2  (:,zlev,2), Nproj/2, lat_i/DEG, b_zonal  )

    dmass(id_i) = 0d0
    dtemp(id_i) = (mean_m(id_i) + mass(id_i)) * (b_zonal_0 - b_zonal) / tau_nudge
  end subroutine nudging_scalars

  subroutine nudging_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN, idNE

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    dvelo(EDGE*id+RT+1) = proj_nudge (vel_nudge, zlev, dom%node%elts(id+1),  dom%node%elts(idE+1))
    dvelo(EDGE*id+DG+1) = proj_nudge (vel_nudge, zlev, dom%node%elts(idN+1), dom%node%elts(id+1))
    dvelo(EDGE*id+UP+1) = proj_nudge (vel_nudge, zlev, dom%node%elts(id+1),  dom%node%elts(idNE+1))
  end subroutine nudging_velo

  subroutine vel_nudge (lon_i, lat_i, k, u, v)
    ! Zonal and meridional components of velocity trend for nudging 
    implicit none
    integer :: k
    real(8) :: lon_i, lat_i, u, v

    real(8) :: uzonal, uzonal_0, vmerid, vmerid_0

    call splint (lat, zonal_0(:,k,3), y2_0(:,k,3), Nproj/2, lat_i/DEG, uzonal_0) 
    call splint (lat, zonal  (:,k,3), y2  (:,k,3), Nproj/2, lat_i/DEG, uzonal  )

    call splint (lat, zonal_0(:,k,4), y2_0(:,k,4), Nproj/2, lat_i/DEG, vmerid_0) 
    call splint (lat, zonal  (:,k,4), y2  (:,k,4), Nproj/2, lat_i/DEG, vmerid  )

    u = (uzonal_0 - uzonal) / tau_nudge
    v = (vmerid_0 - vmerid) / tau_nudge
  end subroutine vel_nudge

  real(8) function proj_nudge (vel_fun, k, ep1, ep2)
    ! Finds velocity in direction from points ep1 to ep2 at mid-point of this vector
    ! given a function for zonal u and meridional v velocities as a function of longitude and latitude
    implicit none
    integer     :: k
    type(Coord) :: ep1, ep2
    external    :: vel_fun

    type(Coord) :: co, e_zonal, e_merid, vel
    real(8)     :: lon, lat, u_zonal, v_merid

    co = mid_pt (ep1, ep2)

    ! Find longitude and latitude coordinates of point co
    call cart2sph (co, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),             0.0_8) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    ! Function returning zonal and meridional velocities given longitude and latitude
    call vel_fun (lon, lat, k, u_zonal, v_merid)

    ! Velocity vector in Cartesian coordinates
    vel = vec_plus (vec_scale(u_zonal,e_zonal), vec_scale(v_merid,e_merid))

    ! Project velocity vector on direction given by points ep1, ep2
    proj_nudge = inner (direction(ep1, ep2), vel)
  end function proj_nudge
end module test_case_mod
