Module test_case_mod
  ! Module file for upwelling test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use equation_of_state_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, npts_penal, resume_init
  real(8)                              :: dt_cfl, tau_diffusion, total_cpu_time
  real(8)                              :: dPdim, R_ddim, specvoldim, dTempdim

  ! Local variables
  real(8)                              :: beta, beta0, drho, f0, Rd, ref_temp
  real(8)                              :: r_max, r_max_loc
  real(8)                              :: n_smth_N, n_smth_S, width_N, width_S
  real(8)                              :: u_wbc
  real(8)                              :: lat_c, lat_width, tau_0, Tcline, width
  real(8), target                      :: bottom_friction_case
  real(4), allocatable, dimension(:,:) :: topo_data
  character(255)                       :: coords

  real(8), parameter :: slope = 1.75d-4 ! slope parameter (larger value -> steeper slope)
  real(8), parameter :: shift = 8
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
       physics_scalar_flux_case = 0d0
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

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    ! Additional physics for the source term of the velocity trend
    ! wind stress and bottom friction are included as surface fluxes in the split eddy viscosity split step
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id
    real(8), dimension(1:EDGE) :: diffusion, penal

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    ! Horizontal diffusion
    if (implicit_diff_divu) then
       diffusion = - visc_rotu * curl_rotu()
    else
       diffusion = (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())
    end if
    
    ! Penalization 
    penal = - penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)/dt * velo(EDGE*id+RT+1:EDGE*id+UP+1)

    physics_velo_source_case = diffusion + penal
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
       write (6,'(A,L1)')     "mode_split                     = ", mode_split
       write (6,'(A,L1)')     "penalize                       = ", penalize
       write (6,'(A,i3)')     "npts_penal                     = ", npts_penal
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
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/bottom_friction_case / DAY
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

  subroutine apply_initial_conditions_case
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

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initial perturbation to mean for an entire vertical column 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    eta = init_free_surface (dom%node%elts(id_i))
    z_s = topography%data(d)%elts(id_i)

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
          sol(S_MASS,zlevels)%data(d)%elts(id_i) = rho * eta
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = 0d0
       end if
       sol(S_TEMP,k)%data(d)%elts(id_i)                      = rho * dz(k) * buoyancy_init (z_k)
       sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do

    if (mode_split) then
       phi = phi_node (d, id_i, zlevels)
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i) = phi * eta ! free surface perturbation
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0d0
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

    eta = 0d0
    z_s = topography%data(d)%elts(id_i)

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

    b_vert(0) = 1d0 ; b_vert(zlevels) = 0d0
    if (trim (coords) == "uniform") then
       do k = 1, zlevels-1
          b_vert(k) = 1d0 - dble(k)/dble(zlevels)
       end do
    elseif (trim(coords) == "croco") then
       p = (/ -5.7831d0,  18.9754d0, -24.6521d0,  16.1698d0, -5.7092d0, 0.9972d0 /)
       do k = 1, zlevels-1
          z = dble(k)/dble(zlevels)
          b_vert(k) = p(1)*z**5 + p(2)*z**4 + p(3)*z**3 + p(4)*z**2 + p(5)*z + p(6)
       end do
    end if
    a_vert = 1d0 - b_vert

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

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential: postive if greater than mean seafloor                                                                                        
    ! Gives minimum depth of approximately 24 m                                                                                                          
    implicit none
    integer       :: d, id
    
    real(8)     :: amp, b_max, lat, lon, y, w_N, w_S, ns_S, ns_N
    type(Coord) :: p

    p = grid(d)%node%elts(id)
    
    call cart2sph (p, lon, lat)

    b_max = abs (max_depth - min_depth)
    lat = lat / DEG

    if (abs(lat-lat_c) <= lat_width/2d0) then
       amp = b_max / (1d0 - tanh (-slope*width/shift))
       y = (lat - (lat_c - lat_width/2d0))/180d0 * MATH_PI*radius ! y = 0 at low latitude boundary of channel                                               
       surf_geopot_case = amp * (1d0 - tanh (slope * (f(y) - width/shift)))
    else
       surf_geopot_case = b_max
    end if
    Surf_geopot_case = grav_accel * surf_geopot_case
  end function surf_geopot_case

  real(8) function f (y)
    implicit none
    real(8) :: y

    if (y <= width/2d0) then
       f = y
    else
       f = width - y
    end if
  end function f

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0d0
  end function init_free_surface

  real(8) function buoyancy_init (z)
    ! Initial buoyancy at depth z
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8) :: z

    real(8) :: rho

    rho = density_eos (S_ref, temp_init (z), z)
    buoyancy_init = (ref_density - rho) / ref_density
  end function buoyancy_init

  real(8) function temp_init (z)
    implicit none
    real(8) :: z

    real(8)            :: hz, strat, z0, z1
    real(8), parameter :: h0 = 150d0, hz_0 = 6.5d0, T0 = 14d0, z0_0 = -35d0, z1_0 = -75d0

    strat = abs(max_depth)
    hz = hz_0 * abs(max_depth/h0)
    z0 = z0_0 * abs(max_depth/h0)
    z1 = z1_0 * abs(max_depth/h0)

    temp_init = T_ref + 4d0*tanh ((z - z0) / hz) + (z - z1) / strat
  end function temp_init

  subroutine print_density
    implicit none
    integer                       :: k
    real(8)                       :: bv, c_k, c1, drho, dz_l, eta, rho, rho_above, z_k, z_s, z_above
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

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
            k, z_k, dz(k), ref_density * (1d0 - buoyancy_init (z_k))
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

       rho_above = ref_density * (1d0 - buoyancy_init (z_above))
       rho  = ref_density * (1d0 - buoyancy_init (z_k))
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
       call cal_lnorm (order)
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
       lnorm(S_TEMP,k) = drho * dz(k)
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C, C_b, C_divu, C_mu, C_rotu, dlat, tau_b, tau_divu, tau_mu, tau_rotu, tau_sclr

    area = 4d0*MATH_PI*radius**2/(20d0*4**max_level) ! average area of a triangle
    dx_min = 0.891d0 * sqrt (4d0/sqrt(3d0) * area) ! edge length of average triangle

    area = 4d0*MATH_PI*radius**2/(20d0*4**min_level)
    dx_max = sqrt (4d0/sqrt(3d0) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, 1.4d0*dx_min/u_wbc, 1.2d0*dx_min/c1)
    dt_init = dt_cfl

    C = 5d-3 ! <= 1/2 if explicit
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
       write (6,'(4(a,es8.2),/)') "C_mu = ", C_mu,  " C_b = ", C_mu, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if

    ! Penalization parameterss
    dlat = 0.5d0*dble(npts_penal+2) * (dx_max/radius) / DEG ! widen channel to account for boundary smoothing

    width_S = 90d0 + (lat_c - (lat_width/2 + dlat))
    width_N = 90d0 - (lat_c + (lat_width/2 + dlat))

    ! Smoothing exponent for land mass
    n_smth_S = 4d0*radius * width_S*DEG / (dx_max * dble(npts_penal+2))
    n_smth_N = 4d0*radius * width_N*DEG / (dx_max * dble(npts_penal+2))
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_upwelling (dom, i, j, zlev, offs, dims, 'bathymetry')
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
       call topo_upwelling (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0       
    end if
  end subroutine set_penal

  subroutine topo_upwelling (dom, i, j, zlev, offs, dims, itype)
    ! Returns penalization mask for land penal and bathymetry coordinate topo 
    ! uses radial basis function for smoothing (if specified)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer                        :: d, e, id, id_e, id_i, l, n_coarse, nsmth
    real(8)                        :: dx
    type(Coord)                    :: p
    type(Coord), dimension(1:EDGE) :: q

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1

    p = dom%node%elts(id_i)
    l = dom%level%elts(id_i)

    select case (itype)
    case ("bathymetry")
       n_coarse = 8 ! ensure r_max is small enough
       nsmth = n_coarse * 2**(l - min_level) 
       dx    = dx_max / 2**(l - min_level) 
       topography%data(d)%elts(id_i) = max_depth + smooth (surf_geopot_case, d, id_i, dx, nsmth) / grav_accel
    case ("penalize") ! analytic land mass with smoothing
       nsmth = npts_penal * 2**(l - min_level) 
       dx    = dx_max / 2**(l - min_level) 
       penal_node(zlev)%data(d)%elts(id_i) = smooth (mask, d, id_i, dx, nsmth)

       q(RT+1) = dom%node%elts(idx(i+1, j,   offs, dims)+1)
       q(DG+1) = dom%node%elts(idx(i+1, j+1, offs, dims)+1) 
       q(UP+1) = dom%node%elts(idx(i,   j+1, offs, dims)+1)
       do e = 1, EDGE
          id_e = EDGE*id + e
          penal_edge(zlev)%data(d)%elts(id_e) = smooth (mask, d, id_i, dx, nsmth)
       end do
    end select
  end subroutine topo_upwelling

  real(8) function smooth (fun, d, id, dx, npts)
    ! Smooth a function using radial basis functions
    implicit none
    integer      :: d, id, npts
    real(8)      :: dx

    integer     :: ii, jj
    real(8)     :: dtheta, lat, lat0, lon, lon0, nrm, r, rbf, wgt
    type(Coord) :: p, q

    interface
       real(8) function fun (d, id)
         use domain_mod
         integer :: d, id
       end function fun
    end interface

    p = grid(d)%node%elts(id)
    
    if (npts == 0) then
       smooth = fun (d, id)
    else
       dtheta = dx/radius
       call cart2sph (p, lon0, lat0)
       nrm = 0d0
       rbf = 0d0
       do ii = -npts, npts
          lat = lat0 + dtheta * ii
          do jj = -npts, npts
             lon = lon0 + dtheta * jj

             q = sph2cart(lon, lat)
             r = norm (vector(p, q))
             wgt = radial_basis_fun ()

             nrm = nrm + wgt
             rbf = rbf + wgt * fun (d, id)
          end do
       end do
       smooth = rbf /nrm
    end if
  contains
    real(8) function radial_basis_fun ()
      ! Radial basis function for smoothing topography
      implicit none

      radial_basis_fun = exp (-(r/(dble(npts)*dx/2d0))**2)
    end function radial_basis_fun
  end function smooth

  real(8) function mask (d, id)
    implicit none
    integer      :: d, id
    
    type(Coord) :: p

    real(8) :: lat, lon

    p = grid(d)%node%elts(id)

    call cart2sph (p, lon, lat)

    mask = exp__flush (- abs((lat/DEG+90d0)/width_S)**n_smth_S) + exp__flush (- abs((lat/DEG-90d0)/width_N)**n_smth_N)
  end function mask

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid

    if (time/DAY <= 2d0) then
       tau_zonal = tau_0 * sin (MATH_PI/4d0 * time/DAY)
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
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
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

    integer :: d, id, idE, idN, idNE
    real(8) :: dz0, dz_e, r_loc

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    d    = dom%id + 1

    if (dom%mask_n%elts(id+1) >= ADJZONE .and. penal_node(zlev)%data(d)%elts(id+1) < 1d-3) then
       dz0  = (sol(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(id+1)) &
            / porous_density (d, id+1, zlev)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)) &
            / porous_density (d, idE+1, zlev)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idNE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1)) &
            / porous_density (d, idNE+1, zlev)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)
       
       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idN+1)  + sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)) &
            / porous_density(d, idN+1, zlev)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc

  real(8) function bottom_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Bottom boundary condition for vertical diffusion of buoyancy (e.g. heat source)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    bottom_buoy_flux_case= 0d0
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
end module test_case_mod

