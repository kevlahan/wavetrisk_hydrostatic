Module test_case_mod
  ! Module file for seamount test case
  use domain_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  implicit none

  ! Standard variables
  integer                              :: bathy_per_deg, CP_EVERY, resume_init
  real(8)                              :: dt_cfl, k_T, tau_diffusion, total_cpu_time
  real(8)                              :: dPdim, R_ddim, specvoldim, dTempdim

  ! Local variables
  real(8)                              :: beta, bu, bv, drho, drho_dz, f0, Rb, Rd, Rey, Ro
  real(8)                              :: delta, h0, lat_c, lon_c, r_max, r_max_loc, width
  real(8)                              :: radius_earth, omega_earth, scale, tau_0, tke_sea, visc
  real(8),                      target :: bottom_friction_case 
  real(4), allocatable, dimension(:,:) :: topo_data

  character(255)                       :: coords, stratification                   
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

    integer                    :: d, id, id_i
    real(8), dimension(1:EDGE) :: diffusion

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    diffusion =  (-1)**(Laplace_order-1) * (visc_divu * grad_divu() - visc_rotu * curl_rotu())

    physics_velo_source_case = diffusion
  contains
    function grad_divu()
      implicit none
      real(8), dimension(3) :: grad_divu

      integer :: idE, idN, idNE

      idE  = idx (i+1, j,   offs, dims)
      idNE = idx (i+1, j+1, offs, dims)
      idN  = idx (i,   j+1, offs, dims)

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
    read (fid,*) varname, tol
    read (fid,*) varname, cfl_num
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    read (fid,*) varname, Bu
    read (fid,*) varname, stratification   
    read (fid,*) varname, coords
    close(fid)

    drho    = - ref_density * (Bu * f0 * width/wave_speed)**2            ! density difference based on Burger number Bu
    drho_dz = drho/max_depth                                             ! approximate density gradient
    bv      = sqrt (grav_accel * drho_dz/ref_density)                    ! Brunt-Vaisala frequency
    c1      = bv * sqrt (abs(max_depth)/grav_accel)/MATH_PI * wave_speed ! first baroclinic mode speed for linear stratification
    Rb      = bv * abs(max_depth) / (MATH_PI*f0)                         ! first baroclinic Rossby radius of deformation

    ! Vertical level to save
    save_zlev = zlevels 

    press_save = 0d0
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    Rey  = Udim * delta / visc_rotu ! Reynolds number 

    call set_save_level_case

    call cal_r_max

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case                      = ", trim (test_case)
       write (6,'(A,A)')      "run_id                         = ", trim (run_id)
       write (6,'(A,es10.4)') "scale                          = ", scale
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
       write (6,'(A,L1)')     "remap                          = ", remap
       write (6,'(a,a)')      "remapscalar_type               = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type                 = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap                         = ", iremap
       write (6,'(A,es8.2)')  "theta1                         = ", theta1
       write (6,'(A,es8.2)')  "theta2                         = ", theta2
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
       write (6,'(A,es10.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", abs (min_depth)
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", abs (max_depth)
       write (6,'(a,a)')      "vertical coordinates           = ", trim (coords)
       write (6,'(a,a)')      "stratification                 = ", trim (stratification)
       write (6,'(A,es11.4)') "density difference   [kg/m^3]  = ", drho
       write (6,'(A,es11.4)') "Brunt-Vaisala freq      [1/s]  = ", bv
       write (6,'(A,es11.4)') "Burger number                  = ", bu
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed           [m/s]  = ", c1
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "alpha (porosity)               = ", alpha
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/bottom_friction_case / DAY
       write (6,'(A,es11.4)') "seamount latitude       [deg]  = ", lat_c / DEG
       write (6,'(A,es11.4)') "f0 at seamount        [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at seamount   [rad/ms]  = ", beta
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(A,es11.4)') "barotropic Rossby radius [km]  = ", Rd / KM
       write (6,'(A,es11.4,/)') "baroclinic Rossby radius [km]  = ", Rb / KM
       write (6,'(A,es11.4)') "Rossby number                  = ", Ro
       write (6,'(A,es11.4)') "Re (delta_I u_wbc / nu)        = ", Rey
       write (6,'(a,es11.4)') "r_max                          = ", r_max
       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'

       call print_density_pert
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    use lnorms_mod
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, maxv, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

    call cal_lnorm ("inf")
    maxv = 100 * maxval (lnorm(S_VELO,:))

    if (rank == 0) then
       write (6,'(a,es12.6,6(a,es8.2),a,i2,a,i12,4(a,es9.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', threshold(S_MASS,zlevels), &
            ' temp tol = ', threshold(S_TEMP,zlevels), &
            ' velo tol = ', threshold(S_VELO,zlevels), &
            ' maxv = ', maxv, &
            ' tke = ', tke_sea, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(7(es15.9,1x),i2,1x,i12,1x,4(es15.9,1x))')  time/DAY, dt, &
            threshold(S_MASS,zlevels), threshold(S_TEMP,zlevels), threshold(S_VELO,zlevels), maxv, tke_sea, &
            level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

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

    integer     :: d, id, id_i
    real (8)    :: dz, eta_surf, phi, porous_density, z, z_s
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)
    eta_surf = init_free_surface (x_i)
    z_s = topography%data(d)%elts(id_i)

    if (zlev == zlevels+1) then ! 2D barotropic mode
       phi = 1d0 + (alpha - 1d0) * penal_node(zlevels)%data(d)%elts(id_i)

       sol(S_MASS,zlev)%data(d)%elts(id_i) = phi * eta_surf ! free surface perturbation
       sol(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
    else ! 3D layers
       dz = a_vert_mass(zlev) * eta_surf + b_vert_mass(zlev) * z_s

       porous_density = ref_density * (1d0 + (alpha - 1d0) * penal_node(zlev)%data(d)%elts(id_i))

       if (zlev == zlevels) then
          sol(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * eta_surf
       else
          sol(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
       end if
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
    real (8)    :: dz, eta_surf, porous_density, z_s
    type(Coord) :: x_i

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1
    x_i  = dom%node%elts(id_i)
    eta_surf = init_free_surface (x_i)
    z_s = topography%data(d)%elts(id_i)

    if (zlev == zlevels+1) then
       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0d0
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0d0
    else
       dz = a_vert_mass(zlev) * eta_surf + b_vert_mass(zlev) * z_s

       porous_density = ref_density * (1d0 + (alpha - 1d0) * penal_node(zlev)%data(d)%elts(id_i))

       sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = porous_density * dz
       sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy_init (z_s, x_i, zlev)
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

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential: postive if greater than mean seafloor
    implicit none
    integer :: d, id
    
    real(8)     :: lon, lat, rgrc
    type(Coord) :: x_i

    x_i = grid(d)%node%elts(id)
    
    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)

    rgrc = radius*acos(sin(lat_c)*sin(lat)+cos(lat_c)*cos(lat)*cos(lon-lon_c))

    surf_geopot_case = grav_accel*h0 * exp__flush (-(rgrc/width)**2)
  end function surf_geopot_case

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0d0
  end function init_free_surface

  real(8) function buoyancy_init (z_s, x_i, zlev)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    integer     :: id_i, zlev
    real(8)     :: z_s
    type(Coord) :: x_i

    real(8) :: rho, z1, z2
    real(8) :: eta_surf = 0d0

    z1 = a_vert(zlev-1) * eta_surf + b_vert(zlev-1) * z_s
    z2 = a_vert(zlev)   * eta_surf + b_vert(zlev)   * z_s

    rho = 0.5 * (density_init (x_i, z1) + density_init (x_i, z2))

    buoyancy_init = (ref_density - rho) / ref_density 
  end function buoyancy_init

  real(8) function density_init (x_i, z)
    implicit none
    real(8)     :: z
    type(Coord) :: x_i

    if (trim(stratification) == "linear") then 
       density_init = ref_density + drho * (max_depth-z)/max_depth
    elseif (trim(stratification) == "exponential") then
       density_init = ref_density + drho * exp__flush (z/delta)
    end if
  end function density_init

  subroutine print_density_pert
    implicit none
    integer     :: k
    real(8)     :: eta_surf, lat, lon, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord) :: p 

    eta_surf = 0d0

!!$    p = sph2cart (lon_c, lat_c) ! centre of seamount
!!$    p = sph2cart (lon_c + width/radius / DEG, lat_c) ! edge of seamount
    p = sph2cart (lon_c + width/radius / DEG, lat_c + width/radius / DEG) ! mid ocean

    p%x = radius * p%x ; p%y = radius * p%y ; p%z = radius * p%z  

    z_s = max_depth 

    if (sigma_z) then
       z = z_coords_case (eta_surf, z_s)
    else
       z = a_vert * eta_surf + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)


    write (6,'(a)') " Layer      z         dz         rho "
    do k = 1, zlevels
       z_k = interp (z(k-1), z(k))
       write (6, '(2x, i2, 4x, 2(es9.2, 1x), es12.5)') &
            k, z_k, dz(k), ref_density * (1d0 - buoyancy_init (z_s, p, k))
    end do

    write (6,'(/,a)') " Interface     z"
    do k = 0, zlevels
       write (6, '(3x, i3, 5x, es9.2)') k, z(k)
    end do

    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density_pert

  subroutine set_thresholds_case
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
       call cal_lnorm (order)
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
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer     :: k
    real(8)     :: dz, eta_surf, z
    type(Coord) :: x_i

    x_i = Coord (radius, 0d0, 0d0)
    eta_surf = 0d0
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta_surf + b_vert_mass(k) * max_depth
       z = 0.5 * ((a_vert(k)+a_vert(k-1)) * eta_surf + (b_vert(k)+b_vert(k-1)) * max_depth)

       lnorm(S_MASS,k) = ref_density*dz

!!$       lnorm(S_TEMP,k) = ref_density*dz * buoyancy_init (max_depth, x_i, k)
       lnorm(S_TEMP,k) = 1d16
       if (lnorm(S_TEMP,k) == 0d0) lnorm(S_TEMP,k) = 1d16

       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, tau_divu, tau_rotu, tau_sclr

    area = 4d0*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3d0) * area)         ! edge length of average triangle

    area = 4d0*MATH_PI*radius**2/(20*4**min_level)
    dx_max = sqrt (4/sqrt(3d0) * area)

    ! Initial CFL limit for time step
    dt_cfl = min (cfl_num*dx_min/wave_speed, dx_min/c1)
    dt_init = dt_cfl

    if (Laplace_order_init == 0) then
       visc_sclr = 0d0
       visc_divu = 0d0
       visc_rotu = 0d0
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr = visc
       visc_divu = visc
       visc_rotu = visc
    elseif (Laplace_order_init >= 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 supported)'
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
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id + 1
    id  = idx (i, j, offs, dims)  + 1

    topography%data(d)%elts(id) = surf_geopot_case (d, id) / grav_accel
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

    penal_node(zlev)%data(d)%elts(id_i)                      = 0d0
    penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0       
  end subroutine set_penal

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    ! Idealized zonally and temporally averaged zonal and meridional wind stresses
    ! (based on Figure 4 from Ferreira et al J Climate 24, 992-1012 (2011) and Figure 4 from Gille J Atmos Ocean Tech 22, 1353-1372 respectively)

    implicit none
    real(8) :: lat, lon, peak, tau_zonal, tau_merid

    logical, parameter :: merid_stress = .false.

    peak = (abs(lat)*180/MATH_PI - 35d0) / 20
    tau_zonal = -tau_0 * 1.2 * exp (-peak**2) * sin (abs(lat)*6) - 5d-3*exp(-(lat*180/MATH_PI/10)**2)

    if (merid_stress) then
       peak = lat*180/MATH_PI / 15
       tau_merid = -tau_0 * 2.5 * exp (-peak**2) * sin (2*lat) * peak**2
    else
       tau_merid = 0d0
    end if
  end subroutine wind_stress

  subroutine set_save_level_case
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level_case

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k, kp

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (trim (coords) == "uniform") then 
       do k = 0, zlevels
          b_vert(k) = 1d0 - dble(k)/dble(zlevels)
       end do
    elseif (trim (coords) == "chebyshev") then
       do k = 0, zlevels
          b_vert(k) = (1d0 + cos (dble(k)/dble(zlevels) * MATH_PI)) / 2d0
       end do
    elseif (trim (coords) == "chebyshev_half") then
       do k = 0, zlevels
          b_vert(k) = 1d0 - sin (dble(k)/dble(zlevels) * MATH_PI/2d0)
       end do
    end if
    a_vert = 1d0 - b_vert
    
    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

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
    real(8) :: r_loc

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    d    = dom%id + 1

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       r_loc = abs (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) - sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)) / &
            (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idE+1))
       r_max_loc = max (r_max_loc, r_loc)

       r_loc = abs (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) - sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1)) / &
            (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1))
       r_max_loc = max (r_max_loc, r_loc)

       r_loc = abs (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) - sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)) / &
            (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idN+1))
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc

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

!!$    hc = min (abs(min_depth), abs(Tcline)) ! if specify position of thermocline
    hc = abs(min_depth)

    cff1 = 1d0 / sinh (theta_s)
    cff2 = 0.5d0 / tanh (0.50 * theta_s)

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
end module test_case_mod
