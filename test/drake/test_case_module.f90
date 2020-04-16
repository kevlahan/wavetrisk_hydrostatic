Module test_case_mod
  ! Module file for Drake passage test case
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
  integer                              :: npts_penal, n_smooth
  real(8)                              :: beta, bv, c1, delta_I, delta_M, delta_S, delta_sm, drho, f0, L_R, Rey, Ro
  real(8)                              :: friction_coeff, max_depth, min_depth, scale, tau, top_layer, u_wbc
  real(8)                              :: resolution
  logical                              :: drag, mean_split
contains
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
          sol(S_TEMP,zlev)%data(d)%elts(id_i) = sol(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy (x_i, zlev)
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
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) * buoyancy (x_i, zlev)
       else
          sol_mean(S_MASS,zlev)%data(d)%elts(id_i) = 0.0_8
          sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) = 0.0_8
       end if
    end if
    sol_mean(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
  end subroutine init_mean
  
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

  real(8) function buoyancy (x_i, zlev)
    ! Buoyancy profile
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    integer     :: zlev
    type(Coord) :: x_i

    if (zlevels /= 2) then
       buoyancy = 0.0_8
    else
       if (zlev == 2) then ! less dense fluid in upper layer
          buoyancy = -drho/ref_density
       else
          buoyancy = 0.0_8
       end if
    end if
  end function buoyancy

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
       threshold_new = lnorm
       ! Correct for zero velocity case
       do k = 1, zmax
          if (threshold_new(S_MASS,k) == 0.0_8) threshold_new(S_MASS,k) = 1d16
          if (threshold_new(S_TEMP,k) == 0.0_8) threshold_new(S_TEMP,k) = 1d16
          if (threshold_new(S_VELO,k) == 0.0_8) threshold_new(S_VELO,k) = 1d16 
       end do
       threshold_new = tol * threshold_new
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
    real(8) :: dz, eta_surf

    allocate (threshold(1:N_VARIABLE,1:zmax));     threshold     = 0.0_8
    allocate (threshold_def(1:N_VARIABLE,1:zmax)); threshold_def = 0.0_8

    eta_surf = 0.0_8
    do k = 1, zlevels
       dz = a_vert_mass(k) * eta_surf + b_vert_mass(k) * max_depth
       lnorm(S_MASS,k) = ref_density*dz
       lnorm(S_TEMP,k) = ref_density*dz
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used
    if (adapt_trend) lnorm = lnorm/Tdim

    threshold_def = tol * lnorm
  end subroutine initialize_thresholds

  subroutine initialize_dt_viscosity 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, C_visc, tau_divu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle

    ! CFL limit for time step
    if (mode_split) then ! geostrophic time scale
       dt_cfl = cfl_num*dx_min/Udim 
    else ! external wave time scale
       dt_cfl = cfl_num*dx_min/(wave_speed + Udim) 
    end if
    dt_init = dt_cfl

    ! Permeability penalization parameter
    eta = dt_cfl
      
    ! Diffusion constants
    !C_visc = max (dt_cfl * beta * dx_min * resolution**3, 2.5d-4)  ! to ensure that Munk layer is resolved with resolution grid points
    C_visc = dt_cfl * beta * dx_min * resolution**3 ! to ensure that Munk layer is resolved with "resolution" grid points
    C_rotu = C_visc
    C_divu = C_visc 
    C_sclr = C_visc 

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
    if (.not. mean_split) visc_sclr = 0.0_8 ! cannot diffuse on scalars if mean_split = .false.
   
    if (rank == 0) then
       write (6,'(/,3(a,es8.2),a,/)') "dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
          " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
          " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
       write (6,'(a,es10.4,a)') "eta = ", eta," [s]"
    end if
  end subroutine initialize_dt_viscosity

  subroutine set_bathymetry (dom, i, j, z_null, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call analytic_topography (dom, i, j, z_null, offs, dims, "bathymetry")
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
       call analytic_topography (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
    end if
  end subroutine set_penal

  subroutine analytic_topography (dom, i, j, zlev, offs, dims, itype)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer     :: d, id, id_i, idW, idSW, idS
    real(8)     :: mask, strip, width, xmin, zmin
    type(Coord) :: p

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    p = dom%node%elts(id_i)

    select case (itype)
    case ("bathymetry")
       dom%topo%elts(id_i) = max_depth + surf_geopot (p) / grav_accel
    case ("penalize") ! Smoothed strip land mass
       if (p%x > dx_min) then
          xmin = 2*n_smooth*dx_min
          zmin = -radius * sin (35 * DEG)
          strip = 2*MATH_PI*radius / 20 ! width of land strip
          width = n_smooth * dx_min

          mask = (tanh ((p%y+strip/2)/width) - tanh ((p%y-strip/2)/width)) * &
               (1.0_8+tanh ((p%z-zmin)/width)) * (1.0_8+tanh ((p%x-xmin)/width)) / 8

          penal_node(zlev)%data(d)%elts(id_i) = mask
          
          ! Set edges to same mask value as associated node
          idW  = idx (i-1, j,   offs, dims)
          idSW = idx (i-1, j-1, offs, dims)
          idS  = idx (i,   j-1, offs, dims)

          penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1)   = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1))
          penal_edge(zlev)%data(d)%elts(EDGE*id+DG+1)   = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*id+DG+1))
          penal_edge(zlev)%data(d)%elts(EDGE*id+UP+1)   = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*id+UP+1))
          penal_edge(zlev)%data(d)%elts(EDGE*idW +RT+1) = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*idW+RT+1))
          penal_edge(zlev)%data(d)%elts(EDGE*idSW+DG+1) = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*idSW+DG+1))
          penal_edge(zlev)%data(d)%elts(EDGE*idS +UP+1) = max (mask, penal_edge(zlev)%data(d)%elts(EDGE*idS+UP+1))
       end if
    end select
  end subroutine analytic_topography

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    ! Idealized zonally and temporally averaged zonal and meridional wind stresses
    ! (based on Figure 4 from Ferreira et al J Climate 24, 992-1012 (2011) and Figure 4 from Gille J Atmos Ocean Tech 22, 1353-1372 respectively)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid
    real(8) :: peak

    peak = (abs(lat)*180/MATH_PI - 35.0_8) / 20
    tau_zonal = -0.15 * exp (-peak**2) * sin (abs(lat)*6) - 5d-3

    !peak = lat*180/MATH_PI / 15
    !tau_merid = -0.25 * exp (-peak**2) * sin (2*lat) * peak**2
    tau_merid = 0.0_8
  end subroutine wind_stress

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
    read (fid,*) varname, mode_split
    read (fid,*) varname, implicit_fs
    read (fid,*) varname, penalize
    read (fid,*) varname, max_level
    read (fid,*) varname, level_fill
    read (fid,*) varname, zlevels
    read (fid,*) varname, remap
    read (fid,*) varname, iremap
    read (fid,*) varname, adapt_trend
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, tol
    read (fid,*) varname, cfl_num
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, timeint_type
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    read (fid,*) varname, alpha
    read (fid,*) varname, drag
    close(fid)

    press_save = 0.0_8
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    call set_save_level
    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    delta_M = (visc_rotu/beta)**(1.0_8/3.0_8) * METRE ! Munk layer

    if (drag) then
       !friction_coeff =  1d-3                                 ! quadratic bottom friction coefficient (nemo)
       friction_coeff =  beta*abs(max_depth) * delta_M/4
       tau = abs(max_depth) / friction_coeff
       !friction_coeff = abs(max_depth)/(110*DAY) * METRE/SECOND ! linear bottom friction coefficient (nemo is 4e-4 for depth = 4000m)
    else
       friction_coeff = 0.0_8
    end if

    delta_S = (friction_coeff/abs(max_depth))/beta * METRE ! Stommel layer (want delta_S = delta_M/4)
    Rey     = u_wbc * delta_I / visc_rotu                  ! Reynolds number of western boundary current
    Ro      = u_wbc / (delta_M*f0)                         ! Rossby number (based on boundary current)
    
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
       write (6,'(A,I1)')     "implicit_fs          = ", implicit_fs
       write (6,'(A,L1)')     "penalize             = ", penalize
       write (6,'(A,i3)')     "n_smooth             = ", n_smooth
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
       write (6,'(A,i1)')     "n_diffuse            = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write [d]         = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY             = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance            = ", rebalance
       write (6,'(A,es10.4)') "time_end [d]         = ", time_end/DAY
       write (6,'(A,i6)')     "resume               = ", resume_init
       write (6,'(A,L1)')     "bottom drag          = ", drag
        
       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius             [km]     = ", radius / KM
       write (6,'(A,es10.4)') "omega              [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density        [kg/m^3] = ", ref_density
       write (6,'(A,es10.4)') "grav accel         [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "min_depth          [m]      = ", abs (min_depth)
       write (6,'(A,es11.4)') "max_depth          [m]      = ", abs (max_depth)
       write (6,'(A,es11.4)') "top_layer          [m]      = ", abs (top_layer)
       write (6,'(A,es11.4)') "density difference [kg/m^3] = ", drho
       write (6,'(A,es11.4)') "Brunt-Vaisala freq [1/s]    = ", bv
       write (6,'(A,es11.4)') "c0 wave speed      [m/s]    = ", wave_speed
       write (6,'(A,es11.4)') "c1 wave speed      [m/s]    = ", c1
       write (6,'(A,es11.4)') "eta (permeability) [s]      = ", eta
       write (6,'(A,es11.4)') "alpha (porosity)            = ", alpha
       write (6,'(A,es11.4)') "bottom friction    [m/s]    = ", friction_coeff
       write (6,'(A,es11.4)') "friction decay     [d]      = ", tau / DAY
       write (6,'(A,es11.4)') "f0 at 30 deg       [rad/s]  = ", f0
       write (6,'(A,es11.4,/)') "beta at 30 deg     [rad/ms] = ", beta
       write (6,'(A,es11.4)') "dx_min             [km]     = ", dx_min   / KM
       write (6,'(A,es11.4)') "L_R at 30 deg      [km]     = ", L_R      / KM
       write (6,'(A,es11.4)') "Inertial layer     [km]     = ", delta_I  / KM
       write (6,'(A,es11.4)') "Munk layer         [km]     = ", delta_M  / KM
       write (6,'(A,es11.4)') "Stommel layer      [km]     = ", delta_S  / KM
       write (6,'(A,es11.4,/)') "submesoscale       [km]     = ", delta_sm / KM
       write (6,'(A,es11.4)') "Rossby number               = ", Ro
       write (6,'(A,es11.4)') "Resolution of Munk layer    = ", resolution
       write (6,'(A,es11.4)') "Reynolds number             = ", Rey 
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

            write (12,'(5(es15.9,1x),i2,1x,i12,1x,4(es15.9,1x))')  time/HOUR, dt, &
                 threshold(S_MASS,zlevels), threshold(S_TEMP,zlevels), threshold(S_VELO,zlevels), &
                 level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine apply_initial_conditions
    use wavelet_mod
    implicit none
    integer :: k, l
    
    do l = level_start, level_end
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zlevels
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

  subroutine update
    ! Update means, bathymetry and penalization mask
    use wavelet_mod
    implicit none
    integer :: d, k, l, p
    
    do d = 1, size(grid)
!!$       do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches (causes memory leak?)
       do p = 3, grid(d)%patch%length ! update all patches
          call apply_onescale_to_patch (set_bathymetry, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          do k = 1, zmax
             call apply_onescale_to_patch (set_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do
    call barrier
    
    do k = 1, zmax
       do d = 1, size(grid)
!!$       do p = n_patch_old(d)+1, grid(d)%patch%length ! only update new patches (causes memory leak?)
          do p = 3, grid(d)%patch%length ! update all patches
             call apply_onescale_to_patch (init_mean, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end do
  end subroutine update

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
  
  subroutine set_save_level
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_zlev = zlevels
    save_height = 0.0_8

    if (rank==0) write (6,'(/,A,i2,A,es10.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level

  subroutine initialize_a_b_vert
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
    
    if (zlevels == 2) then ! special two layer case
       a_vert(0) = 0.0_8; a_vert(1) = 0.0_8;                    a_vert(2) = 1.0_8
       b_vert(0) = 1.0_8; b_vert(1) = top_layer/max_depth; b_vert(2) = 0.0_8
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
