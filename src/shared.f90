module shared_mod
  use param_mod
  implicit none

  type Coord
     real(8) :: x, y, z
  end type Coord

  type Areas
     real(8), dimension(6) :: part
     real(8)               :: hex_inv
  end type Areas

  type(Coord), parameter :: ORIGIN = Coord (0d0, 0d0, 0d0)

  ! Numbers of triangles and edges per grid element
  integer, parameter :: TRIAG = 2, EDGE = 3

  ! Indices for edges
  integer, parameter :: RT = 0, DG = 1, UP = 2

  ! Indices for triangles
  integer, parameter :: LORT = 0, UPLT = 1

  ! Indices for nodes and edges of type Float_Field variables (e.g. sol, wave_coeff)
  integer, parameter :: AT_NODE = 1, AT_EDGE = 2

  ! Shifts on regular (i,j) grid
  integer, parameter :: JPLUS = 1
  integer, parameter :: IPLUS = 2
  integer, parameter :: JMINUS = 3
  integer, parameter :: IMINUS = 4
  integer, parameter :: IJPLUS = 5
  integer, parameter :: IPLUSJMINUS = 6
  integer, parameter :: IJMINUS = 7
  integer, parameter :: IMINUSJPLUS = 8

  ! Neighbouring patch indices for use in index arrays offs and dims 
  integer, parameter :: NORTH     = 1
  integer, parameter :: EAST      = 2
  integer, parameter :: SOUTH     = 3
  integer, parameter :: WEST      = 4
  integer, parameter :: NORTHEAST = 5
  integer, parameter :: SOUTHEAST = 6
  integer, parameter :: SOUTHWEST = 7
  integer, parameter :: NORTHWEST = 8

  ! Number of children nodes associated to each parent node
  integer, parameter :: N_CHDRN = 4 

  ! Domain parameters
  integer, parameter :: N_BDRY            = 8                            ! number of boundary patches associated to each patch
  integer, parameter :: N_ICOSAH_LOZENGE  = 10                           ! number of lozenges (coarse regular domains) in icosahedron
  integer, parameter :: N_SUB_DOM_PER_DIM = 2**DOMAIN_LEVEL              ! number of subdomains per lozenge in each direction
  integer, parameter :: N_SUB_DOM         = N_SUB_DOM_PER_DIM**2         ! total number of sub-domains per lozenge
  integer, parameter :: N_GLO_DOMAIN      = N_ICOSAH_LOZENGE * N_SUB_DOM ! total number of domains at coarsest level (number of cores must be <= N_GLO_DOMAIN)
  integer, parameter :: PATCH_LEVEL       = MIN_LEVEL - DOMAIN_LEVEL - 1 ! patch level: MIN_LEVEL = DOMAIN_LEVEL + PATCH_LEVEL + 1
  integer, dimension(:), allocatable :: n_domain                         ! number of subdomains on each processor

  ! Thickness of boundary overlaps between lozenges (ghost points or halo)
  integer, parameter :: BDRY_THICKNESS = 2

  integer, parameter :: FROZEN = 32

  ! Label for active nodes
  integer, parameter :: TOLRNZ = 16

  ! Label for adjacent nodes  in position (space) only 
  integer, parameter :: ADJSPACE = 14

  ! Label for nodes whose flux can be obtained by restriction from fine level
  integer, parameter :: RESTRCT = 12

  ! Label for nodes where trend is uniformly accurate
  integer, parameter :: TRND = 10

  ! Label for adjacent zone nodes in either position (space) or scale
  integer, parameter :: ADJZONE = 8

  ! Label for nodes added for consistency between adaptive velocity (edge) and mass (hexagon) nodes
  integer, parameter :: CONSIST = 4

  ! Label nodes needed for trisk operators 
  integer, parameter :: TRSK = 2

  integer, parameter :: NODE = 3

  integer, parameter :: ZERO =  0 
  integer, parameter :: NONE = -1
  integer, parameter :: POLE = -2 ! label for two pole points 

  ! logical integer parameters
  integer, parameter :: FALSE = 0
  integer, parameter :: TRUE  = 1

  integer, parameter :: ON_LINE  = 2
  integer, parameter :: INSIDE   = 0
  integer, parameter :: OUTER1   = 1
  integer, parameter :: OUTER2   = 2
  integer, parameter :: COINSIDE = 3

  ! Nearest two neighbour flux/velocity interpolation points U, V, W (i.e. RT,UP,DG)
  ! 
  ! Z = zero shift
  ! P = positive shift
  ! M = negative shift

  ! Note that there are 16 flux locations but only 14 distinct weights
  !
  ! First neighbours Uij, Vij, Wij where i and j can be any of (M,Z,P)
  integer, parameter :: UZM = 0
  integer, parameter :: UPZ = 1
  integer, parameter :: UMZ = 2
  integer, parameter :: UZP = 3

  integer, parameter :: VMM = 4
  integer, parameter :: VPM = 5
  integer, parameter :: VMP = 6
  integer, parameter :: VPP = 7

  integer, parameter :: WMM = 8
  integer, parameter :: WPM = 9
  integer, parameter :: WMP = 10
  integer, parameter :: WPP = 11

  ! Indices used by two flux locations (C is centre: no shift in x direction)
  integer, parameter :: CMM = 12 ! same as WMMM
  integer, parameter :: CPP = 13 ! same as WPPP

  ! Second neighbours Wijj, Vijj
  integer, parameter :: WMMM = 12
  integer, parameter :: WPPP = 13

  integer, parameter :: VMPP = 14
  integer, parameter :: VPMM = 15

  ! First diagonal neighbours of hexagon points 
  integer, parameter :: MP = 16
  integer, parameter :: PP = 17
  integer, parameter :: PM = 18
  integer, parameter :: MM = 19

  ! Weights for various interpolation schemes
  integer, dimension(3)     :: hex_s_offs
  integer, dimension(2,10)  :: nghb_pt
  integer, dimension(3,10)  :: hex_sides, no_adj_tri
  integer, dimension(2,2,3) :: end_pt, opp_no
  integer, dimension(2,4,3) :: bfly_no, bfly_no2
  integer, dimension(3,2,3) :: adj_tri
  integer, dimension(3,4,3) :: bfly_tri

  ! Used in grid smoothing routine
  integer, dimension(2,3) :: O2 
  data O2 /2,3, 3,1, 1,2/ 

  ! Indices of prognostic variables in sol, trend etc
  integer, parameter    :: S_VELO = 1, S_MASS = 2, S_TEMP = 3
  integer, parameter    :: S_DIVU = 1, S_ROTU = 2
  integer               :: N_VECTOR, N_SCALAR, N_VARIABLE
  integer, dimension(2) :: scalars

  ! Number of each variable per grid element (at hexagon nodes, triangle nodes, or edges) 
  integer, dimension(:), allocatable :: MULT, POSIT

  ! Grid optimization choices
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, HR_GRID = 2

  ! Define land and sea regions
  real(8), parameter :: LAND = 1, SEA = 0

  ! Basic grid parameters
  integer, parameter :: z_null = -1 ! place holder argument for functions not currently using z levels
  integer :: max_level              ! maximum grid refinement levels in pseudo-horizontal directions
  integer :: level_fill             ! make all grid points active for scales l <= level_fill
  integer :: zlevels                ! number of levels in vertical direction
  integer :: zmax                   ! zmax=zlevels+1 for a separate free surface layer, zmax=zlevels otherwise
  integer :: save_levels            ! number of vertical levels to save
  integer :: level_start, level_end, level_save, optimize_grid
  
  integer, dimension(AT_NODE:AT_EDGE) :: n_active ! number of active points at grid locations (node and edge)
  
  real(8) :: tol ! relative tolerance for all variables

  ! Basic constants (uses MKS system of units)

  ! Math
  real(8), parameter :: MATH_PI = acos (-1d0)

  ! Length
  real(8), parameter :: METRE   = 1
  real(8), parameter :: KM      = 1000 * METRE

  ! Mass
  real(8), parameter :: KG      = 1
  real(8), parameter :: GRAM    = KG / 1000

  ! Time
  real(8), parameter :: SECOND  = 1d0
  real(8), parameter :: MINUTE  = 60d0  * SECOND
  real(8), parameter :: HOUR    = 60d0  * MINUTE
  real(8), parameter :: DAY     = 24d0  * HOUR
  real(8), parameter :: WEEK    =   7d0 * DAY
  real(8), parameter :: YEAR    = 365d0 * DAY
   
  ! Angle
  real(8), parameter :: RAD     = 1d0
  real(8), parameter :: DEG     = MATH_PI / 180d0

  ! Force
  real(8), parameter :: NEWTON  = KG * METRE / SECOND**2

  ! Pressure
  real(8), parameter :: Pa      = NEWTON / METRE**2
  real(8), parameter :: hPa     =  100d0 * Pa
  real(8), parameter :: kPa     = 1000d0 * Pa

  ! Heat and energy
  real(8), parameter :: KELVIN  = 1d0
  real(8), parameter :: CELSIUS = KELVIN
  real(8), parameter :: JOULE   = KG * METRE**2 / SECOND**2
  real(8), parameter :: WATT    = JOULE / SECOND
  
  ! Simulation variables
  integer                                       :: coarse_iter, cp_idx, err_restart, fine_iter
  integer                                       :: iadapt, ibin, irebalance, iremap, istep, istep_cumul, iwrite
  integer                                       :: n_diffuse, nbins, nstep_init, save_zlev
  integer                                       :: resume, Laplace_order, Laplace_order_init
  integer(8)                                    :: itime
  integer, parameter                            :: nvar_zonal = 9   ! number of zonal statistics to calculate
  integer, dimension(:), allocatable            :: n_node_old, n_patch_old
  integer, dimension(:,:), allocatable          :: Nstats, Nstats_glo

  real(8)                                       :: alpha, a_0, b_0, lambda_1, lambda_2, mu_1, mu_2, nu_0, T_ref, S_ref
  real(8)                                       :: dbin, dt, dt_init, dt_write, dx_min, dx_max, time_end, time
  real(8)                                       :: omega, radius, grav_accel, cfl_adv, cfl_bar, cfl_num, kmax, Q_sr, ref_density
  real(8)                                       :: initotalmass, mass_error, max_depth, min_depth, min_mass, totalmass
  real(8)                                       :: e_min, Kt_const, Kt_0, Kv_0, Kv_bottom, rb_0
  real(8)                                       :: theta1, theta2, coarse_tol, fine_tol, visc_divu, visc_rotu
  real(8)                                       :: c1, c_p, c_s, c_v, gamma, H_rho, kappa, p_0, p_top, R_d, wave_speed
  real(8)                                       :: hex_int
  real(8), dimension(:),         allocatable    :: bounds, C_visc, pressure_save, visc_sclr
  real(8), dimension(:),         allocatable    :: a_vert, b_vert, a_vert_mass, b_vert_mass
  real(8), dimension(:,:),       allocatable    :: lnorm, threshold, threshold_def
  real(8), dimension(:,:,:),     allocatable    :: zonal_avg, zonal_avg_glo
  real(8), dimension(3)                         :: L_diffusion
  real(8), dimension (10*2**(2*DOMAIN_LEVEL),3) :: nonunique_pent_locs
  real(8), dimension (12,3)                     :: unique_pent_locs

  character(255)                                :: run_id, test_case, remapscalar_type, remapvelo_type, timeint_type
  
  logical :: adapt_dt, compressible, default_thresholds, eos_nl, fill, implicit_diff_sclr, implicit_diff_divu
  logical :: log_iter, log_mass, match_time, mode_split, penalize, rebalance, remap, uniform, vert_diffuse
  logical :: sigma_z, tke_closure
contains
  subroutine init_shared_mod
    logical :: initialized = .false.

    if (initialized) return ! Initialize only once
    initialized = .true.

    ! Initialize variable indices and arrays
    N_VECTOR = 1
    N_SCALAR = 2
    N_VARIABLE = N_VECTOR + N_SCALAR
    scalars = (/ N_VECTOR+1, N_VARIABLE /)
    allocate (MULT(1:N_VARIABLE), POSIT(1:N_VARIABLE))
    allocate (visc_sclr(scalars(1):scalars(2)))
    allocate (C_visc(1:N_VARIABLE))
       
    ! Specify the multiplicity per grid element of each quantity
    MULT(S_VELO) = EDGE
    MULT(scalars(1):scalars(2)) = 1

    ! Specify the position on the grid of each quantity
    POSIT(S_VELO) = AT_EDGE
    POSIT(scalars(1):scalars(2)) = AT_NODE

    ! i                        nghb_pt
    ! 1    1    1    0   -1   -1    0    1    1    0   -1
    ! 2    0    1    1    0   -1   -1    0    1    1    0
    nghb_pt  = reshape ((/ 1, 0, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 0, 1, 1, 0, 1, -1, 0 /), (/2, 10/))
    
    ! i  j       end_pt
    ! 1  1    0    1    0
    ! 1  2    1    0    0
    ! 2  1    0    1    0
    ! 2  2    0    0    1
    end_pt = reshape ((/ 0,  0, 1, 0, 1, 1, 0, 0, 0, 0,  0, 1 /), (/2, 2, 3/))

    ! i  j       opp_no
    ! 1  1    0    0    1
    ! 1  2    1    1   -1
    ! 2  1   -1    1    1
    ! 2  2    1    0    0
    opp_no = reshape ((/ 0, -1, 1, 1, 0, 1, 1, 0, 1, 1, -1, 0 /), (/2, 2, 3/))

    ! i  j      adj_tri
    ! 1  1    0    0    0
    ! 1  2    0    0   -1
    ! 2  1   -1    0    0
    ! 2  2    0    0    0
    ! 3  1    1    1    1
    ! 3  2    0    0    0
    adj_tri  = reshape ((/ 0, -1, UPLT, 0, 0, LORT, 0, 0, UPLT, 0, 0, LORT, 0, 0, UPLT, -1, 0, LORT /), (/3, 2, 3/))
    
    ! i                     no_adj_tri
    ! 1     0    0   -1   -1   -1    0    0    0   -1   -1
    ! 2     0   -1   -1   -1    0    0    0   -1   -1   -1
    ! 3     0    1    0    1    0    1    0    1    0    1
    no_adj_tri = reshape ((/0, 0, LORT, 0, -1, UPLT, -1, -1, LORT, -1, -1, UPLT, -1, 0, LORT, 0, 0, UPLT, 0, 0, LORT, &
         0, -1, UPLT, -1, -1, LORT, -1, -1, UPLT/), (/3, 10/))

    ! i                     hex_sides
    ! 1     0    0    0   -1   -1    0    0    0    0   -1
    ! 2     0    0    0    0   -1   -1    0    0    0    0
    ! 3     0    1    2    0    1    2    0    1    2    0
    hex_sides = reshape ((/ 0, 0, RT, 0, 0, DG, 0, 0, UP, -1, 0, RT, -1, -1, DG, 0, -1, UP, 0, 0, RT, 0, 0, DG, 0, 0, UP, &
         -1, 0, RT /), (/3, 10/))

    hex_s_offs = (/ 2, 0, 4 /)

    bfly_tri = reshape ((/ -1, -1, LORT, 0, -1, LORT, 1, 0, UPLT, 0, 0, UPLT, 0, 1, LORT, -1, 0, LORT, 0, -1, UPLT, 1, 0, &
         UPLT, 0, 0, LORT, 0, 1, LORT, -1, 0, UPLT, -1, -1, UPLT /), (/3, 4, 3/))
    
    bfly_no  = reshape ((/ -1, -1, 1, -1, 2, 1,  0, 1, 1, 2, -1,  0,  0, -1, 2, 1, 1, 0, 1, 2,  -1, 1, -1, -1 /), (/2, 4, 3/))
    bfly_no2 = reshape ((/ -3, -2, 1, -2, 3, 2, -1, 2, 1, 3, -3, -1, -1, -3, 3, 1, 2, -1, 2, 3, -2, 1, -2, -3 /), (/2, 4, 3/))

    ! Initialize values
    ! (these parameters may be reset in the test case file, but are needed for compilation)
    resume              = NONE
    cp_idx              = NONE
    err_restart         = 0
    istep               = 0
    istep_cumul         = 0
    iwrite              = 0
    time                = 0d0
    max_level           = MIN_LEVEL
    level_start         = MIN_LEVEL
    level_end           = level_start
    level_fill          = MIN_LEVEL
    nstep_init          = -1                                ! nstep_init gradually increasing small time steps after restart
    
    ! Default logical switches, most are reset in the input file
    adapt_dt            = .true.                            ! dynamically adapt time step (T) or use time step based on initial conditions (F) 
    compressible        = .true.                            ! compressible equations (T) or Boussinesq incompressible (F)
    default_thresholds  = .true.                            ! use default thresholds (T) or calculate dynamically (F)
    fill                = .false.                           ! fill up grid to level j_fill if true (T)
    log_iter            = .false.                           ! print residual error in elliptic solver (T)
    log_mass            = .true.                            ! compute minimum mass and mass conservation, relatively expensive (T)
    match_time          = .false.                           ! match time exactly for data saving (T)
    mode_split          = .false.                           ! calculate barotropic free surface mode separately (T)
    penalize            = .false.                           ! include penalization of topography (T)
    rebalance           = .true.                            ! rebalance computational load at each checkpoint if T
    sigma_z             = .false.                           ! use Schepetkin/CROCO type sigma-z vertical coordinates (T) or A/B hybrid coordinates (F)
    remap               = .true.                            ! remap Lagrangian coordinates (T) or no remapping (F)
    tke_closure         = .false.                           ! use TKE closure for eddy viscosity (T) or analytic form (F)
    uniform             = .true.                            ! uniform vertical grid in pressure (T) or hybrid (F)
    vert_diffuse        = .false.                           ! include vertical diffusion in ocean models (T)

    ! Default numerical method values
    alpha               = 1d0                               ! porosity
    cfl_adv             = 1.4d0                             ! advective CFL number in mode split case
    cfl_bar             = 1d0                               ! baroclinic CFL number in mode split case
    cfl_num             = 1d0                               ! CFL number (barotropic CFL in mode split case)
    C_visc              = 5d-4                              ! dimensionless diffusion coefficients
    iadapt              = 1                                 ! adapt horizontal grid every iadapt time step
    irebalance          = 5                                 ! interval for checking rebalance (only active if using AMPI)
    iremap              = 10                                ! remap every iremap time steps
    level_save          = level_start                       ! level to save
    Laplace_order_init  = 0                                 ! 0 = no diffusion, 1 = Laplacian diffusion, 2 = second-order iterated Laplacian hyperdiffusion
    n_diffuse           = 1                                 ! include diffusion every n_diffuse steps
    optimize_grid       = HR_GRID                           ! type of optimization of coarse grid
    remapscalar_type    = "PPR"                             ! remapping scheme for scalars
    remapvelo_type      = "PPR"                             ! remapping scheme for velocity
    save_levels         = 1                                 ! vertical level to save
    timeint_type        = "RK45"                            ! time integration scheme (RK3 is default for incompressible case)
    coarse_iter         = 30                                ! maximum number of coarse scale bicgstab iterations for elliptic solver
    fine_iter           = 200                               ! maximum number of fine scale jacobi iterations for elliptic solver
    tol                 = 5d-3                              ! relative tolerance for adaptivity
    coarse_tol          = 1d-9                              ! tolerance for coarse scale bicgstab elliptic solver
    fine_tol            = 1d-3                              ! tolerance for fine scale jacobi iterations
    zlevels             = 20                                ! number of vertical levels
    
    ! Default physical parameters
    ! (these parameters are typically reset in test case file, but are needed for compilation)
    c_p                 = 1004.64d0   * JOULE / (KG*KELVIN)   ! specific heat at constant pressure for air (= 3991.87 for seawater)
    c_v                 = 717.6d0     * JOULE / (KG*KELVIN)   ! specfic heat at constant volume c_v = R_d - c_p
    grav_accel          = 9.80616d0   * METRE / SECOND**2     ! gravitational acceleration
    p_top               = 0d0         * hPa                   ! pressure at upper interface of top vertical layer (should be non-zero for Lin remapping)
    R_d                 = 287d0       * JOULE / (KG*KELVIN)   ! ideal gas constant for dry air in joules per kilogram Kelvin
    ref_density         = 1000d0      * KG                    ! reference density for incompressible case
    omega               = 7.292d-05   * RAD / SECOND          ! rotation rate of Earth
    radius              = 6371.22d0   * KM                    ! radius of Earth
    p_0                 = 1000d0      * hPA                   ! standard pressure
    visc_sclr           = 0d0         * METRE**2 / SECOND     ! kinematic viscosity of scalars 
    visc_divu           = 0d0         * METRE**2 / SECOND     ! kinematic viscosity of divergence of velocity 
    visc_rotu           = 0d0         * METRE**2 / SECOND     ! kinematic viscosity of vorticity 

    kappa               = R_d/c_p                             ! heat capacity ratio

    ! Parameters for ocean (incompressible) model
    c1                  = 1d-16     * METRE / SECOND              ! value for internal wave speed (used for incompressible cases)
    c_s                 = 1500d0    * METRE / SECOND              ! sound speed for seawater
    e_min               = 0d0       * METRE**2 / SECOND**2        ! minimum TKE for vertical diffusion
    Kt_0                = 1.2d-5    * METRE**2 / SECOND           ! NEMO value for minimum/initial eddy diffusion
    Kv_0                = 1.2d-4    * METRE**2 / SECOND           ! NEMO value for minimum/initial eddy viscosity
    Kt_const            = 1d-6      * METRE**2 / SECOND           ! analytic value for eddy diffusion (tke_closure = .false.)
    Kv_bottom           = 2d-3      * METRE**2 / SECOND           ! analytic value for eddy viscosity (tke_closure = .false.)
    max_depth           = 4d0       * KM                          ! maximum depth 
    min_depth           = max_depth                               ! minimum depth 
    Q_sr                = 0d0       * WATT / METRE**2             ! penetrative part of solar short wave radiation
    rb_0                = 4d-4      * METRE /SECOND               ! NEMO value for bottom friction
    H_rho               = c_s**2 / grav_accel                     ! density scale height

    ! Equation of state parameters for ocean model
    eos_nl              = .false.                                 ! nonlinear equation of state if true
    a_0                 = 1.6550d-1 * KG / METRE**3 / CELSIUS     ! linear coefficient of thermal expansion for seawater 
    b_0                 = 7.6554d-1 * KG / METRE**3 / (GRAM / KG) ! linear haline expansion coefficient for seawater
    lambda_1            = 5.9520d-2                               ! cabbeling coefficient in T^2
    lambda_2            = 5.4914d-4                               ! cabbeling coefficient in S^2
    mu_1                = 1.4970d-4 / METRE                       ! thermobaric coefficient in temperature (pressure effect)
    mu_2                = 1.1090d-5 / METRE                       ! thermobaric coefficient in salinity (pressure effect)
    nu_0                = 2.4341d-3                               ! cabbeling coefficient in temperature S, T
    T_ref               = 10d0      * CELSIUS                     ! reference temperature
    S_ref               = 35d0      * GRAM / KG                   ! reference salinity

    ! Theta parameters for barotropic-baroclinic mode splitting
    ! theta1 > 0.75 and theta2 > 0.75 stable for all wavenumbers (otherwise unstable over a small interval of small wavenumbers)
    theta1              = 0.8d0                                   ! external pressure gradient in barotropic-baroclinic splitting (1 = fully implicit, 0.5 = Crank-Nicolson)
    theta2              = 0.8d0                                   ! barotropic flow divergence in barotropic-baroclinic splitting (1 = fully implicit, 0.5 = Crank-Nicolson)
  end subroutine init_shared_mod

  real(8) function eps ()
    eps = radius * 1d-13
  end function eps

  integer function max_nodes_per_level (lev, entity)
    integer           :: lev
    integer, optional :: entity
    
    max_nodes_per_level = 10*4**lev
    if (present (entity)) then 
       max_nodes_per_level = max_nodes_per_level*entity
    else ! mass node
       max_nodes_per_level = max_nodes_per_level+2
    end if
  end function max_nodes_per_level

  real(8) function exp__flush(x)
    real(8) :: x
    
    if (x > -1.0d2) then
       exp__flush = exp (x)
    else
       exp__flush = 0d0
    end if
  end function exp__flush
end module shared_mod
