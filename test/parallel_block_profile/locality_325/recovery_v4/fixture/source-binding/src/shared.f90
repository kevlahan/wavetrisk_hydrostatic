module shared_mod
  
  use kind_mod,  only : sp, dp
  use param_mod, only : DOMAIN_LEVEL, MIN_LEVEL
  
  implicit none
  
  type Coord
     real(dp) :: x, y, z
  end type Coord

  type Coord_r4
     real(sp):: x, y, z
  end type Coord_r4

  type Areas
     real(dp) :: part(6)
     real(dp) :: hex_inv
  end type Areas
  
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Parameters
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  type(coord), parameter :: penta_node_std(12) = [ &   ! coordinates of pentagons (unrotated grid)
       Coord ( 0.0,  0.0,  1.0), &
       Coord ( 0.0,  0.0, -1.0), &

       Coord ( 7.23606797749978936e-01,  5.25731112119133592e-01,  4.47213595499957983e-01), &
       Coord (-2.76393202250021008e-01,  8.50650808352039989e-01,  4.47213595499957983e-01), &
       Coord (-8.94427190999915855e-01,  0.0,                      4.47213595499957983e-01), &
       Coord (-2.76393202250021119e-01, -8.50650808352039878e-01,  4.47213595499957983e-01), &       
       Coord ( 7.23606797749978825e-01, -5.25731112119133814e-01,  4.47213595499957983e-01), &
       
       Coord ( 8.94427190999915855e-01,  0.0,                     -4.47213595499957983e-01), &
       Coord ( 2.76393202250021119e-01,  8.50650808352039878e-01, -4.47213595499957983e-01), &
       Coord (-7.23606797749978825e-01,  5.25731112119133703e-01, -4.47213595499957983e-01), &
       Coord (-7.23606797749979047e-01, -5.25731112119133481e-01, -4.47213595499957983e-01), &
       Coord ( 2.76393202250020786e-01, -8.50650808352039989e-01, -4.47213595499957983e-01) ]

  ! Domain parameters
  integer, parameter :: N_BDRY            =  8                             ! number of boundary patches associated to each patch
  integer, parameter :: POLE              = -2                             ! label for two pole points
  integer, parameter :: N_ICOSAH_LOZENGE  = 10                             ! number of lozenges (coarse regular domains) in icosahedron
  integer, parameter :: BDRY_THICKNESS    =  3                             ! width of halo/ghost cells boundary (DO NOT MODIFY)
  integer, parameter :: N_CHDRN           =  4                             ! number of children nodes associated to each parent node

  integer, parameter :: N_SUB_DOM_PER_DIM = 2**DOMAIN_LEVEL                ! number of subdomains per lozenge in each direction
  integer, parameter :: N_SUB_DOM         = N_SUB_DOM_PER_DIM**2           ! total number of sub-domains per lozenge
  integer, parameter :: N_GLO_DOMAIN      = N_ICOSAH_LOZENGE * N_SUB_DOM   ! total number of domains at coarsest level 
  integer, parameter :: PATCH_LEVEL       = MIN_LEVEL - (DOMAIN_LEVEL + 1) ! patch level: MIN_LEVEL = DOMAIN_LEVEL+1 + PATCH_LEVEL
  
  ! Shifts on regular (i,j) grid
  integer, parameter :: JPLUS             = 1
  integer, parameter :: IPLUS             = 2
  integer, parameter :: JMINUS            = 3
  integer, parameter :: IMINUS            = 4
  integer, parameter :: IJPLUS            = 5
  integer, parameter :: IPLUSJMINUS       = 6
  integer, parameter :: IJMINUS           = 7
  integer, parameter :: IMINUSJPLUS       = 8

  ! Neighbouring patch indices for use in index arrays offs and dims 
  integer, parameter :: NORTH             = 1
  integer, parameter :: EAST              = 2
  integer, parameter :: SOUTH             = 3
  integer, parameter :: WEST              = 4
  integer, parameter :: NORTHEAST         = 5
  integer, parameter :: SOUTHEAST         = 6
  integer, parameter :: SOUTHWEST         = 7
  integer, parameter :: NORTHWEST         = 8 

  ! Mask values
  integer, parameter :: FROZEN            = 32  ! nodes that should not be modified
  integer, parameter :: TOLRNZ            = 16  ! active nodes
  integer, parameter :: ADJSPACE          = 14  ! adjacent nodes in position (space) only 
  integer, parameter :: RESTRCT           = 12  ! nodes whose flux can be obtained by restriction from fine level
  integer, parameter :: ADJZONE           = 8   ! adjacent zone nodes in either position (space) or scale
  integer, parameter :: TRSK              = 2   ! nodes needed for TRISK operators
  
  integer, parameter :: INSIDE            = 0
  integer, parameter :: OUTER1            = 1
  integer, parameter :: OUTER2            = 2

  ! Logical integer parameters
  integer, parameter :: FALSE             = 0
  integer, parameter :: TRUE              = 1

  integer, parameter :: ZERO              =  0 
  integer, parameter :: NONE              = -1

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

   ! Numbers of triangles and edges per grid element
  integer, parameter :: TRIAG = 2, EDGE = 3

  ! Indices for edges
  integer, parameter :: RT = 0, DG = 1, UP = 2

  ! Indices for triangles
  integer, parameter :: LORT = 0, UPLT = 1

  ! Label for node
  integer, parameter :: NODE = 3

  ! Indices for nodes and edges of type Float_Field variables (e.g. sol, wave_coeff)
  integer, parameter :: AT_NODE = 1, AT_EDGE = 2

  ! Indices for longitude and latitude components
  integer, parameter :: LON_x = 1, LAT_y = 2

  ! Indices of prognostic variables in sol, trend etc
  integer, parameter :: S_VELO = 1, S_MASS = 2, S_TEMP = 3, S_DIVU = 4, S_ROTU = 5
  integer, parameter :: N_VECTOR = 1  ! number of vector variables
  integer, parameter :: N_SCALAR = 2  ! number of scalar variables
  integer, parameter :: N_VARIABLE = N_VECTOR + N_SCALAR

  ! Indices of scalar variables
  integer, parameter :: scalars(2) = [ N_VECTOR+1, N_VARIABLE ]

  ! Number of each variable per grid element (at hexagon nodes, triangle nodes, or edges)
  integer, parameter :: MULT(1:N_VARIABLE) = [EDGE, 1, 1]

  ! Position on the grid of each quantity
  integer, parameter :: POSIT(1:N_VARIABLE) = [AT_EDGE, AT_NODE, AT_NODE]

  ! Indices of sub scale orography model (SSO)
  integer, parameter :: S_MU = 1, S_THETA = 2, S_GAMMA = 3, S_SIGMA = 4

  ! Grid optimization choices
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, DATA_GRID = 2

  ! Define land and sea regions
  real(dp), parameter :: LAND = 1, SEA = 0

  ! Basic grid parameters
  integer, parameter :: z_null = -1 ! place holder argument for functions not currently using z levels
  
  type(Coord), parameter :: ORIGIN = Coord (0.0_dp, 0.0_dp, 0.0_dp)

  ! Weights for various interpolation schemes

  ! i  j       end_pt
  ! 1  1    0    1    0
  ! 1  2    1    0    0
  ! 2  1    0    1    0
  ! 2  2    0    0    1
  integer, parameter :: end_pt(2,2,3) = reshape( &
       [0,0, 1,0, 1,1, 0,0, 0,0, 0,1], [2,2,3])

  ! i  j       opp_no
  ! 1  1    0    0    1
  ! 1  2    1    1   -1
  ! 2  1   -1    1    1
  ! 2  2    1    0    0
  integer, parameter :: opp_no(2,2,3) = reshape( &
       [0,-1, 1,1, 0,1, 1,0, 1,1, -1,0], [2,2,3])

  ! i  j      adj_tri
  ! 1  1    0    0    0
  ! 1  2    0    0   -1
  ! 2  1   -1    0    0
  ! 2  2    0    0    0
  ! 3  1    1    1    1
  ! 3  2    0    0    0
  integer, parameter :: adj_tri(3,2,3) = reshape( &
       [0,-1,UPLT, 0,0,LORT, 0,0,UPLT, 0,0,LORT, 0,0,UPLT, -1,0,LORT], &
       [3,2,3])

  ! i                     no_adj_tri
  ! 1     0    0   -1   -1   -1    0    0    0   -1   -1
  ! 2     0   -1   -1   -1    0    0    0   -1   -1   -1
  ! 3     0    1    0    1    0    1    0    1    0    1
  integer, parameter :: no_adj_tri(3,10) = reshape( &
       [0,0,LORT, 0,-1,UPLT, -1,-1,LORT, -1,-1,UPLT, -1,0,LORT, &
       0,0,UPLT, 0,0,LORT, 0,-1,UPLT, -1,-1,LORT, -1,-1,UPLT], &
       [3,10])

  ! i                     hex_sides
  ! 1     0    0    0   -1   -1    0    0    0    0   -1
  ! 2     0    0    0    0   -1   -1    0    0    0    0
  ! 3     0    1    2    0    1    2    0    1    2    0
  integer, parameter :: hex_sides(3,10) = reshape( &
       [0,0,RT, 0,0,DG, 0,0,UP, -1,0,RT, -1,-1,DG, 0,-1,UP, &
       0,0,RT, 0,0,DG, 0,0,UP, -1,0,RT], &
       [3,10])

  integer, parameter :: hex_s_offs(3) = [2,0,4]

  integer, parameter :: bfly_tri(3,4,3) = reshape( &
       [-1,-1,LORT, 0,-1,LORT, 1,0,UPLT, 0,0,UPLT, &
       0, 1,LORT, -1,0,LORT, 0,-1,UPLT, 1,0,UPLT, &
       0, 0,LORT, 0,1,LORT, -1,0,UPLT, -1,-1,UPLT], &
       [3,4,3])

  integer, parameter :: bfly_no(2,4,3) = reshape( &
       [-1,-1, 1,-1, 2,1, 0,1, 1,2,-1,0, &
       0,-1, 2,1, 1,0, 1,2,-1,1,-1,-1], &
       [2,4,3])

  integer, parameter :: bfly_no2(2,4,3) = reshape( &
       [-3,-2, 1,-2, 3,2, -1,2, 1,3,-3,-1, &
       -1,-3, 3,1, 2,-1, 2,3,-2,1,-2,-3], &
       [2,4,3])

  ! nghb_pt are offsets of neighbours of a node (i,j)
  !                 E    NE   N    W    SW   S    E    NE   N    S
  !   offset in i   1    1    0   -1   -1    0    1    1    0   -1
  !   offset in j   0    1    1    0   -1   -1    0    1    1    0
  integer, parameter :: nghb_pt(2,10) = reshape( &
       [1,0, 1,1, 0,1, -1,0, -1,-1, 0,-1, 1,0, 1,1, 0,1, -1,0], &
       [2,10] )

  ! Used in grid smoothing routine
  integer,  parameter :: O2(2,3) = reshape([2,3, 3,1, 1,2], [2,3])

  ! Basic constants (uses MKS system of units)

  ! Math
  real(dp), parameter :: MATH_PI = acos (-1.0_dp)

  ! Length
  real(dp), parameter :: METRE   = 1.0_dp
  real(dp), parameter :: KM      = 1000 * METRE

  ! Mass
  real(dp), parameter :: KG      = 1.0_dp
  real(dp), parameter :: GRAM    = KG / 1000

  ! Time
  real(dp), parameter :: SECOND  = 1.0_dp
  real(dp), parameter :: MINUTE  = 60  * SECOND
  real(dp), parameter :: HOUR    = 60  * MINUTE
  real(dp), parameter :: DAY     = 24  * HOUR
  real(dp), parameter :: WEEK    =   7 * DAY
  real(dp), parameter :: YEAR    = 365 * DAY

  ! Angle
  real(dp), parameter :: RAD     = 1.0_dp
  real(dp), parameter :: DEG     = MATH_PI / 180

  ! Force
  real(dp), parameter :: NEWTON  = KG * METRE / SECOND**2

  ! Pressure
  real(dp), parameter :: Pa      = NEWTON / METRE**2
  real(dp), parameter :: hPa     =  100 * Pa
  real(dp), parameter :: kPa     = 1000 * Pa

  ! Heat and energy
  real(dp), parameter :: KELVIN  = 1.0_dp
  real(dp), parameter :: CELSIUS = KELVIN
  real(dp), parameter :: JOULE   = KG * METRE**2 / SECOND**2
  real(dp), parameter :: WATT    = JOULE / SECOND

  ! Line feed character
  character(1), parameter :: lf  = char(10)                            
  

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Initialized variables
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Names
  character(255)          :: run_id       = "test"                              ! name of run
  character(255)          :: test_case    = "climate"                           ! name of test case
  character(255)          :: topo_file    = " "                                 ! name of file containing multiscale NCAR topography data
  character(255)          :: physics_type = "Held_Suarez"                       ! physics model used if physics_model is T

  ! Grid variables
  character(255) :: grid_type                 = "HR95JT"                        ! type of coarse data grid
  integer        :: optimize_grid             = DATA_GRID                       ! type of optimization of coarse grid
  integer        :: n_active(AT_NODE:AT_EDGE) = 0                               ! number of active points at grid locations (node and edge)
  real(dp)       :: theta_grid                = 0.0_dp                          ! rotation of standard icosahedron around pole
  type(coord)    :: penta_node(12)            = penta_node_std                  ! coordinates of rotated pentagon

  ! Simulation variables
  
  ! Initialize values
  ! (these parameters may be reset in the test case file, but are needed for compilation)
  integer :: resume                  = NONE
  integer :: cp_every                = 1
  integer :: cp_idx                  = NONE
  integer :: cp_init                 = 0
  integer :: err_restart             = 0
  integer :: istep                   = 0
  integer :: istep_cumul             = 0
  integer :: iwrite                  = 0
  integer :: iwrite_init             = 0
  integer :: max_level               = MIN_LEVEL                                ! maximum grid refinement levels in pseudo-horizontal directions
  integer :: level_start             = MIN_LEVEL
  integer :: level_end               = MIN_LEVEL
  integer :: level_fill              = MIN_LEVEL                                ! make all grid points active for scales l <= level_fill
  integer :: topo_min_level          = MIN_LEVEL                                ! minimum level of NCAR topography data 
  integer :: topo_max_level          = MIN_LEVEL                                ! maximum level of NCAR topography data 
  
  ! Logical switches
  logical :: adapt_dt                = .false.                                  ! dynamically adapt time step (T) or use time step based on initial conditions (F) 
  logical :: compressible            = .true.                                   ! compressible equations (T) or Boussinesq incompressible (F)
  logical :: default_thresholds      = .true.                                   ! use default thresholds (T) or calculate dynamically (F)
  logical :: fill                    = .false.                                  ! fill up grid to level j_fill if true (T)
  logical :: log_iter                = .false.                                  ! print residual error in elliptic solver (T)
  logical :: log_min_mass            = .false.                                  ! compute minimum mass, relatively expensive (T)
  logical :: log_total_mass          = .false.                                  ! mass conservation, relatively expensive (F)
  logical :: match_time              = .false.                                  ! exactly match time interval dt_write when saving data (F)
  logical :: mode_split              = .false.                                  ! calculate barotropic free surface mode separately (T)
  logical :: NCAR_topo               = .false.                                  ! use NCAR topography (requires pre-computation of float field topography) (T)
  logical :: penalize                = .false.                                  ! include penalization of topography (T)
  logical :: physics_model           = .false.                                  ! use physics model sub-step for compressible cases (T)
  logical :: remap                   = .true.                                   ! remap Lagrangian coordinates (T) or no remapping (F)
  logical :: sigma_z                 = .false.                                  ! use Schepetkin/CROCO type sigma-z vertical coordinates (T) or A/B hybrid coordinates (F)
  logical :: split_mean_perturbation = .false.
  logical :: sso                     = .false.                                  ! SSO (subgrid scale orography) parameterization for atmosphere
  logical :: tke_closure             = .false.                                  ! use TKE closure for eddy viscosity (T) or analytic form (F)
  logical :: uniform                 = .true.                                   ! uniform vertical grid in pressure (T) or hybrid (F)
  logical :: vert_diffuse            = .false.                                  ! include vertical diffusion in ocean models (T)
  
  ! Numerical method
  integer(8) :: itime                = 0
  integer    :: nstep_init           = -1                                       ! nstep_init gradually increasing small time steps after restart
  integer    :: iadapt               = 1                                        ! adapt horizontal grid every iadapt time step
  integer    :: irebalance           = 5                                        ! interval for checking rebalance (only active if using AMPI)
  integer    :: iremap               = 1                                        ! remap counter
  integer    :: iremap_max           = 5                                        ! maximum remap interval (every iremap_max dt)
  real(dp)   :: time                 = 0   * DAY
  real(dp)   :: time_end             = 1   * DAY
  real(dp)   :: dt                   = 100 * SECOND
  real(dp)   :: dt_init              = 100 * SECOND
  real(dp)   :: cfl_safety           = 0.9_dp                                   ! safety factor for maximum stable advectice time step
  real(dp)   :: dt_phys              = 30 * MINUTE                              ! interval for physics split step
  real(dp)   :: dt_write             = 5  * DAY                                 ! interval for writing data
  real(dp)   :: min_mass_remap       = 0.9_dp                                   ! minimum relative layer mass compared to initial value at which to remap
  real(dp)   :: porosity             = 1e-2_dp                                  ! porosity for solid regions when using penalization

  ! Order of Laplacian diffusion  0 = no diffusion, 1 = Laplacian diffusion, 2 = second-order iterated Laplacian hyperdiffusion
  integer  :: Laplace_sclr            = 2                                       ! scalars
  integer  :: Laplace_divu            = 2                                       ! div u
  integer  :: Laplace_rotu            = 2                                       ! rot u 
  integer  :: n_diffuse               = 1                                       ! include diffusion every n_diffuse steps

  character(255) :: remap_type        = "PPR"                                   ! remapping scheme for scalars
  character(255) :: remapscalar_type  = "PPR"                                   ! remapping scheme for scalars
  character(255) :: remapvelo_type    = "PPR"                                   ! remapping scheme for velocity

  character(255) :: timeint_type      = "RK4"                                   ! time integration scheme
  real(dp)       :: r_adv             = 2 * sqrt (2.0_dp)                       ! advective stability constant
  real(dp)       :: r_dif             = 2.77_dp                                 ! diffusive stability constant
  
  character(3)   :: linear_solver     = "FMG"                                   ! type of linear solver (for barotropic-barocline mode splitting)      
  character(3)   :: vtk_grid          = "tri"                                   ! type of grid to save when exporting data
  
  real(dp) :: tol                     = 0.0_dp                                  ! relative tolerance for adaptivity (default is non-adaptive)
  integer  :: zlevels                 = 30                                      ! number of vertical layers
  integer  :: zmin                    = 1                                       ! index of lowest vertical level,  1 for atmosphere, 0 for simple phys surf temp or -Nsoil for soil mod
  integer  :: zmax                    = 30                                      ! zmax=zlevels+1 for a separate free surface layer, zmax=zlevels otherwise
  integer  :: Nsoil                   = 0                                       ! number of soil layers in vertical direction: k = [zmin 0], zmin <= 0

  ! Physical parameters
  real(dp) :: c_p                     = 1004.64      * JOULE / (KG*KELVIN)      ! specific heat at constant pressure for air (= 3991.87 for seawater)
  real(dp) :: c_v                     = 717.6        * JOULE / (KG*KELVIN)      ! specific heat at constant volume c_v = R_d - c_p
  real(dp) :: gamma                   = 1004.64_dp / 717.6_dp                   ! heat capacity ratio
  real(dp) :: grav_accel              = 9.80616      * METRE / SECOND**2        ! gravitational acceleration
  real(dp) :: p_top                   = 0            * hPa                      ! pressure at upper interface of top vertical layer (should be non-zero for Lin remapping)
  real(dp) :: R_d                     = 287          * JOULE / (KG*KELVIN)      ! ideal gas constant for dry air in joules per kilogram Kelvin
  real(dp) :: kappa                   = 287.0_dp/1004.64_dp                     ! heat capacity ratio
  real(dp) :: ref_density             = 0            * KG / METRE**3            ! set ref_density to correct default value below if not set in test case
  real(dp) :: ref_density_air         = 1.225        * KG / METRE**3            ! reference density (compressible case: atmosphere)
  real(dp) :: ref_density_water       = 1028.0       * KG / METRE**3            ! reference density (incompressible case: seawater)
  real(dp) :: omega                   = 7.292e-05    * RAD / SECOND             ! rotation rate of Earth
  real(dp) :: radius                  = 6371.22      * KM                       ! radius of Earth
  real(dp) :: p_0                     = 1000         * hPA                      ! standard pressure
  real(dp) :: visc_sclr(scalars(1):scalars(2)) = 0   * METRE**2 / SECOND        ! kinematic viscosity of scalars
  real(dp) :: visc_divu               = 0            * METRE**2 / SECOND        ! kinematic viscosity of divergence of velocity 
  real(dp) :: visc_rotu               = 0            * METRE**2 / SECOND        ! kinematic viscosity of vorticity 


  ! Parameters for atmosphere (compressible) model
  real(dp) :: c_g                     = 287          * METRE / SECOND           ! gravity wave speed for atmosphere ( 198 m/s for ocean)
  real(dp) :: c_s                     = 340          * METRE / SECOND           ! sound speed for atmosphere        (1500 m/s for ocean)
  real(dp) :: wave_speed              = 340          * METRE / SECOND           ! wave speed used for CFL number (use c_g for ocean)

  ! Parameters for ocean (incompressible) model
  real(dp) :: c1                      = 1e-16        * METRE / SECOND           ! value for internal wave speed (used for incompressible cases)
  real(dp) :: max_depth               = 4            * KM                       ! maximum depth 
  real(dp) :: min_depth               = 4            * KM                       ! minimum depth 
  real(dp) :: H_rho                   = 340.0/ 9.80616                          ! density scale height

  ! Equation of state parameters for ocean model
  logical  :: eos_nl                  = .false.                                 ! nonlinear equation of state if true
  real(dp) :: a_0                     = 1.6550e-1 * KG / METRE**3 / CELSIUS     ! linear coefficient of thermal expansion for seawater 
  real(dp) :: b_0                     = 7.6554e-1 * KG / METRE**3 / (GRAM / KG) ! linear haline expansion coefficient for seawater
  real(dp) :: lambda_1                = 5.9520e-2                               ! cabbeling coefficient in T^2
  real(dp) :: lambda_2                = 5.4914e-4                               ! cabbeling coefficient in S^2
  real(dp) :: mu_1                    = 1.4970e-4 / METRE                       ! thermobaric coefficient in temperature (pressure effect)
  real(dp) :: mu_2                    = 1.1090e-5 / METRE                       ! thermobaric coefficient in salinity (pressure effect)
  real(dp) :: nu_0                    = 2.4341e-3                               ! cabbeling coefficient in temperature S, T
  real(dp) :: T_ref                   = 10        * CELSIUS                     ! reference temperature
  real(dp) :: S_ref                   = 35        * GRAM / KG                   ! reference salinity

  ! Nondimensionalization scales
  real(dp) :: Hdim                    = 4         * KM  
  real(dp) :: Ldim                    = 6371.22   * KM 
  real(dp) :: Mudim                   = 1028 * (4e3 / 30)
  real(dp) :: Pdim                    = 1000      * hPA 
  real(dp) :: Tdim                    = 1         * DAY
  real(dp) :: Tempdim                 = 200       * KELVIN
  real(dp) :: Thetadim                = 1028 * (4e3 / 30 )* 200
  real(dp) :: Udim                    = 1.0       * METRE / SECOND

  ! Theta parameters for barotropic-baroclinic mode splitting: 1 = fully implicit, 0.5 = Crank-Nicolson
  ! (avoid theta = 0.5 due to free surface standing wave instability)
  ! NEMO uses theta = 0.55 to ensure stability while avoiding over-damping  
  real(dp) :: theta1                  = 0.55_dp                                 ! external pressure gradient
  real(dp) :: theta2                  = 0.55_dp                                 ! barotropic flow divergence

  ! Diagnostic values
  real(dp) :: mass_error              = 0.0_dp
  real(dp) :: min_mass                = 1.0e16_dp


  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Allocatable variables
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  integer,   allocatable  :: n_domain(:), n_node_old(:), n_patch_old(:)
  real(dp),  allocatable  :: Area_avg(:), dx_avg(:), pressure_save(:)
  real(dp),  allocatable  :: a_vert(:), b_vert(:), a_vert_mass(:), b_vert_mass(:)
  real(dp),  allocatable  :: C_visc(:,:), lnorm(:,:), threshold(:,:), threshold_def(:,:)

  
contains
  
  
  function eps (scale) result(val)
    ! Minimum absolute error with respect to scale

    implicit none

    real(dp) :: val

    real(dp) :: scale
    val = abs(scale) * epsilon (1.0_dp)
  end function eps


  function exp__flush (x) result(val)

    implicit none

    real(dp) :: x
    real(dp) :: val

    if (x > -1.0d2) then
       val = exp (x)
    else
       val = 0.0_dp
    end if
  end function exp__flush

  
end module shared_mod
