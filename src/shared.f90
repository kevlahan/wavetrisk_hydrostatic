module shared_mod
  use param_mod
  implicit none
  logical :: DEB_ON

  ! numbers of triangles and edges per grid element
  integer, parameter :: TRIAG = 2, EDGE = 3

  ! indices for edges
  integer, parameter :: RT = 0, DG = 1, UP = 2

  ! indices for triangles
  integer, parameter :: LORT = 0, UPLT = 1

  ! find nodes and edges of type Float_Field variables (e.g. sol, wave_coeff)
  integer, parameter :: AT_NODE = 1, AT_EDGE = 2

  ! shifts on regular (i,j) grid
  integer, parameter :: JPLUS = 1
  integer, parameter :: IPLUS = 2
  integer, parameter :: JMINUS = 3
  integer, parameter :: IMINUS = 4
  integer, parameter :: IJPLUS = 5
  integer, parameter :: IPLUSJMINUS = 6
  integer, parameter :: IJMINUS = 7
  integer, parameter :: IMINUSJPLUS = 8

  ! neighbouring patch indices for use in index arrays offs and dims 
  integer, parameter :: NORTH     = 1
  integer, parameter :: EAST      = 2
  integer, parameter :: SOUTH     = 3
  integer, parameter :: WEST      = 4
  integer, parameter :: NORTHEAST = 5
  integer, parameter :: SOUTHEAST = 6
  integer, parameter :: SOUTHWEST = 7
  integer, parameter :: NORTHWEST = 8

  ! number of children nodes associated to each parent node
  integer, parameter :: N_CHDRN = 4 

  ! domain parameters
  integer, parameter :: N_ICOSAH_LOZENGE = 10               ! number of lozenges (coarse regular domains) in icosahedron
  integer, parameter :: N_SUB_DOM_PER_DIM = 2**DOMAIN_LEVEL ! number of subdomains per lozenge in each direction
  integer, parameter :: N_SUB_DOM = N_SUB_DOM_PER_DIM**2    ! total number of sub-domains per lozenge
  integer, parameter :: N_BDRY = 8                          ! number of boundary patches associated to each patch
  integer, dimension(:), allocatable :: n_domain            ! number of subdomains on each processor

  ! thickness of boundary overlaps between lozenges (ghost points or halo)
  integer, parameter :: BDRY_THICKNESS = 2

  integer, parameter :: FROZEN = 32

  ! label for active nodes
  integer, parameter :: TOLRNZ = 16

   ! label for adjacent nodes  in position (space) only 
  integer, parameter :: ADJSPACE = 14

  ! label for nodes whose flux can be obtained by restriction from fine level
  integer, parameter :: RESTRCT = 12

  ! label for adjacent zone nodes in either position (space) or scale
  integer, parameter :: ADJZONE = 8

  ! label for nodes added for consistency between adaptive velocity (edge) and mass (hexagon) nodes
  integer, parameter :: CONSIST = 4

  ! label nodes needed for trisk operators 
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

  ! note that there are 16 flux locations but only 14 distinct weights
  !
  ! first neighbours Uij, Vij, Wij where i and j can be any of (M,Z,P)
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

  ! indices used by two flux locations (C is centre: no shift in x direction)
  integer, parameter :: CMM = 12 ! same as WMMM
  integer, parameter :: CPP = 13 ! same as WPPP

  ! second neighbours Wijj, Vijj
  integer, parameter :: WMMM = 12
  integer, parameter :: WPPP = 13

  integer, parameter :: VMPP = 14
  integer, parameter :: VPMM = 15

  ! first diagonal neighbours of hexagon points 
  integer, parameter :: MP = 16
  integer, parameter :: PP = 17
  integer, parameter :: PM = 18
  integer, parameter :: MM = 19

  ! weights for various interpolation schemes
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
  integer, parameter :: S_MASS = 1, S_TEMP = 2, S_PENL = 3 ! scalar variable labels
  integer, parameter :: N_SCALARS = 3                      ! number of scalar variables (first N_SCALAR elements of sol are scalars)

  integer, parameter :: S_VELO = N_SCALARS+1               ! vector variable labels
  integer, parameter :: N_VECTORS = 1                      ! number of vector variables
  
  integer, parameter :: N_VARS = N_SCALARS + N_VECTORS     ! total number of variables
  
  integer, parameter :: S_DIVU = 1, S_ROTU = 2

  ! Number of each variable per grid element (at hexagon nodes, triangle nodes, or edges) 
  integer, dimension(1:N_VARS) :: MULT

  ! Location of each variable on the grid (at an edge or at a node)
  integer, dimension(1:N_VARS) :: POSIT

  ! Grid optimization choices
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, HR_GRID = 2

  ! Define land and sea regions
  real(8), parameter :: LAND = 1, SEA = 0

  ! Basic grid parameters
  integer, parameter :: z_null = -1 ! place holder argument for functions not currently using z levels
  integer :: min_level, max_level   ! minimum and maximum grid refinement levels in pseudo-horizontal directions
  integer :: zlevels                ! number of levels in vertical direction
  integer :: save_levels            ! number of vertical levels to save
  integer :: level_start, level_end, level_save, optimize_grid
  
  integer, dimension(AT_NODE:AT_EDGE) :: n_active ! number of active points at grid locations (node and edge)
  
  real(8) :: tol ! relative tolerance for all variables

  ! Basic constants (assumes basic unit of time is seconds)
  integer, parameter :: SECOND = 1
  integer, parameter :: MINUTE = 60*SECOND
  integer, parameter :: HOUR = 60*MINUTE
  integer, parameter :: DAY = 24*HOUR
  integer, parameter :: WEEK = 7*DAY
  real(8), parameter :: METRE = 1
  real(8), parameter :: KM = 1000*METRE
  real(8), parameter :: MATH_PI = acos(-1.0_8)
  
  ! Simulation variables
  integer                                       :: cp_idx, err_restart, ibin, iremap, istep, istep_cumul, iwrite, n_diffuse, nbins
  integer                                       :: resume, Laplace_order, Laplace_order_init
  integer(8)                                    :: itime
  integer, parameter                            :: nvar_zonal = 9   ! number of zonal statistics to calculate
  integer, dimension(:,:), allocatable          :: Nstats, Nstats_glo
  
  real(8)                                       :: C_visc, dbin, dt, dt_init, dt_write, dx_min, dx_max, time_end, time
  real(8)                                       :: omega, radius, grav_accel, cfl_num, kmax, ref_density
  real(8)                                       :: visc_divu, visc_rotu
  real(8)                                       :: alpha, eta
  real(8), dimension(S_MASS:S_TEMP)             :: visc_sclr
  real(8)                                       :: c_p, c_v, gamma, gk, kappa, mean_depth, p_0, p_top, R_d, wave_speed
  real(8)                                       :: hex_int
  real(8)                                       :: min_mass, min_allowed_mass
  real(8), dimension(:),         allocatable    :: pressure_save, bounds
  real(8), dimension(:),         allocatable    :: a_vert, b_vert, a_vert_mass, b_vert_mass
  real(8), dimension(:,:),       allocatable    :: lnorm, threshold
  real(8), dimension(:,:,:),     allocatable    :: zonal_avg, zonal_avg_glo
  real(8), dimension(3)                         :: L_diffusion
  real(8), dimension (10*2**(2*DOMAIN_LEVEL),3) :: nonunique_pent_locs
  real(8), dimension (12,3)                     :: unique_pent_locs

  character(255)                                :: run_id, test_case, remapscalar_type, remapvelo_type, timeint_type
  
  logical :: adapt_dt, adapt_trend, compressible, default_thresholds, penalize, perfect, rebalance, remap, uniform
contains
  subroutine init_shared_mod
    integer :: v
    logical :: initialized = .false.

    if (initialized) return ! Initialize only once
    initialized = .true.
    
    ! Specify the position and multiplicity on the grid of each variable
    do v = 1, N_SCALARS
       MULT(v)  = 1
       POSIT(v) = AT_NODE
    end do

    do v = N_SCALARS+1, N_VARS
       MULT(v)  = EDGE
       POSIT(v) = AT_EDGE
    end do

    end_pt = reshape ((/0,  0, 1, 0, 1, 1, 0, 0, 0, 0,  0, 1/), (/2, 2, 3/))
    opp_no = reshape ((/0, -1, 1, 1, 0, 1, 1, 0, 1, 1, -1, 0/), (/2, 2, 3/))

    no_adj_tri = reshape ((/0, 0, LORT, 0, -1, UPLT, -1, -1, LORT, -1, -1, UPLT, -1, 0, LORT, 0, 0, UPLT, 0, 0, LORT, &
         0, -1, UPLT, -1, -1, LORT, -1, -1, UPLT/), (/3, 10/))

    hex_sides = reshape ((/0, 0, RT, 0, 0, DG, 0, 0, UP, -1, 0, RT, -1, -1, DG, 0, -1, UP, 0, 0, RT, 0, 0, DG, 0, 0, UP, &
         -1, 0, RT/), (/3, 10/))

    hex_s_offs = (/2, 0, 4/)

    bfly_tri = reshape ((/-1, -1, LORT, 0, -1, LORT, 1, 0, UPLT, 0, 0, UPLT, 0, 1, LORT, -1, 0, LORT, 0, -1, UPLT, 1, 0, &
         UPLT, 0, 0, LORT, 0, 1, LORT, -1, 0, UPLT, -1, -1, UPLT/), (/3, 4, 3/))
    bfly_no  = reshape ((/-1, -1, 1, -1, 2, 1,  0, 1, 1, 2, -1,  0,  0, -1, 2, 1, 1, 0, 1, 2,  -1, 1, -1, -1/), (/2, 4, 3/))
    bfly_no2 = reshape ((/-3, -2, 1, -2, 3, 2, -1, 2, 1, 3, -3, -1, -1, -3, 3, 1, 2, -1, 2, 3, -2, 1, -2, -3/), (/2, 4, 3/))
    
    adj_tri  = reshape ((/0, -1, UPLT, 0, 0, LORT, 0, 0, UPLT, 0, 0, LORT, 0, 0, UPLT, -1, 0, LORT/), (/3, 2, 3/))

    nghb_pt  = reshape ((/1, 0, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 0, 1, 1, 0, 1, -1, 0/), (/2, 10/))

    ! Initialize values
    resume              = NONE
    cp_idx              = NONE
    err_restart         = 0
    istep               = 0
    istep_cumul         = 0
    iwrite              = 0
    time                = 0.0_8
    min_level           = DOMAIN_LEVEL+PATCH_LEVEL+1
    max_level           = min_level
    level_start         = min_level
    level_end           = level_start
    
    ! Default logical switches, most are reset in the input file
    adapt_dt            = .true.  ! dynamically adapt time step (T) or use time step based on initial conditions (F) 
    adapt_trend         = .false. ! adapt on trend (T) or on solution (F)
    compressible        = .true.  ! compressible equations (T) or Boussinesq incompressible (F)
    default_thresholds  = .true.  ! use default thresholds (T) or calculate dynamically (F)
    perfect             = .false. ! use perfect reconstruction criteria for wavelets and exact TRiSK operators (T) or less conservative wavetrisk version (F)
    rebalance           = .true.  ! rebalance computational load at each checkpoint if T
    remap               = .true.  ! remap Lagrangian coordinates (T) or no remapping (F)
    penalize            = .false. ! include penalization of topography if T
    uniform             = .true.  ! Uniform vertical grid in pressure (T) or hybrid (F)

    ! Default run values
    ! these parameters are typically reset in the input file, but are needed for compilation
    alpha               = 1d-4
    eta                 = 1d-2
    cfl_num             = 1.0_8
    C_visc              = 1d-2
    level_save          = level_start
    Laplace_order_init  = 0 ! 0 = no diffusion, 1 = Laplacian diffusion, 2 = second-order iterated Laplacian hyperdiffusion
    remapscalar_type    = "0" 
    remapvelo_type      = "0"
    timeint_type        = "RK45"
    iremap              = 1
    min_allowed_mass    = 1.0_8
    n_diffuse           = 1
    optimize_grid       = HR_GRID
    save_levels         = 1
    tol                 = 0.0_8
    zlevels             = 20
    
    ! Default physical parameters
    ! these parameters are typically reset in test case file, but are needed for compilation
    c_p            = 1004.64_8                   ! specific heat at constant pressure in joules per kilogram Kelvin
    c_v            = 717.6_8                     ! specfic heat at constant volume c_v = R_d - c_p
    grav_accel     = 9.80616_8
    p_top          = 0.0_8                     ! pressure at upper interface of top vertical layer (should be non-zero for Lin remapping)
    R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
    ref_density    = 1.0d3
    kappa          = R_d/c_p
    omega          = 7.292d-05
    radius         = 6371.22*KM
    p_0            = 1000.0d2

    visc_sclr = 0.0_8
    visc_divu = 0.0_8
    visc_rotu = 0.0_8
  end subroutine init_shared_mod

  real(8) function eps()
    eps = radius*1d-13
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
       exp__flush = 0.0_8
    end if
  end function exp__flush
end module shared_mod
