!
! BASICS OF CODE
!
!
! GRID STRUCTURE 
!
! The icosahedron is divided into a network of ten regular lozenge grids (with the exception of the POLES),
! each lozenge is then divided into N_SUB_DOM = 4**DOMAIN_LEVEL regular sub-domains with the number of sub-domains on each processor given by n_domain(rank+1).
!
! There are two special nodes called POLEs. One connects the five lozenge vertices at the top of the network and the other connects the 
! five lozenge grid vertices at the bottom of the network. The network is oriented so the POLEs are at the geographic North and South poles.  
!
! The coarsest level Jmin = DOMAIN_LEVEL + PATCH_LEVEL + 1 with PATCH_LEVEL>=2. The geometry of this coarsest level may be optimized by reading 
! in a Heikes & Randall (1995) grid of resolution Jmin-1 or using the Xu (2006) smoothing algorithm.
! The size of patches is 2**PATCH_LEVEL (e.g. if PATCH_LEVEL = 2 then the patches are 4x4).
! The larger PATCH_LEVEL is the fewer levels of tree structure must be traversed before reaching a uniform grid (which will in general contain inactive nodes).
! Each computation patch element is made up of one node for scalars and three edges for vectors U, V,W (see also GRID ELEMENTS below).
!
! PARALLEL EXECUTION
!
! If compiled with ARCH = mpi or mpi-lb (load balancing) the code can be run in parallel using mpi.  With option mpi-lb the computational load is re-balanced
! each time a checkpoint is saved.
! 
! The maximum number of computational cores must be less than or equal to the number of domains, i.e. Ncore <= 10*4**DOMAIN_LEVEL.
! Since DOMAIN_LEVEL = Jmin - (PATCH_LEVEL+1) larger Jmin or smaller PATCH_LEVEL allows more cores to be used.  For inhomogeneous problems
! (i.e. unbalanced) it is best to set PATCH_LEVEL=2 and use a larger Jmin while for homogeneous (i.e. well-balanced) problems it is more efficient to choose
! a larger PATCH_LEVEL. For example, if Jmin = 7 and PATCH_LEVEL = 4 then up to 160 cores may be used, while if PATCH_LEVEL = 2 then up to 2560 cores may be used.
!
! BASIC GRID DATA TYPES
!
! Type(Coord_Array) has components elts(:)%length (where elts is Type(Coord) array of x%y%z coordinates on the sphere for each grid element and
! length is an integer giving the size of elts).
!
! Type(Float_Array) has components elts(:)%length (where elts is a real array of variable values on each grid element on each sub-domain and
! and length is an integer giving the size of elts). It is used for physical variables such as coriolis, divu, vort etc.
!
! Type(Int_Array) has components elts(:)%length (where elts is an integer array of element indices on each sub-domain and
! and length is an integer giving the size of elts). It is used for masks, levels and parallel communication variables.
!
! Type(Float_Field) has components data(:)%bdry_uptodate%pos (where data is a Float_Array, bdry_uptodate is a logical and pos is an integer)
! and is used for equation variables such as sol, trend, wav_coeff, sca_coeff.
!
! Type(Domain) has many components defining the grid and is the size of the total number of sub-domains on each processor n_domain(rank+1).
! It is used for the variable grid(:).
!
! Various other dynamical data types and subroutines for allocating them are defined in dyn_array.f90.


! GRID ELEMENTS
!
! The triangular primal grid element (one hexagonal node H, three edges UP, DG, RT and two triangular nodes UPLT, LORT) is

!             ------------ 
!             \           / \ 
!              \  UPLT   /   \ 
!               \       /     \
!               UP     DG      \   
!                 \   /   LORT  \
!                  \ /           \
!                   H ------RT---- 


! The hexagonal dual grid element: one hexagonal node H, three adjacent edges UP, DG, RT and two triangular nodes O
! (rotated clockwise 30 degrees wrt primal grid above) is

!
!
!             -----UP---- UPLT 
!           /               \ 
!          /                 \ 
!         /                  DG
!        /                     \   
!       /                       \
!                   H           LORT
!       \                       /
!        \                     /
!         \                   RT
!          \                 / 
!           \               /
!             -------------


! Patch neighour directions (based on regular coordinates i,j).
! Note that within a patch similar notation is used for shifts in i and j (e.g. shift i-1, j+1 is denoted idNW).
!                   
! ------------- ------------- ------------- 
!\              \             \            \
! \              \             \            \
!  \              \             \            \
!   \  NORTHWEST   \    NORTH    \ NORTHEAST  \  
!    \              \             \            \
!     \              \             \            \
!       -------------  ------------- ------------- 
!       \              \             \            \
!        \              \             \            \
!         \              \             \            \  
!          \     WEST     \      0      \    EAST    \ 
!           \              \             \            \
!            \              \             \            \
!              -------------  ------------- ------------- 
!              \              \             \            \
!               \              \             \            \ 
!                \              \             \            \
!                 \   SOUTHWEST  \  SOUTH      \  SOUTHEAST \
!                  \              \             \            \
!                   \              \             \            \
!                     -------------   ------------ ------------ 
!

!
!
! INDEXING OF GRID ELEMENTS AND NEIGHBOURS
!
! Quantities (e.g. mass, velocities) are all stored in a single array, whose elements are organized in patches.
! Each patch has regular array coordinates (i,j).
!
! Patch offset array offs(N_BDRY+1) contains the starting index in the single array for the current patch as offs(0). The starting 
! indices for neighbouring patches are given by offs(NORTH) etc.
!
! Patch dimension array dims(2, N_BDRY+1) gives the dimensions of the current patch as dims(2,0) and neighbouring patches as dims(2, NORTH) etc.
!
! Function id = idx(i,j,offs,dims) returns the element index for the hexagon node with coordinates (i,j) on the patch selected by offs with dimensions dim. 
!
! The components of the grid elements are then found as:
!
! %elts(id+1)         - the one grid element hexagon node (e.g. masses)
! %elts(EDGE*id+e+1)  - the three grid element edges RT, DG, UP, where e = 0,1,2 (e.g. velocities, fluxes)
! %elts(TRIAG*id+t+1) - the two grid element triangles LORT, UPLT, where t = 0,1 (e.g. circulation)
!
! Wavelet coefficients are stored at the SAME nodes/edges as the nodes/edges they characterize.

! DATA OUTPUT
!
! Data is written for plotting by routines io.f90/write_primal and io.f90/write_dual. 
! The full state of the simulation is saved in io.f90/dump_adapt_mpi and read in again by io.f90/load_adapt_mpi.

! CALCULATIONS ON ADAPTED GRID
!
! By default fields are calculated and operators are applied on the ENTIRE grid (including at nodes where the result can be obtained by  restriction 
! indicated by mask=12 and adjacent zone nodes indicated by mask=8) and the results are then over-written by correct values.  
! (Note that the solution in the adjacent zone at fine scale j+1 are found from values at coarse scale j so the values calculated at scale j+1 are not actually used.)
! This means that some operations on the entire grid could produce intermediate overflows, inf/NaN, or invalid indices ue to incorrect values at these nodes or their neighbours.
! Functions and subroutines should take this into account.  Similarly, circulation, vorticity and qe are first computed incorrectly at pentagon points
! in step1 and then corrected in post_step1.

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
  integer, parameter :: S_MASS = 1, S_TEMP = 2, S_VELO = 3

  ! Number of each variable per grid element (at hexagon nodes, triangle nodes, or edges) 
  integer, dimension(S_MASS:S_VELO) :: MULT

  ! Location of each variable on the grid (at an edge or at a node)
  integer, dimension(S_MASS:S_VELO) :: POSIT

  ! Grid optimization choices
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, HR_GRID = 2

  ! Basic grid parameters
  integer, parameter :: z_null = -1 ! place holder argument for functions not currently using z levels
  integer :: min_level, max_level ! minimum and maximum grid refinement levels in pseudo-horizontal directions
  integer :: zlevels ! number of levels in vertical direction
  integer :: save_levels ! number of vertical levels to save
  integer :: level_start, level_end, level_save

  real(8) :: tol ! relative tolerance for all variables

  integer, dimension(AT_NODE:AT_EDGE) :: n_active ! number of active points at grid locations (node and edge)

  integer :: optimize_grid

  ! Basic constants (assumes basic unit of time is seconds)
  integer, parameter :: SECOND = 1
  integer, parameter :: MINUTE = 60*SECOND
  integer, parameter :: HOUR = 60*MINUTE
  integer, parameter :: DAY = 24*HOUR
  integer, parameter :: WEEK = 7*DAY
  real(8), parameter :: MATH_PI = acos(-1.0_8)
  
  ! Simulation variables
  integer                                       :: cp_idx, err_restart, istep, istep_cumul, iwrite, n_diffuse
  integer                                       :: resume, Laplace_order, Laplace_order_init
  integer(8)                                    :: itime
  integer, parameter                            :: max_restart = 10 ! Maximum allowed number of restarts after negative mass event
  
  real(8)                                       :: dt, dt_init, dt_write, dx_min, dx_max, time_end, time
  real(8)                                       :: omega, radius, grav_accel, cfl_num, kmax, ref_density
  real(8)                                       :: viscosity_rotu, viscosity_mass, viscosity_temp
  real(8)                                       :: p_0, p_top, gamma, gk, kappa, c_p, c_v, R_d, wave_speed
  real(8)                                       :: min_mass, min_allowed_mass
  real(8), dimension(:),         allocatable    :: pressure_save
  real(8), dimension(:),         allocatable    :: a_vert, b_vert, a_vert_mass, b_vert_mass, viscosity_divu
  real(8), dimension(:,:),       allocatable    :: threshold
  real(8), dimension(3)                         :: L_diffusion
  real(8), dimension (10*2**(2*DOMAIN_LEVEL),3) :: nonunique_pent_locs
  real(8), dimension (12,3)                     :: unique_pent_locs

  character(255)                                :: run_id, test_case
  
  logical :: adapt_dt, adapt_trend, compressible, default_thresholds, lagrangian_vertical, perfect, remap, uniform
contains
  subroutine init_shared_mod
    logical :: initialized = .false.

    if (initialized) return ! Initialize only once
    initialized = .true.
    
    ! Specify the multiplicity per grid element of each quantity
    MULT(S_MASS) = 1
    MULT(S_TEMP) = 1
    MULT(S_VELO) = EDGE

    ! Specify the position on the grid of each quantity
    POSIT(S_MASS) = AT_NODE
    POSIT(S_TEMP) = AT_NODE
    POSIT(S_VELO) = AT_EDGE

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
    adapt_dt            = .true. ! Dynamically adapt time step (T) or use time step based on initial conditions (F) 
    adapt_trend         = .true. ! Adapt on trend (T) or on solution (F)
    compressible        = .true. ! Compressible equations (T) or Boussinesq incompressible (F)
    default_thresholds  = .true. ! Use default thresholds (T) or calculate dynamically (F)
    perfect             = .true. ! Use perfect reconstruction criteria for wavelets and exact TRiSK operators (T) or less conservative wavetrisk version (F)
    remap               = .true. ! Remap Lagrangian coordinates (T) or no remapping (F)
    lagrangian_vertical = .true. ! Lagrangian or mass based vertical coordinates (only option implement is T, mass-based coordinates not implemented)
    uniform             = .true. ! Uniform vertical grid in pressure (T) or hybrid (F)

    ! Default run values
    ! these parameters are typically reset in the input file, but are needed for compilation
    cfl_num             = 1.5_8
    level_save          = level_start
    Laplace_order_init  = 0                      ! 0 = no diffusion, 1 = Laplacian diffusion, 2 = second-order iterated Laplacian hyperdiffusion
    min_allowed_mass    = 0.3_8
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
    p_top          = 0.0_8                       ! pressure at upper interface of top vertical layer
    R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
    ref_density    = 1.0_8
    kappa          = R_d/c_p
    omega          = 7.292d-05
    radius         = 6371220.0_8
    p_0            = 1000.0d2
    viscosity_rotu = 0.0_8
    viscosity_mass = 0.0_8
    viscosity_temp = 0.0_8
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
