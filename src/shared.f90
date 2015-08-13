!
! BASICS OF CODE


! GRID STRUCTURE 
!
! The icosahedron is divided into ten regular lozange grids (with the exception of the poles), each lozange is then divided
! into N_SUB_DOM = 2**(2*DOMAIN_LEVEL) regular sub-domains with the number of sub-domains on each processor given by n_domain(rank+1)).


! BASIC GRID DATA TYPES
!
! Type(Coord_Array) has components elts(:)%length (where elts is Type(Coord) array of x%y%z coordinates on the sphere for each grid element and
! length is an integer giving the size of elts).
!
! Type(Float_Array) has components elts(:)%length (where elts is a real array of variable values on each grid element on each sub-domain and
! and length is an integer giving the size of elts). It is used for physical variables such as coriolis, kin_energy, divu, vort etc.
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


! Patch neighour directions
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
! Patch offset array offs(N_BDRY+1) contains the starting index for the current patch as offs(0) and the starting 
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


! DATA OUTPUT
!
! Data is written for plotting by routines io.f90/write_primal and io.f90/write_dual. 
! The full state of the simulation is saved in io.f90/dump_adapt_mpi and read in again by io.f90/load_adapt_mpi.

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

  ! neighbour indices for use in patch index arrays offs and dims 
  integer, parameter :: NORTH = 1
  integer, parameter :: EAST = 2
  integer, parameter :: SOUTH = 3
  integer, parameter :: WEST = 4
  integer, parameter :: NORTHEAST = 5
  integer, parameter :: SOUTHEAST = 6
  integer, parameter :: SOUTHWEST = 7
  integer, parameter :: NORTHWEST = 8

  integer, parameter :: N_BDRY = 8

  ! number of children nodes associated to each parent node
  integer, parameter :: N_CHDRN = 4 

  ! domain parameters
  integer, parameter :: N_ICOSAH_LOZANGE = 10 ! number of lozanges (coarse regular domains) in icosahedron
  integer, parameter :: N_SUB_DOM_PER_DIM = 2**DOMAIN_LEVEL ! number of subdomains per lozange in each direction
  integer, parameter :: N_SUB_DOM = N_SUB_DOM_PER_DIM**2 ! total number of sub-domains per lozange
  integer, allocatable :: n_domain(:) ! number of subdomains on each processor

  ! thickness of boundary overlaps between lozanges
  integer, parameter :: BDRY_THICKNESS = 3

  integer, parameter :: FROZEN = 32

  ! label for active points
  integer, parameter :: TOLRNZ = 16

  ! label for points whose flux can be obtained by restriction from fine level
  integer, parameter :: RESTRCT = 12

  ! label for adjacent zone points in either position (space) or scale
  integer, parameter :: ADJZONE = 8

  ! label for adjacent points in position (space) only 
  integer, parameter :: ADJSPACE = 14

  integer, parameter :: CONSIST = 4

  integer, parameter :: NODE = 3

  integer, parameter :: ZERO = 0 
  integer, parameter :: NONE = -1
  integer, parameter :: POLE = -2

  integer, parameter :: FALSE = 0
  integer, parameter :: TRUE = 1

  integer, parameter :: ON_LINE = 2
  integer, parameter :: INSIDE = 0
  integer, parameter :: OUTER1 = 1
  integer, parameter :: OUTER2 = 2
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

  ! first diagonal neighbours of hexagon points ??
  integer, parameter :: MP = 16
  integer, parameter :: PP = 17
  integer, parameter :: PM = 18
  integer, parameter :: MM = 19

  ! weights for various interpolation schemes
  integer, dimension(2,2,3) :: end_pt
  integer, dimension(2,2,3) :: opp_no
  integer, dimension(3,10)  :: no_adj_tri
  integer, dimension(3,10)  :: hex_sides
  integer, dimension(3)     :: hex_s_offs
  integer, dimension(3,4,3) :: bfly_tri
  integer, dimension(3,2,3) :: adj_tri
  integer, dimension(2,10)  :: nghb_pt
  integer, dimension(2,4,3) :: bfly_no
  integer, dimension(2,4,3) :: bfly_no2

  ! used in grid smoothing routine
  integer :: O2(2,3)
  data O2 /2,3, 3,1, 1,2/ 

  ! indices of prognostic variables in sol, trend etc
  integer, parameter :: S_MASS = 1, S_VELO = 2

  ! number of each variable per grid element (at hexagon nodes, triangle nodes, or edges) 
  integer, dimension(S_MASS:S_VELO) :: MULT 

  ! grid optimization choices
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, HR_GRID = 2

  ! basic grid parameters
  integer, parameter :: z_null = -1 ! place holder argument for functions not currently using z levels
  integer min_level, max_level ! minimum and maximum grid refinement levels in pseudo-horizontal directions
  integer zlevels ! number of levels in vertical direction
  integer level_start, level_end 

  real(8) threshold ! threshold level on wavelet coefficients for grid adaptation
  
  integer n_active(S_MASS:S_VELO) ! number of active points in each variable
  
  logical dynamic_adapt
  integer optimize_grid

  ! basic constants
  real(8), parameter :: MATH_PI = 3.14159265359_8
  integer, parameter :: DAY = 24*60*60

  ! simulation variables
  real(8) dt_write, time_end, time
  real(8) viscosity, friction_coeff, omega, radius, grav_accel, cfl_num, kmax
  integer istep, resume
  logical advect_only, wind_stress, bottom_friction

  ! for penalization boundary condition
  logical penalize
  real(8) alpha_m1, ieta

contains
  
  subroutine init_shared_mod()
    logical :: initialized = .False.

    if (initialized) return ! initialize only once

    MULT(S_MASS) = 1
    MULT(S_VELO) = EDGE

    end_pt = reshape((/0, 0, 1, 0, 1, 1, 0, 0, 0, 0, 0, 1/), (/2, 2, 3/))

    opp_no = reshape((/0, -1, 1, 1, 0, 1, 1, 0, 1, 1, -1, 0/), (/2, 2, 3/))

    no_adj_tri = reshape((/0, 0, LORT, 0, -1, UPLT, -1, -1, LORT, -1, -1, &
         UPLT, -1, 0, LORT, 0, 0, UPLT, 0, 0, LORT, 0, -1, UPLT, -1, -1, &
         LORT, -1, -1, UPLT/), (/3, 10/))

    hex_sides = reshape((/0, 0, RT, 0, 0, DG, 0, 0, UP, -1, 0, RT, -1, -1, &
         DG, 0, -1, UP, 0, 0, RT, 0, 0, DG, 0, 0, UP, -1, 0, RT/), (/3, &
         10/))

    hex_s_offs = (/2, 0, 4/)

    bfly_tri = reshape((/-1, -1, LORT, 0, -1, LORT, 1, 0, UPLT, 0, 0, UPLT, &
         0, 1, LORT, -1, 0, LORT, 0, -1, UPLT, 1, 0, UPLT, 0, 0, LORT, 0, &
         1, LORT, -1, 0, UPLT, -1, -1, UPLT/), (/3, 4, 3/))

    adj_tri = reshape((/0, -1, UPLT, 0, 0, LORT, 0, 0, UPLT, 0, 0, LORT, 0, &
         0, UPLT, -1, 0, LORT/), (/3, 2, 3/))

    nghb_pt = reshape((/1, 0, 1, 1, 0, 1, -1, 0, -1, -1, 0, -1, 1, 0, 1, 1, &
         0, 1, -1, 0/), (/2, 10/))

    bfly_no = reshape((/-1, -1, 1, -1, 2, 1, 0, 1, 1, 2, -1, 0, 0, -1, 2, 1, &
         1, 0, 1, 2, -1, 1, -1, -1/), (/2, 4, 3/))

    bfly_no2 = reshape((/-3, -2, 1, -2, 3, 2, -1, 2, 1, 3, -3, -1, -1, -3, 3, &
         1, 2, -1, 2, 3, -2, 1, -2, -3/), (/2, 4, 3/))

    ! earth parameters
    omega = 7.292e-05_8
    grav_accel = 9.80616_8
    radius = 6371220.0_8

    ! default values
    threshold = 0.0_8
    cfl_num = 1.0_8
    min_level = DOMAIN_LEVEL+PATCH_LEVEL+1
    max_level = min_level
    level_start = min_level
    level_end = level_start

    resume = NONE
    istep = 0
    time = 0.0_8
    viscosity = 0
    friction_coeff = 0
    wind_stress = .False.
    bottom_friction = .False.
    advect_only = .False.
    dynamic_adapt = .True.
    optimize_grid = NO_OPTIM

    ! for penalization boundary condition
    ieta = 1.0_8 / 0.1_8
    alpha_m1 = sqrt(1.0_8/ieta) - 1.0_8
    penalize = .False.

    initialized = .True.
  end subroutine init_shared_mod

  real(8) function eps()
    eps = radius*1e-13_8
  end function eps

  integer function max_nodes_per_level(lev, entity)
    integer lev
    integer, optional :: entity
    max_nodes_per_level = 10*4**lev
    if (present(entity)) then 
       max_nodes_per_level = max_nodes_per_level*entity
    else ! mass node
       max_nodes_per_level = max_nodes_per_level+2
    end if
  end function max_nodes_per_level

  real(8) function exp__flush(x)
    real(8) :: x
    if (x .gt. -100.0_8) then
       exp__flush = exp(x)
    else
       exp__flush = 0.0_8
    end if
  end function exp__flush
  
end module shared_mod
