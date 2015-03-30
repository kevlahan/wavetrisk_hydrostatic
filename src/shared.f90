module shared_mod
  use param_mod
  implicit none
  logical :: DEB_ON
  integer, parameter :: NORTH = 1
  integer, parameter :: EAST = 2
  integer, parameter :: SOUTH = 3
  integer, parameter :: WEST = 4
  integer, parameter :: NORTHEAST = 5
  integer, parameter :: SOUTHEAST = 6
  integer, parameter :: SOUTHWEST = 7
  integer, parameter :: NORTHWEST = 8
  integer, parameter :: JPLUS = 1
  integer, parameter :: IPLUS = 2
  integer, parameter :: JMINUS = 3
  integer, parameter :: IMINUS = 4
  integer, parameter :: IJPLUS = 5
  integer, parameter :: IPLUSJMINUS = 6
  integer, parameter :: IJMINUS = 7
  integer, parameter :: IMINUSJPLUS = 8
  integer, parameter :: N_BDRY = 8
  integer, parameter :: N_CHDRN = 4
  real(8), parameter :: MATH_PI = 3.14159265359_8
  integer, parameter :: N_ICOSAH_LOZANGE = 10
  integer, parameter :: BDRY_THICKNESS = 3
  integer, parameter :: EDGE = 3
  integer, parameter :: TRIAG = 2
  integer, parameter :: FROZEN = 32
  integer, parameter :: TOLLRNZ = 16
  integer, parameter :: ADJSPACE = 14
  integer, parameter :: RESTRCT = 12
  integer, parameter :: ADJZONE = 8
  integer, parameter :: CONSIST = 4
  integer, parameter :: RT = 0
  integer, parameter :: DG = 1
  integer, parameter :: UP = 2
  integer, parameter :: NODE = 3
  integer, parameter :: LORT = 0
  integer, parameter :: UPLT = 1
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
  integer, parameter :: CMM = 12
  integer, parameter :: CPP = 13
  integer, parameter :: WMMM = 12
  integer, parameter :: WPPP = 13
  integer, parameter :: VMPP = 14
  integer, parameter :: VPMM = 15
  integer, parameter :: MP = 16
  integer, parameter :: PP = 17
  integer, parameter :: PM = 18
  integer, parameter :: MM = 19
  integer, parameter :: S_HEIGHT = 1, AT_NODE = 1
  integer, parameter :: S_VELO = 2, AT_EDGE = 2
  integer, parameter :: NO_OPTIM = 0, XU_GRID = 1, HR_GRID = 2
  integer, dimension(S_HEIGHT:S_VELO) :: MULT
  integer, dimension(2,2,3) :: end_pt
  integer, dimension(2,2,3) :: opp_no
  integer, dimension(3,10) :: no_adj_tri
  integer, dimension(3,10) :: hex_sides
  integer, dimension(3) :: hex_s_offs
  integer, dimension(3,4,3) :: bfly_tri
  integer, dimension(3,2,3) :: adj_tri
  integer, dimension(2,10) :: nghb_pt
  integer, dimension(2,4,3) :: bfly_no
  integer, dimension(2,4,3) :: bfly_no2
  integer, parameter :: N_SUB_DOM_PER_DIM = 2**DOMAIN_LEVEL
  integer, parameter :: N_SUB_DOM = N_SUB_DOM_PER_DIM**2
  integer, parameter :: DAY = 24*60*60
  integer :: O2(2,3)
  integer, allocatable :: n_domain(:)
  integer n_active(S_HEIGHT:S_VELO)
  integer min_level, max_level
  integer level_start, level_end
  real(8) threshold
  real(8) time_end
  real(8) viscosity
  real(8) omega, radius, grav_accel, cfl_num, kmax
  real(8) dt_write
  real(8) friction_coeff
  integer resume
  integer istep
  real(8) time
  logical wind_stress
  logical bottom_friction
  logical advect_only
  logical dynamic_adapt
  integer optimize_grid
  ! for penalization boundary condition
  logical penalize
  real(8) alpha_m1, ieta

  data O2 /2,3, 3,1, 1,2/ 

contains
  subroutine init_shared_mod()
      logical :: initialized = .False.
      if (initialized) return ! initialize only once
      MULT(S_HEIGHT) = 1
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
      threshold = 0.0_8

      ! earth
      omega = 7.292e-05_8
      grav_accel = 9.80616_8
      radius = 6371220.0_8

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
  end function

  integer function max_nodes_per_level(lev, entity)
      integer lev
      integer, optional :: entity
      max_nodes_per_level = 10*4**lev
      if (present(entity)) then 
          max_nodes_per_level = max_nodes_per_level*entity
      else ! height node
          max_nodes_per_level = max_nodes_per_level+2
      end if
  end function

  real(8) function exp__flush(x)
      real(8) :: x
      if (x .gt. -100.0_8) then
          exp__flush = exp(x)
      else
          exp__flush = 0.0_8
      end if
  end function
end module shared_mod
