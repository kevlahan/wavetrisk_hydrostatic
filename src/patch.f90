module patch_mod
  use param_mod
  use shared_mod
  implicit none

  integer, parameter :: PATCH_SIZE = 2**PATCH_LEVEL
  integer, parameter :: LAST = PATCH_SIZE - 1
  integer, parameter :: LAST_BDRY = BDRY_THICKNESS - 1

  type Patch
     integer elts_start
     integer level
     integer, dimension(N_CHDRN) :: children
     integer, dimension(N_BDRY) :: neigh
     integer active
     logical deleted
  end type Patch

  type Bdry_Patch
     integer elts_start
     integer side
     integer neigh
  end type Bdry_Patch

  type Overl_Area
     real(8), dimension(4) :: a
     real(4), dimension(2) :: split  
  end type Overl_Area

  type Iu_Wgt
     real(4), dimension(9) :: enc  
  end type Iu_Wgt

  type RF_Wgt
     real(4), dimension(3) :: enc
  end type RF_Wgt
contains
  subroutine init_patch_mod
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_shared_mod
    initialized = .True.
  end subroutine init_patch_mod

  subroutine init_Patch (self, elts_start, level, chdrn, neigh)
    type(Patch) :: self
    integer     :: elts_start, chdrn, level, neigh

    self%elts_start = elts_start
    self%level = level
    self%deleted = .False.
  end subroutine init_Patch

  subroutine init_Bdry_Patch (self, elts_start, side, neigh)
    type(Bdry_Patch) :: self
    integer          :: elts_start, neigh, side

    self%elts_start = elts_start
    self%side = side
  end subroutine init_Bdry_Patch

  subroutine init_Overl_Area (self, areas)
    type(Overl_Area)      :: self
    real(8), dimension(8) :: areas

    self%a = areas(1:4)
    self%split = areas(5:6)
  end subroutine init_Overl_Area

  subroutine init_Iu_Wgt (self, wgt)
    type(Iu_Wgt)          :: self
    real(8), dimension(9) :: wgt

    self%enc = wgt
  end subroutine init_Iu_Wgt

  subroutine init_RF_Wgt (self, wgt)
    type(RF_Wgt) self
    real(4), dimension(3) :: wgt

    self%enc = wgt
  end subroutine init_RF_Wgt
end module patch_mod
