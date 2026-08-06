module patch_mod
  
  use kind_mod,   only : dp
  use shared_mod, only : BDRY_THICKNESS, N_BDRY, N_CHDRN, PATCH_LEVEL, init_shared_mod
  
  implicit none

  private
  public :: LAST, LAST_BDRY, PATCH_SIZE
  public :: Bdry_Patch, Overl_Area, Iu_Wgt, Patch, RF_Wgt
  public :: init_Iu_Wgt, init_Overl_Area, init_Bdry_Patch, init_Patch, init_patch_mod, init_RF_Wgt

  integer, parameter :: PATCH_SIZE = 2**PATCH_LEVEL
  integer, parameter :: LAST       = PATCH_SIZE - 1
  integer, parameter :: LAST_BDRY  = BDRY_THICKNESS - 1

  type Patch
     integer                     :: elts_start
     integer                     :: level
     integer, dimension(N_CHDRN) :: children
     integer, dimension(N_BDRY)  :: neigh
     integer                     :: active
     logical                     :: deleted
  end type Patch

  type Bdry_Patch
     integer :: elts_start
     integer :: side
     integer :: neigh
  end type Bdry_Patch

  type Overl_Area
     real(dp), dimension(4) :: a
     real(dp), dimension(2) :: split  
  end type Overl_Area

  type Iu_Wgt
     real(dp), dimension(9) :: enc  
  end type Iu_Wgt

  type RF_Wgt
     real(dp), dimension(3) :: enc
  end type RF_Wgt
contains
  subroutine init_patch_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_shared_mod
    initialized = .True.
  end subroutine init_patch_mod

  subroutine init_Patch (self, elts_start, level, chdrn, neigh)
    ! Initializes new patch
    implicit none
    type(Patch), intent(inout) :: self
    integer,     intent(in)    :: elts_start, chdrn, level, neigh

    self%elts_start = elts_start
    self%level      = level
    self%deleted     = .false.
  end subroutine init_Patch

  subroutine init_Bdry_Patch (self, elts_start, side, neigh)
    implicit none
    type(Bdry_Patch), intent(inout) :: self
    integer,          intent(in)    :: elts_start, neigh, side

    self%elts_start = elts_start
    self%side = side
  end subroutine init_Bdry_Patch

  subroutine init_Overl_Area (self, areas)
    implicit none
    type(Overl_Area),       intent(inout) :: self
    real(dp), dimension(8), intent(in)    :: areas

    self%a = areas(1:4)
    self%split = areas(5:6)
  end subroutine init_Overl_Area

  subroutine init_Iu_Wgt (self, wgt)
    implicit none
    type(Iu_Wgt),           intent(inout) :: self
    real(dp), dimension(9), intent(in)    :: wgt

    self%enc = wgt
  end subroutine init_Iu_Wgt
  
  subroutine init_RF_Wgt (self, wgt)
    implicit none
    type(RF_Wgt),           intent(inout) :: self
    real(dp), dimension(3), intent(in)    :: wgt

    self%enc = wgt
  end subroutine init_RF_Wgt
end module patch_mod
