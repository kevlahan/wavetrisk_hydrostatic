module dyn_arrays
  use kind_mod,   only : dp
  use shared_mod, only : Areas, Coord, N_BDRY, N_CHDRN, ORIGIN
  
  use patch_mod, only : BDRY_patch, Iu_wgt, RF_wgt, Overl_Area, patch
  
  implicit none

  private
  public :: append, extend, init
  public :: Areas_Array, Bdry_Patch_Array, Coord_Array, Float_Array, Int_Array, Iu_Wgt_Array
  public :: Logical_Array, Overl_Area_Array, Patch_Array, RF_Wgt_Array, Topo_Array
  
  type Int_Array
     integer, dimension(:), allocatable :: elts
     integer                            :: length
  end type Int_Array
  
  type Float_Array
     real(dp), dimension(:), allocatable :: elts
     integer                             :: length
  end type Float_Array

  type Logical_Array
     logical, dimension(:), allocatable :: elts
     integer                            :: length
  end type Logical_Array

  type Topo_Array
     real(dp),    dimension(:), allocatable :: elts
     type(Coord), dimension(:), allocatable :: node
  end type Topo_Array

  type Coord_Array
     type(Coord), dimension(:), allocatable :: elts
     integer                                :: length
  end type Coord_Array

  type Areas_Array
     type(Areas), dimension(:), allocatable :: elts
     integer                                :: length
  end type Areas_Array

  type Overl_Area_Array
     type(Overl_Area), dimension(:), allocatable :: elts
     integer                                     :: length
  end type Overl_Area_Array

  type Iu_Wgt_Array
     type(Iu_Wgt), dimension(:), allocatable :: elts
     integer                                 :: length
  end type Iu_Wgt_Array

  type RF_Wgt_Array
     type(RF_Wgt), dimension(:), allocatable :: elts
     integer                                 :: length
  end type RF_Wgt_Array

  type Patch_Array
     type(Patch), dimension(:), allocatable :: elts
     integer                                :: length
  end type Patch_Array

  type Bdry_Patch_Array
     type(Bdry_Patch), dimension(:), allocatable :: elts
     integer                                     :: length
  end type Bdry_Patch_Array

  interface init
     module procedure init_Float_Array, init_Int_Array, init_Logical_Array, &
           init_Coord_Array, init_Areas_Array, &
          init_Overl_Area_Array, init_Iu_Wgt_Array, init_RF_Wgt_Array, &
          init_Patch_Array, init_Bdry_Patch_Array
  end interface init

  interface append
     module procedure append_Int_Array, append_Float_Array, append_Coord_Array, &
          append_Areas_Array, &
          append_Overl_Area_Array, append_Iu_Wgt_Array, append_RF_Wgt_Array, &
          append_Patch_Array, append_Bdry_Patch_Array
  end interface append

  interface extend
     module procedure extend_Int_Array, extend_Float_Array, extend_Coord_Array, &
          extend_Areas_Array, &
          extend_Overl_Area_Array, extend_Iu_Wgt_Array, extend_RF_Wgt_Array, &
          extend_Patch_Array, extend_Bdry_Patch_Array
  end interface extend

  interface dbl_alloc
     module procedure dbl_alloc_Int_Array, dbl_alloc_Float_Array, &
          dbl_alloc_Coord_Array, &
          dbl_alloc_Areas_Array, dbl_alloc_Overl_Area_Array, &
          dbl_alloc_Iu_Wgt_Array, dbl_alloc_RF_Wgt_Array, &
          dbl_alloc_Patch_Array, dbl_alloc_Bdry_Patch_Array
  end interface dbl_alloc
  
contains


  subroutine init_Float_Array (arr, N)
    implicit none
    type(Float_Array), intent(out) :: arr
    integer,           intent(in)    :: N

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts = 0.0_dp
  end subroutine init_Float_Array


  subroutine init_Int_Array (arr, N)
    implicit none
    type(Int_Array), intent(out) :: arr
    integer,         intent(in)  :: N

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts = 0
  end subroutine init_Int_Array


  subroutine init_Logical_Array (arr, N)
    implicit none
    type(Logical_Array), intent(out) :: arr
    integer,             intent(in)  :: N

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts = .false.
  end subroutine init_Logical_Array


  subroutine init_Coord_Array (arr, N)
    implicit none
    type(Coord_Array), intent(out) :: arr
    integer,           intent(in)  :: N

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts = ORIGIN
  end subroutine init_Coord_Array


  subroutine init_Areas_Array (arr, N)
    type(Areas_Array), intent(out) :: arr
    integer,           intent(in)  :: N

    integer :: i

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    do i = 1, 6
       arr%elts%part(i) = 0.0_dp
    end do
    arr%elts%hex_inv = 0.0_dp
  end subroutine init_Areas_Array


  subroutine init_Overl_Area_Array (arr, N)
    implicit none
    type(Overl_Area_Array), intent(out) :: arr
    integer,                intent(in)  :: N

    integer :: i

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    do i = 1, 4
       arr%elts%a(i) = 0.0_dp
    end do
    do i = 1, 2
       arr%elts%split(i) = 0.0_dp
    end do
  end subroutine init_Overl_Area_Array


  subroutine init_Iu_Wgt_Array (arr, N)
    implicit none
    type(Iu_Wgt_Array), intent(out) :: arr
    integer,            intent(in)  :: N

    integer :: i

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    do i = 1, 9
       arr%elts%enc(i) = 0.0_dp
    end do
  end subroutine init_Iu_Wgt_Array


  subroutine init_RF_Wgt_Array (arr, N)
    implicit none
    type(RF_Wgt_Array), intent(out) :: arr
    integer,            intent(in)  :: N

    integer :: i

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    do i = 1, 3
       arr%elts%enc(i) = 0.0_dp
    end do
  end subroutine init_RF_Wgt_Array


  subroutine init_Patch_Array (arr, N)
    implicit none
    type(Patch_Array), intent(out) :: arr
    integer,           intent(in)  :: N

    integer :: i

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts%elts_start = 0
    arr%elts%level      = 0
    do i = 1, N_CHDRN
       arr%elts%children(i) = 0
    end do
    do i = 1, N_BDRY
       arr%elts%neigh(i) = 0
    end do
    arr%elts%active  = 0
    arr%elts%deleted = .false.
  end subroutine init_Patch_Array


  subroutine init_Bdry_Patch_Array (arr, N)
    implicit none
    type(Bdry_Patch_Array), intent(out) :: arr
    integer,                intent(in)  :: N

    arr%length = N
    allocate(arr%elts(max(N,1))) ! min. 1 -> no 0 alloc
    arr%elts%elts_start = 0
    arr%elts%side       = 0
    arr%elts%neigh      = 0
  end subroutine init_Bdry_Patch_Array


  subroutine append_Int_Array (arr, item)
    type(Int_Array), intent(inout) :: arr
    integer,         intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Int_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Int_Array


  subroutine append_Float_Array (arr, item)
    implicit none
    type(Float_Array), intent(inout) :: arr
    real(dp),          intent(in)    :: item

    if (arr%length ==  size(arr%elts)) call dbl_alloc_Float_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Float_Array


  subroutine append_Coord_Array (arr, item)
    implicit none
    type(Coord_Array), intent(inout) :: arr
    type(Coord),       intent(in)    ::  item

    if (arr%length == size(arr%elts)) call dbl_alloc_Coord_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Coord_Array


  subroutine append_Areas_Array (arr, item)
    implicit none
    type(Areas_Array), intent(inout) :: arr
    type(Areas),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Areas_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Areas_Array


  subroutine append_Overl_Area_Array (arr, item)
    implicit none
    type(Overl_Area_Array), intent(inout) :: arr
    type(Overl_Area),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Overl_Area_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Overl_Area_Array


  subroutine append_Iu_Wgt_Array (arr, item)
    implicit none
    type(Iu_Wgt_Array), intent(inout) :: arr
    type(Iu_Wgt),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Iu_Wgt_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Iu_Wgt_Array


  subroutine append_RF_Wgt_Array (arr, item)
    implicit none
    type(RF_Wgt_Array), intent(inout) :: arr
    type(RF_Wgt),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_RF_Wgt_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_RF_Wgt_Array


  subroutine append_Patch_Array (arr, item)
    implicit none
    type(Patch_Array), intent(inout) :: arr
    type(Patch),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Patch_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Patch_Array


  subroutine append_Bdry_Patch_Array(arr, item)
    implicit none
    type(Bdry_Patch_Array), intent(inout) :: arr
    type(Bdry_Patch),       intent(in)    :: item

    if (arr%length == size(arr%elts)) call dbl_alloc_Bdry_Patch_Array (arr, arr%length + 1)

    arr%length = arr%length + 1
    arr%elts(arr%length) = item
  end subroutine append_Bdry_Patch_Array


  subroutine extend_Int_Array (arr, N, items)
    implicit none
    type(Int_Array), intent(inout) :: arr
    integer,         intent(in)    :: items, N

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Int_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Int_Array


  subroutine extend_Float_Array (arr, N, items)
    type(Float_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    real(dp),          intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Float_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Float_Array


  subroutine extend_Coord_Array (arr, N, items)
    implicit none
    type(Coord_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    type(Coord),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Coord_Array(arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Coord_Array


  subroutine extend_Areas_Array (arr, N, items)
    implicit none
    type(Areas_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    type(Areas),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Areas_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Areas_Array


  subroutine extend_Overl_Area_Array (arr, N, items)
    implicit none
    type(Overl_Area_Array), intent(inout) :: arr
    integer,                intent(in)    :: N
    type(Overl_Area),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Overl_Area_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Overl_Area_Array


  subroutine extend_Iu_Wgt_Array (arr, N, items)
    implicit none
    type(Iu_Wgt_Array), intent(inout) :: arr
    integer,            intent(in)    :: N
    type(Iu_Wgt),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Iu_Wgt_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Iu_Wgt_Array


  subroutine extend_RF_Wgt_Array (arr, N, items)
    implicit none
    type(RF_Wgt_Array), intent(inout) :: arr
    integer,            intent(in)    :: N
    type(RF_Wgt),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_RF_Wgt_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_RF_Wgt_Array


  subroutine extend_Patch_Array (arr, N, items)
    implicit none
    type(Patch_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    type(Patch),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Patch_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Patch_Array


  subroutine extend_Bdry_Patch_Array (arr, N, items)
    implicit none
    type(Bdry_Patch_Array), intent(inout) :: arr
    integer,                intent(in)    :: N
    type(Bdry_Patch),       intent(in)    :: items

    if (arr%length + N > size(arr%elts)) call dbl_alloc_Bdry_Patch_Array (arr, arr%length + N)

    arr%elts(arr%length+1:arr%length+N) = items
    arr%length = arr%length + N
  end subroutine extend_Bdry_Patch_Array


  subroutine dbl_alloc_Int_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Int_Array), intent(inout) :: arr
    integer,         intent(in)    :: N

    integer                            :: ierr
    integer, dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N),stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    tmparr = 0
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    arr%elts = 0
    
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Int_Array


  subroutine dbl_alloc_Float_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Float_Array), intent(inout) :: arr
    integer,           intent(in)    :: N

    integer                             :: ierr
    real(dp), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    tmparr = 0.0_dp
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,*) "ERROR: not enough memory"
       stop
    endif
    arr%elts = 0.0_dp
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Float_Array


  subroutine dbl_alloc_Coord_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Coord_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    
    integer                                :: ierr
    type(Coord), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,*) "ERROR: not enough memory"
       stop
    endif
    tmparr = ORIGIN
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    arr%elts = ORIGIN
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Coord_Array


  subroutine dbl_alloc_Areas_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Areas_Array), intent(inout) :: arr
    integer,           intent(in)    :: N

    integer                                :: i, ierr
    type(Areas), dimension(:), allocatable :: tmparr
    
    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 6
      tmparr%part(i) = 0.0_dp
    end do
    tmparr%hex_inv = 0.0_dp
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 6
       arr%elts%part(i) = 0.0_dp
    end do
    arr%elts%hex_inv = 0.0_dp
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Areas_Array


  subroutine dbl_alloc_Overl_Area_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Overl_Area_Array), intent(inout) :: arr
    integer,                intent(in)    :: N
    
    integer                                     :: i, ierr
    type(Overl_Area), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 4
       tmparr%a(i) = 0.0_dp
    end do
    do i = 1, 2
       tmparr%split(i) = 0.0_dp
    end do
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 4
       arr%elts%a(i) = 0.0_dp
    end do
    do i = 1, 2
       arr%elts%split(i) = 0.0_dp
    end do
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Overl_Area_Array


  subroutine dbl_alloc_Iu_Wgt_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Iu_Wgt_Array), intent(inout) :: arr
    integer,            intent(in)    :: N

    integer                                 :: i, ierr
    type(Iu_Wgt), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 9
       tmparr%enc(i) = 0.0_dp
    end do
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 9
       arr%elts%enc(i) = 0.0_dp
    end do
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Iu_Wgt_Array


  subroutine dbl_alloc_RF_Wgt_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(RF_Wgt_Array), intent(inout) :: arr
    integer,            intent(in)    :: N
    
    integer                                 :: i, ierr
    type(RF_Wgt), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 3
       tmparr%enc(i) = 0.0_dp
    end do
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    do i = 1, 3
       arr%elts%enc(i) = 0.0_dp
    end do
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_RF_Wgt_Array


  subroutine dbl_alloc_Patch_Array(arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Patch_Array), intent(inout) :: arr
    integer,           intent(in)    :: N
    
    integer                                :: i, ierr
    type(Patch), dimension(:), allocatable :: tmparr

    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    tmparr%elts_start = 0
    tmparr%level      = 0
    do i = 1, N_CHDRN
       tmparr%children(i) = 0
    end do
    do i = 1, N_BDRY
       tmparr%neigh(i) = 0
    end do
    tmparr%active  = 0
    tmparr%deleted = .false.
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    arr%elts%elts_start = 0
    arr%elts%level      = 0
    do i = 1, N_CHDRN
       arr%elts%children(i) = 0
    end do
    do i = 1, N_BDRY
       arr%elts%neigh(i) = 0
    end do
    arr%elts%active  = 0
    arr%elts%deleted = .false.
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Patch_Array


  subroutine dbl_alloc_Bdry_Patch_Array (arr, N)
    ! double allocated memory to avoid frequent reallocation
    implicit none
    type(Bdry_Patch_Array), intent(inout) :: arr
    integer,                intent(in)    :: N

    integer                                     :: ierr
    type(Bdry_Patch), dimension(:), allocatable :: tmparr
    
    allocate (tmparr(2*N), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif
    tmparr%elts_start = 0
    tmparr%side       = 0
    tmparr%neigh      = 0
    tmparr(1:size(arr%elts)) = arr%elts

    if (allocated (arr%elts)) deallocate (arr%elts)
    allocate (arr%elts(size(tmparr)), stat=ierr)
    if (ierr /= 0) then
       write(0,'(A)') "ERROR: not enough memory"
       stop
    endif 
    arr%elts%elts_start = 0
    arr%elts%side       = 0
    arr%elts%neigh      = 0
    arr%elts(1:size(arr%elts)) = tmparr(1:size(arr%elts))
    deallocate (tmparr)
  end subroutine dbl_alloc_Bdry_Patch_Array


end module dyn_arrays
