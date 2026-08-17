module parallel_block_mod

  use, intrinsic :: iso_fortran_env, only : int8, int64

  use kind_mod,   only : dp
  use shared_mod, only : Coord, MULT, N_BDRY, S_VELO, scalars, zmin, zmax
  use patch_mod,  only : Patch, PATCH_SIZE

  implicit none

  private

  integer, parameter, public :: STORE_NONE  = 0
  integer, parameter, public :: STORE_PATCH = 1
  integer, parameter, public :: STORE_BDRY  = 2
  integer, parameter, public :: STORE_GHOST = 3

  integer, parameter, public :: NGB_INTERNAL = 1
  integer, parameter, public :: NGB_BLOCK    = 2
  integer, parameter, public :: NGB_DOMAIN   = 3
  integer, parameter, public :: NGB_ADAPT    = 4
  integer, parameter, public :: NGB_OTHER    = 5

  integer, parameter :: BLOCK_PACK_MAGIC = &
       int(z'54424C4B')
  integer, parameter :: BLOCK_PACK_VERSION = 5
  integer, parameter :: BLOCK_PACK_HEADER_SIZE = 36


  type, public :: Block_Bdry_Storage
     integer :: source_bdry = -1
     integer :: elts_start  = -1
     integer :: dims(2)     = 0
     integer :: n_node      = 0
     integer :: local_start = -1
  end type Block_Bdry_Storage


  type, public :: Block_Bdry_Link
     integer :: patch       = -1
     integer :: side        = -1
     integer :: class       = -1

     integer :: root_domain = -1

     integer :: neigh_patch = -1
     integer :: source_bdry = -1

     integer :: elts_start  = -1
     integer :: bdry_side   = 0
     integer :: bdry_neigh  = 0

     integer :: dims(2)     = 0
     integer :: n_node      = 0

     integer :: storage_id  = -1

     integer :: source_block    = -1
     integer :: source_block_id = -1
     integer :: source_owner    = -1
     integer :: ghost_id        = -1
  end type Block_Bdry_Link


  type, public :: Block_Ghost_Storage
     integer :: source_patch = -1
     integer :: elts_start   = -1
     integer :: local_start  = -1
     integer :: n_node       = 0

     integer :: source_domain   = -1
     integer :: source_block    = -1
     integer :: source_block_id = -1
     integer :: source_owner    = -1
  end type Block_Ghost_Storage


  type, public :: Block_Stencil_Address
     integer :: storage = STORE_NONE
     integer :: id      = -1
     integer :: offset  = 0
     integer :: dims(2) = 0
  end type Block_Stencil_Address


  type, public :: Block_Data
     integer :: id          = -1
     integer :: root_domain = -1
     integer :: root_patch  = -1
     integer :: level       = -1

     integer :: scalar_variable = -1
     integer :: n_scalar_variable = 0
     integer :: vector_variable = -1
     integer :: field_level     = -1
     integer :: n_field_level   = 0
     integer :: scalar_mult     = 0
     integer :: vector_mult     = 0

     type(Patch), allocatable :: patch(:)
     type(Coord), allocatable :: node(:)
     type(Coord), allocatable :: bdry_node(:)

     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
     real(dp), allocatable :: scalar_mean(:)
     real(dp), allocatable :: vector_mean(:)
     real(dp), allocatable :: bdry_scalar(:)
     real(dp), allocatable :: bdry_vector(:)
     real(dp), allocatable :: bdry_scalar_mean(:)
     real(dp), allocatable :: bdry_vector_mean(:)

     integer, allocatable :: neigh_class(:,:)

     type(Block_Bdry_Link), allocatable :: block_bdry(:)
     type(Block_Bdry_Storage), allocatable :: bdry_storage(:)
     type(Block_Stencil_Address), allocatable :: stencil(:,:)
     type(Block_Ghost_Storage), allocatable :: ghost_storage(:)

     type(Coord), allocatable :: ghost_node(:)

     real(dp), allocatable :: ghost_scalar(:)
     real(dp), allocatable :: ghost_vector(:)
     real(dp), allocatable :: ghost_scalar_mean(:)
     real(dp), allocatable :: ghost_vector_mean(:)
  end type Block_Data


  type(Block_Data), allocatable, public :: block_source(:)
  type(Block_Data), allocatable, public :: block_received(:)
  type(Block_Data), allocatable :: block_local(:)

  integer, allocatable, public :: block_source_catalog_index(:)
  integer, allocatable, public :: block_retained_source_index(:)
  integer, allocatable, public :: block_migrating_source_index(:)
  integer, allocatable, public :: block_received_catalog_index(:)
  integer, allocatable :: block_local_catalog_index(:)
  integer, allocatable :: block_catalog_local_index(:)

  logical :: block_store_ready = .false.

  public :: packed_block_nbyte
  public :: pack_block
  public :: unpack_block
  public :: check_block_storage
  public :: clear_block_staging
  public :: clear_local_blocks
  public :: local_block_store_ready
  public :: n_local_blocks
  public :: local_block_catalog
  public :: catalog_local_block
  public :: get_local_block_identity
  public :: get_block_field_layout
  public :: get_local_block_field_layout
  public :: check_local_block_storage
  public :: local_block_field_statistics
  public :: local_block_mean_field_statistics
  public :: install_local_blocks

contains


subroutine get_block_field_layout ( &
     scalar_variable,n_scalar_variable,vector_variable,field_level, &
     n_field_level,scalar_mult,vector_mult)
  ! Return the field layout represented by this Block_Data format.
  ! Centralizing the descriptor prevents builders and consumers from
  ! independently hard-coding the selected variables and sol levels.

  implicit none

  integer, intent(out) :: scalar_variable
  integer, intent(out) :: n_scalar_variable
  integer, intent(out) :: vector_variable
  integer, intent(out) :: field_level
  integer, intent(out) :: n_field_level
  integer, intent(out) :: scalar_mult
  integer, intent(out) :: vector_mult

  scalar_variable = scalars(1)
  n_scalar_variable = scalars(2) - scalars(1) + 1
  vector_variable = S_VELO
  field_level     = zmin
  n_field_level   = zmax - zmin + 1
  scalar_mult     = MULT(scalar_variable)
  vector_mult     = MULT(vector_variable)

end subroutine get_block_field_layout


subroutine check_block_storage (block,check_serialization)
  ! Check the allocation and extents of one block without referring to
  ! ownership, catalogues or MPI state. Optionally verify an exact
  ! pack/unpack/pack round trip.

  implicit none

  type(Block_Data), intent(in) :: block
  logical, optional, intent(in) :: check_serialization

  integer :: n_bdry_node
  integer :: n_ghost_node
  integer :: n_node
  integer :: scalar_variable
  integer :: n_scalar_variable
  integer :: vector_variable
  integer :: field_level
  integer :: n_field_level
  integer :: scalar_mult
  integer :: vector_mult

  integer(int8), allocatable :: buffer_copy(:)
  integer(int8), allocatable :: buffer_source(:)

  type(Block_Data) :: block_copy

  logical :: serialize

  serialize = .false.
  if (present(check_serialization)) serialize = check_serialization

  call get_block_field_layout( &
       scalar_variable,n_scalar_variable,vector_variable,field_level, &
       n_field_level,scalar_mult,vector_mult)

  if (block%scalar_variable /= scalar_variable .or. &
       block%n_scalar_variable /= n_scalar_variable .or. &
       block%vector_variable /= vector_variable .or. &
       block%field_level /= field_level .or. &
       block%n_field_level /= n_field_level .or. &
       block%scalar_mult /= scalar_mult .or. &
       block%vector_mult /= vector_mult) then

     error stop "check_block_storage: field layout mismatch"

  end if

  if (.not. allocated(block%patch) .or. &
       .not. allocated(block%node) .or. &
       .not. allocated(block%scalar) .or. &
       .not. allocated(block%vector) .or. &
       .not. allocated(block%scalar_mean) .or. &
       .not. allocated(block%vector_mean) .or. &
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%bdry_scalar_mean) .or. &
       .not. allocated(block%bdry_vector_mean) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector) .or. &
       .not. allocated(block%ghost_scalar_mean) .or. &
       .not. allocated(block%ghost_vector_mean)) then

     error stop "check_block_storage: unallocated component"

  end if

  n_node = size(block%patch) * PATCH_SIZE**2

  if (size(block%node) /= n_node .or. &
       size(block%scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_node .or. &
       size(block%vector) /= &
       block%n_field_level*block%vector_mult*n_node .or. &
       size(block%scalar_mean) /= size(block%scalar) .or. &
       size(block%vector_mean) /= size(block%vector)) then

     error stop "check_block_storage: interior extent mismatch"

  end if

  if (size(block%neigh_class,1) /= N_BDRY .or. &
       size(block%neigh_class,2) /= size(block%patch) .or. &
       size(block%stencil,1) /= N_BDRY .or. &
       size(block%stencil,2) /= size(block%patch)) then

     error stop "check_block_storage: topology extent mismatch"

  end if

  n_bdry_node = sum(block%bdry_storage%n_node)

  if (size(block%bdry_node) /= n_bdry_node .or. &
       size(block%bdry_scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_bdry_node .or. &
       size(block%bdry_vector) /= &
       block%n_field_level*block%vector_mult*n_bdry_node .or. &
       size(block%bdry_scalar_mean) /= size(block%bdry_scalar) .or. &
       size(block%bdry_vector_mean) /= size(block%bdry_vector)) then

     error stop "check_block_storage: boundary extent mismatch"

  end if

  n_ghost_node = sum(block%ghost_storage%n_node)

  if (size(block%ghost_node) /= n_ghost_node .or. &
       size(block%ghost_scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_ghost_node .or. &
       size(block%ghost_vector) /= &
       block%n_field_level*block%vector_mult*n_ghost_node .or. &
       size(block%ghost_scalar_mean) /= size(block%ghost_scalar) .or. &
       size(block%ghost_vector_mean) /= size(block%ghost_vector)) then

     error stop "check_block_storage: ghost extent mismatch"

  end if

  if (serialize) then

     call pack_block(block,buffer_source)
     call unpack_block(buffer_source,block_copy)
     call pack_block(block_copy,buffer_copy)

     if (size(buffer_copy) /= size(buffer_source)) then
        error stop &
             "check_block_storage: packed size mismatch"
     end if

     if (any(buffer_copy /= buffer_source)) then
        error stop &
             "check_block_storage: serialization mismatch"
     end if

  end if

end subroutine check_block_storage


subroutine install_local_blocks (n_catalog,local_seen)
  ! Build the final-owner local store from deep copies of the retained
  ! source blocks and validated received blocks. Catalogue ownership is
  ! deliberately checked by the caller, outside this non-MPI module.

  implicit none

  integer, intent(in) :: n_catalog
  integer, allocatable, intent(out) :: local_seen(:)

  integer :: b
  integer :: i
  integer :: ib
  integer :: ilocal
  integer :: n_local

  integer(int8), allocatable :: buffer_local(:)
  integer(int8), allocatable :: buffer_reference(:)

  if (n_catalog < 1) then
     error stop "install_local_blocks: invalid catalogue size"
  end if

  if (.not. allocated(block_source) .or. &
       .not. allocated(block_retained_source_index) .or. &
       .not. allocated(block_source_catalog_index)) then

     error stop "install_local_blocks: retained store missing"

  end if

  if (size(block_source_catalog_index) /= size(block_source)) then
     error stop "install_local_blocks: source map size mismatch"
  end if

  if (.not. allocated(block_received) .or. &
       .not. allocated(block_received_catalog_index)) then

     error stop "install_local_blocks: received store missing"

  end if

  if (size(block_received) /= &
       size(block_received_catalog_index)) then
     error stop "install_local_blocks: received map size mismatch"
  end if

  n_local = size(block_retained_source_index) + &
       size(block_received)

  call clear_local_blocks

  allocate(block_local(n_local))
  allocate(block_local_catalog_index(n_local))
  allocate(block_catalog_local_index(n_catalog))
  allocate(local_seen(n_catalog))

  block_local_catalog_index = -1
  block_catalog_local_index = 0
  local_seen = 0
  ilocal = 0

  do i = 1, size(block_retained_source_index)

     ib = block_retained_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop "install_local_blocks: invalid retained index"
     end if

     b = block_source_catalog_index(ib)

     if (b < 1 .or. b > n_catalog) then
        error stop &
             "install_local_blocks: invalid retained catalogue index"
     end if

     if (local_seen(b) /= 0) then
        error stop "install_local_blocks: duplicate retained block"
     end if

     ilocal = ilocal + 1
     block_local(ilocal) = block_source(ib)
     block_local_catalog_index(ilocal) = b
     block_catalog_local_index(b) = ilocal
     local_seen(b) = 1

     call pack_block(block_local(ilocal),buffer_local)
     call pack_block(block_source(ib),buffer_reference)

     if (size(buffer_local) /= size(buffer_reference)) then
        error stop &
             "install_local_blocks: retained copy size mismatch"
     end if

     if (any(buffer_local /= buffer_reference)) then
        error stop &
             "install_local_blocks: retained deep-copy mismatch"
     end if

  end do

  do i = 1, size(block_received)

     b = block_received_catalog_index(i)

     if (b < 1 .or. b > n_catalog) then
        error stop &
             "install_local_blocks: invalid received catalogue index"
     end if

     if (local_seen(b) /= 0) then
        error stop "install_local_blocks: duplicate received block"
     end if

     ilocal = ilocal + 1
     block_local(ilocal) = block_received(i)
     block_local_catalog_index(ilocal) = b
     block_catalog_local_index(b) = ilocal
     local_seen(b) = 1

     call pack_block(block_local(ilocal),buffer_local)
     call pack_block(block_received(i),buffer_reference)

     if (size(buffer_local) /= size(buffer_reference)) then
        error stop &
             "install_local_blocks: received copy size mismatch"
     end if

     if (any(buffer_local /= buffer_reference)) then
        error stop &
             "install_local_blocks: received deep-copy mismatch"
     end if

  end do

  if (ilocal /= n_local) then
     error stop "install_local_blocks: local fill count mismatch"
  end if

  if (count(local_seen /= 0) /= n_local) then
     error stop "install_local_blocks: local inventory mismatch"
  end if

  if (count(block_catalog_local_index > 0) /= n_local) then
     error stop "install_local_blocks: inverse map count mismatch"
  end if

  do ilocal = 1, n_local

     b = block_local_catalog_index(ilocal)

     if (block_catalog_local_index(b) /= ilocal) then
        error stop "install_local_blocks: inverse map mismatch"
     end if

  end do

  block_store_ready = .true.

end subroutine install_local_blocks


logical function local_block_store_ready () result(ready)
  ! Report whether a complete final-owner local store has been
  ! installed. Allocation checks make a partially modified store
  ! unavailable even if its readiness flag was set previously.

  implicit none

  ready = .false.

  if (.not. block_store_ready) return

  if (.not. allocated(block_local)) return
  if (.not. allocated(block_local_catalog_index)) return
  if (.not. allocated(block_catalog_local_index)) return

  if (size(block_local) /= &
       size(block_local_catalog_index)) return

  if (size(block_catalog_local_index) < 1) return

  ready = .true.

end function local_block_store_ready


integer function n_local_blocks () result(n_block)
  ! Return the number of blocks in the ready final-owner store.

  implicit none

  if (.not. local_block_store_ready()) then
     error stop "n_local_blocks: local store is not ready"
  end if

  n_block = size(block_local)

end function n_local_blocks


integer function local_block_catalog (local_index) result(catalog_index)
  ! Map a valid local block index to its replicated catalogue index.

  implicit none

  integer, intent(in) :: local_index

  if (.not. local_block_store_ready()) then
     error stop "local_block_catalog: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "local_block_catalog: invalid local index"
  end if

  catalog_index = block_local_catalog_index(local_index)

end function local_block_catalog


integer function catalog_local_block (catalog_index) result(local_index)
  ! Map a replicated catalogue index to its local block index. Return
  ! zero when the catalogue block belongs to another rank.

  implicit none

  integer, intent(in) :: catalog_index

  if (.not. local_block_store_ready()) then
     error stop "catalog_local_block: local store is not ready"
  end if

  if (catalog_index < 1 .or. &
       catalog_index > size(block_catalog_local_index)) then
     error stop "catalog_local_block: invalid catalogue index"
  end if

  local_index = block_catalog_local_index(catalog_index)

end function catalog_local_block


subroutine get_local_block_identity ( &
     local_index,id,root_domain,root_patch,level)
  ! Return immutable identity metadata for one local block.

  implicit none

  integer, intent(in)  :: local_index
  integer, intent(out) :: id
  integer, intent(out) :: root_domain
  integer, intent(out) :: root_patch
  integer, intent(out) :: level

  if (.not. local_block_store_ready()) then
     error stop "get_local_block_identity: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "get_local_block_identity: invalid local index"
  end if

  id          = block_local(local_index)%id
  root_domain = block_local(local_index)%root_domain
  root_patch  = block_local(local_index)%root_patch
  level       = block_local(local_index)%level

end subroutine get_local_block_identity


subroutine get_local_block_field_layout ( &
     local_index,scalar_variable,n_scalar_variable, &
     vector_variable,field_level, &
     n_field_level,scalar_mult,vector_mult)
  ! Return the serialized field descriptor for one installed block.

  implicit none

  integer, intent(in)  :: local_index
  integer, intent(out) :: scalar_variable
  integer, intent(out) :: n_scalar_variable
  integer, intent(out) :: vector_variable
  integer, intent(out) :: field_level
  integer, intent(out) :: n_field_level
  integer, intent(out) :: scalar_mult
  integer, intent(out) :: vector_mult

  if (.not. local_block_store_ready()) then
     error stop &
          "get_local_block_field_layout: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop &
          "get_local_block_field_layout: invalid local index"
  end if

  scalar_variable = block_local(local_index)%scalar_variable
  n_scalar_variable = block_local(local_index)%n_scalar_variable
  vector_variable = block_local(local_index)%vector_variable
  field_level     = block_local(local_index)%field_level
  n_field_level   = block_local(local_index)%n_field_level
  scalar_mult     = block_local(local_index)%scalar_mult
  vector_mult     = block_local(local_index)%vector_mult

end subroutine get_local_block_field_layout


subroutine check_local_block_storage (local_index,check_serialization)
  ! Validate one local block without exposing the private store.

  implicit none

  integer, intent(in) :: local_index
  logical, optional, intent(in) :: check_serialization

  logical :: serialize

  if (.not. local_block_store_ready()) then
     error stop "check_local_block_storage: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "check_local_block_storage: invalid local index"
  end if

  serialize = .false.
  if (present(check_serialization)) serialize = check_serialization

  call check_block_storage(block_local(local_index),serialize)

end subroutine check_local_block_storage


subroutine local_block_field_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Compute order-independent field inventory statistics over the
  ! ready final-owner block store. This is the first read-only field
  ! consumer of the private store; no block data are exposed or
  ! modified.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count

  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_field_statistics: local store is not ready"
  end if

  scalar_count  = 0_int64
  vector_count  = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%scalar) .or. &
          .not. allocated(block_local(i)%vector)) then
        error stop &
             "local_block_field_statistics: field storage missing"
     end if

     scalar_count = scalar_count + &
          int(size(block_local(i)%scalar),int64)
     vector_count = vector_count + &
          int(size(block_local(i)%vector),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_local(i)%scalar)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_local(i)%scalar))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_local(i)%scalar**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_local(i)%vector)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_local(i)%vector))
     vector_moment(3) = vector_moment(3) + &
          sum(block_local(i)%vector**2)

  end do

end subroutine local_block_field_statistics


subroutine local_block_mean_field_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Compute order-independent sol_mean inventory statistics over the
  ! ready final-owner block store.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count

  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_mean_field_statistics: store is not ready"
  end if

  scalar_count  = 0_int64
  vector_count  = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%scalar_mean) .or. &
          .not. allocated(block_local(i)%vector_mean)) then
        error stop &
             "local_block_mean_field_statistics: storage missing"
     end if

     scalar_count = scalar_count + &
          int(size(block_local(i)%scalar_mean),int64)
     vector_count = vector_count + &
          int(size(block_local(i)%vector_mean),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_local(i)%scalar_mean)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_local(i)%scalar_mean))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_local(i)%scalar_mean**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_local(i)%vector_mean)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_local(i)%vector_mean))
     vector_moment(3) = vector_moment(3) + &
          sum(block_local(i)%vector_mean**2)

  end do

end subroutine local_block_mean_field_statistics


subroutine clear_local_blocks
  ! Invalidate and release the persistent final-owner local store.
  ! This routine is deliberately idempotent so it is safe before the
  ! first restart and before every subsequent checkpoint restart.

  implicit none

  block_store_ready = .false.

  if (allocated(block_local)) then
     deallocate(block_local)
  end if

  if (allocated(block_local_catalog_index)) then
     deallocate(block_local_catalog_index)
  end if

  if (allocated(block_catalog_local_index)) then
     deallocate(block_catalog_local_index)
  end if

  if (allocated(block_local) .or. &
       allocated(block_local_catalog_index) .or. &
       allocated(block_catalog_local_index)) then
     error stop "clear_local_blocks: cleanup failed"
  end if

end subroutine clear_local_blocks


subroutine clear_block_staging
  ! Release only temporary source/receive migration stores. The
  ! final-owner block_local store and its catalogue are retained.

  implicit none

  if (allocated(block_source)) then
     deallocate(block_source)
  end if
  if (allocated(block_source_catalog_index)) then
     deallocate(block_source_catalog_index)
  end if
  if (allocated(block_retained_source_index)) then
     deallocate(block_retained_source_index)
  end if
  if (allocated(block_migrating_source_index)) then
     deallocate(block_migrating_source_index)
  end if
  if (allocated(block_received)) then
     deallocate(block_received)
  end if
  if (allocated(block_received_catalog_index)) then
     deallocate(block_received_catalog_index)
  end if

  if (allocated(block_source) .or. &
       allocated(block_source_catalog_index) .or. &
       allocated(block_retained_source_index) .or. &
       allocated(block_migrating_source_index) .or. &
       allocated(block_received) .or. &
       allocated(block_received_catalog_index)) then

     error stop "clear_block_staging: cleanup failed"

  end if

end subroutine clear_block_staging

integer function packed_block_nbyte (block) result(nbyte)
  ! Return the exact byte count used by pack_block.

  implicit none

  type(Block_Data), intent(in) :: block

  integer :: nbyte_integer

  if (.not. allocated(block%patch) .or. &
       .not. allocated(block%node) .or. &
       .not. allocated(block%scalar) .or. &
       .not. allocated(block%vector) .or. &
       .not. allocated(block%scalar_mean) .or. &
       .not. allocated(block%vector_mean) .or. &
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%bdry_scalar_mean) .or. &
       .not. allocated(block%bdry_vector_mean) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector) .or. &
       .not. allocated(block%ghost_scalar_mean) .or. &
       .not. allocated(block%ghost_vector_mean)) then

     error stop "packed_block_nbyte: unallocated component"

  end if

  nbyte_integer = storage_size(0) / 8

  nbyte = BLOCK_PACK_HEADER_SIZE * nbyte_integer

  nbyte = nbyte + &
       size(block%patch) * storage_size(block%patch) / 8

  nbyte = nbyte + &
       size(block%node) * storage_size(block%node) / 8

  nbyte = nbyte + &
       size(block%scalar) * storage_size(block%scalar) / 8

  nbyte = nbyte + &
       size(block%vector) * storage_size(block%vector) / 8

  nbyte = nbyte + &
       size(block%scalar_mean) * storage_size(block%scalar_mean) / 8

  nbyte = nbyte + &
       size(block%vector_mean) * storage_size(block%vector_mean) / 8

  nbyte = nbyte + size(block%neigh_class) * nbyte_integer

  nbyte = nbyte + &
       size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8

  nbyte = nbyte + &
       size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8

  nbyte = nbyte + &
       size(block%stencil) * storage_size(block%stencil) / 8

  nbyte = nbyte + &
       size(block%bdry_node) * storage_size(block%bdry_node) / 8

  nbyte = nbyte + &
       size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8

  nbyte = nbyte + &
       size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8

  nbyte = nbyte + &
       size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8

  nbyte = nbyte + &
       size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8

  nbyte = nbyte + &
       size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8

  nbyte = nbyte + &
       size(block%ghost_node) * storage_size(block%ghost_node) / 8

  nbyte = nbyte + &
       size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8

  nbyte = nbyte + &
       size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8

  nbyte = nbyte + &
       size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8

  nbyte = nbyte + &
       size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8

end function packed_block_nbyte


subroutine pack_block (block,buffer)
  ! Serialize a Block_Data into a versioned contiguous byte buffer.
  ! This raw representation is intended for homogeneous MPI jobs.

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int8), allocatable, intent(out) :: buffer(:)

  integer :: header(BLOCK_PACK_HEADER_SIZE)
  integer :: n
  integer :: nbyte
  integer :: pos

  nbyte = packed_block_nbyte(block)

  allocate(buffer(nbyte))

  header = [ &
       BLOCK_PACK_MAGIC, &
       BLOCK_PACK_VERSION, &
       block%id, &
       block%root_domain, &
       block%root_patch, &
       block%level, &
       block%scalar_variable, &
       block%n_scalar_variable, &
       block%vector_variable, &
       block%field_level, &
       block%n_field_level, &
       block%scalar_mult, &
       block%vector_mult, &
       size(block%patch), &
       size(block%node), &
       size(block%scalar), &
       size(block%vector), &
       size(block%scalar_mean), &
       size(block%vector_mean), &
       size(block%neigh_class,1), &
       size(block%neigh_class,2), &
       size(block%block_bdry), &
       size(block%bdry_storage), &
       size(block%stencil,1), &
       size(block%stencil,2), &
       size(block%bdry_node), &
       size(block%bdry_scalar), &
       size(block%bdry_vector), &
       size(block%bdry_scalar_mean), &
       size(block%bdry_vector_mean), &
       size(block%ghost_storage), &
       size(block%ghost_node), &
       size(block%ghost_scalar), &
       size(block%ghost_vector), &
       size(block%ghost_scalar_mean), &
       size(block%ghost_vector_mean) ]

  pos = 0

  n = size(header) * storage_size(header) / 8
  buffer(pos+1:pos+n) = transfer(header,0_int8,n)
  pos = pos + n

  n = size(block%patch) * storage_size(block%patch) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%patch,0_int8,n)
     pos = pos + n
  end if

  n = size(block%node) * storage_size(block%node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%scalar) * storage_size(block%scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%vector) * storage_size(block%vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%scalar_mean) * &
       storage_size(block%scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%vector_mean) * &
       storage_size(block%vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%vector_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%neigh_class) * storage_size(0) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%neigh_class,0_int8,n)
     pos = pos + n
  end if

  n = size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%block_bdry,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_storage,0_int8,n)
     pos = pos + n
  end if

  n = size(block%stencil) * storage_size(block%stencil) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%stencil,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_node) * storage_size(block%bdry_node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_vector_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_storage,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_node) * storage_size(block%ghost_node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_vector_mean,0_int8,n)
     pos = pos + n
  end if

  if (pos /= size(buffer)) then
     error stop "pack_block: final byte count mismatch"
  end if

end subroutine pack_block


subroutine unpack_block (buffer,block)
  ! Reconstruct a Block_Data from pack_block's byte stream.

  implicit none

  integer(int8), intent(in) :: buffer(:)
  type(Block_Data), intent(out) :: block

  integer :: header(BLOCK_PACK_HEADER_SIZE)
  integer :: n
  integer :: pos
  integer :: scalar_variable
  integer :: n_scalar_variable
  integer :: vector_variable
  integer :: field_level
  integer :: n_field_level
  integer :: scalar_mult
  integer :: vector_mult

  n = size(header) * storage_size(header) / 8

  if (size(buffer) < n) then
     error stop "unpack_block: truncated header"
  end if

  header = transfer(buffer(1:n),header,size(header))
  pos = n

  if (header(1) /= BLOCK_PACK_MAGIC) then
     error stop "unpack_block: invalid pack magic"
  end if

  if (header(2) /= BLOCK_PACK_VERSION) then
     error stop "unpack_block: unsupported pack version"
  end if

  if (any(header(14:BLOCK_PACK_HEADER_SIZE) < 0)) then
     error stop "unpack_block: negative component extent"
  end if

  call get_block_field_layout( &
       scalar_variable,n_scalar_variable,vector_variable,field_level, &
       n_field_level,scalar_mult,vector_mult)

  if (header(7) /= scalar_variable .or. &
       header(8) /= n_scalar_variable .or. &
       header(9) /= vector_variable .or. &
       header(10) /= field_level .or. &
       header(11) /= n_field_level .or. &
       header(12) /= scalar_mult .or. &
       header(13) /= vector_mult) then

     error stop "unpack_block: unsupported field layout"

  end if

  if (header(20) /= N_BDRY .or. &
       header(21) /= header(14) .or. &
       header(24) /= N_BDRY .or. &
       header(25) /= header(14)) then

     error stop "unpack_block: invalid topology extents"

  end if

  if (header(15) /= header(14)*PATCH_SIZE**2) then
     error stop "unpack_block: invalid interior node extent"
  end if

  block%id          = header(3)
  block%root_domain = header(4)
  block%root_patch  = header(5)
  block%level       = header(6)
  block%scalar_variable = header(7)
  block%n_scalar_variable = header(8)
  block%vector_variable = header(9)
  block%field_level     = header(10)
  block%n_field_level   = header(11)
  block%scalar_mult     = header(12)
  block%vector_mult     = header(13)

  allocate(block%patch(header(14)))
  allocate(block%node(header(15)))
  allocate(block%scalar(header(16)))
  allocate(block%vector(header(17)))
  allocate(block%scalar_mean(header(18)))
  allocate(block%vector_mean(header(19)))
  allocate(block%neigh_class(header(20),header(21)))
  allocate(block%block_bdry(header(22)))
  allocate(block%bdry_storage(header(23)))
  allocate(block%stencil(header(24),header(25)))
  allocate(block%bdry_node(header(26)))
  allocate(block%bdry_scalar(header(27)))
  allocate(block%bdry_vector(header(28)))
  allocate(block%bdry_scalar_mean(header(29)))
  allocate(block%bdry_vector_mean(header(30)))
  allocate(block%ghost_storage(header(31)))
  allocate(block%ghost_node(header(32)))
  allocate(block%ghost_scalar(header(33)))
  allocate(block%ghost_vector(header(34)))
  allocate(block%ghost_scalar_mean(header(35)))
  allocate(block%ghost_vector_mean(header(36)))

  n = size(block%patch) * storage_size(block%patch) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated patch data"
  end if
  if (n > 0) then
     block%patch = transfer( &
          buffer(pos+1:pos+n),block%patch,size(block%patch))
     pos = pos + n
  end if

  n = size(block%node) * storage_size(block%node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated node data"
  end if
  if (n > 0) then
     block%node = transfer( &
          buffer(pos+1:pos+n),block%node,size(block%node))
     pos = pos + n
  end if

  n = size(block%scalar) * storage_size(block%scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated scalar data"
  end if
  if (n > 0) then
     block%scalar = transfer( &
          buffer(pos+1:pos+n),block%scalar,size(block%scalar))
     pos = pos + n
  end if

  n = size(block%vector) * storage_size(block%vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated vector data"
  end if
  if (n > 0) then
     block%vector = transfer( &
          buffer(pos+1:pos+n),block%vector,size(block%vector))
     pos = pos + n
  end if

  n = size(block%scalar_mean) * &
       storage_size(block%scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean scalar data"
  end if
  if (n > 0) then
     block%scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%scalar_mean, &
          size(block%scalar_mean))
     pos = pos + n
  end if

  n = size(block%vector_mean) * &
       storage_size(block%vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean vector data"
  end if
  if (n > 0) then
     block%vector_mean = transfer( &
          buffer(pos+1:pos+n),block%vector_mean, &
          size(block%vector_mean))
     pos = pos + n
  end if

  n = size(block%neigh_class) * storage_size(0) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated neighbour data"
  end if
  if (n > 0) then
     block%neigh_class = reshape( &
          transfer(buffer(pos+1:pos+n),0, &
          size(block%neigh_class)),shape(block%neigh_class))
     pos = pos + n
  end if

  n = size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated block boundary data"
  end if
  if (n > 0) then
     block%block_bdry = transfer( &
          buffer(pos+1:pos+n),block%block_bdry, &
          size(block%block_bdry))
     pos = pos + n
  end if

  n = size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary catalogue"
  end if
  if (n > 0) then
     block%bdry_storage = transfer( &
          buffer(pos+1:pos+n),block%bdry_storage, &
          size(block%bdry_storage))
     pos = pos + n
  end if

  n = size(block%stencil) * storage_size(block%stencil) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated stencil data"
  end if
  if (n > 0) then
     block%stencil = reshape( &
          transfer(buffer(pos+1:pos+n),block%stencil(1,1), &
          size(block%stencil)),shape(block%stencil))
     pos = pos + n
  end if

  n = size(block%bdry_node) * storage_size(block%bdry_node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary nodes"
  end if
  if (n > 0) then
     block%bdry_node = transfer( &
          buffer(pos+1:pos+n),block%bdry_node, &
          size(block%bdry_node))
     pos = pos + n
  end if

  n = size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary scalar"
  end if
  if (n > 0) then
     block%bdry_scalar = transfer( &
          buffer(pos+1:pos+n),block%bdry_scalar, &
          size(block%bdry_scalar))
     pos = pos + n
  end if

  n = size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary vector"
  end if
  if (n > 0) then
     block%bdry_vector = transfer( &
          buffer(pos+1:pos+n),block%bdry_vector, &
          size(block%bdry_vector))
     pos = pos + n
  end if

  n = size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean boundary scalar"
  end if
  if (n > 0) then
     block%bdry_scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%bdry_scalar_mean, &
          size(block%bdry_scalar_mean))
     pos = pos + n
  end if

  n = size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean boundary vector"
  end if
  if (n > 0) then
     block%bdry_vector_mean = transfer( &
          buffer(pos+1:pos+n),block%bdry_vector_mean, &
          size(block%bdry_vector_mean))
     pos = pos + n
  end if

  n = size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost catalogue"
  end if
  if (n > 0) then
     block%ghost_storage = transfer( &
          buffer(pos+1:pos+n),block%ghost_storage, &
          size(block%ghost_storage))
     pos = pos + n
  end if

  n = size(block%ghost_node) * storage_size(block%ghost_node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost nodes"
  end if
  if (n > 0) then
     block%ghost_node = transfer( &
          buffer(pos+1:pos+n),block%ghost_node, &
          size(block%ghost_node))
     pos = pos + n
  end if

  n = size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost scalar"
  end if
  if (n > 0) then
     block%ghost_scalar = transfer( &
          buffer(pos+1:pos+n),block%ghost_scalar, &
          size(block%ghost_scalar))
     pos = pos + n
  end if

  n = size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost vector"
  end if
  if (n > 0) then
     block%ghost_vector = transfer( &
          buffer(pos+1:pos+n),block%ghost_vector, &
          size(block%ghost_vector))
     pos = pos + n
  end if

  n = size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean ghost scalar"
  end if
  if (n > 0) then
     block%ghost_scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%ghost_scalar_mean, &
          size(block%ghost_scalar_mean))
     pos = pos + n
  end if

  n = size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean ghost vector"
  end if
  if (n > 0) then
     block%ghost_vector_mean = transfer( &
          buffer(pos+1:pos+n),block%ghost_vector_mean, &
          size(block%ghost_vector_mean))
     pos = pos + n
  end if

  if (pos /= size(buffer)) then
     error stop "unpack_block: final byte count mismatch"
  end if

end subroutine unpack_block

end module parallel_block_mod
