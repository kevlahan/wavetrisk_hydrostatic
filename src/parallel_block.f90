module parallel_block_mod

  use, intrinsic :: iso_fortran_env, only : int8

  use kind_mod,   only : dp
  use shared_mod, only : Coord, MULT, N_BDRY, S_VELO, scalars
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
  integer, parameter :: BLOCK_PACK_VERSION = 1
  integer, parameter :: BLOCK_PACK_HEADER_SIZE = 23


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

     type(Patch), allocatable :: patch(:)
     type(Coord), allocatable :: node(:)
     type(Coord), allocatable :: bdry_node(:)

     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
     real(dp), allocatable :: bdry_scalar(:)
     real(dp), allocatable :: bdry_vector(:)

     integer, allocatable :: neigh_class(:,:)

     type(Block_Bdry_Link), allocatable :: block_bdry(:)
     type(Block_Bdry_Storage), allocatable :: bdry_storage(:)
     type(Block_Stencil_Address), allocatable :: stencil(:,:)
     type(Block_Ghost_Storage), allocatable :: ghost_storage(:)

     type(Coord), allocatable :: ghost_node(:)

     real(dp), allocatable :: ghost_scalar(:)
     real(dp), allocatable :: ghost_vector(:)
  end type Block_Data


  type(Block_Data), allocatable, public :: block_source(:)
  type(Block_Data), allocatable, public :: block_received(:)
  type(Block_Data), allocatable, public :: block_local(:)

  integer, allocatable, public :: block_source_catalog_index(:)
  integer, allocatable, public :: block_retained_source_index(:)
  integer, allocatable, public :: block_migrating_source_index(:)
  integer, allocatable, public :: block_received_catalog_index(:)
  integer, allocatable, public :: block_local_catalog_index(:)

  logical :: block_store_ready = .false.

  public :: packed_block_nbyte
  public :: pack_block
  public :: unpack_block
  public :: check_block_storage
  public :: clear_block_staging
  public :: clear_local_blocks
  public :: local_block_store_ready
  public :: install_local_blocks

contains


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

  integer(int8), allocatable :: buffer_copy(:)
  integer(int8), allocatable :: buffer_source(:)

  type(Block_Data) :: block_copy

  logical :: serialize

  serialize = .false.
  if (present(check_serialization)) serialize = check_serialization

  if (.not. allocated(block%patch) .or. &
       .not. allocated(block%node) .or. &
       .not. allocated(block%scalar) .or. &
       .not. allocated(block%vector) .or. &
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector)) then

     error stop "check_block_storage: unallocated component"

  end if

  n_node = size(block%patch) * PATCH_SIZE**2

  if (size(block%node) /= n_node .or. &
       size(block%scalar) /= MULT(scalars(1))*n_node .or. &
       size(block%vector) /= MULT(S_VELO)*n_node) then

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
       size(block%bdry_scalar) /= MULT(scalars(1))*n_bdry_node .or. &
       size(block%bdry_vector) /= MULT(S_VELO)*n_bdry_node) then

     error stop "check_block_storage: boundary extent mismatch"

  end if

  n_ghost_node = sum(block%ghost_storage%n_node)

  if (size(block%ghost_node) /= n_ghost_node .or. &
       size(block%ghost_scalar) /= MULT(scalars(1))*n_ghost_node .or. &
       size(block%ghost_vector) /= MULT(S_VELO)*n_ghost_node) then

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
  allocate(local_seen(n_catalog))

  block_local_catalog_index = -1
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

  block_store_ready = .true.

end subroutine install_local_blocks


logical function local_block_store_ready () result(ready)
  ! Report whether a complete final-owner local store has been
  ! installed. Allocation checks make a partially modified store
  ! unavailable even if its readiness flag was set previously.

  implicit none

  ready = block_store_ready .and. &
       allocated(block_local) .and. &
       allocated(block_local_catalog_index)

  if (ready) then
     ready = size(block_local) == &
          size(block_local_catalog_index)
  end if

end function local_block_store_ready


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

  if (allocated(block_local) .or. &
       allocated(block_local_catalog_index)) then
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
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector)) then

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
       size(block%patch), &
       size(block%node), &
       size(block%scalar), &
       size(block%vector), &
       size(block%neigh_class,1), &
       size(block%neigh_class,2), &
       size(block%block_bdry), &
       size(block%bdry_storage), &
       size(block%stencil,1), &
       size(block%stencil,2), &
       size(block%bdry_node), &
       size(block%bdry_scalar), &
       size(block%bdry_vector), &
       size(block%ghost_storage), &
       size(block%ghost_node), &
       size(block%ghost_scalar), &
       size(block%ghost_vector) ]

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

  if (any(header(7:BLOCK_PACK_HEADER_SIZE) < 0)) then
     error stop "unpack_block: negative component extent"
  end if

  if (header(11) /= N_BDRY .or. &
       header(12) /= header(7) .or. &
       header(15) /= N_BDRY .or. &
       header(16) /= header(7)) then

     error stop "unpack_block: invalid topology extents"

  end if

  if (header(8) /= header(7)*PATCH_SIZE**2) then
     error stop "unpack_block: invalid interior node extent"
  end if

  block%id          = header(3)
  block%root_domain = header(4)
  block%root_patch  = header(5)
  block%level       = header(6)

  allocate(block%patch(header(7)))
  allocate(block%node(header(8)))
  allocate(block%scalar(header(9)))
  allocate(block%vector(header(10)))
  allocate(block%neigh_class(header(11),header(12)))
  allocate(block%block_bdry(header(13)))
  allocate(block%bdry_storage(header(14)))
  allocate(block%stencil(header(15),header(16)))
  allocate(block%bdry_node(header(17)))
  allocate(block%bdry_scalar(header(18)))
  allocate(block%bdry_vector(header(19)))
  allocate(block%ghost_storage(header(20)))
  allocate(block%ghost_node(header(21)))
  allocate(block%ghost_scalar(header(22)))
  allocate(block%ghost_vector(header(23)))

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

  if (pos /= size(buffer)) then
     error stop "unpack_block: final byte count mismatch"
  end if

end subroutine unpack_block

end module parallel_block_mod
