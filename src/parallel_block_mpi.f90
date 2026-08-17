module parallel_block_mpi_mod

  use iso_fortran_env, only : error_unit, int8, int64
  use mpi_f08,        only : MPI_Allgather, MPI_Allgatherv, MPI_Allreduce, &
       MPI_Alltoall, MPI_Alltoallv, MPI_Exscan, MPI_BYTE, MPI_INTEGER, &
       MPI_INTEGER8, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_SUCCESS, MPI_SUM

  use kind_mod,   only : dp
  use shared_mod, only : N_CHDRN, N_GLO_DOMAIN

  use domain_mod, only : grid, sol, sol_mean, wav_coeff, &
       subtree_weight_Domain

  use patch_mod, only : PATCH_SIZE

  use arch_mod, only : abort_run, block_catalog, comm, glo_id, loc_id, &
       n_process, owner, Parallel_Block, rank

  use parallel_block_mod, only : block_source, block_received, &
       block_source_catalog_index, &
       block_migrating_source_index, block_received_catalog_index, &
       packed_block_nbyte, pack_block, unpack_block, &
       check_block_storage, install_local_blocks, clear_block_staging, &
       clear_local_blocks, local_block_store_ready, n_local_blocks, &
       local_block_catalog, catalog_local_block, &
       get_local_block_identity, check_local_block_storage, &
       get_block_field_layout, get_local_block_field_layout, &
       local_block_field_statistics, local_block_wavelet_statistics, &
       local_block_mean_field_statistics

  implicit none

  private

  type, public :: Block_Migration_Manifest
     integer :: n_send = 0
     integer :: n_recv = 0
     integer, allocatable :: send_count(:)
     integer, allocatable :: recv_count(:)
     integer, allocatable :: send_displ(:)
     integer, allocatable :: recv_displ(:)
     integer, allocatable :: send_block(:)
     integer, allocatable :: recv_block(:)
     integer, allocatable :: send_nbyte(:)
     integer, allocatable :: recv_nbyte(:)
     integer, allocatable :: send_byte_count(:)
     integer, allocatable :: recv_byte_count(:)
     integer, allocatable :: send_byte_displ(:)
     integer, allocatable :: recv_byte_displ(:)
     integer(int8), allocatable :: recv_payload(:)
     integer(int64) :: total_send_nbyte = 0_int64
     integer(int64) :: total_recv_nbyte = 0_int64
     integer :: max_send_nbyte = 0
     integer :: max_recv_nbyte = 0
     logical :: validated = .false.
     logical :: sizes_validated = .false.
     logical :: payload_validated = .false.
  end type Block_Migration_Manifest

  public :: build_block_migration_manifest
  public :: check_block_migration_manifest
  public :: exchange_block_migration_sizes
  public :: exchange_block_migration_payloads
  public :: clear_block_migration_manifest
  public :: build_parallel_block_catalog
  public :: clear_parallel_block_state
  public :: migrate_blocks
  public :: check_local_blocks
  public :: check_block_field_inventory

contains


  subroutine clear_parallel_block_state
    ! Release persistent and staging block data before a checkpoint
    ! restart invalidates the geometric-domain state from which those
    ! blocks were constructed. Safe to call when nothing is allocated.

    implicit none

    call clear_local_blocks
    call clear_block_staging

    if (allocated(block_catalog)) deallocate(block_catalog)

    if (local_block_store_ready()) then
       call fail("local block store remained ready after reset")
    end if

    if (allocated(block_catalog)) then
       call fail("block catalogue remained allocated after reset")
    end if

  end subroutine clear_parallel_block_state


  subroutine build_parallel_block_catalog
  ! Construct the parallel-block decomposition of all locally owned
  ! geometric root domains, build the replicated global block catalogue,
  ! and assign the target block owners using load balancing.
  !
  ! Rule:
  !   - if root patch 1 has all four children, use those four child
  !     subtrees as candidate blocks;
  !   - otherwise keep the whole root domain as one candidate block.
  !
  ! This routine does not modify the active geometric-domain storage.

  implicit none

  integer, parameter :: N_BLOCK_DATA = 6

  integer :: b, c, d, i, k, p, p_chd, r
  integer :: n_block_local, n_block_before, n_block_total
  integer :: n_data_local, n_data_total
  integer :: n_assigned, n_changed
  integer :: ierr

  integer, allocatable :: block_count(:)
  integer, allocatable :: recv_count(:)
  integer, allocatable :: recv_disp(:)
  integer, allocatable :: send_data(:)
  integer, allocatable :: recv_data(:)
  integer, allocatable :: domain_count(:)

  integer, allocatable :: block_owner_new(:)
  integer, allocatable :: load_current(:)
  integer, allocatable :: load_proposed(:)

  real(dp) :: balanced_weight
  real(dp) :: imbalance_goal

  real(dp), parameter :: init_goal = 0.05_dp
  real(dp), parameter :: incr_goal = 1.20_dp

  type(Parallel_Block), allocatable :: block_loc(:)

  if (local_block_store_ready()) then
     call fail("existing local block store was not reset")
  end if

  !
  ! Count local candidate blocks.
  !
  n_block_local = 0

  do d = 1, size(grid)

     p = 1

     if (all(grid(d)%patch%elts(p+1)%children > 0)) then
        n_block_local = n_block_local + N_CHDRN
     else
        n_block_local = n_block_local + 1
     end if

  end do

  allocate(block_loc(n_block_local))

  !
  ! Construct local candidate block descriptors.
  !
  i = 0

  do d = 1, size(grid)

     p = 1

     if (all(grid(d)%patch%elts(p+1)%children > 0)) then

        do c = 1, N_CHDRN

           p_chd = grid(d)%patch%elts(p+1)%children(c)

           i = i + 1

           block_loc(i)%root_domain = glo_id(rank+1,d)
           block_loc(i)%root_patch  = p_chd
           block_loc(i)%level       = grid(d)%patch%elts(p_chd+1)%level
           block_loc(i)%owner       = rank
           block_loc(i)%weight      = &
                subtree_weight_Domain(grid(d), p_chd)

        end do

     else

        i = i + 1

        block_loc(i)%root_domain = glo_id(rank+1,d)
        block_loc(i)%root_patch  = p
        block_loc(i)%level       = grid(d)%patch%elts(p+1)%level
        block_loc(i)%owner       = rank
        block_loc(i)%weight      = &
             subtree_weight_Domain(grid(d), p)

     end if

  end do

  if (i /= n_block_local) then
     error stop &
          "build_parallel_block_catalog: local block count mismatch"
  end if

  !
  ! Determine the number of candidate blocks preceding those on
  ! this rank. This gives globally unique contiguous block IDs.
  !
  n_block_before = 0

  call MPI_Exscan( &
       n_block_local, n_block_before, &
       1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "build_parallel_block_catalog: MPI_Exscan failed"
  end if

  if (rank == 0) n_block_before = 0

  !
  ! Determine the total number of candidate blocks.
  !
  call MPI_Allreduce( &
       n_block_local, n_block_total, &
       1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "build_parallel_block_catalog: MPI_Allreduce failed"
  end if

  !
  ! Assign global candidate-block IDs.
  !
  do i = 1, n_block_local
     block_loc(i)%id = n_block_before + i - 1
  end do

  !
  ! Determine how many candidate blocks are contributed by every rank.
  !
  allocate(block_count(n_process))

  call MPI_Allgather( &
       n_block_local, 1, MPI_INTEGER, &
       block_count, 1, MPI_INTEGER, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "build_parallel_block_catalog: MPI_Allgather failed"
  end if

  if (sum(block_count) /= n_block_total) then
     error stop &
          "build_parallel_block_catalog: inconsistent global block count"
  end if

  !
  ! Pack each local candidate block as six integers:
  !
  !   id, root_domain, root_patch, level, owner, weight
  !
  n_data_local = N_BLOCK_DATA * n_block_local

  allocate(send_data(n_data_local))

  k = 0

  do i = 1, n_block_local

     send_data(k+1:k+N_BLOCK_DATA) = [ &
          block_loc(i)%id, &
          block_loc(i)%root_domain, &
          block_loc(i)%root_patch, &
          block_loc(i)%level, &
          block_loc(i)%owner, &
          block_loc(i)%weight ]

     k = k + N_BLOCK_DATA

  end do

  if (k /= n_data_local) then
     error stop &
          "build_parallel_block_catalog: local packing count mismatch"
  end if

  !
  ! Construct Allgatherv receive counts and displacements.
  !
  allocate(recv_count(n_process))
  allocate(recv_disp(n_process))

  recv_count = N_BLOCK_DATA * block_count

  recv_disp(1) = 0

  do i = 2, n_process
     recv_disp(i) = recv_disp(i-1) + recv_count(i-1)
  end do

  n_data_total = N_BLOCK_DATA * n_block_total

  allocate(recv_data(n_data_total))

  !
  ! Every rank receives the complete packed candidate-block catalogue.
  !
  call MPI_Allgatherv( &
       send_data, n_data_local, MPI_INTEGER, &
       recv_data, recv_count, recv_disp, MPI_INTEGER, &
       comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "build_parallel_block_catalog: MPI_Allgatherv failed"
  end if

  !
  ! Construct the global Parallel_Block catalogue on every rank.
  !
  if (allocated(block_catalog)) deallocate(block_catalog)

  allocate(block_catalog(n_block_total))

  do i = 1, n_block_total

     k = N_BLOCK_DATA * (i-1)

     block_catalog(i)%id          = recv_data(k+1)
     block_catalog(i)%root_domain = recv_data(k+2)
     block_catalog(i)%root_patch  = recv_data(k+3)
     block_catalog(i)%level       = recv_data(k+4)
     block_catalog(i)%owner       = recv_data(k+5)
     block_catalog(i)%weight      = recv_data(k+6)

  end do

  !
  ! Validate global block IDs and owners.
  !
  do i = 1, n_block_total

     if (block_catalog(i)%id /= i-1) then
        error stop &
             "build_parallel_block_catalog: invalid global block ID ordering"
     end if

     if (block_catalog(i)%owner < 0 .or. &
          block_catalog(i)%owner >= n_process) then
        error stop &
             "build_parallel_block_catalog: invalid block owner"
     end if

  end do

  !
  ! Verify that every geometric root domain occurs in the catalogue.
  ! With the current one-level split rule, each domain should occur
  ! either once or four times.
  !
  allocate(domain_count(N_GLO_DOMAIN))
  domain_count = 0

  do i = 1, n_block_total

     if (block_catalog(i)%root_domain < 0 .or. &
          block_catalog(i)%root_domain >= N_GLO_DOMAIN) then
        error stop &
             "build_parallel_block_catalog: invalid root domain"
     end if

     domain_count(block_catalog(i)%root_domain+1) = &
          domain_count(block_catalog(i)%root_domain+1) + 1

  end do

  if (any(domain_count == 0)) then
     error stop &
          "build_parallel_block_catalog: one or more root domains missing"
  end if

  if (any(domain_count /= 1 .and. domain_count /= N_CHDRN)) then
     error stop &
          "build_parallel_block_catalog: invalid number of blocks per root domain"
  end if

  !
  ! ---------------------------------------------------------------
  ! Compute the target load balancing of candidate parallel blocks.
  ! ---------------------------------------------------------------
  !
  ! This does not change actual ownership.
  !
  allocate(block_owner_new(n_block_total))
  allocate(load_current(n_process))
  allocate(load_proposed(n_process))

  block_owner_new = -1
  load_current     = 0
  load_proposed    = 0

  !
  ! Current load implied by the existing domain decomposition.
  !
  do b = 1, n_block_total

     r = block_catalog(b)%owner + 1

     if (r < 1 .or. r > n_process) then
        error stop &
             "build_parallel_block_catalog: invalid current block owner"
     end if

     load_current(r) = load_current(r) + block_catalog(b)%weight

  end do

  !
  ! Ideal average load per MPI rank.
  !
  balanced_weight = &
       real(sum(block_catalog%weight), kind=dp) / &
       real(n_process, kind=dp)

  !
  ! Use the same next-fit strategy as distribute_grid(), but apply it
  ! only to the prospective candidate blocks.
  !
  n_assigned     = 0
  imbalance_goal = init_goal

  do while (n_assigned < n_block_total)

     n_assigned      = 0
     load_proposed   = 0
     block_owner_new = -1

     do r = 1, n_process

        do while ( &
             real(load_proposed(r),kind=dp) < balanced_weight .and. &
             n_block_total - n_assigned > n_process - r)

           block_owner_new(n_assigned+1) = r - 1

           load_proposed(r) = load_proposed(r) + &
                block_catalog(n_assigned+1)%weight

           n_assigned = n_assigned + 1

        end do

        !
        ! If the final block makes this rank too heavy,
        ! move that block to the next rank.
        !
        if (n_assigned > 0) then

           if (real(load_proposed(r),kind=dp) > &
                balanced_weight * (1.0_dp + imbalance_goal)) then

              load_proposed(r) = load_proposed(r) - &
                   block_catalog(n_assigned)%weight

              block_owner_new(n_assigned) = -1

              n_assigned = n_assigned - 1

           end if

        end if

     end do

     !
     ! If all blocks could not be assigned, relax the allowed imbalance.
     !
     if (n_assigned < n_block_total) then
        imbalance_goal = imbalance_goal * incr_goal
     end if

  end do

  if (any(block_owner_new < 0)) then
     error stop &
          "build_parallel_block_catalog: one or more candidate blocks unassigned"
  end if

  !
  ! Independently reconstruct proposed rank loads from the owner map.
  !
  load_proposed = 0

  do b = 1, n_block_total

     r = block_owner_new(b) + 1

     if (r < 1 .or. r > n_process) then
        error stop &
             "build_parallel_block_catalog: invalid proposed block owner"
     end if

     load_proposed(r) = load_proposed(r) + block_catalog(b)%weight

  end do

  n_changed = count(block_owner_new /= block_catalog%owner)

  !
  ! Commit the prospective balanced ownership to the global block
  ! catalogue. The current source-domain owner remains available through
  ! owner(root_domain+1).
  !
  block_catalog%owner = block_owner_new

  !
  ! Validate the committed target owners.
  !
  if (any(block_catalog%owner < 0) .or. &
       any(block_catalog%owner >= n_process)) then

     error stop &
          "build_parallel_block_catalog: invalid committed block owner"

  end if
  
  !
  ! Print compact diagnostic summary.
  !
  if (rank == 0) then

     write(6,'(/,a,i0)') &
          "Total candidate parallel blocks = ", n_block_total

     write(6,'(a,i0)') &
          "Unsplit root domains = ", count(domain_count == 1)

     write(6,'(a,i0)') &
          "Split root domains   = ", count(domain_count == N_CHDRN)

     write(6,'(/,a)') &
          "Current candidate-block load:"

     write(6,'(a,i0)') &
          "  minimum = ", minval(load_current)

     write(6,'(a,f10.2)') &
          "  average = ", balanced_weight

     write(6,'(a,i0)') &
          "  maximum = ", maxval(load_current)

     write(6,'(a,f10.4)') &
          "  max/avg = ", &
          real(maxval(load_current),kind=dp) / balanced_weight

     write(6,'(/,a)') &
          "Prospective balanced block load:"

     write(6,'(a,i0)') &
          "  minimum = ", minval(load_proposed)

     write(6,'(a,f10.2)') &
          "  average = ", balanced_weight

     write(6,'(a,i0)') &
          "  maximum = ", maxval(load_proposed)

     write(6,'(a,f10.4)') &
          "  max/avg = ", &
          real(maxval(load_proposed),kind=dp) / balanced_weight

     write(6,'(a,i0,a,i0,/)') &
          "Blocks changing owner = ", n_changed, &
          " / ", n_block_total

  end if

  !
  ! Cleanup temporary arrays.
  !
  deallocate(block_loc)
  deallocate(block_count)
  deallocate(recv_count)
  deallocate(recv_disp)
  deallocate(send_data)
  deallocate(recv_data)
  deallocate(domain_count)

  deallocate(block_owner_new)
  deallocate(load_current)
  deallocate(load_proposed)

end subroutine build_parallel_block_catalog


  subroutine migrate_blocks (verbose)
    ! Migrate the source-local blocks to their final catalogue owners,
    ! install a self-contained local block store, and release all
    ! source, receive and MPI migration staging storage.

    implicit none

    logical, optional, intent(in) :: verbose

    type(Block_Migration_Manifest) :: manifest

    integer :: b
    integer :: expected_local
    integer :: i
    integer :: ib
    integer :: n_recv_byte
    integer :: n_received
    integer :: n_send_byte
    integer :: n_sent
    integer :: nbyte
    integer :: pos

    integer, allocatable :: local_seen(:)
    integer, allocatable :: seen_catalog(:)
    integer, allocatable :: seen_source(:)
    integer, allocatable :: send_nbyte(:)

    integer(int8), allocatable :: block_buffer(:)
    integer(int8), allocatable :: send_payload(:)

    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. allocated(block_source) .or. &
         .not. allocated(block_source_catalog_index) .or. &
         .not. allocated(block_migrating_source_index)) then
       call fail("source block store is incomplete")
    end if

    if (size(block_source_catalog_index) /= size(block_source)) then
       call fail("source block catalogue map has the wrong extent")
    end if

    call build_block_migration_manifest(manifest)
    call check_block_migration_manifest(manifest,print_local)

    if (manifest%n_send /= size(block_migrating_source_index)) then
       call fail("manifest and source migration counts differ")
    end if

    allocate(send_nbyte(manifest%n_send))
    allocate(seen_source(size(block_source)))

    send_nbyte = 0
    seen_source = 0

    do i = 1, manifest%n_send

       b = manifest%send_block(i)
       ib = findloc(block_source_catalog_index,b,dim=1)

       if (ib < 1 .or. ib > size(block_source)) then
          call fail("manifest source block was not constructed locally")
       end if

       if (findloc(block_migrating_source_index,ib,dim=1) < 1) then
          call fail("manifest source block is not marked for migration")
       end if

       if (seen_source(ib) /= 0) then
          call fail("source block occurs more than once in manifest")
       end if

       send_nbyte(i) = packed_block_nbyte(block_source(ib))

       if (send_nbyte(i) <= 0) then
          call fail("source block has an invalid packed byte count")
       end if

       seen_source(ib) = 1

    end do

    do i = 1, size(block_migrating_source_index)

       ib = block_migrating_source_index(i)

       if (ib < 1 .or. ib > size(block_source)) then
          call fail("migrating source index is invalid")
       end if

       if (seen_source(ib) /= 1) then
          call fail("migrating source block is absent from manifest")
       end if

    end do

    call exchange_block_migration_sizes( &
         manifest,send_nbyte,print_local)

    if (manifest%total_send_nbyte > int(huge(0),int64) .or. &
         manifest%total_recv_nbyte > int(huge(0),int64)) then
       call fail("packed migration size exceeds default integer range")
    end if

    n_send_byte = int(manifest%total_send_nbyte)
    n_recv_byte = int(manifest%total_recv_nbyte)

    allocate(send_payload(max(1,n_send_byte)))
    send_payload = 0_int8
    pos = 0

    do i = 1, manifest%n_send

       b = manifest%send_block(i)
       ib = findloc(block_source_catalog_index,b,dim=1)

       if (ib < 1 .or. ib > size(block_source)) then
          call fail("source block disappeared before packing")
       end if

       call pack_block(block_source(ib),block_buffer)

       if (size(block_buffer) /= manifest%send_nbyte(i)) then
          call fail("source block packed size changed")
       end if

       nbyte = size(block_buffer)

       if (pos+nbyte > n_send_byte) then
          call fail("outgoing packed payload exceeds its buffer")
       end if

       send_payload(pos+1:pos+nbyte) = block_buffer
       pos = pos + nbyte

    end do

    if (pos /= n_send_byte) then
       call fail("outgoing packed payload has the wrong extent")
    end if

    call exchange_block_migration_payloads( &
         manifest,send_payload,print_local)

    if (allocated(block_received)) deallocate(block_received)
    if (allocated(block_received_catalog_index)) then
       deallocate(block_received_catalog_index)
    end if

    allocate(block_received(manifest%n_recv))
    allocate(block_received_catalog_index(manifest%n_recv))
    allocate(seen_catalog(size(block_catalog)))

    block_received_catalog_index = -1
    seen_catalog = 0
    pos = 0

    do i = 1, manifest%n_recv

       b = manifest%recv_block(i)
       nbyte = manifest%recv_nbyte(i)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("received block catalogue index is invalid")
       end if

       if (seen_catalog(b) /= 0) then
          call fail("received block occurs more than once")
       end if

       if (block_catalog(b)%owner /= rank) then
          call fail("received block has the wrong final owner")
       end if

       if (nbyte <= 0 .or. pos+nbyte > n_recv_byte) then
          call fail("received packed block extent is invalid")
       end if

       call unpack_block( &
            manifest%recv_payload(pos+1:pos+nbyte), &
            block_received(i))

       if (block_received(i)%id /= block_catalog(b)%id .or. &
            block_received(i)%root_domain /= &
            block_catalog(b)%root_domain .or. &
            block_received(i)%root_patch /= &
            block_catalog(b)%root_patch .or. &
            block_received(i)%level /= block_catalog(b)%level) then
          call fail("received block identity does not match catalogue")
       end if

       call check_block_storage(block_received(i),.true.)

       block_received_catalog_index(i) = b
       seen_catalog(b) = 1
       pos = pos + nbyte

    end do

    if (pos /= n_recv_byte) then
       call fail("received packed payload has the wrong extent")
    end if

    if (count(seen_catalog /= 0) /= manifest%n_recv) then
       call fail("received block inventory is incomplete")
    end if

    call install_local_blocks(size(block_catalog),local_seen)

    expected_local = 0

    do b = 1, size(block_catalog)
       if (block_catalog(b)%owner == rank) then
          expected_local = expected_local + 1
          if (local_seen(b) /= 1) then
             call fail("owned block is absent from installed store")
          end if
       else
          if (local_seen(b) /= 0) then
             call fail("nonlocal block appears in installed store")
          end if
       end if
    end do

    if (n_local_blocks() /= expected_local) then
       call fail("installed block count does not match final ownership")
    end if

    n_sent     = manifest%n_send
    n_received = manifest%n_recv

    call clear_block_staging
    call clear_block_migration_manifest(manifest)

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Completed block migration for rank ", rank, ":"
       write(6,'(a,i0)') "  blocks sent      = ", n_sent
       write(6,'(a,i0)') "  blocks received  = ", n_received
       write(6,'(a,i0)') &
            "  blocks installed = ", n_local_blocks()
       write(6,'(a,/)') &
            "  migration staging storage released"
    end if

    deallocate(local_seen)
    deallocate(seen_catalog)
    deallocate(seen_source)
    deallocate(send_nbyte)
    deallocate(send_payload)

    call check_local_blocks(print_local)
    call check_block_field_inventory(print_local)

  end subroutine migrate_blocks


  subroutine check_local_blocks (verbose)
    ! Validate the final-owner local block store without referring to
    ! source, receive or migration-manifest staging allocation.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: expected_local
    integer :: expected_field_level
    integer :: expected_n_field_level
    integer :: expected_n_scalar_variable
    integer :: expected_scalar_mult
    integer :: expected_scalar_variable
    integer :: expected_vector_mult
    integer :: expected_vector_variable
    integer :: field_level
    integer :: global_count
    integer :: global_weight
    integer :: i
    integer :: id
    integer :: ierr
    integer :: level
    integer :: local_count
    integer :: local_weight
    integer :: root_domain
    integer :: root_patch
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: scalar_mult
    integer :: scalar_variable
    integer :: vector_mult
    integer :: vector_variable

    integer, allocatable :: global_seen(:)
    integer, allocatable :: local_seen(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("local block store is not ready")
    end if

    allocate(local_seen(size(block_catalog)))
    allocate(global_seen(size(block_catalog)))

    local_seen   = 0
    global_seen  = 0
    local_count  = n_local_blocks()
    local_weight = 0

    call get_block_field_layout( &
         expected_scalar_variable,expected_n_scalar_variable, &
         expected_vector_variable, &
         expected_field_level,expected_n_field_level, &
         expected_scalar_mult, &
         expected_vector_mult)

    do i = 1, local_count

       b = local_block_catalog(i)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("local block has an invalid catalogue index")
       end if

       if (local_seen(b) /= 0) then
          call fail("local block occurs more than once")
       end if

       if (catalog_local_block(b) /= i) then
          call fail("local/catalogue inverse mapping mismatch")
       end if

       if (block_catalog(b)%owner /= rank) then
          call fail("local block has the wrong final owner")
       end if

       call get_local_block_identity( &
            i,id,root_domain,root_patch,level)

       call get_local_block_field_layout( &
            i,scalar_variable,n_scalar_variable, &
            vector_variable,field_level, &
            n_field_level,scalar_mult,vector_mult)

       if (id /= block_catalog(b)%id .or. &
            root_domain /= block_catalog(b)%root_domain .or. &
            root_patch /= block_catalog(b)%root_patch .or. &
            level /= block_catalog(b)%level) then
          call fail("local block identity does not match catalogue")
       end if

       if (scalar_variable /= expected_scalar_variable .or. &
            n_scalar_variable /= expected_n_scalar_variable .or. &
            vector_variable /= expected_vector_variable .or. &
            field_level /= expected_field_level .or. &
            n_field_level /= expected_n_field_level .or. &
            scalar_mult /= expected_scalar_mult .or. &
            vector_mult /= expected_vector_mult) then
          call fail("local block field layout does not match format")
       end if

       call check_local_block_storage(i,.true.)

       local_seen(b) = 1
       local_weight = local_weight + block_catalog(b)%weight

    end do

    expected_local = count(block_catalog%owner == rank)

    if (local_count /= expected_local) then
       call fail("local block count does not match final ownership")
    end if

    call MPI_Allreduce( &
         local_count,global_count,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block count")

    call MPI_Allreduce( &
         local_weight,global_weight,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block weight")

    call MPI_Allreduce( &
         local_seen,global_seen,size(local_seen),MPI_INTEGER, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block inventory")

    if (global_count /= size(block_catalog)) then
       call fail("global block count mismatch")
    end if

    if (global_weight /= sum(block_catalog%weight)) then
       call fail("global block weight mismatch")
    end if

    if (any(global_seen /= 1)) then
       call fail("global block ownership is not unique")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Standalone local block store for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  final-owner blocks = ", local_count
       write(6,'(a,i0)') &
            "  final-owner weight = ", local_weight
       write(6,'(a)') &
            "  component and serialization checks passed"
       write(6,'(a)') &
            "  persistent store readiness check passed"
       write(6,'(a)') &
            "  bidirectional local/catalogue mapping check passed"
       write(6,'(a)') &
            "  self-describing field metadata check passed"
       write(6,'(a,i0)') &
            "  scalar variables represented = ", &
            expected_n_scalar_variable
       write(6,'(a,i0)') &
            "  first sol field level = ", expected_field_level
       write(6,'(a,i0)') &
            "  sol field levels represented = ", &
            expected_n_field_level
       write(6,'(a,/)') &
            "  unique global inventory check passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Standalone global block objects verified = ", &
            global_count
       write(6,'(a,i0)') &
            "Standalone global block weight verified  = ", &
            global_weight
       write(6,'(a,/)') &
            "Final-owner block store is self-contained"
    end if

    deallocate(global_seen)
    deallocate(local_seen)

  end subroutine check_local_blocks


  subroutine check_block_field_inventory (verbose)
    ! Compare sol, sol_mean and wav_coeff interior fields in the
    ! migrated final-owner block store with the still-authoritative
    ! legacy Domain fields covered by the catalogue-rooted subtrees. The
    ! fixed coarse scaffold above those roots is not block storage.
    ! Only global counts and order-independent moments are compared,
    ! because block and Domain ownership differ after migration.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: d
    integer :: ierr
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: v_scalar
    integer :: v_vector
    integer :: k_field

    integer(int64) :: block_count_local(2)
    integer(int64) :: block_count_global(2)
    integer(int64) :: block_mean_count_local(2)
    integer(int64) :: block_mean_count_global(2)
    integer(int64) :: block_wavelet_count_local(2)
    integer(int64) :: block_wavelet_count_global(2)
    integer(int64) :: domain_count_local(2)
    integer(int64) :: domain_count_global(2)
    integer(int64) :: domain_mean_count_local(2)
    integer(int64) :: domain_mean_count_global(2)
    integer(int64) :: domain_wavelet_count_local(2)
    integer(int64) :: domain_wavelet_count_global(2)

    real(dp) :: block_moment_local(3,2)
    real(dp) :: block_moment_global(3,2)
    real(dp) :: block_mean_moment_local(3,2)
    real(dp) :: block_mean_moment_global(3,2)
    real(dp) :: block_wavelet_moment_local(3,2)
    real(dp) :: block_wavelet_moment_global(3,2)
    real(dp) :: domain_moment_local(3,2)
    real(dp) :: domain_moment_global(3,2)
    real(dp) :: domain_mean_moment_local(3,2)
    real(dp) :: domain_mean_moment_global(3,2)
    real(dp) :: domain_wavelet_moment_local(3,2)
    real(dp) :: domain_wavelet_moment_global(3,2)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("field inventory requested before block installation")
    end if

    call local_block_field_statistics( &
         block_count_local(1),block_count_local(2), &
         block_moment_local(:,1),block_moment_local(:,2))

    call local_block_wavelet_statistics( &
         block_wavelet_count_local(1),block_wavelet_count_local(2), &
         block_wavelet_moment_local(:,1), &
         block_wavelet_moment_local(:,2))

    call local_block_mean_field_statistics( &
         block_mean_count_local(1),block_mean_count_local(2), &
         block_mean_moment_local(:,1),block_mean_moment_local(:,2))

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,k_field, &
         n_field_level,mult_scalar,mult_vector)

    domain_count_local  = 0_int64
    domain_moment_local = 0.0_dp
    domain_mean_count_local  = 0_int64
    domain_mean_moment_local = 0.0_dp
    domain_wavelet_count_local  = 0_int64
    domain_wavelet_moment_local = 0.0_dp

    do b = 1, size(block_catalog)

       if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

       d = loc_id(block_catalog(b)%root_domain+1) + 1

       if (d < 1 .or. d > size(grid)) then
          call fail("invalid source Domain in field inventory")
       end if

       call accumulate_domain_subtree_fields( &
            d,block_catalog(b)%root_patch, &
            v_scalar,n_scalar_variable,v_vector,k_field,n_field_level, &
            mult_scalar,mult_vector, &
            domain_count_local,domain_moment_local, &
            domain_mean_count_local,domain_mean_moment_local, &
            domain_wavelet_count_local,domain_wavelet_moment_local)

    end do

    call MPI_Allreduce( &
         block_count_local,block_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block field counts")

    call MPI_Allreduce( &
         block_mean_count_local,block_mean_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block mean-field counts")

    call MPI_Allreduce( &
         block_wavelet_count_local,block_wavelet_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wavelet counts")

    call MPI_Allreduce( &
         domain_count_local,domain_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain field counts")

    call MPI_Allreduce( &
         domain_mean_count_local,domain_mean_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain mean-field counts")

    call MPI_Allreduce( &
         domain_wavelet_count_local,domain_wavelet_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wavelet counts")

    call MPI_Allreduce( &
         block_moment_local,block_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block field moments")

    call MPI_Allreduce( &
         block_mean_moment_local,block_mean_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block mean-field moments")

    call MPI_Allreduce( &
         block_wavelet_moment_local,block_wavelet_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wavelet moments")

    call MPI_Allreduce( &
         domain_moment_local,domain_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain field moments")

    call MPI_Allreduce( &
         domain_mean_moment_local,domain_mean_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain mean-field moments")

    call MPI_Allreduce( &
         domain_wavelet_moment_local,domain_wavelet_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wavelet moments")

    if (any(block_count_global /= domain_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain field-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block counts  = ", block_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain counts = ", domain_count_global
       end if

       call fail("block and Domain field counts differ")
    end if

    if (any(block_mean_count_global /= domain_mean_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean-field-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block mean counts  = ", block_mean_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain mean counts = ", domain_mean_count_global
       end if

       call fail("block and Domain mean-field counts differ")
    end if

    if (any(block_wavelet_count_global /= &
         domain_wavelet_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain wavelet-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block wavelet counts  = ", block_wavelet_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain wavelet counts = ", domain_wavelet_count_global
       end if

       call fail("block and Domain wavelet counts differ")
    end if

    if (.not. field_moments_match( &
         block_moment_global(:,1),domain_moment_global(:,1), &
         block_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain scalar-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_moment_global(:,1)
       end if

       call fail("block and Domain scalar moments differ")
    end if

    if (.not. field_moments_match( &
         block_moment_global(:,2),domain_moment_global(:,2), &
         block_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain vector-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_moment_global(:,2)
       end if

       call fail("block and Domain vector moments differ")
    end if

    if (.not. field_moments_match( &
         block_mean_moment_global(:,1), &
         domain_mean_moment_global(:,1), &
         block_mean_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean scalar-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_mean_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_mean_moment_global(:,1)
       end if

       call fail("block and Domain mean scalar moments differ")
    end if

    if (.not. field_moments_match( &
         block_mean_moment_global(:,2), &
         domain_mean_moment_global(:,2), &
         block_mean_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean vector-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_mean_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_mean_moment_global(:,2)
       end if

       call fail("block and Domain mean vector moments differ")
    end if

    if (.not. field_moments_match( &
         block_wavelet_moment_global(:,1), &
         domain_wavelet_moment_global(:,1), &
         block_wavelet_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain scalar wavelet-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_wavelet_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_wavelet_moment_global(:,1)
       end if

       call fail("block and Domain scalar wavelet moments differ")
    end if

    if (.not. field_moments_match( &
         block_wavelet_moment_global(:,2), &
         domain_wavelet_moment_global(:,2), &
         block_wavelet_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain vector wavelet-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_wavelet_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_wavelet_moment_global(:,2)
       end if

       call fail("block and Domain vector wavelet moments differ")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Read-only block field consumer for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  local block scalar values = ", block_count_local(1)
       write(6,'(a,i0)') &
            "  local block vector values = ", block_count_local(2)
       write(6,'(a,i0)') &
            "  local block wavelet scalar values = ", &
            block_wavelet_count_local(1)
       write(6,'(a,i0)') &
            "  local block wavelet vector values = ", &
            block_wavelet_count_local(2)
       write(6,'(a,i0)') &
            "  local block mean scalar values = ", &
            block_mean_count_local(1)
       write(6,'(a,i0)') &
            "  local block mean vector values = ", &
            block_mean_count_local(2)
       write(6,'(a,/)') &
            "  global sol/sol_mean/wav_coeff inventory checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global scalar interior values verified = ", &
            block_count_global(1)
       write(6,'(a,i0)') &
            "Global vector interior values verified = ", &
            block_count_global(2)
       write(6,'(a,i0)') &
            "Global wavelet scalar values verified = ", &
            block_wavelet_count_global(1)
       write(6,'(a,i0)') &
            "Global wavelet vector values verified = ", &
            block_wavelet_count_global(2)
       write(6,'(a,i0)') &
            "Global mean scalar interior values verified = ", &
            block_mean_count_global(1)
       write(6,'(a,i0)') &
            "Global mean vector interior values verified = ", &
            block_mean_count_global(2)
       write(6,'(a,/)') &
            "Block sol, sol_mean and wav_coeff match legacy Domain data"
    end if

  end subroutine check_block_field_inventory


  logical function field_moments_match (first,second,n_value) &
       result(match)
    ! Allow only the floating-point reassociation error introduced by
    ! summing an identical value inventory in different block/rank
    ! orders.

    implicit none

    real(dp), intent(in) :: first(3)
    real(dp), intent(in) :: second(3)
    integer(int64), intent(in) :: n_value

    real(dp) :: factor
    real(dp) :: scale

    factor = 256.0_dp * epsilon(1.0_dp) * &
         real(max(1_int64,n_value),dp)

    scale = max(1.0_dp,first(2),second(2))

    match = abs(first(1)-second(1)) <= factor*scale
    if (.not. match) return

    match = abs(first(2)-second(2)) <= factor*scale
    if (.not. match) return

    scale = max(1.0_dp,first(3),second(3))
    match = abs(first(3)-second(3)) <= factor*scale

  end function field_moments_match


  recursive subroutine accumulate_domain_subtree_fields ( &
       d,p,v_scalar,n_scalar_variable,v_vector,k_field,n_field_level, &
       mult_scalar,mult_vector,field_count,field_moment, &
       mean_field_count,mean_field_moment, &
       wavelet_field_count,wavelet_field_moment)
    ! Accumulate one catalogue-rooted subtree from the authoritative
    ! Domain representation using the same patch coverage copied into
    ! Block_Data.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: v_scalar
    integer, intent(in) :: n_scalar_variable
    integer, intent(in) :: v_vector
    integer, intent(in) :: k_field
    integer, intent(in) :: n_field_level
    integer, intent(in) :: mult_scalar
    integer, intent(in) :: mult_vector

    integer(int64), intent(inout) :: field_count(2)
    real(dp), intent(inout) :: field_moment(3,2)
    integer(int64), intent(inout) :: mean_field_count(2)
    real(dp), intent(inout) :: mean_field_moment(3,2)
    integer(int64), intent(inout) :: wavelet_field_count(2)
    real(dp), intent(inout) :: wavelet_field_moment(3,2)

    integer :: c
    integer :: field_level
    integer :: level_slot
    integer :: n_value
    integer :: p_child
    integer :: scalar_id
    integer :: scalar_slot
    integer :: start

    if (p < 0 .or. p >= grid(d)%patch%length) then
       call fail("invalid patch in Domain subtree field inventory")
    end if

    if (grid(d)%patch%elts(p+1)%deleted) return

    do scalar_slot = 1, n_scalar_variable

       scalar_id = v_scalar + scalar_slot - 1

       do level_slot = 1, n_field_level

          field_level = k_field + level_slot - 1
          start = mult_scalar * grid(d)%patch%elts(p+1)%elts_start
          n_value = mult_scalar * PATCH_SIZE**2

          if (start < 0 .or. &
               start+n_value > &
               size(sol(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy scalar patch extent is invalid")
          end if

          field_count(1) = field_count(1) + int(n_value,int64)
          field_moment(1,1) = field_moment(1,1) + &
               sum(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          field_moment(2,1) = field_moment(2,1) + &
               sum(abs(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          field_moment(3,1) = field_moment(3,1) + &
               sum(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

          if (start+n_value > &
               size(wav_coeff(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy scalar wavelet patch extent is invalid")
          end if

          wavelet_field_count(1) = &
               wavelet_field_count(1) + int(n_value,int64)
          wavelet_field_moment(1,1) = wavelet_field_moment(1,1) + &
               sum(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          wavelet_field_moment(2,1) = wavelet_field_moment(2,1) + &
               sum(abs(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          wavelet_field_moment(3,1) = wavelet_field_moment(3,1) + &
               sum(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

          if (start+n_value > &
               size(sol_mean(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy mean scalar patch extent is invalid")
          end if

          mean_field_count(1) = &
               mean_field_count(1) + int(n_value,int64)
          mean_field_moment(1,1) = mean_field_moment(1,1) + &
               sum(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          mean_field_moment(2,1) = mean_field_moment(2,1) + &
               sum(abs(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          mean_field_moment(3,1) = mean_field_moment(3,1) + &
               sum(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

       end do

    end do

    do level_slot = 1, n_field_level

       field_level = k_field + level_slot - 1
       start = mult_vector * grid(d)%patch%elts(p+1)%elts_start
       n_value = mult_vector * PATCH_SIZE**2

       if (start < 0 .or. &
            start+n_value > &
            size(sol(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy vector patch extent is invalid")
       end if

       field_count(2) = field_count(2) + int(n_value,int64)
       field_moment(1,2) = field_moment(1,2) + &
            sum(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       field_moment(2,2) = field_moment(2,2) + &
            sum(abs(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       field_moment(3,2) = field_moment(3,2) + &
            sum(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

       if (start+n_value > &
            size(wav_coeff(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy vector wavelet patch extent is invalid")
       end if

       wavelet_field_count(2) = &
            wavelet_field_count(2) + int(n_value,int64)
       wavelet_field_moment(1,2) = wavelet_field_moment(1,2) + &
            sum(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       wavelet_field_moment(2,2) = wavelet_field_moment(2,2) + &
            sum(abs(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       wavelet_field_moment(3,2) = wavelet_field_moment(3,2) + &
            sum(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

       if (start+n_value > &
            size(sol_mean(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy mean vector patch extent is invalid")
       end if

       mean_field_count(2) = &
            mean_field_count(2) + int(n_value,int64)
       mean_field_moment(1,2) = mean_field_moment(1,2) + &
            sum(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       mean_field_moment(2,2) = mean_field_moment(2,2) + &
            sum(abs(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       mean_field_moment(3,2) = mean_field_moment(3,2) + &
            sum(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

    end do

    do c = 1, N_CHDRN

       p_child = grid(d)%patch%elts(p+1)%children(c)
       if (p_child == 0) cycle

       call accumulate_domain_subtree_fields( &
            d,p_child,v_scalar,n_scalar_variable,v_vector,k_field, &
            n_field_level,mult_scalar,mult_vector, &
            field_count,field_moment, &
            mean_field_count,mean_field_moment, &
            wavelet_field_count,wavelet_field_moment)

    end do

  end subroutine accumulate_domain_subtree_fields


  subroutine build_block_migration_manifest (manifest)
    ! Build and exchange a migration manifest containing only indices
    ! into the replicated block_catalog. No block payload is moved and
    ! neither grid nor any persistent block object is modified.

    implicit none

    type(Block_Migration_Manifest), intent(out) :: manifest

    integer :: b
    integer :: destination
    integer :: ierr
    integer :: r
    integer :: slot
    integer :: source
    integer, allocatable :: cursor(:)

    if (.not. allocated(block_catalog)) then
       call fail("block_catalog is not allocated")
    end if

    if (.not. allocated(owner)) then
       call fail("domain owner map is not allocated")
    end if

    allocate(manifest%send_count(n_process))
    allocate(manifest%recv_count(n_process))
    allocate(manifest%send_displ(n_process))
    allocate(manifest%recv_displ(n_process))

    manifest%send_count = 0
    manifest%recv_count = 0
    manifest%send_displ = 0
    manifest%recv_displ = 0
    manifest%validated  = .false.

    do b = 1, size(block_catalog)
       source      = source_rank(b)
       destination = block_catalog(b)%owner

       call check_rank(destination,"destination owner")

       if (source == rank .and. destination /= rank) then
          manifest%send_count(destination+1) = &
               manifest%send_count(destination+1) + 1
       end if
    end do

    call MPI_Alltoall( &
         manifest%send_count, 1, MPI_INTEGER, &
         manifest%recv_count, 1, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoall migration counts")

    do r = 2, n_process
       manifest%send_displ(r) = manifest%send_displ(r-1) + &
            manifest%send_count(r-1)
       manifest%recv_displ(r) = manifest%recv_displ(r-1) + &
            manifest%recv_count(r-1)
    end do

    manifest%n_send = sum(manifest%send_count)
    manifest%n_recv = sum(manifest%recv_count)

    ! A one-element allocation avoids passing an unallocated or
    ! zero-sized buffer to MPI implementations with older Fortran
    ! choice-buffer handling. n_send/n_recv remain the true lengths.
    allocate(manifest%send_block(max(1,manifest%n_send)))
    allocate(manifest%recv_block(max(1,manifest%n_recv)))

    manifest%send_block = 0
    manifest%recv_block = 0

    allocate(cursor(n_process))
    cursor = manifest%send_displ

    do b = 1, size(block_catalog)
       source      = source_rank(b)
       destination = block_catalog(b)%owner

       if (source /= rank .or. destination == rank) cycle

       slot = cursor(destination+1) + 1

       if (slot < 1 .or. slot > manifest%n_send) then
          call fail("send-manifest fill position is invalid")
       end if

       manifest%send_block(slot) = b
       cursor(destination+1) = cursor(destination+1) + 1
    end do

    do r = 1, n_process
       if (cursor(r) /= manifest%send_displ(r) + &
            manifest%send_count(r)) then
          call fail("send-manifest count mismatch")
       end if
    end do

    deallocate(cursor)

    call MPI_Alltoallv( &
         manifest%send_block, manifest%send_count, &
         manifest%send_displ, MPI_INTEGER, &
         manifest%recv_block, manifest%recv_count, &
         manifest%recv_displ, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv migration manifest")

  end subroutine build_block_migration_manifest


  subroutine check_block_migration_manifest (manifest, verbose)
    ! Verify that every manifest entry has the expected source rank and
    ! this rank as its prospective destination. Global counts and a
    ! catalogue-index checksum must agree with the replicated catalogue.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    logical, optional, intent(in)                   :: verbose

    integer :: b
    integer :: destination
    integer :: expected_changed
    integer :: expected_checksum
    integer :: global_recv
    integer :: global_recv_checksum
    integer :: global_send
    integer :: global_send_checksum
    integer :: ierr
    integer :: local_recv_checksum
    integer :: local_send_checksum
    integer :: pos
    integer :: r
    integer :: source
    logical :: print_local
    logical, allocatable :: seen(:)

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. allocated(manifest%send_count) .or. &
         .not. allocated(manifest%recv_count) .or. &
         .not. allocated(manifest%send_displ) .or. &
         .not. allocated(manifest%recv_displ) .or. &
         .not. allocated(manifest%send_block) .or. &
         .not. allocated(manifest%recv_block)) then
       call fail("migration manifest is incomplete")
    end if

    if (size(manifest%send_count) /= n_process .or. &
         size(manifest%recv_count) /= n_process .or. &
         size(manifest%send_displ) /= n_process .or. &
         size(manifest%recv_displ) /= n_process) then
       call fail("migration manifest has an invalid rank dimension")
    end if

    if (manifest%n_send /= sum(manifest%send_count) .or. &
         manifest%n_recv /= sum(manifest%recv_count)) then
       call fail("migration manifest has an invalid local total")
    end if

    allocate(seen(size(block_catalog)))
    seen = .false.

    do r = 0, n_process-1
       do pos = manifest%recv_displ(r+1) + 1, &
            manifest%recv_displ(r+1) + manifest%recv_count(r+1)

          b = manifest%recv_block(pos)

          if (b < 1 .or. b > size(block_catalog)) then
             call fail("received block-catalog index is invalid")
          end if

          if (seen(b)) then
             call fail("duplicate block in received migration manifest")
          end if
          seen(b) = .true.

          source      = source_rank(b)
          destination = block_catalog(b)%owner

          if (source /= r) then
             call fail("received block has the wrong source rank")
          end if

          if (destination /= rank) then
             call fail("received block has the wrong destination rank")
          end if

          if (source == destination) then
             call fail("retained block appears in migration manifest")
          end if
       end do
    end do

    deallocate(seen)

    do pos = 1, manifest%n_send
       b = manifest%send_block(pos)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("sent block-catalog index is invalid")
       end if

       if (source_rank(b) /= rank) then
          call fail("sent block does not originate on this rank")
       end if

       if (block_catalog(b)%owner == rank) then
          call fail("retained block appears in send manifest")
       end if
    end do

    local_send_checksum = 0
    local_recv_checksum = 0

    if (manifest%n_send > 0) then
       local_send_checksum = &
            sum(manifest%send_block(1:manifest%n_send))
    end if

    if (manifest%n_recv > 0) then
       local_recv_checksum = &
            sum(manifest%recv_block(1:manifest%n_recv))
    end if

    call MPI_Allreduce( &
         manifest%n_send, global_send, 1, MPI_INTEGER, MPI_SUM, &
         comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce sent migration count")

    call MPI_Allreduce( &
         manifest%n_recv, global_recv, 1, MPI_INTEGER, MPI_SUM, &
         comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce received migration count")

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce sent manifest checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce received manifest checksum")

    expected_changed  = 0
    expected_checksum = 0

    do b = 1, size(block_catalog)
       if (source_rank(b) == block_catalog(b)%owner) cycle
       expected_changed  = expected_changed + 1
       expected_checksum = expected_checksum + b
    end do

    if (global_send /= expected_changed .or. &
         global_recv /= expected_changed) then
       call fail("global migration count does not match block_catalog")
    end if

    if (global_send_checksum /= expected_checksum .or. &
         global_recv_checksum /= expected_checksum) then
       call fail("global migration checksum does not match block_catalog")
    end if

    manifest%validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration manifest for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  blocks sent     = ", manifest%n_send
       write(6,'(a,i0)') &
            "  blocks received = ", manifest%n_recv
       write(6,'(a)') &
            "  source/destination checks passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global migrating blocks verified = ", expected_changed
       write(6,'(a,/)') &
            "Block migration manifest exchange passed"
    end if

  end subroutine check_block_migration_manifest


  subroutine exchange_block_migration_sizes (manifest, send_nbyte, verbose)
    ! Exchange only the packed byte count associated with each entry in
    ! the already validated migration manifest. Block payloads are not
    ! communicated and persistent block storage is not modified.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    integer, intent(in)                            :: send_nbyte(:)
    logical, optional, intent(in)                  :: verbose

    integer :: global_max_recv
    integer :: global_max_send
    integer :: ierr
    integer :: pos
    integer(int64) :: global_recv_checksum
    integer(int64) :: global_recv_nbyte
    integer(int64) :: global_send_checksum
    integer(int64) :: global_send_nbyte
    integer(int64) :: local_recv_checksum
    integer(int64) :: local_send_checksum
    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. manifest%validated) then
       call fail("migration manifest must be validated before sizes")
    end if

    if (size(send_nbyte) /= manifest%n_send) then
       call fail("outgoing byte-count array has the wrong extent")
    end if

    if (manifest%n_send > 0) then
       if (any(send_nbyte <= 0)) then
          call fail("outgoing packed block has a nonpositive byte count")
       end if
    end if

    if (allocated(manifest%send_nbyte)) then
       deallocate(manifest%send_nbyte)
    end if

    if (allocated(manifest%recv_nbyte)) then
       deallocate(manifest%recv_nbyte)
    end if

    allocate(manifest%send_nbyte(max(1,manifest%n_send)))
    allocate(manifest%recv_nbyte(max(1,manifest%n_recv)))

    manifest%send_nbyte = 0
    manifest%recv_nbyte = 0

    if (manifest%n_send > 0) then
       manifest%send_nbyte(1:manifest%n_send) = send_nbyte
    end if

    call MPI_Alltoallv( &
         manifest%send_nbyte, manifest%send_count, &
         manifest%send_displ, MPI_INTEGER, &
         manifest%recv_nbyte, manifest%recv_count, &
         manifest%recv_displ, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv packed block sizes")

    if (manifest%n_recv > 0) then
       if (any(manifest%recv_nbyte(1:manifest%n_recv) <= 0)) then
          call fail("received packed block has a nonpositive byte count")
       end if
    end if

    manifest%total_send_nbyte = 0_int64
    manifest%total_recv_nbyte = 0_int64
    manifest%max_send_nbyte   = 0
    manifest%max_recv_nbyte   = 0

    if (manifest%n_send > 0) then
       manifest%total_send_nbyte = &
            sum(int(manifest%send_nbyte(1:manifest%n_send),int64))
       manifest%max_send_nbyte = &
            maxval(manifest%send_nbyte(1:manifest%n_send))
    end if

    if (manifest%n_recv > 0) then
       manifest%total_recv_nbyte = &
            sum(int(manifest%recv_nbyte(1:manifest%n_recv),int64))
       manifest%max_recv_nbyte = &
            maxval(manifest%recv_nbyte(1:manifest%n_recv))
    end if

    local_send_checksum = 0_int64
    local_recv_checksum = 0_int64

    do pos = 1, manifest%n_send
       local_send_checksum = local_send_checksum + &
            int(manifest%send_block(pos),int64) * &
            int(manifest%send_nbyte(pos),int64)
    end do

    do pos = 1, manifest%n_recv
       local_recv_checksum = local_recv_checksum + &
            int(manifest%recv_block(pos),int64) * &
            int(manifest%recv_nbyte(pos),int64)
    end do

    call MPI_Allreduce( &
         manifest%total_send_nbyte, global_send_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing packed bytes")

    call MPI_Allreduce( &
         manifest%total_recv_nbyte, global_recv_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming packed bytes")

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing size checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming size checksum")

    call MPI_Allreduce( &
         manifest%max_send_nbyte, global_max_send, 1, &
         MPI_INTEGER, MPI_MAX, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce maximum outgoing block size")

    call MPI_Allreduce( &
         manifest%max_recv_nbyte, global_max_recv, 1, &
         MPI_INTEGER, MPI_MAX, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce maximum incoming block size")

    if (global_send_nbyte /= global_recv_nbyte) then
       call fail("global outgoing and incoming byte totals differ")
    end if

    if (global_send_checksum /= global_recv_checksum) then
       call fail("global block/size checksums differ")
    end if

    if (global_max_send /= global_max_recv) then
       call fail("global maximum outgoing and incoming sizes differ")
    end if

    manifest%sizes_validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration sizes for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  outgoing blocks       = ", manifest%n_send
       write(6,'(a,i0)') &
            "  outgoing packed bytes = ", manifest%total_send_nbyte
       write(6,'(a,i0)') &
            "  incoming blocks       = ", manifest%n_recv
       write(6,'(a,i0)') &
            "  incoming packed bytes = ", manifest%total_recv_nbyte
       write(6,'(a)') &
            "  packed-size routing checks passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global packed migration bytes verified = ", &
            global_send_nbyte
       write(6,'(a,i0)') &
            "Maximum packed migrating block bytes   = ", &
            global_max_send
       write(6,'(a,/)') &
            "Block migration size exchange passed"
    end if

  end subroutine exchange_block_migration_sizes


  subroutine exchange_block_migration_payloads ( &
       manifest, send_payload, verbose)
    ! Exchange the serialized payloads using the already validated
    ! per-block byte counts. The received byte stream remains separate
    ! from all source and installed block storage.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    integer(int8), intent(in)                      :: send_payload(:)
    logical, optional, intent(in)                  :: verbose

    integer :: first
    integer :: ierr
    integer :: last
    integer :: n_recv_byte
    integer :: n_send_byte
    integer :: r
    integer(int64) :: global_recv_checksum
    integer(int64) :: global_send_checksum
    integer(int64) :: global_send_nbyte
    integer(int64) :: local_recv_checksum
    integer(int64) :: local_send_checksum
    integer(int64) :: rank_nbyte
    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. manifest%validated .or. &
         .not. manifest%sizes_validated) then
       call fail("manifest and sizes must precede payload exchange")
    end if

    if (manifest%total_send_nbyte > int(huge(0),int64) .or. &
         manifest%total_recv_nbyte > int(huge(0),int64)) then
       call fail("packed byte total exceeds MPI default-count range")
    end if

    n_send_byte = int(manifest%total_send_nbyte)
    n_recv_byte = int(manifest%total_recv_nbyte)

    if (n_send_byte > 0) then
       if (size(send_payload) /= n_send_byte) then
          call fail("outgoing payload buffer has the wrong extent")
       end if
    else
       if (size(send_payload) < 1) then
          call fail("zero-byte sender requires a one-element MPI buffer")
       end if
    end if

    if (allocated(manifest%send_byte_count)) then
       deallocate(manifest%send_byte_count)
    end if
    if (allocated(manifest%recv_byte_count)) then
       deallocate(manifest%recv_byte_count)
    end if
    if (allocated(manifest%send_byte_displ)) then
       deallocate(manifest%send_byte_displ)
    end if
    if (allocated(manifest%recv_byte_displ)) then
       deallocate(manifest%recv_byte_displ)
    end if
    if (allocated(manifest%recv_payload)) then
       deallocate(manifest%recv_payload)
    end if

    allocate(manifest%send_byte_count(n_process))
    allocate(manifest%recv_byte_count(n_process))
    allocate(manifest%send_byte_displ(n_process))
    allocate(manifest%recv_byte_displ(n_process))

    manifest%send_byte_count = 0
    manifest%recv_byte_count = 0
    manifest%send_byte_displ = 0
    manifest%recv_byte_displ = 0
    manifest%payload_validated = .false.

    do r = 1, n_process

       first = manifest%send_displ(r) + 1
       last  = manifest%send_displ(r) + manifest%send_count(r)

       rank_nbyte = 0_int64
       if (last >= first) then
          rank_nbyte = &
               sum(int(manifest%send_nbyte(first:last),int64))
       end if

       if (rank_nbyte > int(huge(0),int64)) then
          call fail("per-rank outgoing payload exceeds MPI count range")
       end if

       manifest%send_byte_count(r) = int(rank_nbyte)

       first = manifest%recv_displ(r) + 1
       last  = manifest%recv_displ(r) + manifest%recv_count(r)

       rank_nbyte = 0_int64
       if (last >= first) then
          rank_nbyte = &
               sum(int(manifest%recv_nbyte(first:last),int64))
       end if

       if (rank_nbyte > int(huge(0),int64)) then
          call fail("per-rank incoming payload exceeds MPI count range")
       end if

       manifest%recv_byte_count(r) = int(rank_nbyte)

    end do

    do r = 2, n_process
       manifest%send_byte_displ(r) = &
            manifest%send_byte_displ(r-1) + &
            manifest%send_byte_count(r-1)
       manifest%recv_byte_displ(r) = &
            manifest%recv_byte_displ(r-1) + &
            manifest%recv_byte_count(r-1)
    end do

    if (sum(int(manifest%send_byte_count,int64)) /= &
         manifest%total_send_nbyte) then
       call fail("outgoing per-rank byte counts do not sum correctly")
    end if

    if (sum(int(manifest%recv_byte_count,int64)) /= &
         manifest%total_recv_nbyte) then
       call fail("incoming per-rank byte counts do not sum correctly")
    end if

    allocate(manifest%recv_payload(max(1,n_recv_byte)))
    manifest%recv_payload = 0_int8

    call MPI_Alltoallv( &
         send_payload, manifest%send_byte_count, &
         manifest%send_byte_displ, MPI_BYTE, &
         manifest%recv_payload, manifest%recv_byte_count, &
         manifest%recv_byte_displ, MPI_BYTE, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv packed block payloads")

    local_send_checksum = 0_int64
    local_recv_checksum = 0_int64

    if (n_send_byte > 0) then
       local_send_checksum = &
            sum(int(send_payload(1:n_send_byte),int64))
    end if

    if (n_recv_byte > 0) then
       local_recv_checksum = &
            sum(int(manifest%recv_payload(1:n_recv_byte),int64))
    end if

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing payload checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming payload checksum")

    call MPI_Allreduce( &
         manifest%total_send_nbyte, global_send_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce transferred payload bytes")

    if (global_send_checksum /= global_recv_checksum) then
       call fail("global payload byte checksums differ")
    end if

    manifest%payload_validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration payload transport for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  outgoing payload bytes = ", &
            manifest%total_send_nbyte
       write(6,'(a,i0)') &
            "  incoming payload bytes = ", &
            manifest%total_recv_nbyte
       write(6,'(a)') &
            "  payload transport checksum passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global packed payload bytes transferred = ", &
            global_send_nbyte
       write(6,'(a,/)') &
            "Block migration payload transport passed"
    end if

  end subroutine exchange_block_migration_payloads


  subroutine clear_block_migration_manifest (manifest)
    ! Release all count, routing, size and payload staging storage after
    ! the received blocks have been installed and independently checked.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest

    if (allocated(manifest%send_count)) then
       deallocate(manifest%send_count)
    end if
    if (allocated(manifest%recv_count)) then
       deallocate(manifest%recv_count)
    end if
    if (allocated(manifest%send_displ)) then
       deallocate(manifest%send_displ)
    end if
    if (allocated(manifest%recv_displ)) then
       deallocate(manifest%recv_displ)
    end if
    if (allocated(manifest%send_block)) then
       deallocate(manifest%send_block)
    end if
    if (allocated(manifest%recv_block)) then
       deallocate(manifest%recv_block)
    end if
    if (allocated(manifest%send_nbyte)) then
       deallocate(manifest%send_nbyte)
    end if
    if (allocated(manifest%recv_nbyte)) then
       deallocate(manifest%recv_nbyte)
    end if
    if (allocated(manifest%send_byte_count)) then
       deallocate(manifest%send_byte_count)
    end if
    if (allocated(manifest%recv_byte_count)) then
       deallocate(manifest%recv_byte_count)
    end if
    if (allocated(manifest%send_byte_displ)) then
       deallocate(manifest%send_byte_displ)
    end if
    if (allocated(manifest%recv_byte_displ)) then
       deallocate(manifest%recv_byte_displ)
    end if
    if (allocated(manifest%recv_payload)) then
       deallocate(manifest%recv_payload)
    end if

    manifest%n_send = 0
    manifest%n_recv = 0
    manifest%total_send_nbyte = 0_int64
    manifest%total_recv_nbyte = 0_int64
    manifest%max_send_nbyte = 0
    manifest%max_recv_nbyte = 0
    manifest%validated = .false.
    manifest%sizes_validated = .false.
    manifest%payload_validated = .false.

  end subroutine clear_block_migration_manifest


  integer function source_rank (b) result(source)

    implicit none

    integer, intent(in) :: b
    integer             :: root_domain

    if (b < 1 .or. b > size(block_catalog)) then
       call fail("block-catalog index is invalid")
    end if

    root_domain = block_catalog(b)%root_domain

    if (root_domain < 0 .or. root_domain >= size(owner)) then
       call fail("block root-domain index is invalid")
    end if

    source = owner(root_domain+1)
    call check_rank(source,"source domain owner")

  end function source_rank


  subroutine check_rank (r, description)

    implicit none

    integer, intent(in)          :: r
    character(*), intent(in)     :: description

    if (r < 0 .or. r >= n_process) then
       call fail(trim(description)//" is outside the MPI rank range")
    end if

  end subroutine check_rank


  subroutine check_mpi (ierr, operation)

    implicit none

    integer, intent(in)      :: ierr
    character(*), intent(in) :: operation

    if (ierr /= MPI_SUCCESS) then
       write(error_unit,'(a,i0,2a,i0)') &
            "Rank ", rank, ": ", trim(operation)// &
            " failed with MPI error ", ierr
       call abort_run
    end if

  end subroutine check_mpi


  subroutine fail (message)

    implicit none

    character(*), intent(in) :: message

    write(error_unit,'(a,i0,2a)') &
         "Rank ", rank, ": parallel_block_mpi_mod: ", trim(message)
    call abort_run

  end subroutine fail


end module parallel_block_mpi_mod
