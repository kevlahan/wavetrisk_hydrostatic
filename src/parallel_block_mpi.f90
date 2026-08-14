module parallel_block_mpi_mod

  use iso_fortran_env, only : error_unit, int8, int64
  use mpi_f08,        only : MPI_Allreduce, MPI_Alltoall, MPI_Alltoallv, &
       MPI_BYTE, MPI_INTEGER, MPI_INTEGER8, MPI_MAX, MPI_SUCCESS, MPI_SUM

  use arch_mod, only : abort_run, block_catalog, comm, n_process, &
       owner, rank

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

contains


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
