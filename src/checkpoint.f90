module checkpoint_mod

  ! Provides checkpointing routines

  use mpi_f08
  use iso_fortran_env, only: error_unit, int8, int64

  use kind_mod,   only : dp
  use shared_mod, only : N_BDRY,  N_CHDRN, N_VARIABLE, MULT, cp_idx, istep_cumul, itime, iwrite, level_end, &
     min_level, max_level,  n_glo_domain, run_id, scalars, threshold, time, vert_diffuse, zmin, zmax, zlevels, NONE, z_null
  
  use arch_mod,         only : barrier, comm, glo_id, n_process, rank, cp_load
  use comm_mod,         only : domain_load
  use comm_mpi_mod,     only : update_bdry
  use domain_mod,       only : Domain, Float_Field, grid, idx, sol, tke, scalar, wav_coeff, wav_tke, wc_s
  use domain_ops_mod,   only : apply_interscale_d, apply_to_pole_d 
  use init_mod,         only : dump, load
  use patch_mod,        only : PATCH_SIZE
  use refine_patch_mod, only : check_child_required, post_refine, refine_patch1, refine_patch2
  use wavelet_mod,      only : Restrict_scalar

  implicit none 

  private
  public :: dump_adapt_mpi, load_adapt_mpi, read_checkpoint_directory


  ! Magic number identifying a WAVETRISK checkpoint file.
  ! (ASCII encoding of the first eight characters of WAVETRISK)
  integer(int64), parameter :: CP_MAGIC = &
       int(z'5741564554524953', int64)

  ! Version number of the checkpoint file layout.
  integer(int64), parameter :: CP_VERSION = 1_int64

  ! Absolute byte offset of each global domain within the checkpoint file.
  integer(int64), allocatable :: cp_offset(:)

  ! Serialized size (bytes) of each global domain.
  integer(int64), allocatable :: cp_nbytes(:)

  ! True if the restart .bin file was temporarily decompressed
  ! from the persistent .bin.zst checkpoint.
  logical :: cp_temporary_bin = .false.

  ! Size of the fixed checkpoint header (magic number, version,
  ! number of global domains).
  integer(int64), parameter :: CP_HEADER_BYTES = 24_int64

  ! Byte offset of the domain load array in the checkpoint file.
  integer(int64), parameter :: CP_LOAD_POS = &
       CP_HEADER_BYTES

  ! Byte offset of the domain offset array in the checkpoint file.
  integer(int64), parameter :: CP_OFFSET_POS = &
       CP_LOAD_POS + &
       int(storage_size(0)/8, int64) * int(N_GLO_DOMAIN, int64)

  ! Byte offset of the domain size array in the checkpoint file.
  integer(int64), parameter :: CP_NBYTES_POS = &
       CP_OFFSET_POS + &
       8_int64 * int(N_GLO_DOMAIN, int64)

  ! Byte offset of the serialized domain data.
  integer(int64), parameter :: CP_DATA_POS = &
       CP_NBYTES_POS + &
       8_int64 * int(N_GLO_DOMAIN, int64)

  
contains


  subroutine require_mpi_transfer_count ( &
       status,datatype,expected_count,context)
    ! MPI-IO may return MPI_SUCCESS after a short read. Convert that
    ! otherwise-latent truncation into an error at the read site.

    implicit none

    type(MPI_Status), intent(in) :: status
    type(MPI_Datatype), intent(in) :: datatype
    integer, intent(in) :: expected_count
    character(len=*), intent(in) :: context

    integer :: actual_count
    integer :: ierr

    call MPI_Get_count(status,datatype,actual_count,ierr)
    if (ierr /= MPI_SUCCESS) then
       write(error_unit,'(a)') trim(context)// &
            ": unable to query MPI transfer count"
       error stop "checkpoint MPI transfer count query failed"
    end if
    if (actual_count /= expected_count) then
       write(error_unit,'(2a,2(a,i0))') trim(context), &
            ": short MPI transfer", &
            ", expected count = ",expected_count, &
            ", actual count = ",actual_count
       error stop "checkpoint MPI transfer was incomplete"
    end if

  end subroutine require_mpi_transfer_count


  subroutine validate_checkpoint_file_extent (fh,context)
    ! Require every directory record to lie within an exactly sized file.
    ! This catches stale decompressions and inconsistent metadata before any
    ! rank-local reconstruction begins.

    implicit none

    type(MPI_File), intent(in) :: fh
    character(len=*), intent(in) :: context

    integer :: gid
    integer :: ierr
    integer(MPI_OFFSET_KIND) :: mpi_file_size
    integer(int64) :: expected_size
    integer(int64) :: file_size

    if (.not. allocated(cp_offset) .or. .not. allocated(cp_nbytes)) &
         error stop "checkpoint extent validation has no directory"
    if (size(cp_offset) /= N_GLO_DOMAIN .or. &
         size(cp_nbytes) /= N_GLO_DOMAIN) &
         error stop "checkpoint extent validation directory is invalid"

    expected_size = CP_DATA_POS
    do gid = 1,N_GLO_DOMAIN
       if (cp_offset(gid) < CP_DATA_POS .or. &
            cp_nbytes(gid) <= 0_int64) then
          write(error_unit,'(2a,a,i0)') trim(context), &
               ": invalid checkpoint directory record", &
               ", global domain = ",gid-1
          error stop "checkpoint directory record is invalid"
       end if
       if (cp_offset(gid) > huge(0_int64)-cp_nbytes(gid)) &
            error stop "checkpoint directory extent overflows int64"
       expected_size = max( &
            expected_size,cp_offset(gid)+cp_nbytes(gid))
    end do

    call MPI_File_get_size(fh,mpi_file_size,ierr)
    if (ierr /= MPI_SUCCESS) &
         error stop "checkpoint MPI file-size query failed"
    file_size = int(mpi_file_size,int64)
    if (file_size /= expected_size) then
       write(error_unit,'(2a,2(a,i0))') trim(context), &
            ": checkpoint file extent differs", &
            ", expected bytes = ",expected_size, &
            ", actual bytes = ",file_size
       error stop "checkpoint file is truncated or inconsistent"
    end if

  end subroutine validate_checkpoint_file_extent


  subroutine dump_adapt_mpi (id)
    ! Save adaptive checkpoint data to one shared MPI-IO file.
    !
    ! Each MPI rank serializes all locally owned domains into one
    ! temporary stream.  The rank-local streams are then written
    ! collectively into one shared checkpoint file.
    !
    ! The checkpoint directory contains, for every global domain:
    !
    !   cp_load   : computational load
    !   cp_offset : absolute byte offset in checkpoint file
    !   cp_nbytes : serialized domain size in bytes
    !
    ! After the shared checkpoint has been closed collectively,
    ! rank zero compresses it using zstd -1 and removes the
    ! uncompressed .bin file.
    !
    ! NOTE: Checkpoint generation modifies the adaptive grid structure
    ! by deleting patches that are not required in the adjacent zone.

    implicit none

    integer, intent(in) :: id

    integer :: d, k, v
    integer :: fid
    integer :: ierr
    integer :: cmdstat, exitstat

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: p0, p1
    integer(int64) :: rank_bytes
    integer(int64) :: rank_disp

    integer,        allocatable :: gid_loc(:)
    integer,        allocatable :: load_loc(:)
    integer(int64), allocatable :: off_loc(:)
    integer(int64), allocatable :: len_loc(:)

    integer(int8), allocatable :: buf(:)

    character(4)              :: cp4
    character(:), allocatable :: filename
    character(:), allocatable :: cmd


    ! ---------------------------------------------------------------
    ! Prepare checkpoint data
    ! ---------------------------------------------------------------

    call update_bdry ( &
         wav_coeff(scalars(1):scalars(2),zmin:zmax), &
         NONE, 964)

    if (vert_diffuse) then
       call update_bdry (wav_tke, NONE, 965)
    end if


    ! Restrict physical solution to coarsest scale.

    do d = 1, size(grid)

       do k = zmin, zmax
          do v = scalars(1), scalars(2)

             scalar => sol(v,k)%data(d)%elts
             wc_s   => wav_coeff(v,k)%data(d)%elts

             call apply_interscale_d ( &
                  Restrict_scalar, &
                  grid(d), &
                  min_level-1, &
                  k, 0, 1)

             nullify (scalar, wc_s)

          end do
       end do

       if (vert_diffuse) then

          do k = 1, zlevels

             scalar => tke(k)%data(d)%elts
             wc_s   => wav_tke(k)%data(d)%elts

             call apply_interscale_d ( &
                  Restrict_scalar, &
                  grid(d), &
                  min_level-1, &
                  k, 0, 1)

             nullify (scalar, wc_s)

          end do

       end if

    end do


    ! ---------------------------------------------------------------
    ! Serialize all locally owned domains into one scratch stream
    ! ---------------------------------------------------------------

    allocate(gid_loc(size(grid)))
    allocate(load_loc(size(grid)))
    allocate(off_loc(size(grid)))
    allocate(len_loc(size(grid)))

    open( &
         newunit=fid, &
         status='scratch', &
         form='unformatted', &
         access='stream', &
         action='readwrite')


    do d = 1, size(grid)

       gid_loc(d)  = glo_id(rank+1,d)
       load_loc(d) = domain_load(grid(d))

       ! Position immediately before this domain.
       inquire(unit=fid, pos=p0)

       ! Complete stage-2 domain serialization.
       call dump_domain(fid, d)

       ! Position immediately after this domain.
       inquire(unit=fid, pos=p1)

       ! Fortran stream positions are one based.
       off_loc(d) = p0 - 1_int64
       len_loc(d) = p1 - p0

    end do


    ! Total rank-local serialized payload size.

    inquire(unit=fid, pos=p1)

    rank_bytes = p1 - 1_int64


    ! ---------------------------------------------------------------
    ! Determine global payload offsets
    ! ---------------------------------------------------------------

    rank_disp = 0_int64

    call MPI_Exscan( &
         rank_bytes, rank_disp, &
         1, MPI_INTEGER8, MPI_SUM, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_Exscan failed"
    end if

    ! MPI_Exscan receive value is undefined on rank zero.
    if (rank == 0) rank_disp = 0_int64


    ! Convert domain offsets from rank-local stream offsets to
    ! absolute byte offsets in the shared checkpoint file.

    off_loc = CP_DATA_POS + rank_disp + off_loc


    ! ---------------------------------------------------------------
    ! Construct global checkpoint directory
    ! ---------------------------------------------------------------

    call build_checkpoint_directory( &
         gid_loc, load_loc, off_loc, len_loc)


    ! build_checkpoint_directory allocates cp_* on every rank but
    ! fills them only on rank zero.  Broadcast the completed
    ! directory so distribute_grid/load_adapt_mpi can use it
    ! immediately without rereading the checkpoint.

    call MPI_Bcast( &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_load broadcast failed"
    end if

    call MPI_Bcast( &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_offset broadcast failed"
    end if

    call MPI_Bcast( &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: cp_nbytes broadcast failed"
    end if


    ! ---------------------------------------------------------------
    ! Copy rank-local scratch stream into memory
    ! ---------------------------------------------------------------

    if (rank_bytes > int(huge(0),int64)) then
       error stop "dump_adapt_mpi: rank checkpoint buffer too large"
    end if

    allocate(buf(int(rank_bytes)))

    if (rank_bytes > 0_int64) then

       flush(fid)

       read(fid, pos=1) buf

    end if

    close(fid)


    ! ---------------------------------------------------------------
    ! Create shared checkpoint file
    ! ---------------------------------------------------------------

    write(cp4,'(i4.4)') id

    filename = &
         trim(run_id) // "_checkpoint_" // cp4 // ".bin"


    call MPI_File_open( &
         comm, &
         filename, &
         MPI_MODE_CREATE + MPI_MODE_WRONLY, &
         MPI_INFO_NULL, &
         fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_open failed"
    end if


    ! Truncate an existing checkpoint with the same name.
    ! MPI_File_set_size is collective.

    call MPI_File_set_size( &
         fh, &
         int(0, MPI_OFFSET_KIND), &
         ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_set_size failed"
    end if


    ! Rank zero writes header and directory.

    if (rank == 0) then
       call write_checkpoint_directory(fh)
    end if


    ! Make sure directory writing has completed before payload writes.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_Barrier failed"
    end if


    ! ---------------------------------------------------------------
    ! Write rank-local payloads collectively
    ! ---------------------------------------------------------------

    call MPI_File_write_at_all( &
         fh, &
         int(CP_DATA_POS + rank_disp, MPI_OFFSET_KIND), &
         buf, &
         int(rank_bytes), &
         MPI_BYTE, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_write_at_all failed"
    end if


    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: MPI_File_close failed"
    end if


    ! ---------------------------------------------------------------
    ! Compress completed checkpoint
    ! ---------------------------------------------------------------
    !
    ! All MPI ranks must have closed the shared file before rank zero
    ! starts compression.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: pre-zstd barrier failed"
    end if


    if (rank == 0) then

       cmd = 'zstd -q -1 --rm -f "' // trim(filename) // '"'

       call execute_command_line( &
            cmd, &
            exitstat=exitstat, &
            cmdstat=cmdstat)

       if (cmdstat /= 0) then
          error stop "dump_adapt_mpi: failed to execute zstd"
       end if

       if (exitstat /= 0) then
          error stop "dump_adapt_mpi: zstd compression failed"
       end if

    end if


    ! Ensure compression has completed before any rank proceeds.
    ! The checkpoint now exists as filename//".zst".

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "dump_adapt_mpi: post-zstd barrier failed"
    end if


    ! ---------------------------------------------------------------
    ! Cleanup
    ! ---------------------------------------------------------------

    deallocate(buf)

    deallocate(gid_loc)
    deallocate(load_loc)
    deallocate(off_loc)
    deallocate(len_loc)

  end subroutine dump_adapt_mpi


  subroutine build_checkpoint_directory (gid_loc, load_loc, off_loc, len_loc)

    ! Gather global-domain IDs, loads, offsets and serialized sizes
    ! onto rank zero and construct checkpoint directory arrays in
    ! global-domain order.
    !
    ! cp_load, cp_offset and cp_nbytes are allocated on ALL ranks
    ! here so that dump_adapt_mpi can broadcast them immediately
    ! after this routine returns.

    implicit none

    integer,        intent(in) :: gid_loc(:)
    integer,        intent(in) :: load_loc(:)
    integer(int64), intent(in) :: off_loc(:)
    integer(int64), intent(in) :: len_loc(:)

    integer :: sz
    integer :: n_tot
    integer :: ierr
    integer :: i, r, gid

    integer :: rcounts(n_process)
    integer :: displs(n_process)

    integer, allocatable :: gid_glo(:)
    integer, allocatable :: load_recv(:)

    integer(int64), allocatable :: off_recv(:)
    integer(int64), allocatable :: len_recv(:)

    logical, allocatable :: seen(:)


    ! ---------------------------------------------------------------
    ! Allocate checkpoint directory arrays on every rank
    ! ---------------------------------------------------------------

    if (allocated(cp_load))   deallocate(cp_load)
    if (allocated(cp_offset)) deallocate(cp_offset)
    if (allocated(cp_nbytes)) deallocate(cp_nbytes)

    allocate(cp_load(N_GLO_DOMAIN))
    allocate(cp_offset(N_GLO_DOMAIN))
    allocate(cp_nbytes(N_GLO_DOMAIN))

    ! Initialize for easier debugging.  Only rank zero fills them
    ! before the subsequent broadcasts.

    cp_load   = 0
    cp_offset = 0_int64
    cp_nbytes = 0_int64


    ! ---------------------------------------------------------------
    ! Gather number of domains contributed by each rank
    ! ---------------------------------------------------------------

    sz = size(gid_loc)

    call MPI_Gather( &
         sz, 1, MPI_INTEGER, &
         rcounts, 1, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "build_checkpoint_directory: MPI_Gather failed"
    end if


    ! ---------------------------------------------------------------
    ! Construct Gatherv receive layout
    ! ---------------------------------------------------------------

    if (rank == 0) then

       displs(1) = 0

       do r = 2, n_process
          displs(r) = displs(r-1) + rcounts(r-1)
       end do

       n_tot = sum(rcounts)

       if (n_tot /= N_GLO_DOMAIN) then
          write(*,'(A,I0,A,I0)') &
               "build_checkpoint_directory: expected ", &
               N_GLO_DOMAIN, &
               " domains, received ", &
               n_tot

          error stop
       end if

       allocate(gid_glo(n_tot))
       allocate(load_recv(n_tot))
       allocate(off_recv(n_tot))
       allocate(len_recv(n_tot))

    else

       ! Receive arguments are ignored on non-root ranks.
       rcounts = 0
       displs  = 0

       allocate(gid_glo(0))
       allocate(load_recv(0))
       allocate(off_recv(0))
       allocate(len_recv(0))

    end if


    ! ---------------------------------------------------------------
    ! Gather global-domain IDs
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         gid_loc, sz, MPI_INTEGER, &
         gid_glo, rcounts, displs, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: gid MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather domain loads
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         load_loc, sz, MPI_INTEGER, &
         load_recv, rcounts, displs, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: load MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather absolute checkpoint offsets
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         off_loc, sz, MPI_INTEGER8, &
         off_recv, rcounts, displs, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: offset MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Gather serialized byte lengths
    ! ---------------------------------------------------------------

    call MPI_Gatherv( &
         len_loc, sz, MPI_INTEGER8, &
         len_recv, rcounts, displs, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "build_checkpoint_directory: length MPI_Gatherv failed"
    end if


    ! ---------------------------------------------------------------
    ! Construct global-domain ordered directory on rank zero
    ! ---------------------------------------------------------------

    if (rank == 0) then

       allocate(seen(N_GLO_DOMAIN))

       seen = .false.

       do i = 1, n_tot

          gid = gid_glo(i)

          if (gid < 0 .or. gid >= N_GLO_DOMAIN) then

             write(*,'(A,I0)') &
                  "build_checkpoint_directory: invalid gid ", gid

             error stop

          end if


          if (seen(gid+1)) then

             write(*,'(A,I0)') &
                  "build_checkpoint_directory: duplicate gid ", gid

             error stop

          end if


          if (len_recv(i) <= 0_int64) then

             write(*,'(A,I0,A,I0)') &
                  "build_checkpoint_directory: domain ", gid, &
                  " has invalid byte count ", len_recv(i)

             error stop

          end if


          cp_load(gid+1)   = load_recv(i)
          cp_offset(gid+1) = off_recv(i)
          cp_nbytes(gid+1) = len_recv(i)

          seen(gid+1) = .true.

       end do


       if (.not. all(seen)) then
          error stop &
               "build_checkpoint_directory: one or more domains missing"
       end if

       deallocate(seen)

    end if


    ! ---------------------------------------------------------------
    ! Cleanup temporary gather arrays
    ! ---------------------------------------------------------------

    deallocate(gid_glo)
    deallocate(load_recv)
    deallocate(off_recv)
    deallocate(len_recv)

  end subroutine build_checkpoint_directory


  subroutine write_checkpoint_directory (fh)
    implicit none

    type(MPI_File), intent(in) :: fh

    integer :: ierr
    type(MPI_Status) :: status

    integer(int64) :: header(3)

    header(1) = CP_MAGIC
    header(2) = CP_VERSION
    header(3) = int(N_GLO_DOMAIN, int64)

    call MPI_File_write_at( &
         fh, &
         int(0, MPI_OFFSET_KIND), &
         header, 3, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: header write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_LOAD_POS, MPI_OFFSET_KIND), &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: load write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_OFFSET_POS, MPI_OFFSET_KIND), &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: offset write failed"
    end if

    call MPI_File_write_at( &
         fh, &
         int(CP_NBYTES_POS, MPI_OFFSET_KIND), &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         status, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_checkpoint_directory: nbytes write failed"
    end if

  end subroutine write_checkpoint_directory


  subroutine prepare_checkpoint_file (id, filename)
    ! Ensure that the uncompressed checkpoint file exists.
    !
    ! A persistent filename.zst is canonical.  On a new process, rank zero
    ! recreates filename from it even if a stale filename survived a failed
    ! restart.  Repeated calls in one restart retain the prepared temporary
    ! file until load_adapt_mpi completes.
    !
    ! This routine is safe to call more than once for the same restart.

    implicit none

    integer, intent(in) :: id

    character(:), allocatable, intent(out) :: filename

    integer :: ierr
    integer :: cmdstat, exitstat
    integer :: decompress_status

    logical :: bin_exists
    logical :: zst_exists
    logical :: local_exists
    logical :: all_exist

    character(4) :: cp4
    character(:), allocatable :: zst_filename
    character(:), allocatable :: cmd


    write(cp4,'(i4.4)') id

    filename = &
         trim(run_id) // "_checkpoint_" // cp4 // ".bin"

    zst_filename = filename // ".zst"

    decompress_status = 0


    ! ---------------------------------------------------------------
    ! Rank zero ensures that the uncompressed checkpoint exists.
    ! ---------------------------------------------------------------

    if (rank == 0) then

       inquire(file=filename,     exist=bin_exists)
       inquire(file=zst_filename, exist=zst_exists)

       if (zst_exists .and. &
            (.not. bin_exists .or. .not. cp_temporary_bin)) then

          cmd = &
               'zstd -q -d -k -f "' // trim(zst_filename) // '"'

          call execute_command_line( &
               cmd, &
               exitstat=exitstat, &
               cmdstat=cmdstat)

          if (cmdstat /= 0) then

             decompress_status = 2

          elseif (exitstat /= 0) then

             decompress_status = 3

          else

             inquire(file=filename, exist=bin_exists)

             if (.not. bin_exists) then
                decompress_status = 4
             else
                cp_temporary_bin = .true.
             end if

          end if

       elseif (.not. bin_exists) then

          if (.not. zst_exists) then

             decompress_status = 1
          end if

       end if

    end if


    ! ---------------------------------------------------------------
    ! Broadcast preparation state.
    ! ---------------------------------------------------------------

    call MPI_Bcast( &
         decompress_status, 1, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "prepare_checkpoint_file: status broadcast failed"
    end if

    if (decompress_status /= 0) then

       select case (decompress_status)

       case (1)
          error stop &
               "prepare_checkpoint_file: checkpoint .bin and .bin.zst not found"

       case (2)
          error stop &
               "prepare_checkpoint_file: failed to execute zstd"

       case (3)
          error stop &
               "prepare_checkpoint_file: zstd decompression failed"

       case (4)
          error stop &
               "prepare_checkpoint_file: decompressed checkpoint not found"

       case default
          error stop &
               "prepare_checkpoint_file: checkpoint preparation failed"

       end select

    end if

    call MPI_Bcast( &
         cp_temporary_bin, 1, MPI_LOGICAL, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "prepare_checkpoint_file: temporary flag broadcast failed"
    end if


    ! Rank zero has completed decompression before this collective.
    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "prepare_checkpoint_file: pre-open barrier failed"
    end if


    ! Verify that every rank can see the uncompressed file before
    ! returning to a collective MPI_File_open.
    inquire(file=filename, exist=local_exists)

    call MPI_Allreduce( &
         local_exists, all_exist, 1, MPI_LOGICAL, MPI_LAND, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop &
            "prepare_checkpoint_file: visibility reduction failed"
    end if

    if (.not. all_exist) then
       if (.not. local_exists) then
          write(*,'(A,I0,A,A,A)') &
               'rank ', rank, ': checkpoint not visible: "', &
               trim(filename), '"'
       end if
       error stop &
            "prepare_checkpoint_file: checkpoint .bin not visible on all ranks"
    end if

  end subroutine prepare_checkpoint_file


  subroutine read_checkpoint_directory (id)
    implicit none

    integer, intent(in) :: id

    integer :: ierr
    integer :: error_len

    character(len=MPI_MAX_ERROR_STRING) :: error_string

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: header(3)

    character(:), allocatable :: filename


    ! Ensure that the uncompressed checkpoint exists.  If only the
    ! .bin.zst file exists, rank zero creates a temporary .bin file.

    call prepare_checkpoint_file(id, filename)


    ! ---------------------------------------------------------------
    ! Open checkpoint
    ! ---------------------------------------------------------------

    call MPI_File_open( &
         comm, filename, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       call MPI_Error_string(ierr, error_string, error_len)
       write(*,'(A,I0,A,A)') &
            'rank ', rank, ': read_checkpoint_directory MPI_File_open: ', &
            error_string(1:error_len)
       write(*,'(A,I0,A,A,A)') &
            'rank ', rank, ': filename = "', trim(filename), '"'
       error stop "read_checkpoint_directory: MPI_File_open failed"
    end if


    ! ---------------------------------------------------------------
    ! Allocate checkpoint directory
    ! ---------------------------------------------------------------

    if (allocated(cp_load))   deallocate(cp_load)
    if (allocated(cp_offset)) deallocate(cp_offset)
    if (allocated(cp_nbytes)) deallocate(cp_nbytes)

    allocate(cp_load(N_GLO_DOMAIN))
    allocate(cp_offset(N_GLO_DOMAIN))
    allocate(cp_nbytes(N_GLO_DOMAIN))


    ! ---------------------------------------------------------------
    ! Rank zero reads checkpoint header and directory
    ! ---------------------------------------------------------------

    if (rank == 0) then

       call MPI_File_read_at( &
            fh, &
            int(0, MPI_OFFSET_KIND), &
            header, 3, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: header read failed"
       end if
       call require_mpi_transfer_count( &
            status,MPI_INTEGER8,3, &
            "read_checkpoint_directory header")

       if (header(1) /= CP_MAGIC) then
          error stop "read_checkpoint_directory: bad checkpoint magic"
       end if

       if (header(2) /= CP_VERSION) then
          error stop "read_checkpoint_directory: unsupported version"
       end if

       if (header(3) /= int(N_GLO_DOMAIN, int64)) then
          error stop "read_checkpoint_directory: wrong domain count"
       end if

       call MPI_File_read_at( &
            fh, &
            int(CP_LOAD_POS, MPI_OFFSET_KIND), &
            cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: load read failed"
       end if
       call require_mpi_transfer_count( &
            status,MPI_INTEGER,N_GLO_DOMAIN, &
            "read_checkpoint_directory load")

       call MPI_File_read_at( &
            fh, &
            int(CP_OFFSET_POS, MPI_OFFSET_KIND), &
            cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: offset read failed"
       end if
       call require_mpi_transfer_count( &
            status,MPI_INTEGER8,N_GLO_DOMAIN, &
            "read_checkpoint_directory offset")

       call MPI_File_read_at( &
            fh, &
            int(CP_NBYTES_POS, MPI_OFFSET_KIND), &
            cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "read_checkpoint_directory: nbytes read failed"
       end if
       call require_mpi_transfer_count( &
            status,MPI_INTEGER8,N_GLO_DOMAIN, &
            "read_checkpoint_directory nbytes")

       call validate_checkpoint_file_extent( &
            fh,"read_checkpoint_directory")

    end if


    ! ---------------------------------------------------------------
    ! Close checkpoint
    ! ---------------------------------------------------------------

    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: MPI_File_close failed"
    end if


    ! ---------------------------------------------------------------
    ! Broadcast checkpoint directory
    ! ---------------------------------------------------------------

    call MPI_Bcast( &
         cp_load, N_GLO_DOMAIN, MPI_INTEGER, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: cp_load broadcast failed"
    end if

    call MPI_Bcast( &
         cp_offset, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: cp_offset broadcast failed"
    end if

    call MPI_Bcast( &
         cp_nbytes, N_GLO_DOMAIN, MPI_INTEGER8, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "read_checkpoint_directory: cp_nbytes broadcast failed"
    end if

  end subroutine read_checkpoint_directory


  subroutine dump_domain (fid, d)
    implicit none

    integer, intent(in) :: fid
    integer, intent(in) :: d

    integer :: c, ibeg, iend, j, k
    integer :: l, p_chd, p_lev, p_par, v

    logical, dimension(1:N_CHDRN) :: required

    ! Checkpoint-global/model state.

    write(fid) istep_cumul
    write(fid) time
    write(fid) itime
    write(fid) iwrite
    write(fid) threshold

    call dump(fid)

    ! Coarsest-scale solution.
    call apply_to_pole_d ( &
         write_scalar, grid(d), min_level-1, z_null, fid, .true.)

    p_par = 1

    do k = zmin, zmax
       do v = 1, N_VARIABLE

          ibeg = MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1
          iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1

          write(fid) sol(v,k)%data(d)%elts(ibeg:iend)

       end do
    end do

    if (vert_diffuse) then

       do k = 1, zlevels

          ibeg = grid(d)%patch%elts(p_par+1)%elts_start + 1
          iend = ibeg + PATCH_SIZE**2 - 1

          write(fid) tke(k)%data(d)%elts(ibeg:iend)

       end do

    end if

    ! Finer scales.
    do l = min_level, level_end

       p_lev = 0

       do j = 1, grid(d)%lev(l)%length

          p_par = grid(d)%lev(l)%elts(j)

          ! Do not save children of a deleted parent.
          if (grid(d)%patch%elts(p_par+1)%deleted) then

             do c = 1, N_CHDRN

                p_chd = grid(d)%patch%elts(p_par+1)%children(c)

                if (p_chd > 0) then
                   grid(d)%patch%elts(p_chd+1)%deleted = .true.
                end if

             end do

             cycle

          end if

          ! Wavelet coefficients.
          do k = zmin, zmax
             do v = 1, N_VARIABLE

                ibeg = &
                     MULT(v)*grid(d)%patch%elts(p_par+1)%elts_start + 1

                iend = ibeg + MULT(v)*PATCH_SIZE**2 - 1

                write(fid) &
                     wav_coeff(v,k)%data(d)%elts(ibeg:iend)

             end do
          end do

          if (vert_diffuse) then

             do k = 1, zlevels

                ibeg = &
                     grid(d)%patch%elts(p_par+1)%elts_start + 1

                iend = ibeg + PATCH_SIZE**2 - 1

                write(fid) wav_tke(k)%data(d)%elts(ibeg:iend)

             end do

          end if

          ! Record which children are required.
          do c = 1, N_CHDRN

             p_chd = grid(d)%patch%elts(p_par+1)%children(c)

             if (p_chd > 0) then

                required(c) = &
                     check_child_required(grid(d), p_par, c-1)

                grid(d)%patch%elts(p_chd+1)%deleted = &
                     .not. required(c)

                if (required(c) .and. &
                     p_lev+1 <= size(grid(d)%lev(l+1)%elts)) then

                   p_lev = p_lev + 1
                   grid(d)%lev(l+1)%elts(p_lev) = p_chd

                end if

             else

                required(c) = .false.

             end if

          end do

          ! This used to be written to fid_gr.
          ! It now immediately follows this patch's field data.
          write(fid) required

       end do

       if (l+1 <= max_level) then
          grid(d)%lev(l+1)%length = p_lev
       end if

    end do

  end subroutine dump_domain


  subroutine write_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: fid, i, j, p, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, k, v

    d = dom%id+1
    id = idx(i, j, offs, dims) + 1

    do k = zmin, zmax
       do v = scalars(1), scalars(2)
          write (fid) sol(v,k)%data(d)%elts(id)
       end do
    end do

    if (vert_diffuse) then
       do k = 1, zlevels
          write (fid) tke(k)%data(d)%elts(id)
       end do
    end if
  end subroutine write_scalar


  subroutine load_adapt_mpi (id)
    ! Read adaptive checkpoint data from one shared MPI-IO file.

    implicit none

    integer, intent(in) :: id

    integer :: d, gid
    integer :: fid
    integer :: ierr
    integer :: error_len

    integer :: delete_unit

    character(len=MPI_MAX_ERROR_STRING) :: error_string

    type(MPI_File)   :: fh
    type(MPI_Status) :: status

    integer(int64) :: nbytes
    integer(int64) :: p

    integer(int64), allocatable :: domain_pos(:)
    integer(int8),  allocatable :: buf(:)

    character(:), allocatable :: filename


    ! ---------------------------------------------------------------
    ! Validate checkpoint directory
    ! ---------------------------------------------------------------

    if (.not. allocated(cp_offset)) then
       error stop "load_adapt_mpi: cp_offset is not allocated"
    end if

    if (.not. allocated(cp_nbytes)) then
       error stop "load_adapt_mpi: cp_nbytes is not allocated"
    end if

    if (size(cp_offset) /= N_GLO_DOMAIN .or. &
         size(cp_nbytes) /= N_GLO_DOMAIN) then

       error stop "load_adapt_mpi: invalid checkpoint directory"

    end if


    ! ---------------------------------------------------------------
    ! Ensure the uncompressed checkpoint exists
    ! ---------------------------------------------------------------
    !
    ! This is required even when read_checkpoint_directory was not
    ! called.  In particular, dump_adapt_mpi leaves the checkpoint
    ! directory in memory but compresses the shared .bin to .bin.zst
    ! before load_adapt_mpi may be called.

    call prepare_checkpoint_file(id, filename)


    ! ---------------------------------------------------------------
    ! Open shared checkpoint file
    ! ---------------------------------------------------------------
   
    call MPI_File_open( &
         comm, filename, MPI_MODE_RDONLY, &
         MPI_INFO_NULL, fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       call MPI_Error_string(ierr, error_string, error_len)
       write(*,'(A,I0,A,A)') &
            'rank ', rank, ': load_adapt_mpi MPI_File_open: ', &
            error_string(1:error_len)
       write(*,'(A,I0,A,A,A)') &
            'rank ', rank, ': filename = "', trim(filename), '"'
       error stop "load_adapt_mpi: MPI_File_open failed"
    end if

    call validate_checkpoint_file_extent(fh,"load_adapt_mpi")


    ! ---------------------------------------------------------------
    ! Create rank-local scratch stream
    ! ---------------------------------------------------------------

    open( &
         newunit=fid, &
         status="scratch", &
         form="unformatted", &
         access="stream", &
         action="readwrite")


    allocate(domain_pos(size(grid)))


    ! ---------------------------------------------------------------
    ! Read locally owned domains from shared checkpoint
    ! ---------------------------------------------------------------

    do d = 1, size(grid)

       gid = glo_id(rank+1,d)


       if (gid < 0 .or. gid >= N_GLO_DOMAIN) then
          error stop "load_adapt_mpi: invalid global domain ID"
       end if


       nbytes = cp_nbytes(gid+1)


       if (nbytes <= 0_int64) then
          error stop "load_adapt_mpi: invalid domain record size"
       end if


       ! buf() has default-integer bounds, so guard the conversion.
       if (nbytes > int(huge(0),int64)) then
          error stop "load_adapt_mpi: domain record too large"
       end if


       allocate(buf(int(nbytes)))


       call MPI_File_read_at( &
            fh, &
            int(cp_offset(gid+1), MPI_OFFSET_KIND), &
            buf, int(nbytes), MPI_BYTE, &
            status, ierr)

       if (ierr /= MPI_SUCCESS) then
          error stop "load_adapt_mpi: domain read failed"
       end if
       call require_mpi_transfer_count( &
            status,MPI_BYTE,int(nbytes), &
            "load_adapt_mpi domain payload")


       inquire(unit=fid, pos=p)

       domain_pos(d) = p


       write(fid) buf


       deallocate(buf)

    end do


    ! ---------------------------------------------------------------
    ! Close shared checkpoint
    ! ---------------------------------------------------------------

    call MPI_File_close(fh, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: MPI_File_close failed"
    end if


    ! ---------------------------------------------------------------
    ! Reconstruct local domains
    ! ---------------------------------------------------------------

    flush(fid)

    call load_domains_stream(fid, domain_pos)

    close(fid)

    deallocate(domain_pos)


    ! ---------------------------------------------------------------
    ! Remove temporary uncompressed checkpoint
    ! ---------------------------------------------------------------
    !
    ! All ranks must have completely finished reading and reconstructing
    ! their domains before rank zero removes the shared .bin file.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: pre-delete barrier failed"
    end if


    if (rank == 0 .and. cp_temporary_bin) then

       ! Delete using Fortran rather than invoking an external rm command.

       open( &
            newunit=delete_unit, &
            file=filename, &
            status='old', &
            action='read', &
            iostat=ierr)

       if (ierr /= 0) then
          error stop "load_adapt_mpi: unable to open temporary checkpoint for deletion"
       end if

       close(delete_unit, status='delete', iostat=ierr)

       if (ierr /= 0) then
          error stop "load_adapt_mpi: unable to delete temporary checkpoint"
       end if

       cp_temporary_bin = .false.

    end if


    ! Ensure deletion is complete before returning.

    call MPI_Barrier(comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "load_adapt_mpi: post-delete barrier failed"
    end if

    ! Keep the module state consistent on all ranks.
    cp_temporary_bin = .false.


  end subroutine load_adapt_mpi

  
  subroutine load_domains_stream (fid, domain_pos)
    ! Reconstruct all locally owned domains from a single
    ! rank-local unformatted stream.
    !
    ! domain_pos(d) is the current position in the serialized
    ! record for local domain d.

    implicit none

    integer,        intent(in)    :: fid
    integer(int64), intent(inout) :: domain_pos(:)

    integer :: c, d
    integer :: ibeg, iend
    integer :: j, k, l
    integer :: old_n_patch
    integer :: p_chd, p_par
    integer :: v

    logical :: first_read

    logical, dimension(N_CHDRN) :: required

    if (size(domain_pos) /= size(grid)) then
       error stop "load_domains_stream: domain_pos has wrong size"
    end if

    !
    ! ------------------------------------------------------------
    ! Coarsest-scale solution
    ! ------------------------------------------------------------
    !
    do d = 1, size(grid)

       !
       ! The first read explicitly seeks to this domain's record.
       ! All subsequent reads for this domain are sequential.
       !
       read(fid, pos=domain_pos(d)) istep_cumul

       read(fid) time
       read(fid) itime
       read(fid) iwrite
       read(fid) threshold

       !
       ! Test-case-specific data. This may be empty.
       !
       call load(fid)

       !
       ! Pole scalar values.
       !
       call apply_to_pole_d ( &
            read_scalar, &
            grid(d), &
            min_level-1, &
            z_null, &
            fid, &
            .true.)

       !
       ! Coarse parent patch.
       !
       p_par = 1

       do k = zmin, zmax
          do v = 1, N_VARIABLE

             ibeg = &
                  MULT(v) * &
                  grid(d)%patch%elts(p_par+1)%elts_start + 1

             iend = &
                  ibeg + MULT(v)*PATCH_SIZE**2 - 1

             read(fid) &
                  sol(v,k)%data(d)%elts(ibeg:iend)

          end do
       end do

       if (vert_diffuse) then

          do k = 1, zlevels

             ibeg = &
                  grid(d)%patch%elts(p_par+1)%elts_start + 1

             iend = ibeg + PATCH_SIZE**2 - 1

             read(fid) &
                  tke(k)%data(d)%elts(ibeg:iend)

          end do

       end if

       !
       ! Save the position of this domain's first finer-scale
       ! wavelet record.
       !
       inquire(unit=fid, pos=domain_pos(d))

    end do

    !
    ! ------------------------------------------------------------
    ! Finer-scale wavelets
    ! ------------------------------------------------------------
    !
    ! Preserve the original level-by-level reconstruction:
    !
    !     process all domains at level l
    !     call post_refine
    !     proceed to newly created level
    !
    l = 1

    do while (level_end > l)

       l = level_end

       if (rank == 0) then
          write(6,'(a,i2)') 'Loading level ', l
       end if

       do d = 1, size(grid)

          old_n_patch = grid(d)%patch%length

          !
          ! The first actual read for this domain at this level
          ! must seek back to domain_pos(d). Thereafter all reads
          ! are sequential.
          !
          first_read = .true.

          do j = 1, grid(d)%lev(l)%length

             p_par = grid(d)%lev(l)%elts(j)

             !
             ! Wavelet coefficients.
             !
             do k = zmin, zmax
                do v = 1, N_VARIABLE

                   ibeg = &
                        MULT(v) * &
                        grid(d)%patch%elts(p_par+1)%elts_start + 1

                   iend = &
                        ibeg + MULT(v)*PATCH_SIZE**2 - 1

                   if (first_read) then

                      read(fid, pos=domain_pos(d)) &
                           wav_coeff(v,k)%data(d)%elts(ibeg:iend)

                      first_read = .false.

                   else

                      read(fid) &
                           wav_coeff(v,k)%data(d)%elts(ibeg:iend)

                   end if

                end do
             end do

             !
             ! Turbulent kinetic energy wavelets.
             !
             if (vert_diffuse) then

                do k = 1, zlevels

                   ibeg = &
                        grid(d)%patch%elts(p_par+1)%elts_start + 1

                   iend = ibeg + PATCH_SIZE**2 - 1

                   !
                   ! Normally first_read is already false because
                   ! N_VARIABLE > 0. This branch keeps the stream
                   ! handling correct even if that ever changes.
                   !
                   if (first_read) then

                      read(fid, pos=domain_pos(d)) &
                           wav_tke(k)%data(d)%elts(ibeg:iend)

                      first_read = .false.

                   else

                      read(fid) &
                           wav_tke(k)%data(d)%elts(ibeg:iend)

                   end if

                end do

             end if

             !
             ! Child-presence information follows this patch's
             ! coefficient data in the stage-2/3 stream.
             !
             if (first_read) then

                read(fid, pos=domain_pos(d)) required
                first_read = .false.

             else

                read(fid) required

             end if

             do c = 1, N_CHDRN

                if (required(c)) then
                   call refine_patch1(grid(d), p_par, c-1)
                end if

             end do

          end do

          !
          ! Only update the domain position if something was
          ! actually read at this level.
          !
          ! This is important for a domain with zero patches at
          ! the current level: its saved stream position must not
          ! accidentally become the position of another domain.
          !
          if (.not. first_read) then
             inquire(unit=fid, pos=domain_pos(d))
          end if

          !
          ! Complete child patch construction exactly as in the
          ! original loader.
          !
          do p_par = 2, old_n_patch

             do c = 1, N_CHDRN

                p_chd = &
                     grid(d)%patch%elts(p_par)%children(c)

                if (p_chd+1 > old_n_patch) then

                   call refine_patch2( &
                        grid(d), &
                        p_par-1, &
                        c-1)

                end if

             end do

          end do

       end do

       !
       ! Must remain outside the domain loop.
       !
       call post_refine

    end do

    !
    ! ------------------------------------------------------------
    ! Boundary state
    ! ------------------------------------------------------------
    !
    sol%bdry_uptodate       = .false.
    wav_coeff%bdry_uptodate = .false.

    if (vert_diffuse) then
       tke%bdry_uptodate     = .false.
       wav_tke%bdry_uptodate = .false.
    end if

  end subroutine load_domains_stream


  subroutine read_scalar (dom, p, i, j, zlev, offs, dims, fid)
    ! For poles

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: fid, i, j, p, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, k, v

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    do k = zmin, zmax
       do v = scalars(1), scalars(2)
          read (fid) sol(v,k)%data(d)%elts(id)
       end do
    end do

    if (vert_diffuse) then
       do k = 1, zlevels
          read (fid) tke(k)%data(d)%elts(id)
       end do
    end if
  end subroutine read_scalar


end module checkpoint_mod
