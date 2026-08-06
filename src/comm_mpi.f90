module comm_mpi_mod

  use mpi_f08

  use kind_mod,   only : dp
  use shared_mod, only : AT_NODE, AT_EDGE, Coord, EDGE, N_BDRY, N_GLO_DOMAIN, NONE, POLE, level_start, level_end, run_id, eps, &
       n_domain

  use arch_mod,       only : abort_run, barrier, comm, MPI_IN, MPI_DP, MPI_SP, glo_id, n_process, owner, rank
  use dyn_arrays,     only : Int_Array, Float_Array, append, extend, init
  use domain_mod,     only : Domain, Float_Field, grid
  use domain_ops_mod, only : apply_onescale, apply_onescale__int, apply_to_pole, sub8

  use comm_mod,       only : coord_get, coord_set, init_comm, get9, set9, sync_val, rot_direction, domain_load, &
       comm_communication, comm_masks, comm_nodes3, comm_nodes9, comm_patch_conn, cp_bdry_inside, unpack, unpack_comm_struct

  implicit none

  private
  public :: cal_load_balance, check_alltoall_lengths, init_comm_mpi, init_comm_mpi_mod
  public :: get_timing, start_timing, stop_timing
  public :: comm_communication_mpi, comm_masks_mpi, comm_nodes3_mpi, comm_nodes9_mpi, comm_patch_conn_mpi
  public :: recv_lengths, send_lengths,recv_offsets, send_offsets, req
  public :: gather_int, gather_vec, sum_int, sync_array, sync_max_real, sync_max_int, sync_min_int, sync_min_real, sum_real
  public :: update_bdry, update_bdry1, update_bdry__finish, update_bdry__start
  public :: write_level_mpi, write_load_conn

  integer                         :: nreq
  type(Int_Array)                 :: recv_buf_i, send_buf_i
  type(Float_Array)               :: recv_buf,   send_buf
  type(MPI_Request),  allocatable :: req(:)
  integer,            allocatable :: recv_lengths(:), recv_offsets(:), send_lengths(:), send_offsets(:)
  real(dp), dimension(2)          :: times

  integer,            parameter   :: TAG_BDRY_S = 1100
  integer,            parameter   :: TAG_BDRY_V = 1101
  integer,            parameter   :: TAG_BDRY_A = 1102
  logical,            parameter   :: deadlock = .false. ! test for communication deadlocks (use for debug only as it is slow)

  interface sum_real
     procedure :: sum_real_0, sum_real_1
  end interface sum_real

  interface sum_int
     procedure :: sum_int_0, sum_int_1
  end interface sum_int

  interface sync_max_real
     procedure :: sync_max_real_0, sync_max_real_1
  end interface sync_max_real

  interface sync_min_real
     procedure :: sync_min_real_0, sync_min_real_1
  end interface sync_min_real

  interface update_bdry
     procedure :: update_bdry_0, update_bdry_1, update_bdry_2
  end interface update_bdry

  interface update_bdry1
     procedure :: update_bdry1_0, update_bdry1_1, update_bdry1_2
  end interface update_bdry1

  interface update_bdry__start
     procedure :: update_bdry__start_0, update_bdry__start_1, update_bdry__start_2
  end interface update_bdry__start

  interface update_bdry__start1 
     procedure :: update_bdry__start1_0, update_bdry__start1_1, update_bdry__start1_2
  end interface update_bdry__start1

  interface update_bdry__finish
     procedure :: update_bdry__finish_0, update_bdry__finish_1, update_bdry__finish_2
  end interface update_bdry__finish

  interface update_bdry__finish1 
     procedure :: update_bdry__finish1_0, update_bdry__finish1_1, update_bdry__finish1_2
  end interface update_bdry__finish1

  interface gather_vec
     procedure :: gatherv_int, gatherv_real4, gatherv_real8
  end interface gather_vec

contains

  
  subroutine init_comm_mpi
    implicit none
    allocate (recv_lengths(n_process), recv_offsets(n_process))
    allocate (send_lengths(n_process), send_offsets(n_process))
    allocate (req(2*n_process))
    recv_lengths = 0
    recv_offsets = 0
    send_lengths = 0
    send_offsets = 0
    req = MPI_REQUEST_NULL  
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  
  subroutine write_load_conn (id)
    implicit none

    integer, intent(in) :: id

    integer, parameter :: fid = 599

    integer :: d, ii, r, sz, n_tot
    integer :: ierr, io_stat
    integer :: displs(n_process)
    integer :: rcounts(n_process)

    integer, allocatable :: n_active_glo(:)
    integer, allocatable :: n_active_loc(:)

    character(len=255) :: filename
    character(len=256) :: io_msg

    ! Number of load entries contributed by this MPI rank.
    sz = size(grid)*(N_GLO_DOMAIN+1)

    ! Only the root process uses rcounts after this collective.
    call MPI_Gather (sz, 1, MPI_IN, rcounts, 1, MPI_IN, 0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_load_conn: MPI_Gather failed"
    end if

    ! Construct local load data.
    allocate(n_active_loc(sz))
    n_active_loc = 0

    ii = 1

    do d = 1, size(grid)
       n_active_loc(ii) = domain_load(grid(d))

       n_active_loc(ii+1:ii+N_GLO_DOMAIN) = ( &
            grid(d)%pack(AT_NODE,:)%length   &
            + grid(d)%pack(AT_EDGE,:)%length &
            + grid(d)%unpk(AT_NODE,:)%length &
            + grid(d)%unpk(AT_EDGE,:)%length) / 2

       ii = ii + N_GLO_DOMAIN + 1
    end do

    ! MPI_Gatherv uses rcounts, displs, and the receive buffer only
    ! on the root process.
    n_tot = 0
    if (rank == 0) then
       displs(1) = 0

       do r = 2, n_process
          displs(r) = displs(r-1) + rcounts(r-1)
       end do

       n_tot = sum (rcounts)
       allocate(n_active_glo(n_tot))
    else
       displs = 0
       rcounts = 0
       allocate(n_active_glo(0))
    end if

    call MPI_Gatherv ( &
         n_active_loc, sz, MPI_IN, &
         n_active_glo, rcounts, displs, MPI_IN, &
         0, comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "write_load_conn: MPI_Gatherv failed"
    end if

    deallocate(n_active_loc)

    if (rank == 0) then
       ! Each global domain should contribute one block containing:
       !
       !   domain load + N_GLO_DOMAIN connection loads
       if (n_tot /= N_GLO_DOMAIN*(N_GLO_DOMAIN+1)) then
          error stop "write_load_conn: unexpected gathered data size"
       end if

       write(filename, '(A,A,I4.4)') trim(run_id), "_conn.", id

       open( &
            unit=fid, file=trim(filename), status="replace", &
            action="write", form="formatted", &
            iostat=io_stat, iomsg=io_msg)

       if (io_stat /= 0) then
          write(*,'(A)') &
               "write_load_conn: unable to open "//trim(filename)
          write(*,'(A)') trim(io_msg)
          error stop
       end if

       ii = 1

       do d = 1, N_GLO_DOMAIN
          write( &
               fid, '(I10,*(1X,I8))', &
               iostat=io_stat, iomsg=io_msg) &
               n_active_glo(ii:ii+N_GLO_DOMAIN)

          if (io_stat /= 0) then
             write(*,'(A)') &
                  "write_load_conn: error writing "//trim(filename)
             write(*,'(A)') trim(io_msg)
             close(fid)
             error stop
          end if

          ii = ii + N_GLO_DOMAIN + 1
       end do

       close(fid, iostat=io_stat, iomsg=io_msg)

       if (io_stat /= 0) then
          write(*,'(A)') &
               "write_load_conn: error closing "//trim(filename)
          write(*,'(A)') trim(io_msg)
          error stop
       end if
    end if

    deallocate (n_active_glo)
  end subroutine write_load_conn

  
  subroutine cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    ! Compute the minimum, average, and maximum processor loads,
    ! together with the maximum-to-average load ratio.
    implicit none

    integer,  intent(out) :: min_load, max_load
    real(dp), intent(out) :: avg_load, rel_imbalance

    call get_load_balance(min_load, avg_load, max_load)

    if (avg_load > 0.0_dp) then
       rel_imbalance = real(max_load, kind=dp)/avg_load
    else
       rel_imbalance = 1.0_dp
    end if
  end subroutine cal_load_balance

  
  subroutine get_load_balance (mini, avg, maxi)
    implicit none

    integer,  intent(out) :: mini, maxi
    real(dp), intent(out) :: avg

    integer :: d
    integer :: load, load_sum

    load = 0

    do d = 1, size(grid)
       load = load + domain_load(grid(d))
    end do

    call MPI_Allreduce (load, maxi,     1, MPI_IN, MPI_MAX, comm)
    call MPI_Allreduce (load, mini,     1, MPI_IN, MPI_MIN, comm)
    call MPI_Allreduce (load, load_sum, 1, MPI_IN, MPI_SUM, comm)

    avg = real (load_sum, kind=dp)/real(n_process, kind=dp)
  end subroutine get_load_balance

  
  subroutine write_level_mpi (routine, l, zlev, eval_pole, filename)
    implicit none

    procedure(sub8)              :: routine
    integer,          intent(in) :: l, zlev
    logical,          intent(in) :: eval_pole
    character(len=*), intent(in) :: filename

    integer, parameter :: funit = 300

    integer :: r
    integer :: io_stat
    character(len=256) :: io_msg

    do r = 1, n_process
       if (r == rank+1) then
          ! The first process replaces any existing file. Subsequent
          ! processes append their records after the preceding process
          ! has closed the file.
          if (r == 1) then
             open( &
                  unit=funit, file=trim(filename), &
                  form="unformatted", status="replace", action="write", &
                  iostat=io_stat, iomsg=io_msg)
          else
             open( &
                  unit=funit, file=trim(filename), &
                  form="unformatted", status="old", position="append", &
                  action="write", iostat=io_stat, iomsg=io_msg)
          end if

          if (io_stat /= 0) then
             write(*,'(A,I0,A)') &
                  "write_level_mpi: rank ", rank, &
                  " could not open "//trim(filename)
             write(*,'(A)') trim(io_msg)
             error stop
          end if

          if (eval_pole) then
             call apply_to_pole (routine, l, zlev, funit, .true.)
          end if

          call apply_onescale__int (routine, l, zlev, 0, 0, funit)

          close(funit, iostat=io_stat, iomsg=io_msg)

          if (io_stat /= 0) then
             write(*,'(A,I0,A)') &
                  "write_level_mpi: rank ", rank, &
                  " could not close "//trim(filename)
             write(*,'(A)') trim(io_msg)
             error stop
          end if
       end if

       ! Every rank reaches exactly one barrier for each value of r.
       ! This serializes access to the shared output file.
       call barrier
    end do
  end subroutine write_level_mpi

  
  subroutine init_comm_mpi_mod()
    implicit none

    call init(send_buf_i, 0)
    call init(recv_buf_i, 0)
    call init(send_buf,   0)
    call init(recv_buf,   0)
  end subroutine init_comm_mpi_mod


  subroutine comm_communication_mpi()
    implicit none

    call alltoall_dom(unpack_comm_struct, 4)
    call comm_communication()
    call recreate_send_patch_lists()
  end subroutine comm_communication_mpi


  subroutine recreate_send_patch_lists()
    implicit none

    integer :: d, d_neigh
    integer :: i, k, l, n, p, r, s, typ

    do l = level_start, level_end
       do d = 1, size(grid)

          do r = 1, n_process
             grid(d)%src_patch(r,l)%length = 0
          end do

          do k = 1, grid(d)%lev(l)%length
             p = grid(d)%lev(l)%elts(k)

             do s = 1, N_BDRY
                n = grid(d)%patch%elts(p+1)%neigh(s)

                if (n >= 0) cycle

                n   = -n
                typ = grid(d)%bdry_patch%elts(n+1)%side

                if (typ < 1) cycle

                d_neigh = grid(d)%neigh(typ)

                select case (d_neigh)
                case (POLE)
                   do i = 1, 2
                      call handle_neigh( &
                           grid(d), grid(d)%neigh_over_pole(i))
                   end do

                case (NONE)
                   cycle

                case default
                   call handle_neigh(grid(d), d_neigh)
                end select
             end do
          end do
       end do
    end do

  contains

    subroutine handle_neigh (dom, d0)
      implicit none

      type(Domain), intent(inout) :: dom
      integer,      intent(in)    :: d0

      integer :: r0

      r0 = owner(d0+1) + 1

      if (r0 == rank+1) return

      if (dom%src_patch(r0,l)%length > 0) then
         if (dom%src_patch(r0,l)%elts( &
              dom%src_patch(r0,l)%length) == p) return
      end if

      call append(dom%src_patch(r0,l), p)
    end subroutine handle_neigh

  end subroutine recreate_send_patch_lists

  
  subroutine check_alltoall_lengths
    implicit none

    integer :: ierr
    integer :: test_recv_len(n_process)

    call MPI_Alltoall ( &
         send_lengths, 1, MPI_INTEGER, &
         test_recv_len, 1, MPI_INTEGER, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "check_alltoall_lengths: MPI_Alltoall failed"
    end if

    write(3000+rank,*) test_recv_len - recv_lengths
    close(3000+rank)

    call MPI_Barrier (comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "check_alltoall_lengths: MPI_Barrier failed"
    end if
  end subroutine check_alltoall_lengths

  
  subroutine alltoall_dom (unpack_rout, n_values)
    implicit none

    procedure(unpack)     :: unpack_rout
    integer, intent(in)   :: n_values

    integer :: d_dest, d_src
    integer :: dest, src
    integer :: i, k, n_record, r_dest, r_src
    integer, allocatable :: values(:)

    if (n_values /= 4) then
       error stop "alltoall_dom: unpack interface requires 4 integers"
    end if

    allocate(values(n_values))

    send_buf_i%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf_i%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             src = d_src

             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                n_record = grid(src)%send_conn(dest)%length

                if (mod(n_record, n_values) /= 0) then
                   error stop &
                        "alltoall_dom: send_conn length not divisible by record size"
                end if

                call append (send_buf_i, n_record)

                do i = 1, n_record
                   call append (send_buf_i, grid(src)%send_conn(dest)%elts(i))
                end do

                grid(src)%send_conn(dest)%length = 0
             end do
          end do
       end if

       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    call alltoall

    i = 1

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src)

          do d_dest = 1, n_domain(rank+1)
             if (i > recv_buf_i%length) then
                error stop "alltoall_dom: missing record length"
             end if

             n_record = recv_buf_i%elts(i)
             i = i + 1

             if (n_record < 0) then
                error stop "alltoall_dom: negative record length"
             end if

             if (mod(n_record, n_values) /= 0) then
                error stop &
                     "alltoall_dom: received length not divisible by record size"
             end if

             if (i+n_record-1 > recv_buf_i%length) then
                error stop "alltoall_dom: receive buffer truncated"
             end if

             do k = 1, n_record, n_values
                values = recv_buf_i%elts( &
                     i+k-1:i+k+n_values-2)

                call unpack_rout ( &
                     grid(d_dest), src, &
                     values(1), values(2), values(3), values(4))
             end do

             i = i + n_record
          end do
       end do
    end do

    if (i /= recv_buf_i%length+1) then
       error stop "alltoall_dom: receive buffer length mismatch"
    end if

    deallocate (values)
  end subroutine alltoall_dom

  
  subroutine alltoall
    implicit none

    integer :: i

    call MPI_Alltoall ( &
         send_lengths, 1, MPI_IN, &
         recv_lengths, 1, MPI_IN, comm)

    recv_offsets(1) = 0

    do i = 2, n_process
       recv_offsets(i) = &
            recv_offsets(i-1) + recv_lengths(i-1)
    end do

    recv_buf_i%length = sum(recv_lengths)

    if (.not. allocated(recv_buf_i%elts)) then
       allocate(recv_buf_i%elts(recv_buf_i%length))

    else if (size(recv_buf_i%elts) < recv_buf_i%length) then
       deallocate(recv_buf_i%elts)
       allocate(recv_buf_i%elts(recv_buf_i%length))
    end if

    if (recv_buf_i%length > 0) recv_buf_i%elts = 0

    call MPI_Alltoallv ( &
         send_buf_i%elts, send_lengths, send_offsets, MPI_IN, &
         recv_buf_i%elts, recv_lengths, recv_offsets, MPI_IN, comm)
  end subroutine alltoall

  
  subroutine comm_masks_mpi (l)
    ! Communicate node and edge mask information between MPI processes.
    implicit none

    integer, intent(in) :: l

    integer :: d_dest, d_src, dest, src
    integer :: i, id, kk
    integer :: r_dest, r_src
    integer :: ierr

    send_buf_i%length = 0

    ! Pack data for remote processes.
    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf_i%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                   id = grid(d_src)%pack(AT_NODE,dest)%elts(i)

                   if (l == NONE .or. &
                        l == grid(d_src)%level%elts(abs(id)+1)) then

                      call append (send_buf_i, grid(d_src)%mask_n%elts(abs(id)+1))
                   end if
                end do

                do i = 1, grid(d_src)%pack(AT_EDGE,dest)%length
                   id = grid(d_src)%pack(AT_EDGE,dest)%elts(i)

                   if (l == NONE .or. &
                        l == grid(d_src)%level%elts(abs(id)/EDGE+1)) then

                      call append (send_buf_i, grid(d_src)%mask_e%elts(abs(id)+1))
                   end if
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    ! Determine receive lengths and offsets.
    recv_buf_i%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf_i%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             src = glo_id(r_src,d_src) + 1

             do d_dest = 1, n_domain(rank+1)
                do i = 1, grid(d_dest)%unpk(AT_NODE,src)%length
                   id = abs(grid(d_dest)%unpk(AT_NODE,src)%elts(i))

                   if (l == NONE .or. &
                        l == grid(d_dest)%level%elts(id+1)) then
                      recv_buf_i%length = recv_buf_i%length + 1
                   end if
                end do

                do i = 1, grid(d_dest)%unpk(AT_EDGE,src)%length
                   id = abs(grid(d_dest)%unpk(AT_EDGE,src)%elts(i))

                   if (l == NONE .or. &
                        l == grid(d_dest)%level%elts(id/EDGE+1)) then
                      recv_buf_i%length = recv_buf_i%length + 1
                   end if
                end do
             end do
          end do
       end if

       recv_lengths(r_src) = &
            recv_buf_i%length - recv_offsets(r_src)
    end do

    ! Ensure the receive buffer is large enough.
    if (.not. allocated(recv_buf_i%elts)) then
       allocate (recv_buf_i%elts(max(1,recv_buf_i%length)))

    else if (size(recv_buf_i%elts) < recv_buf_i%length) then
       deallocate (recv_buf_i%elts)
       allocate (recv_buf_i%elts(max(1,recv_buf_i%length)))
    end if

    if (recv_buf_i%length > 0) then
       recv_buf_i%elts(1:recv_buf_i%length) = 0
    end if

    ! Optional diagnostic:
    ! call check_alltoall_lengths()

    call MPI_Alltoallv ( &
         send_buf_i%elts, send_lengths, send_offsets, MPI_IN, &
         recv_buf_i%elts, recv_lengths, recv_offsets, MPI_IN, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "comm_masks_mpi: MPI_Alltoallv failed"
    end if

    ! Communicate between domains owned by the same process.
    call comm_masks

    ! Unpack data received from remote processes.
    kk = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,src)%length
                id = grid(d_dest)%unpk(AT_NODE,src)%elts(i)

                if (l == NONE .or. &
                     l == grid(d_dest)%level%elts(abs(id)+1)) then

                   kk = kk + 1

                   if (kk > recv_buf_i%length) then
                      error stop &
                           "comm_masks_mpi: node receive buffer overrun"
                   end if

                   grid(d_dest)%mask_n%elts(abs(id)+1) = &
                        recv_buf_i%elts(kk)
                end if
             end do

             do i = 1, grid(d_dest)%unpk(AT_EDGE,src)%length
                id = grid(d_dest)%unpk(AT_EDGE,src)%elts(i)

                if (l == NONE .or. &
                     l == grid(d_dest)%level%elts(abs(id)/EDGE+1)) then

                   kk = kk + 1

                   if (kk > recv_buf_i%length) then
                      error stop &
                           "comm_masks_mpi: edge receive buffer overrun"
                   end if

                   grid(d_dest)%mask_e%elts(abs(id)+1) = recv_buf_i%elts(kk)
                end if
             end do
          end do
       end do
    end do

    if (kk /= recv_buf_i%length) then
       error stop "comm_masks_mpi: receive buffer length mismatch"
    end if
  end subroutine comm_masks_mpi

  
  subroutine deadlock_test (flag)
    implicit none
    integer, intent(in) :: flag

    logical  :: all_done
    real(dp) :: now, t_start
    real(dp), parameter :: timeout_time = 100.0_dp  ! seconds

    ! Nothing to wait for
    if (nreq <= 0) return

    t_start = MPI_Wtime ()

    do
       call MPI_Testall (nreq, req(1:nreq), all_done, MPI_STATUSES_IGNORE)
       if (all_done) exit

       now = MPI_Wtime ()
       if ((now - t_start) >= timeout_time) then
          write (6,'(a,i0,a,i0)') &
               "ERROR: boundary update deadlocked at call ", flag, &
               " on rank ", rank
          call abort_run
       end if
    end do
  end subroutine deadlock_test

  
  subroutine update_bdry1_0 (field, l_start, l_end, flag_in)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l_start, l_end
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    flag = 9999
    if (present(flag_in)) flag = flag_in

    call update_bdry__start1 (field, l_start, l_end)

    if (deadlock) call deadlock_test (flag)

    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_0

  
  subroutine update_bdry1_1 (field, l_start, l_end, flag_in)
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l_start, l_end
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    flag = 9999
    if (present(flag_in)) flag = flag_in

    call update_bdry__start1 (field, l_start, l_end)

    if (deadlock) call deadlock_test (flag)

    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_1

  
  subroutine update_bdry1_2 (field, l_start, l_end, flag_in)
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l_start, l_end
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    flag = 9999
    if (present(flag_in)) flag = flag_in

    call update_bdry__start1 (field, l_start, l_end)

    if (deadlock) call deadlock_test(flag)

    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_2

  
  subroutine update_bdry_0 (field, l, flag_in)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    if (present(flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start (field, l)
    if (deadlock) call deadlock_test(flag)
    call update_bdry__finish( field, l)
  end subroutine update_bdry_0

  
  subroutine update_bdry_1 (field, l, flag_in)
    ! Update a rank-1 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    if (present(flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start (field, l)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry_1

  
  subroutine update_bdry_2 (field, l, flag_in)
    ! Update a rank-2 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l
    integer, optional, intent(in)    :: flag_in

    integer :: flag

    if (present(flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start (field, l)
    if (deadlock) call deadlock_test(flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry_2

  
  subroutine update_bdry__start_0(field, l)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__start1(field, level_start-1, level_end)
    else
       call update_bdry__start1(field, l, l)
    end if
  end subroutine update_bdry__start_0

  
  subroutine update_bdry__start_1 (field, l)
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    end if
  end subroutine update_bdry__start_1

  
  subroutine update_bdry__start_2 (field, l)
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    end if
  end subroutine update_bdry__start_2

  
  subroutine update_bdry__start1_0( field, l_start, l_end)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l_start, l_end

    integer :: d_src, d_dest, dest
    integer :: i, id, lev, multipl
    integer :: r, r_dest, r_src, tag

    if (field%bdry_uptodate) return

    tag = TAG_BDRY_S
    field%bdry_tag = tag

    if (field%pos == AT_EDGE) then
       multipl = EDGE
    else
       multipl = 1
    end if

    send_buf%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i = 1, grid(d_src)%pack(field%pos,dest)%length
                   id  = grid(d_src)%pack(field%pos,dest)%elts(i)
                   lev = grid(d_src)%level%elts(id/multipl+1)

                   if (lev >= l_start .and. lev <= l_end) then
                      call append (send_buf, field%data(d_src)%elts(id+1))
                   end if
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    recv_buf%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             do d_dest = 1, n_domain(rank+1)
                dest = glo_id(r_src,d_src) + 1

                do i = 1, grid(d_dest)%unpk(field%pos,dest)%length
                   id = abs(grid(d_dest)%unpk(field%pos,dest)%elts(i))
                   lev = grid(d_dest)%level%elts(id/multipl+1)

                   if (lev >= l_start .and. lev <= l_end) then
                      recv_buf%length = recv_buf%length + 1
                   end if
                end do
             end do
          end do
       end if

       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (.not. allocated(recv_buf%elts)) then
       allocate(recv_buf%elts(recv_buf%length))
    else if (size(recv_buf%elts) < recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    if (recv_buf%length > 0) recv_buf%elts = 0.0_dp

    nreq = 0

    ! Post all receives before sends.
    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Irecv( &
            recv_buf%elts( &
            recv_offsets(r)+1:recv_offsets(r)+recv_lengths(r)), &
            recv_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Isend( &
            send_buf%elts( &
            send_offsets(r)+1:send_offsets(r)+send_lengths(r)), &
            send_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    call cp_bdry_inside (field)
  end subroutine update_bdry__start1_0

  
  subroutine update_bdry__start1_1 (field, l_start, l_end)
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l_start, l_end

    integer :: d_dest, d_src, dest
    integer :: i, i1, id, lev, multipl, pos
    integer :: r, r_dest, r_src, tag

    if (all(field%bdry_uptodate)) return

    tag = TAG_BDRY_V
    field%bdry_tag = tag

    send_buf%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i1 = 1, size(field)
                   pos = field(i1)%pos

                   if (pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if

                   do i = 1, grid(d_src)%pack(pos,dest)%length
                      id  = grid(d_src)%pack(pos,dest)%elts(i)
                      lev = grid(d_src)%level%elts(id/multipl+1)

                      if (lev >= l_start .and. lev <= l_end) then
                         call append (send_buf, field(i1)%data(d_src)%elts(id+1))
                      end if
                   end do
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    recv_buf%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             dest = glo_id(r_src,d_src) + 1

             do d_dest = 1, n_domain(rank+1)
                do i1 = 1, size(field)
                   pos = field(i1)%pos

                   if (pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if

                   do i = 1, grid(d_dest)%unpk(pos,dest)%length
                      id = abs(grid(d_dest)%unpk(pos,dest)%elts(i))
                      lev = grid(d_dest)%level%elts(id/multipl+1)

                      if (lev >= l_start .and. lev <= l_end) then
                         recv_buf%length = recv_buf%length + 1
                      end if
                   end do
                end do
             end do
          end do
       end if

       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (.not. allocated(recv_buf%elts)) then
       allocate(recv_buf%elts(recv_buf%length))
    else if (size(recv_buf%elts) < recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    if (recv_buf%length > 0) recv_buf%elts = 0.0_dp

    nreq = 0

    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Irecv ( &
            recv_buf%elts( &
            recv_offsets(r)+1:recv_offsets(r)+recv_lengths(r)), &
            recv_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Isend ( &
            send_buf%elts( &
            send_offsets(r)+1:send_offsets(r)+send_lengths(r)), &
            send_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    call cp_bdry_inside (field)
  end subroutine update_bdry__start1_1

  
  subroutine update_bdry__start1_2 (field, l_start, l_end)
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l_start, l_end

    integer :: d_dest, d_src, dest
    integer :: i, i1, i2, id, lev, multipl, pos
    integer :: r, r_dest, r_src, tag

    if (all(field%bdry_uptodate)) return

    tag = TAG_BDRY_A
    field%bdry_tag = tag

    send_buf%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i2 = 1, size(field,2)
                   do i1 = 1, size(field,1)
                      pos = field(i1,i2)%pos

                      if (pos == AT_EDGE) then
                         multipl = EDGE
                      else
                         multipl = 1
                      end if

                      do i = 1, grid(d_src)%pack(pos,dest)%length
                         id  = grid(d_src)%pack(pos,dest)%elts(i)
                         lev = grid(d_src)%level%elts(id/multipl+1)

                         if (lev >= l_start .and. lev <= l_end) then
                            call append ( &
                                 send_buf, &
                                 field(i1,i2)%data(d_src)%elts(id+1))
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    recv_buf%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             dest = glo_id(r_src,d_src) + 1

             do d_dest = 1, n_domain(rank+1)
                do i2 = 1, size(field,2)
                   do i1 = 1, size(field,1)
                      pos = field(i1,i2)%pos

                      if (pos == AT_EDGE) then
                         multipl = EDGE
                      else
                         multipl = 1
                      end if

                      do i = 1, grid(d_dest)%unpk(pos,dest)%length
                         id = abs(grid(d_dest)%unpk(pos,dest)%elts(i))
                         lev = grid(d_dest)%level%elts(id/multipl+1)

                         if (lev >= l_start .and. lev <= l_end) then
                            recv_buf%length = recv_buf%length + 1
                         end if
                      end do
                   end do
                end do
             end do
          end do
       end if

       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (.not. allocated(recv_buf%elts)) then
       allocate(recv_buf%elts(recv_buf%length))
    else if (size(recv_buf%elts) < recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    if (recv_buf%length > 0) recv_buf%elts = 0.0_dp

    nreq = 0

    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Irecv ( &
            recv_buf%elts( &
            recv_offsets(r)+1:recv_offsets(r)+recv_lengths(r)), &
            recv_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle

       nreq = nreq + 1

       call MPI_Isend ( &
            send_buf%elts( &
            send_offsets(r)+1:send_offsets(r)+send_lengths(r)), &
            send_lengths(r), MPI_DP, r-1, tag, comm, req(nreq))
    end do

    call cp_bdry_inside (field)
  end subroutine update_bdry__start1_2

  
  subroutine update_bdry__finish_0(field, l)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    end if
  end subroutine update_bdry__finish_0

  
  subroutine update_bdry__finish_1(field, l)
    ! Finish boundary updates for a rank-1 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    end if
  end subroutine update_bdry__finish_1

  
  subroutine update_bdry__finish_2 (field, l)
    ! Finish boundary updates for a rank-2 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l

    if (l == NONE) then
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    end if
  end subroutine update_bdry__finish_2

  
  subroutine update_bdry__finish1_0 (field, l_start, l_end)
    implicit none

    type(Float_Field), intent(inout) :: field
    integer,           intent(in)    :: l_start, l_end

    integer :: r_src, d_src, d_dest
    integer :: id, i, k, lev, multipl, src

    if (field%bdry_uptodate) return

    if (field%bdry_tag == -1) then
       error stop "update_bdry__finish1_0: finish without start"
    end if

    if (field%pos == AT_EDGE) then
       multipl = EDGE
    else
       multipl = 1
    end if

    if (nreq > 0) then
       call MPI_Waitall (nreq, req(1:nreq), MPI_STATUSES_IGNORE)
    end if

    k = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(field%pos,src)%length
                id = grid(d_dest)%unpk(field%pos,src)%elts(i)
                lev = grid(d_dest)%level%elts(abs(id)/multipl+1)

                if (lev >= l_start .and. lev <= l_end) then
                   k = k + 1

                   field%data(d_dest)%elts(abs(id)+1) = &
                        recv_buf%elts(k)

                   if (id < 0 .and. field%pos == AT_EDGE) then
                      field%data(d_dest)%elts(abs(id)+1) = &
                           -field%data(d_dest)%elts(abs(id)+1)
                   end if
                end if
             end do
          end do
       end do
    end do

    if (k /= recv_buf%length) then
       error stop "update_bdry__finish1_0: receive count mismatch"
    end if

    ! Assumes this is called either for one level or for all levels
    ! whose boundaries are to be considered current.
    if (l_start < l_end) field%bdry_uptodate = .true.

    field%bdry_tag = -1
    nreq = 0
  end subroutine update_bdry__finish1_0

  
  subroutine update_bdry__finish1_1 (field, l_start, l_end)
    ! Complete boundary communication for a rank-1 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:)
    integer,           intent(in)    :: l_start, l_end

    integer :: d_dest, d_src, src
    integer :: id, i, i1, k, lev, multipl, pos, r_src

    if (all(field%bdry_uptodate)) return

    if (.not. all(field%bdry_uptodate .or. field%bdry_tag /= -1)) then
       error stop &
            "update_bdry__finish1_1: finish without matching start"
    end if

    if (nreq > 0) then
       call MPI_Waitall (nreq, req(1:nreq), MPI_STATUSES_IGNORE)
    end if

    k = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i1 = 1, size(field)
                pos = field(i1)%pos

                if (pos == AT_EDGE) then
                   multipl = EDGE
                else
                   multipl = 1
                end if

                do i = 1, grid(d_dest)%unpk(pos,src)%length
                   id = grid(d_dest)%unpk(pos,src)%elts(i)
                   lev = grid(d_dest)%level%elts(abs(id)/multipl+1)

                   if (lev >= l_start .and. lev <= l_end) then
                      k = k + 1

                      field(i1)%data(d_dest)%elts(abs(id)+1) = &
                           recv_buf%elts(k)

                      if (id < 0 .and. pos == AT_EDGE) then
                         field(i1)%data(d_dest)%elts(abs(id)+1) = &
                              -field(i1)%data(d_dest)%elts(abs(id)+1)
                      end if
                   end if
                end do
             end do
          end do
       end do
    end do

    if (k /= recv_buf%length) then
       error stop "update_bdry__finish1_1: receive count mismatch"
    end if

    if (l_start < l_end) field%bdry_uptodate = .true.

    field%bdry_tag = -1
    nreq = 0
  end subroutine update_bdry__finish1_1

  
  subroutine update_bdry__finish1_2 (field, l_start, l_end)
    ! Complete boundary communication for a rank-2 field array.
    implicit none

    type(Float_Field), intent(inout) :: field(:,:)
    integer,           intent(in)    :: l_start, l_end

    integer :: d_dest, d_src, src
    integer :: i, i1, i2, id, k, lev, multipl, pos, r_src

    if (all(field%bdry_uptodate)) return

    if (.not. all(field%bdry_uptodate .or. field%bdry_tag /= -1)) then
       error stop &
            "update_bdry__finish1_2: finish without matching start"
    end if

    if (nreq > 0) then
       call MPI_Waitall (nreq, req(1:nreq), MPI_STATUSES_IGNORE)
    end if

    k = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i2 = 1, size(field,2)
                do i1 = 1, size(field,1)
                   pos = field(i1,i2)%pos

                   if (pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if

                   do i = 1, grid(d_dest)%unpk(pos,src)%length
                      id = grid(d_dest)%unpk(pos,src)%elts(i)
                      lev = grid(d_dest)%level%elts(abs(id)/multipl+1)

                      if (lev >= l_start .and. lev <= l_end) then
                         k = k + 1

                         field(i1,i2)%data(d_dest)%elts(abs(id)+1) = &
                              recv_buf%elts(k)

                         if (id < 0 .and. pos == AT_EDGE) then
                            field(i1,i2)%data(d_dest)%elts(abs(id)+1) = &
                                 -field(i1,i2)%data(d_dest)%elts(abs(id)+1)
                         end if
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

    if (k /= recv_buf%length) then
       error stop "update_bdry__finish1_2: receive count mismatch"
    end if

    if (l_start < l_end) field%bdry_uptodate = .true.

    field%bdry_tag = -1
    nreq = 0
  end subroutine update_bdry__finish1_2

  
  subroutine comm_nodes3_mpi(get, set)
    implicit none

    procedure(coord_get) :: get
    procedure(coord_set) :: set

    integer :: r_dest, r_src
    integer :: d_src, d_dest, dest, src
    integer :: id, i, k, ierr
    type(Coord) :: c

    send_buf%length = 0

    ! Pack data for remote processes.
    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                   id = grid(d_src)%pack(AT_NODE,dest)%elts(i)
                   c  = get(grid(d_src), id)

                   k = send_buf%length
                   call extend(send_buf, 3, 0.0_dp)

                   send_buf%elts(k+1:k+3) = [c%x, c%y, c%z]
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = &
            send_buf%length - send_offsets(r_dest)
    end do

    ! Determine receive lengths and offsets.
    recv_buf%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             src = glo_id(r_src,d_src) + 1

             do d_dest = 1, n_domain(rank+1)
                recv_buf%length = recv_buf%length + &
                     3*grid(d_dest)%unpk(AT_NODE,src)%length
             end do
          end do
       end if

       recv_lengths(r_src) = &
            recv_buf%length - recv_offsets(r_src)
    end do

    if (.not. allocated(recv_buf%elts)) then
       allocate(recv_buf%elts(max(1,recv_buf%length)))

    else if (size(recv_buf%elts) < recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(max(1,recv_buf%length)))
    end if

    if (recv_buf%length > 0) then
       recv_buf%elts(1:recv_buf%length) = 0.0_dp
    end if

    call MPI_Alltoallv( &
         send_buf%elts, send_lengths, send_offsets, MPI_DP, &
         recv_buf%elts, recv_lengths, recv_offsets, MPI_DP, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "comm_nodes3_mpi: MPI_Alltoallv failed"
    end if

    ! Communicate between domains on this MPI process.
    call comm_nodes3(get, set)

    ! Unpack remote data.
    k = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,src)%length
                if (k+3 > recv_buf%length) then
                   error stop "comm_nodes3_mpi: receive buffer overrun"
                end if

                id = grid(d_dest)%unpk(AT_NODE,src)%elts(i)

                c%x = recv_buf%elts(k+1)
                c%y = recv_buf%elts(k+2)
                c%z = recv_buf%elts(k+3)

                call set(grid(d_dest), id, c)

                k = k + 3
             end do
          end do
       end do
    end do

    if (k /= recv_buf%length) then
       error stop "comm_nodes3_mpi: receive buffer length mismatch"
    end if
  end subroutine comm_nodes3_mpi

  
  subroutine comm_nodes9_mpi (get, set)
    implicit none

    procedure(get9)    :: get
    procedure(set9)    :: set

    integer, parameter :: NVAL = 7

    real(dp) :: val(NVAL)

    integer :: r_dest, r_src
    integer :: d_src, d_dest, dest, src
    integer :: id, i, k, ierr

    send_buf%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)
             do d_dest = 1, n_domain(r_dest)
                dest = glo_id(r_dest,d_dest) + 1

                do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                   id = grid(d_src)%pack(AT_NODE,dest)%elts(i)

                   call get(grid(d_src), id, val)

                   k = send_buf%length
                   call extend(send_buf, NVAL, 0.0_dp)

                   send_buf%elts(k+1:k+NVAL) = val
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = &
            send_buf%length - send_offsets(r_dest)
    end do

    recv_buf%length = 0

    do r_src = 1, n_process
       recv_offsets(r_src) = recv_buf%length

       if (r_src /= rank+1) then
          do d_src = 1, n_domain(r_src)
             src = glo_id(r_src,d_src) + 1

             do d_dest = 1, n_domain(rank+1)
                recv_buf%length = recv_buf%length + &
                     NVAL*grid(d_dest)%unpk(AT_NODE,src)%length
             end do
          end do
       end if

       recv_lengths(r_src) = &
            recv_buf%length - recv_offsets(r_src)
    end do

    if (.not. allocated(recv_buf%elts)) then
       allocate(recv_buf%elts(max(1,recv_buf%length)))

    else if (size(recv_buf%elts) < recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(max(1,recv_buf%length)))
    end if

    if (recv_buf%length > 0) then
       recv_buf%elts(1:recv_buf%length) = 0.0_dp
    end if

    call MPI_Alltoallv( &
         send_buf%elts, send_lengths, send_offsets, MPI_DP, &
         recv_buf%elts, recv_lengths, recv_offsets, MPI_DP, &
         comm, ierr)

    if (ierr /= MPI_SUCCESS) then
       error stop "comm_nodes9_mpi: MPI_Alltoallv failed"
    end if

    call comm_nodes9(get, set)

    k = 0

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       do d_src = 1, n_domain(r_src)
          src = glo_id(r_src,d_src) + 1

          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,src)%length
                if (k+NVAL > recv_buf%length) then
                   error stop "comm_nodes9_mpi: receive buffer overrun"
                end if

                id = grid(d_dest)%unpk(AT_NODE,src)%elts(i)

                call set( &
                     grid(d_dest), id, &
                     recv_buf%elts(k+1:k+NVAL))

                k = k + NVAL
             end do
          end do
       end do
    end do

    if (k /= recv_buf%length) then
       error stop "comm_nodes9_mpi: receive buffer length mismatch"
    end if
  end subroutine comm_nodes9_mpi

  
  subroutine comm_patch_conn_mpi
    implicit none

    integer, parameter :: SEND_RECORD_SIZE = 6
    integer, parameter :: RECV_DATA_SIZE   = 4

    integer :: r_src, r_dest
    integer :: d_src, d_dest, d_glo
    integer :: i, b, c, p, s, k
    integer :: rot, rot_shift
    integer :: d, ngh_pa, typ, l_par

    integer :: st(RECV_DATA_SIZE)
    logical :: is_pole

    send_buf_i%length = 0

    do r_dest = 1, n_process
       send_offsets(r_dest) = send_buf_i%length

       if (r_dest /= rank+1) then
          do d_src = 1, n_domain(rank+1)

             if (mod(grid(d_src)%send_pa_all%length, 4) /= 0) then
                error stop &
                     "comm_patch_conn_mpi: send_pa_all has invalid length"
             end if

             do d_dest = 1, n_domain(r_dest)
                do i = 1, grid(d_src)%send_pa_all%length, 4
                   b = grid(d_src)%send_pa_all%elts(i)
                   c = grid(d_src)%send_pa_all%elts(i+1)
                   p = grid(d_src)%send_pa_all%elts(i+2)
                   s = grid(d_src)%send_pa_all%elts(i+3)

                   typ = grid(d_src)%bdry_patch%elts(b+1)%side

                   if (typ < 1) cycle

                   d_glo   = grid(d_src)%neigh(typ)
                   is_pole = d_glo == POLE

                   if (is_pole) then
                      d_glo = grid(d_src)%neigh_over_pole(c+1)
                      l_par = grid(d_src)%patch%elts(p+1)%level - 1

                      if (grid(d_src)%neigh_pa_over_pole%length < &
                           2*l_par+c+1) then
                         ngh_pa = 0
                      else
                         ngh_pa = grid(d_src)%neigh_pa_over_pole%elts( &
                              2*l_par+c+1)
                      end if
                   else
                      ngh_pa = grid(d_src)%bdry_patch%elts(b+1)%neigh
                   end if

                   if (ngh_pa == 0) cycle
                   if (d_glo /= glo_id(r_dest,d_dest)) cycle

                   rot = grid(d_src)%neigh_rot(typ)
                   rot_shift = &
                        (2*rot_direction(grid(d_src),typ)-1)*rot

                   call append (send_buf_i, d_dest)
                   call append (send_buf_i, glo_id(rank+1,d_src))
                   call append (send_buf_i, p)

                   if (is_pole) then
                      call append (send_buf_i, c)
                   else
                      call append (send_buf_i, modulo(c+rot_shift,4))
                   end if

                   call append (send_buf_i, ngh_pa)

                   if (is_pole) then
                      call append(send_buf_i, s)
                   else
                      call append ( &
                           send_buf_i, &
                           modulo(s+rot_shift+2,4) + 4*(s/4))
                   end if
                end do
             end do
          end do
       end if

       send_lengths(r_dest) = &
            send_buf_i%length - send_offsets(r_dest)

       if (mod(send_lengths(r_dest),SEND_RECORD_SIZE) /= 0) then
          error stop &
               "comm_patch_conn_mpi: invalid packed record length"
       end if
    end do

    call alltoall

    call comm_patch_conn

    do r_src = 1, n_process
       if (r_src == rank+1) cycle

       if (mod(recv_lengths(r_src),SEND_RECORD_SIZE) /= 0) then
          error stop &
               "comm_patch_conn_mpi: invalid received record length"
       end if

       do k = recv_offsets(r_src)+1, &
            recv_offsets(r_src)+recv_lengths(r_src), &
            SEND_RECORD_SIZE

          if (k+SEND_RECORD_SIZE-1 > recv_buf_i%length) then
             error stop &
                  "comm_patch_conn_mpi: receive buffer overrun"
          end if

          d     = recv_buf_i%elts(k)
          d_src = recv_buf_i%elts(k+1) + 1
          st    = recv_buf_i%elts(k+2:k+5)

          if (d < 1 .or. d > size(grid)) then
             error stop &
                  "comm_patch_conn_mpi: invalid destination domain"
          end if

          if (d_src < 1 .or. d_src > N_GLO_DOMAIN) then
             error stop &
                  "comm_patch_conn_mpi: invalid source domain"
          end if

          call append (grid(d)%recv_pa(d_src), st(1))
          call append (grid(d)%recv_pa(d_src), st(2))
          call append (grid(d)%recv_pa(d_src), st(3))
          call append (grid(d)%recv_pa(d_src), st(4))
       end do
    end do
  end subroutine comm_patch_conn_mpi

  
  function sync_max_int (val) result(value)
    implicit none
    integer, intent(in) :: val
    integer             :: value

    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_IN, MPI_MAX, comm)
    value = val_glo
  end function sync_max_int


  function sync_min_int (val) result(value)
    implicit none
    integer, intent(in) :: val
    integer             :: value

    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_IN, MPI_MIN, comm)
    value = val_glo
  end function sync_min_int

  
  function sync_max_real_0 (val) result(value)
    implicit none
    real(dp), intent(in) :: val
    real(dp)             :: value

    real(dp) :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_DP, MPI_MAX, comm)
    value = val_glo
  end function sync_max_real_0


  function sync_max_real_1 (val) result(value)
    implicit none
    real(dp), intent(in)  :: val(:)
    real(dp), allocatable :: value(:)

    integer                             :: n
    real(dp), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (value(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_DP, MPI_MAX, comm)
    value = val_glo
  end function sync_max_real_1


  function sync_min_real_0 (val) result(value)
    implicit none
    real(dp), intent(in) :: val
    real(dp)             ::  value

    real(dp) :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_DP, MPI_MIN, comm)
    value = val_glo
  end function sync_min_real_0

  
  function sync_min_real_1 (val) result(value)
    implicit none

    real(dp), intent(in)  :: val(:)
    real(dp), allocatable :: value(:)
    
    integer                             :: n
    real(dp), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (value(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_DP, MPI_MIN, comm)
    value = val_glo
  end function sync_min_real_1
  
  
  function sum_real_0 (val) result(value)
    implicit none
    real(dp), intent(in) :: val
    real(dp) :: value

    real(dp) :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_DP, MPI_SUM, comm)
    value = val_glo
  end function sum_real_0
  
  
  function sum_real_1 (val) result(value)
    implicit none
    real(dp), intent(in)  :: val(:)
    real(dp), allocatable :: value(:)

    integer                             :: n
    real(dp), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (value(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_DP, MPI_SUM, comm)
    value = val_glo
  end function sum_real_1

  function sum_int_0 (val) result(value)
    implicit none
    integer, intent(in) :: val
    integer             :: value

    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_IN, MPI_SUM, comm)
    value = val_glo
  end function sum_int_0

  function sum_int_1 (val) result(value)
    implicit none
    integer, intent(in)  :: val(:)
    integer, allocatable :: value(:)

    integer                            :: n
    integer, dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (value(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_IN, MPI_SUM, comm)
    value = val_glo
  end function sum_int_1

  subroutine start_timing
    implicit none
    times(1) = MPI_Wtime () 
  end subroutine start_timing

  subroutine stop_timing
    implicit none
    times(2) = MPI_Wtime ()  
  end subroutine stop_timing

  real(dp) function get_timing ()
    implicit none
    get_timing = times(2) - times(1)
  end function get_timing

  subroutine sync_array (arr, n)
    implicit none
    integer,  intent(in)    :: n
    real(dp), intent(inout) :: arr(n)

    real(dp), allocatable :: vloc(:), vmax(:)
    integer,  allocatable :: mloc(:), mmax(:)
    integer               :: i

    allocate (vloc(n), vmax(n), mloc(n), mmax(n))

    do i = 1, n
       if (abs(arr(i) - sync_val) > eps(1.0_dp)) then
          mloc(i) = 1
          vloc(i) = arr(i)
       else
          mloc(i) = 0
          vloc(i) = -huge(1.0_dp)   ! neutral element for MAX
       end if
    end do

    call MPI_Allreduce (mloc, mmax, n, MPI_IN, MPI_MAX, comm)
    call MPI_Allreduce (vloc, vmax, n, MPI_DP, MPI_MAX, comm)

    do i = 1, n
       if (mmax(i) == 1) then
          arr(i) = vmax(i)
       else
          arr(i) = sync_val
       end if
    end do

    deallocate (vloc, vmax, mloc, mmax)
  end subroutine sync_array

  subroutine gatherv_int(n_loc, n_glo, vec_loc, vec_glo)
    implicit none

    integer,              intent(in)  :: n_loc
    integer,              intent(out) :: n_glo(n_process)
    integer,              intent(in)  :: vec_loc(n_loc)
    integer, allocatable, intent(out) :: vec_glo(:)

    integer :: displs(n_process)
    integer :: n_tot

    call gather_int(n_loc, n_glo, displs)

    if (rank == 0) then
       n_tot = sum(n_glo)
    else
       n_tot = 0
    end if

    allocate(vec_glo(n_tot))

    call MPI_Gatherv( &
         vec_loc, n_loc,         MPI_IN, &
         vec_glo, n_glo, displs, MPI_IN, &
         0, comm)
  end subroutine gatherv_int

  
  subroutine gatherv_real4 (n_loc, n_glo, vec_loc, vec_glo)
    implicit none

    integer,              intent(in)  :: n_loc
    integer,              intent(out) :: n_glo(n_process)
    real(4),              intent(in)  :: vec_loc(n_loc)
    real(4), allocatable, intent(out) :: vec_glo(:)

    integer :: displs(n_process)
    integer :: n_tot

    call gather_int(n_loc, n_glo, displs)

    if (rank == 0) then
       n_tot = sum(n_glo)
    else
       n_tot = 0
    end if

    allocate(vec_glo(n_tot))

    call MPI_Gatherv( &
         vec_loc, n_loc,         MPI_SP, &
         vec_glo, n_glo, displs, MPI_SP, &
         0, comm)
  end subroutine gatherv_real4


  subroutine gatherv_real8 (n_loc, n_glo, vec_loc, vec_glo)
    implicit none

    integer,              intent(in)  :: n_loc
    integer,              intent(out) :: n_glo(n_process)
    real(dp),             intent(in)  :: vec_loc(n_loc)
    real(dp), allocatable, intent(out) :: vec_glo(:)

    integer :: displs(n_process)
    integer :: n_tot

    call gather_int(n_loc, n_glo, displs)

    if (rank == 0) then
       n_tot = sum(n_glo)
    else
       n_tot = 0
    end if

    allocate(vec_glo(n_tot))

    call MPI_Gatherv( &
         vec_loc, n_loc,         MPI_DP, &
         vec_glo, n_glo, displs, MPI_DP, &
         0, comm)
  end subroutine gatherv_real8

  
  subroutine gather_int(n_loc, n_glo, displs)
    implicit none

    integer, intent(in)  :: n_loc
    integer, intent(out) :: n_glo(n_process)
    integer, intent(out) :: displs(n_process)

    integer :: r

    ! Only rank 0 initially receives valid n_glo from MPI_Gather.
    call MPI_Gather(n_loc, 1, MPI_IN, n_glo, 1, MPI_IN, 0, comm)

    if (rank == 0) then
       displs(1) = 0

       do r = 2, n_process
          displs(r) = displs(r-1) + n_glo(r-1)
       end do
    end if

    ! Define both arrays on every rank.
    call MPI_Bcast(n_glo,  n_process, MPI_IN, 0, comm)
    call MPI_Bcast(displs, n_process, MPI_IN, 0, comm)
  end subroutine gather_int

  
end module comm_mpi_mod
