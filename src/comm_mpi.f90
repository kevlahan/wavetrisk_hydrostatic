module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
  integer                            :: nreq
  type(Int_Array)                    :: recv_buf_i, send_buf_i
  type(Float_Array)                  :: recv_buf, send_buf
  integer, dimension(:), allocatable :: recv_lengths, recv_offsets, req, send_lengths, send_offsets
  real(8), dimension(2)              :: times

  logical, parameter                 :: deadlock = .true. ! test for communication deadlocks
  
  interface sum_real
     procedure :: sum_real_0, sum_real_1
  end interface sum_real

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
    req          = 0
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  subroutine write_load_conn (id, run_id)
    ! Write out load distribution and connectivity for load balancing
    use mpi
    implicit none
    integer      :: id
    character(*) :: run_id
    
    integer                                           :: d, ii, r, sz
    integer, parameter                                :: fid = 599
    integer, dimension(n_process)                     :: displs, rcounts
    integer, dimension(:), allocatable                :: n_active_loc
    integer, dimension(N_GLO_DOMAIN*(N_GLO_DOMAIN+1)) :: n_active_glo
    character(255)                                    :: filename

    ! Size of load data for current rank
    sz = size(grid) * (N_GLO_DOMAIN+1)

    ! Gather rcounts on root
    call MPI_Gather (sz, 1, MPI_INTEGER, rcounts, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    ! Set displacements for contiguous positions
    displs(1) = 0
    do r = 2, n_process
       displs(r) = displs(r-1) + rcounts(r-1)
    end do

    ! Set load data on this rank
    allocate (n_active_loc(sz))
    n_active_loc = 0
    ii = 1 
    do d = 1, size(grid)
       n_active_loc(ii) = domain_load(grid(d))
       
       n_active_loc(ii+1:ii+N_GLO_DOMAIN) = (grid(d)%pack(AT_NODE,:)%length + grid(d)%pack(AT_EDGE,:)%length + &
            grid(d)%unpk(AT_NODE,:)%length + grid(d)%unpk(AT_EDGE,:)%length)/2
       
       ii = ii + 1 + N_GLO_DOMAIN
    end do

    ! Gather all load data onto rank 0
    call MPI_Gatherv (n_active_loc, sz, MPI_INTEGER, n_active_glo, rcounts, displs, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    deallocate (n_active_loc)

    ! Write out load data
    if (rank==0) then
       write (filename, '(A,A,I4.4)')  trim (run_id), "_conn.", id
       open (unit = fid, file = trim(filename), recl = 333333, status = 'REPLACE')

       ii = 1
       do d = 1, N_GLO_DOMAIN
          write (fid,'(I10, 99999(1X,I8))') n_active_glo(ii:ii+N_GLO_DOMAIN)
          ii = ii + 1 + N_GLO_DOMAIN
       end do
       close (fid)
    end if
  end subroutine write_load_conn

  subroutine cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    ! Finds load balance between processors
    implicit none
    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance

    call get_load_balance (min_load, avg_load, max_load)
    rel_imbalance = dble(max_load)/avg_load
  end subroutine cal_load_balance

  subroutine get_load_balance (mini, avg, maxi)
    use mpi
    implicit none
    integer :: mini, maxi
    real(8) :: avg

    integer :: d, load, load_sum

    load = 0
    do d = 1, size(grid)
       load = load + domain_load(grid(d))
    end do

    call MPI_Reduce (load, maxi,     1, MPI_INTEGER, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (load, mini,     1, MPI_INTEGER, MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (load, load_sum, 1, MPI_INTEGER, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    avg = dble(load_sum)/dble(n_process)
  end subroutine get_load_balance

  subroutine write_level_mpi (out_rout, l, zlev, eval_pole, filename)
    use mpi
    implicit none
    external       :: out_rout
    integer        :: l, zlev
    logical        :: eval_pole
    character(*)   :: filename

    integer            :: r
    integer, parameter :: funit = 300
            
    do r = 1, n_process
       if (r /= rank+1) then ! write only if our turn, otherwise only wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if

       if (r == 1) then ! first process opens without append to delete old file if existing
          open (unit=funit, file=trim(filename), form='unformatted', status='replace')
       else
          open (unit=funit, file=trim(filename), form='unformatted', access='append', status='old')
       end if

       if (eval_pole) call apply_to_pole (out_rout, l, zlev, funit, .true.)

       call apply_onescale__int (out_rout, l, zlev, 0, 0, funit)

       close (funit)

       call MPI_Barrier (MPI_Comm_World, ierror)
    end do
  end subroutine write_level_mpi

  subroutine init_comm_mpi_mod
    implicit none
    call init (send_buf_i, 0)
    call init (recv_buf_i, 0) 
    call init (send_buf,   0)
    call init (recv_buf,   0) 
  end subroutine init_comm_mpi_mod

  subroutine comm_communication_mpi
    implicit none
    call alltoall_dom (unpack_comm_struct, 4)
    call comm_communication
    call recreate_send_patch_lists
  end subroutine comm_communication_mpi

  subroutine recreate_send_patch_lists
    implicit none
    integer :: d, d_neigh, i, k, l, n, p, r, s, typ 

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
                n = -n
                typ = grid(d)%bdry_patch%elts(n+1)%side
                if (typ < 1) cycle
                d_neigh = grid(d)%neigh(typ)
                if (d_neigh == POLE) then
                   do i = 1, 2
                      call handle_neigh (grid(d), grid(d)%neigh_over_pole(i))
                   end do
                else if (d_neigh == NONE) then
                   cycle
                else
                   call handle_neigh (grid(d), d_neigh)
                end if
             end do
          end do
       end do
    end do
  contains
    subroutine handle_neigh(dom, d0)
      implicit none
      type(Domain) :: dom
      integer      :: d0
      
      integer :: r0

      r0 = owner(d0+1) + 1
      if (r0 == rank+1) return

      if (dom%src_patch(r0,l)%length > 0) then
         if (dom%src_patch(r0,l)%elts(dom%src_patch(r0,l)%length) == p) return ! skip if just added
      end if

      call append (dom%src_patch(r0,l), p)
    end subroutine handle_neigh
  end subroutine recreate_send_patch_lists

  subroutine alltoall_dom (unpack_rout, N)
    implicit none
    external :: unpack_rout
    integer  :: N
    
    integer               :: d_dest, d_src, dest, i, k, length, r_dest, r_src, src
    integer, dimension(N) :: st

    send_buf_i%length = 0 ! reset
    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf_i%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             src = d_src
             dest = glo_id(r_dest,d_dest) + 1
             call append (send_buf_i, grid(src)%send_conn(dest)%length)
             do i = 1, grid(src)%send_conn(dest)%length
                call append (send_buf_i, grid(src)%send_conn(dest)%elts(i))
             end do
             grid(src)%send_conn(dest)%length = 0
          end do
       end do
       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    call alltoall

    i = 1
    do r_src = 1, n_process
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             length = recv_buf_i%elts(i); i = i + 1
             do k = 1, length, N
                st = recv_buf_i%elts(i+0:i+(N-1))
                call unpack_rout (grid(d_dest), glo_id(r_src,d_src), st(1), st(2), st(3), st(4))
                i = i + N
             end do
          end do
       end do
    end do
  end subroutine alltoall_dom

  subroutine check_alltoall_lengths
    use mpi
    implicit none
    integer, dimension(n_process) :: test_recv_len

    call MPI_Alltoall (send_lengths, 1, MPI_INTEGER, test_recv_len, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)

    write (3000+rank,*) test_recv_len-recv_lengths
    close (3000+rank)
    call MPI_Barrier (MPI_Comm_World, ierror)
  end subroutine check_alltoall_lengths

  subroutine alltoall
    use mpi
    implicit none
    integer :: i

    call MPI_Alltoall (send_lengths, 1, MPI_INTEGER, recv_lengths, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)

    recv_offsets(1) = 0

    do i = 2, n_process
       recv_offsets(i) = recv_offsets(i-1) + recv_lengths(i-1)
    end do

    recv_buf_i%length = sum (recv_lengths)

    if (size(recv_buf_i%elts) < recv_buf_i%length) then
       deallocate (recv_buf_i%elts)
       allocate (recv_buf_i%elts(recv_buf_i%length))
       recv_buf_i%elts = 0
    end if

    call MPI_Alltoallv (send_buf_i%elts, send_lengths, send_offsets, MPI_INTEGER, recv_buf_i%elts, &
         recv_lengths, recv_offsets, MPI_INTEGER, MPI_COMM_WORLD, ierror)
  end subroutine alltoall

  subroutine comm_masks_mpi (l)
    ! Communication of mask information in a subdomain between different processes
    use mpi
    implicit none
    integer :: l
    
    integer :: d_dest, d_src, dest, i, id, kk, r_dest, r_src

    send_buf_i%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf_i%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                id = grid(d_src)%pack(AT_NODE,dest)%elts(i)
                if (l == NONE .or. l == grid(d_src)%level%elts(id+1)) then
                   call append (send_buf_i, grid(d_src)%mask_n%elts(abs(id)+1))
                end if
             end do
             do i = 1, grid(d_src)%pack(AT_EDGE,dest)%length
                id = grid(d_src)%pack(AT_EDGE,dest)%elts(i)
                if (l == NONE .or. l == grid(d_src)%level%elts(id/EDGE+1)) then
                   call append (send_buf_i, grid(d_src)%mask_e%elts(abs(id)+1))
                end if
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    ! Determine recv buff lengths
    recv_buf_i%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf_i%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = abs(grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i))
                if (l == NONE .or. l == grid(d_dest)%level%elts(id+1)) recv_buf_i%length = recv_buf_i%length + 1
             end do
             do i = 1, grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%length
                id = abs (grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%elts(i))
                if (l == NONE .or. l == grid(d_dest)%level%elts(id/EDGE+1)) recv_buf_i%length = recv_buf_i%length + 1
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf_i%length - recv_offsets(r_src)
    end do
    if (size(recv_buf_i%elts) < recv_buf_i%length) then
       deallocate (recv_buf_i%elts)
       allocate (recv_buf_i%elts(recv_buf_i%length))
       recv_buf_i%elts = 0
    end if

    ! Call check_alltoall_lengths()
    call MPI_Alltoallv (send_buf_i%elts, send_lengths, send_offsets, MPI_INTEGER, &
                        recv_buf_i%elts, recv_lengths, recv_offsets, MPI_INTEGER, &
                        MPI_COMM_WORLD, ierror)

    ! Communicate inside domain
    call comm_masks

    kk = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                if (l == NONE .or. l == grid(d_dest)%level%elts(abs(id)+1)) then
                   kk = kk + 1
                   grid(d_dest)%mask_n%elts(abs(id)+1) = recv_buf_i%elts(kk)
                end if
             end do
             do i = 1, grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%elts(i)
                if (l == NONE .or. l == grid(d_dest)%level%elts(abs(id)/EDGE+1)) then
                   kk = kk + 1
                   grid(d_dest)%mask_e%elts(abs(id)+1) = recv_buf_i%elts(kk)
                end if
             end do
          end do
       end do
    end do
  end subroutine comm_masks_mpi

  subroutine deadlock_test (flag)
    ! Detects deadlock and aborts
    implicit none
    integer :: flag

    logical            :: got_data
    real(8)            :: t_start
    real(8), parameter :: timeout_time = 1d2

    t_start = MPI_Wtime () 
    call MPI_Test (req, got_data, MPI_STATUS_IGNORE, ierror)
    do while (.not. got_data .and. MPI_Wtime()-t_start < timeout_time) 
       call MPI_Test (req, got_data, MPI_STATUS_IGNORE, ierror)
    end do
    if (.not. got_data) then
       write (6,'(a,i4,a,i5)') "ERROR: boundary update deadlocked at call ", flag, " on rank ", rank
       call abort
    end if
  end subroutine deadlock_test

  subroutine update_bdry1_0 (field, l_start, l_end, flag_in)
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    integer, optional :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start1  (field, l_start, l_end)
    if (deadlock) call deadlock_test (flag_in)
    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_0

  subroutine update_bdry1_1 (field, l_start, l_end, flag_in)
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l_start, l_end
    integer, optional               :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start1  (field, l_start, l_end)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_1

  subroutine update_bdry1_2 (field, l_start, l_end, flag_in)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l_start, l_end
    integer, optional                 :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start1  (field, l_start, l_end)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1_2

  subroutine update_bdry_0 (field, l, flag_in)
    implicit none
    type(Float_Field) :: field
    integer           :: l
    integer, optional :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start  (field, l)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry_0

  subroutine update_bdry_1 (field, l, flag_in)
    implicit none
    ! Updates field array
    type(Float_Field), dimension(:) :: field
    integer                         :: l
    integer, optional               :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start  (field, l)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry_1
  
  subroutine update_bdry_2 (field, l, flag_in)
    ! Updates field array
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l
    integer, optional                 :: flag_in

    integer :: flag

    if (present (flag_in)) then
       flag = flag_in
    else
       flag = 9999
    end if

    call update_bdry__start  (field, l)
    if (deadlock) call deadlock_test (flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry_2

  subroutine update_bdry__start_0 (field, l)
    implicit none
    type(Float_Field) :: field
    integer           ::  l

    if (l == NONE) then 
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    endif
  end subroutine update_bdry__start_0
  
  subroutine update_bdry__start_1 (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l == NONE) then 
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    endif
  end subroutine update_bdry__start_1
  
  subroutine update_bdry__start_2 (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l == NONE) then 
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    endif
  end subroutine update_bdry__start_2

  subroutine update_bdry__start1_0 (field, l_start, l_end)
    use mpi
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: d_src, d_dest, dest, i, id, k, lev, multipl, r, r_dest, r_src

    if (field%bdry_uptodate) return

    if (field%pos == AT_EDGE) then
       multipl = EDGE
    else
       multipl = 1
    end if

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest) + 1
             do i = 1, grid(d_src)%pack(field%pos,dest)%length
                id = grid(d_src)%pack(field%pos,dest)%elts(i)
                lev = grid(d_src)%level%elts(id/multipl+1)
                if (lev >= l_start .and. lev <= l_end) call append (send_buf, field%data(d_src)%elts(id+1))
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    ! Determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%length
                id = abs (grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%elts(i))
                lev = grid(d_dest)%level%elts(id/multipl+1)
                if (lev >= l_start .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
       recv_buf%elts = 0d0
    end if

    ! Post all receives first so buffer is available
    nreq = 0
    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Irecv (recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! Communicate inside domain
    call cp_bdry_inside (field)
  end subroutine update_bdry__start1_0

  subroutine update_bdry__start1_1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    use mpi
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l_start, l_end
    
    integer :: d_dest, d_src, dest, i, i1, id, k, lev, multipl, pos, r, r_dest, r_src
    logical :: ret

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i1 = 1, size(field)
       if (.not. field(i1)%bdry_uptodate) ret=.false.
    end do
    if (ret) return

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest) + 1
             ! Loop over each element of field array
             do i1 = 1, size(field)
                pos = field(i1)%pos
                if (pos == AT_EDGE) then
                   multipl = EDGE
                else
                   multipl = 1
                end if
                do i = 1, grid(d_src)%pack(pos,dest)%length
                   id = grid(d_src)%pack(pos,dest)%elts(i)
                   lev = grid(d_src)%level%elts(id/multipl+1)
                   if (lev >= l_start .and. lev <= l_end) call append (send_buf, field(i1)%data(d_src)%elts(id+1))
                end do
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    ! Determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             ! Loop over each element of field array
             do i1 = 1, size(field)
                pos = field(i1)%pos
                if (pos == AT_EDGE) then
                   multipl = EDGE
                else
                   multipl = 1
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))
                   lev = grid(d_dest)%level%elts(id/multipl+1)
                   if (lev >= l_start .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
                end do
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
       recv_buf%elts = 0d0
    end if

    ! Post all receives first so buffer is available
    nreq = 0
    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Irecv (recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! Communicate inside domain
    call cp_bdry_inside_vector (field)
  end subroutine update_bdry__start1_1

  subroutine update_bdry__start1_2 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    use mpi
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           ::  l_start, l_end
    
    integer :: d_dest, d_src, dest, i, i1, i2, id, k, lev, multipl, pos, r, r_dest, r_src
    logical :: ret

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i2 = 1, size(field,2)
       do i1 = 1, size(field,1)
          if (.not. field(i1,i2)%bdry_uptodate) ret = .false.
       end do
    end do
    if (ret) return

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest) + 1
             ! Loop over each element of field array
             do i2 = 1, size(field,2)
                do i1 = 1, size(field,1)
                   pos = field(i1,i2)%pos
                   if (pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if
                   do i = 1, grid(d_src)%pack(pos,dest)%length
                      id = grid(d_src)%pack(pos,dest)%elts(i)
                      lev = grid(d_src)%level%elts(id/multipl+1)
                      if (lev >= l_start .and. lev <= l_end) call append (send_buf, field(i1,i2)%data(d_src)%elts(id+1))
                   end do
                end do
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    ! Determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             ! Loop over each element of field array
             do i2 = 1, size(field,2)
                do i1 = 1, size(field,1)
                   pos = field(i1,i2)%pos
                   if (field(i1,i2)%pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))
                      lev = grid(d_dest)%level%elts(id/multipl+1)
                      if (lev >= l_start .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
                   end do
                end do
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
       recv_buf%elts = 0d0
    end if

    ! Post all receives first so buffer is available
    nreq = 0
    do r = 1, n_process
       if (r == rank+1 .or. recv_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Irecv (recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r == rank+1 .or. send_lengths(r) == 0) cycle
       nreq = nreq + 1
       call MPI_Isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, r-1, 1, &
            MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! Communicate inside domain
    call cp_bdry_inside_array (field)
  end subroutine update_bdry__start1_2

  subroutine update_bdry__finish_0 (field, l)
    implicit none
    type(Float_Field) :: field
    integer           :: l

    if (l == NONE) then 
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    endif
  end subroutine update_bdry__finish_0
  
  subroutine update_bdry__finish_1 (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l == NONE) then 
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    endif
  end subroutine update_bdry__finish_1
  
  subroutine update_bdry__finish_2 (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l == NONE) then 
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    endif
  end subroutine update_bdry__finish_2

  subroutine update_bdry__finish1_0 (field, l_start, l_end)
    use mpi
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: r_dest, r_src, d_src, d_dest, dest, id, i, k, multipl, lev

    if (field%bdry_uptodate) return

    if (field%pos == AT_EDGE) then
       multipl = EDGE
    else
       multipl = 1
    end if

    call MPI_Waitall (nreq, req, MPI_STATUSES_IGNORE, ierror)

    k = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%elts(i)
                lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                if (lev >= l_start .and. lev <= l_end) then
                   k = k + 1
                   field%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                   if (id < 0 .and. field%pos == AT_EDGE) field%data(d_dest)%elts(abs(id)+1) = -field%data(d_dest)%elts(abs(id)+1)
                end if
             end do
          end do
       end do
    end do

    ! Assumes routine is either called for one level, or all levels ever to be updated
    if (l_start < l_end) field%bdry_uptodate = .true.
  end subroutine update_bdry__finish1_0

  subroutine update_bdry__finish1_1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    use mpi
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l_start, l_end
    
    integer :: d_dest, d_src, dest, id, i, i1, k, lev, multipl, pos, r_dest, r_src
    logical :: ret

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i1 = 1, size(field)
       if (.not. field(i1)%bdry_uptodate) ret=.false.
    end do
    if (ret) return

    call MPI_Waitall (nreq, req, MPI_STATUSES_IGNORE, ierror)
    
    k = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i1 = 1, size(field)
                pos = field(i1)%pos 
                if (pos == AT_EDGE) then
                   multipl = EDGE
                else
                   multipl = 1
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)
                   lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                   if (lev >= l_start .and. lev <= l_end) then
                      k = k + 1
                      field(i1)%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                      if (id < 0 .and. pos == AT_EDGE) &
                           field(i1)%data(d_dest)%elts(abs(id)+1) = -field(i1)%data(d_dest)%elts(abs(id)+1)
                   end if
                end do
             end do
          end do
       end do
    end do

    ! Assumes routine is either called for one level, or all levels ever to be updated
    if (l_start < l_end) field%bdry_uptodate = .true.
  end subroutine update_bdry__finish1_1
  
  subroutine update_bdry__finish1_2 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    use mpi
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l_start, l_end
    
    integer :: d_dest, d_src, dest, i, i1, i2, id, k, lev, multipl, pos, r_dest, r_src
    logical :: ret

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i2 = 1, size(field,2)
       do i1 = 1, size(field,1)
          if (.not. field(i1,i2)%bdry_uptodate) ret=.false.
       end do
    end do
    if (ret) return

    call MPI_Waitall (nreq, req, MPI_STATUSES_IGNORE, ierror)

    k = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i2 = 1, size(field,2)
                do i1 = 1, size(field,1)
                   pos = field(i1,i2)%pos
                   if (pos == AT_EDGE) then
                      multipl = EDGE
                   else
                      multipl = 1
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)
                      lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                      if (lev >= l_start .and. lev <= l_end) then
                         k = k + 1
                         field(i1,i2)%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                         if (id < 0 .and. pos == AT_EDGE) &
                              field(i1,i2)%data(d_dest)%elts(abs(id)+1) = -field(i1,i2)%data(d_dest)%elts(abs(id)+1)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

    ! Assumes routine is either called for one level, or all levels ever to be updated
    if (l_start < l_end) field%bdry_uptodate = .true.
  end subroutine update_bdry__finish1_2
  
  subroutine comm_nodes9_mpi (get, set, l)
    use mpi
    implicit none
    external :: get, set
    integer :: l
    
    real(8), dimension(7) :: val
    integer               :: r_dest, r_src, d_src, d_dest, dest, id, i, k

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle 
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                id = grid(d_src)%pack(AT_NODE,dest)%elts(i)
                k = send_buf%length
                call extend (send_buf, 7, 0d0)
                call get (grid(d_src), id, val)
                send_buf%elts(k+1:k+7) = val
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    ! determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                recv_buf%length = recv_buf%length + 7
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
       recv_buf%elts = 0d0
    end if

    call MPI_Alltoallv (send_buf%elts, send_lengths, send_offsets, MPI_DOUBLE_PRECISION, &
                        recv_buf%elts, recv_lengths, recv_offsets, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, ierror)

    call comm_nodes9 (get, set) ! communicate inside domain

    k = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                call set (grid(d_dest), id, recv_buf%elts(k+1:k+7))
                k = k + 7
             end do
          end do
       end do
    end do
  end subroutine comm_nodes9_mpi

  subroutine comm_nodes3_mpi (get, set, l)
    use mpi
    implicit none
    external    :: get, set
    type(Coord) :: get
    integer     :: l
    
    integer     :: r_dest, r_src, d_src, d_dest, dest, id, i, k
    type(Coord) :: c

    send_buf%length = 0 ! reset
    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                id = grid(d_src)%pack(AT_NODE,dest)%elts(i)
                c = get(grid(d_src), id)
                k = send_buf%length
                call extend (send_buf, 3, 0d0)
                send_buf%elts(k+1:k+3) = (/c%x, c%y, c%z/)
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf%length - send_offsets(r_dest)
    end do

    ! Determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src == rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                recv_buf%length = recv_buf%length + 3
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
       recv_buf%elts = 0d0
    end if

    call MPI_Alltoallv (send_buf%elts, send_lengths, send_offsets, MPI_DOUBLE_PRECISION, &
                        recv_buf%elts, recv_lengths, recv_offsets, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, ierror)

    call comm_nodes3 (get, set) ! communicate inside domain

    k = 0
    do r_src = 1, n_process 
       if (r_src == rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                c%x = recv_buf%elts(k+1) 
                c%y = recv_buf%elts(k+2) 
                c%z = recv_buf%elts(k+3) 
                call set (grid(d_dest), id, c)
                k = k + 3
             end do
          end do
       end do
    end do
  end subroutine comm_nodes3_mpi

  subroutine comm_patch_conn_mpi
    use mpi
    implicit none
    integer               :: r_src, r_dest, d_src, d_dest, i, b, c, p, s, d_glo, k, rot, d, ngh_pa, typ, l_par, rot_shift
    integer, dimension(4) :: st
    logical               :: is_pole
    
    send_buf_i%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf_i%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) exit ! inside domain
          do d_dest = 1, n_domain(r_dest)
             do i = 1, grid(d_src)%send_pa_all%length, 4
                b = grid(d_src)%send_pa_all%elts(0+i)
                c = grid(d_src)%send_pa_all%elts(1+i)
                p = grid(d_src)%send_pa_all%elts(2+i)
                s = grid(d_src)%send_pa_all%elts(3+i)
                typ = grid(d_src)%bdry_patch%elts(b+1)%side
                d_glo = grid(d_src)%neigh(typ) ! 0 ...
                is_pole = d_glo == POLE
                if (is_pole) then
                   d_glo = grid(d_src)%neigh_over_pole(c+1)
                   l_par = grid(d_src)%patch%elts(p+1)%level - 1
                   if (grid(d_src)%neigh_pa_over_pole%length < l_par*2 + c + 1) then
                      ngh_pa = 0
                   else
                      ngh_pa = grid(d_src)%neigh_pa_over_pole%elts(l_par*2 + c + 1)
                   end if
                else
                   ngh_pa = grid(d_src)%bdry_patch%elts(b+1)%neigh
                end if
                ! Also skips if dest == 0
                if (ngh_pa /= 0 .and. d_glo == glo_id(r_dest,d_dest)) then 
                   rot = grid(d_src)%neigh_rot(typ)
                   rot_shift = (rot_direction(grid(d_src), typ)*2 - 1)*rot
                   call append (send_buf_i, d_dest)
                   call append (send_buf_i, glo_id(rank+1,d_src))
                   call append (send_buf_i, p)
                   if (is_pole) then
                      call append (send_buf_i, c)
                   else
                      call append (send_buf_i, modulo(c + rot_shift, 4))
                   end if
                   call append (send_buf_i, ngh_pa)
                   if (is_pole) then
                      call append (send_buf_i, s)
                   else
                      call append (send_buf_i, modulo(s + rot_shift + 2, 4) + 4*(s/4))
                   end if
                end if
             end do
          end do
       end do
       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    call alltoall

    call comm_patch_conn

    do r_src = 1, n_process
       if (r_src == rank+1) cycle ! inside domain
       do k = recv_offsets(r_src) + 1, recv_offsets(r_src) + recv_lengths(r_src), 6
          d = recv_buf_i%elts(k)
          d_src = recv_buf_i%elts(k+1)+1
          st = recv_buf_i%elts(k+2:k+5)
          call append (grid(d)%recv_pa(d_src), st(1))
          call append (grid(d)%recv_pa(d_src), st(2))
          call append (grid(d)%recv_pa(d_src), st(3))
          call append (grid(d)%recv_pa(d_src), st(4))
       end do
    end do
  end subroutine comm_patch_conn_mpi

  integer function sync_max_int (val)
    use mpi
    implicit none
    integer :: val
    
    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_int = val_glo
  end function sync_max_int

  integer function sync_min_int (val)
    use mpi
    implicit none
    integer :: val
    
    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_MIN, MPI_COMM_WORLD, ierror)
    sync_min_int = val_glo
  end function sync_min_int

  real(8) function sync_max_real_0 (val)
    use mpi
    implicit none
    real(8) :: val

    real(8) :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_real_0 = val_glo
  end function sync_max_real_0

  function sync_max_real_1 (val)
    use mpi
    implicit none
    real(8), dimension(:), allocatable :: sync_max_real_1
    real(8), dimension(:)              :: val

    integer                            :: n
    real(8), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (sync_max_real_1(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_real_1 = val_glo
  end function sync_max_real_1

  real(8) function sync_min_real_0 (val)
    use mpi
    implicit none
    real(8) :: val

    real(8) :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    sync_min_real_0 = val_glo
  end function sync_min_real_0

  function sync_min_real_1 (val)
    use mpi
    implicit none
    real(8), dimension(:), allocatable :: sync_min_real_1
    real(8), dimension(:)              :: val

    integer                            :: n
    real(8), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (sync_min_real_1(n), val_glo(n))

    call MPI_Allreduce (val, val_glo, n, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    sync_min_real_1 = val_glo
  end function sync_min_real_1

  real(8) function sum_real_0 (val)
    use mpi
    implicit none
    real(8) :: val

    real(8) :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_real_0 = val_glo
  end function sum_real_0

  function sum_real_1 (val)
    use mpi
    implicit none
    real(8), dimension(:), allocatable :: sum_real_1
    real(8), dimension(:)              :: val

    integer                            :: n
    real(8), dimension(:), allocatable :: val_glo

    n = size (val,1)
    allocate (sum_real_1(n), val_glo(n))
    
    call MPI_Allreduce (val, val_glo, n, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_real_1 = val_glo
  end function sum_real_1

  integer function sum_int (val)
    use mpi
    implicit none
    integer :: val

    integer :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_int = val_glo
  end function sum_int

  function sum_int_vector (val, n)
    use mpi
    implicit none
    integer, dimension(n) :: sum_int_vector
    integer               :: n
    integer, dimension(n) :: val

    integer, dimension(n) :: val_glo
    
    call MPI_Allreduce (val, val_glo, n, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_int_vector = val_glo
  end function sum_int_vector

  subroutine start_timing
    use mpi
    implicit none
    times(1) = MPI_Wtime()  
  end subroutine start_timing

  subroutine stop_timing
    use mpi
    implicit none
    times(2) = MPI_Wtime()  
  end subroutine stop_timing

  real(8) function get_timing()
    use mpi
    implicit none
    get_timing = times(2) - times(1)
  end function get_timing

  subroutine sync_array (arr, N)
    ! Synchronizes arr, with subsections computed on different ranks
    ! the complete array is broadcast to all ranks
    ! (see project_field_onto_plane for example of use)
    use mpi
    implicit none
    real(8), dimension(N) :: arr
    integer               :: N
    
    integer               :: myop
    real(8), dimension(N) :: garr

    call MPI_Op_create (sync, .true., myop, ierror)  
    call MPI_Reduce (arr, garr, N, MPI_DOUBLE_PRECISION, myop, 0, MPI_COMM_WORLD, ierror)
    if (rank == 0) arr(1:N) = garr(1:N)
    call MPI_Bcast (arr, N, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror) 
  end subroutine sync_array

  subroutine sync (in, inout, len, type)
    use mpi
    implicit none
    real(8), dimension(len) :: in, inout
    integer                 :: len, type
  
    where (in /= sync_val) inout = in
  end subroutine sync
end module comm_mpi_mod
