module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
  integer                              :: nreq
  type(Int_Array)                      :: recv_buf_i, send_buf_i
  type(Float_Array)                    :: recv_buf, send_buf
  integer, dimension(:),   allocatable :: recv_lengths, recv_offsets, req, send_lengths, send_offsets
  real(8), dimension(2)                :: times
contains
  subroutine init_comm_mpi
    implicit none
    allocate (send_lengths(n_process), send_offsets(n_process))
    allocate (recv_lengths(n_process), recv_offsets(n_process))
    allocate (req(2*n_process))
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  integer function write_active_per_level ()
    ! Write out distribution of active nodes over levels
    implicit none
    integer                                         :: l, n_full, fillin, n_lev_cur, recommended_level_start
    integer, dimension(2*(level_end-level_start+1)) :: n_active_all_loc, n_active_all_glo
    integer, dimension(level_start:level_end)       :: n_active_per_lev
    real(8)                                         :: dt

    dt = cpt_dt() ! to set n_active_*

    n_lev_cur = level_end - level_start + 1

    n_active_all_loc = (/n_active_nodes(level_start:level_end), n_active_edges(level_start:level_end)/)
    
    ! Sum n_active_all_loc up across all processes and distribute result n_active_all_glo among all processes
    call MPI_Allreduce (n_active_all_loc, n_active_all_glo, n_lev_cur*2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    
    n_active_nodes(level_start:level_end) = n_active_all_glo(1:n_lev_cur)
    n_active_edges(level_start:level_end) = n_active_all_glo(n_lev_cur+1:n_lev_cur*2)
    n_active_per_lev = n_active_edges(level_start:level_end) + n_active_nodes(level_start:level_end)

    if (rank == 0) write (6,'(6X,A,A,3(1X,A))') '   N_p   ', '   N_u   ','of all active', 'of full level', 'fill-in'

    recommended_level_start = level_start

    do l = level_start, level_end
       n_full = max_nodes_per_level(l) + max_nodes_per_level(l,EDGE)

       ! Fill-in: additional nodes on level `l` if it'd become lowest level 
       ! minus the nodes on lower levels which would be removed
       fillin = n_full-n_active_per_lev(l)-sum(n_active_per_lev(level_start:l-1))

       if (rank == 0) then
          write (6,'(A,I2,I9,I9,2(1X,F9.1,A),1X,I9,1X,F9.1,A)') &
               'lev', l, n_active_nodes(l), n_active_edges(l), &
               float(n_active_per_lev(l))/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%', &
               float(n_active_per_lev(l))/float(n_full)*100.0, '%', &
               fillin, float(fillin)/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%'
       end if

       if (fillin <= 0) recommended_level_start = l
    end do

    if (rank == 0) then
       write (6,'(A,I9,I9,2(1X,F9.1,A),9X,I9)') 'total', n_active(AT_NODE:AT_EDGE), 100.0, '%', &
            float(sum(n_active(AT_NODE:AT_EDGE)))/float(n_full)*100.0, '%', &
            n_full/sum(n_active(AT_NODE:AT_EDGE))
    end if

    write_active_per_level = recommended_level_start
  end function write_active_per_level

  subroutine write_load_conn (id, run_id)
    ! Write out load distribution and connectivity for load balancing
    implicit none
    integer      :: id
    character(*) :: run_id
    
    integer        :: r, fid
    character(255) :: filename

    fid = 599
    write (filename, '(A,A,I4.4)')  trim (run_id), "_conn.", id

    do r = 1, n_process
       if (r /= rank+1) then ! write only if our turn, otherwise only wait at Barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if

       if (r == 1) then ! first process opens without append to delete old file if existing
          open (unit=fid, file=trim(filename), recl=333333)
       else
          open (unit=fid, file=trim(filename), recl=333333, access='APPEND')
       end if

       call write_load_conn1 (fid)
       close(fid)
       call MPI_Barrier (MPI_Comm_World, ierror)
    end do
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

  subroutine write_level_mpi (out_rout, fid, l, zlev, eval_pole, run_id)
    implicit none
    external       :: out_rout
    integer        :: fid, l, zlev
    logical        :: eval_pole
    character(*)   :: run_id

    integer            :: r
    integer, parameter :: funit = 300
    character(7)       :: var_file
    character(255)     :: filename
        
    do r = 1, n_process
       if (r /= rank+1) then ! Write only if our turn, otherwise only wait at barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if

       write (var_file, '(I7)')  fid
       filename = trim(run_id)//'.'//var_file
       if (r == 1) then ! first process opens without append to delete old file if existing
          open (unit=funit, file=filename)
       else
          open (unit=funit, file=filename, access='APPEND')
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
             dest = glo_id(r_dest,d_dest)+1
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
    implicit none
    integer, dimension(n_process) :: test_recv_len

    call MPI_Alltoall (send_lengths, 1, MPI_INTEGER, test_recv_len, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)

    write (3000+rank,*) test_recv_len-recv_lengths
    close (3000+rank)
    call MPI_Barrier (MPI_Comm_World, ierror)
  end subroutine check_alltoall_lengths

  subroutine alltoall
    implicit none
    integer :: i

    call MPI_Alltoall (send_lengths, 1, MPI_INTEGER, recv_lengths, 1, MPI_INTEGER, MPI_COMM_WORLD, ierror)

    recv_offsets(1) = 0

    do i = 2, n_process
       recv_offsets(i) = recv_offsets(i-1) + recv_lengths(i-1)
    end do

    recv_buf_i%length = sum(recv_lengths)

    if (size(recv_buf_i%elts) < recv_buf_i%length) then
       deallocate (recv_buf_i%elts)
       allocate (recv_buf_i%elts(recv_buf_i%length))
    end if

    call MPI_Alltoallv (send_buf_i%elts, send_lengths, send_offsets, MPI_INTEGER, recv_buf_i%elts, &
         recv_lengths, recv_offsets, MPI_INTEGER, MPI_COMM_WORLD, ierror)
  end subroutine alltoall

  subroutine comm_masks_mpi (l)
    ! Communication of mask information in a subdomain between different processes
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
    integer                         :: flag

    logical            :: got_data
    real(8)            :: t_start
    real(8), parameter :: timeout_time = 1.0d2

    t_start = MPI_Wtime() 
    call MPI_Test (req, got_data, MPI_STATUS_IGNORE, ierror)
    do while (.not. got_data .and. MPI_Wtime()-t_start < timeout_time) 
       call MPI_Test (req, got_data, MPI_STATUS_IGNORE, ierror)
    end do
    if (.not. got_data) then
       write (6,'(A,i2,A,i5)') "ERROR: boundary update deadlocked at call ", flag, " on rank ", rank
       call abort
    end if
  end subroutine deadlock_test

  subroutine update_bdry1 (field, l_start, l_end, flag)
    implicit none
    type(Float_Field) :: field
    integer           :: flag, l_start, l_end

    call update_bdry__start1  (field, l_start, l_end)
    call deadlock_test (flag)
    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1

  subroutine update_array_bdry1 (field, l_start, l_end, flag)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: flag, l_start, l_end

    call update_array_bdry__start1  (field, l_start, l_end)
    call deadlock_test (flag)
    call update_array_bdry__finish1 (field, l_start, l_end)
  end subroutine update_array_bdry1

  subroutine update_bdry (field, l, flag)
    implicit none
    type(Float_Field) :: field
    integer           :: flag, l

    call update_bdry__start  (field, l)
    call deadlock_test (flag)
    call update_bdry__finish (field, l)
  end subroutine update_bdry

  subroutine update_vector_bdry (field, l, flag)
    implicit none
    ! Updates field array
    type(Float_Field), dimension(:) :: field
    integer                         :: flag, l

    call update_vector_bdry__start  (field, l)
    call deadlock_test (flag)
    call update_vector_bdry__finish (field, l)
  end subroutine update_vector_bdry
  
  subroutine update_array_bdry (field, l, flag)
    ! Updates field array
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: flag, l

    call update_array_bdry__start  (field, l)
    call deadlock_test (flag)
    call update_array_bdry__finish (field, l)
  end subroutine update_array_bdry

  subroutine update_bdry__start (field, l)
    implicit none
    type(Float_Field) :: field
    integer           ::  l

    if (l == NONE) then 
       call update_bdry__start1 (field, level_start-1, level_end)
    else
       call update_bdry__start1 (field, l, l)
    endif
  end subroutine update_bdry__start
  
  subroutine update_vector_bdry__start (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l == NONE) then 
       call update_vector_bdry__start1 (field, level_start-1, level_end)
    else
       call update_vector_bdry__start1 (field, l, l)
    endif
  end subroutine update_vector_bdry__start
  
  subroutine update_array_bdry__start (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l == NONE) then 
       call update_array_bdry__start1 (field, level_start-1, level_end)
    else
       call update_array_bdry__start1 (field, l, l)
    endif
  end subroutine update_array_bdry__start

  subroutine update_bdry__start1 (field, l_start, l_end)
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: d_src, d_dest, dest, i, id, k, lev, multipl, r, r_dest, r_src

    if (field%bdry_uptodate) return

    if (field%pos == AT_NODE) then
       multipl = 1
    else
       multipl = EDGE
    end if

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             do i = 1, grid(d_src)%pack(field%pos,dest)%length
                id = grid(d_src)%pack(field%pos,dest)%elts(i)
                lev = grid(d_src)%level%elts(id/multipl+1)
                if (l_start <= lev .and. lev <= l_end) call append (send_buf, field%data(d_src)%elts(id+1))
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
                if (l_start <= lev .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
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
  end subroutine update_bdry__start1

  subroutine update_vector_bdry__start1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
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
             dest = glo_id(r_dest,d_dest)+1
             ! Loop over each element of field array
             do i1 = 1, size(field)
                pos = field(i1)%pos
                if (pos == AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_src)%pack(pos,dest)%length
                   id = grid(d_src)%pack(pos,dest)%elts(i)
                   lev = grid(d_src)%level%elts(id/multipl+1)
                   if (l_start <= lev .and. lev <= l_end) call append (send_buf, field(i1)%data(d_src)%elts(id+1))
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
                if (pos == AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))
                   lev = grid(d_dest)%level%elts(id/multipl+1)
                   if (l_start <= lev .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
                end do
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) < recv_buf%length) then
       deallocate (recv_buf%elts)
       allocate (recv_buf%elts(recv_buf%length))
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
  end subroutine update_vector_bdry__start1

  subroutine update_array_bdry__start1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           ::  l_start, l_end
    
    integer :: d_dest, d_src, dest, i, i1, i2, id, k, lev, multipl, pos, r, r_dest, r_src
    logical :: ret

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i2 = 1, size(field,2)
       do i1 = 1, size(field,1)
          if (.not. field(i1,i2)%bdry_uptodate) ret=.false.
       end do
    end do
    if (ret) return

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest == rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             ! Loop over each element of field array
             do i2 = 1, size(field,2)
                do i1 = 1, size(field,1)
                   pos = field(i1,i2)%pos
                   if (pos == AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_src)%pack(pos,dest)%length
                      id = grid(d_src)%pack(pos,dest)%elts(i)
                      lev = grid(d_src)%level%elts(id/multipl+1)
                      if (l_start <= lev .and. lev <= l_end) call append (send_buf, field(i1,i2)%data(d_src)%elts(id+1))
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
                   if (field(i1,i2)%pos == AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))
                      lev = grid(d_dest)%level%elts(id/multipl+1)
                      if (l_start <= lev .and. lev <= l_end) recv_buf%length = recv_buf%length + 1
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
  end subroutine update_array_bdry__start1

  subroutine update_bdry__finish (field, l)
    implicit none
    type(Float_Field) :: field
    integer           :: l

    if (l == NONE) then 
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    endif
  end subroutine update_bdry__finish
  
  subroutine update_vector_bdry__finish (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l == NONE) then 
       call update_vector_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_vector_bdry__finish1 (field, l, l)
    endif
  end subroutine update_vector_bdry__finish
  
  subroutine update_array_bdry__finish (field, l)
    ! Finishes boundary update for field arrays
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l == NONE) then 
       call update_array_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_array_bdry__finish1 (field, l, l)
    endif
  end subroutine update_array_bdry__finish

  subroutine update_bdry__finish1 (field, l_start, l_end)
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: r_dest, r_src, d_src, d_dest, dest, id, i, k, multipl, lev

    if (field%bdry_uptodate) return

    if (field%pos == AT_NODE) then
       multipl = 1
    else
       multipl = EDGE
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
                if (l_start <= lev .and. lev <= l_end) then
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
  end subroutine update_bdry__finish1

  subroutine update_vector_bdry__finish1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
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
                if (pos == AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)
                   lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                   if (l_start <= lev .and. lev <= l_end) then
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
  end subroutine update_vector_bdry__finish1
  
  subroutine update_array_bdry__finish1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
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
                   if (pos == AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)
                      lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                      if (l_start <= lev .and. lev <= l_end) then
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
  end subroutine update_array_bdry__finish1
  
  subroutine comm_nodes9_mpi (get, set, l)
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
                call extend (send_buf, 7, 0._8)
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
                call extend (send_buf, 3, 0._8)
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

  real(8) function cpt_dt ()
    ! Calculates time step, minimum relative mass and active nodes and edges
    implicit none
    integer               :: l, ierror
    integer, dimension(2) :: n_active_loc
    
    if (adapt_dt) then
       dt_loc = 1d16
    else
       dt_loc = dt_init
    end if
    n_active_nodes = 0
    n_active_edges = 0

    ! Calculate minimum time step, number of active nodes and edges
    do l = level_start, level_end
       call apply_onescale (min_dt, l, z_null, 0, 0)
    end do

    ! Time step
    if (adapt_dt) then
       call MPI_Allreduce (dt_loc, cpt_dt, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    else
       cpt_dt = dt_loc
    end if

    ! Active nodes and edges
    n_active_loc = (/ sum(n_active_nodes(level_start:level_end)), sum(n_active_edges(level_start:level_end)) /)
    call MPI_Allreduce (n_active_loc, n_active,  2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    call MPI_Allreduce (MPI_IN_PLACE, level_end, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
  end function cpt_dt

  real(8) function cpt_min_mass ()
    ! Calculates minimum relative mass
    implicit none
    integer :: ierror, l

    min_mass_loc = 1d16
    do l = level_start, level_end
       call apply_onescale (cal_min_mass, l, z_null, 0, 0)
    end do

    call MPI_Allreduce (min_mass_loc, cpt_min_mass, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
  end function cpt_min_mass

  integer function sync_max_int (val)
    implicit none
    integer :: val
    
    integer :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_int = val_glo
  end function sync_max_int

  real(8) function sync_max_real (val)
    implicit none
    real(8) :: val

    real(8) :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_real = val_glo
  end function sync_max_real

  real(8) function sum_real (val)
    implicit none
    real(8) :: val

    real(8) :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_real = val_glo
  end function sum_real

  integer function sum_int (val)
    implicit none
    integer :: val

    integer :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_int = val_glo
  end function sum_int

  subroutine start_timing
    implicit none
    times(1) = MPI_Wtime()  
  end subroutine start_timing

  subroutine stop_timing
    implicit none
    times(2) = MPI_Wtime()  
  end subroutine stop_timing

  real(8) function get_timing()
    implicit none
    get_timing = times(2) - times(1)
  end function get_timing

  subroutine sync (in, inout, len, type)
    implicit none
    real, dimension(len) :: in, inout
    integer              :: len, type

    where (in /= sync_val) inout = in
  end subroutine sync

  subroutine sync_array (arr, N)
    implicit none
    real, dimension(N) :: arr
    integer            :: N
    
    integer            :: myop
    real, dimension(N) :: garr

    call MPI_Op_create (sync, .true., myop, ierror)  
    call MPI_Reduce (arr, garr, N, MPI_REAL, myop, 0, MPI_COMM_WORLD, ierror)
    if (rank == 0) arr(1:N) = garr(1:N)
  end subroutine sync_array

  subroutine stop_and_record_timings (id)
    ! Use like:
    ! call stop_and_record_timings(6500); call start_timing()
    ! call stop_and_record_timings(6501); call start_timing()
    implicit none
    integer :: id

    real(8) :: time_loc, time_max, time_min, time_sum

    call stop_timing

    time_loc = get_timing()

    call MPI_Reduce (time_loc, time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (time_loc, time_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (time_loc, time_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    if (rank == 0) write (id,'(3(es8.2,1x))') time_max, time_min, time_sum
  end subroutine stop_and_record_timings

  subroutine combine_stats
    ! Uses Chan, Golub and LeVeque (1983) algorithm for partitioned data sets to combine zonal average results from each rank
    !   T.F. Chan, G.H. Golub & R.J. LeVeque (1983):
    !   "Algorithms for computing the sample variance: Analysis and recommendations." The American Statistician 37: 242–247.
    implicit none
    integer                                      :: bin, ivar, k, r
    integer, dimension(MPI_STATUS_SIZE)          :: status
    integer, dimension(zlevels,nbins)            :: Nstats_loc
    real(8), dimension(zlevels,nbins,nvar_zonal) :: zonal_avg_loc

    ! Initialize to values on rank 0
    if (rank == 0) then
       Nstats_glo    = Nstats
       zonal_avg_glo = zonal_avg
    end if
       
    do r = 1, n_process - 1
       if (rank == r) then
          call MPI_Send (Nstats,     zlevels*nbins,            MPI_INT,              0, 0, MPI_COMM_WORLD, ierror)
          call MPI_Send (zonal_avg,  zlevels*nbins*nvar_zonal, MPI_DOUBLE_PRECISION, 0, 0, MPI_COMM_WORLD, ierror)
       elseif (rank == 0) then
          call MPI_Recv (Nstats_loc,     zlevels*nbins,            MPI_INT,              r, 0, MPI_COMM_WORLD, status, ierror)
          call MPI_Recv (zonal_avg_loc,  zlevels*nbins*nvar_zonal, MPI_DOUBLE_PRECISION, r, 0, MPI_COMM_WORLD, status, ierror)

          do k = 1, zlevels
             do bin = 1, nbins
                if (Nstats_loc(k,bin) /= 0) call combine_var
             end do
          end do
       end if
    end do
  contains
    subroutine combine_var
      integer :: nA, nB, nAB
      real(8) :: delta_KE, delta_T, delta_U, delta_V
      real(8) :: mA_T, mB_T, mA_U, mB_U, mA_V, mB_V, mA_zonal, mB_zonal, mA_merid, mB_merid, mA_VT, mB_VT, mA_UV, mB_UV

      nA = Nstats_glo(k,bin)
      nB = Nstats_loc(k,bin)
      nAB = nA + nB

      delta_T  = zonal_avg_loc(k,bin,1) - zonal_avg_glo(k,bin,1)
      delta_U  = zonal_avg_loc(k,bin,3) - zonal_avg_glo(k,bin,3)
      delta_V  = zonal_avg_loc(k,bin,4) - zonal_avg_glo(k,bin,4)
      delta_KE = zonal_avg_loc(k,bin,5) - zonal_avg_glo(k,bin,5)

      mA_T   = zonal_avg_glo(k,bin,2)
      mB_T   = zonal_avg_loc(k,bin,2)

      mA_UV  = zonal_avg_glo(k,bin,6)
      mB_UV  = zonal_avg_loc(k,bin,6)

      mA_zonal = zonal_avg_glo(k,bin,7)
      mB_zonal = zonal_avg_loc(k,bin,7)

      mA_merid = zonal_avg_glo(k,bin,8)
      mB_merid = zonal_avg_loc(k,bin,8)

      mA_VT  = zonal_avg_glo(k,bin,9)
      mB_VT  = zonal_avg_loc(k,bin,9)

      ! Combine means
      zonal_avg_glo(k,bin,1) = zonal_avg_glo(k,bin,1) + delta_T  * nB/nAB
      zonal_avg_glo(k,bin,3) = zonal_avg_glo(k,bin,3) + delta_U  * nB/nAB
      zonal_avg_glo(k,bin,4) = zonal_avg_glo(k,bin,4) + delta_V  * nB/nAB
      zonal_avg_glo(k,bin,5) = zonal_avg_glo(k,bin,5) + delta_KE * nB/nAB

      ! Combine sums of squares (for variances)
      zonal_avg_glo(k,bin,2) = mA_T     + mB_T     + delta_T**2      * nA*nB/nAB ! temperature variance
      zonal_avg_glo(k,bin,6) = mA_UV    + mB_UV    + delta_U*delta_V * nA*nB/nAB ! velocity covariance
      zonal_avg_glo(k,bin,7) = mA_zonal + mB_zonal + delta_U**2      * nA*nB/nAB ! zonal wind variance
      zonal_avg_glo(k,bin,8) = mA_merid + mB_merid + delta_V**2      * nA*nB/nAB ! meridional wind variance
      zonal_avg_glo(k,bin,9) = mA_VT    + mB_VT    + delta_V*delta_T * nA*nB/nAB ! V-T covariance (eddy heat flux)

      ! Update total number of data points
      Nstats_glo(k,bin) = nAB
    end subroutine combine_var
  end subroutine combine_stats
end module comm_mpi_mod
