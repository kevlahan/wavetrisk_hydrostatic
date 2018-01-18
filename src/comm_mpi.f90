module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
  type(int_Array) recv_buf_i, send_buf_i
  type(float_Array) recv_buf, send_buf
  integer, allocatable :: send_lengths(:), send_offsets(:)
  integer, allocatable :: recv_lengths(:), recv_offsets(:)
  real(8) times(2)
  integer nreq
  integer, allocatable :: req(:)
  integer, allocatable :: stat_ray(:,:)

contains

  subroutine init_comm_mpi()
    allocate(send_lengths(n_process), send_offsets(n_process))
    allocate(recv_lengths(n_process), recv_offsets(n_process))
    allocate(req(2*n_process))
    allocate(stat_ray(MPI_STATUS_SIZE,2*n_process))
    call init_comm()
    call comm_communication_mpi()
  end subroutine init_comm_mpi

  integer function write_active_per_level()
    ! write out distribution of active nodes over levels
    integer l, n_full, fillin, n_lev_cur
    integer, dimension(2*(level_end-level_start+1)) :: n_active_all_loc, n_active_all_glo
    integer, dimension(level_start:level_end) :: n_active_per_lev
    integer recommended_level_start
    real(8) dt

    dt = cpt_dt_mpi() ! to set n_active_*

    n_lev_cur = level_end - level_start + 1

    n_active_all_loc = (/n_active_nodes(level_start:level_end), n_active_edges(level_start:level_end)/)

    call MPI_Allreduce(n_active_all_loc, n_active_all_glo, n_lev_cur*2, &
         MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror) !sum n_active_all_loc up across all processes and distribute result n_active_all_glo among all processes

    n_active_nodes(level_start:level_end) = n_active_all_glo(1:n_lev_cur)
    n_active_edges(level_start:level_end) = n_active_all_glo(n_lev_cur+1:n_lev_cur*2)
    n_active_per_lev = n_active_edges(level_start:level_end) + n_active_nodes(level_start:level_end)

    if (rank .eq. 0) write(*,'(6X,A,A,3(1X,A))') '   N_p   ', '   N_u   ','of all active', 'of full level', 'fill-in'

    recommended_level_start = level_start

    do l = level_start, level_end
       n_full = max_nodes_per_level(l) + max_nodes_per_level(l,EDGE)

       ! fill-in: additional nodes on level `l` if it'd become lowest level 
       ! minus the nodes on lower levels which would be removed
       fillin = n_full-n_active_per_lev(l)-sum(n_active_per_lev(level_start:l-1))

       if (rank .eq. 0) then
          write(*,'(A,I2,I9,I9,2(1X,F9.1,A),1X,I9,1X,F9.1,A)') &
               'lev', l, n_active_nodes(l), n_active_edges(l), &
               float(n_active_per_lev(l))/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%', &
               float(n_active_per_lev(l))/float(n_full)*100.0, '%', &
               fillin, float(fillin)/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%'
       end if

       if (fillin .le. 0) recommended_level_start = l
    end do

    if (rank .eq. 0) then
       write(*,'(A,I9,I9,2(1X,F9.1,A),9X,I9)') 'total', n_active(AT_NODE:AT_EDGE), 100.0, '%', &
            float(sum(n_active(AT_NODE:AT_EDGE)))/float(n_full)*100.0, '%', &
            n_full/sum(n_active(AT_NODE:AT_EDGE))
    end if

    write_active_per_level = recommended_level_start

  end function write_active_per_level

  subroutine write_load_conn (id)
    ! Write out load distribution and connectivity for load balancing
    integer :: id
    
    integer        :: r, fid
    character(5+4) :: filename

    fid = 599
    write(filename, '(A,I4.4)')  "conn.", id

    do r = 1, n_process
       if (r .ne. rank+1) then ! write only if our turn, otherwise only wait at Barrier
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if

       if (r .eq. 1) then ! first process opens without append to delete old file if existing
          open (unit=fid, file=filename, recl=333333)
       else
          open (unit=fid, file=filename, recl=333333, access='APPEND')
       end if

       call write_load_conn1 (fid)
       close(fid)
       call MPI_Barrier (MPI_Comm_World, ierror)
    end do
  end subroutine write_load_conn

  subroutine get_load_balance (mini, avg, maxi)
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

  subroutine print_load_balance
    ! Prints out load balance between processors
    integer :: load_min, load_max
    real(8) :: rel_imbalance, load_avg

    call get_load_balance (load_min, load_avg, load_max)
    rel_imbalance = dble(load_max)/load_avg

    if (rank .eq. 0) then
       write(6,'(A,1x,i9,1x,f10.1,1x,i9)') 'min load, average load, max load:', load_min, load_avg, load_max
       write(6,'(A,1x,f10.2)') 'relative imbalance (1=perfect balance)', rel_imbalance
    end if
  end subroutine print_load_balance

  subroutine write_level_mpi (out_rout, fid, l, zlev, eval_pole)
    external       :: out_rout
    integer        :: fid, l, zlev
    character(5+6) :: filename
    logical        :: eval_pole

    integer        :: r
    
    write (filename, '(A,I6)')  "fort.", fid
    
    do r = 1, n_process
       if (r .ne. rank+1) then ! Write only if our turn, otherwise only wait at Barrier0
          call MPI_Barrier (MPI_Comm_World, ierror)
          cycle 
       end if

       if (r .eq. 1) then ! first process opens without append to delete old file if existing
          open (unit=fid, file=filename)
       else
          open (unit=fid, file=filename, access='APPEND')
       end if

       if (eval_pole) call apply_to_pole (out_rout, l, zlev, fid, .true.)

       call apply_onescale__int (out_rout, l, zlev, 0, 0, fid)

       close (fid)

       call MPI_Barrier (MPI_Comm_World, ierror)
    end do
  end subroutine write_level_mpi

  subroutine init_comm_mpi_mod
    call init (send_buf_i, 0)
    call init (recv_buf_i, 0) 
    call init (send_buf, 0)
    call init (recv_buf, 0) 
  end subroutine init_comm_mpi_mod

  subroutine comm_communication_mpi
    call alltoall_dom (unpack_comm_struct, 4)
    call comm_communication
    call recreate_send_patch_lists
  end subroutine comm_communication_mpi

  subroutine recreate_send_patch_lists
    integer :: l, r, d, k, p, s, n, typ, d_neigh, i

    do l = level_start, level_end
       do d = 1, size(grid)
          do r = 1, n_process
             grid(d)%src_patch(r,l)%length = 0
          end do
          do k = 1, grid(d)%lev(l)%length
             p = grid(d)%lev(l)%elts(k)
             do s = 1, N_BDRY
                n = grid(d)%patch%elts(p+1)%neigh(s)
                if (n .ge. 0) cycle
                n = -n
                typ = grid(d)%bdry_patch%elts(n+1)%side
                if (typ .lt. 1) cycle
                d_neigh = grid(d)%neigh(typ)
                if (d_neigh .eq. POLE) then
                   do i = 1, 2
                      call handle_neigh(grid(d), grid(d)%neigh_over_pole(i))
                   end do
                else if (d_neigh .eq. NONE) then
                   cycle
                else
                   call handle_neigh(grid(d), d_neigh)
                end if
             end do
          end do
       end do
    end do
  contains
    subroutine handle_neigh(dom, d0)
      type(Domain) :: dom
      integer      :: d0
      
      integer :: r0

      r0 = owner(d0+1) + 1
      if (r0 .eq. rank+1) return

      if (dom%src_patch(r0,l)%length .gt. 0) then
         ! skip if just added
         if (dom%src_patch(r0,l)%elts(dom%src_patch(r0,l)%length) .eq. p) return
      end if

      call append(dom%src_patch(r0,l), p)
    end subroutine handle_neigh

  end subroutine recreate_send_patch_lists

  subroutine alltoall_dom (unpack_rout, N)
    external :: unpack_rout
    integer :: N
    
    integer :: src, dest, r_src, r_dest, d_src, d_dest, i, k, length
    integer, dimension(N) :: st

    send_buf_i%length = 0 ! reset
    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf_i%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             src = d_src
             dest = glo_id(r_dest,d_dest)+1
             call append(send_buf_i, grid(src)%send_conn(dest)%length)
             do i = 1, grid(src)%send_conn(dest)%length
                call append(send_buf_i, grid(src)%send_conn(dest)%elts(i))
             end do
             grid(src)%send_conn(dest)%length = 0
          end do
       end do
       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    call alltoall

    i = 1
    do r_src = 1, n_process
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             length = recv_buf_i%elts(i); i = i + 1
             do k = 1, length, N
                st = recv_buf_i%elts(i+0:i+(N-1))
                call unpack_rout(grid(d_dest), glo_id(r_src,d_src), st(1), st(2), st(3), st(4))
                i = i + N
             end do
          end do
       end do
    end do
  end subroutine alltoall_dom

  subroutine check_alltoall_lengths
    integer :: test_recv_len(n_process)

    call MPI_Alltoall(send_lengths, 1, MPI_INTEGER, &
         test_recv_len, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierror)

    write(3000+rank,*) test_recv_len-recv_lengths

    close(3000+rank)
    call MPI_Barrier(MPI_Comm_World, ierror)
  end subroutine check_alltoall_lengths

  subroutine alltoall
    integer :: i

    call MPI_Alltoall(send_lengths, 1, MPI_INTEGER, &
         recv_lengths, 1, MPI_INTEGER, &
         MPI_COMM_WORLD, ierror)

    recv_offsets(1) = 0

    do i = 2, n_process
       recv_offsets(i) = recv_offsets(i-1) + recv_lengths(i-1)
    end do

    recv_buf_i%length = sum(recv_lengths)

    if (size(recv_buf_i%elts) .lt. recv_buf_i%length) then
       deallocate(recv_buf_i%elts)
       allocate(recv_buf_i%elts(recv_buf_i%length))
    end if

    call MPI_Alltoallv(send_buf_i%elts, send_lengths, send_offsets, MPI_INTEGER, &
         recv_buf_i%elts, recv_lengths, recv_offsets, MPI_INTEGER, &
         MPI_COMM_WORLD, ierror)
  end subroutine alltoall

  subroutine comm_masks_mpi (l)
    !communication of mask information in a subdomain between different processes
    integer :: l
    
    integer :: r_dest, r_src, d_src, d_dest, dest, id, i, kk

    send_buf_i%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication

       send_offsets(r_dest) = send_buf_i%length

       do d_src = 1, n_domain(rank+1)

          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain

          do d_dest = 1, n_domain(r_dest)

             dest = glo_id(r_dest,d_dest)+1

             do i = 1, grid(d_src)%pack(AT_NODE,dest)%length
                id = grid(d_src)%pack(AT_NODE,dest)%elts(i)
                if (l .eq. NONE .or. l .eq. grid(d_src)%level%elts(id+1)) then
                   call append(send_buf_i, grid(d_src)%mask_n%elts(abs(id)+1))
                end if
             end do

             do i = 1, grid(d_src)%pack(AT_EDGE,dest)%length
                id = grid(d_src)%pack(AT_EDGE,dest)%elts(i)
                if (l .eq. NONE .or. l .eq. grid(d_src)%level%elts(id/EDGE+1)) then
                   call append(send_buf_i, grid(d_src)%mask_e%elts(abs(id)+1))
                end if
             end do

          end do
       end do
       send_lengths(r_dest) = send_buf_i%length - send_offsets(r_dest)
    end do

    ! determine recv buff lengths
    recv_buf_i%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf_i%length
       do d_src = 1, n_domain(r_src)
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = abs(grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i))
                if (l .eq. NONE .or. l .eq. grid(d_dest)%level%elts(id+1)) &
                     recv_buf_i%length = recv_buf_i%length + 1
             end do
             do i = 1, grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%length
                id = abs(grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%elts(i))
                if (l .eq. NONE .or. l .eq. grid(d_dest)%level%elts(id/EDGE+1)) &
                     recv_buf_i%length = recv_buf_i%length + 1
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf_i%length - recv_offsets(r_src)
    end do
    if (size(recv_buf_i%elts) .lt. recv_buf_i%length) then
       deallocate(recv_buf_i%elts)
       allocate(recv_buf_i%elts(recv_buf_i%length))
    end if

    !     call check_alltoall_lengths()
    call MPI_Alltoallv(send_buf_i%elts, send_lengths, send_offsets, MPI_INTEGER, &
         recv_buf_i%elts, recv_lengths, recv_offsets, MPI_INTEGER, &
         MPI_COMM_WORLD, ierror)

    ! communicate inside domain
    call comm_masks()

    kk = 0
    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                if (l .eq. NONE .or. l .eq. grid(d_dest)%level%elts(abs(id)+1)) then
                   kk = kk + 1
                   grid(d_dest)%mask_n%elts(abs(id)+1) = recv_buf_i%elts(kk)
                end if
             end do
             do i = 1, grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_EDGE,glo_id(r_src,d_src)+1)%elts(i)
                if (l .eq. NONE .or. l .eq. grid(d_dest)%level%elts(abs(id)/EDGE+1)) then
                   kk = kk + 1
                   grid(d_dest)%mask_e%elts(abs(id)+1) = recv_buf_i%elts(kk)
                end if
             end do
          end do
       end do
    end do
  end subroutine comm_masks_mpi

  subroutine update_bdry1 (field, l_start, l_end)
    type(Float_Field) :: field
    integer l_start, l_end

    call update_bdry__start1 (field, l_start, l_end)
    call update_bdry__finish1 (field, l_start, l_end)
  end subroutine update_bdry1

  subroutine update_array_bdry1 (field, l_start, l_end)
    type(Float_Field), dimension(:,:) :: field
    integer l_start, l_end

    call update_array_bdry__start1 (field, l_start, l_end)
    call update_array_bdry__finish1(field, l_start, l_end)
  end subroutine update_array_bdry1

  subroutine update_bdry (field, l)
    type(Float_Field) :: field
    integer           :: l

    call update_bdry__start (field, l)
    call update_bdry__finish (field, l)
  end subroutine update_bdry

  subroutine update_vector_bdry (field, l)
    ! Updates field array
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    call update_vector_bdry__start (field, l)
    call update_vector_bdry__finish(field, l)
  end subroutine update_vector_bdry
  
  subroutine update_array_bdry (field, l)
    ! Updates field array
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    call update_array_bdry__start  (field, l)
    call update_array_bdry__finish (field, l)
  end subroutine update_array_bdry

  subroutine update_bdry__start (field, l)
    type(Float_Field) :: field
    integer           ::  l

    if (l .eq. NONE) then 
       call update_bdry__start1(field, level_start-1, level_end)
    else
       call update_bdry__start1(field, l, l)
    endif
  end subroutine update_bdry__start
  
  subroutine update_vector_bdry__start (field, l)
    ! Finishes boundary update for field arrays
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l .eq. NONE) then 
       call update_vector_bdry__start1 (field, level_start-1, level_end)
    else
       call update_vector_bdry__start1 (field, l, l)
    endif
  end subroutine update_vector_bdry__start
  
  subroutine update_array_bdry__start (field, l)
    ! Finishes boundary update for field arrays
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l .eq. NONE) then 
       call update_array_bdry__start1 (field, level_start-1, level_end)
    else
       call update_array_bdry__start1 (field, l, l)
    endif
  end subroutine update_array_bdry__start

  subroutine update_bdry__start1 (field, l_start, l_end)
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: r_dest, r_src, d_src, d_dest, dest, id, i, k, r, multipl, lev

    if (field%bdry_uptodate) return

    if (field%pos .eq. AT_NODE) then
       multipl = 1
    else
       multipl = EDGE
    end if

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             do i = 1, grid(d_src)%pack(field%pos,dest)%length
                id = grid(d_src)%pack(field%pos,dest)%elts(i)
                lev = grid(d_src)%level%elts(id/multipl+1)
                if (l_start .le. lev .and. lev .le. l_end) then
                   call append (send_buf, field%data(d_src)%elts(id+1))
                end if
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
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%length
                id = abs(grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%elts(i))
                lev = grid(d_dest)%level%elts(id/multipl+1)
                if (l_start .le. lev .and. lev .le. l_end) then
                   recv_buf%length = recv_buf%length + 1
                end if
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) .lt. recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    ! Post all receives first
    nreq = 0
    do r = 1, n_process
       if (r .eq. rank+1 .or. recv_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_irecv (recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r .eq. rank+1 .or. send_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! communicate inside domain
    call cp_bdry_inside (field)
  end subroutine update_bdry__start1

  subroutine update_vector_bdry__start1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    type(Float_Field), dimension(:) :: field
    integer                         :: l_start, l_end
    
    integer :: i1, sz, r_dest, r_src, d_src, d_dest, dest, id, i, k, r, multipl, lev, pos
    logical :: ret

    ! Find shape of field
    sz = size(field)
    
    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i1 = 1, sz
       if (.not. field(i1)%bdry_uptodate) ret=.false.
    end do
    if (ret) return

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1
             ! Loop over each element of field array
             do i1 = 1, sz
                pos = field(i1)%pos
                if (pos .eq. AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_src)%pack(pos,dest)%length
                   id = grid(d_src)%pack(pos,dest)%elts(i)

                   lev = grid(d_src)%level%elts(id/multipl+1)
                   if (l_start .le. lev .and. lev .le. l_end) call append (send_buf, field(i1)%data(d_src)%elts(id+1))
                end do
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
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             ! Loop over each element of field array
             do i1 = 1, sz
                pos = field(i1)%pos
                if (pos .eq. AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))

                   lev = grid(d_dest)%level%elts(id/multipl+1)
                   if (l_start .le. lev .and. lev .le. l_end) recv_buf%length = recv_buf%length + 1
                end do
             end do

          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) .lt. recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    ! Post all receives first
    nreq = 0
    do r = 1, n_process
       if (r .eq. rank+1 .or. recv_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_irecv(recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r .eq. rank+1 .or. send_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! communicate inside domain
    call cp_bdry_inside_vector (field)
  end subroutine update_vector_bdry__start1

  subroutine update_array_bdry__start1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    type(Float_Field), dimension(:,:) :: field
    integer                           ::  l_start, l_end
    
    integer               :: i1, i2, r_dest, r_src, d_src, d_dest, dest, id, i, k, r, multipl, lev, pos
    integer, dimension(2) :: sz
    logical               :: ret

    ! Find shape of field
    sz = shape(field)

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i2 = 1, sz(2)
       do i1 = 1, sz(1)
          if (.not. field(i1,i2)%bdry_uptodate) ret=.false.
       end do
    end do
    if (ret) return

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain
          do d_dest = 1, n_domain(r_dest)
             dest = glo_id(r_dest,d_dest)+1

             ! Loop over each element of field array
             do i2 = 1, sz(2)
                do i1 = 1, sz(1)
                   pos = field(i1,i2)%pos
                   if (pos .eq. AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_src)%pack(pos,dest)%length
                      id = grid(d_src)%pack(pos,dest)%elts(i)

                      lev = grid(d_src)%level%elts(id/multipl+1)
                      if (l_start .le. lev .and. lev .le. l_end) call append (send_buf, field(i1,i2)%data(d_src)%elts(id+1))
                   end do
                end do
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
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)

             ! Loop over each element of field array
             do i2 = 1, sz(2)
                do i1 = 1, sz(1)
                   pos = field(i1,i2)%pos
                   if (field(i1,i2)%pos .eq. AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = abs(grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i))

                      lev = grid(d_dest)%level%elts(id/multipl+1)
                      if (l_start .le. lev .and. lev .le. l_end) recv_buf%length = recv_buf%length + 1
                   end do
                end do
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) .lt. recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    ! Post all receives first
    nreq = 0
    do r = 1, n_process
       if (r .eq. rank+1 .or. recv_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_irecv (recv_buf%elts(recv_offsets(r)+1), recv_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    do r = 1, n_process
       if (r .eq. rank+1 .or. send_lengths(r) .eq. 0) cycle
       nreq = nreq + 1
       call MPI_isend (send_buf%elts(send_offsets(r)+1), send_lengths(r), MPI_DOUBLE_PRECISION, &
            r-1, 1, MPI_COMM_WORLD, req(nreq), ierror)
    end do

    ! communicate inside domain
    call cp_bdry_inside_array (field)
  end subroutine update_array_bdry__start1

  subroutine update_bdry__finish (field, l)
    type(Float_Field) :: field
    integer           :: l

    if (l .eq. NONE) then 
       call update_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_bdry__finish1 (field, l, l)
    endif
  end subroutine update_bdry__finish
  
  subroutine update_vector_bdry__finish (field, l)
    ! Finishes boundary update for field arrays
    type(Float_Field), dimension(:) :: field
    integer                         :: l

    if (l .eq. NONE) then 
       call update_vector_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_vector_bdry__finish1 (field, l, l)
    endif
  end subroutine update_vector_bdry__finish
  
  subroutine update_array_bdry__finish (field, l)
    ! Finishes boundary update for field arrays
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    if (l .eq. NONE) then 
       call update_array_bdry__finish1 (field, level_start-1, level_end)
    else
       call update_array_bdry__finish1 (field, l, l)
    endif
  end subroutine update_array_bdry__finish

  subroutine update_bdry__finish1 (field, l_start, l_end)
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    integer :: r_dest, r_src, d_src, d_dest, dest, id, i, k, multipl, lev

    if (field%bdry_uptodate) return

    if (field%pos .eq. AT_NODE) then
       multipl = 1
    else
       multipl = EDGE
    end if

    k = 0
    call MPI_Waitall (nreq, req, stat_ray, ierror)
    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(field%pos,glo_id(r_src,d_src)+1)%elts(i)
                lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                if (l_start .le. lev .and. lev .le. l_end) then
                   k = k + 1
                   field%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                   if (id .lt. 0 .and. field%pos .eq. AT_EDGE) &
                        field%data(d_dest)%elts(abs(id)+1) = &
                        -field%data(d_dest)%elts(abs(id)+1)
                end if
             end do
          end do
       end do
    end do

    ! assumes routine is either called for one level, or all levels ever to be updated
    if (l_start .lt. l_end) field%bdry_uptodate = .True.
  end subroutine update_bdry__finish1

  subroutine update_vector_bdry__finish1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    type(Float_Field), dimension(:) :: field
    integer                         :: l_start, l_end
    
    integer :: i1, sz, r_dest, r_src, d_src, d_dest, dest, id, i, k, multipl, lev, pos
    logical :: ret

    ! Find shape of field
    sz = size(field)

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i1 = 1, sz
       if (.not. field(i1)%bdry_uptodate) ret=.false.
    end do
    if (ret) return

    k = 0
    call MPI_Waitall (nreq, req, stat_ray, ierror)
    
    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i1 = 1, sz
                pos = field(i1)%pos 
                if (pos .eq. AT_NODE) then
                   multipl = 1
                else
                   multipl = EDGE
                end if
                do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                   id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)

                   lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                   if (l_start .le. lev .and. lev .le. l_end) then
                      k = k + 1
                      field(i1)%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                      if (id .lt. 0 .and. pos .eq. AT_EDGE) &
                           field(i1)%data(d_dest)%elts(abs(id)+1) = -field(i1)%data(d_dest)%elts(abs(id)+1)
                   end if
                end do
             end do
          end do
       end do
    end do

    ! assumes routine is either called for one level, or all levels ever to be updated
    if (l_start .lt. l_end) field%bdry_uptodate = .True.
  end subroutine update_vector_bdry__finish1
  
  subroutine update_array_bdry__finish1 (field, l_start, l_end)
    ! Communicates boundary data in field, where fields is a Float_Field array
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l_start, l_end
    
    integer               ::  i1, i2, r_dest, r_src, d_src, d_dest, dest, id, i, k, multipl, lev, pos
    integer, dimension(2) :: sz
    logical               :: ret

    ! Find shape of field
    sz = shape(field)

    ! Check if boundaries of all field elements are up to date
    ret = .true.
    do i2 = 1, sz(2)
       do i1 = 1, sz(1)
          if (.not. field(i1,i2)%bdry_uptodate) ret=.false.
       end do
    end do
    if (ret) return

    k = 0
    call MPI_Waitall (nreq, req, stat_ray, ierror)

    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i2 = 1, sz(2)
                do i1 = 1, sz(1)
                   pos = field(i1,i2)%pos
                   if (pos .eq. AT_NODE) then
                      multipl = 1
                   else
                      multipl = EDGE
                   end if
                   do i = 1, grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%length
                      id = grid(d_dest)%unpk(pos,glo_id(r_src,d_src)+1)%elts(i)
                      lev = grid(d_dest)%level%elts(abs(id)/multipl+1)
                      if (l_start .le. lev .and. lev .le. l_end) then
                         k = k + 1
                         field(i1,i2)%data(d_dest)%elts(abs(id)+1) = recv_buf%elts(k)
                         if (id .lt. 0 .and. pos .eq. AT_EDGE) &
                              field(i1,i2)%data(d_dest)%elts(abs(id)+1) = -field(i1,i2)%data(d_dest)%elts(abs(id)+1)
                      end if
                   end do
                end do
             end do
          end do
       end do
    end do

    ! assumes routine is either called for one level, or all levels ever to be updated
    if (l_start .lt. l_end) field%bdry_uptodate = .True.
  end subroutine update_array_bdry__finish1
  
  subroutine comm_nodes9_mpi (get, set, l)
    external :: get, set
    integer :: l
    
    real(8), dimension(7) :: val
    integer               :: r_dest, r_src, d_src, d_dest, dest, id, i, k

    send_buf%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle 
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
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                recv_buf%length = recv_buf%length + 7
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) .lt. recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    call MPI_Alltoallv (send_buf%elts, send_lengths, send_offsets, MPI_DOUBLE_PRECISION, &
         recv_buf%elts, recv_lengths, recv_offsets, MPI_DOUBLE_PRECISION, &
         MPI_COMM_WORLD, ierror)

    call comm_nodes9 (get, set) ! communicate inside domain

    k = 0
    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                call set(grid(d_dest), id, recv_buf%elts(k+1:k+7))
                k = k + 7
             end do
          end do
       end do
    end do
  end subroutine comm_nodes9_mpi

  subroutine comm_nodes3_mpi (get, set, l)
    external    :: get, set
    type(Coord) :: get
    integer     :: l
    
    integer     :: r_dest, r_src, d_src, d_dest, dest, id, i, k
    type(Coord) :: c

    send_buf%length = 0 ! reset
    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) cycle ! TODO communicate inside domain
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

    ! determine recv buff lengths
    recv_buf%length = 0
    do r_src = 1, n_process 
       recv_offsets(r_src) = recv_buf%length
       do d_src = 1, n_domain(r_src)
          if (r_src .eq. rank+1) cycle 
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                recv_buf%length = recv_buf%length + 3
             end do
          end do
       end do
       recv_lengths(r_src) = recv_buf%length - recv_offsets(r_src)
    end do

    if (size(recv_buf%elts) .lt. recv_buf%length) then
       deallocate(recv_buf%elts)
       allocate(recv_buf%elts(recv_buf%length))
    end if

    call MPI_Alltoallv (send_buf%elts, send_lengths, send_offsets, MPI_DOUBLE_PRECISION, &
         recv_buf%elts, recv_lengths, recv_offsets, MPI_DOUBLE_PRECISION, &
         MPI_COMM_WORLD, ierror)

    call comm_nodes3 (get, set) ! communicate inside domain

    k = 0
    do r_src = 1, n_process 
       if (r_src .eq. rank+1) cycle ! inside domain
       do d_src = 1, n_domain(r_src)
          do d_dest = 1, n_domain(rank+1)
             do i = 1, grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%length
                id = grid(d_dest)%unpk(AT_NODE,glo_id(r_src,d_src)+1)%elts(i)
                c%x = recv_buf%elts(k+1) 
                c%y = recv_buf%elts(k+2) 
                c%z = recv_buf%elts(k+3) 
                call set(grid(d_dest), id, c)
                k = k + 3
             end do
          end do
       end do
    end do
  end subroutine comm_nodes3_mpi

  subroutine comm_patch_conn_mpi
    integer               :: r_src, r_dest, d_src, d_dest, i, b, c, p, s, d_glo, k, rot, d, ngh_pa, typ, l_par, rot_shift
    integer, dimension(4) :: st
    logical               :: is_pole
    
    send_buf_i%length = 0 ! reset

    do r_dest = 1, n_process ! destination for inter process communication
       send_offsets(r_dest) = send_buf_i%length
       do d_src = 1, n_domain(rank+1)
          if (r_dest .eq. rank+1) exit ! inside domain
          do d_dest = 1, n_domain(r_dest)
             do i = 1, grid(d_src)%send_pa_all%length, 4
                b = grid(d_src)%send_pa_all%elts(0+i)
                c = grid(d_src)%send_pa_all%elts(1+i)
                p = grid(d_src)%send_pa_all%elts(2+i)
                s = grid(d_src)%send_pa_all%elts(3+i)
                typ = grid(d_src)%bdry_patch%elts(b+1)%side
                d_glo = grid(d_src)%neigh(typ) ! 0 ...
                is_pole = d_glo .eq. POLE
                if (is_pole) then
                   d_glo = grid(d_src)%neigh_over_pole(c+1)
                   l_par = grid(d_src)%patch%elts(p+1)%level - 1
                   if (grid(d_src)%neigh_pa_over_pole%length .lt. l_par*2 + c + 1) then
                      ngh_pa = 0
                   else
                      ngh_pa = grid(d_src)%neigh_pa_over_pole%elts(l_par*2 + c + 1)
                   end if
                else
                   ngh_pa = grid(d_src)%bdry_patch%elts(b+1)%neigh
                end if
                ! also skips if dest .eq. 0
                if (ngh_pa .ne. 0 .and. d_glo .eq. glo_id(r_dest,d_dest)) then 
                   rot = grid(d_src)%neigh_rot(typ)
                   rot_shift = (rot_direction(grid(d_src), typ)*2 - 1)*rot
                   call append(send_buf_i, d_dest)
                   call append(send_buf_i, glo_id(rank+1,d_src))
                   call append(send_buf_i, p)
                   if (is_pole) then
                      call append(send_buf_i, c)
                   else
                      call append(send_buf_i, modulo(c + rot_shift, 4))
                   end if
                   call append(send_buf_i, ngh_pa)
                   if (is_pole) then
                      call append(send_buf_i, s)
                   else
                      call append(send_buf_i, modulo(s + rot_shift + 2, 4) + 4*(s/4))
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
       if (r_src .eq. rank+1) cycle ! inside domain
       do k = recv_offsets(r_src) + 1, recv_offsets(r_src) + recv_lengths(r_src), 6
          d = recv_buf_i%elts(k)
          d_src = recv_buf_i%elts(k+1)+1
          st = recv_buf_i%elts(k+2:k+5)
          call append(grid(d)%recv_pa(d_src), st(1))
          call append(grid(d)%recv_pa(d_src), st(2))
          call append(grid(d)%recv_pa(d_src), st(3))
          call append(grid(d)%recv_pa(d_src), st(4))
       end do
    end do
  end subroutine comm_patch_conn_mpi

  function cpt_dt_mpi()
    real(8) :: cpt_dt_mpi
    
    integer               :: l, ierror, n_level_glo
    integer, dimension(2) :: n_active_loc, n_active_glo
    real(8)               :: loc_min, glo_min

    if (adapt_dt) then
       dt_loc = 1d16
    else
       dt_loc = dt_init
    end if
    n_active_nodes = 0
    n_active_edges = 0

    min_mass = 1d16
    loc_min = min_mass
    call MPI_Allreduce (loc_min, glo_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    min_mass = glo_min

    do l = level_start, level_end
       call apply_onescale (min_dt, l, z_null, 0, 0)
    end do

    loc_min = dt_loc
    call MPI_Allreduce (loc_min, glo_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, MPI_COMM_WORLD, ierror)
    cpt_dt_mpi = glo_min
    
    n_active_loc = (/sum(n_active_nodes(level_start:level_end)), sum(n_active_edges(level_start:level_end))/)

    call MPI_Allreduce (n_active_loc, n_active_glo, 2, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    n_active = n_active_glo

    call MPI_Allreduce (level_end, n_level_glo, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
    level_end = n_level_glo
  end function cpt_dt_mpi

  function sync_max (val)
    integer :: sync_max
    integer :: val

    integer :: val_glo

    call MPI_Allreduce(val, val_glo, 1, MPI_INTEGER, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max = val_glo
  end function sync_max

  function sync_max_d (val)
    real(8) :: sync_max_d
    real(8) :: val

    real(8) :: val_glo

    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
    sync_max_d = val_glo
  end function sync_max_d

  function sum_real (val)
    real(8) :: sum_real
    real(8) :: val

    real(8) :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_real = val_glo
  end function sum_real

  function sum_int (val)
    integer :: sum_int
    integer :: val

    integer :: val_glo
    
    call MPI_Allreduce (val, val_glo, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierror)
    sum_int = val_glo
  end function sum_int

  subroutine start_timing
    times(1) = MPI_Wtime()  
  end subroutine start_timing

  subroutine stop_timing
    times(2) = MPI_Wtime()  
  end subroutine stop_timing

  function get_timing()
    real(8) :: get_timing
    get_timing = times(2) - times(1)
  end function get_timing

  subroutine sync (in, inout, len, type)
    real, dimension(len) :: in, inout
    integer              :: len, type

    where (in .ne. sync_val) inout = in
  end subroutine sync

  subroutine sync_array (arr, N)
    real, dimension(N) :: arr
    integer            :: N
    
    integer            :: myop
    real, dimension(N) :: garr

    call MPI_Op_create (sync, .True., myop, ierror)  
    call MPI_Reduce (arr, garr, N, MPI_REAL, myop, 0, MPI_COMM_WORLD, ierror)
    if (rank .eq. 0) arr(1:N) = garr(1:N)
  end subroutine sync_array

  subroutine stop_and_record_timings (id)
    ! use like:
    !call stop_and_record_timings(6500); call start_timing()
    !call stop_and_record_timings(6501); call start_timing()
    integer :: id

    real(8) :: time_loc, time_max, time_min, time_sum

    call stop_timing

    time_loc = get_timing()

    call MPI_Reduce (time_loc, time_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (time_loc, time_min, 1, MPI_DOUBLE_PRECISION, MPI_MIN, 0, MPI_COMM_WORLD, ierror)
    call MPI_Reduce (time_loc, time_sum, 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierror)

    if (rank .eq. 0) write(id,*) time_max, time_min, time_sum
  end subroutine stop_and_record_timings
end module comm_mpi_mod
