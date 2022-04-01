module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
  integer, dimension(:), allocatable :: recv_lengths, recv_offsets, req, send_lengths, send_offsets
  real(8), dimension(2)              :: times
contains
  subroutine init_comm_mpi
    ! Needed for compatibility with mpi code (not actually used in serial case)
    allocate(send_lengths(n_process), send_offsets(n_process))
    allocate(recv_lengths(n_process), recv_offsets(n_process))
    allocate(req(2*n_process))
    recv_lengths = 0
    recv_offsets = 0
    send_lengths = 0
    send_offsets = 0
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  subroutine cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    implicit none
    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance

    min_load = 1; max_load = 1; avg_load = 1d0; rel_imbalance = 1d0
  end subroutine cal_load_balance

  subroutine write_level_mpi (out_rout, fid, l, zlev, eval_pole, run_id)
    implicit none
    external     :: out_rout
    integer      :: fid, l, zlev
    logical      :: eval_pole
    character(*) :: run_id

    integer, parameter :: funit = 300
    character(7)       :: var_file
    character(255)     :: filename
    
    write (var_file, '(I7)')  fid
    filename = trim(run_id)//'.'//var_file
    open (unit=funit, file=filename, status='REPLACE')
    if (eval_pole) call apply_to_pole (out_rout, l, zlev, funit, .false.)
    call apply_onescale__int (out_rout, l, zlev, 0, 0, funit)
    close (funit)
  end subroutine write_level_mpi

  subroutine write_load_conn (id, run_id)
    ! write out load distribution and connectivity for load balancing
    implicit none
    integer        :: id
    character(*)   :: run_id
    
    character(255) :: filename
    integer        :: d, fid, n_active_d

    fid = 599

    write (filename,'(A,A,I4.4)') trim (run_id), "_conn.", id
    open (unit=fid, file=trim(filename), status='REPLACE')
    do d = 1, size(grid)
       ! The following includes load for boundaries, but that seem just fair
       n_active_d = domain_load(grid(d))
       write(fid,'(I10, 99999(1X,I8))') n_active_d, ( &
            grid(d)%pack(AT_NODE,:)%length + grid(d)%pack(AT_EDGE,:)%length + &
            grid(d)%unpk(AT_NODE,:)%length + grid(d)%unpk(AT_EDGE,:)%length)/2
    end do
    close (fid)
  end subroutine write_load_conn

  subroutine init_comm_mpi_mod
  end subroutine init_comm_mpi_mod

  subroutine comm_communication_mpi
    call comm_communication
  end subroutine comm_communication_mpi

  subroutine comm_masks_mpi (l)
    implicit none
    integer :: l
    call comm_masks
  end subroutine comm_masks_mpi

  subroutine update_bdry1 (field, l_start, l_end, flag)
    implicit none
    type(Float_Field) :: field
    integer           :: flag, l_start, l_end
    
    call cp_bdry_inside(field)
  end subroutine update_bdry1

  subroutine update_vector_bdry1 (field, l_start, l_end, flag)
     implicit none
     integer                         :: flag, l_start, l_end
     type(Float_Field), dimension(:) :: field
    
    integer :: l, i1, sz

    sz = size(field)

    do i1 = 1, sz
       call cp_bdry_inside (field(i1))
    end do
  end subroutine update_vector_bdry1

  subroutine update_array_bdry1 (field, l_start, l_end, flag)
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l_start, l_end
    
    integer               :: flag, i1, i2
    integer, dimension(2) :: sz

    sz = shape(field)

    do i2 = 1, sz(2)
       do i1 = 1, sz(1)
          call cp_bdry_inside (field(i1,i2))
       end do
    end do
  end subroutine update_array_bdry1

  subroutine comm_nodes9_mpi (get, set, l)
    implicit none
    external :: get, set
    integer  :: l
    
    call comm_nodes9 (get, set) ! communicate inside domain
  end subroutine comm_nodes9_mpi

  subroutine comm_nodes3_mpi (get, set, l)
    implicit none
    external    :: get, set
    integer     :: l
    type(Coord) :: get
    
    call comm_nodes3 (get, set) ! communicate inside domain
  end subroutine comm_nodes3_mpi

  subroutine comm_patch_conn_mpi
    implicit none
    call comm_patch_conn
  end subroutine comm_patch_conn_mpi

  subroutine update_bdry (field, l, flag)
    implicit none
    type(Float_Field) :: field
    integer           :: flag, l
    call cp_bdry_inside (field)
  end subroutine update_bdry

  subroutine update_vector_bdry (field, l, flag)
    implicit none
    type(Float_Field), dimension(:) :: field
    
    integer :: flag, l, i1, sz

    sz = size(field)

    do i1 = 1, sz
       call cp_bdry_inside (field(i1))
    end do
  end subroutine update_vector_bdry

  subroutine update_array_bdry (field, l, flag)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l
    
    integer                :: flag, i1, i2
    integer, dimension (2) :: sz

    sz = shape(field)

    do i2 = 1, sz(2)
       do i1 = 1, sz(1)
          call cp_bdry_inside (field(i1,i2))
       end do
    end do
  end subroutine update_array_bdry

  subroutine update_bdry__start (field, l)
    implicit none
    type(Float_Field) :: field
    integer           :: l
    call cp_bdry_inside (field)
  end subroutine update_bdry__start
  
  subroutine update_vector_bdry__start (field, l)
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l
    
    integer :: i1, sz
    
    sz = size(field)

    do i1 = 1, sz
       call cp_bdry_inside (field(i1))
    end do
  end subroutine update_vector_bdry__start

  subroutine update_array_bdry__start (field, l)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l

    integer                :: i1, i2
    integer, dimension (2) :: sz

    sz = shape(field)

    do i2 = 1, sz(2)
       do i1 = 1, sz(1)
          call cp_bdry_inside(field(i1,i2))
       end do
    end do
  end subroutine update_array_bdry__start

  subroutine update_bdry__finish (field, l)
    implicit none
    type(Float_Field) :: field
    integer           :: l
  end subroutine update_bdry__finish

  subroutine update_vector_bdry__finish (field, l)
    implicit none
    type(Float_Field), dimension(:) :: field
    integer                         :: l
  end subroutine update_vector_bdry__finish

  subroutine update_array_bdry__finish (field, l)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l
  end subroutine update_array_bdry__finish

  integer function sync_max_int (val)
    implicit none
    integer :: val
    
    sync_max_int = val
  end function sync_max_int

  real(8) function sync_max_real (val)
    implicit none
    real(8) :: val
    
    sync_max_real = val
  end function sync_max_real

  real(8) function sum_real (val)
    implicit none
    real(8) :: val
    
    sum_real = val
  end function sum_real

  integer function sum_int (val)
    integer :: val

    sum_int = val
  end function sum_int

  subroutine start_timing
    implicit none
    call cpu_time(times(1))
  end subroutine start_timing

  subroutine stop_timing
    implicit none
    call cpu_time(times(2))
  end subroutine stop_timing

 real(8) function get_timing()
    implicit none
    get_timing = times(2) - times(1)
  end function get_timing  

  subroutine sync_array (arr, N)
    implicit none
    integer            :: N
    real(8), dimension(N) :: arr
  end subroutine sync_array
  
  subroutine combine_stats
    implicit none
    Nstats_glo    = Nstats
    zonal_avg_glo = zonal_avg
  end subroutine combine_stats
end module comm_mpi_mod
