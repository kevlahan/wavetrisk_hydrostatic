module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
contains
  subroutine init_comm_mpi
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  function write_active_per_level ()
    ! write out distribution of active nodes over levels
    implicit none
    integer :: write_active_per_level
    integer :: l, recommended_level_start

    recommended_level_start = level_start
    do l = level_start, level_end
       if (rank == 0) write(*,'(A,I2,I9)') 'lev', l, n_active_nodes(l), n_active_edges(l)
    end do

    if (rank == 0) write(*,'(A,I9)') 'total', sum(n_active(AT_NODE:AT_EDGE))

    write_active_per_level = recommended_level_start
  end function write_active_per_level

  subroutine print_load_balance
    
  end subroutine print_load_balance

  subroutine write_level_mpi (out_rout, fid, l, zlev, eval_pole, test_case)
    implicit none
    external     :: out_rout
    integer      :: fid, l, zlev
    logical      :: eval_pole
    character(*) :: test_case

    integer, parameter :: funit = 300
    character(7)       :: var_file
    character(255)     :: filename
    
    write (var_file, '(I7)')  fid
    filename = trim(test_case)//'.'//var_file
    open (unit=funit, file=filename)
    if (eval_pole) call apply_to_pole (out_rout, l, zlev, funit, .False.)
    call apply_onescale__int (out_rout, l, zlev, 0, 0, funit)
    close (funit)
  end subroutine write_level_mpi

  subroutine write_load_conn (id)
    ! write out load distribution and connectivity for load balancing
    implicit none
    integer        :: id
    character(5+4) :: filename
    integer        :: fid

    fid = 599

    write(filename,'(A,I4.4)')  "conn.", id
    open(unit=fid, file=filename)
    call write_load_conn1(fid)
    close(fid)
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

  subroutine update_bdry1 (field, l_start, l_end)
    implicit none
    type(Float_Field) :: field
    integer           :: l_start, l_end
    
    call cp_bdry_inside(field)
  end subroutine update_bdry1

  subroutine update_array_bdry1 (field, l_start, l_end)
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l_start, l_end
    
    integer               :: i1, i2
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

  subroutine update_bdry (field, l)
    implicit none
    type(Float_Field) :: field
    integer           :: l
    call cp_bdry_inside (field)
  end subroutine update_bdry

  subroutine update_vector_bdry (field, l)
    implicit none
    type(Float_Field), dimension(:) :: field
    
    integer :: l, i1, i2, sz

    sz = size(field)

    do i1 = 1, sz
       call cp_bdry_inside (field(i1))
    end do
  end subroutine update_vector_bdry

  subroutine update_array_bdry (field, l)
    implicit none
    type(Float_Field), dimension(:,:) :: field
    integer                           :: l
    
    integer                :: i1, i2
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

  function cpt_dt_mpi()
    implicit none
    real(8) :: cpt_dt_mpi

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: d, i, ierror, j, l, p

    if (adapt_dt) then
       dt_loc = 1d16
    else
       dt_loc = dt_init
    end if
    
    n_active_nodes = 0
    n_active_edges = 0
    do l = level_start, level_end
       call apply_onescale (min_dt, l, z_null, 0, 0)
    end do
    change_mass = change_loc
    n_active = (/ sum(n_active_nodes), sum(n_active_edges) /)
    cpt_dt_mpi = dt_loc
  end function cpt_dt_mpi

  function sync_max (val)
    implicit none
    integer :: sync_max
    integer :: val
    
    sync_max = val
  end function sync_max

  function sync_max_d (val)
    implicit none
    real(8) :: sync_max_d
    real(8) :: val
    
    sync_max_d = val
  end function sync_max_d

  function sum_real (val)
    implicit none
    real(8) :: sum_real
    real(8) :: val
    
    sum_real = val
  end function sum_real

  function sum_int (val)
    integer :: sum_int
    integer :: val

    sum_int = val
  end function sum_int

  subroutine start_timing
  end subroutine start_timing

  subroutine stop_timing
  end subroutine stop_timing

  function get_timing ()
    implicit none
    real(8) :: get_timing
    get_timing = 0.0_8
  end function get_timing

  subroutine sync_array (arr, N)
    implicit none
    integer            :: N
    real, dimension(N) :: arr
  end subroutine sync_array
end module comm_mpi_mod
