module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none
  integer, dimension(:), allocatable   :: recv_lengths, recv_offsets, req, send_lengths, send_offsets
  integer, dimension(:,:), allocatable :: stat_ray
contains
  subroutine init_comm_mpi
    ! Needed for compatibility with mpi code (not actually used in serial case)
    allocate(send_lengths(n_process), send_offsets(n_process))
    allocate(recv_lengths(n_process), recv_offsets(n_process))
    allocate(req(2*n_process))
    allocate(stat_ray(1,2*n_process))
    
    call init_comm
    call comm_communication_mpi
  end subroutine init_comm_mpi

  integer function write_active_per_level()
    ! write out distribution of active nodes over levels
    implicit none

    integer :: l, recommended_level_start

    recommended_level_start = level_start
    do l = level_start, level_end
       if (rank == 0) write(*,'(A,I2,I9)') 'lev', l, n_active_nodes(l), n_active_edges(l)
    end do

    if (rank == 0) write(*,'(A,I9)') 'total', sum(n_active(AT_NODE:AT_EDGE))

    write_active_per_level = recommended_level_start
  end function write_active_per_level

  subroutine cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    implicit none
    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance

    min_load = 1; max_load = 1; avg_load = 1.0_8; rel_imbalance = 1.0_8
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
    open (unit=funit, file=filename)
    if (eval_pole) call apply_to_pole (out_rout, l, zlev, funit, .False.)
    call apply_onescale__int (out_rout, l, zlev, 0, 0, funit)
    close (funit)
  end subroutine write_level_mpi

  subroutine write_load_conn (id, run_id)
    ! write out load distribution and connectivity for load balancing
    implicit none
    integer        :: id
    character(*)   :: run_id
    
    character(255) :: filename
    integer        :: fid

    fid = 599

    write (filename,'(A,A,I4.4)') trim (run_id), "_conn.", id
    open (unit=fid, file=trim(filename))
    call write_load_conn1 (fid)
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

  real(8) function cpt_dt ()
    implicit none

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: l

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
    n_active = (/ sum(n_active_nodes), sum(n_active_edges) /)
    cpt_dt = dt_loc
  end function cpt_dt

  real(8) function cpt_min_mass ()
    implicit none

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: l

    min_mass_loc = 1.0d16
    do l = level_start, level_end
       call apply_onescale (cal_min_mass, l, z_null, 0, 0)
    end do
    cpt_min_mass = min_mass_loc
  end function cpt_min_mass

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
  end subroutine start_timing

  subroutine stop_timing
  end subroutine stop_timing

  real(8) function get_timing()
    implicit none
    get_timing = 0.0_8
  end function get_timing

  subroutine sync_array (arr, N)
    implicit none
    integer            :: N
    real, dimension(N) :: arr
  end subroutine sync_array
  
  subroutine combine_stats
    implicit none
    Nstats_glo    = Nstats
    zonal_avg_glo = zonal_avg
  end subroutine combine_stats
end module comm_mpi_mod
