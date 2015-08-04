module comm_mpi_mod
  use domain_mod
  use comm_mod
  implicit none

contains
  subroutine init_comm_mpi()
      call init_comm()
      call comm_communication_mpi()
  end subroutine

  integer function write_active_per_level()
  ! write out distribution of active nodes over levels
      integer l
      integer recommended_level_start
      recommended_level_start = level_start
      do l = level_start, level_end
          if (rank .eq. 0) write(*,'(A,I2,I9)') 'lev', l, n_active_mass(l), n_active_velo(l)
      end do
      if (rank .eq. 0) write(*,'(A,I9)') 'total', sum(n_active(S_MASS:S_VELO))
      write_active_per_level = recommended_level_start
  end function

  subroutine print_load_balance()
  end subroutine

  subroutine write_level_mpi(out_rout, fid, l, eval_pole)
      external out_rout
      integer fid, l
      character(5+6) filename
      logical eval_pole
      write(filename,   '(A,I6)')  "fort.", fid
      open(unit=fid, file=filename)
      if (eval_pole) call apply_to_pole(out_rout, l, fid, .False.)
      call apply_onescale__int(out_rout, l, 0, 0, fid)
      close(fid)
  end subroutine

  subroutine write_load_conn(id)
      ! write out load distribution and connectivity for load balancing
      integer id
      character(5+4) filename
      integer fid
      write(filename, '(A,I4.4)')  "conn.", id
      open(unit=fid, file=filename)
      call write_load_conn1(fid)
      close(fid)
  end subroutine

  subroutine init_comm_mpi_mod()
  end subroutine init_comm_mpi_mod

  subroutine comm_communication_mpi()
    call comm_communication()
  end subroutine

  subroutine comm_masks_mpi(l)
      integer l
      call comm_masks()
  end subroutine

  subroutine update_bdry1(field, l_start, l_end)
      type(Float_Field) :: field
      integer l_start, l_end
      call cp_bdry_inside(field)
  end subroutine

  subroutine comm_nodes9_mpi(get, set, l)
      external get, set
      integer l
      call comm_nodes9(get, set) ! communicate inside domain
  end subroutine

  subroutine comm_nodes3_mpi(get, set, l)
      external get, set
      integer l
      type(Coord) get
      call comm_nodes3(get, set) ! communicate inside domain
  end subroutine

  subroutine comm_patch_conn_mpi()
      call comm_patch_conn()
  end subroutine

  subroutine update_bdry(field, l)
      type(Float_Field) field
      integer l
      call cp_bdry_inside(field)
  end subroutine

  subroutine update_bdry__start(field, l)
      type(Float_Field) field
      integer l
      call cp_bdry_inside(field)
  end subroutine

  subroutine update_bdry__finish(field, l)
      type(Float_Field) field
      integer l
  end subroutine

  real(8) function cpt_dt_mpi()
    integer l, ierror
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer i, j
    integer d, p
    dt = 1.0e16_8
    n_active_mass = 0
    n_active_velo = 0
    do l = level_start, level_end
        call apply_onescale(min_dt, l, 0, 0)
    end do
!   TODO FIXME
!   do while (n_active_mass(level_end) .eq. 0 .and. &
!             n_active_velo(level_end) .eq. 0 )
!       level_end = level_end - 1
!   end do
    n_active = (/sum(n_active_mass), sum(n_active_velo)/)
  end function

  integer function sync_max(val)
    integer val
    sync_max = val
  end function

  real(8) function sync_max_d(val)
    real(8) val
    sync_max_d = val
  end function

  real(8) function sum_real(val)
     real(8) val
     sum_real = val
  end function

  subroutine start_timing()
  end subroutine

  subroutine stop_timing()
  end subroutine

  real(8) function get_timing()
    get_timing = 0.0_8
  end function

  subroutine sync_array(arr, N)
    real arr(N)
    integer N
  end subroutine
end module
