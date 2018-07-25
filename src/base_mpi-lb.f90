module arch_mod
  use shared_mod
  use geom_mod
  implicit none
  include 'mpif.h'
  integer                              :: ierror, n_process, rank
  integer, dimension(N_GLO_DOMAIN)     :: loc_id, owner
  integer, dimension(:,:), allocatable :: glo_id
contains
  subroutine distribute_grid (k)
    implicit none
    integer :: k

    integer                          :: i, d, r, d_ngb, n_domain_floor, total_wgt
    integer, dimension(N_GLO_DOMAIN) :: adj_line, vwgt
    integer, parameter               :: fid = 599
    character(5+4)                   :: filename
    real(8)                          :: wgt_per_rank, wgt_cur_rank, accepted_imbalance
    
    write(filename, '(A,I4.4)')  "conn.", k

    if (k >= 0 .and. n_process > 1) then

       open (unit=fid, file=filename)
       do d = 1, N_GLO_DOMAIN
          read (fid,*) vwgt(d), adj_line
       end do
       close (fid)

       total_wgt = sum (vwgt)
       wgt_per_rank = dble(total_wgt)/dble(n_process)
       d = 0

       ! Goals:
       !  - every rank has at least one domain
       !  - every domain is assigned to a rank

       accepted_imbalance = 0.1_8

       do while (d < N_GLO_DOMAIN) ! increase accepted_imbalance until all domains fit
          d = 0
          do r = 1, n_process
             wgt_cur_rank = 0
             do while (wgt_cur_rank < wgt_per_rank .and. n_process - r < N_GLO_DOMAIN - d)
                owner(d+1) = r-1
                wgt_cur_rank = wgt_cur_rank + vwgt(d+1)
                d = d + 1
             end do
             ! If load too much, keep last item for next rank
             if (wgt_cur_rank > dble(wgt_per_rank)*(1.0_8 + accepted_imbalance)) d = d - 1
          end do
          ! Did not find enough room for all domains > accepted_imbalance was too tight
          accepted_imbalance = accepted_imbalance*2.0_8
       end do
    else
       n_domain_floor = N_GLO_DOMAIN/n_process
       d = 0
       do r = 1, n_process
          owner(d+1:d+n_domain_floor) = r-1
          d = d + n_domain_floor
          if (r <= N_GLO_DOMAIN - n_process*n_domain_floor) then
             owner(d+1) = r-1
             d = d + 1
          end if
       end do
    end if

    n_domain = 0
    do d = 1, N_GLO_DOMAIN
       r = owner(d)
       loc_id(d) = n_domain(r+1)
       n_domain(r+1) = n_domain(r+1) + 1
    end do

    if (allocated(glo_id)) deallocate (glo_id)
    allocate (glo_id(n_process,maxval(n_domain)))

    do d = 1, N_GLO_DOMAIN
       glo_id(owner(d)+1,loc_id(d)+1) = d - 1
    end do
  end subroutine distribute_grid

  subroutine init_arch_mod
    implicit none
    logical :: initialized = .False.

    if (initialized) return ! initialize only once

    call init_shared_mod

    call MPI_Init (ierror)
    call MPI_Comm_Size (MPI_COMM_WORLD, n_process, ierror)
    call MPI_Comm_Rank (MPI_COMM_WORLD, rank, ierror)

    allocate (n_domain(n_process))

    initialized = .True.
  end subroutine init_arch_mod

  subroutine finalize
    call MPI_Finalize (ierror)
  end subroutine finalize

  subroutine barrier
    call MPI_Barrier (MPI_Comm_World, ierror)
  end subroutine barrier

  subroutine abort
    call MPI_Abort (MPI_Comm_World, 1, ierror)
  end subroutine abort
end module arch_mod
