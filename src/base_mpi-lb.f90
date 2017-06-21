module arch_mod
  use shared_mod
  use geom_mod
  implicit none
  include 'mpif.h'
  integer rank
  integer n_process
  integer ierror
  integer, dimension(N_GLO_DOMAIN) :: loc_id
  integer, dimension(N_GLO_DOMAIN) :: owner
  integer, allocatable :: glo_id(:,:)

contains

  subroutine distribute_grid(k)
    integer k

    integer vwgt(N_GLO_DOMAIN)
    integer adj_line(N_GLO_DOMAIN)
    integer i, d, r, d_ngb
    integer n_domain_floor
    character(5+4) filename
    integer :: fid = 599
    integer total_wgt
    real(8) wgt_per_rank, wgt_cur_rank, accepted_inbalance

    write(filename, '(A,I4.4)')  "conn.", k

    if (k .ge. 0 .and. n_process .gt. 1) then

       open(unit=fid, file=filename)
       do d = 1, N_GLO_DOMAIN
          read(fid,*) vwgt(d), adj_line
       end do
       close(fid)

       total_wgt = sum(vwgt)
       wgt_per_rank = dble(total_wgt)/dble(n_process)
       d = 0

       ! Goals:
       !  - every rank has at least one domain
       !  - every domain is assigned to a rank

       accepted_inbalance = 0.1_8

       do while (d .lt. N_GLO_DOMAIN) ! increase accepted_inbalance until all domains fit
          d = 0
          do r = 1, n_process
             wgt_cur_rank = 0
             do while (wgt_cur_rank .lt. wgt_per_rank .and. n_process - r .lt. N_GLO_DOMAIN - d)
                owner(d+1) = r-1
                wgt_cur_rank = wgt_cur_rank + vwgt(d+1)
                d = d + 1
             end do
             ! if load too much, keep last item for next rank
             if (wgt_cur_rank .gt. dble(wgt_per_rank)*(1.0_8 + accepted_inbalance)) then
                d = d - 1
             end if
          end do
          ! did not find enough room for all domains > accepted_inbalance was too tight
          accepted_inbalance = accepted_inbalance*2.0_8
       end do
       if (rank .eq. 0) write(*,'(A,es12.4)') 'Accepted imbalance:', accepted_inbalance / 2.0_8
    else
       n_domain_floor = N_GLO_DOMAIN/n_process
       d = 0
       do r = 1, n_process
          owner(d+1:d+n_domain_floor) = r-1
          d = d + n_domain_floor
          if (r .le. N_GLO_DOMAIN - n_process*n_domain_floor) then
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

    if (allocated(glo_id)) deallocate(glo_id)
    allocate(glo_id(n_process,maxval(n_domain)))

    do d = 1, N_GLO_DOMAIN
       glo_id(owner(d)+1,loc_id(d)+1) = d - 1
    end do

  end subroutine distribute_grid

  subroutine init_arch_mod()
    logical :: initialized = .False.

    if (initialized) return ! initialize only once

    call init_shared_mod()

    call MPI_Init(ierror)
    call MPI_Comm_Size(MPI_COMM_WORLD, n_process, ierror)
    call MPI_Comm_Rank(MPI_COMM_WORLD, rank, ierror)

    allocate(n_domain(n_process))

    initialized = .True.

  end subroutine init_arch_mod

  subroutine finalize()
    call MPI_Finalize(ierror)
  end subroutine finalize

  subroutine barrier()
    call MPI_Barrier(MPI_Comm_World, ierror)
  end subroutine barrier

  subroutine abort()
    call MPI_Abort(MPI_Comm_World, 1, ierror)
  end subroutine abort

end module arch_mod
