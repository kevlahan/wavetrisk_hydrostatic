module arch_mod
  use shared_mod
  use geom_mod
  implicit none
  include 'mpif.h'
  integer                              :: ierror, n_process, rank
  integer, dimension(N_GLO_DOMAIN)     :: loc_id, owner
  integer, dimension(:,:), allocatable :: glo_id
contains
  subroutine distribute_grid (cp_idx, run_id)
    ! Uses simple next-fit algorithm to allocate each domain to a processor (does not use adjacency information)
    ! Attempts to balance the total load using load data for each domain from checkpoint
    implicit none
    integer        :: cp_idx
    character(255) :: run_id

    integer                                       :: d, r, n_domain_floor
    integer, dimension(N_GLO_DOMAIN)              :: wgt_d
    integer, dimension(N_GLO_DOMAIN,N_GLO_DOMAIN) :: adj_line
    
    integer, parameter                            :: fid = 599
    
    real(8)                                       :: balanced_wgt, imbalance_goal
    real(8), dimension(n_process)                 :: wgt_cur_rank
    
    real(8), parameter                            :: init_goal = 0.05d0 ! starting goal for maximum imbalance 
    real(8), parameter                            :: incr_goal = 1.20d0 ! factor to increase goal by each iteration until domains fit

    character(255)                                :: filename
    
    if (rank == 0) then
       if (cp_idx >= 0 .and. n_process > 1 .and. n_process < N_GLO_DOMAIN) then
          write (filename, '(a,a,i4.4)') trim (run_id), "_conn.", cp_idx
          open (unit=fid, file=trim(filename), status='OLD')
          do d = 1, N_GLO_DOMAIN
             read (fid,*) wgt_d(d), adj_line(d,:)
          end do
          close (fid)
          balanced_wgt = dble(sum(wgt_d)) / dble(n_process) ! average load per rank (perfect balance)

          ! Goals: use a variant of next-fit algorithm to maximize balance with the constraints that
          !  - every rank has at least one domain
          !  - every domain is assigned to a rank
          d = 0
          imbalance_goal = init_goal ! initial imbalance goal is 1 + imbalance_goal
          do while (d < N_GLO_DOMAIN)
             d = 0
             wgt_cur_rank = 0
             do r = 1, n_process
                do while (wgt_cur_rank(r) < balanced_wgt .and. N_GLO_DOMAIN - d > n_process - r)
                   owner(d+1) = r - 1
                   wgt_cur_rank(r) = wgt_cur_rank(r) + wgt_d(d+1)
                   d = d + 1
                end do
                if (wgt_cur_rank(r) > balanced_wgt * (1d0 + imbalance_goal)) then ! last domain unbalanced current rank -> put it on next rank
                   wgt_cur_rank(r) = wgt_cur_rank(r) - wgt_d(d)
                   d = d - 1 
                end if
             end do
             ! Not enough room for all domains -> increase imbalance_goal and try again
             imbalance_goal = imbalance_goal * incr_goal
          end do

          if (rank == 0) write (6,'(a,es8.2,/)') 'New maximum load imbalance = ', maxval (wgt_cur_rank) / balanced_wgt
       else ! distribute domains equally 
          n_domain_floor = N_GLO_DOMAIN / n_process
          d = 0
          do r = 1, n_process
             owner(d+1:d+n_domain_floor) = r - 1
             d = d + n_domain_floor
             if (r <= N_GLO_DOMAIN - n_process * n_domain_floor) then
                owner(d+1) = r - 1
                d = d + 1
             end if
          end do
       end if
    end if
    call MPI_Bcast (owner, N_GLO_DOMAIN, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)

    n_domain = 0
    do d = 1, N_GLO_DOMAIN
       r = owner(d)
       loc_id(d) = n_domain(r+1)
       n_domain(r+1) = n_domain(r+1) + 1
    end do

    if (allocated (glo_id)) deallocate (glo_id)
    allocate (glo_id(n_process, maxval (n_domain)))
    glo_id = 0

    do d = 1, N_GLO_DOMAIN
       glo_id(owner(d)+1,loc_id(d)+1) = d - 1
    end do
  end subroutine distribute_grid

  subroutine init_arch_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_shared_mod

    call MPI_Init (ierror)
    call MPI_Comm_Size (MPI_COMM_WORLD, n_process, ierror)
    call MPI_Comm_Rank (MPI_COMM_WORLD, rank,      ierror)

    allocate (n_domain(n_process))
    n_domain = 0
    
    initialized = .true.

    if (n_process > N_GLO_DOMAIN) then
       if (rank == 0) then
          write (6,'(/,a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(2(a,i4),a)') "!!          Number of cores ", n_process, " > number of domains ", N_GLO_DOMAIN, &
               " ... aborting             !!"
          write (6,'(a,/)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       end if
       call abort
    end if
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
