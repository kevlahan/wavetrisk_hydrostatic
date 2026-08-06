module arch_mod
  use mpi_f08

  use kind_mod,   only : dp
  use shared_mod, only : N_GLO_DOMAIN, init_shared_mod, n_domain, run_id

  implicit none

  private
  public :: MPI_IN, MPI_DP, MPI_SP
  public :: abort_run, barrier, comm, distribute_grid, glo_id, finalize, init_arch_mod, loc_id, n_process, owner, rank

  type(MPI_Comm)                :: comm
  type(MPI_Datatype), parameter :: MPI_IN = MPI_INTEGER
  type(MPI_Datatype), parameter :: MPI_DP = MPI_DOUBLE_PRECISION
  type(MPI_Datatype), parameter :: MPI_SP = MPI_REAL

  integer                              :: n_process, rank
  integer, dimension(N_GLO_DOMAIN)     :: loc_id, owner
  integer, dimension(:,:), allocatable :: glo_id

  
contains

  
  subroutine distribute_grid (cp_idx)
    ! Uses simple next-fit algorithm to allocate each domain to a processor (does not use adjacency information)
    ! Attempts to balance the total load using load data for each domain from checkpoint
    implicit none
    integer, intent(in) :: cp_idx

    integer                              :: d, r, n_domain_floor
    integer, dimension(:),   allocatable :: wgt_d
    integer, dimension(:,:), allocatable :: adj_line
    
    integer, parameter                   :: fid = 599
    
    real(dp)                             :: balanced_wgt, imbalance_goal
    real(dp), dimension(n_process)       :: wgt_cur_rank
    
    real(dp), parameter                  :: init_goal = 0.05_dp ! starting goal for maximum imbalance 
    real(dp), parameter                  :: incr_goal = 1.20_dp ! factor to increase goal by each iteration until domains fit

    character(255)                       :: filename

    allocate (wgt_d(N_GLO_DOMAIN), adj_line(N_GLO_DOMAIN,N_GLO_DOMAIN))
    
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
                if (wgt_cur_rank(r) > balanced_wgt * (1.0_dp + imbalance_goal)) then ! last domain unbalanced current rank -> put it on next rank
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
    call MPI_Bcast (owner, N_GLO_DOMAIN, MPI_IN, 0, comm)

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
    logical, save :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_shared_mod

    call MPI_Init 

    comm = MPI_COMM_WORLD
    
    call MPI_Comm_Size (comm, n_process)
    call MPI_Comm_Rank (comm, rank)

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
       call abort_run
    end if
  end subroutine init_arch_mod
  

  subroutine finalize
    call MPI_Finalize
  end subroutine finalize
  

  subroutine barrier
    call MPI_Barrier (comm)
  end subroutine barrier
  

  subroutine abort_run
    call MPI_Abort (comm, 1)
  end subroutine abort_run

  
end module arch_mod
