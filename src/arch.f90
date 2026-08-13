module arch_mod

  use mpi_f08

  use kind_mod,   only : dp
  use shared_mod, only : N_GLO_DOMAIN, n_domain, run_id

  implicit none

  private
  public :: MPI_IN, MPI_DP, MPI_SP
  public :: abort_run, barrier, comm, cp_load, distribute_grid, glo_id, finalize, init_arch_mod, loc_id, n_glo_block, &
       n_process, owner, rank

  type(MPI_Comm)                :: comm
  type(MPI_Datatype), parameter :: MPI_IN = MPI_INTEGER
  type(MPI_Datatype), parameter :: MPI_DP = MPI_DOUBLE_PRECISION
  type(MPI_Datatype), parameter :: MPI_SP = MPI_REAL

  integer                       :: n_process, rank
  integer                       :: n_glo_block
  integer,          allocatable :: loc_id(:), owner(:)
  integer,          allocatable :: glo_id(:,:)
  integer,          allocatable :: cp_load(:)
  

  type, public :: Parallel_Block
     integer :: id          = -1   ! global parallel-block ID
     integer :: root_domain = -1   ! original global geometric domain ID
     integer :: root_patch  = -1   ! root patch of subtree
     integer :: level       = -1   ! level of root patch
     integer :: owner       = -1   ! current MPI rank
     integer :: weight      = 0    ! provisional subtree weight
  end type Parallel_Block

  type(Parallel_Block), allocatable, public :: block_catalog(:)

  
contains


  subroutine init_arch_mod

    implicit none

    logical, save :: initialized = .false.

    if (initialized) return ! initialize only once

    call MPI_Init

    comm = MPI_COMM_WORLD

    call MPI_Comm_Size (comm, n_process)
    call MPI_Comm_Rank (comm, rank)

    allocate (n_domain(n_process))
    n_domain = 0

    ! Initially each geometric root domain is one parallel block.
    n_glo_block = N_GLO_DOMAIN

    allocate (loc_id(n_glo_block), owner(n_glo_block))

    loc_id = 0
    owner  = 0

    initialized = .true.

    if (n_process > n_glo_block) then
       if (rank == 0) then
          write (6,'(/,a)') &
               "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(2(a,i4),a)') &
               "!!          Number of cores ", n_process, &
               " > number of parallel blocks ", n_glo_block, &
               " ... aborting             !!"
          write (6,'(a,/)') &
               "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       end if
       call abort_run
    end if

  end subroutine init_arch_mod


  subroutine distribute_grid (cp_idx)
    ! Uses simple next-fit algorithm to allocate each parallel block to a
    ! processor. Attempts to balance the total load using load data stored
    ! in the checkpoint directory.
    !
    ! For cp_idx < 0, or when balancing is not appropriate, blocks are
    ! distributed equally by number.
    !
    ! Initially n_glo_block = N_GLO_DOMAIN, so each parallel block is
    ! identical to one geometric root domain.

    implicit none

    integer, intent(in) :: cp_idx

    integer                            :: d, r, n_domain_floor
    integer, dimension(:), allocatable :: wgt_d

    real(dp)                           :: balanced_wgt
    real(dp)                           :: imbalance_goal
    real(dp), dimension(n_process)     :: wgt_cur_rank

    real(dp), parameter :: init_goal = 0.05_dp
    real(dp), parameter :: incr_goal = 1.20_dp

    allocate (wgt_d(n_glo_block))

    if (rank == 0) then

       if (cp_idx >= 0 .and. &
            n_process > 1 .and. &
            n_process < n_glo_block) then

          !
          ! Load values have already been read from the checkpoint
          ! directory by read_checkpoint_directory().
          !
          if (.not. allocated(cp_load)) then
             error stop "distribute_grid: cp_load is not allocated"
          end if

          if (size(cp_load) /= n_glo_block) then
             error stop "distribute_grid: wrong size for cp_load"
          end if

          wgt_d = cp_load

          balanced_wgt = &
               real(sum(wgt_d), kind=dp) / real(n_process, kind=dp)

          !
          ! Use a variant of next-fit to maximize balance subject to:
          !
          !   - every rank receives at least one block
          !   - every block is assigned to a rank
          !
          d = 0
          imbalance_goal = init_goal

          do while (d < n_glo_block)

             d = 0
             wgt_cur_rank = 0.0_dp

             do r = 1, n_process

                do while ( &
                     wgt_cur_rank(r) < balanced_wgt .and. &
                     n_glo_block - d > n_process - r)

                   owner(d+1) = r - 1

                   wgt_cur_rank(r) = &
                        wgt_cur_rank(r) + real(wgt_d(d+1), kind=dp)

                   d = d + 1

                end do

                !
                ! If the final block makes this rank too heavy,
                ! move that block to the next rank.
                !
                if (d > 0) then
                   if (wgt_cur_rank(r) > &
                        balanced_wgt * (1.0_dp + imbalance_goal)) then

                      wgt_cur_rank(r) = &
                           wgt_cur_rank(r) - real(wgt_d(d), kind=dp)

                      d = d - 1

                   end if
                end if

             end do

             !
             ! If all blocks could not be fitted, relax the
             ! imbalance tolerance and try again.
             !
             if (d < n_glo_block) then
                imbalance_goal = imbalance_goal * incr_goal
             end if

          end do

          write (6,'(a,es8.2,/)') &
               'New maximum load imbalance = ', &
               maxval(wgt_cur_rank) / balanced_wgt

       else

          !
          ! Equal distribution by number of blocks.
          !
          n_domain_floor = n_glo_block / n_process

          d = 0

          do r = 1, n_process

             if (n_domain_floor > 0) then
                owner(d+1:d+n_domain_floor) = r - 1
                d = d + n_domain_floor
             end if

             if (r <= n_glo_block - n_process*n_domain_floor) then
                owner(d+1) = r - 1
                d = d + 1
             end if

          end do

       end if

    end if

    !
    ! Every rank now receives the new ownership map.
    !
    call MPI_Bcast ( &
         owner, n_glo_block, MPI_IN, 0, comm)

    !
    ! Construct local block numbering for the new distribution.
    !
    n_domain = 0

    do d = 1, n_glo_block

       r = owner(d)

       loc_id(d) = n_domain(r+1)

       n_domain(r+1) = n_domain(r+1) + 1

    end do

    !
    ! Construct global-block ID list for each rank.
    !
    if (allocated(glo_id)) deallocate (glo_id)

    allocate (glo_id(n_process, maxval(n_domain)))

    glo_id = 0

    do d = 1, n_glo_block

       glo_id( &
            owner(d)+1, &
            loc_id(d)+1) = d - 1

    end do

    deallocate (wgt_d)

  end subroutine distribute_grid


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
