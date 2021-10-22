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
    ! Uses METIS library to allocate each domain to a processor
    ! Attempts to balance the total load using load and adjacency data for each domain from checkpoint
    implicit none
    integer        :: cp_idx
    character(255) :: run_id
    
    integer                                 :: d, d_adj, edgecut, i, n_domain_floor, r
    integer, parameter                      :: fid = 599

    integer, dimension(N_GLO_DOMAIN)        :: adj_line, wgt_d
    integer, dimension(N_GLO_DOMAIN+1)      :: xadj
    integer, dimension(N_GLO_DOMAIN*N_BDRY) :: adjncy, wgt_adj
    
    character(255)                          :: filename

    ! METIS options
    integer, parameter                     :: METIS_NOPTIONS = 40
    integer, dimension(0:METIS_NOPTIONS-1) :: options
    integer, parameter :: METIS_OPTION_OBJTYPE   = 1  !! Objective (0 = edge cut, 1 = total communication, 2 = node)
    integer, parameter :: METIS_OPTION_CTYPE     = 2  !! Matching scheme for coarsening (0 = random, 1 = sorted heavy edge).
    integer, parameter :: METIS_OPTION_IPTYPE    = 3  !! Algorithm for  initial partitioning (0 = greedy, 1 = random, 2 = from edge cut, 3 = greedy node-based).
    integer, parameter :: METIS_OPTION_RTYPE     = 4  !! Algorithm for refinement (0 = FM-based, 1 = greedy, 2 = 2-sided node FM, 3 = 1-sided node FM).
    integer, parameter :: METIS_OPTION_DBGLVL    = 5  !! Amount of progress/debugging information.
    integer, parameter :: METIS_OPTION_NITER     = 6  !! Number of iterations for refinement algorithm (default = 10).
    integer, parameter :: METIS_OPTION_NCUTS     = 7  !! Number of different partitionings computed (default = 1).
    integer, parameter :: METIS_OPTION_SEED      = 8  !! Seed for the random number generator.
    integer, parameter :: METIS_OPTION_NO2HOP    = 9  !! Use 2–hop matchings (0 = T, 1 = F).
    integer, parameter :: METIS_OPTION_MINCONN   = 10 !! Minimize maximum connnectivity  (0 = F, 1 = T).
    integer, parameter :: METIS_OPTION_CONTIG    = 11 !! Contiguous partitions (0 = T, 1 = F).
    integer, parameter :: METIS_OPTION_UFACTOR   = 16 !! Maximum allowed load imbalance (1 + x)/1000 (default is 1 for recursive, 30 for k-way).
    
    logical, parameter :: recursive = .false.         !! T = recursive partitioning, F = k-way partitioning.

    call METIS_SetDefaultOptions (options)

    if (recursive) then
       options(METIS_OPTION_CTYPE) = 1
       options(METIS_OPTION_NCUTS) = 4
    end if
    
    write (filename, '(a,a,i4.4)') trim (run_id), "_conn.", cp_idx
    
    if (cp_idx >= 0 .and. n_process > 1) then
       open (unit=fid, file=filename)
       i = 0
       xadj(1) = 0
       do d = 1, N_GLO_DOMAIN
          read(fid,*) wgt_d(d), adj_line
          do d_adj = 1, N_GLO_DOMAIN ! adjacency information for domain d for METIS
             if (adj_line(d_adj) > 0) then ! only include the N_BDRY edges for each domain
                i = i + 1              
                if (i > size(adjncy)) then
                   write(0,*) "ERROR: adjncy to short"
                   call MPI_Finalize (ierror)
                   stop
                end if
                adjncy(i)  = d_adj - 1
                wgt_adj(i) = adj_line(d_adj)
             end if
          end do
          xadj(d+1) = i
       end do
       close(fid)
       if (recursive) then
          call METIS_PartGraphRecursive (N_GLO_DOMAIN, 1, xadj, adjncy, wgt_d, , wgt_adj, n_process, , , options, edgecut, owner)
       else
          call METIS_PartGraphKway      (N_GLO_DOMAIN, 1, xadj, adjncy, wgt_d, , wgt_adj, n_process, , , options, edgecut, owner)
       end if
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

    n_domain = 0
    do d = 1, N_GLO_DOMAIN
       r = owner(d)
       loc_id(d) = n_domain(r+1)
       n_domain(r+1) = n_domain(r+1) + 1
    end do

    if (allocated (glo_id)) deallocate (glo_id)
    allocate (glo_id(n_process, maxval (n_domain)))

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
