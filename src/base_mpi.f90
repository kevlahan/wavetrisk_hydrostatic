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
    ! Allocates each domain to a processor
    implicit none
    integer      :: cp_idx
    character(*) :: run_id
    
    integer :: i, d, r, d_ngb, n_domain_floor

    if (rank == 0) write (6,'(/,a,/)') 'Distributing each domain to a processor'
    
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

    n_domain = 0
    do d = 1, N_GLO_DOMAIN
       r = owner(d)
       loc_id(d) = n_domain(r+1)
       n_domain(r+1) = n_domain(r+1) + 1
    end do

    if (.not. allocated (glo_id)) allocate (glo_id(n_process, maxval (n_domain)))

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
    call MPI_Barrier(MPI_Comm_World, ierror)
  end subroutine barrier

  subroutine abort
    call MPI_Abort(MPI_Comm_World, 1, ierror)
  end subroutine abort
end module arch_mod
