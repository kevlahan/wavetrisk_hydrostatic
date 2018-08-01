module arch_mod
  use shared_mod
  use geom_mod
  implicit none
  integer                              :: n_process, rank
  integer, dimension(N_GLO_DOMAIN)     :: loc_id
  integer, dimension(N_GLO_DOMAIN)     :: owner
  integer, dimension(:,:), allocatable :: glo_id
contains
  subroutine init_arch_mod
    implicit none
    
    integer :: p, d
    logical :: initialized = .false.

    if (initialized) return ! initialize only once

    call init_shared_mod

    rank = 0
    n_process = 1
    allocate (n_domain(n_process))
    n_domain = N_GLO_DOMAIN
    initialized = .true.
  end subroutine init_arch_mod

  subroutine distribute_grid (cp_idx)
    implicit none
    integer :: cp_idx

    integer :: p, d

    loc_id = 0
    owner = 0

    if (.not. allocated (glo_id)) allocate (glo_id(1,n_domain(rank+1)))

    glo_id = 0
    p = 1
    do d = 1, n_domain(rank+1)
       owner(d) = p - 1 ! serial execution so one owner only
       loc_id(d) = d - 1 ! the local id is counting up
       glo_id(p,d) = d - 1 ! single owner so global id simply counting up as well
    end do
  end subroutine distribute_grid

  subroutine finalize()
  end subroutine finalize

  subroutine barrier()
  end subroutine barrier

  subroutine abort()
    stop
  end subroutine abort
end module arch_mod
