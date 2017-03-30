module arch_mod
  use shared_mod
  use geom_mod
  implicit none
  integer rank
  integer n_process
  integer, dimension(N_GLO_DOMAIN) :: loc_id
  integer, dimension(N_GLO_DOMAIN) :: owner
  integer, allocatable :: glo_id(:,:)

contains

  subroutine init_arch_mod()
    integer p, d
    logical :: initialized = .False.

    if (initialized) return ! initialize only once

    call init_shared_mod()

    rank = 0
    n_process = 1

    allocate(n_domain(n_process))

    n_domain = N_GLO_DOMAIN

    initialized = .True.

  end subroutine init_arch_mod

  subroutine distribute_grid(dummy)
    integer p, d, dummy

    loc_id = 0
    owner = 0

    if (.not. allocated(glo_id)) then
       allocate(glo_id(1,n_domain(rank+1)))
    end if

    glo_id = 0
    p = 1
    do d = 1, n_domain(rank+1)
       owner(d) = p - 1 !serial execution so one owner only
       loc_id(d) = d - 1 !the local id is counting up
       glo_id(p,d) = d - 1 !single owner so global id simply counting up as well
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
