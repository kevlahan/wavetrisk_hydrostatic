module lnorms_mod
  use domain_mod
  use comm_mpi_mod
  implicit none
contains
  subroutine cal_lnorm_sol (scaling, order)
    ! Calculates l norm of a float_field
    implicit none

    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: scaling
    character(*)                                              :: order

    integer :: k, l, v

    call update_array_bdry (scaling, NONE)

    lnorm = 0d0
    do k = 1, zmax
       do l = level_start, level_end
          select case (order)
          case ("1")
             call apply_onescale (l1_scalar_sol, l, k, 0, 1)
             call apply_onescale (l1_velo_sol,   l, k, 0, 0)
          case ("2")
             call apply_onescale (l2_scalar_sol, l, k, 0, 1)
             call apply_onescale (l2_velo_sol,   l, k, 0, 0)
          case ("inf")
             call apply_onescale (linf_scalar_sol, l, k, 0, 1)
             call apply_onescale (linf_velo_sol,   l, k, 0, 0)
          end select
       end do
       select case (order)
       case ("1", "2")
          do v = scalars(1), scalars(2)
             lnorm(v,k) = sum_real (lnorm(v,k))
          end do
          lnorm(S_VELO,k) = sum_real (lnorm(S_VELO,k))
       case ("inf")
          do v = scalars(1), scalars(2)
             lnorm(v,k) = sync_max_real (lnorm(v,k))
          end do
          lnorm(S_VELO,k) = sync_max_real (lnorm(S_VELO,k))
       end select
    end do
  end subroutine cal_lnorm_sol

  subroutine l1_scalar_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = lnorm(v,zlev) + abs (sol(v,zlev)%data(d)%elts(id+1))
       end do
    endif
  end subroutine l1_scalar_sol

  subroutine l1_velo_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) then
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + abs (sol(S_VELO,zlev)%data(d)%elts(id_e))
       end if
    end do
  end subroutine l1_velo_sol

  subroutine l2_scalar_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = lnorm(v,zlev) + sol(v,zlev)%data(d)%elts(id+1)**2
       end do
    end if
  end subroutine l2_scalar_sol

  subroutine l2_velo_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) then
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + sol(S_VELO,zlev)%data(d)%elts(id_e)**2
       end if
    end do
  end subroutine l2_velo_sol

  subroutine linf_scalar_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, v

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = max (lnorm(v,zlev), abs (sol(v,zlev)%data(d)%elts(id+1)))
       end do
    end if
  end subroutine linf_scalar_sol

  subroutine linf_velo_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_e, e

    d = dom%id+1
    id = idx(i, j, offs, dims)

    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) then
          lnorm(S_VELO,zlev) = max (lnorm(S_VELO,zlev), abs (sol(S_VELO,zlev)%data(d)%elts(id_e)))
       end if
    end do
  end subroutine linf_velo_sol
end module lnorms_mod
