module lnorms_mod
  use domain_mod
  use comm_mpi_mod
  implicit none
  type(Float_Field), dimension(:,:), allocatable, target :: scaling
contains
  subroutine cal_lnorm (q, order)
    ! Calculates l norm of a float_field
    implicit none

    type(Float_Field), dimension(1:N_VARIABLE,zmin:zmax), target :: q
    character(*)                                                 :: order

    integer :: k, l, v

    allocate (scaling(1:N_VARIABLE,zmin:zmax))
    scaling = q
    call update_array_bdry (scaling, NONE)
    
    lnorm = 0d0
    do k = zmin, zmax
       do l = level_start, level_end
          select case (order)
          case ("1")
             call apply_onescale (l1_scalar, l, k, 0, 1)
             call apply_onescale (l1_velo,   l, k, 0, 0)
          case ("2")
             call apply_onescale (l2_scalar, l, k, 0, 1)
             call apply_onescale (l2_velo,   l, k, 0, 0)
          case ("inf")
             call apply_onescale (linf_scalar, l, k, 0, 1)
             call apply_onescale (linf_velo,   l, k, 0, 0)
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
    
    deallocate (scaling)
  end subroutine cal_lnorm

  subroutine l1_scalar (dom, i, j, zlev, offs, dims)
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
          lnorm(v,zlev) = lnorm(v,zlev) + abs (scaling(v,zlev)%data(d)%elts(id+1))
       end do
    endif
  end subroutine l1_scalar

  subroutine l1_velo (dom, i, j, zlev, offs, dims)
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
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + abs (scaling(S_VELO,zlev)%data(d)%elts(id_e))
       end if
    end do
  end subroutine l1_velo

  subroutine l2_scalar (dom, i, j, zlev, offs, dims)
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
          lnorm(v,zlev) = lnorm(v,zlev) + scaling(v,zlev)%data(d)%elts(id+1)**2
       end do
    end if
  end subroutine l2_scalar

  subroutine l2_velo (dom, i, j, zlev, offs, dims)
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
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + scaling(S_VELO,zlev)%data(d)%elts(id_e)**2
       end if
    end do
  end subroutine l2_velo

  subroutine linf_scalar (dom, i, j, zlev, offs, dims)
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
          lnorm(v,zlev) = max (lnorm(v,zlev), abs (scaling(v,zlev)%data(d)%elts(id+1)))
       end do
    end if
  end subroutine linf_scalar

  subroutine linf_velo (dom, i, j, zlev, offs, dims)
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
          lnorm(S_VELO,zlev) = max (lnorm(S_VELO,zlev), abs (scaling(S_VELO,zlev)%data(d)%elts(id_e)))
       end if
    end do
  end subroutine linf_velo
end module lnorms_mod
