module lnorms_mod
  use domain_mod
  use comm_mpi_mod
  implicit none
  integer :: n_norm
contains
  subroutine cal_lnorm (order)
    ! Calculates order = "l1", "l2", "linf", "lmin" norm of all prognostic variables
    ! result is divided by number of active cells for l1 and l2 norms to get average 
    implicit none
    character(*) :: order

    integer :: k, v

    ! Count number of active cells
    n_norm = 0
    call apply (count_norm, z_null)
    n_norm = sum_int (n_norm)

    ! Initialize norms
    select case (order)
    case ("1")
       lnorm = 0d0
    case ("2")
       lnorm = 0d0
    case ("inf")
       lnorm = 0d0
    case ("min")
       lnorm = 1d16
    end select

    do k = zmin, zmax
       select case (order)
       case ("1")
          call apply_bdry (l1_scalar, k, 0, 1)
          call apply_bdry (l1_velo,   k, 0, 0)
       case ("2")
          call apply_bdry (l2_scalar, k, 0, 1)
          call apply_bdry (l2_velo,   k, 0, 0)
       case ("inf")
          call apply_bdry (linf_scalar, k, 0, 1)
          call apply_bdry (linf_velo,   k, 0, 0)
       case ("min")
          call apply_bdry (lmin_scalar, k, 0, 1)
          call apply_bdry (lmin_velo,   k, 0, 0)
       end select
       select case (order)
       case ("1")
          lnorm(:,k) = sum_real (lnorm(:,k)) / dble (n_norm)
       case ("2")
          lnorm(:,k) = sqrt (sum_real (lnorm(:,k)) / dble (n_norm))
       case ("inf")
          lnorm(:,k) = sync_max_real (lnorm(:,k))
       case ("min")
          lnorm(:,k) = sync_min_real (lnorm(:,k))
       end select
    end do
  end subroutine cal_lnorm

  subroutine count_norm (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims)
    
    if (dom%mask_n%elts(id+1) >= ADJZONE) n_norm = n_norm + 1
  end subroutine count_norm

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
          lnorm(v,zlev) = lnorm(v,zlev) + abs (sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1))
       end do
    end if
  end subroutine l1_scalar

  subroutine l1_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e

    d  = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do e = 1, EDGE
          id_e = EDGE*id+e
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + abs (sol(S_VELO,zlev)%data(d)%elts(id_e))
       end do
    end if
  end subroutine l1_velo

  subroutine l2_scalar (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, v
    
    d  = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then 
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = lnorm(v,zlev) + (sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1))**2
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

    if (dom%mask_n%elts(id+1) >= ADJZONE) then 
       do e = 1, EDGE
          id_e = EDGE*id+e
          lnorm(S_VELO,zlev) = lnorm(S_VELO,zlev) + sol(S_VELO,zlev)%data(d)%elts(id_e)**2
       end do
    end if
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
          lnorm(v,zlev) = max (lnorm(v,zlev), abs (sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1)))
       end do
    end if
  end subroutine linf_scalar

  subroutine lmin_scalar (dom, i, j, zlev, offs, dims)
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
          lnorm(v,zlev) = min (lnorm(v,zlev), sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1))
       end do
    end if
  end subroutine lmin_scalar

  subroutine linf_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_e, e

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then 
       do e = 1, EDGE
          id_e = EDGE*id+e
          lnorm(S_VELO,zlev) = max (lnorm(S_VELO,zlev), abs (sol(S_VELO,zlev)%data(d)%elts(id_e)))
       end do
    end if
  end subroutine linf_velo

   subroutine lmin_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_e, e

    d = dom%id+1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then 
       do e = 1, EDGE
          id_e = EDGE*id+e
          lnorm(S_VELO,zlev) = min (lnorm(S_VELO,zlev), sol(S_VELO,zlev)%data(d)%elts(id_e))
       end do
    end if
  end subroutine lmin_velo
end module lnorms_mod
