module lnorms_mod

  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, EDGE, N_BDRY, lnorm, S_TEMP, S_MASS, S_VELO, scalars
  use domain_mod, only : Domain, sol, sol_mean, idx
 
  implicit none

  private
  public :: cal_lnorm, lnorm

contains

  subroutine cal_lnorm (order)
    ! Computes order = "1", "2", "inf", or "min" norm of prognostic variables.
    !
    ! The l1 and l2 norms are computed as area integrals over the adaptive triangular grid.
    
    use shared_mod,     only : MATH_PI, radius, zmin, zmax
    use comm_mpi_mod,   only : sum_int, sum_real, sync_max_real, sync_min_real
    use domain_ops_mod, only : apply, apply_bdry
    use integrate_mod,  only : integrate_tri
    
    implicit none

    character(len=*), intent(in) :: order

    integer  :: k
    real(dp) :: Area
    
    Area = 4.0_dp*MATH_PI*radius**2

    ! Initialize norms.
    select case (order)
    case ("1", "2", "inf")
       lnorm = 0.0_dp
    case ("min")
       lnorm = 1.0e16_dp
    end select

    do k = zmin, zmax
       select case (order)
       case ("1")
          lnorm(S_MASS, k) = integrate_tri (l1_mass_fun, k) / Area
          lnorm(S_TEMP, k) = integrate_tri (l1_temp_fun, k) / Area
          lnorm(S_VELO, k) = integrate_tri (l1_velo_fun, k) / Area
       case ("2")
          lnorm(S_MASS, k) = sqrt (integrate_tri (l2_mass_fun, k) / Area)
          lnorm(S_TEMP, k) = sqrt (integrate_tri (l2_temp_fun, k) / Area)
          lnorm(S_VELO, k) = sqrt (integrate_tri (l2_velo_fun, k) / Area)
       case ("inf")
          call apply_bdry (linf_scalar, k, 0, 1)
          call apply_bdry (linf_velo,   k, 0, 0)
          
          lnorm(:,k) = sync_max_real(lnorm(:,k))
       case ("min")
          call apply_bdry (lmin_scalar, k, 0, 1)
          call apply_bdry (lmin_velo,   k, 0, 0)
          
          lnorm(:,k) = sync_min_real(lnorm(:,k))
       end select
    end do
  end subroutine cal_lnorm


  function l1_mass_fun (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    val = abs(sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol(S_MASS,zlev)%data(d)%elts(id+1))
  end function l1_mass_fun

  
  function l1_temp_fun (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    val = abs(sol_mean(S_TEMP,zlev)%data(d)%elts(id+1) + sol(S_TEMP,zlev)%data(d)%elts(id+1))
  end function l1_temp_fun


  function l1_velo_fun (dom, i, j, zlev, offs, dims) result(val)
    use diagnostics_mod, only : u_mag
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: id
    
    id = idx(i, j, offs, dims)

    val = u_mag (dom, i, j, zlev, offs, dims)
  end function l1_velo_fun
  
  
  function l2_mass_fun (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    val = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol(S_MASS,zlev)%data(d)%elts(id+1))**2 
  end function l2_mass_fun

  
  function l2_temp_fun (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    val = (sol_mean(S_TEMP,zlev)%data(d)%elts(id+1) + sol(S_TEMP,zlev)%data(d)%elts(id+1))**2 
  end function l2_temp_fun


  function l2_velo_fun (dom, i, j, zlev, offs, dims) result(val)
    use diagnostics_mod, only : u_mag
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: id

    id = idx(i, j, offs, dims)

    val = u_mag (dom, i, j, zlev, offs, dims)**2
  end function l2_velo_fun


  subroutine linf_scalar (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, v

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = max ( &
               lnorm(v,zlev), &
               abs (sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1)))
       end do
    end if
  end subroutine linf_scalar


  subroutine lmin_scalar (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, v

    d  = dom%id + 1
    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       do v = scalars(1), scalars(2)
          lnorm(v,zlev) = min ( &
               lnorm(v,zlev), &
               abs(sol_mean(v,zlev)%data(d)%elts(id+1) + sol(v,zlev)%data(d)%elts(id+1)))
       end do
    end if
  end subroutine lmin_scalar


  subroutine linf_velo (dom, i, j, zlev, offs, dims)
    use diagnostics_mod, only : u_mag
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) &
         lnorm(S_VELO,zlev) = max (lnorm(S_VELO,zlev), u_mag (dom, i, j, zlev, offs, dims))
  end subroutine linf_velo


  subroutine lmin_velo (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    
    integer :: id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) &
         lnorm(S_VELO,zlev) = min (lnorm(S_VELO,zlev), u_mag (dom, i, j, zlev, offs, dims))
  end subroutine lmin_velo

  
end module lnorms_mod

