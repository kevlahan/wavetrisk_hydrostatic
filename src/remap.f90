module remap_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use comm_mpi_mod
  use ops_mod
  use wavelet_mod
  use time_integr_mod
  implicit none
  integer                             :: order
  integer, dimension (:), allocatable :: stencil
contains
  subroutine remap_vertical_coordinates 
    ! Remap the Lagrangian layers to target vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! Conserves mass, heat and momentum flux
    integer            :: l
    integer, parameter :: order_default = 7 ! order must be odd

    ! Set order of Newton interpolation
    order = min (zlevels+1, order_default)
    if (allocated (stencil)) deallocate (stencil)
    allocate (stencil(1:order))
    if (zlevels+1 < 3) then
       write (6,'(A)') "Cannot remap fewer than 3 vertical levels"
       stop
    end if

    ! Ensure boundary values are up to date
    call update_array_bdry (sol, NONE)

    do l = level_start, level_end
       call apply_onescale (remap_scalars,  l, z_null, 0, 1)
       call apply_onescale (remap_velocity, l, z_null, 0, 0)
    end do

    ! Re-adapt grid after remapping
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine remap_vertical_coordinates

  subroutine remap_scalars (dom, i, j, zlev, offs, dims)
    ! Remap scalars to target grid
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, id, id_i, k, kb
    real(8)                          :: new_p, p_s
    real(8), dimension (zlevels+1)   :: cumul_temp, new_cumul_temp, p

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    
    call cal_p (p, p_s, d, id_i)

    ! Integrate full full mass-weighted potential temperature vertically downward from the top
    ! All quantities located at interfaces
    cumul_temp(1) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 
       cumul_temp(kb) = cumul_temp(kb-1) + sol(S_TEMP,k)%data(d)%elts(id_i)
       trend(S_MASS,k)%data(d)%elts(id_i) = sol(S_MASS,k)%data(d)%elts(id_i) ! Save current mass for momentum interpolation
    end do

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    new_cumul_temp(1) = 0.0_8
    do kb = 2, zlevels
       k = zlevels-kb+2
       new_p = a_vert(k) + b_vert(k)*p_s
       call find_stencil (p, new_p)
       new_cumul_temp(kb) = Newton_interp (p(stencil), cumul_temp(stencil), new_p)
    end do
    new_cumul_temp(zlevels+1) = cumul_temp(zlevels+1)

    ! Variables on new vertical grid
    do k = 1, zlevels
       kb = zlevels-k+1
       sol(S_MASS,k)%data(d)%elts(id_i) = a_vert_mass(k) + b_vert_mass(k)*p_s/grav_accel
       sol(S_TEMP,k)%data(d)%elts(id_i) = new_cumul_temp(kb+1) - new_cumul_temp(kb)
    end do
  end subroutine remap_scalars

  subroutine remap_velocity (dom, i, j, zlev, offs, dims)
    ! Remap velocity to target vertical grid
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                               :: d, e, id, id_i, k, kb
    integer, dimension(1:EDGE)            :: idr

    real(8)                               :: new_p, p_s
    real(8), dimension (zlevels+1)        :: p
    real(8), dimension (zlevels+1,1:EDGE) :: cumul_flux, new_cumul_flux

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    idr(RT+1) = idx (i,   j+1, offs, dims) + 1
    idr(DG+1) = idx (i+1, j,   offs, dims) + 1
    idr(UP+1) = idx (i+1, j+1, offs, dims) + 1

    call cal_p (p, p_s, d, id_i)

    ! Integrate full momentum flux vertically downward from the top
    ! All quantities located at interfaces
    cumul_flux(1,:) = 0.0_8
    do kb = 2, zlevels+1
       k = zlevels-kb+2 ! Actual zlevel
       do e = 1, EDGE
          cumul_flux(kb,e) = cumul_flux(kb-1,e) + sol(S_VELO,k)%data(d)%elts(EDGE*id+e) &
               * (trend(S_MASS,k)%data(d)%elts(id_i) + trend(S_MASS,k)%data(d)%elts(idr(e)))
       end do
    end do

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    new_cumul_flux(1,:) = 0.0_8
    do kb = 2, zlevels
       k = zlevels-kb+2
       new_p = a_vert(k) + b_vert(k)*p_s
       call find_stencil (p, new_p)
       do e = 1, EDGE
          new_cumul_flux(kb,e) = Newton_interp (p(stencil), cumul_flux(stencil,e), new_p)
       end do
    end do
    new_cumul_flux(zlevels+1,:) = cumul_flux(zlevels+1,:)

    ! Find velocity on new grid from mass flux
    do k = 1, zlevels
       kb = zlevels-k+1
       do e = 1, EDGE
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (new_cumul_flux(kb+1,e) - new_cumul_flux(kb,e)) &
               / (sol(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(idr(e)))
       end do
    end do
  end subroutine remap_velocity

  subroutine cal_p (p, p_s, d, id_i)
    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    implicit none
    integer                        :: d, id_i
    real(8)                        :: p_s
    real(8), dimension (zlevels+1) :: p

    integer :: k, kb
    
    p(1) = p_top
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       p(kb) = p(kb-1) + grav_accel * sol(S_MASS,k)%data(d)%elts(id_i)
    end do
    p_s = p(zlevels+1)
  end subroutine cal_p

  subroutine find_stencil (p, new_p)
    ! Find stencil associated with new pressure level
    implicit none
    real(8)                        :: new_p
    real(8), dimension (zlevels+1) :: p

    integer                        :: kc, kk, m
    real(8)                        :: diff, dmin
    
    dmin = 1d16
    do kk = 1, zlevels+1
       diff = abs (p(kk)-new_p)
       if (diff < dmin) then
          kc = kk
          dmin = diff
       end if
    end do

    if (kc < (order-1)/2+1) then
       stencil = (/ (m, m = 1, order) /)
    else if (kc > zlevels+1-(order-1)/2) then
       stencil = (/ (m, m = zlevels+1-(order-1), zlevels+1) /)
    else
       stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
    end if
  end subroutine find_stencil

  function Newton_interp (xv, yv, xd)
    ! Order point Newton form polynomial interpolation scheme as in Yang (2001)
    ! Uses the p values xv and yv to interpolate the value yd at given point xd
    ! order must be odd >= 3

    real(8)                               :: Newton_interp
    real(8), dimension(order), intent(in) :: xv, yv
    real(8),                   intent(in) :: xd

    real(8), dimension(order,order) :: interp_diff
    integer                         :: mi, ni

    ! Construct interpolating polynomial by calculating Newton differences
    interp_diff(1,1:order) = yv ! zeroth order finite differences for x_0 thru x_6
    do mi = 1, order-1
       do ni = 0, (order-1)-mi
          interp_diff(mi+1,ni+1) = (interp_diff(mi,ni+2)-interp_diff(mi,ni+1))/(xv(ni+mi+1)-xv(ni+1))
       end do
    end do

    ! Evaluate the polynomial using Horner's algorithm
    Newton_interp = interp_diff(order,1)
    do mi = order-2, 0, -1
       Newton_interp = interp_diff(mi+1,1)+(xd-xv(mi+1)) * Newton_interp
    end do
  end function Newton_interp
end module remap_mod
