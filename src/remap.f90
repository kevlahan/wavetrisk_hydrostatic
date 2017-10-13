module remap_mod
  use geom_mod
  use domain_mod
  use arch_mod
  use comm_mpi_mod
  implicit none
  logical :: wasprinted1
  integer :: p
  integer, dimension (:), allocatable :: stencil
  integer, dimension(6) :: cnt
  real(8) :: val
contains
  subroutine remap_vertical_coordinates()
    ! Remap the Lagrangian coordinate given a_vert and b_vert vertical coordinate parameters 
    ! interpolate on the full (not the perturbation) quantities
    ! Conserves mass, heat and momentum
    !
    ! Assumes mean velocity is zero
    
    integer :: l, zlev
    integer, parameter :: p_default=7 ! p must be odd

    ! Set order of Newton interpolation
    p = min(zlevels+1, p_default)
    if (allocated(stencil)) deallocate(stencil)
    allocate(stencil(1:p))
    if (zlevels+1.lt.3) then
       write(6,*) "Cannot remap fewer than 3 vertical levels"
       stop
    end if

    ! Ensure boundary values are up to date
    call update_array_bdry (sol, NONE)

    do l = level_start, level_end
       ! Find mass, mass-weighted potential temperature and momentum on new vertical grid
       ! including nearest neighbours in N and E directions needed for momentum calculation
       call apply_onescale (remap_solution, l, z_null, 0, 1)
    end do

    ! Update boundary values
    call update_array_bdry (sol, NONE)
  end subroutine remap_vertical_coordinates

  subroutine remap_solution (dom, i, j, zlev, offs, dims)
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer                        :: d, e, id, id_i, idN, idE, idNE, k, kb, kc, kk, m
    real(8)                        :: layer_pressure, p_surf, diff, dmin, vel_old
    real(8), dimension (zlevels+1) :: integrated_mass, integrated_temp, new_mass, new_temp, pressure
    real(8), dimension (3)           :: mu
    real(8), dimension (zlevels+1,3) :: integrated_momentum, new_momentum

    d = dom%id + 1
    id = idx(i, j, offs, dims)
    id_i = id+1

    idN  = idx(i,     j + 1, offs, dims) + 1
    idE  = idx(i + 1, j,     offs, dims) + 1
    idNE = idx(i + 1, j + 1, offs, dims) + 1

    ! Integrate full mass, full mass-weighted potential temperature and full momentum vertically downward from the top
    ! All quantities located at interfaces
    integrated_mass(1)       = 0.0_8
    integrated_temp(1)       = 0.0_8
    integrated_momentum(1,:) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 ! Actual zlevel

       integrated_mass(kb) = integrated_mass(kb-1) + sol(S_MASS,k)%data(d)%elts(id_i) + mean(S_MASS,k)
       integrated_temp(kb) = integrated_temp(kb-1) + sol(S_TEMP,k)%data(d)%elts(id_i) + mean(S_TEMP,k)
 
       ! Interpolate mass to edges
       mu(RT+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idE))
       mu(DG+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idNE))
       mu(UP+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idN))
       mu = mu + mean(S_MASS,k)
       ! Integrate momentum at vertical edges
       do e = 1, 3
          integrated_momentum(kb,e) = integrated_momentum(kb-1,e) + sol(S_VELO,k)%data(d)%elts(EDGE*id+e)*mu(e)
       end do
    end do

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    ! pressure(zlevels+1) is surface pressure
    pressure(1) = press_infty
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel * (sol(S_MASS,k)%data(d)%elts(id_i) + mean(S_MASS,k))
    end do
    p_surf = pressure(zlevels+1)
    
    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    do kb = 1, zlevels+1
       k = zlevels-kb+2
       ! Pressure at top interface of cells of new vertical grid
       layer_pressure = a_vert(k)*ref_press + b_vert(k)*p_surf

       ! Find index of pressure on old vertical grid closest to layer_pressure on new grid
       dmin = 1e16_8
       do kk = 1, zlevels+1
          diff = abs(pressure(kk)-layer_pressure)
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do
       
       ! Set interpolation stencil based on layer_pressure
       if (kc .lt. (p-1)/2+1) then
          stencil = (/ (m, m = 1, p) /)
       else if (kc .gt. zlevels+1-(p-1)/2) then
          stencil = (/ (m, m = zlevels-p+2, zlevels+1) /)
       else
          stencil = (/ (m, m = kc-(p-1)/2, kc+(p-1)/2) /)
       end if
   
       ! Interpolate mass, integrated temperature and momentum at top interfaces of new vertical grid
       new_mass(kb) = Newton_interp(pressure(stencil), integrated_mass(stencil), layer_pressure)
       new_temp(kb) = Newton_interp(pressure(stencil), integrated_temp(stencil), layer_pressure)
       do e = 1, 3
          new_momentum(kb,e) = Newton_interp(pressure(stencil), integrated_momentum(stencil,e), layer_pressure)
       end do
    end do

    ! Cell values on new vertical grid from integrated values at interfaces calculated bottom to top
    do k = 1, zlevels
       sol(S_MASS,k)%data(d)%elts(id_i) = (new_mass(zlevels-k+2) - new_mass(zlevels-k+1)) - mean(S_MASS,k)
       sol(S_TEMP,k)%data(d)%elts(id_i) = (new_temp(zlevels-k+2) - new_temp(zlevels-k+1)) - mean(S_TEMP,k)
       
       ! Interpolate new masses to edges
       mu(RT+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idE))
       mu(DG+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idNE))
       mu(UP+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id_i)+sol(S_MASS,k)%data(d)%elts(idN))
       mu = mu + mean(S_MASS,k)
       do e = 1, 3
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (new_momentum(zlevels-k+2,e) - new_momentum(zlevels-k+1,e))/mu(e)
       end do
    end do
  end subroutine remap_solution
  
  function Newton_interp (xv, yv, xd)
    ! p-point Newton form polynomial interpolation scheme as in Yang (2001)
    ! Uses the p values xv and yv to interpolate the value yd at given point xd
    ! p must be odd >= 3
    
    real(8)                           :: Newton_interp
    real(8), dimension(p), intent(in) :: xv, yv
    real(8),               intent(in) :: xd

    real(8), dimension(p,p) :: interp_diff
    integer                 ::  mi, ni

    ! Construct interpolating polynomial by calculating Newton differences
    interp_diff(1,:) = yv ! zeroth order finite differences for x_0 thru x_6
    do mi = 1, p-1
       do ni = 0, (p-1)-mi
          interp_diff(mi+1,ni+1) = (interp_diff(mi,ni+2)-interp_diff(mi,ni+1))/(xv(ni+mi+1)-xv(ni+1))
       end do
    end do

    ! Evaluate the polynomial using Horner's algorithm
    Newton_interp = interp_diff(p,1)
    do mi = p-2, 0, -1
       Newton_interp = interp_diff(mi+1,1)+(xd-xv(mi+1)) * Newton_interp
    end do
  end function Newton_interp
end module remap_mod
