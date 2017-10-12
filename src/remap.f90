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
    !trend = sol
    
    do l = level_start, level_end
       ! Find mass, mass-weighted potential temperature and momentum on new vertical grid
       ! including nearest neighbours in N and E directions needed for momentum calculation
       call apply_onescale (remap_scalars, l, z_null, 0, 1)

       ! Find velocity on new vertical grid by interpolating momentum (no boundary points included)
       call apply_onescale (remap_velocity, l, z_null, 0, 0)
    end do

    ! Update boundary values
    call update_array_bdry (sol, NONE)

    ! cnt = 0
    ! do l = level_start, level_end
    !    call apply_onescale (compare, l, z_null, 0, 0)
    !    do zlev = 1, zlevels
    !       call apply_onescale (write_example, l, zlev, 0, 0)
    !    end do
    !  end do
    ! write(6,'(A,5(i5,1x))') "Incorrect at ", cnt(1:5)
  end subroutine remap_vertical_coordinates

  subroutine check_bdry()
    integer :: l, zlev

    zlev = 5

    val = 25.0_8
    do l = level_start, level_end
       call apply_onescale (change_value, l, zlev, 0, 0)
    end do
    
    call update_array_bdry (sol, NONE)

    do l = level_start, level_end
       cnt = 0
       call apply_onescale (check_mass, l, zlev, 0, 0)
       !call apply_onescale (check_velo, l, zlev, 0, 0)
       
       write (6,'(7(A,i4))') "Count = ", sum(cnt), " IdN = ", cnt(1), " IdE = ", cnt(2), " IdNE = ", cnt(3), &
            " IdS = ", cnt(4), " IdSW = ", cnt(5), " IdW = ", cnt(6)
    end do
  end subroutine check_bdry

   subroutine change_value (dom, i, j, zlev, offs, dims)
    ! Finds relative changes before and after interpolation
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer :: d, e, id, idN, idE, idNE, kb
    real(8), dimension(3) :: u

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    sol(S_MASS,zlev)%data(d)%elts(id+1) = val
  end subroutine change_value

 subroutine check_mass (dom, i, j, zlev, offs, dims)
    ! Finds relative changes before and after interpolation
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer :: d, e, id, idN, idE, idNE, idS, idSW, idW, kb, p
    real(8), dimension(3) :: u

    d = dom%id + 1
    
    id   = idx(i, j, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)

    
    if(sol(S_MASS,zlev)%data(d)%elts(idN+1).ne.val) then
       write(6,'(A, es11.4)') "idN", sol(S_MASS,zlev)%data(d)%elts(idN+1)
       cnt(1) = cnt(1)+1
    end if
    if(sol(S_MASS,zlev)%data(d)%elts(idE+1).ne.val) then
       write(6,'(A, es11.4)') "idE", sol(S_MASS,zlev)%data(d)%elts(idE+1)
       cnt(2)=cnt(2)+1
    end if
    if(sol(S_MASS,zlev)%data(d)%elts(idNE+1).ne.val) then
       write(6,'(A, es11.4)') "idNE", sol(S_MASS,zlev)%data(d)%elts(idNE+1)
       cnt(3)=cnt(3)+1
    end if
    if(sol(S_MASS,zlev)%data(d)%elts(idS+1).ne.val) then
       write(6,'(A, es11.4)') "idS", sol(S_MASS,zlev)%data(d)%elts(idS+1)
       cnt(4) = cnt(4)+1
    end if
    if(sol(S_MASS,zlev)%data(d)%elts(idSW+1).ne.val) then
       write(6,'(A, es11.4)') "idSW", sol(S_MASS,zlev)%data(d)%elts(idSW+1)
       cnt(5) = cnt(5)+1
    end if
    if(sol(S_MASS,zlev)%data(d)%elts(idW+1).ne.val) then
       write(6,'(A, es11.4)') "idW", sol(S_MASS,zlev)%data(d)%elts(idW+1)
       cnt(6) = cnt(6)+1
    end if
  end subroutine check_mass

   subroutine check_velo (dom, i, j, zlev, offs, dims)
    ! Finds relative changes before and after interpolation
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer :: d, e, id, idN, idE, idNE, idS, idSW, idW, kb, p
    real(8), dimension(3) :: u

    d = dom%id + 1
    id = idx(i, j, offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)

    do e = 1, 3
       if(abs(sol(S_VELO,zlev)%data(d)%elts(EDGE*id+e)).ne.5.0_8) cnt(e) = cnt(e)+1
    end do
    
    if(abs(sol(S_VELO,zlev)%data(d)%elts(EDGE*idW+RT+1)).ne.5.0_8) cnt(4) = cnt(4)+1
    if(abs(sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1)).ne.5.0_8) cnt(5) = cnt(5)+1
    if(abs(sol(S_VELO,zlev)%data(d)%elts(EDGE*idS+UP+1)).ne.5.0_8) cnt(6) = cnt(6)+1
  end subroutine check_velo
  
  subroutine compare (dom, i, j, zlev, offs, dims)
    ! Finds relative changes before and after interpolation
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer :: d, e, id, kb, k
    real(8), dimension (5) :: err

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    ! Calculate relative changes
    do kb = 1, zlevels
       err(1) = abs(sol(S_MASS,kb)%data(d)%elts(id+1) - trend(S_MASS,kb)%data(d)%elts(id+1))&
            /abs(trend(S_MASS,kb)%data(d)%elts(id+1))
       err(2) = abs(sol(S_TEMP,kb)%data(d)%elts(id+1) - trend(S_TEMP,kb)%data(d)%elts(id+1))&
            /abs(trend(S_TEMP,kb)%data(d)%elts(id+1))
       do e = 1, 3
          err(e+2) = abs(abs(sol(S_VELO,kb)%data(d)%elts(EDGE*id+e)) - abs(trend(S_VELO,kb)%data(d)%elts(EDGE*id+e)))
       end do
       do k = 1, 5
          if (err(k).gt.1e-5) cnt(k)=cnt(k)+1
       end do
          
          ! write (6,'(i2,1x,10(es11.4,1x))') kb, &
          !      trend(S_MASS,kb)%data(d)%elts(id+1), sol(S_MASS,kb)%data(d)%elts(id+1),&
          !      trend(S_TEMP,kb)%data(d)%elts(id+1), sol(S_TEMP,kb)%data(d)%elts(id+1),&
          !   trend(S_VELO,kb)%data(d)%elts(EDGE*id+RT+1), sol(S_VELO,kb)%data(d)%elts(EDGE*id+RT+1),&
          !   trend(S_VELO,kb)%data(d)%elts(EDGE*id+DG+1), sol(S_VELO,kb)%data(d)%elts(EDGE*id+DG+1),&
          !   trend(S_VELO,kb)%data(d)%elts(EDGE*id+UP+1), sol(S_VELO,kb)%data(d)%elts(EDGE*id+UP+1)
    end do
  end subroutine compare

   subroutine write_example (dom, i, j, zlev, offs, dims)
    ! Finds relative changes before and after interpolation
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer :: d, e, id, kb
    real(8), dimension (5) :: err

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    if (rank.eq.0 .and. d .eq. 7 .and.  id.eq. 110) then
       write(6,'(A,i4,1x,A,I2,1x,5(A,es11.4))') &
            ' id = ', id, &
            ' zlev = ', zlev, &
            ' mu =', sol(S_MASS,zlev)%data(d)%elts(id+1), &
            ' Theta =', sol(S_TEMP,zlev)%data(d)%elts(id+1), &
            ' U = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1), &
            ' V = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1), &
            ' W = ', sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1)
    end if
  end subroutine write_example

  subroutine remap_scalars (dom, i, j, zlev, offs, dims)
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer                        :: d, e, id, k, kb, kc, kk, m
    real(8)                        :: layer_pressure, p_surf, diff, dmin
    real(8), dimension (zlevels+1) :: integrated_mass, integrated_temp, new_mass, new_temp, pressure 

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    ! Integrate full mass and full mass-weighted potential temperature  vertically downward from the top
    ! All quantities located at interfaces
    integrated_mass(1) = 0.0_8
    integrated_temp(1) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 ! Actual zlevel
       
       if (k.eq.zlevels) then
          layer_pressure = press_infty + 0.5_8*grav_accel*sol(S_MASS,k)%data(d)%elts(id+1)
       else ! Interpolate mass to lower interface
          layer_pressure = layer_pressure + 0.5_8*grav_accel*(sol(S_MASS,k+1)%data(d)%elts(id+1) &
               + sol(S_MASS,k)%data(d)%elts(id+1))
       end if
       ! Integrate mass-weighted temperature
       integrated_mass(kb) = integrated_mass(kb-1) + sol(S_MASS,k)%data(d)%elts(id+1) + mean(S_MASS,k)
       integrated_temp(kb) = integrated_temp(kb-1) + sol(S_TEMP,k)%data(d)%elts(id+1)*(layer_pressure/ref_press)**kappa &
            + mean(S_TEMP,k)
    end do

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    ! pressure(zlevels+1) is surface pressure
    pressure(1) = press_infty
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel * (sol(S_MASS,k)%data(d)%elts(id+1) + mean(S_MASS,k))
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
   
       ! Interpolate integrated temperature and momentum at top interfaces of new vertical grid
       new_mass(kb) = Newton_interp(pressure(stencil), integrated_mass(stencil), layer_pressure)
       new_temp(kb) = Newton_interp(pressure(stencil), integrated_temp(stencil), layer_pressure)
    end do

    ! Cell values on new vertical grid from integrated values at interfaces calculated bottom to top
    do k = 1, zlevels
       sol(S_MASS,k)%data(dom%id+1)%elts(id+1) = (new_mass(zlevels-k+2) - new_mass(zlevels-k+1)) - mean(S_MASS,k)

       layer_pressure = 0.5_8*(a_vert(k)+a_vert(k+1))*ref_press + 0.5_8*(b_vert(k)+b_vert(k+1))*p_surf
       sol(S_TEMP,k)%data(d)%elts(id+1) = (new_temp(zlevels-k+2) - new_temp(zlevels-k+1))*(layer_pressure/ref_press)**(-kappa)&
            - mean(S_TEMP,k)
    end do
  end subroutine remap_scalars

  subroutine remap_velocity (dom, i, j, zlev, offs, dims)
    type (Domain)                 :: dom
    integer                       :: i, j, zlev
    integer, dimension (N_BDRY+1) :: offs
    integer, dimension (2,9)      :: dims

    integer                          :: d, e, id, idN, idE, idNE, k, kb, kc, kk, m
    real(8)                          :: layer_pressure, p_surf, diff, dmin
    real(8), dimension (3)           :: mu
    real(8), dimension (zlevels+1)   :: pressure 
    real(8), dimension (zlevels+1,3) :: integrated_momentum, new_momentum
    
    d = dom%id + 1
 
    id   = idx(i,     j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    
    ! Integrate full momentum vertically downward from the top
    ! Located at interfaces
    integrated_momentum(1,:) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 ! Actual zlevel

       ! Interpolate mass to edges (use masses on new grid)
       mu(RT+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id+1)+sol(S_MASS,k)%data(d)%elts(idE+1))
       mu(DG+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id+1)+sol(S_MASS,k)%data(d)%elts(idNE+1))
       mu(UP+1) = 0.5_8*(sol(S_MASS,k)%data(d)%elts(id+1)+sol(S_MASS,k)%data(d)%elts(idN+1))
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
       pressure(kb) = pressure(kb-1) + grav_accel * (sol(S_MASS,k)%data(d)%elts(id+1) + mean(S_MASS,k))
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

       do e = 1, 3
          new_momentum(kb,e) = Newton_interp(pressure(stencil), integrated_momentum(stencil,e), layer_pressure)
       end do
    end do

    ! Velocity from integrated momentum at interfaces calculated bottom to top
    do k = 1, zlevels
       do e = 1, 3
          sol(S_VELO,k)%data(dom%id+1)%elts(EDGE*id+e) = (new_momentum(zlevels-k+2,e) - new_momentum(zlevels-k+1,e))/mu(e)
       end do
    end do
  end subroutine remap_velocity
  
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
