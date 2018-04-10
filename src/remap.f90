module remap_mod
  use adapt_mod
  use ops_mod
  use time_integr_mod
  use wavelet_mod
  implicit none
  
  integer                             :: order
  integer, dimension (:), allocatable :: stencil
  real(8), parameter                  :: ex_val = -1.0d13
contains
  subroutine remap_vertical_coordinates (set_thresholds)
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! interpolate on the full (not the perturbation) quantities
    ! Conserves mass, heat and momentum flux

    external :: set_thresholds
    
    integer   :: d, j, k, l

    if (rank.eq.0) write(6,*) "Remapping vertical coordinates"

    ! Ensure boundary values are up to date
    call update_array_bdry (sol, NONE)
    
    do k = 1, zlevels
       do d = 1, size(grid)
          exner_fun(k)%data(d)%elts = ex_val
       end do
    end do

    ! Remap on finest level
    call apply_onescale (remap_scalars, level_end, z_null, 0, 1)

    ! Remap scalars at coarser levels
    do l = level_end-1, level_start-1, -1
       call update_array_bdry (sol, l+1)

       ! Compute scalar wavelet coefficients
       do d = 1, size(grid)
          do k = 1, zlevels
             mass => sol(S_MASS,k)%data(d)%elts
             temp => sol(S_TEMP,k)%data(d)%elts
             wc_m => wav_coeff(S_MASS,k)%data(d)%elts
             wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             nullify (wc_m, wc_t, mass, temp)
          end do
       end do
       call update_array_bdry (wav_coeff(S_MASS:S_TEMP,:), l+1)

       ! Remap at level l (over-written if value available from restriction)
       call apply_onescale (remap_scalars, l, z_null, 0, 1)

       ! Restrict scalars (sub-sample and lift) and velocity (average) to coarser grid
       do d = 1, size(grid)
          do k = 1, zlevels
             mass => sol(S_MASS,k)%data(d)%elts
             temp => sol(S_TEMP,k)%data(d)%elts
             wc_m => wav_coeff(S_MASS,k)%data(d)%elts
             wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (restrict_scalar, grid(d), l, k, 0, 0)
             nullify (mass, temp, wc_m, wc_t)
          end do
       end do
    end do

    ! Remap velocity
    call apply_onescale (remap_velo,    level_end, z_null, 0, 0)
    do l = level_end-1, level_start-1, -1
       call update_array_bdry (sol, l+1)

       ! Remap at level l (over-written if value available from restriction)
       call apply_onescale (remap_velo,    l, z_null, 0, 0)

       ! Restrict scalars (sub-sample and lift) and velocity (average) to coarser grid
       do d = 1, size(grid)
          do k = 1, zlevels
             velo => sol(S_VELO,k)%data(d)%elts
             call apply_interscale_d (restrict_velo, grid(d), l, k, 0, 0)
             nullify (velo)
          end do
       end do
    end do
    sol%bdry_uptodate       = .False.
    wav_coeff%bdry_uptodate = .False.
        
    ! Re-adapt grid after remapping
    call WT_after_step (sol, wav_coeff, level_start-1)
    call adapt_grid (set_thresholds)
  end subroutine remap_vertical_coordinates

  subroutine remap_save 
    ! Remap the Lagrangian layers to pressure levels given by pressure_save array
    integer            :: d, j, k, l, p
    integer, parameter :: order_default = 7 ! order must be odd

    if (rank.eq.0) write(6,*) "Remapping vertical coordinates for saving"
    
    ! Set order of Newton interpolation
    order = min(zlevels+1, order_default)
    if (allocated(stencil)) deallocate(stencil)
    allocate(stencil(1:order))
    if (zlevels+1.lt.3) then
       write(6,*) "Cannot remap fewer than 3 vertical levels"
       stop
    end if

    ! Ensure boundary values are up to date
    call update_array_bdry (sol, NONE)
    sol_save = sol

    ! Remap at save level
    call apply_onescale (remap_vars_save, level_save, z_null, 0, 1)
    
    ! Remap poles
    do d = 1, size(grid)
       call apply_to_pole_d (remap_vars_save, grid(d), min_level-1, z_null, z_null, .True.)
    end do
    
    ! Ensure boundary values are up to date
    call update_array_bdry (sol_save, NONE)
  end subroutine remap_save
  
  subroutine remap_scalars (dom, i, j, z_null, offs, dims)
    ! Remap mass and vertical variable eta onto grid defined by A, B coefficients
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: current_zlev, d, id, id_i, k, zlev
    real(8)                          :: column_mass, cumul_mass_zlev, cumul_mass_target, new_cumul_mass, cumul_mass_upper, X
    real(8)                          :: new_mass
    real(8), dimension (zlevels+1)   :: cumul_mass, cumul_temp, new_cumul_temp

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    ! Calculate cumulative mass and total mass of column
    cumul_mass(1) = 0.0_8
    cumul_temp(1) = 0.0_8
    do k = 1, zlevels
       cumul_mass(k+1) = cumul_mass(k) + sol(S_MASS,k)%data(d)%elts(id_i)
       cumul_temp(k+1) = cumul_temp(k) + sol(S_TEMP,k)%data(d)%elts(id_i)
    end do
    column_mass = cumul_mass(zlevels+1)
    
    do k = 1, zlevels
      trend(S_MASS,k)%data(d)%elts(id_i) = sol(S_MASS,k)%data(d)%elts(id_i)
       new_mass = a_vert_mass(k) + b_vert_mass(k) * column_mass
       if (new_mass<0.0_8) then
          return ! Do not try to remap pole except at coarsest level
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = new_mass
       end if
    end do

    current_zlev = 1
    exner_fun(1)%data(d)%elts(id_i) = 1.0_8
    new_cumul_mass = 0.0_8
    do k = 1, zlevels
       cumul_mass_target = new_cumul_mass + sol(S_MASS,k)%data(d)%elts(id_i)

       do zlev = current_zlev, zlevels
          cumul_mass_upper = cumul_mass(zlev+1)
          if (cumul_mass_target <= cumul_mass_upper) exit
       end do

       if (zlev > zlevels) zlev = zlevels
       cumul_mass_zlev = cumul_mass(zlev)

       current_zlev = zlev
       new_cumul_mass = cumul_mass_target

       ! New vertical coordinate (saved in exner_fun)
       exner_fun(k+1)%data(d)%elts(id_i) = zlev + (cumul_mass_target - cumul_mass_zlev)/(cumul_mass_upper - cumul_mass_zlev)
    end do

    ! Remap theta
    do k = 1, zlevels+1
       X = exner_fun(k)%data(d)%elts(id_i)
       zlev = min(zlevels,floor(X))
       X = X - zlev
       new_cumul_temp(k) = cumul_temp(zlev) + X*sol(S_TEMP,zlev)%data(d)%elts(id_i)
    end do

    do k = 1, zlevels
       sol(S_TEMP,k)%data(d)%elts(id_i) = new_cumul_temp(k+1) - new_cumul_temp(k)
    end do
  end subroutine remap_scalars

  subroutine remap_velo (dom, i, j, z_null, offs, dims)
    ! Remap velocity onto original vertical grid by linear interpolation of mass flux
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: d, e, id, idE, idN, idNE, id_i, k, zlev
    real(8), dimension(EDGE) :: mass_e, X
    real(8), dimension(zlevels+1,1:EDGE) :: massflux_cumul, massflux, new_massflux_cumul

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    idN  = idx(i,   j+1, offs, dims) + 1
    idE  = idx(i+1, j,   offs, dims) + 1
    idNE = idx(i+1, j+1, offs, dims) + 1

    massflux_cumul(1,:) = 0.0_8
    do k = 1, zlevels
       ! Interpolate old masses (stored in trend)
       mass_e(RT+1) = trend(S_MASS,k)%data(d)%elts(id_i) + trend(S_MASS,k)%data(d)%elts(idE)
       mass_e(DG+1) = trend(S_MASS,k)%data(d)%elts(id_i) + trend(S_MASS,k)%data(d)%elts(idNE)
       mass_e(UP+1) = trend(S_MASS,k)%data(d)%elts(id_i) + trend(S_MASS,k)%data(d)%elts(idN)

       if (exner_fun(k)%data(d)%elts(idE).eq.ex_val)  mass_e(RT+1) = trend(S_MASS,k)%data(d)%elts(id_i)
       if (exner_fun(k)%data(d)%elts(idNE).eq.ex_val) mass_e(DG+1) = trend(S_MASS,k)%data(d)%elts(id_i)
       if (exner_fun(k)%data(d)%elts(idN).eq.ex_val)  mass_e(UP+1) = trend(S_MASS,k)%data(d)%elts(id_i)
          
       do e = 1, EDGE
          massflux(k,e) = sol(S_VELO,k)%data(d)%elts(EDGE*id+e) * mass_e(e)
          massflux_cumul(k+1,e) = massflux_cumul(k,e) + massflux(k,e)
       end do
    end do

    do k = 1, zlevels+1
       X(RT+1) = 0.5_8*(exner_fun(k)%data(d)%elts(id_i) + exner_fun(k)%data(d)%elts(idE))
       X(DG+1) = 0.5_8*(exner_fun(k)%data(d)%elts(id_i) + exner_fun(k)%data(d)%elts(idNE))
       X(UP+1) = 0.5_8*(exner_fun(k)%data(d)%elts(id_i) + exner_fun(k)%data(d)%elts(idN))

       if (exner_fun(k)%data(d)%elts(idE).eq.ex_val)  X(RT+1) = exner_fun(k)%data(d)%elts(id_i) 
       if (exner_fun(k)%data(d)%elts(idNE).eq.ex_val) X(DG+1) = exner_fun(k)%data(d)%elts(id_i) 
       if (exner_fun(k)%data(d)%elts(idN).eq.ex_val)  X(UP+1) = exner_fun(k)%data(d)%elts(id_i) 
       
       do e = 1, EDGE
          zlev = min(zlevels,floor(X(e)))
          X(e) = X(e) - zlev
          new_massflux_cumul(k,e) = massflux_cumul(zlev,e) + X(e)*massflux(zlev,e)
       end do
    end do

    do k = 1, zlevels
       ! Interpolate new masses
       mass_e(RT+1) = sol(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(idE)
       mass_e(DG+1) = sol(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(idNE)
       mass_e(UP+1) = sol(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(idN)
       
       if (exner_fun(k)%data(d)%elts(idE).eq.ex_val)  mass_e(RT+1) = sol(S_MASS,k)%data(d)%elts(id_i)
       if (exner_fun(k)%data(d)%elts(idNE).eq.ex_val) mass_e(DG+1) = sol(S_MASS,k)%data(d)%elts(id_i)
       if (exner_fun(k)%data(d)%elts(idN).eq.ex_val)  mass_e(UP+1) = sol(S_MASS,k)%data(d)%elts(id_i)
       
       do e = 1, EDGE
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (new_massflux_cumul(k+1,e) - new_massflux_cumul(k,e)) / mass_e(e)
       end do
    end do
  end subroutine remap_velo

  subroutine remap_vars_save (dom, i, j, zlev, offs, dims)
    ! Non-conservative interpolation of solution to pressure_save pressure levels
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                           :: d, e, id, id_i, k, kc, kk, m
    real(8)                           :: diff, dmin
    real(8), dimension (zlevels)      :: pressure, mass_interp, temp_interp
    real(8), dimension (EDGE,zlevels) :: velo_interp

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    do k = 1, zlevels
       mass_interp(k) = sol(S_MASS,k)%data(d)%elts(id_i)
       temp_interp(k) = sol(S_TEMP,k)%data(d)%elts(id_i)
       do e = 1, EDGE
          velo_interp(e,k) = sol(S_VELO,k)%data(d)%elts(EDGE*id+e) 
       end do
    end do

    ! Calculate pressure at each vertical level
    pressure(1) = dom%surf_press%elts(id_i) - 0.5_8*grav_accel*sol(S_MASS,1)%data(d)%elts(id_i)
    do k = 2, zlevels
       pressure(k) = pressure(k-1) - grav_accel*interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k-1)%data(d)%elts(id_i))
    end do

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    do k = 1, save_levels
       ! Find index of pressure on old vertical grid closest to pressure_save on new grid
       dmin = 1d16
       do kk = 1, zlevels
          diff = abs(pressure(kk)-pressure_save(k))
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do

       ! Set interpolation stencil based on pressure_save
       if (kc .lt. (order-1)/2+1) then
          stencil = (/ (m, m = 1, order) /)
       else if (kc .gt. zlevels-(order-1)/2) then
          stencil = (/ (m, m = zlevels-(order-1), zlevels) /)
       else
          stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
       end if

       ! Interpolate mass and temperature
       sol_save(S_MASS,k)%data(d)%elts(id_i) = Newton_interp(pressure(stencil), mass_interp(stencil), pressure_save(k))
       sol_save(S_TEMP,k)%data(d)%elts(id_i) = Newton_interp(pressure(stencil), temp_interp(stencil), pressure_save(k))
       
       ! Interpolate velocity
       do e = 1, EDGE
          sol_save(S_VELO,k)%data(d)%elts(EDGE*id+e) = Newton_interp(pressure(stencil), velo_interp(e,stencil), pressure_save(k))
       end do
    end do
  end subroutine remap_vars_save

  function Newton_interp (xv, yv, xd)
    ! Order point Newton form polynomial interpolation scheme as in Yang (2001)
    ! Uses the p values xv and yv to interpolate the value yd at given point xd
    ! order must be odd >= 3

    real(8)                               :: Newton_interp
    real(8), dimension(order), intent(in) :: xv, yv
    real(8),                   intent(in) :: xd

    real(8), dimension(order,order) :: interp_diff
    integer                         ::  mi, ni

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
