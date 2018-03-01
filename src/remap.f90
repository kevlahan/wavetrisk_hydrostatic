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
  subroutine remap_vertical_coordinates (set_thresholds)
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! interpolate on the full (not the perturbation) quantities
    ! Conserves mass, heat and momentum flux
    external           :: set_thresholds
    integer            :: d, j, k, l, p
    integer, parameter :: order_default = 7 ! order must be odd

    if (rank.eq.0) write(6,*) "Remapping vertical coordinates"

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

    ! Remap on finest scale
    call apply_onescale (remap_velocity_fv, level_end, z_null, 0, 0)
    call apply_onescale (remap_scalars_fv,  level_end, z_null, 0, 0)

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
       call apply_onescale (remap_velocity_fv, l, z_null, 0, 0)
       call apply_onescale (remap_scalars_fv,  l, z_null, 0, 0)

       ! Restrict scalars (sub-sample and lift) and velocity (average) to coarser grid
       do d = 1, size(grid)
          do k = 1, zlevels
             mass => sol(S_MASS,k)%data(d)%elts
             temp => sol(S_TEMP,k)%data(d)%elts
             velo => sol(S_VELO,k)%data(d)%elts
             wc_m => wav_coeff(S_MASS,k)%data(d)%elts
             wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (restrict_scalar, grid(d), l, k, 0, 0)
             call apply_interscale_d (restrict_velo,   grid(d), l, k, 0, 0)
             nullify (mass, temp, velo, wc_m, wc_t)
          end do
       end do
    end do

    sol%bdry_uptodate       = .False.
    wav_coeff%bdry_uptodate = .False.
    
    ! Remap poles
    do d = 1, size(grid)
       call apply_to_pole_d (remap_scalars, grid(d), min_level-1, z_null, z_null, .True.)
    end do
    
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

  subroutine remap_velocity (dom, i, j, zlev, offs, dims)
    ! Remap velocity to original vertical grid
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, idE, idN, idNE, k, kb, kc, kk, m
    real(8)                          :: new_pressure, p_surf, p_surf_E, p_surf_N, p_surf_NE, diff, dmin, vel_old
    real(8)                          :: mass_id, mass_idE, mass_idN, mass_idNE
    real(8), dimension (3)           :: mass_e
    real(8), dimension (zlevels+1)   :: pressure
    real(8), dimension (zlevels+1,3) :: old_flux, new_flux

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    idN  = idx(i,   j+1, offs, dims) + 1
    idE  = idx(i+1, j,   offs, dims) + 1
    idNE = idx(i+1, j+1, offs, dims) + 1

    ! Integrate full momentum flux vertically downward from the top
    ! All quantities located at interfaces
    old_flux(1,:) = 0.0_8
    do kb = 2, zlevels+1
       k = zlevels-kb+2 ! Actual zlevel

       ! Interpolate mass on current vertical grid to edges 
       mass_e(RT+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idE))
       mass_e(DG+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idNE))
       mass_e(UP+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idN))
       ! Mass fluxes
       do e = 1, EDGE
          old_flux(kb,e) = old_flux(kb-1,e) + sol(S_VELO,k)%data(d)%elts(EDGE*id+e) * mass_e(e)
       end do
    end do

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    ! also calculate surface pressure at adjacent nodes needed to interpolate mass to edges
    pressure(1) = press_infty
    p_surf_E    = press_infty
    p_surf_NE   = press_infty
    p_surf_N    = press_infty
    do kb = 2, zlevels+1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel*sol(S_MASS,k)%data(d)%elts(id_i)
       p_surf_E  = p_surf_E  + grav_accel * sol(S_MASS,k)%data(d)%elts(idE) 
       p_surf_NE = p_surf_NE + grav_accel * sol(S_MASS,k)%data(d)%elts(idNE)
       p_surf_N  = p_surf_N  + grav_accel * sol(S_MASS,k)%data(d)%elts(idN) 
    end do
    p_surf = pressure(zlevels+1)
    
    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    new_flux(1,:) = 0.0_8
    do kb = 2, zlevels
       k = zlevels-kb+2
       ! Pressure at bottom interface of cells of new vertical grid
       new_pressure = a_vert(k)*ref_press + b_vert(k)*p_surf

       ! Find index of pressure on old vertical grid closest to new_pressure on new grid
       dmin = 1.0d16
       do kk = 1, zlevels+1
          diff = abs(pressure(kk)-new_pressure)
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do

       ! Set interpolation stencil based on new_pressure
       if (kc .lt. (order-1)/2+1) then
          stencil = (/ (m, m = 1, order) /)
       else if (kc .gt. zlevels+1-(order-1)/2) then
          stencil = (/ (m, m = zlevels+1-(order-1), zlevels+1) /)
       else
          stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
       end if

       ! Interpolate integrated mass flux at top interfaces of new vertical grid
       do e = 1, EDGE
          new_flux(kb,e) = Newton_interp(pressure(stencil), old_flux(stencil,e), new_pressure)
       end do
    end do
    new_flux(zlevels+1,:) = old_flux(zlevels+1,:)

    ! Variables on new vertical grid
    do k = 1, zlevels
       kb = zlevels-k+1

       ! Remapped masses needed to interpolate mass at edges
       mass_id   = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf)   /grav_accel
       mass_idE  = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_E) /grav_accel
       mass_idNE = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_NE)/grav_accel
       mass_idN  = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_N) /grav_accel

       ! Interpolate remapped masses to edges
       mass_e(RT+1) = interp(mass_id, mass_idE)
       mass_e(DG+1) = interp(mass_id, mass_idNE)
       mass_e(UP+1) = interp(mass_id, mass_idN)

       ! Find velocity on new grid from mass flux
       do e = 1, EDGE
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (new_flux(kb+1,e) - new_flux(kb,e)) / mass_e(e)
       end do
    end do
  end subroutine remap_velocity

  subroutine remap_scalars (dom, i, j, zlev, offs, dims)
    ! Remap scalars to original grid
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, k, kb, kc, kk, m
    real(8)                          :: new_pressure, p_surf, diff, dmin
    real(8), dimension (zlevels+1)   :: integrated_temp, new_temp, pressure

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    ! Integrate full full mass-weighted potential temperature vertically downward from the top
    ! All quantities located at interfaces
    integrated_temp(1) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 ! Actual zlevel
       integrated_temp(kb) = integrated_temp(kb-1) + sol(S_TEMP,k)%data(d)%elts(id_i)
    end do

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    pressure(1) = press_infty
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel * sol(S_MASS,k)%data(d)%elts(id_i)
    end do
    p_surf = pressure(zlevels+1)

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    new_temp(1) = 0.0_8
    do kb = 2, zlevels
       k = zlevels-kb+2
       ! Pressure at bottom interface of cells of new vertical grid
       new_pressure = a_vert(k)*ref_press + b_vert(k)*p_surf

       ! Find index of pressure on old vertical grid closest to new_pressure on new grid
       dmin = 1.0d16
       do kk = 1, zlevels+1
          diff = abs(pressure(kk)-new_pressure)
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do

       ! Set interpolation stencil based on new_pressure
       if (kc .lt. (order-1)/2+1) then
          stencil = (/ (m, m = 1, order) /)
       else if (kc .gt. zlevels+1-(order-1)/2) then
          stencil = (/ (m, m = zlevels+1-(order-1), zlevels+1) /)
       else
          stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
       end if

       ! Interpolate integrated temperature at top interfaces of new vertical grid
       new_temp(kb) = Newton_interp(pressure(stencil), integrated_temp(stencil), new_pressure)
    end do
    new_temp(zlevels+1) = integrated_temp(zlevels+1)

    ! Variables on new vertical grid
    do k = 1, zlevels
       kb = zlevels-k+1
       ! Remapped mass-weighted potential temperature from integrated value interpolated to new grid
       sol(S_TEMP,k)%data(d)%elts(id_i) = new_temp(kb+1) - new_temp(kb)
       sol(S_MASS,k)%data(d)%elts(id_i) = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf)/grav_accel
    end do
  end subroutine remap_scalars

   subroutine remap_velocity_fv (dom, i, j, zlev, offs, dims)
    ! Remap velocity to original grid using low order finite volume scheme
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, idE, idN, idNE, k, kb
    real(8)                          :: mass_id, mass_idE, mass_idN, mass_idNE
    real(8)                          :: new_pressure, p_surf, p_surf_E, p_surf_NE, p_surf_N
    real(8), dimension (3)           :: mass_e
    real(8), dimension (zlevels+1)   :: pressure
    real(8), dimension (zlevels+1,3) :: integrated_flux

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    idN  = idx(i,   j+1, offs, dims) + 1
    idE  = idx(i+1, j,   offs, dims) + 1
    idNE = idx(i+1, j+1, offs, dims) + 1

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    ! also calculate surface pressure at adjacent nodes needed to interpolate mass to edges
    pressure(1) = press_infty
    p_surf_E  = press_infty
    p_surf_NE = press_infty
    p_surf_N  = press_infty
    do kb = 2, zlevels+1
       pressure(kb) = pressure(kb-1) + grav_accel * sol(S_MASS,zlevels-kb+2)%data(d)%elts(id_i)
       p_surf_E     = p_surf_E       + grav_accel * sol(S_MASS,zlevels-kb+2)%data(d)%elts(idE) 
       p_surf_NE    = p_surf_NE      + grav_accel * sol(S_MASS,zlevels-kb+2)%data(d)%elts(idNE)
       p_surf_N     = p_surf_N       + grav_accel * sol(S_MASS,zlevels-kb+2)%data(d)%elts(idN) 
    end do
    p_surf = pressure(zlevels+1)

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    integrated_flux = 0.0_8
    do kb = 2, zlevels ! kb is lower interface on new grid
       ! Pressure at lower interface on new vertical grid
       new_pressure = a_vert(zlevels-kb+2)*ref_press + b_vert(zlevels-kb+2)*p_surf

       ! Integrate temperature downwards to new pressure level
       k = 2 ! k is lower interface on old grid
       do while (pressure(k) .lt. new_pressure)
          ! Interpolate mass on current vertical grid to edges 
          mass_e(RT+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idE))
          mass_e(DG+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idNE))
          mass_e(UP+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idN))

          ! Mass fluxes
          do e = 1, EDGE
             integrated_flux(kb,e) = integrated_flux(kb,e) + mass_e(e) * sol(S_VELO,zlevels-k+2)%data(d)%elts(EDGE*id+e)
          end do
          k = k+1
       end do
       ! Add additional contribution from partial level from level kb to layer pressure by linearly interpolating the integrated flux
       mass_e(RT+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idE))
       mass_e(DG+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idNE))
       mass_e(UP+1) = interp(sol(S_MASS,zlevels-k+2)%data(d)%elts(id_i), sol(S_MASS,zlevels-k+2)%data(d)%elts(idN))
       do e = 1, EDGE
          integrated_flux(kb,e) = integrated_flux(kb,e) &
               + (new_pressure-pressure(k-1))/(pressure(k)-pressure(k-1)) &
               * mass_e(e) * sol(S_VELO,zlevels-k+2)%data(d)%elts(EDGE*id+e)
       end do
    end do

    ! Flux integrated over all vertical levels (does not change)
    do k = 1, zlevels
       mass_e(RT+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idE))
       mass_e(DG+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idNE))
       mass_e(UP+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idN))
       do e = 1, EDGE
          integrated_flux(zlevels+1,e) = integrated_flux(zlevels+1,e) + mass_e(e) * sol(S_VELO,k)%data(d)%elts(EDGE*id+e)
       end do
    end do
      
    ! Velocity on new vertical grid
    do k = 1, zlevels
       ! Remapped masses needed to interpolate mass at edges
       mass_id   = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf)   /grav_accel
       mass_idE  = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_E) /grav_accel
       mass_idNE = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_NE)/grav_accel
       mass_idN  = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf_N) /grav_accel

       ! Interpolate remapped masses to edges
       mass_e(RT+1) = interp(mass_id, mass_idE)
       mass_e(DG+1) = interp(mass_id, mass_idNE)
       mass_e(UP+1) = interp(mass_id, mass_idN)

       ! Find velocity on new grid from mass flux
       do e = 1, EDGE
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (integrated_flux(zlevels-k+2,e) - integrated_flux(zlevels-k+1,e)) / mass_e(e)
       end do
    end do
  end subroutine remap_velocity_fv

  subroutine remap_scalars_fv (dom, i, j, zlev, offs, dims)
    ! Remap scalars to original grid using low order finite volume scheme
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, k, kb
    real(8)                          :: new_pressure, p_surf, m, P1, P2
    real(8), dimension (zlevels+1)   :: integrated_temp, pressure

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    ! Calculate pressure at interfaces of current vertical grid integrating from top down
    pressure(1) = press_infty
    do kb = 2, zlevels+1
       pressure(kb) = pressure(kb-1) + grav_accel * sol(S_MASS,zlevels-kb+2)%data(d)%elts(id_i)
    end do
    p_surf = pressure(zlevels+1)
    
    ! Compute integrated temperature at each interface downward from top
    integrated_temp = 0.0_8
    do kb = 2, zlevels ! kb is lower interface on new grid
       ! Pressure at lower interface on new vertical grid
       new_pressure = a_vert(zlevels-kb+2)*ref_press + b_vert(zlevels-kb+2)*p_surf

       ! Integrate temperature downwards to new pressure level
       k = 2 ! k is lower interface on old grid
       do while (pressure(k) .lt. new_pressure)
          integrated_temp(kb) = integrated_temp(kb) + sol(S_TEMP,zlevels-k+2)%data(d)%elts(id_i)
          k = k+1
       end do
       ! Add additional contribution from partial level from level kb to layer pressure by linearly interpolating integrated temperature
       integrated_temp(kb) = integrated_temp(kb) &
            + (new_pressure-pressure(k-1))/(pressure(k)-pressure(k-1)) * sol(S_TEMP,zlevels-k+2)%data(d)%elts(id_i)
    end do

    ! Temperature integrated over all vertical levels (does not change)
    do k = 1, zlevels
       integrated_temp(zlevels+1) = integrated_temp(zlevels+1) + sol(S_TEMP,k)%data(d)%elts(id_i)
    end do

    ! Scalars on new vertical grid
    do k = 1, zlevels
       sol(S_TEMP,k)%data(d)%elts(id_i) = integrated_temp(zlevels-k+2) - integrated_temp(zlevels-k+1)
       sol(S_MASS,k)%data(d)%elts(id_i) = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf)/grav_accel
    end do
  end subroutine remap_scalars_fv

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
