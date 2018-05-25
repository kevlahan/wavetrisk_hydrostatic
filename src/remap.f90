module remap_mod
  use adapt_mod
  use ops_mod
  use time_integr_mod
  use wavelet_mod
  implicit none
contains
  subroutine remap_vertical_coordinates (set_thresholds)
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! interpolate on the full (not the perturbation) quantities
    ! Conserves mass, heat and momentum flux

    external :: set_thresholds
    
    integer :: d, j, k, l

    if (rank == 0) write(6,*) "Remapping vertical coordinates"

    ! Ensure boundary values are up to date
    call update_array_bdry (sol, NONE)
    
    ! Remap on finest level
    call apply_onescale (remap_scalars, level_end, z_null, 0, 1)
    call apply_onescale (remap_velo,    level_end, z_null, 0, 0)

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
       call apply_onescale (remap_velo,    l, z_null, 0, 0)

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
        
    ! Re-adapt grid after remapping
    call WT_after_step (sol, wav_coeff, level_start-1)
    if (adapt_trend) then ! Find trend wavelets on new vertical grid
       call trend_ml (sol, trend)
       call WT_after_step (trend, trend_wav_coeff, level_start-1)
    end if
    call adapt_grid (set_thresholds)
  end subroutine remap_vertical_coordinates

  subroutine remap_scalars (dom, i, j, z_null, offs, dims)
    ! Remap mass and vertical variable eta onto grid defined by A, B coefficients
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: current_zlev, d, id, k, zlev
    real(8)                          :: column_mass, cumul_mass_zlev, cumul_mass_target, new_cumul_mass, cumul_mass_upper, X
    real(8)                          :: new_mass
    real(8), dimension (zlevels+1)   :: cumul_mass, cumul_temp, new_cumul_temp

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1)<TRSK) return
    
    ! Calculate cumulative mass and total mass of column
    cumul_mass(1) = 0.0_8
    cumul_temp(1) = 0.0_8
    do k = 1, zlevels
       cumul_mass(k+1) = cumul_mass(k) + sol(S_MASS,k)%data(d)%elts(id+1)
       cumul_temp(k+1) = cumul_temp(k) + sol(S_TEMP,k)%data(d)%elts(id+1)
    end do
    column_mass = cumul_mass(zlevels+1)

    current_zlev = 1
    exner_fun(1)%data(d)%elts(id+1) = 1.0_8
    new_cumul_mass = 0.0_8
    do k = 1, zlevels
       trend(S_MASS,k)%data(d)%elts(id+1) = sol(S_MASS,k)%data(d)%elts(id+1)
       
       sol(S_MASS,k)%data(d)%elts(id+1) = a_vert_mass(k) + b_vert_mass(k) * column_mass
       
       cumul_mass_target = new_cumul_mass + sol(S_MASS,k)%data(d)%elts(id+1)

       do zlev = current_zlev, zlevels
          cumul_mass_upper = cumul_mass(zlev+1)
          if (cumul_mass_target <= cumul_mass_upper) exit
       end do

       if (zlev > zlevels) zlev = zlevels
       cumul_mass_zlev = cumul_mass(zlev)

       current_zlev = zlev
       new_cumul_mass = cumul_mass_target

       ! New vertical coordinate (saved in exner_fun)
       exner_fun(k+1)%data(d)%elts(id+1) = zlev + (cumul_mass_target - cumul_mass_zlev)/(cumul_mass_upper - cumul_mass_zlev)
    end do

    ! Remap theta
    do k = 1, zlevels+1
       X = exner_fun(k)%data(d)%elts(id+1)
       zlev = min(zlevels,floor(X)) ; if (zlev<1) return
       X = X - zlev
       new_cumul_temp(k) = cumul_temp(zlev) + X*sol(S_TEMP,zlev)%data(d)%elts(id+1)
    end do

    do k = 1, zlevels
       sol(S_TEMP,k)%data(d)%elts(id+1) = new_cumul_temp(k+1) - new_cumul_temp(k)
    end do
  end subroutine remap_scalars

  subroutine remap_velo (dom, i, j, z_null, offs, dims)
    ! Remap velocity onto original vertical grid by linear interpolation of mass flux
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                              :: d, e, id, idE, idN, idNE, k, zlev
    real(8), dimension(EDGE)             :: mass_e, X
    real(8), dimension(zlevels+1,1:EDGE) :: massflux_cumul, massflux, new_massflux_cumul

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)

    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    
    ! Do not remap inactive nodes
    if (dom%mask_n%elts(id+1)<TRSK  .or. &
         dom%mask_e%elts(EDGE*id+RT+1)<TRSK .or. dom%mask_e%elts(EDGE*id+DG+1)<TRSK .or. dom%mask_e%elts(EDGE*id+UP+1)<TRSK) return 

    massflux_cumul(1,:) = 0.0_8
    do k = 1, zlevels
       ! Interpolate old masses (stored in trend)
       mass_e(RT+1) = interp(trend(S_MASS,k)%data(d)%elts(id+1), trend(S_MASS,k)%data(d)%elts(idE+1))
       mass_e(DG+1) = interp(trend(S_MASS,k)%data(d)%elts(id+1), trend(S_MASS,k)%data(d)%elts(idNE+1))
       mass_e(UP+1) = interp(trend(S_MASS,k)%data(d)%elts(id+1), trend(S_MASS,k)%data(d)%elts(idN+1))

       do e = 1, EDGE
          massflux(k,e) = sol(S_VELO,k)%data(d)%elts(EDGE*id+e) * mass_e(e)
          massflux_cumul(k+1,e) = massflux_cumul(k,e) + massflux(k,e)
       end do
    end do

    do k = 1, zlevels+1
       X(RT+1) = interp(exner_fun(k)%data(d)%elts(id+1), exner_fun(k)%data(d)%elts(idE+1))
       X(DG+1) = interp(exner_fun(k)%data(d)%elts(id+1), exner_fun(k)%data(d)%elts(idNE+1))
       X(UP+1) = interp(exner_fun(k)%data(d)%elts(id+1), exner_fun(k)%data(d)%elts(idN+1))

       do e = 1, EDGE
          zlev = min(zlevels,floor(X(e))) ; if (zlev<1) return
          X(e) = X(e) - zlev
          new_massflux_cumul(k,e) = massflux_cumul(zlev,e) + X(e)*massflux(zlev,e)
       end do
    end do

    do k = 1, zlevels
       ! Interpolate new masses
       mass_e(RT+1) = interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k)%data(d)%elts(idE+1))
       mass_e(DG+1) = interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k)%data(d)%elts(idNE+1))
       mass_e(UP+1) = interp(sol(S_MASS,k)%data(d)%elts(id+1), sol(S_MASS,k)%data(d)%elts(idN+1))

       do e = 1, EDGE
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = (new_massflux_cumul(k+1,e) - new_massflux_cumul(k,e)) / mass_e(e)
       end do
    end do
  end subroutine remap_velo
end module remap_mod
