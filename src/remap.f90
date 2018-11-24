module remap_mod
  use time_integr_mod
  use wavelet_mod
  implicit none
contains
  subroutine remap_vertical_coordinates
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! Conserves mass, heat and momentum flux
    integer :: l

    ! Save old mass in trend
    trend = sol

    ! Remap
    do l = level_start, level_end
       call apply_onescale (remap_scalars, l, z_null, 0, 1)
       call apply_onescale (remap_velo_RT, l, z_null, 0, 0)
       call apply_onescale (remap_velo_DG, l, z_null, 0, 0)
       call apply_onescale (remap_velo_UP, l, z_null, 0, 0)
    end do
    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE)
        
    ! Wavelet transform and interpolate back onto adapted grid
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine remap_vertical_coordinates

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

    if (dom%mask_n%elts(id_i) < TRSK) return
    
    ! Calculate cumulative mass and total mass of column
    cumul_mass(1) = 0.0_8
    cumul_temp(1) = 0.0_8
    do k = 1, zlevels
       cumul_mass(k+1) = cumul_mass(k) + sol(S_MASS,k)%data(d)%elts(id_i)
       cumul_temp(k+1) = cumul_temp(k) + sol(S_TEMP,k)%data(d)%elts(id_i)
    end do
    column_mass = cumul_mass(zlevels+1)

    current_zlev = 1
    exner_fun(1)%data(d)%elts(id_i) = 1.0_8
    new_cumul_mass = 0.0_8
    do k = 1, zlevels
       sol(S_MASS,k)%data(d)%elts(id_i) = a_vert_mass(k) + b_vert_mass(k) * column_mass
       
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
       zlev = min (zlevels, floor (X)) ; if (zlev < 1) return
       X = X - zlev
       new_cumul_temp(k) = cumul_temp(zlev) + X * sol(S_TEMP,zlev)%data(d)%elts(id_i)
    end do
    
    do k = 1, zlevels
       sol(S_TEMP,k)%data(d)%elts(id_i) = new_cumul_temp(k+1) - new_cumul_temp(k)
    end do
  end subroutine remap_scalars

  subroutine remap_velo_RT (dom, i, j, z_null, offs, dims)
    ! Remap RT velocity onto original vertical grid by linear interpolation of mass flux
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, idE_i, k, zlev
    real(8)                       :: mass_e, X_e
    real(8), dimension(zlevels)   :: massflux
    real(8), dimension(zlevels+1) :: massflux_cumul, new_massflux_cumul

    d     = dom%id + 1
    id    = idx (i, j, offs, dims)
    id_i  = id + 1
    id_e  = EDGE*id + RT + 1
    idE_i = idx(i+1, j,   offs, dims) + 1
    
    ! Do not remap inactive nodes
    if (dom%mask_n%elts(id_i) < TRSK .or. dom%mask_e%elts(id_e) < TRSK) return 

    massflux_cumul(1) = 0.0_8
    do k = 1, zlevels
       mass_e = interp (trend(S_MASS,k)%data(d)%elts(id_i), trend(S_MASS,k)%data(d)%elts(idE_i))
       massflux(k) = sol(S_VELO,k)%data(d)%elts(id_e) * mass_e 
       massflux_cumul(k+1) = massflux_cumul(k) + massflux(k)
    end do
      
    do k = 1, zlevels+1
       X_e = interp (exner_fun(k)%data(d)%elts(id_i), exner_fun(k)%data(d)%elts(idE_i))
       do e = 1, EDGE
          zlev = min (zlevels, floor (X_e)) ; if (zlev < 1) return
          X_e = X_e - zlev
       end do
       new_massflux_cumul(k) = massflux_cumul(zlev) + X_e * massflux(zlev)
    end do
    
    do k = 1, zlevels
       mass_e = interp (sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idE_i))
       sol(S_VELO,k)%data(d)%elts(id_e) = (new_massflux_cumul(k+1) - new_massflux_cumul(k)) / mass_e
    end do
  end subroutine remap_velo_RT

   subroutine remap_velo_DG (dom, i, j, z_null, offs, dims)
    ! Remap velocity onto original vertical grid by linear interpolation of mass flux
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, idNE_i, k, zlev
    real(8)                       :: mass_e, X_e
    real(8), dimension(zlevels)   :: massflux
    real(8), dimension(zlevels+1) :: massflux_cumul, new_massflux_cumul

    d      = dom%id + 1
    id     = idx (i, j, offs, dims)
    id_i   = id + 1
    id_e   = EDGE*id + DG + 1
    idNE_i = idx (i+1, j+1, offs, dims) + 1
    
    ! Do not remap inactive nodes
    if (dom%mask_n%elts(id_i) < TRSK .or. dom%mask_e%elts(id_e) < TRSK) return 

    massflux_cumul(1) = 0.0_8
    do k = 1, zlevels
       mass_e = interp (trend(S_MASS,k)%data(d)%elts(id_i), trend(S_MASS,k)%data(d)%elts(idNE_i))
       massflux(k) = sol(S_VELO,k)%data(d)%elts(id_e) * mass_e 
       massflux_cumul(k+1) = massflux_cumul(k) + massflux(k)
    end do
      
    do k = 1, zlevels+1
       X_e = interp (exner_fun(k)%data(d)%elts(id_i), exner_fun(k)%data(d)%elts(idNE_i))
       do e = 1, EDGE
          zlev = min (zlevels, floor (X_e)) ; if (zlev < 1) return
          X_e = X_e - zlev
       end do
       new_massflux_cumul(k) = massflux_cumul(zlev) + X_e * massflux(zlev)
    end do
    
    do k = 1, zlevels
       mass_e = interp (sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idNE_i))
       sol(S_VELO,k)%data(d)%elts(id_e) = (new_massflux_cumul(k+1) - new_massflux_cumul(k)) / mass_e
    end do
  end subroutine remap_velo_DG

   subroutine remap_velo_UP (dom, i, j, z_null, offs, dims)
    ! Remap UP velocity onto original vertical grid by linear interpolation of mass flux
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, idN_i, k, zlev
    real(8)                       :: mass_e, X_e
    real(8), dimension(zlevels)   :: massflux
    real(8), dimension(zlevels+1) :: massflux_cumul, new_massflux_cumul

    d     = dom%id + 1
    id    = idx (i, j, offs, dims)
    id_i  = id + 1
    id_e  = EDGE*id + UP + 1
    idN_i = idx (i, j+1, offs, dims) + 1
    
    ! Do not remap inactive nodes
    if (dom%mask_n%elts(id_i) < TRSK .or. dom%mask_e%elts(id_e) < TRSK) return 

    massflux_cumul(1) = 0.0_8
    do k = 1, zlevels
       mass_e = interp (trend(S_MASS,k)%data(d)%elts(id_i), trend(S_MASS,k)%data(d)%elts(idN_i))
       massflux(k) = sol(S_VELO,k)%data(d)%elts(id_e) * mass_e 
       massflux_cumul(k+1) = massflux_cumul(k) + massflux(k)
    end do
      
    do k = 1, zlevels+1
       X_e = interp (exner_fun(k)%data(d)%elts(id_i), exner_fun(k)%data(d)%elts(idN_i))
       do e = 1, EDGE
          zlev = min (zlevels, floor (X_e)) ; if (zlev < 1) return
          X_e = X_e - zlev
       end do
       new_massflux_cumul(k) = massflux_cumul(zlev) + X_e * massflux(zlev)
    end do
    
    do k = 1, zlevels
       mass_e = interp (sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idN_i))
       sol(S_VELO,k)%data(d)%elts(id_e) = (new_massflux_cumul(k+1) - new_massflux_cumul(k)) / mass_e
    end do
  end subroutine remap_velo_UP
 end module remap_mod
