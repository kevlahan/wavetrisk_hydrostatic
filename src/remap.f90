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
  subroutine remap_vertical_coordinates()
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! interpolate on the full (not the perturbation) quantities
    ! Conserves mass, heat and momentum flux
    integer            :: d, j, k, l
    integer, parameter :: order_default = 7 ! order must be odd

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

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Remap variables on finest scale !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    ! Remap velocity first to use only values of mass on current grid to find old momentum
    do l = level_start, level_end
       ! Remap velocity first to use only values of mass on current grid to find old momentum
       call apply_onescale (remap_velocity, l, z_null, 0, 0)
       ! Remap mass and mass-weighted temperature
       call apply_onescale (remap_scalars,  l, z_null, 0, 1)
    end do
    call WT_after_step (sol, wav_coeff, level_start-1)

    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Remap at coarser scales !
    !!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! do l = level_end-1, level_start, -1
    !    do d = 1, size(grid)
    !       call cpt_or_restr_velo (grid(d), l)
    !    end do
    !    do d = 1, size(grid)
    !       do j = 1, grid(d)%lev(l)%length
    !          call apply_onescale_to_patch (remap_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
    !       end do
    !    end do

    !    do d = 1, size(grid)
    !       do j = 1, grid(d)%lev(l)%length
    !          call apply_onescale_to_patch (remap_scalars, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
    !       end do
    !    end do

    !    do d = 1, size(grid)
    !       do k = 1, zlevels
    !          mass => sol(S_MASS,k)%data(d)%elts
    !          temp => sol(S_TEMP,k)%data(d)%elts
    !          call cpt_or_restr_scalar (grid(d), l)
    !          nullify (mass, temp)
    !       end do
    !     end do
    ! end do

     ! Ensure boundary values are up to date
   ! call update_array_bdry (sol, NONE)
  end subroutine remap_vertical_coordinates

  subroutine cpt_or_restr_velo (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer :: c, j, p_par, p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) call apply_onescale_to_patch (remap_velocity, dom, p_par, z_null, -1, 0)
       end do
       call apply_interscale_to_patch (velo_cpt_restr, dom, dom%lev(l)%elts(j), z_null, -1, 0)
    end do
  end subroutine cpt_or_restr_velo

  subroutine velo_cpt_restr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    integer :: d, id_par, id_chd, idE_chd, idNE_chd, idN_chd, k
    real (8) :: u_prim_RT, u_prim_RT_E, u_prim_DG, u_prim_DG_NE, u_prim_UP, u_prim_UP_N

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)
    d = dom%id + 1

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd + RT + 1:EDGE*id_chd + UP + 1)) .lt. ADJZONE) then
       call remap_velocity (dom, i_par, j_par, z_null, offs_par, dims_par)
    end if

    do k = 1, zlevels
       ! Restriction is defined for edge integrated velocity (i.e. flux)
       u_prim_RT    = sol(S_VELO,k)%data(d)%elts(EDGE*id_chd+RT+1)*dom%len%elts(EDGE*id_chd+RT+1)
       u_prim_RT_E  = sol(S_VELO,k)%data(d)%elts(EDGE*idE_chd+RT+1)*dom%len%elts(EDGE*idE_chd+RT+1)
       u_prim_DG    = sol(S_VELO,k)%data(d)%elts(EDGE*id_chd+DG+1)*dom%len%elts(EDGE*id_chd+DG+1)
       u_prim_DG_NE = sol(S_VELO,k)%data(d)%elts(EDGE*idNE_chd+DG+1)*dom%len%elts(EDGE*idNE_chd+DG+1)
       u_prim_UP    = sol(S_VELO,k)%data(d)%elts(EDGE*id_chd+UP+1)*dom%len%elts(EDGE*id_chd+UP+1)
       u_prim_UP_N  = sol(S_VELO,k)%data(d)%elts(EDGE*idN_chd+UP+1)*dom%len%elts(EDGE*idN_chd+UP+1)

       if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) &
            sol(S_VELO,k)%data(d)%elts(EDGE*id_par+RT+1) = (u_prim_RT + u_prim_RT_E)/dom%len%elts(EDGE*id_par+RT+1)
       
       if (dom%mask_e%elts(DG+EDGE*id_chd+1) .ge. ADJZONE) &
            sol(S_VELO,k)%data(d)%elts(EDGE*id_par+DG+1) = (u_prim_DG + u_prim_DG_NE)/dom%len%elts(EDGE*id_par+DG+1)
       
       if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) &
            sol(S_VELO,k)%data(d)%elts(EDGE*id_par+UP+1) = (u_prim_UP + u_prim_UP_N)/dom%len%elts(EDGE*id_par+UP+1)
    end do
  end subroutine velo_cpt_restr

  subroutine cpt_or_restr_scalar (dom, l)
    type(Domain) :: dom
    integer      :: l

    integer :: j, p_par, c, p_chd
    logical :: restrict(N_CHDRN)

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .False.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .gt. 0) restrict(c) = .True.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) call apply_interscale_to_patch3 (scalar_cpt_restr, dom, p_par, c, z_null, 0, 1)
       end do
    end do
  end subroutine cpt_or_restr_scalar

  subroutine scalar_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute or restrict scalars
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par

    id_par = idx(i_par, j_par, offs_par, dims_par)
    
    if (dom%mask_n%elts(id_par+1) .ge. RESTRCT) then
       mass(id_par+1) = scalar_restr (mass)
       temp(id_par+1) = scalar_restr (temp)
    end if
  contains
    function scalar_restr (scalar)
      ! Mass conserving restriction
      real(8)                        :: scalar_restr
      real(8), dimension(:)          :: scalar
     
      integer :: id_chd, idE_chd, idNE_chd, idN_chd, idW_chd, idSW_chd, idS_chd

      id_chd   = idx(i_chd,   j_chd,   offs_chd, dims_chd)
      idE_chd  = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
      idNE_chd = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
      idN_chd  = idx(i_chd,   j_chd,   offs_chd, dims_chd)
      idW_chd  = idx(i_chd-1, j_chd,   offs_chd, dims_chd)
      idSW_chd = idx(i_chd-1, j_chd-1, offs_chd, dims_chd)
      idS_chd  = idx(i_chd,   j_chd-1, offs_chd, dims_chd)

      scalar_restr = (scalar(id_chd+1)/dom%areas%elts(id_chd+1)%hex_inv + 0.5_8 * ( &
           scalar(idE_chd+1) /dom%areas%elts(idE_chd+1)%hex_inv  + &
           scalar(idNE_chd+1)/dom%areas%elts(idNE_chd+1)%hex_inv + &
           scalar(idN_chd+1) /dom%areas%elts(idN_chd+1)%hex_inv  + &
           scalar(idW_chd+1) /dom%areas%elts(idW_chd+1)%hex_inv  + &
           scalar(idSW_chd+1)/dom%areas%elts(idSW_chd+1)%hex_inv + &
           scalar(idS_chd+1) /dom%areas%elts(idS_chd+1)%hex_inv)) * dom%areas%elts(id_par+1)%hex_inv
    end function scalar_restr
  end subroutine scalar_cpt_restr

  subroutine remap_velocity (dom, i, j, zlev, offs, dims)
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, idE, idN, idNE, k, kb, kc, kk, m
    real(8)                          :: layer_pressure, p_surf, p_surf_E, p_surf_N, p_surf_NE, diff, dmin, vel_old
    real(8)                          :: mass_id, mass_idE, mass_idN, mass_idNE
    real(8), dimension (3)           :: mass_e
    real(8), dimension (zlevels+1)   :: pressure
    real(8), dimension (zlevels+1,3) :: integrated_flux, new_flux

    d    = dom%id + 1
    id   = idx(i, j, offs, dims)
    id_i = id + 1

    idN  = idx(i,     j + 1, offs, dims) + 1
    idE  = idx(i + 1, j,     offs, dims) + 1
    idNE = idx(i + 1, j + 1, offs, dims) + 1

    ! Integrate full momentum flux vertically downward from the top
    ! All quantities located at interfaces
    integrated_flux(1,:) = 0.0_8
    do kb = 2, zlevels + 1
       k = zlevels-kb+2 ! Actual zlevel

       ! Interpolate mass on current vertical grid to edges 
       mass_e(RT+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idE))
       mass_e(DG+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idNE))
       mass_e(UP+1) = interp(sol(S_MASS,k)%data(d)%elts(id_i), sol(S_MASS,k)%data(d)%elts(idN))
       ! Mass fluxes
       do e = 1, EDGE
          integrated_flux(kb,e) = integrated_flux(kb-1,e) + mass_e(e) * sol(S_VELO,k)%data(d)%elts(EDGE*id+e)
       end do
    end do

    ! Calculate pressure at interfaces of current vertical grid, used as independent coordinate
    ! also calculate surface pressure at adjacent nodes needed to interpolate mass to edges
    pressure(1) = press_infty
    p_surf_E    = press_infty
    p_surf_NE   = press_infty
    p_surf_N    = press_infty
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel*sol(S_MASS,k)%data(d)%elts(id_i)
       p_surf_E  = p_surf_E  + grav_accel * sol(S_MASS,k)%data(d)%elts(idE) 
       p_surf_NE = p_surf_NE + grav_accel * sol(S_MASS,k)%data(d)%elts(idNE)
       p_surf_N  = p_surf_N  + grav_accel * sol(S_MASS,k)%data(d)%elts(idN) 
    end do
    p_surf = pressure(zlevels+1)

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    do kb = 1, zlevels+1
       k = zlevels-kb+2
       ! Pressure at top interface of cells of new vertical grid
       layer_pressure = a_vert(k)*ref_press + b_vert(k)*p_surf

       ! Find index of pressure on old vertical grid closest to layer_pressure on new grid
       dmin = 1d16
       do kk = 1, zlevels+1
          diff = abs(pressure(kk)-layer_pressure)
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do

       ! Set interpolation stencil based on layer_pressure
       if (kc .lt. (order-1)/2+1) then
          stencil = (/ (m, m = 1, order) /)
       else if (kc .gt. zlevels+1-(order-1)/2) then
          stencil = (/ (m, m = zlevels-order+2, zlevels+1) /)
       else
          stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
       end if

       ! Interpolate integrated mass flux at top interfaces of new vertical grid
       do e = 1, EDGE
          new_flux(kb,e) = Newton_interp(pressure(stencil), integrated_flux(stencil,e), layer_pressure)
       end do
    end do

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
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                          :: d, e, id, id_i, k, kb, kc, kk, m
    real(8)                          :: layer_pressure, p_surf, diff, dmin, vel_old
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
    ! also calculate surface pressure at adjacent nodes needed to interpolate mass to edges
    pressure(1) = press_infty
    do kb = 2, zlevels + 1
       k = zlevels-kb+2
       pressure(kb) = pressure(kb-1) + grav_accel * sol(S_MASS,k)%data(d)%elts(id_i)
    end do
    p_surf = pressure(zlevels+1)

    ! Interpolate using the moving stencil centred at each interpolation point computed downward from top
    do kb = 1, zlevels+1
       k = zlevels-kb+2
       ! Pressure at top interface of cells of new vertical grid
       layer_pressure = a_vert(k)*ref_press + b_vert(k)*p_surf

       ! Find index of pressure on old vertical grid closest to layer_pressure on new grid
       dmin = 1d16
       do kk = 1, zlevels+1
          diff = abs(pressure(kk)-layer_pressure)
          if (diff.lt.dmin) then
             kc = kk
             dmin = diff
          end if
       end do

       ! Set interpolation stencil based on layer_pressure
       if (kc .lt. (order-1)/2+1) then
          stencil = (/ (m, m = 1, order) /)
       else if (kc .gt. zlevels+1-(order-1)/2) then
          stencil = (/ (m, m = zlevels-order+2, zlevels+1) /)
       else
          stencil = (/ (m, m = kc-(order-1)/2, kc+(order-1)/2) /)
       end if

       ! Interpolate integrated temperature and integrated mass flux at top interfaces of new vertical grid
       new_temp(kb) = Newton_interp(pressure(stencil), integrated_temp(stencil), layer_pressure)
    end do

    ! Variables on new vertical grid
    do k = 1, zlevels
       kb = zlevels-k+1
       ! Remapped mass-weighted potential temperature from integrated value interpolated to new grid
       sol(S_TEMP,k)%data(d)%elts(id_i) = new_temp(kb+1) - new_temp(kb)
       sol(S_MASS,k)%data(d)%elts(id_i) = ((a_vert(k)-a_vert(k+1))*ref_press + (b_vert(k)-b_vert(k+1))*p_surf)/grav_accel
    end do
  end subroutine remap_scalars

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
