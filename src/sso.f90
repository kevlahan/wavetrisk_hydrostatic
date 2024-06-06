module sso_mod
  ! SSO (subgrid scale orography) parameterization
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use utils_mod
  implicit none
  logical :: circle        = .true. ! assume circular mountain (not enough SSO statistics)
  logical :: blocking_drag = .true. 
  logical :: wave_drag     = .false.
contains
  function sso_drag (dom, i, j, z_null, offs, dims)
    ! SSO block and wave drag at edges
    ! uses version from Japanese Meteorological Agency (2019) report
    use coord_arithmetic_mod
    implicit none
    type(Domain)                         :: dom
    integer                              :: i, j, z_null
    integer, dimension(N_BDRY+1)         :: offs
    integer, dimension(2,N_BDRY+1)       :: dims
    real(8), dimension(1:zlevels,1:EDGE) :: sso_drag

    integer                           :: d, id, idE, idNE, idN, k, nlev
    integer, dimension(1:EDGE)        :: id_e

    real(8)                           :: mu, gamma, sigma, p, theta, N1, N2
    real(8)                           :: B, C, H, H_eff, H_env, H_peak, r, Z_block
    real(8)                           :: N_above, N_below, N_av, phi_av, psi_av, rho_av, rho_dz
    real(8)                           :: lat, lon
    real(8), dimension(1:EDGE)        :: u_av
    real(8), dimension(2)             :: vel
    real(8), dimension(1:zlevels)     :: N_bv, phi, psi, rho, drag_block, drag_wave, umag
    real(8), dimension(0:zlevels)     :: z

    type(Coord)                       :: e_zonal, e_merid, e_Uav, e_U, e_V, e_W

    real(8), parameter                :: C_d    = 2.00d0
    real(8), parameter                :: H_crit = 0.50d0
    real(8), parameter                :: G      = 0.25d0
    real(8), parameter                :: n_eff  = 2.4d0

    sso_drag   = 0d0
    drag_block = 0d0
    drag_wave  = 0d0

    d    = dom%id + 1
    
    id   = idx (i, j, offs, dims)
    id_e = id_edge (id)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    ! SSO parameters
    mu    =    sso_param(S_MU)%data(d)%elts(id+1)    
    theta = sso_param(S_THETA)%data(d)%elts(id+1) 
    gamma = sso_param(S_GAMMA)%data(d)%elts(id+1) 
    sigma = sso_param(S_SIGMA)%data(d)%elts(id+1) 

    H_peak = 3d0 * mu ! maximum SSO peak height
    H_env  = 2d0 * mu ! SSO envelope

    B = 1d0 - 0.18d0 * gamma - 0.04d0 * gamma**2
    C =       0.48d0 * gamma + 0.30d0 * gamma**2

    ! Unit vectors
    call cart2sph (dom%node%elts(id+1), lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),               0d0) 
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))

    e_U = direction (dom%node%elts(id+1),   dom%node%elts(idE+1))
    e_V = direction (dom%node%elts(idNE+1), dom%node%elts(id+1))
    e_W = direction (dom%node%elts(id+1),   dom%node%elts(idN+1))

    ! Brunt-Vaisala frequency at layer centres
    N_bv(1) = N_i (dom, i, j, 1, offs, dims); N_below = N_bv(1)
    do k = 2, zlevels-1
       N_above = N_i (dom, i, j, k, offs, dims)
       N_bv(k) = interp (N_below, N_above)
       N_below = N_above
    end do
    N_bv(zlevels) = N_below   

    ! Save layer values
    z(0) = 0d0
    do k = 1, zlevels
       rho(k)  = density_i (dom, i, j, k, offs, dims, sol)     ! density
       vel     = uvw2zonal_merid (dom, i, j, k, offs, dims)    ! zonal and meridional velocities at node
       umag(k) = sqrt (sum (vel**2))                           ! velocity magnitude
       z(k)    = z(k-1) + dz_i (dom, i, j, k, offs, dims, sol) ! height of upper interface above topography

       if (.not. circle) then
          phi(k) = atan2 (vel(2), vel(1))                      ! angle of incident flow
          psi(k) = theta - phi(k)                              ! angle of incident flow with respect to principal axis of ellipse
       end if
    end do

    ! Find blocking height
    if (blocking_drag) then
       H = 0d0
       Z_block = 0d0
       do k = zlevels, 1, -1
          if (z(k) <= H_peak) then
             H = H + N_bv(k) * dz_i (dom, i, j, k, offs, dims, sol) / umag(k) ! non-dimensional height
             if (H >= H_crit) then
                Z_block = z(k)
                exit
             end if
          end if
       end do
    end if

    ! Compute mean values vertical layers mu <= z - z_s <= 2 mu
    N_av    = 0d0
    u_av    = 0d0
    rho_av  = 0d0
    psi_av  = 0d0
    nlev = 0
    do k = 1, zlevels
       if (z(k-1) > mu .and. z(k) < H_env) then ! average over layers *between* mu and 2 mu
          nlev = nlev + 1
          N_av   = N_av   + N_bv(k)    
          rho_av = rho_av + rho(k)  
          u_av   = u_av   + sol(S_VELO,k)%data(d)%elts(id_e)
          if (.not. circle) phi_av = phi_av + phi(k)
       elseif (z(k) >= H_env) then
          exit
       end if
    end do

    if (nlev > 0) then
       N_av   = N_av   / dble (nlev)
       phi_av = phi_av / dble (nlev)
       rho_av = rho_av / dble (nlev)
       u_av   = u_av   / dble (nlev)
    end if
    if (.not. circle) then
       psi_av = theta - phi_av
       e_Uav = cos(phi_av) * e_zonal + sin(phi_av) * e_merid
    end if

    ! Compute gravity wave drag drag (non-zero in lowest layer only)
    if (wave_drag) then
       if (blocking_drag) then
          H_eff = n_eff * (H_peak - Z_block)
       else
          H_eff = H_peak
       end if
       rho_dz = sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1)

       if (circle) then ! circular mountain
          drag_wave(1) = - rho_av * N_av * (H_eff/3d0)**2 * sigma/mu * G  * 0.78d0 / rho_dz
       else
          drag_wave(1) = - rho_av * N_av * (H_eff/3d0)**2 * sigma/mu * G  &
               * ((B*cos(psi_av)**2 + C*sin(psi_av)**2) + (B-C)*sin(psi_av)*cos(psi_av)) / rho_dz
       end if
    end if

    ! Compute blocking drag magnitude
    if (blocking_drag) then
       do k = 1, zlevels
          if (z(k) <= Z_block) then
             if (circle) then ! circular mountain
                drag_block(k) = - 0.5d0 * C_d * 4d0  * sigma/H_env * sqrt ((Z_block - z(k))/(z(k) + mu)) * umag(k) * 0.78d0  
             else
                r = sqrt ( (cos(psi(k))**2 + gamma**2 * sin(psi(k))**2) / (gamma**2 * cos(psi(k))**2 + sin(psi(k))**2) )

                drag_block(k) = - 0.5d0 * C_d * max (5d0 - 1d0/r**3, 0d0)  * sigma/H_env &
                     * sqrt ( (Z_block - z(k))/(z(k) + mu) ) * umag(k) * (B * cos(psi(k))**2  + C * sin(psi(k))**2)
             end if
          else 
             exit
          end if
       end do
    end if

    ! Complete drag
    do k = 1, zlevels
       if (blocking_drag) sso_drag(k,:) = drag_block(k) * sol(S_VELO,k)%data(d)%elts(id_e)

       if (wave_drag) then
          if (circle) then
             sso_drag(k,:) = sso_drag(k,:) + drag_wave(k) * u_av
          else
             sso_drag(k,RT+1) = sso_drag(k,RT+1) + inner (drag_wave(k)*e_Uav, e_U) 
             sso_drag(k,DG+1) = sso_drag(k,DG+1) + inner (drag_wave(k)*e_Uav, e_V) 
             sso_drag(k,UP+1) = sso_drag(k,UP+1) + inner (drag_wave(k)*e_Uav, e_W) 
          end if
       end if
    end do
  end function sso_drag

  subroutine cal_sso_param (dom, i, j, zlev, offs, dims)
    ! mu (standard deviation) of Subgrid Scale Orography (SSO) compared to Grid Scale Orography (GSO)
    ! Use area-weighted integral including approximate overlapping coarse-fine hexagonal cells for levels < max_level-1
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, ii, jj, n_topo
    real(8) :: distance, dx, h, hx, hy, h_sq, hx_sq, hy_sq, hxhy, total_area
    real(8) :: K, L, M
    real(8) :: gamma, mu, sigma, theta, topo_Area_min

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    ! Include all finest grid topography in a disk of radius dx
    dx  = dom%len%elts(EDGE*id+RT+1) / sqrt(3d0)                       ! current cell radius

    topo_Area_min = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**topo_max_level) ! topography cell area

    h_sq       = 0d0
    hx_sq      = 0d0
    hy_sq      = 0d0
    hxhy       = 0d0

    total_area = 0d0
    n_topo     = size (topography_data(topo_max_level,d)%node)

    ! Integrate finest resolution topography over current cell size
    do ii = 1, n_topo
       distance = dist (dom%node%elts(id+1), topography_data(topo_max_level,d)%node(ii))
       if (distance <= dx) then
          jj = 3*(ii-1) + 1

          h  = topography_data(topo_max_level,d)%elts(jj) - topography%data(d)%elts(id+1)

          hx = topography_data(topo_max_level,d)%elts(jj+1)  ! topography gradient in longitude direction
          hy = topography_data(topo_max_level,d)%elts(jj+2)  ! topography gradient in latitude direction

          h_sq  = h_sq  + h**2
          hx_sq = hx_sq + hx**2
          hy_sq = hy_sq + hy**2
          hxhy  = hxhy  + hx * hy

          total_area = total_area + topo_Area_min
       end if
    end do

    hx_sq = hx_sq * topo_Area_min / total_area
    hy_sq = hy_sq * topo_Area_min / total_area
    M     = hxhy  * topo_Area_min / total_area

    K = 0.5d0 * (hx_sq + hy_sq)
    L = 0.5d0 * (hx_sq - hy_sq) 

    ! SSO parameters
    sso_param(S_MU)%data(d)%elts(id+1) = sqrt (h_sq * topo_Area_min / total_area)
    sso_param(S_THETA)%data(d)%elts(id+1) = 0.5d0 * atan2 (M, L)
    sso_param(S_GAMMA)%data(d)%elts(id+1) = sqrt ( (K - sqrt (L**2 + M**2)) / (K + sqrt (L**2 + M**2)) )
    if (circle) then
       sso_param(S_SIGMA)%data(d)%elts(id+1) = sqrt (hx_sq + hy_sq)
    else
       sso_param(S_SIGMA)%data(d)%elts(id+1) = sqrt (hx_sq * cos (theta) + hy_sq * sin (theta))
    end if
  end subroutine cal_sso_param
end module sso_mod
