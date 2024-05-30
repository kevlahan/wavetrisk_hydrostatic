module sso_mod
  ! SSO (subgrid scale orography) parameterization
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use utils_mod
  use io_mod
  implicit none
contains
  function sso_wave_drag (dom, i, j, zlev, offs, dims)
    ! Gravity wave drag at edges
    ! uses version from Japanese Meteorological Agency (2019) report
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: sso_wave_drag

    integer               :: d, id, idE, idNE, idN, k, nlev

    real(8)               :: mu, gamma, sigma, p, theta
    real(8)               :: B, C, H, H_eff, H_mount, N, phi, psi, rho, tau_mag, z, U, Z_block
    real(8)               :: lat, lon
    real(8), dimension(2) :: vel

    type(Coord)           :: e_zonal, e_merid, e_U, e_V, e_W, tau

    real(8), parameter    :: H_crit   = 0.5d0
    real(8), parameter    :: G        = 0.25d0
    logical, parameter    :: blocking = .false. ! reduce wave drag by including blocking effect

    id = idx (i, j, offs, dims)

    if (dom%level%elts(id+1) < max_level .and. zlev == 1) then ! apply surface stress only in bottom layer
       d = dom%id + 1

       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)

       call cart2sph (dom%node%elts(id+1), lon, lat)

       e_zonal = Coord (-sin(lon),           cos(lon),               0d0) 
       e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))

       e_U = direction (dom%node%elts(id+1),   dom%node%elts(idE+1))
       e_V = direction (dom%node%elts(idNE+1), dom%node%elts(id+1))
       e_W = direction (dom%node%elts(id+1),   dom%node%elts(idN+1))

       ! SSO parameters
       mu    = sso_param(S_MU)%data(d)%elts(id+1)    
       theta = sso_param(S_THETA)%data(d)%elts(id+1) 
       gamma = sso_param(S_GAMMA)%data(d)%elts(id+1) 
       sigma = sso_param(S_SIGMA)%data(d)%elts(id+1) 

       H_mount = 3d0 * mu ! assume entire flow passes over mountain

       B = 1d0 - 0.18d0 * gamma -0.04d0 * gamma**2
       C = 0.48d0 * gamma + 0.30d0 * gamma*2

       ! Find blocking height
       H = 0d0
       Z_block = 0d0
       do k = zlevels, 1, -1
          z = zl_i (dom, i, j, k, offs, dims, sol, 1) - topography%data(d)%elts(id+1)
          if (z <= H_mount) then
             ! Non-dimensional height
             H = H + N_i (dom, i, j, k, offs, dims) * dz_i (dom, i, j, zlev, offs, dims, sol) / u_mag (dom, i, j, k, offs, dims)
             if (H >= H_crit) then
                Z_block = z
                exit
             end if
          end if
       end do

       ! Compute mean values vertical layers mu <= z - z_s <= 2 mu
       N    = 0d0
       U    = 0d0
       rho  = 0d0
       phi  = 0d0
       nlev = 0
       do k = 1, zlevels
          z =  zl_i (dom, i, j, k, offs, dims, sol, 1) - topography%data(d)%elts(id+1) ! height of upper interface above topography
          if (z >= mu .and. z <= 2d0*mu) then
             nlev = nlev + 1
             N    = N + N_i (dom, i, j, k, offs, dims)              ! Brunt-Vaisala frequency
             U    = U + u_mag (dom, i, j, k, offs, dims)            ! velocity magnitude
             rho  = rho + density_i (dom, i, j, k, offs, dims, sol) ! density

             vel = uvw2zonal_merid (dom, i, j, zlev, offs, dims)    ! zonal and meridional velocities at node
             phi = phi + atan2 (vel(2), vel(1))                     ! incident flow angle
          elseif (z > 2d0*mu) then
             exit
          end if
       end do
       if (nlev > 0) then
          N   = N   / dble (nlev)
          U   = U   / dble (nlev)
          rho = rho / dble (nlev)
          phi = phi / dble (nlev)
       end if

       ! Effective mountain height 
       if (blocking) then
          H_eff = H_mount - Z_block
       else
          H_eff = H_mount
       end if

       psi = theta - psi ! angle between incident flow and principal axis of elliptical mountain

       tau_mag = rho *  U * N * (H_eff/3d0)**2 * sigma / mu * G

       tau = tau_mag * ( (cos(psi)**2 + C*sin(psi)**2) * e_zonal + (B-C)*sin(psi)*cos(psi) * e_merid )

       sso_wave_drag(RT+1) = inner (tau, e_U)
       sso_wave_drag(DG+1) = inner (tau, e_V)
       sso_wave_drag(UP+1) = inner (tau, e_W)
    else
       sso_wave_drag = 0d0
    end if
  end function sso_wave_drag

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
    id = idx (i, j, offs, dims) + 1

    dx  = dom%len%elts(EDGE*id+RT+1) * 2d0/sqrt(3d0)                       ! current cell radius

    topo_Area_min = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**topo_max_level) ! topography cell area
    ! Include all finest grid topography in a disk of radius dx
    h_sq       = 0d0
    hx_sq      = 0d0
    hy_sq      = 0d0
    hxhy       = 0d0

    total_area = 0d0
    n_topo     = size (topography_data(topo_max_level,d)%node)

    ! Integrate finest resolution topography over current cell size
    do ii = 1, n_topo
       distance = dist (dom%node%elts(id), topography_data(topo_max_level,d)%node(ii))
       if (distance <= dx) then
          jj = 3*(ii-1) + 1

          h  = topography_data(topo_max_level,d)%elts(jj) - topography%data(d)%elts(id)

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
    sso_param(S_MU)%data(d)%elts(id)    = sqrt (h_sq * topo_Area_min / total_area)
    sso_param(S_THETA)%data(d)%elts(id) = 0.5d0 * atan2 (M, L)
    sso_param(S_GAMMA)%data(d)%elts(id) = sqrt ( (K - sqrt (L**2 + M**2)) / (K + sqrt (L**2 + M**2)) )
    sso_param(S_SIGMA)%data(d)%elts(id) = sqrt (hx_sq * cos (theta) + hy_sq * sin (theta))
  end subroutine cal_sso_param
end module sso_mod
