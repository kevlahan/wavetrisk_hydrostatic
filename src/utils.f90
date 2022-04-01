module utils_mod
  ! Basic functions
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  implicit none
contains
  real(8) function z_i (dom, i, j, zlev, offs, dims)
    ! Position of vertical level zlev at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, k
    real(8) :: dz, dz_below, z_s

    id = idx (i, j, offs, dims)

    z_s = dom%topo%elts(id+1) 
    
    dz_below = dz_i (dom, i, j, 1, offs, dims)
    z_i = z_s + dz_below / 2
    do k = 2, zlev
       dz = dz_i (dom, i, j, k, offs, dims)
       z_i = z_i + interp (dz, dz_below)
       dz_below = dz
    end do
  end function z_i

  real(8) function dz_i (dom, i, j, zlev, offs, dims)
    ! Thickness of layer zlev at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    dz_i = (sol_mean(S_MASS,zlev)%data(d)%elts(id) + sol(S_MASS,zlev)%data(d)%elts(id)) / porous_density (d, id, zlev)
  end function dz_i
  
  function dz_e (dom, i, j, zlev, offs, dims)
    ! Thickness of layer zlev at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: dz_e

    integer                    :: d, id, idE, idN, idNE
    real(8), dimension(0:EDGE) :: dz

    d = dom%id + 1

    id   = idx (i, j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
   
    dz(0) = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1)   + sol(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (d, id+1, zlev)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)  + sol(S_MASS,zlev)%data(d)%elts(idE+1)) &
         / porous_density (d, idE+1, zlev)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1) + sol(S_MASS,zlev)%data(d)%elts(idNE+1)) &
         / porous_density (d, idNE+1, zlev)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)  + sol(S_MASS,zlev)%data(d)%elts(idN+1)) &
         / porous_density (d, idN+1, zlev)

    dz_e = 0.5d0 * (dz(0) + dz(1:EDGE))
  end function dz_e

  function dz_SW_e (dom, i, j, zlev, offs, dims)
    ! Thickness of layer zlev at edges (SW edges)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: dz_SW_e

    integer                    :: d, id, idW, idSW, idS
    real(8), dimension(0:EDGE) :: dz

    d = dom%id + 1

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
   
    dz(0) = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1)   + sol(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (d, id+1, zlev)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idW+1)  + sol(S_MASS,zlev)%data(d)%elts(idW+1)) &
         / porous_density (d, idW+1, zlev)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idSW+1) + sol(S_MASS,zlev)%data(d)%elts(idSW+1)) &
         / porous_density (d, idSW+1, zlev)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idS+1)  + sol(S_MASS,zlev)%data(d)%elts(idS+1)) &
         / porous_density (d, idS+1, zlev)

    dz_SW_e = 0.5d0 * (dz(0) + dz(1:EDGE))
  end function dz_SW_e
  
  real(8) function zl_i (dom, i, j, zlev, offs, dims, l)
    ! Position of interface below (l=-1) or above (l=1) vertical level zlev nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, k

    id = idx (i, j, offs, dims)

    zl_i = dom%topo%elts(id+1)
    do k = 1, zlev+l
       zl_i = zl_i + dz_i (dom, i, j, k, offs, dims)
    end do
  end function zl_i

  function zl_e (dom, i, j, zlev, offs, dims, l)
    ! Position of interface below (l=-1) or above (l=1) vertical level zlev at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: zl_e

    integer :: id, idE, idN, idNE, k

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    zl_e(RT+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idE+1))  
    zl_e(DG+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idNE+1)) 
    zl_e(UP+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idN+1))  
    do k = 1, zlev+l
       zl_e = zl_e + dz_e (dom, i, j, k, offs, dims)
    end do
  end function zl_e

  function eta_e (dom, i, j, zlev, offs, dims)
    ! Free surface perturbation at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: eta_e

    integer :: d, id, idE, idN, idNE
    real(8) :: eta0

    d = dom%id + 1
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (mode_split) then
       eta_e(RT+1) = interp (sol(S_MASS,zlevels+1)%data(d)%elts(id+1), sol(S_MASS,zlevels+1)%data(d)%elts(idE+1))
       eta_e(DG+1) = interp (sol(S_MASS,zlevels+1)%data(d)%elts(id+1), sol(S_MASS,zlevels+1)%data(d)%elts(idNE+1))
       eta_e(UP+1) = interp (sol(S_MASS,zlevels+1)%data(d)%elts(id+1), sol(S_MASS,zlevels+1)%data(d)%elts(idN+1))
    else
       eta0 = free_surface (dom, i, j, zlev, offs, dims)
       eta_e(RT+1) = interp (eta0, free_surface (dom, i+1, j,   zlev, offs, dims))
       eta_e(DG+1) = interp (eta0, free_surface (dom, i+1, j+1, zlev, offs, dims))
       eta_e(UP+1) = interp (eta0, free_surface (dom, i,   j+1, zlev, offs, dims))
    end if
  end function eta_e

  real(8) function porous_density (d, id_i, zlev)
    ! Porous density at nodes
    implicit none
    integer :: d, id_i, zlev
    
    porous_density = ref_density * (1d0 + (alpha - 1d0) * penal_node(zlev)%data(d)%elts(id_i))
  end function porous_density

  function porous_density_edge (d, id, zlev)
    ! Porous density at edges
    implicit none
    integer                        :: d, id, zlev
    real(8), dimension(1:EDGE)     :: porous_density_edge
    
    porous_density_edge = ref_density * (1d0 + (alpha - 1d0) * penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1))
  end function porous_density_edge

  real(8) function phi_node (d, id_i, zlev)
    ! Returns porosity at node given by (d, id_i, zlev)
    implicit none
    integer :: d, id_i, zlev

    phi_node = 1d0 + (alpha - 1d0) * penal_node(zlev)%data(d)%elts(id_i)
  end function phi_node

  real(8) function phi_edge (d, id_e, zlev)
    ! Returns porosity at edge given by (d, id_e, zlev)
    implicit none
    integer :: d, id_e, zlev

    phi_edge = 1d0 + (alpha - 1d0) * penal_edge(zlev)%data(d)%elts(id_e)
  end function phi_edge

  real(8) function potential_density (dom, i, j, zlev, offs, dims)
    ! Potential density at nodes (neglect free surface perturbation)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i
    real(8) :: z
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims)

    potential_density = density (dom, i, j, zlev, offs, dims) - ref_density * z / H_rho
  end function potential_density

   real(8) function potential_buoyancy (dom, i, j, zlev, offs, dims)
    ! Potential buoyancy at nodes (neglect free surface perturbation)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i
    real(8) :: z
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims)

    potential_buoyancy = buoyancy (dom, i, j, zlev, offs, dims) + z / H_rho
  end function potential_buoyancy

  real(8) function buoyancy (dom, i, j, zlev, offs, dims)
    ! at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i
    real(8) :: full_mass, full_theta
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    full_mass  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)
    full_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)

    buoyancy = full_theta / full_mass
  end function buoyancy

  subroutine cal_buoyancy (dom, i, j, zlev, offs, dims)
    ! at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i
    real(8) :: full_mass, full_theta
    
    id_i = idx (i, j, offs, dims) + 1

    full_mass  = mean_m(id_i) + mass(id_i)
    full_theta = mean_t(id_i) + temp(id_i)

    scalar(id_i) = full_theta / full_mass
  end subroutine cal_buoyancy

  real(8) function density (dom, i, j, zlev, offs, dims)
    ! Density at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    density = ref_density * (1d0 - buoyancy (dom, i, j, zlev, offs, dims))
  end function density

  real(8) function free_surface (dom, i, j, zlev, offs, dims)
    ! Computes free surface perturbations
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i, k
    real(8) :: full_mass, total_depth

    d = dom%id + 1
    id_i  = idx (i, j, offs, dims) + 1
    
    total_depth = 0d0
    do k = 1, zlevels
       full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
       total_depth = total_depth + full_mass /  porous_density (d, id_i, k)
    end do
    free_surface = total_depth + dom%topo%elts(id_i)
  end function free_surface

  real(8) function interp (e1, e2)
    ! Centred average interpolation of quantities e1 and e2
    implicit none
    real(8) :: e1, e2

    interp = 0.5d0 * (e1 + e2)
  end function interp

  function interp_e (e1, e2)
    ! Centred average interpolation of edge quantities e1 and e2
    implicit none
    real(8), dimension(1:EDGE) :: interp_e
    real(8), dimension(1:EDGE) :: e1, e2

    interp_e = 0.5d0 * (e1 + e2)
  end function interp_e
  
  real(8) function u_mag (dom, i, j, zlev, offs, dims)
    ! Velocity magnitude using data from a single element
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id
    real(8), dimension(1:EDGE) :: prim_dual, u

    id = idx (i, j, offs, dims)

    prim_dual = dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1)
   
    u = sol(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
   
    u_mag = sqrt (sum (u**2 * prim_dual) * dom%areas%elts(id+1)%hex_inv)
  end function u_mag

  subroutine interp_node_edge (dom, i, j, zlev, offs, dims, uvw)
    ! Interpolate from zonal, meridional velocity components at nodes to U, V, W velocity components at edges
    ! (assumes that dom%u_zonal and dom%v_merid have been set over all grid points)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    real(8), dimension (1:EDGE)     :: uvw

    integer     :: id, idE, idN, idNE
    real(8)     :: u_zonal, v_merid
    type(Coord) :: x_e, x_i

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    
    x_i = dom%node%elts(id+1)

    x_e = dom%node%elts(idE+1)
    u_zonal = interp (dom%u_zonal%elts(id+1), dom%u_zonal%elts(idE+1))
    v_merid = interp (dom%v_merid%elts(id+1), dom%v_merid%elts(idE+1))
    uvw(RT+1) = proj_edge ()

    x_e = dom%node%elts(idNE+1)
    u_zonal = interp (dom%u_zonal%elts(id+1), dom%u_zonal%elts(idNE+1))
    v_merid = interp (dom%v_merid%elts(id+1), dom%v_merid%elts(idNE+1))
    uvw(DG+1) = - proj_edge ()

    x_e = dom%node%elts(idN+1)
    u_zonal = interp (dom%u_zonal%elts(id+1), dom%u_zonal%elts(idN+1))
    v_merid = interp (dom%v_merid%elts(id+1), dom%v_merid%elts(idN+1))
    uvw(UP+1) = proj_edge ()
  contains
    real(8) function proj_edge ()
      implicit none
      real(8)     :: lon, lat
      type(Coord) :: co, e_merid, e_zonal, uvw_dir, vel

      ! Find longitude and latitude coordinates
      co = mid_pt (x_i, x_e)
      call cart2sph (co, lon, lat)
      e_zonal = Coord (-sin(lon),           cos(lon),               0d0) ! Zonal direction
      e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

      ! Velocity vector in Cartesian coordinates
      vel = vec_plus (vec_scale (u_zonal, e_zonal), vec_scale (v_merid, e_merid))

      ! Project velocity vector on direction given by points ep1, ep2
      uvw_dir = direction (x_i, x_e)
      proj_edge = inner (uvw_dir, vel)
    end function proj_edge
  end subroutine interp_node_edge

  subroutine interp_edge_node (dom, i, j, zlev, offs, dims)
    ! Interpolate velocity from U, V, W velocity components at edges to zonal, meridional velocity components at nodes
    ! (uses Perot formula as also used for kinetic energy: 
    ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_center) ) / cell_area)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    type(Coord) :: vel, x_e, x_i
    type(Coord) :: e_zonal, e_merid
    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(8)     :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    ! Fluxes normal to hexagon edges
    u_dual_RT    =  velo(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG    = -velo(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    =  velo(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)

    u_dual_RT_W  = -velo(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW =  velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_UP_S  = -velo(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)

    ! Coordinate of hexagon centre (circumcentre)
    x_i = dom%node%elts(id+1)

    ! Sum over 6 hexagon edges
    vel = Coord (0d0, 0d0, 0d0)

    x_e = dom%midpt%elts(EDGE*id+RT+1)
    vel = vec_plus (vel, vec_scale (u_dual_RT,    vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*idW+RT+1)
    vel = vec_plus (vel, vec_scale (u_dual_RT_W,  vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*id+DG+1)
    vel = vec_plus (vel, vec_scale (u_dual_DG,    vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*idSW+DG+1)
    vel = vec_plus (vel, vec_scale (u_dual_DG_SW, vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*id+UP+1)
    vel = vec_plus (vel, vec_scale (u_dual_UP,    vec_minus(x_e, x_i)))

    x_e = dom%midpt%elts(EDGE*idS+UP+1)
    vel = vec_plus (vel, vec_scale (u_dual_UP_S,  vec_minus(x_e, x_i)))

    vel = vec_scale (dom%areas%elts(id+1)%hex_inv, vel) ! construct velocity at hexagonal node

    ! Project velocity onto zonal and meridional directions
    call cart2sph (x_i, lon, lat)

    e_zonal = Coord (-sin(lon),           cos(lon),               0d0) ! Zonal direction
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat)) ! Meridional direction

    velo1(id+1) = inner (vel, e_zonal)
    velo2(id+1) = inner (vel, e_merid)
  end subroutine interp_edge_node

  subroutine vort_triag_to_hex (dom, i, j, zlev, offs, dims)
    ! Approximate vorticity at hexagon points
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idW, idSW, idS, d

    d = dom%id + 1
    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    vort(id+1) = ( &
         dom%areas%elts(id+1)%part(1)*dom%vort%elts(TRIAG*id+LORT+1)   + &
         dom%areas%elts(id+1)%part(2)*dom%vort%elts(TRIAG*id+UPLT+1)   + &
         dom%areas%elts(id+1)%part(3)*dom%vort%elts(TRIAG*idW+LORT+1)  + &
         dom%areas%elts(id+1)%part(4)*dom%vort%elts(TRIAG*idSW+UPLT+1) + &
         dom%areas%elts(id+1)%part(5)*dom%vort%elts(TRIAG*idSW+LORT+1) + &
         dom%areas%elts(id+1)%part(6)*dom%vort%elts(TRIAG*idS+UPLT+1)    &
         ) * dom%areas%elts(id+1)%hex_inv
  end subroutine vort_triag_to_hex
  
  function f_coriolis_edge (dom, i, j, zlev, offs, dims)
    ! Coriolis parameter at edges
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: f_coriolis_edge

    integer :: id

    id = idx (i, j, offs, dims)
    
    f_coriolis_edge = dom%midpt%elts(EDGE*id+RT+1:EDGE*id+UP+1)%z / radius * 2d0*omega
  end function f_coriolis_edge

  real(8) function f_coriolis_node (dom, i, j, zlev, offs, dims)
    ! Coriolis parameter at nodes
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims)

    f_coriolis_node = dom%node%elts(id+1)%z / radius * 2d0*omega
  end function f_coriolis_node

  real(8) function integrate_hex (fun, l, zlev)
    ! Integrate function defined by fun over hexagons
    implicit none
    integer  :: l, zlev

    integer                        :: d, ll, p, i, j, c, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s, val

    interface
       real(8) function fun (dom, i, j, zlev, offs, dims)
         import
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    s = 0d0
    do d = 1, size(grid)
       do ll = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(ll)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx (i-1, j-1, offs, dims)
                val = fun (grid(d), i-1, j-1, zlev, offs, dims)
                s = s + val/grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do

       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%pole_master(c/2-2) .or. .not. grid(d)%penta(c)) cycle
          p = 1
          do while (grid(d)%patch%elts(p+1)%level < l)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p == 0) then
                write (6,'(A, i4, A)') "ERROR(rank = ", rank, "): integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)
          if (c == NORTHWEST) then
             id = idx (0, PATCH_SIZE, offs, dims)
             val = fun (grid(d), 0, PATCH_SIZE, zlev, offs, dims)
             s = s + val/grid(d)%areas%elts(id+1)%hex_inv
          else
             id = idx (PATCH_SIZE, 0, offs, dims)
             val = fun (grid(d), PATCH_SIZE, 0, zlev, offs, dims)
             s = s + val/grid(d)%areas%elts(id+1)%hex_inv
          end if
       end do
    end do
    integrate_hex = sum_real (s)
  end function integrate_hex

  real(8) function integrate_tri (fun, k)
    implicit none
    integer  :: k

    integer                        :: d, l, ll, p, i, j, t, id
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8)                        :: s

    interface
       real(8) function fun (dom, i, j, zlev, t, offs, dims)
         import 
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, t, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    s = 0d0
    do l = level_start, level_end
       do d = 1, size(grid)
          do ll = 1, grid(d)%lev(l)%length
             p = grid(d)%lev(l)%elts(ll)
             call get_offs_Domain (grid(d), p, offs, dims)
             do j = 1, PATCH_SIZE
                do i = 1, PATCH_SIZE
                   id = idx (i-1, j-1, offs, dims)
                   do t = LORT, UPLT
                      s = s + fun (grid(d), i-1, j-1, k, t, offs, dims) * grid(d)%triarea%elts(TRIAG*id+t+1)
                   end do
                end do
             end do
          end do
       end do
    end do
    integrate_tri = sum_real (s)
  end function integrate_tri

  real(8) function integrate_adaptive (routine, zlev)
    ! Integrates value define in routine over adaptive grid
    implicit none
    integer           :: zlev
    real(8), external :: routine

    integer ::  l

    hex_int = 0d0
    call fine_hex_area (routine, level_start, zlev)
    do l = level_start+1, level_end-1
       call fine_hex_area (routine, l, zlev)
       call coarse_hex_area (routine, l, zlev)
    end do
    call fine_hex_area (routine, level_end, zlev)
    integrate_adaptive = sum_real (hex_int)
  end function integrate_adaptive

  subroutine fine_hex_area (routine, l, zlev)
    implicit none
    integer           :: l, zlev
    real(8), external :: routine

    integer                          :: d, i, id, j, jj, p
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(jj)
          call get_offs_Domain (grid(d), p, offs, dims, inner_bdry)
          do j = 1, PATCH_SIZE 
             do i = 1, PATCH_SIZE
                id = idx (i, j, offs, dims) + 1
                hex_int = hex_int + routine (grid(d), i, j, zlev, offs, dims) / grid(d)%areas%elts(id)%hex_inv
             end do
          end do
       end do
    end do
  end subroutine fine_hex_area

  subroutine coarse_hex_area (routine, l, zlev)
    ! Remove cells that are not at locally finest scale
    implicit none
    integer           :: l, zlev
    real(8), external :: routine

    integer                          :: c, d, i, id, i_par, j, j_par, jj, p, p_chd, p_par, s
    integer, dimension(N_BDRY+1)     :: offs, offs_chd, offs_par
    integer, dimension(2,N_BDRY+1)   :: dims, dims_chd, dims_par
    integer, dimension(JPlUS:IMINUS) :: bdry
    logical, dimension(JPlUS:IMINUS) :: inner_bdry

    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          p_par = grid(d)%lev(l)%elts(jj)
          call get_offs_Domain (grid(d), p_par, offs_par, dims_par)
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd == 0) cycle
             call get_offs_Domain (grid(d), p_chd, offs_chd, dims_chd, inner_bdry)
             do j = 1, PATCH_SIZE/2
                j_par = j-1 + chd_offs(2,c)
                do i = 1, PATCH_SIZE/2
                   i_par = i-1 + chd_offs(1,c)
                   id = idx (i_par, j_par, offs_par, dims_par) + 1
                   hex_int = hex_int - routine (grid(d), i_par, j_par, zlev, offs_par, dims_par) &
                        / grid(d)%areas%elts(id)%hex_inv
                end do
             end do
          end do
       end do
    end do
  end subroutine coarse_hex_area

  subroutine set_areas (dom, id, val)
    implicit none
    type(Domain)          :: dom
    integer               :: id
    real(8), dimension(7) :: val

    real(8), dimension(4) :: area

    area = val(1:4)
    if (id < 0) area = (/area(2), area(1), area(4), area(3)/)

    dom%overl_areas%elts(abs(id) + 1)%a     = area
    dom%overl_areas%elts(abs(id) + 1)%split = val(5:6)
    dom%areas%elts(abs(id) + 1)%hex_inv     = val(7)
  end subroutine set_areas

  subroutine get_areas (dom, id, val)
    implicit none
    real(8), dimension(7), intent(out) :: val
    type(Domain)                       :: dom
    integer                            :: id

    real(8), dimension(7) :: area

    area = 0d0
    area(1:4) = dom%overl_areas%elts(id+1)%a
    area(5:6) = dom%overl_areas%elts(id+1)%split
    area(7)   = dom%areas%elts(id+1)%hex_inv
    val       = area
    return
  end subroutine get_areas

  subroutine area_post_comm (dom, p, c, offs, dims, zlev)
    implicit none
    type(Domain)                   :: dom
    integer                        :: c, p, zlev
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    integer, dimension(N_BDRY+1)   :: offs

    if (c == IPLUSJMINUS) then
       id = idx(PATCH_SIZE, -1, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0d0
    end if
    if (c == IMINUSJPLUS) then
       id = idx(-1, PATCH_SIZE, offs, dims)
       dom%overl_areas%elts(id+1)%a = 0d0
    end if
  end subroutine area_post_comm
end module utils_mod
