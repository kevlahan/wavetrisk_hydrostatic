module utils_mod
  ! Basic functions
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  implicit none
  real(8) :: integral

  abstract interface
     real(8) function routine_hex (dom, i, j, zlev, offs, dims)
       use domain_mod
       implicit none
       type (Domain)                  :: dom
       integer                        :: i, j, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function routine_hex
     real(8) function routine_tri (dom, i, j, t, zlev, offs, dims)
       use domain_mod
       implicit none
       type (Domain)                  :: dom
       integer                        :: i, j, t, zlev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function routine_tri
  end interface
  procedure (routine_hex), pointer :: integrand_hex => null ()
  procedure (routine_tri), pointer :: integrand_tri => null ()
contains
  real(8) function z_i (dom, i, j, zlev, offs, dims, q)
    ! Position of vertical level zlev at nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: id, k
    real(8) :: dz, dz_below, z_s

    id = idx (i, j, offs, dims)

    z_s = dom%topo%elts(id+1) 
    
    dz_below = dz_i (dom, i, j, 1, offs, dims, q)
    z_i = z_s + dz_below / 2
    do k = 2, zlev
       dz = dz_i (dom, i, j, k, offs, dims, q)
       z_i = z_i + interp (dz, dz_below)
       dz_below = dz
    end do
  end function z_i

  real(8) function dz_i (dom, i, j, zlev, offs, dims, q)
    ! Thickness of layer zlev at nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    dz_i = (sol_mean(S_MASS,zlev)%data(d)%elts(id) + q(S_MASS,zlev)%data(d)%elts(id)) / porous_density (d, id, zlev)
  end function dz_i
  
  function dz_e (dom, i, j, zlev, offs, dims, q)
    ! Thickness of layer zlev at edges
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(8), dimension(1:EDGE)                :: dz_e

    integer                    :: d, id, idE, idN, idNE
    real(8), dimension(0:EDGE) :: dz

    d = dom%id + 1

    id   = idx (i, j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
   
    dz(0) = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1)   + q(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (d, id+1, zlev)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)  + q(S_MASS,zlev)%data(d)%elts(idE+1)) &
         / porous_density (d, idE+1, zlev)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1) + q(S_MASS,zlev)%data(d)%elts(idNE+1)) &
         / porous_density (d, idNE+1, zlev)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)  + q(S_MASS,zlev)%data(d)%elts(idN+1)) &
         / porous_density (d, idN+1, zlev)

    dz_e = 0.5d0 * (dz(0) + dz(1:EDGE))
  end function dz_e

  function dz_SW_e (dom, i, j, zlev, offs, dims, q)
    ! Thickness of layer zlev at edges (SW edges)
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(8), dimension(1:EDGE)                :: dz_SW_e

    integer                    :: d, id, idW, idSW, idS
    real(8), dimension(0:EDGE) :: dz

    d = dom%id + 1

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
   
    dz(0) = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + q(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (d, id+1, zlev)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idW+1)  + q(S_MASS,zlev)%data(d)%elts(idW+1)) &
         / porous_density (d, idW+1, zlev)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idSW+1) + q(S_MASS,zlev)%data(d)%elts(idSW+1)) &
         / porous_density (d, idSW+1, zlev)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idS+1)  + q(S_MASS,zlev)%data(d)%elts(idS+1)) &
         / porous_density (d, idS+1, zlev)

    dz_SW_e = 0.5d0 * (dz(0) + dz(1:EDGE))
  end function dz_SW_e
  
  real(8) function zl_i (dom, i, j, zlev, offs, dims, q, l)
    ! Position of interface below (l=-1) or above (l=1) vertical level zlev nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, l, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: id, k

    id = idx (i, j, offs, dims)

    zl_i = dom%topo%elts(id+1)
    do k = 1, zlev+l
       zl_i = zl_i + dz_i (dom, i, j, k, offs, dims, q)
    end do
  end function zl_i

  function zl_e (dom, i, j, zlev, offs, dims, q, l)
    ! Position of interface below (l=-1) or above (l=1) vertical level zlev at edges
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer                                   :: l
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(8), dimension(1:EDGE)                :: zl_e

    integer :: id, idE, idN, idNE, k

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    zl_e(RT+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idE+1))  
    zl_e(DG+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idNE+1)) 
    zl_e(UP+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idN+1))  
    do k = 1, zlev+l
       zl_e = zl_e + dz_e (dom, i, j, k, offs, dims, q)
    end do
  end function zl_e

  real(8) function dz_l (dom, i, j, zlev, offs, dims, q)
    ! Thickness of layer associated with interface between layers zlev and zlev+1: z_(zlev+1) - z_zlev
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    
    real(8) :: dZ, dZ_above

    dZ        = dz_i (dom, i, j, zlevels,   offs, dims, q)
    dZ_above  = dz_i (dom, i, j, zlevels+1, offs, dims, q)

    dz_l = 0.5d0 * (dZ + dZ_above)
  end function dz_l

  function eta_e (dom, i, j, zlev, offs, dims, q)
    ! Free surface perturbation at edges
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(8), dimension(1:EDGE)                :: eta_e

    integer :: d, id, idE, idN, idNE
    real(8) :: eta0

    d = dom%id + 1
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (mode_split) then
       eta_e(RT+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idE+1))
       eta_e(DG+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idNE+1))
       eta_e(UP+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idN+1))
    else
       eta0 = free_surface (dom, i, j, zlev, offs, dims, q)
       eta_e(RT+1) = interp (eta0, free_surface (dom, i+1, j,   zlev, offs, dims, q))
       eta_e(DG+1) = interp (eta0, free_surface (dom, i+1, j+1, zlev, offs, dims, q))
       eta_e(UP+1) = interp (eta0, free_surface (dom, i,   j+1, zlev, offs, dims, q))
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

  real(8) function potential_density (dom, i, j, zlev, offs, dims, q)
    ! Potential density at nodes (neglect free surface perturbation)
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, id_i
    real(8) :: z
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims, q)

    potential_density = density (dom, i, j, zlev, offs, dims, q) - ref_density * z / H_rho
  end function potential_density

   real(8) function potential_buoyancy (dom, i, j, zlev, offs, dims, q)
    ! Potential buoyancy at nodes (neglect free surface perturbation)
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, id_i
    real(8) :: z
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims, q)

    potential_buoyancy = buoyancy (dom, i, j, zlev, offs, dims, q) + z / H_rho
  end function potential_buoyancy

  real(8) function buoyancy (dom, i, j, zlev, offs, dims, q)
    ! at nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, id_i
    real(8) :: full_mass, full_theta
    
    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    full_mass  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + q(S_MASS,zlev)%data(d)%elts(id_i)
    full_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + q(S_TEMP,zlev)%data(d)%elts(id_i)

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

  real(8) function density (dom, i, j, zlev, offs, dims, q)
    ! Density at nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    density = ref_density * (1d0 - buoyancy (dom, i, j, zlev, offs, dims, q))
  end function density

  real(8) function free_surface (dom, i, j, zlev, offs, dims, q)
    ! Computes free surface perturbations
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, id_i, k
    real(8) :: full_mass, total_depth

    d = dom%id + 1
    id_i  = idx (i, j, offs, dims) + 1
    
    total_depth = 0d0
    do k = 1, zlevels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id_i) + q(S_MASS,k)%data(d)%elts(id_i)
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

  real(8) function integrate_tri (fun, zlev, coarse_only)
    ! Integrate over adaptive triangles, where the integrand is defined by the routine fun.
    ! If optional variable coarse_only = .true. the integration is carried out over level_start only.
    implicit none
    integer           :: zlev
    logical, optional :: coarse_only
    
    integer :: d, j, l
    logical :: finer
    
    interface
       real(8) function fun (dom, i, j, t, zlev, offs, dims)
         use domain_mod
         implicit none
         type (Domain)                  :: dom
         integer                        :: i, j, t, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    if (present(coarse_only)) then
       finer = .not. coarse_only
    else
       finer = .true.
    end if
    integrand_tri => fun

    ! Coarsest level
    integral = 0d0
    call apply_onescale (integrate_tri_coarse, level_start, zlev, 0, 0)

    ! Finer levels
    if (finer) then
       do l = level_start, level_end-1
          do d = 1, size(grid)
             do j = 1, grid(d)%lev(l)%length
                call apply_interscale_to_patch (integrate_tri_fine, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
             end do
          end do
       end do
    end if
    nullify (integrand_tri)
    
    integrate_tri = sum_real (integral)
  end function integrate_tri
    
  subroutine integrate_tri_coarse (dom, i, j, zlev, offs, dims)
    ! Integral over a single level.
    implicit none
    type(Domain),                   intent(in) :: dom
    integer,                        intent(in) :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer :: id, t

    id = idx (i, j, offs, dims)

    do t = LORT, UPLT
       integral = integral + integrand_tri (dom, i, j, t, zlev, offs, dims) * dom%triarea%elts(TRIAG*id+t+1)
    end do
  end subroutine integrate_tri_coarse

  subroutine integrate_tri_fine (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Integrate over active finer levels.
    implicit none
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    integer                       :: id_chd, idE_chd, idNE_chd, idN_chd, t
    real(8), dimension(LORT:UPLT) :: parent_integrand

    id_chd   = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) then
       do t = LORT, UPLT
          parent_integrand(t) = integrand_tri (dom, i_par, j_par, t, zlev, offs_par, dims_par)
       end do

       do t = LORT, UPLT
          integral = integral + (integrand_tri (dom, i_chd, j_chd, t, zlev, offs_chd, dims_chd) - parent_integrand(t)) &
               * dom%triarea%elts(TRIAG*id_chd+t+1)

          integral = integral + (integrand_tri (dom, i_chd+1, j_chd, t, zlev, offs_chd, dims_chd) - parent_integrand(LORT)) &
               * dom%triarea%elts(TRIAG*idE_chd+t+1)

          integral = integral + (integrand_tri (dom, i_chd, j_chd+1, t, zlev, offs_chd, dims_chd) - parent_integrand(UPLT)) &
               * dom%triarea%elts(TRIAG*idN_chd+t+1)

          integral = integral + (integrand_tri (dom, i_chd+1, j_chd+1, t, zlev, offs_chd, dims_chd) - parent_integrand(t)) &
               * dom%triarea%elts(TRIAG*idNE_chd+t+1)
       end do
    end if
  end subroutine integrate_tri_fine
  
  real(8) function integrate_hex (fun, zlev, coarse_only)
    ! Integrate over adaptive hexagons, where the integrand is defined by the routine fun.
    ! If optional variable coarse_only = .true. the integration is carried out over level_start only.
    implicit none
    integer           :: zlev
    logical, optional :: coarse_only

    integer :: d, j, l
    logical :: finer

    interface
       real(8) function fun (dom, i, j, zlev, offs, dims)
         use domain_mod
         implicit none
         type (Domain)                  :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function fun
    end interface

    if (present(coarse_only)) then
       finer = .not. coarse_only
    else
       finer = .true.
    end if
    integrand_hex => fun

    ! Coarsest level
    integral = 0d0
    call integrate_hex_coarse (zlev)
    
    ! Finer levels
    if (finer) then
       do l = level_start, level_end-1
          do d = 1, size(grid)
             do j = 1, grid(d)%lev(l)%length
                call apply_interscale_to_patch (integrate_hex_fine, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
             end do
          end do
       end do
    end if
    nullify (integrand_hex)
    
    integrate_hex = sum_real (integral)
  end function integrate_hex

  subroutine integrate_hex_coarse (zlev)
    ! Integrate function pointer integrand_hex over hexagons at coarsest scale.
    implicit none
    integer :: zlev
    
    integer                        :: c, d, i, id, j, jj, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       ! Regular hexagons/pentagons
       do jj = 1, grid(d)%lev(level_start)%length
          p = grid(d)%lev(level_start)%elts(jj)
          call get_offs_Domain (grid(d), p, offs, dims)
          do j = 1, PATCH_SIZE
             do i = 1, PATCH_SIZE
                id = idx (i-1, j-1, offs, dims)
                integral = integral + integrand_hex (grid(d), i-1, j-1, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
             end do
          end do
       end do

       ! Check domain d to see if its SOUTHEAST or NORTHWEST corners have associated poles
       do c = SOUTHEAST, NORTHWEST, 2
          if (.not. grid(d)%pole_master(c/2-2) .or. .not. grid(d)%penta(c)) cycle
          p = 1
          do while (grid(d)%patch%elts(p+1)%level < level_start)
             p = grid(d)%patch%elts(p+1)%children(c-4)
             if (p == 0) then
                write (6,'(A, i4, A)') "ERROR(rank = ", rank, "): integrate_hex: level incomplete"
                return
             end if
          end do
          call get_offs_Domain (grid(d), p, offs, dims)

          if (c == NORTHWEST) then     ! north pole
             id = idx (0, PATCH_SIZE, offs, dims)
             integral = integral + integrand_hex (grid(d), 0, PATCH_SIZE, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
          elseif (c == SOUTHEAST) then ! south pole
             id = idx (PATCH_SIZE, 0, offs, dims)
             integral = integral + integrand_hex (grid(d), PATCH_SIZE, 0, zlev, offs, dims) / grid(d)%areas%elts(id+1)%hex_inv
          end if

       end do
    end do
  end subroutine integrate_hex_coarse

  subroutine integrate_hex_fine (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Integrate over active finer levels.
    implicit none
    type(Domain)                     :: dom
    integer                          :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_par, dims_chd

    real(8), dimension(LORT:UPLT) :: parent_integrand
    integer :: id_par, id_chd, idE_chd, idNE_chd, idN_chd, t

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) then
       do t = LORT, UPLT
          parent_integrand(t) = hex2tri (dom, i_par, j_par, t, offs_par, dims_par, zlev)
       end do
       
       do t = LORT, UPLT
          integral = integral + (hex2tri (dom, i_chd, j_chd, t, offs_chd, dims_chd, zlev) - parent_integrand(t)) &
               * dom%triarea%elts(TRIAG*id_chd+t+1)

          integral = integral + (hex2tri (dom, i_chd+1, j_chd, t, offs_chd, dims_chd, zlev) - parent_integrand(LORT)) &
               * dom%triarea%elts(TRIAG*idE_chd+t+1)

          integral = integral + &
               (hex2tri (dom, i_chd, j_chd+1, t, offs_chd, dims_chd, zlev) - parent_integrand(UPLT)) &
               * dom%triarea%elts(TRIAG*idN_chd+t+1)

          integral = integral + (hex2tri (dom, i_chd+1, j_chd+1, t, offs_chd, dims_chd, zlev) - parent_integrand(t)) &
               * dom%triarea%elts(TRIAG*idNE_chd+t+1)
       end do
    end if
  end subroutine integrate_hex_fine

  real(8) function hex2tri (dom, i, j, t, offs, dims, zlev)
    ! Integrand at triangles associated with node (i,j) computed from integral over hexagons
    implicit none
    integer :: i, j, t, zlev
    type(Domain)                     :: dom
    integer, dimension(N_BDRY + 1)   :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id, idE, idNE, idN

    id   = idx (i,   j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    if (t == LORT) then
       idE = idx (i+1, j, offs, dims)
       hex2tri = &
            integrand_hex (dom, i,   j,   zlev, offs, dims) * dom%areas%elts(id+1)%part(1)  + &
            integrand_hex (dom, i+1, j+1, zlev, offs, dims) * dom%areas%elts(idNE+1)%part(5) + &
            integrand_hex (dom, i+1, j,   zlev, offs, dims) * dom%areas%elts(idE+1)%part(3)
    elseif (t == UPLT) then
       idN = idx (i, j+1, offs, dims)
       hex2tri = &
            integrand_hex (dom, i,   j,   zlev, offs, dims) * dom%areas%elts(id+1)%part(2)  + &
            integrand_hex (dom, i+1, j+1, zlev, offs, dims) * dom%areas%elts(idNE+1)%part(4) + &
            integrand_hex (dom, i,   j+1, zlev, offs, dims) * dom%areas%elts(idN+1)%part(6)
    end if

    hex2tri = hex2tri / dom%triarea%elts(TRIAG*id+t+1)
  end function hex2tri
end module utils_mod
