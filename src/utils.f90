module utils_mod
  ! Basic functions
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  implicit none
  real(8)                        :: integral
  real(8), dimension(:), pointer :: val1, val2

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

  interface zero_float
     procedure :: zero_float_0, zero_float_1, zero_float_2
  end interface zero_float
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
    z_i = z_s + dz_below / 2d0
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

    integer :: id, k, kmax

    id = idx (i, j, offs, dims)

    if (l == -1) then
       kmax = zlev - 1
    else
       kmax = zlev
    end if

    zl_i = dom%topo%elts(id+1)
    do k = 1, kmax
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

    integer :: id, idE, idN, idNE, k, kmax

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    if (l == -1) then
       kmax = zlev - 1
    else
       kmax = zlev
    end if

    zl_e(RT+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idE+1))  
    zl_e(DG+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idNE+1)) 
    zl_e(UP+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idN+1))  
    do k = 1, kmax
       zl_e = zl_e + dz_e (dom, i, j, k, offs, dims, q)
    end do
  end function zl_e

  real(8) function dz_l (dom, i, j, zlev, offs, dims, q)
    ! Thickness of layer associated with interface between layers zlev and zlev+1: z_(zlev+1) - z(zlev)
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    
    real(8) :: dZ, dZ_above

    dZ        = dz_i (dom, i, j, zlev,   offs, dims, q)
    dZ_above  = dz_i (dom, i, j, zlev+1, offs, dims, q)

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
    
    porous_density = ref_density * phi_node (d, id_i, zlev)
  end function porous_density

  function porous_density_edge (d, id, zlev)
    ! Porous density at edges
    implicit none
    integer                    :: d, id, zlev
    real(8), dimension(1:EDGE) :: porous_density_edge
    
    porous_density_edge = ref_density * phi_edge (d, id, zlev)
  end function porous_density_edge

  real(8) function phi_node (d, id_i, zlev)
    ! Returns porosity at node given by (d, id_i, zlev)
    implicit none
    integer :: d, id_i, zlev

    phi_node = 1d0 + (alpha - 1d0) * penal_node(zlev)%data(d)%elts(id_i)
  end function phi_node

  function phi_edge (d, id, zlev)
    ! Returns porosity at edges associated to node given by (d, id_i, zlev)
    implicit none
    integer                    :: d, e, id, id_e, zlev
    real(8), dimension(1:EDGE) :: phi_edge

    do e = 1, EDGE
       id_e = EDGE*id+e
       phi_edge(e) = 1d0 + (alpha - 1d0) * penal_edge(zlev)%data(d)%elts(id_e)
    end do
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

  subroutine interp_latlon_UVW (dom, i, j, zlev, offs, dims, uvw)
    ! Interpolate from zonal, meridional velocity components at nodes to U, V, W velocity components at edges
    ! (assumes that dom%u_zonal and dom%v_merid have been set over all grid points)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    real(8), dimension (1:EDGE)     :: uvw

    integer     :: id, idE, idN, idNE
    type(Coord) :: vel0

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    vel0 = vel (id)

    uvw(RT+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idE+1)),  0.5d0*(vel0 + vel(idE)))
    uvw(DG+1) = inner (direction (dom%node%elts(idNE+1), dom%node%elts(id+1)) ,  0.5d0*(vel0 + vel(idNE)))
    uvw(UP+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idN+1)),  0.5d0*(vel0 + vel(idN)))
  contains
    type(Coord) function vel (id)
      ! Computes velocity at node id from its latitude and longitude components
      integer :: id

      real(8)     :: lat, lon
      type(Coord) :: e_merid, e_zonal

      call cart2sph (dom%node%elts(id+1), lon, lat)

      e_zonal = Coord (-sin(lon),           cos(lon),               0d0) 
      e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))

      vel = dom%u_zonal%elts(id+1) * e_zonal + dom%v_merid%elts(id+1) * e_merid
    end function vel
  end subroutine interp_latlon_UVW

  subroutine interp_UVW_latlon (dom, i, j, zlev, offs, dims)
    ! Interpolate velocity from U, V, W velocity components at edges to zonal, meridional velocity components at nodes
    ! Perot reconstruction based on Gauss theorem:
    !
    ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_centroid) ) / cell_area
    !
    ! also used for kinetic energy
    !
    ! Output is in pointer arrays velo1 (u_zonal) and velo2 (u_merid)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    
    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(8)     :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW
    type(Coord) :: cent, e_zonal, e_merid, vel

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

    ! Compute hexagon centroid from its vertices
    cent = centroid (                                                                 &
         (/ dom%ccentre%elts(TRIAG*id+LORT+1),   dom%ccentre%elts(TRIAG*id+UPLT+1),   &
            dom%ccentre%elts(TRIAG*idW+LORT+1),  dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS+UPLT+1) /), 6)

    ! Velocity at node from Perot formula
    vel = dom%areas%elts(id+1)%hex_inv * ( &
         u_dual_RT    * (dom%midpt%elts(EDGE*id+RT+1)   - cent) + &
         u_dual_DG    * (dom%midpt%elts(EDGE*id+DG+1)   - cent) + &
         u_dual_UP    * (dom%midpt%elts(EDGE*id+UP+1)   - cent) + &
         u_dual_RT_W  * (dom%midpt%elts(EDGE*idW+RT+1)  - cent) + &
         u_dual_DG_SW * (dom%midpt%elts(EDGE*idSW+DG+1) - cent) + &
         u_dual_UP_S  * (dom%midpt%elts(EDGE*idS+UP+1)  - cent))

    ! Coordinate of hexagon centre (circumcentre)
    call cart2sph (dom%node%elts(id+1), lon, lat)

    ! Zonal and meridional directions
    e_zonal = Coord (-sin(lon),           cos(lon),               0d0) 
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))
    
    ! Project velocity at node onto zonal and meridional directions
    velo1(id+1) = inner (vel, e_zonal)
    velo2(id+1) = inner (vel, e_merid)
  end subroutine interp_UVW_latlon

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
  
  real(8) function integrate_hex (fun, zlev, level)
    ! Integrate over adaptive hexagons, where the integrand is defined by the routine fun.
    ! If optional variable coarse_only = .true. the integration is carried out over level_start only.
    implicit none
    integer           :: zlev
    integer, optional :: level

    integer :: j, d, l

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

    integrand_hex => fun
    
    integral = 0d0
    
    if (present (level)) then ! integrate over a single scales
       call integrate_hex_scale (level, zlev)
    else                      ! integrate over all scales
       call integrate_hex_scale (level_start, zlev)
       do l = level_start, level_end-1
          do d = 1, size(grid)
             do j = 1, grid(d)%lev(l)%length
                call apply_interscale_to_patch (integrate_hex_fine, grid(d), grid(d)%lev(l)%elts(j), zlev, 0, 0)
             end do
          end do
       end do
    end if
    
    integrate_hex = sum_real (integral)
    
    nullify (integrand_hex)
  end function integrate_hex

  subroutine integrate_hex_scale (l, zlev)
    ! Integrate function pointer integrand_hex over hexagons at a single scale l
    implicit none
    integer :: l, zlev
    
    integer                        :: c, d, i, id, j, jj, p
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    do d = 1, size(grid)
       ! Regular hexagons/pentagons
       do jj = 1, grid(d)%lev(l)%length
          p = grid(d)%lev(l)%elts(jj)
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
          do while (grid(d)%patch%elts(p+1)%level < l)
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
  end subroutine integrate_hex_scale

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
  
  subroutine zero_float_0 (q)
    ! Initializes a float field to zero
    implicit none
    type(Float_Field) :: q

    integer :: d

    do d = 1, size(grid)
       q%data(d)%elts = 0d0
    end do
  end subroutine zero_float_0

  subroutine zero_float_1 (q)
    ! Initializes a float vector to zero
    implicit none
    type(Float_Field), dimension(:) :: q

    integer :: d, j

    do j = 1, size(q,1)
       do d = 1, size(grid)
          q(j)%data(d)%elts = 0d0
       end do
    end do
  end subroutine zero_float_1

   subroutine zero_float_2 (q)
    ! Initializes a float array to zero
    implicit none
    type(Float_Field), dimension(:,:) :: q

    integer :: d, j1, j2

    do j1 = 1, size(q,1)
       do j2 = 1, size(q,2)
          do d = 1, size(grid)
             q(j1,j2)%data(d)%elts = 0d0
          end do
       end do
    end do
  end subroutine zero_float_2
  
  subroutine zero_float_vector (q)
    ! Initializes a float vector to zero
    implicit none
    type(Float_Field), dimension(:) :: q

    integer :: d, j

    do j = 1, size(q,1)
       do d = 1, size(grid)
          q(j)%data(d)%elts = 0d0
       end do
    end do
  end subroutine zero_float_vector

  subroutine zero_float_array (q)
    ! Initializes a float array to zero
    implicit none
    type(Float_Field), dimension(:,:) :: q

    integer :: d, j1, j2

    do j1 = 1, size(q,1)
       do j2 = 1, size(q,2)
          do d = 1, size(grid)
             q(j1,j2)%data(d)%elts = 0d0
          end do
       end do
    end do
  end subroutine zero_float_array

  subroutine zero_float_field (q, itype, lmin_in, lmax_in)
    ! Set float field to zero for scales:
    ! (lmin,lmax) if both lmin and lmax are present
    ! lmin if only lmin is present
    ! (level_start, level_end) if lmin is not present
    ! itype = S_MASS or S_VELO
    implicit none
    integer                                  :: itype
    integer, optional                        :: lmin_in, lmax_in
    type(Float_Field), target, intent(inout) :: q

    integer :: d, j, l, lmin, lmax
    
    if (present(lmin_in)) then
       lmin = lmin_in
       if (present(lmax_in)) then
          lmax = lmax_in
       else
          lmax = lmin_in
       end if
    else
       lmin = level_start
       lmax = level_end
    end if

    if (itype == AT_NODE) then
       do l = lmin, lmax
          do d = 1, size(grid)
             val1 => q%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_zero_node, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (val1)
          end do
       end do
    elseif (itype == AT_EDGE) then
       do l = lmin, lmax
          do d = 1, size(grid)
             val1 => q%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_zero_edge, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
             end do
             nullify (val1)
          end do
       end do
    else
       if (rank == 0) write (6,'(a)') "Unsupported type for zero_float_field ... aborting"
       call abort
    end if
    q%bdry_uptodate = .false.
    call update_bdry1 (q, lmin, lmax)
  end subroutine zero_float_field

  subroutine cal_zero_node (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    val1(id) = 0d0
  end subroutine cal_zero_node

  subroutine cal_zero_edge (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx (i, j, offs, dims)

    do e = 1, EDGE
       val1(EDGE*id+e) = 0d0
    end do
  end subroutine cal_zero_edge

  subroutine equals_float_field (q1, q2, itype, lmin_in, lmax_in)
    ! Set elements of float field q1 = q2
    !
    ! itype = S_MASS or S_VELO
    ! if scale l is present, compute only for scale l
    implicit none
    integer                                  :: itype
    integer, optional                        :: lmin_in, lmax_in
    type(Float_Field), target, intent(in)    :: q2
    type(Float_Field), target, intent(inout) :: q1

    integer :: d, j, l, lmin, lmax
    
    if (present(lmin_in)) then
       lmin = lmin_in
       if (present(lmax_in)) then
          lmax = lmax_in
       else
          lmax = lmin_in
       end if
    else
       lmin = level_start
       lmax = level_end
    end if

    if (itype == AT_NODE) then
       do l = lmin, lmax
          do d = 1, size(grid)
             val1 => q1%data(d)%elts
             val2 => q2%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_equals_node, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (val1, val2)
          end do
       end do
    elseif (itype == AT_EDGE) then
       do l = lmin, lmax
          do d = 1, size(grid)
             val1 => q1%data(d)%elts
             val2 => q2%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_equals_edge, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
             end do
             nullify (val1, val2)
          end do
       end do
    else
       if (rank == 0) write (6,'(a)') "Unsupported type for zero_float_field ... aborting"
       call abort
    end if
    q1%bdry_uptodate = .false.
    call update_bdry1 (q1, lmin, lmax)
  end subroutine equals_float_field

   subroutine cal_equals_node (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    val1(id) = val2(id)
  end subroutine cal_equals_node

  subroutine cal_equals_edge (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx (i, j, offs, dims) 

    do e = 1, EDGE
       val1(EDGE*id+e) = val2(EDGE*id+e)
    end do
  end subroutine cal_equals_edge
  
  subroutine smoothing_rbf (dx, npts, nsmth, data)
    ! Smooths data(lon,lat) over neighbouring region using radial basis functions
    implicit none
    integer                 :: npts, nsmth
    real(8)                 :: dx
    real(8), dimension(:,:) :: data

    integer                              :: i, ismth, i0, ii, j, j0, jj, nx, ny
    real(8)                              :: r, M_topo, sw_topo, topo_sum, wgt
    real(8), allocatable, dimension(:,:) :: data_old

    nx = size(data,1); ny = size(data,2)
    allocate (data_old(nx,ny))

    data_old = data

    do ismth = 1, nsmth
       do i = 1, nx
          do j = 1, ny
             sw_topo  = 0d0
             topo_sum = 0d0
             do ii = -npts, npts
                do jj = -npts, npts

                   r = sqrt (dble(ii**2 + jj**2)) * MATH_PI * radius / dble (ny)

                   wgt = radial_basis_fun (dx, npts, r)

                   call wrapij (i+ii, j+jj, nx, ny, i0, j0)

                   M_topo = data_old(i0, j0)

                   topo_sum = topo_sum + wgt * M_topo
                   sw_topo  = sw_topo  + wgt
                end do
             end do
             data(i,j) = topo_sum / sw_topo
          end do
       end do
       data_old = data
    end do
    deallocate (data_old)
  end subroutine smoothing_rbf

  subroutine smoothing_shapiro (nsmth, data)
    ! Smooths data(lon,lat) over neighbouring region using shapiro (diffusion) filter
    implicit none
    integer                 :: nsmth
    real(8), dimension(:,:) :: data

    integer                              :: i, ii, ismth, jj, nx, ny
    real(8), allocatable, dimension(:,:) :: data_old

    nx = size(data,1); ny = size(data,2)
    allocate (data_old(nx,ny))

    data_old = data

    do ismth = 1, nsmth
       ! Smooth in x-direction
       do jj = 1, ny
          data(1,jj) = 0.25d0 * data_old(nx,jj) + 0.5d0 * data_old(1,jj) + 0.25d0 * data_old(2,jj)
          do ii = 2, nx-1
             data(ii,jj) = 0.25d0 * data_old(ii-1,jj) + 0.5d0 * data_old(ii,jj) + 0.25d0 * data_old(ii+1,jj)
          end do
          data(nx,jj) = 0.25d0 * data_old(nx-1,jj) + 0.5d0 * data_old(nx,jj) + 0.25d0 * data_old(1,jj)
       end do
       data_old = data
       
       ! Smooth in y-direction
       do ii = 1, nx
          i = ii + int(dble(ny)/dble(2))
          if (i > nx) i = mod(i,nx)
          data(ii,1) = 0.25d0*data_old(i,1) + 0.5d0*data_old(i,1) + 0.25d0*data_old(i,2)
          do jj = 2, ny-1
             data(ii,jj) = 0.25d0*data_old(ii,jj-1) + 0.5d0*data_old(ii,jj) + 0.25d0*data_old(ii,jj+1)
          end do
          i = ii + int(dble(ny)/dble(2))
          if (i > nx) i = mod(i,nx)
          data(ii,ny) = 0.25d0*data_old(i,ny-1) + 0.5d0*data_old(i,ny) + 0.25d0*data_old(i,ny)
       end do
       data_old = data
    end do
    
    deallocate (data_old)
  end subroutine smoothing_shapiro

  subroutine wrapij (i, j, nx, ny, i0, j0)
    implicit none
    integer :: i, i0, j, j0, nx, ny

    if (j > ny) then
       j = ny - mod (j, ny)
       i = i + int (dble(nx)/dble(2))
    elseif (j < 1) then 
       j = 1 - j
       i = i + int (dble(ny)/dble(2))
    end if
    
    if (i > nx) then
       i = mod (i, nx)
    elseif (i < 1) then
       i = nx + mod (i, nx)
    end if

    i0 = i
    j0 = j
  end subroutine wrapij

  real(8) function radial_basis_fun (dx, npts, r)
    ! Radial basis function for smoothing topography
    implicit none
    integer :: npts
    real(8) :: dx, r

    real(8) :: alph

    alph = 1d0/(dx_max * dble(npts))

    radial_basis_fun = exp__flush (-(alph*r)**2)
  end function radial_basis_fun
end module utils_mod
