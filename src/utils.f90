module utils_mod
  ! Basic functions and routines for evaluating grid-based properties and variables
  
  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, Coord, EDGE, N_BDRY, LORT, UPLT, S_DIVU, S_ROTU, S_MASS, S_TEMP, RT, DG, UP, compressible, &
       z_null, zlevels, min_level, max_level, grav_accel, MATH_PI, C_visc, radius, dt, Laplace_divu, Laplace_rotu, Laplace_sclr, &
       r_dif, level_start, level_end, AT_EDGE, AT_NODE, TRIAG, omega, S_VELO, p_0, p_top, kappa, &
       ref_density, porosity, c_p, a_vert, b_vert, dx_avg, exp__flush
  
  use arch_mod,       only : abort_run, rank
  use comm_mpi_mod,   only : sync_max_real, update_bdry1
  use domain_mod,     only : exner_fun, trend
  use domain_ops_mod, only : apply_onescale, apply_onescale_to_patch, fun3
  use geom_mod,       only : cart2sph, centroid, direction, inner

  use coord_arithmetic_mod
  
  use domain_mod, only : grid, Domain, Float_Field, penal_edge, penal_node, sol, sol_mean, topography, &
       id_edge, idx, velo, velo1, velo2, vort
    
  implicit none
  
  private
  public :: zero_float, zero_float_field
  public :: z_i, zl_i, zl_e, dz_i, dz_e, dz_SW_e, dz_l
  public :: porous_density, porous_density_edge, phi_node, phi_edge
  public :: pressure_i, dA_i, interp, interp_e, interp_latlon_UVW, latlon2uvw
  public :: interp_UVW_latlon, uvw2zonal_merid, vort_triag_to_hex
  public :: hex_pedlen, hex_len, hex_dx, hex_perim, hex2tri, hex2tri2, tri_perim
  public :: equals_float_field, nu_scale, smoothing_rbf, wrapij, radial_basis_fun
  public :: cal_rx0_max, cal_rx1_max

  real(dp)                                  :: rx0_max_loc, rx1_max_loc
  real(dp), dimension(:), pointer           :: val1, val2
  
  interface zero_float
     procedure :: zero_float_0, zero_float_1, zero_float_2
  end interface zero_float


  procedure (fun3), pointer :: hex_fun => null ()
  
contains

  
  function z_i (dom, i, j, zlev, offs, dims, q) result(val)
    ! Position of vertical level zlev at nodes
    ! *** compressible case requires dom%surf_press ***
    
    implicit none
    
    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val

    integer                     :: d, id, l
    real(dp)                    :: dz, exner, rho_dz, rho_dz_theta, p
    real(dp), dimension(0:zlev) :: z

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    z(0) = topography%data(d)%elts(id)
    p    =     dom%surf_press%elts(id)
    do l = 1, zlev
       rho_dz = sol_mean(S_MASS,l)%data(d)%elts(id) + q(S_MASS,l)%data(d)%elts(id) 
       if (compressible) then
          rho_dz_theta = sol_mean(S_TEMP,l)%data(d)%elts(id) + q(S_TEMP,l)%data(d)%elts(id)

          p     = p - grav_accel * rho_dz
          exner = c_p * (p/p_0)**kappa
          
          dz = kappa * rho_dz_theta * exner / p
       else 
          dz = rho_dz / porous_density (d, id, l)
       end if
       z(l) = z(l-1) + dz
    end do

    val = interp (z(zlev-1), z(zlev))
  end function z_i

  
  function zl_i (dom, i, j, zlev, offs, dims, q, pos) result(val)
    ! Position of interface below (l=-1) or above (l=1) vertical level zlev nodes
    ! *** compressible case requires dom%surf_press ***
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, pos, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val
   
    integer  :: d, id, l, lmax
    real(dp) :: dz, exner, rho_dz, rho_dz_theta, p

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    if (pos == -1) then
       lmax = zlev - 1
    else
       lmax = zlev
    end if
    
    val = topography%data(d)%elts(id)
    p    =     dom%surf_press%elts(id)
    do l = 1, lmax
       rho_dz = sol_mean(S_MASS,l)%data(d)%elts(id) + q(S_MASS,l)%data(d)%elts(id) 
       if (compressible) then
          rho_dz_theta = sol_mean(S_TEMP,l)%data(d)%elts(id) + q(S_TEMP,l)%data(d)%elts(id)

          p     = p - grav_accel * rho_dz
          exner = c_p * (p/p_0)**kappa
          
          dz = kappa * rho_dz_theta * exner / p
       else 
          dz = rho_dz / porous_density (d, id, l)
       end if
       val = val + dz
    end do
  end function zl_i

  
  function zl_e (dom, i, j, zlev, offs, dims, q, pos) result(val)
    ! Position of interface below (pos=-1) or above (pos=1) vertical level zlev at edges
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, pos, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val(EDGE)

    real(dp) :: z(0:EDGE)

    z(0)    = zl_i (dom, i,   j,   zlev, offs, dims, q, pos)
    z(RT+1) = zl_i (dom, i+1, j,   zlev, offs, dims, q, pos)
    z(DG+1) = zl_i (dom, i+1, j+1, zlev, offs, dims, q, pos)
    z(UP+1) = zl_i (dom, i,   j+1, zlev, offs, dims, q, pos)

    val = 0.5 * (z(0) + z(1:EDGE))
  end function zl_e

  
  function dz_i (dom, i, j, zlev, offs, dims, q) result(val)
    ! Thickness of layer zlev at nodes
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val

    integer  :: d, id
    real(dp) :: exner, rho_dz, rho_dz_theta, p

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    if (compressible) then ! dz = mu alpha = mu (kappa theta pi) / p = kappa Theta pi / p
       p         = pressure_i (dom, i, j, zlev, offs, dims, sol)
       exner     = c_p * (p/p_0)**kappa
       rho_dz_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id) + q(S_TEMP,zlev)%data(d)%elts(id)

       val = kappa * rho_dz_theta * exner / p
    else                   ! dz = mu/ref_density
       rho_dz = sol_mean(S_MASS,zlev)%data(d)%elts(id) + q(S_MASS,zlev)%data(d)%elts(id)

       val = rho_dz / porous_density (d, id, zlev)
    end if
  end function dz_i

  
  function dz_e (dom, i, j, zlev, offs, dims, q) result(val)
    ! Thickness of layer zlev at edges
    
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1)
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val(EDGE) 

    real(dp) :: dz(0:EDGE)

    dz(0)    = dz_i (dom, i,   j,   zlev, offs, dims, q)
    dz(RT+1) = dz_i (dom, i+1, j,   zlev, offs, dims, q)
    dz(DG+1) = dz_i (dom, i+1, j+1, zlev, offs, dims, q)
    dz(UP+1) = dz_i (dom, i  , j+1, zlev, offs, dims, q)

    val = 0.5 * (dz(0) + dz(1:EDGE))
  end function dz_e

  
  function dz_SW_e (dom, i, j, zlev, offs, dims, q) result(val)
    ! Thickness of layer zlev at edges (SW edges)
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1)
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val(EDGE) 

    integer                     :: d, id, idW, idSW, idS
    real(dp), dimension(0:EDGE) :: dz

    d = dom%id + 1

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    dz(0)    = dz_i (dom, i,   j,   zlev, offs, dims, q)
    dz(RT+1) = dz_i (dom, i-1, j,   zlev, offs, dims, q)
    dz(DG+1) = dz_i (dom, i-1, j-1, zlev, offs, dims, q)
    dz(UP+1) = dz_i (dom, i  , j-1, zlev, offs, dims, q)

    val = 0.5 * (dz(0) + dz(1:EDGE))
  end function dz_SW_e

  
  function dz_l (dom, i, j, zlev, offs, dims, q) result(val)
    ! Thickness of layer associated with interface between layers zlev and zlev+1: z_(zlev+1) - z(zlev)
    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1)
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val
  
    real(dp) :: dz, dz_above

    dz        = dz_i (dom, i, j, zlev,   offs, dims, q)
    dz_above  = dz_i (dom, i, j, zlev+1, offs, dims, q)

    val = 0.5 * (dz + dz_above)
  end function dz_l

  
  function porous_density (d, id_i, zlev) result(val)
    ! Porous density at nodes for incompressible case
    implicit none
    integer, intent(in) :: d, id_i, zlev
    real(dp)            :: val

    val = ref_density * phi_node (d, id_i, zlev)
  end function porous_density

  
  function porous_density_edge (d, id, zlev) result(val)
    ! Porous density at edges
    implicit none
    integer , intent(in) :: d, id, zlev
    real(dp)             :: val(EDGE) 

    val = ref_density * phi_edge (d, id, zlev)
  end function porous_density_edge

  
  function phi_node (d, id_i, zlev) result(val)
    ! Returns porosity at node given by (d, id_i, zlev)
    implicit none
    integer, intent(in) :: d, id_i, zlev
    real(dp)            :: val

    val = 1.0_dp + (porosity - 1.0_dp) * penal_node(zlev)%data(d)%elts(id_i)
  end function phi_node

  
  function phi_edge (d, id, zlev) result(val)
    ! Returns porosity at edges associated to node given by (d, id_i, zlev)
    implicit none
    integer, intent(in) :: d, id, zlev
    real(dp)            :: val(EDGE) 

    integer :: e, id_e
    
    do e = 1, EDGE
       id_e = EDGE*id+e
       val(e) = 1.0_dp + (porosity - 1.0_dp) * penal_edge(zlev)%data(d)%elts(id_e)
    end do
  end function phi_edge

  
  function pressure_i (dom, i, j, zlev, offs, dims, q) result(val)
    ! Pressure at layer zlev computed by integrated down

    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val
 
    integer                             :: d, id, k, l
    real(dp)                            :: rho_dz, rho_dz_theta
    real(dp), dimension(zlev-1:zlevels) :: p

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    p(zlevels) = p_top
    do l = zlevels-1, zlev-1, -1
       k = l + 1 ! layer index
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + q(S_MASS,k)%data(d)%elts(id) 

       ! Pressure at interface l
       if (compressible) then
          p(l) = p(l+1) + grav_accel * rho_dz               
       else 
          rho_dz_theta = sol_mean(S_TEMP,k)%data(d)%elts(id) + q(S_TEMP,k)%data(d)%elts(id)
          p(l) = p(l+1) + grav_accel * (rho_dz - rho_dz_theta) 
       end if
    end do
   val = interp (p(zlev-1), p(zlev)) ! pressure at layer zlev
  end function pressure_i

  
  function dA_i (dom, i, j, zlev, offs, dims) result(val)
    ! For checking areas

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val
   
    val = 1.0_dp
  end function dA_i

  
  function interp (e1, e2) result(val)
    ! Centred average interpolation of quantities e1 and e2
    implicit none
    real(dp), intent(in) :: e1, e2
    real(dp)             :: val

    val = 0.5 * (e1 + e2)
  end function interp

  
  function interp_e (e1, e2) result(val)
    ! Centred average interpolation of edge quantities e1 and e2
    implicit none
    real(dp), intent(in) :: e1(EDGE), e2(EDGE)
    real(dp)             :: val(EDGE)
    
    val = 0.5 * (e1 + e2)
  end function interp_e

  
  subroutine interp_latlon_UVW (dom, i, j, zlev, offs, dims, uvw)
    ! Interpolate from zonal, meridional velocity components at nodes to U, V, W velocity components at edges
    ! (assumes that dom%u_zonal and dom%v_merid have been set over all grid points)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: uvw(1:EDGE)

    integer     :: id, idE, idN, idNE
    type(Coord) :: vel0

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    vel0 = vel (id)

    uvw(RT+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idE+1)),  0.5_dp * (vel0 + vel (idE)))
    uvw(DG+1) = inner (direction (dom%node%elts(idNE+1), dom%node%elts(id+1)),   0.5_dp * (vel0 + vel (idNE)))
    uvw(UP+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idN+1)),  0.5_dp * (vel0 + vel (idN)))
  contains
    function vel (id) result(val)
      ! Computes velocity at node id from its latitude and longitude components
      integer, intent(in) :: id
      type(Coord)         :: val

      real(dp)    :: lat, lon
      type(Coord) :: e_merid, e_zonal

      call cart2sph (dom%node%elts(id+1), lon, lat)

      e_zonal = Coord (-sin(lon),             cos(lon),                 0.0_dp) 
      e_merid = Coord (-cos(lon) * sin(lat), -sin(lon) * sin(lat), cos(lat))

      val = dom%u_zonal%elts(id+1) * e_zonal + dom%v_merid%elts(id+1) * e_merid
    end function vel
  end subroutine interp_latlon_UVW

  
  function latlon2uvw (dom, i, j, zlev, offs, dims) result(val)
    ! Interpolate from zonal, meridional velocity components at nodes to U, V, W velocity components at edges
    ! (assumes that dom%u_zonal and dom%v_merid have been set over all grid points)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val(1:EDGE)

    integer     :: id, idE, idN, idNE
    type(Coord) :: vel0

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    vel0 = vel (id)

    val(RT+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idE+1)),  0.5_dp * (vel0 + vel (idE)))
    val(DG+1) = inner (direction (dom%node%elts(idNE+1), dom%node%elts(id+1)),   0.5_dp * (vel0 + vel (idNE)))
    val(UP+1) = inner (direction (dom%node%elts(id+1),   dom%node%elts(idN+1)),  0.5_dp * (vel0 + vel (idN)))
  contains
    function vel (id) result(val)
      ! Computes velocity at node id from its latitude and longitude components
      integer, intent(in) :: id
      type(Coord)         :: val

      real(dp) :: lat, lon
      
      type(Coord) :: e_merid, e_zonal

      call cart2sph (dom%node%elts(id+1), lon, lat)

      e_zonal = Coord (-sin(lon),             cos(lon),                 0.0_dp) 
      e_merid = Coord (-cos(lon) * sin(lat), -sin(lon) * sin(lat), cos(lat))

      val = dom%u_zonal%elts(id+1) * e_zonal + dom%v_merid%elts(id+1) * e_merid
    end function vel
  end function latlon2uvw

  
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
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer     :: id, idN, idE, idNE, idS, idSW, idW
    real(dp)    :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW
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
         [ &
         dom%ccentre%elts(TRIAG*id+LORT+1),   dom%ccentre%elts(TRIAG*id+UPLT+1),   &
         dom%ccentre%elts(TRIAG*idW+LORT+1),  dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
         dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS+UPLT+1) ], 6)

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
    e_zonal = Coord (-sin(lon),           cos(lon),               0.0_dp) 
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))

    ! Project velocity at node onto zonal and meridional directions
    velo1(id+1) = inner (vel, e_zonal)
    velo2(id+1) = inner (vel, e_merid)
  end subroutine interp_UVW_latlon

  
  function uvw2zonal_merid (dom, i, j, zlev, offs, dims) result(val)
    ! Interpolate velocity from U, V, W velocity components at edges to zonal, meridional velocity components at nodes
    ! uses sol(S_VELO,zlev)
    ! Perot reconstruction based on Gauss theorem:
    !
    ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_centroid) ) / cell_area
    !
    ! also used for kinetic energy
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
       
    real(dp)                 :: val(2) 

    integer     :: d, id, idN, idE, idNE, idS, idSW, idW
    real(dp)    :: lon, lat, u_dual_RT, u_dual_UP, u_dual_DG, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW
    type(Coord) :: cent, e_zonal, e_merid, vel

    d = dom%id + 1
    
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    ! Fluxes normal to hexagon edges
    u_dual_RT    =  sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1)   * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG    = -sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1)   * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    =  sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1)   * dom%pedlen%elts(EDGE*id+UP+1)

    u_dual_RT_W  = -sol(S_VELO,zlev)%data(d)%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
    u_dual_DG_SW =  sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
    u_dual_UP_S  = -sol(S_VELO,zlev)%data(d)%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)

    ! Compute hexagon centroid from its vertices
    cent = centroid (                                                                 &
         [ dom%ccentre%elts(TRIAG*id+LORT+1),   dom%ccentre%elts(TRIAG*id+UPLT+1),   &
            dom%ccentre%elts(TRIAG*idW+LORT+1),  dom%ccentre%elts(TRIAG*idSW+UPLT+1), &
            dom%ccentre%elts(TRIAG*idSW+LORT+1), dom%ccentre%elts(TRIAG*idS+UPLT+1) ], 6)

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
    e_zonal = Coord (-sin(lon),           cos(lon),               0.0_dp) 
    e_merid = Coord (-cos(lon)*sin(lat), -sin(lon)*sin(lat), cos(lat))

    ! Project velocity at node onto zonal and meridional directions
    val(1) = inner (vel, e_zonal)
    val(2) = inner (vel, e_merid)
  end function uvw2zonal_merid

  
  subroutine vort_triag_to_hex (dom, i, j, zlev, offs, dims)
    ! Approximate vorticity at hexagon points
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

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

  
  function hex2tri (dom, i, j, t, offs, dims, zlev) result(val)
    ! Float array hex_fun at triangles associated with node (i,j) computed from integral over hexagons
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, t, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: id, idE, idNE, idN

    id   = idx (i,   j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    select case (t)
    case (LORT)
       idE = idx (i+1, j, offs, dims)
       val = &
            hex_fun (dom, i,   j,   zlev, offs, dims) * dom%areas%elts(id+1)%part(1)  + &
            hex_fun (dom, i+1, j+1, zlev, offs, dims) * dom%areas%elts(idNE+1)%part(5) + &
            hex_fun (dom, i+1, j,   zlev, offs, dims) * dom%areas%elts(idE+1)%part(3)
    case (UPLT)
       idN = idx (i, j+1, offs, dims)
       val = &
            hex_fun (dom, i,   j,   zlev, offs, dims) * dom%areas%elts(id+1)%part(2)  + &
            hex_fun (dom, i+1, j+1, zlev, offs, dims) * dom%areas%elts(idNE+1)%part(4) + &
            hex_fun (dom, i,   j+1, zlev, offs, dims) * dom%areas%elts(idN+1)%part(6)
    case default
       val = 0.0_dp
    end select

    val = val / dom%triarea%elts(TRIAG*id+t+1)
  end function hex2tri

  
  function hex2tri2 (sclr, hex_area, tri_area, t) result(val)
    ! Interpolates sclr given at hexagons to triangles
    
    implicit none
    
    integer, intent(in) :: t
    real(4), intent(in) :: tri_area
    real(4), intent(in) :: sclr(0:EDGE)
    real(4), intent(in) :: hex_area(2*EDGE)
    real(4)             :: val

    if (t == LORT) then
       val = (sclr(0) * hex_area(1) + sclr(1) * hex_area(3) + sclr(2) * hex_area(5)) / tri_area
    else
       val = (sclr(0) * hex_area(2) + sclr(2) * hex_area(4) + sclr(3) * hex_area(6)) / tri_area
    end if
  end function hex2tri2

  
  subroutine zero_float_0 (q)
    ! Initializes a float field to zero
    implicit none
    type(Float_Field), intent(inout) :: q

    integer :: d

    do d = 1, size(grid)
       q%data(d)%elts = 0.0_dp
    end do
  end subroutine zero_float_0

  
  subroutine zero_float_1 (q)
    ! Initializes a float vector to zero
    implicit none
    type(Float_Field), intent(inout) :: q(:)

    integer :: d, j

    do j = 1, size(q,1)
       do d = 1, size(grid)
          q(j)%data(d)%elts = 0.0_dp
       end do
    end do
  end subroutine zero_float_1

  
  subroutine zero_float_2 (q)
    ! Initializes a float array to zero
    implicit none
    type(Float_Field), intent(inout) :: q(:,:)

    integer :: d, j1, j2

    do j1 = 1, size(q,1)
       do j2 = 1, size(q,2)
          do d = 1, size(grid)
             q(j1,j2)%data(d)%elts = 0.0_dp
          end do
       end do
    end do
  end subroutine zero_float_2

  
  subroutine zero_float_field (q, itype, lmin_in, lmax_in)
    ! Set float field to zero for scales:
    ! (lmin,lmax) if both lmin and lmax are present
    ! lmin if only lmin is present
    ! (level_start, level_end) if lmin is not present
    ! itype = S_MASS or S_VELO
    implicit none
    integer,           intent(in)              :: itype
    integer,           intent(in),    optional :: lmin_in, lmax_in
    type(Float_Field), intent(inout), target   :: q

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
       call abort_run
    end if

    q%bdry_uptodate = .false.
    call update_bdry1 (q, lmin, lmax, 909)
  end subroutine zero_float_field

  
  subroutine cal_zero_node (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
   
    integer :: id

    id = idx (i, j, offs, dims) + 1

    val1(id) = 0.0_dp
  end subroutine cal_zero_node

  
  subroutine cal_zero_edge (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer :: id

    id = idx (i, j, offs, dims)

    val1(id_edge(id)) = 0.0_dp
  end subroutine cal_zero_edge

  
  subroutine equals_float_field (q1, q2, itype, lmin_in, lmax_in)
    ! Set elements of float field q1 = q2
    !
    ! itype = S_MASS or S_VELO
    ! if scale l is present, compute only for scale l
    implicit none
    integer,           intent(in)              :: itype
    integer,           intent(in),    optional :: lmin_in, lmax_in
    type(Float_Field), intent(in),    target   :: q2
    type(Float_Field), intent(inout), target   :: q1

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
       call abort_run
    end if

    q1%bdry_uptodate = .false.
    call update_bdry1 (q1, lmin, lmax, 910)
  end subroutine equals_float_field

  
  subroutine cal_equals_node (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    val1(id) = val2(id)
  end subroutine cal_equals_node

  
  subroutine cal_equals_edge (dom, i, j, zlev, offs, dims)
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer :: id

    id = idx (i, j, offs, dims) 
    
    val1(id_edge(id)) = val2(id_edge(id))
  end subroutine cal_equals_edge

  
  function nu_scale (type, zlev, dx_loc) result(val)
    ! Viscosity for diffusion on hexagonal C-grids
    ! include dom, id for scale-aware viscosity
    ! (conservative Gershgorin estimate to be consisistent with pentagons and irregular grid
    ! uses 1/8 factor instead of 1/6 for regular hexagonal grid)
    implicit none
    integer,  intent(in)           :: type, zlev
    real(dp), intent(in), optional :: dx_loc
    real(dp)                       :: val

    real(dp) :: dx

    ! Correction factors for irregular grid
    ! (determined by experiment to satisfy local stability on adaptive grids)
    real(dp), parameter :: rho_sclr = 1.10_dp
    real(dp), parameter :: rho_divu = 1.10_dp 
    real(dp), parameter :: rho_rotu = 1.65_dp 

    if (present(dx_loc)) then
       dx = dx_loc
    else
       dx = dx_avg(max_level)
    end if
    
    select case (type)
    case (S_MASS)
       val = C_visc(S_MASS,zlev) * r_dif/dt * (dx**2/8  / rho_sclr)**Laplace_sclr
    case (S_TEMP)
       val = C_visc(S_TEMP,zlev) * r_dif/dt * (dx**2/8  / rho_sclr)**Laplace_sclr
    case (S_DIVU)
       val = C_visc(S_DIVU,zlev) * r_dif/dt * (dx**2/8  / rho_divu)**Laplace_divu
    case (S_ROTU)
       val = C_visc(S_ROTU,zlev) * r_dif/dt * (dx**2/24 / rho_rotu)**Laplace_rotu
    case default
       val = 0.0_dp
    end select
  end function nu_scale

  
  subroutine smoothing_rbf (dx, npts, nsmth, data)
    ! Smooths data(lon,lat) over neighbouring region using radial basis functions
    implicit none
    integer,  intent(in)    :: npts, nsmth
    real(dp), intent(in)    :: dx
    real(dp), intent(inout) :: data(:,:)

    integer                               :: i, ismth, i0, ii, iw, j, j0, jj, jw, nx, ny
    real(dp)                              :: r, M_topo, sw_topo, topo_sum, wgt
    real(dp), allocatable, dimension(:,:) :: data_old

    nx = size(data,1); ny = size(data,2)
    allocate (data_old(nx,ny))

    data_old = data

    do ismth = 1, nsmth
       do i = 1, nx
          do j = 1, ny
             sw_topo  = 0.0_dp
             topo_sum = 0.0_dp
             do ii = -npts, npts
                do jj = -npts, npts

                   r = sqrt (dble(ii**2 + jj**2)) * MATH_PI * radius / dble (ny)

                   wgt = radial_basis_fun (dx, npts, r)

                   iw = i + ii
                   jw = j + jj
                   call wrapij (iw, jw, nx, ny, i0, j0)

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

  
  subroutine wrapij (i, j, nx, ny, i0, j0)
    implicit none
    integer, intent(in)    :: nx, ny
    integer, intent(inout) :: i, j
    integer, intent(out)   :: i0, j0

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

  
  function radial_basis_fun (dx, npts, r) result(val)
    ! Radial basis function for smoothing topography
    implicit none
    integer, intent(in) :: npts
    real(dp),intent(in) :: dx, r
    real(dp)            :: val

    real(dp) :: alph

    alph = 1 / (dx_avg(min_level) * dble(npts))

    val = exp__flush (-(alph*r)**2)
  end function radial_basis_fun

  
  subroutine cal_rx0_max (l, rx0_max)
    ! Maximum topographic stiffness ratio rx0 < 0.2
    ! (also called the Beckman-Haidvogel number)
    ! the compressible version uses surface pressure instead of topography height
    
    implicit none
    
    integer,  intent(in)  :: l
    real(dp), intent(out) :: rx0_max

    rx0_max_loc = 0.0_dp
    rx0_max     = 0.0_dp

    if (compressible) then
       call apply_onescale (cal_rx0_loc_P, l, z_null, 0, 1)
    else
       call apply_onescale (cal_rx0_loc_Z, l, z_null, 0, 1)
    end if

    rx0_max = sync_max_real (rx0_max_loc)
  end subroutine cal_rx0_max

  
  subroutine cal_rx1_max (l, rx1_max)
    ! Hydrostatic maximum instability number rx1 < 1 (also called Haney number)
    ! (Haney 1991, Shchepetkin and McWilliams 2003)
    ! note that rx1 < 1 is almost impossible to achieve and rx1 <= 5 is usually okay in oceanographic simulations
    ! compute only over lowest layer (most unstable)
    
    implicit none
    
    integer,  intent(in)  :: l
    real(dp), intent(out) :: rx1_max

    integer :: k

    rx1_max_loc = 0.0_dp
    rx1_max     = 0.0_dp

    do k = 1, zlevels
       if (compressible) then
          call apply_onescale (cal_rx1_loc_P, l, k, 0, 1)
       else
          call apply_onescale (cal_rx1_loc_Z, l, k, 0, 1)
       end if
    end do

    rx1_max = sync_max_real (rx1_max_loc)
  end subroutine cal_rx1_max
  

  subroutine cal_rx0_loc_P (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    ! uses surface pressure
    
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer  :: id, idE, idNE, idN
    real(dp) :: h0, h1

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       h0 = dom%surf_press%elts(id+1)

       h1 = dom%surf_press%elts(idE+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())

       h1 = dom%surf_press%elts(idNE+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())

       h1 = dom%surf_press%elts(idN+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())
    end if
    
  contains

    function rx0 () result(val)
      implicit none
      real(dp) :: val

      val = abs (h0 - h1) / (h0 + h1)
    end function rx0
    
  end subroutine cal_rx0_loc_P


  subroutine cal_rx1_loc_P (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    ! uses pressure vertical coordinates

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    
    integer  :: id, idE, idNE, idN
    real(dp) :: z1, z2, z3, z4

    id  = idx (i, j, offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       z1 = p(id, zlev+1)
       z2 = p(id, zlev)

       z3 = p(idE, zlev+1)
       z4 = p(idE, zlev)
       rx1_max_loc = max (rx1_max_loc, rx1 ())

       z3 = p(idNE, zlev+1)
       z4 = p(idNE, zlev)
       rx1_max_loc = max (rx1_max_loc, rx1 ())

       z3 = p(idN, zlev+1)
       z4 = p(idN, zlev)
       rx1_max_loc = max (rx1_max_loc, rx1 ())
    end if
    
  contains
    
   function p (id, k) result(val)
      implicit none
      integer, intent(in) :: id, k
      real(dp)            :: val

      val = a_vert(k) + b_vert(k) * dom%surf_press%elts(id+1)
    end function p

    function rx1 () result(val)
      implicit none
      real(dp) :: val

      val = abs ( (z4 - z2 + z3 - z1) / (z4 + z2 - z3 - z1) )
    end function rx1
    
  end subroutine cal_rx1_loc_P

  
  subroutine cal_rx0_loc_Z (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    ! uses z vertical coordinates
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer  :: d, id, idE, idN, idNE
    real(dp) :: h0, h1

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    d    = dom%id + 1

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       h0 = topography%data(d)%elts(id+1)

       h1 = topography%data(d)%elts(idE+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())

       h1 = topography%data(d)%elts(idNE+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())

       h1 = topography%data(d)%elts(idN+1)
       rx0_max_loc = max (rx0_max_loc, rx0 ())
    end if
    
  contains
    
    function rx0 () result(val)
      implicit none
      real(dp) :: val

      val = abs (h0 - h1) / (h0 + h1)
    end function rx0
    
  end subroutine cal_rx0_loc_Z

  
  subroutine cal_rx1_loc_Z (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    ! uses z vertical coordinates
    
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer  :: id
    real(dp) :: z1, z2, z3, z4

    id  = idx (i, j, offs, dims)

    !if (dom%mask_n%elts(id+1) >= ADJZONE) then
    z1 = zl_i (dom, i, j, zlev, offs, dims, sol, -1)
    z2 = zl_i (dom, i, j, zlev, offs, dims, sol,  1)

    z3 = zl_i (dom, i+1, j, zlev, offs, dims, sol, -1)
    z4 = zl_i (dom, i+1, j, zlev, offs, dims, sol,  1)
    rx1_max_loc = max (rx1_max_loc, rx1 ())

    z3 = zl_i (dom, i+1, j+1, zlev, offs, dims, sol, -1)
    z4 = zl_i (dom, i+1, j+1, zlev, offs, dims, sol,  1)
    rx1_max_loc = max (rx1_max_loc, rx1 ())

    z3 = zl_i (dom, i, j+1, zlev, offs, dims, sol, -1)
    z4 = zl_i (dom, i, j+1, zlev, offs, dims, sol,  1)
    rx1_max_loc = max (rx1_max_loc, rx1 ())
    !end if
  contains
    function rx1 () result(val)
      implicit none
      real(dp) :: val

      val = abs ( (z4 - z2 + z3 - z1) / (z4 + z2 - z3 - z1) )
    end function rx1
  end subroutine cal_rx1_loc_Z


  function hex_perim (dom, i, j, offs, dims) result(val)
    ! Perimeter of hexagon associated to i, j
    use domain_mod
    implicit none
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    type(Domain)                   :: dom
    real(dp)                       :: val

    integer :: id, idW, idSW, idS
    
    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    val = &
         dom%pedlen%elts(EDGE*id +RT+1) + dom%len%elts(EDGE*id  +DG+1) + dom%len%elts(EDGE*id +UP+1) + &
         dom%pedlen%elts(EDGE*idW+RT+1) + dom%len%elts(EDGE*idSW+DG+1) + dom%len%elts(EDGE*idS+UP+1)
  end function hex_perim

  
  function tri_perim (dom, i, j, t, offs, dims) result(val)
    ! Perimeter of triangles associated to i, j
    implicit none
    integer                        :: i, j, t
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    type(Domain)                   :: dom
    real(dp)                       :: val

    integer :: id, idE, idN

    id  = idx (i,   j,   offs, dims)
    idE = idx (i+1, j,   offs, dims)
    idN = idx (i,   j+1, offs, dims)

    select case (t)
    case (LORT) 
       val = dom%len%elts(EDGE*id+RT+1) + dom%len%elts(EDGE*idE+UP+1) + dom%len%elts(EDGE*id+DG+1) 
    case (UPLT) 
       val = dom%len%elts(EDGE*id+DG+1) + dom%len%elts(EDGE*id +UP+1) + dom%len%elts(EDGE*idN+RT+1)
    case default
       val = 0.0_dp
    end select
  end function tri_perim

  
  function hex_pedlen (dom, i, j, offs, dims) result(val)
    ! The six primary grid edges (hexagon edges) associated to hexagon i, j
    
    implicit none
    
    integer,      intent(in) :: i, j
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    type(Domain), intent(in) :: dom
    real(dp)                 :: val(1:2*EDGE)

    integer                       :: id, idW, idSW, idS
    integer,  dimension(1:2*EDGE) :: ide
    
    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    ide = [ id_edge(id), EDGE*idW + RT + 1, EDGE*idSW + DG + 1, EDGE*idS + UP + 1 ]

    val = dom%pedlen%elts(ide)
  end function hex_pedlen

  
  function hex_len (dom, i, j, offs, dims) result(val)
    ! The six dual grid edges (distances to neighbour hexagons) associated to hexagon i, j
    
    implicit none
    
    integer,      intent(in) :: i, j
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    type(Domain), intent(in) :: dom
    real(dp)                 :: val(1:2*EDGE)
   

    integer :: id, idW, idSW, idS
    integer,  dimension(1:2*EDGE) :: ide

    id   = idx (i, j, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    
    ide = [ id_edge(id), EDGE*idW + RT + 1, EDGE*idSW + DG + 1, EDGE*idS + UP + 1 ]

    val = dom%len%elts(ide)
  end function hex_len

  
  function hex_dx (dom, i, j, offs, dims) result(val)
    ! Equivalent dx for a hexagon

    implicit none

    integer,      intent(in) :: i, j
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    type(Domain), intent(in) :: dom
    real(dp)                 :: val

    integer  :: id
    real(dp) :: Area, Perimeter
    
    id = idx (i, j, offs, dims)

    Perimeter = sum (hex_pedlen (dom, i, j, offs, dims))
    Area      = 1 / dom%areas%elts(id+1)%hex_inv
    
    val = 4 * Area / Perimeter
  end function hex_dx
  

end module utils_mod
