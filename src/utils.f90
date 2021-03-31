module utils_mod
  ! Basic functions used primarily by test_case_module
  use domain_mod
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
    id = idx (i, j, offs, dims)

    dz_i = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (dom, i, j, zlev, offs, dims)
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
         / porous_density (dom, i, j, zlev, offs, dims)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)  + sol(S_MASS,zlev)%data(d)%elts(idE+1)) &
         / porous_density (dom, i+1, j, zlev, offs, dims)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1) + sol(S_MASS,zlev)%data(d)%elts(idNE+1)) &
         / porous_density (dom, i+1, j+1, zlev, offs, dims)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)  + sol(S_MASS,zlev)%data(d)%elts(idN+1)) &
         / porous_density (dom, i, j+1, zlev, offs, dims)

    dz_e = 0.5 * (dz(0) + dz(1:EDGE))
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

    id   = idx (i, j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i,   j-1, offs, dims)
    idS  = idx (i-1, j-1, offs, dims)
   
    dz(0) = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1)   + sol(S_MASS,zlev)%data(d)%elts(id+1)) &
         / porous_density (dom, i, j, zlev, offs, dims)
    
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idW+1)  + sol(S_MASS,zlev)%data(d)%elts(idW+1)) &
         / porous_density (dom, i-1, j, zlev, offs, dims)
    
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idSW+1) + sol(S_MASS,zlev)%data(d)%elts(idSW+1)) &
         / porous_density (dom, i-1, j-1, zlev, offs, dims)
    
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idS+1)  + sol(S_MASS,zlev)%data(d)%elts(idS+1)) &
         / porous_density (dom, i, j-1, zlev, offs, dims)

    dz_SW_e = 0.5 * (dz(0) + dz(1:EDGE))
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

  real(8) function porous_density (dom, i, j, zlev, offs, dims)
    ! Porous density at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    
    d = dom%id + 1
    id = idx (i, j, offs, dims)
    
    porous_density = ref_density * (1.0_8 + (alpha - 1.0_8) * penal_node(zlev)%data(d)%elts(id+1))
  end function porous_density

  real(8) function phi_node (d, id_i, zlev)
    ! Returns porosity at node given by (d, id_i, zlev)
    implicit none
    integer :: d, id_i, zlev

    phi_node = 1.0_8 + (alpha - 1.0_8) * penal_node(zlev)%data(d)%elts(id_i)
  end function phi_node

  real(8) function phi_edge (d, id_e, zlev)
    ! Returns porosity at edge given by (d, id_e, zlev)
    implicit none
    integer :: d, id_e, zlev

    phi_edge = 1.0_8 + (alpha - 1.0_8) * penal_edge(zlev)%data(d)%elts(id_e)
  end function phi_edge

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
    
    total_depth = 0.0_8
    do k = 1, zlevels
       full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
       total_depth = total_depth + full_mass /  porous_density (dom, i, j, k, offs, dims)
    end do
    free_surface = total_depth + dom%topo%elts(id_i)
  end function free_surface

  real(8) function interp (e1, e2)
    ! Centred average interpolation of quantities e1 and e2
    implicit none
    real(8) :: e1, e2

    interp = 0.5 * (e1 + e2)
  end function interp

  function interp_e (e1, e2)
    ! Centred average interpolation of edge quantities e1 and e2
    implicit none
    real(8), dimension(1:EDGE) :: interp_e
    real(8), dimension(1:EDGE) :: e1, e2

    interp_e = 0.5 * (e1 + e2)
  end function interp_e
end module utils_mod
