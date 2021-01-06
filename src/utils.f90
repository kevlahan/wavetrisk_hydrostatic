module utils_mod
  ! Basic functions used primarily by test_case_module
  use domain_mod
  implicit none
contains
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

    dz_i = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1) + sol(S_MASS,zlev)%data(d)%elts(id+1)) / porous_density (d, id, zlev)
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
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    dz(0)    = (sol_mean(S_MASS,zlev)%data(d)%elts(id+1)   + sol(S_MASS,zlev)%data(d)%elts(id+1))  /porous_density (d, id,   zlev)
    dz(RT+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)  + sol(S_MASS,zlev)%data(d)%elts(idE+1)) /porous_density (d, idE,  zlev)
    dz(DG+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1) + sol(S_MASS,zlev)%data(d)%elts(idNE+1))/porous_density (d, idNE, zlev)
    dz(UP+1) = (sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)  + sol(S_MASS,zlev)%data(d)%elts(idN+1)) /porous_density (d, idN,  zlev)

    dz_e = 0.5 * (dz(0) + dz(1:EDGE))
  end function dz_e
  
  real(8) function z_i (dom, i, j, zlev, offs, dims)
    ! Position of upper interface of level zlev  at nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, k

    id = idx (i, j, offs, dims)

    z_i = dom%topo%elts(id+1)
    do k = 1, zlev
       z_i = z_i + dz_i (dom, i, j, k, offs, dims)
    end do
  end function z_i

  function z_e (dom, i, j, zlev, offs, dims)
    ! Position of upper interface of level zlev at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: z_e

    integer :: id, idE, idN, idNE, k

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    z_e(RT+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idE+1))  
    z_e(DG+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idNE+1)) 
    z_e(UP+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idN+1))  
    do k = 1, zlev
       z_e = z_e + dz_e (dom, i, j, k, offs, dims)
    end do
  end function z_e

  function eta_e (dom, i, j, zlev, offs, dims)
    ! Free surface perturbation at edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: eta_e

    integer :: id, idE, idN, idNE

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    eta_e(RT+1) = 0.5 * (scalar(id+1) + scalar(idE+1))
    eta_e(DG+1) = 0.5 * (scalar(id+1) + scalar(idNE+1))
    eta_e(UP+1) = 0.5 * (scalar(id+1) + scalar(idN+1))
  end function eta_e

  real(8) function porous_density (d, id, k)
    ! Porous density at nodes
    implicit none
    integer :: d, id, k
    
    porous_density = ref_density * (1.0_8 + (alpha - 1.0_8) * penal_node(k)%data(d)%elts(id+1))
  end function porous_density
end module utils_mod
