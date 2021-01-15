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

    zl_e(RT+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idE+1))  
    zl_e(DG+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idNE+1)) 
    zl_e(UP+1) = 0.5 * (dom%topo%elts(id+1) + dom%topo%elts(idN+1))  
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

    d = dom%id + 1
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    eta_e(RT+1) = 0.5 * (sol(S_MASS,zlevels+1)%data(d)%elts(id+1) + sol(S_MASS,zlevels+1)%data(d)%elts(idE+1))
    eta_e(DG+1) = 0.5 * (sol(S_MASS,zlevels+1)%data(d)%elts(id+1) + sol(S_MASS,zlevels+1)%data(d)%elts(idNE+1))
    eta_e(UP+1) = 0.5 * (sol(S_MASS,zlevels+1)%data(d)%elts(id+1) + sol(S_MASS,zlevels+1)%data(d)%elts(idN+1))
  end function eta_e

  real(8) function richardson (dom, i, j, zlev, offs, dims, l)
    ! Richardson number at interface below (l=-1) or above (l=1) level zlev
    ! (defined at nodes)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id, id_i
    real(8)                    :: drho, dz_l, mass_0, mass_l, rho_0, rho_l, temp_0, temp_l, u_0, u_l, v_mag
    real(8), dimension(1:EDGE) :: du_e

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    mass_0 = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)
    temp_0 = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)
    rho_0 = (1.0_8 - temp_0 / mass_0) * porous_density (dom, i, j, zlev, offs, dims)

    mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
    temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
    rho_l = (1.0_8 - temp_l / mass_l) * porous_density (dom, i, j, zlev+l, offs, dims)

    drho = rho_l - rho_0

    du_e = sol(S_VELO,zlev+l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    v_mag = sum(du_e**2)/3
    
    dz_l = 0.5 * (dz_i(dom, i, j, zlev, offs, dims) + dz_i(dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

    if (drho == 0.0_8) then
       richardson = 0.0_8
    elseif (v_mag == 0.0_8) then
       richardson = 1d10
    else
       richardson = - grav_accel * l * drho/ref_density * dz_l / v_mag
    end if
  end function richardson

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
end module utils_mod
