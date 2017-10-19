module viscous_mod
  use shared_mod
  use domain_mod
  use arch_mod
  implicit none

contains
  subroutine divu(dom, i, j, zlev, offs, dims)
    type(Domain) :: dom
    integer :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    
    integer :: id, idS, idW, idSW

    id   = idx(i,     j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)

    dom%divu%elts(id+1) = (&
           velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)     - velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1) &
         + velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)     - velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1) &
         + velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1) - velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1))* &
         dom%areas%elts(id+1)%hex_inv
  end subroutine divu

  subroutine flux_grad_scalar(dom, i, j, zlev, offs, dims)
    ! Diffuse scalars
    type(Domain) :: dom
    integer :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    
    integer :: id, idE, idN, idNE
    
    id = idx(i, j, offs, dims)
    idE = idx(i + 1, j, offs, dims)
    idN = idx(i, j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    ! Mass
    h_mflux(EDGE*id+RT+1) = h_mflux(EDGE*id+RT+1) - viscosity*dom%pedlen%elts(EDGE*id+RT+1)*(mass(idE+1) - mass(id+1)) &
         /dom%len%elts(EDGE*id+RT+1)
    
    h_mflux(EDGE*id+DG+1) = h_mflux(EDGE*id+DG+1) - viscosity*dom%pedlen%elts(EDGE*id+DG+1)*(mass(id+1)  - mass(idNE+1)) &
         /dom%len%elts(EDGE*id+DG+1)
    
    h_mflux(EDGE*id+UP+1) = h_mflux(EDGE*id+UP+1) - viscosity*dom%pedlen%elts(EDGE*id+UP+1)*(mass(idN+1) - mass(id+1)) &
         /dom%len%elts(EDGE*id+UP+1)

    ! Temperature
    h_tflux(EDGE*id+RT+1) = h_tflux(EDGE*id+RT+1) - viscosity*dom%pedlen%elts(EDGE*id+RT+1)*(temp(idE+1) - temp(id+1)) &
         /dom%len%elts(EDGE*id+RT+1)
    
    h_tflux(EDGE*id+DG+1) = h_tflux(EDGE*id+DG+1) - viscosity*dom%pedlen%elts(EDGE*id+DG+1)*(temp(id+1)  - temp(idNE+1)) &
         /dom%len%elts(EDGE*id+DG+1)
    
    h_tflux(EDGE*id+UP+1) = h_tflux(EDGE*id+UP+1) - viscosity*dom%pedlen%elts(EDGE*id+UP+1)*(temp(idN+1) - temp(id+1)) &
         /dom%len%elts(EDGE*id+UP+1)
  end subroutine flux_grad_scalar

 
end module viscous_mod
