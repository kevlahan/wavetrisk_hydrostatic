module viscous_mod
  use shared_mod
  use domain_mod
  use arch_mod
  use ops_mod
  implicit none

contains

  subroutine cal_divu (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW

    id   = idx(i,     j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)

    dom%divu%elts(id+1) = ( &
           velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)     - velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1)  &
         + velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)     - velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1)  &
         + velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1) - velo(EDGE*id +DG+1)*dom%pedlen%elts(EDGE*id +DG+1)) &
         * dom%areas%elts(id+1)%hex_inv
  end subroutine cal_divu

  subroutine flux_grad_scalar (dom, i, j, zlev, offs, dims)
    ! Diffuse scalars
    type(Domain) :: dom
    integer :: i, j, zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims

    integer :: id, idE, idN, idNE

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    ! Temperature
    h_tflux(EDGE*id+RT+1) = h_tflux(EDGE*id+RT+1) &
         - viscosity_temp * dom%pedlen%elts(EDGE*id+RT+1) * (temp(idE+1) - temp(id+1))/dom%len%elts(id*EDGE+RT+1)  

    h_tflux(EDGE*id+DG+1) = h_tflux(EDGE*id+DG+1) &
         - viscosity_temp * dom%pedlen%elts(EDGE*id+DG+1) * (temp(id+1)  - temp(idNE+1))/dom%len%elts(id*EDGE+DG+1)

    h_tflux(EDGE*id+UP+1) = h_tflux(EDGE*id+UP+1) &
         - viscosity_temp * dom%pedlen%elts(EDGE*id+UP+1) * (temp(idN+1) - temp(id+1))/dom%len%elts(id*EDGE+UP+1) 
  end subroutine flux_grad_scalar

  subroutine diffuse_momentum (dom, i, j, zlev, offs, dims)
    ! Calculate grad(divu) + grad(vort) to diffuse momentum
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idN, idW, idE, idNE
    real(8), dimension(3) :: grad_div

    id   = idx(i,     j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if(dom%pedlen%elts(id*EDGE+RT+1).eq.0.0_8) write(6,*) 'RT'
    if(dom%pedlen%elts(id*EDGE+DG+1).eq.0.0_8) write(6,*) 'DG'
    if(dom%pedlen%elts(id*EDGE+UP+1).eq.0.0_8) write(6,*) 'UP'

    grad_div = grad(divu)
    
    dvelo(id*EDGE+RT+1) = dvelo(id*EDGE+RT+1) &
         + viscosity_velo * (&
           (vort(id*TRIAG+LORT+1)-vort(idS*TRIAG+UPLT+1)) &
         + (divu(idE+1)-divu(id+1))/dom%len%elts(id*EDGE+RT+1))
        ! + (divu(idE+1)-divu(id+1))*dom%len%elts(id*EDGE+RT+1)/dom%pedlen%elts(id*EDGE+RT+1))

    dvelo(id*EDGE+DG+1) = dvelo(id*EDGE+DG+1) &
         + viscosity_velo * (&
          (vort(id*TRIAG+LORT+1)-vort(id*TRIAG+UPLT+1)) &
          + (divu(id+1)-divu(idNE+1))/dom%len%elts(id*EDGE+DG+1))
        ! + (divu(id+1)-divu(idNE+1))*dom%len%elts(id*EDGE+DG+1)/dom%pedlen%elts(id*EDGE+DG+1))

    dvelo(id*EDGE+UP+1) = dvelo(id*EDGE+UP+1) &
         + viscosity_velo * (&
         (vort(idW*TRIAG+LORT+1)-vort(id*TRIAG+UPLT+1)) &
          + (divu(idN+1)-divu(id+1))/dom%len%elts(id*EDGE+UP+1))
        ! + (divu(idN+1)-divu(id+1))*dom%len%elts(id*EDGE+UP+1)/dom%pedlen%elts(id*EDGE+UP+1))
  end subroutine diffuse_momentum
end module viscous_mod
