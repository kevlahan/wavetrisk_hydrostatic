module viscous_mod
  use shared_mod
  use domain_mod
  use arch_mod
  implicit none

contains
  subroutine divu(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idS
    integer idW
    integer idSW

    id   = idx(i,     j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)

    dom%divu%elts(id+1) = (&
         velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1)     - velo(EDGE*id+DG+1)*dom%pedlen%elts(EDGE*id+DG+1) &
         + velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1)     - velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1) &
         + velo(EDGE*idSW+DG+1)*dom%pedlen%elts(EDGE*idSW+DG+1) - velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1))*&
         dom%areas%elts(id+1)%hex_inv
  end subroutine divu

  subroutine flux_gradP(dom, i, j, offs, dims)
    type(Domain) dom
    integer i, j
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idE
    integer idN
    integer idNE
    id = idx(i, j, offs, dims)
    idE = idx(i + 1, j, offs, dims)
    idN = idx(i, j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    h_mflux(EDGE*id+RT+1) = h_mflux(EDGE*id+RT+1) - &
         viscosity*dom%pedlen%elts(EDGE*id+RT+1)* &
         (mass(idE+1) - mass(id+1))/dom%len%elts(EDGE*id+RT+1)
    h_mflux(EDGE*id+DG+1) = h_mflux(EDGE*id+DG+1) - &
         viscosity*dom%pedlen%elts(EDGE*id+DG+1)* &
         (mass(id+1) - mass(idNE+1))/dom%len%elts(EDGE*id+DG+1)
    h_mflux(EDGE*id+UP+1) = h_mflux(EDGE*id+UP+1) - &
         viscosity*dom%pedlen%elts(EDGE*id+UP+1)* &
         (mass(idN+1) - mass(id+1))/dom%len%elts(EDGE*id+UP+1)
  end subroutine flux_gradP

  subroutine diff_mom(dom, i, j, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY + 1) :: dims
    integer id
    integer idS, idN
    integer idW, idE
    integer idNE
    id = idx(i, j, offs, dims)
    idS = idx(i, j - 1, offs, dims)
    idW = idx(i - 1, j, offs, dims)
    idN = idx(i, j + 1, offs, dims)
    idE = idx(i + 1, j, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)
    dvelo(id*EDGE+RT+1) = dvelo(id*EDGE+RT+1) + viscosity*( &
         (dom%divu%elts(idE+1) - dom%divu%elts(id+1)) + &
         (dom%vort%elts(id*TRIAG+LORT+1) - dom%vort%elts(idS*TRIAG+UPLT+1)) &
         /dom%pedlen%elts(id*EDGE+RT+1)*dom%len%elts(id*EDGE+RT+1))
    dvelo(id*EDGE+DG+1) = dvelo(id*EDGE+DG+1) + viscosity*( &
         (dom%divu%elts(id+1) - dom%divu%elts(idNE+1)) + &
         (dom%vort%elts(id*TRIAG+LORT+1) - dom%vort%elts(id*TRIAG+UPLT+1)) &
         /dom%pedlen%elts(id*EDGE+DG+1)*dom%len%elts(id*EDGE+DG+1))
    dvelo(id*EDGE+UP+1) = dvelo(id*EDGE+UP+1) + viscosity*( &
         (dom%divu%elts(idN+1) - dom%divu%elts(id+1)) + &
         (dom%vort%elts(idW*TRIAG+LORT+1) - dom%vort%elts(id*TRIAG+UPLT+1)) &
         /dom%pedlen%elts(id*EDGE+UP+1)*dom%len%elts(id*EDGE+UP+1))
  end subroutine diff_mom
end module viscous_mod
