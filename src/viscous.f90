module viscous_mod
  use shared_mod
  use domain_mod
  use arch_mod
  implicit none

contains
  subroutine divu(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      integer idS
      integer idW
      integer idSW
      id = idx(i, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      dom%divu%elts(id+1) = (velo(EDGE*id+UP+1)*dom%pedlen%elts(EDGE*id+UP+1) &
                  - velo(DG+EDGE*id+1)*dom%pedlen%elts(DG+EDGE*id+1) &
                  + velo(EDGE*id+RT+1)*dom%pedlen%elts(EDGE*id+RT+1) &
                  - velo(EDGE*idS+UP+1)*dom%pedlen%elts(EDGE*idS+UP+1) &
                  + velo(DG+EDGE*idSW+1)*dom%pedlen%elts(DG+EDGE*idSW+1) &
                  - velo(EDGE*idW+RT+1)*dom%pedlen%elts(EDGE*idW+RT+1))*dom%areas%elts(id+1)%hex_inv
  end subroutine

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
      tflux(EDGE*id+RT+1) = tflux(EDGE*id+RT+1) - &
              viscosity*dom%pedlen%elts(EDGE*id+RT+1)* &
              (height(idE+1) - height(id+1))/dom%len%elts(EDGE*id+RT+1)
      tflux(DG+EDGE*id+1) = tflux(DG+EDGE*id+1) - &
              viscosity*dom%pedlen%elts(EDGE*id+DG+1)* &
              (height(id+1) - height(idNE+1))/dom%len%elts(DG+EDGE*id+1)
      tflux(EDGE*id+UP+1) = tflux(EDGE*id+UP+1) - &
              viscosity*dom%pedlen%elts(EDGE*id+UP+1)* &
              (height(idN+1) - height(id+1))/dom%len%elts(EDGE*id+UP+1)
  end subroutine

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
  end subroutine
end module
