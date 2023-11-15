module mask_mod
  use domain_mod
  use comm_mpi_mod
contains
  subroutine set_masks (dom, p, i, j, zlev, offs, dims, mask)
    ! Sets all nodes and edges to value mask 
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, mask
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i

    id = idx (i, j, offs, dims)
    id_i = id + 1

    dom%mask_n%elts(id_i)              = mask
    dom%mask_e%elts(EDGE*id:EDGE*id_i) = mask
  end subroutine set_masks
 
  subroutine mask_active
    ! Determine active grid points
    implicit none
    integer :: d, j, k, l

    ! Set active grid at finest scale
    call apply_onescale (mask_tol_vars, level_end, z_null, -1, 2)
    call comm_masks_mpi (level_end)
       
    ! Set active grid at coarser scales
    do l = level_end-1, level_start, -1
       call apply_onescale (mask_tol_vars, l, z_null, -1, 2)

       ! Make parents active
       call apply_interscale (mask_parent_nodes, l, z_null,  0, 1)
       call apply_interscale (mask_parent_edges, l, z_null, -1, 1)
       call comm_masks_mpi (l)
    end do
    call comm_masks_mpi (NONE)
  end subroutine mask_active
  
  subroutine mask_tol_vars (dom, i, j, zlev, offs, dims)
    ! Set active wavelets (determines which grid points are active at adjacent finer scale)
    use utils_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer            :: d, e, id, id_e, id_i, k, l, v
    logical            :: active
    logical, parameter :: penta_adjust = .false.

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1
    l = dom%level%elts(id_i)

    if (dom%mask_n%elts(id_i) == FROZEN) return
   
    ! Scalars
    active = .false.
    do k = zmin, zlevels
       do v = scalars(1), scalars(2)
          if (abs (wav_coeff(v,k)%data(d)%elts(id_i)) >= threshold(v,k) .or. l < level_fill) active = .true.
       end do
    end do

    if (active) then
       dom%mask_n%elts(id_i) = TOLRNZ
    else
       if (dom%mask_n%elts(id_i) > ADJZONE) dom%mask_n%elts(id_i) = ADJZONE
    end if

    ! Vectors
    do e = 1, EDGE
       id_e = EDGE*id + e

       active = .false.
       do k = zmin, zlevels
          if (abs (wav_coeff(S_VELO,k)%data(d)%elts(id_e)) >= threshold(S_VELO,k) .or. l < level_fill) active = .true.
       end do

       if (active) then
          dom%mask_e%elts(id_e) = TOLRNZ
       else
          if (dom%mask_e%elts(id_e) > ADJZONE) dom%mask_e%elts(id_e) = ADJZONE
       end if
    end do
  end subroutine mask_tol_vars

  subroutine mask_parent_nodes (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Make parent node active if any child is active, also make child active if any child neighbours are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idN, idE, idNE, idSW, idS, idW

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idN    = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE    = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE   = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idSW   = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS    = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idW    = idx (i_chd-1, j_chd,   offs_chd, dims_chd)

    if ( dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ) then

       call set_at_least (dom%mask_n%elts(id_par+1), TOLRNZ)
       call set_at_least (dom%mask_n%elts(id_chd+1), TOLRNZ)
    end if
  end subroutine mask_parent_nodes

  subroutine mask_parent_edges (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Nearest neighbours of active edges at coarser scale
    !
    ! Make parent edge active if at least one of the four neighbouring child edges is active
    ! (check two additional child neighbour edges for each parent edge compared to wavetrisk)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: e, id_chd,  id_par, idE, idN, idNE, idNW, idS, idSE, idW

    id_par = idx (i_par, j_par, offs_par, dims_par) ! parent node
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) ! child node

    ! Three neighbours of child node
    idE  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idW  = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idN  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idNE = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idNW = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idS  = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idSE = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    ! Check six child edges neighbouring each parent edge to see if at least one is active
    ! if true make parent edge active
    if ( dom%mask_e%elts(EDGE*id_chd+RT+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1)    == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), TOLRNZ)

    if ( dom%mask_e%elts(EDGE*id_chd+DG+1)   == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)      == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)      == TOLRNZ) &
         call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), TOLRNZ)

    if ( dom%mask_e%elts(EDGE*id_chd+UP+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1)    == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)     == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), TOLRNZ)
  end subroutine mask_parent_edges

  subroutine mask_node_trsk (dom, i, j, zlev, offs, dims)
    ! Add additional TRISK operator stencils needed for nodes
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id  = idx (i, j, offs, dims)
   
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       ! Stencil for divergence of fluxes and Laplacian diffusion
       call div_grad_stencil (dom, i, j, offs, dims)

       ! Hyperdiffusion
       if (Laplace_order == 2) then
          call div_grad_stencil (dom, i+1, j,   offs, dims)
          call div_grad_stencil (dom, i+1, j+1, offs, dims)
          call div_grad_stencil (dom, i,   j+1, offs, dims)
          call div_grad_stencil (dom, i-1, j,   offs, dims)
          call div_grad_stencil (dom, i-1, j-1, offs, dims)
          call div_grad_stencil (dom, i,   j-1, offs, dims)
       end if
    end if
  end subroutine mask_node_trsk

  subroutine mask_edge_trsk (dom, i, j, zlev, offs, dims)
    ! Add additional TRISK operator stencils needed for edges
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_i, idE, idNE, idN

    id  = idx (i, j, offs, dims)
    id_i = id + 1

    if (maxval (dom%mask_e%elts(EDGE*id+1:EDGE*id_i)) >= ADJZONE) then
       ! Mask for remap
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       call set_at_least (dom%mask_n%elts(idE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
       call set_at_least (dom%mask_n%elts(idN+1),  TRSK)
    
       ! Kinetic energy stencil
       call mask_ke_trsk (dom, i,   j,   zlev, offs, dims)
       call mask_ke_trsk (dom, i+1, j,   zlev, offs, dims)
       call mask_ke_trsk (dom, i+1, j+1, zlev, offs, dims)
       call mask_ke_trsk (dom, i,   j+1, zlev, offs, dims)
       
       ! Stencil for gradients of Bernoulli function
       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), TRSK)

       ! Qperp stencil
       call div_grad_stencil (dom, i+1, j,   offs, dims)
       call div_grad_stencil (dom, i+1, j+1, offs, dims)
       call div_grad_stencil (dom, i,   j+1, offs, dims)

       call qe_stencil (dom, i,   j,   offs, dims)
       call qe_stencil (dom, i+1, j,   offs, dims)
       call qe_stencil (dom, i+1, j+1, offs, dims)
       call qe_stencil (dom, i,   j+1, offs, dims)
       call qe_stencil (dom, i-1, j,   offs, dims) 
       call qe_stencil (dom, i-1, j-1, offs, dims) 
       call qe_stencil (dom, i,   j-1, offs, dims) 
       call qe_stencil (dom, i+1, j-1, offs, dims) 
       call qe_stencil (dom, i-1, j+1, offs, dims)

       ! Diffusion
       if (Laplace_order /= 0) then
          call Laplacian_u_stencil (dom, i, j, offs, dims)
          if (Laplace_order == 2) then
             call Laplacian_u_stencil (dom, i+1, j,   offs, dims)
             call Laplacian_u_stencil (dom, i+1, j+1, offs, dims)
             call Laplacian_u_stencil (dom, i+1, j-1, offs, dims)
             call Laplacian_u_stencil (dom, i,   j+1, offs, dims)
             call Laplacian_u_stencil (dom, i-1, j,   offs, dims)
             call Laplacian_u_stencil (dom, i-1, j-1, offs, dims)
             call Laplacian_u_stencil (dom, i,   j-1, offs, dims)
          end if
       end if
    end if
  end subroutine mask_edge_trsk

  subroutine div_grad_stencil (dom, i, j, offs, dims)
    ! Stencil for flux-divergence operator
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idSW, idS

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idSW   = idx (i-1, j-1, offs, dims)
    idS    = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_n%elts(id+1),   TRSK)
    call set_at_least (dom%mask_n%elts(idE+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    call set_at_least (dom%mask_n%elts(idN+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idW+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idSW+1), TRSK)
    call set_at_least (dom%mask_n%elts(idS+1),  TRSK)

    call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK) 
    call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
  end subroutine div_grad_stencil

  subroutine mask_ke_trsk (dom, i, j, zlev, offs, dims)
    ! Kinetic energy stencil
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idSW, idW
    
    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
  end subroutine mask_ke_trsk

  subroutine qe_stencil (dom, i, j, offs, dims)
    ! Stencil for qe 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idS

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idS    = idx (i,   j-1, offs, dims)

    ! Circulation stencil
    call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1), TRSK)

    ! Potential vorticity stencil
    call set_at_least (dom%mask_n%elts(id+1),   TRSK)
    call set_at_least (dom%mask_n%elts(idE+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    call set_at_least (dom%mask_n%elts(idN+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idW+1),  TRSK)
    call set_at_least (dom%mask_n%elts(idS+1),  TRSK)
  end subroutine qe_stencil

  subroutine Laplacian_u_stencil (dom, i, j, offs, dims)
    ! Stencil for Laplacian(u) operators
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call divu_stencil (dom, i,   j,   offs, dims)
    call divu_stencil (dom, i+1, j,   offs, dims)
    call divu_stencil (dom, i+1, j+1, offs, dims)
    call divu_stencil (dom, i,   j+1, offs, dims)
  end subroutine Laplacian_u_stencil

  subroutine divu_stencil (dom, i, j, offs, dims)
    ! Stencil for divu operator
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idW, idS, idSW

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idSW   = idx (i-1, j-1, offs, dims)
    idS    = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK) 
    call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
  end subroutine divu_stencil

  subroutine mask_adj_same_scale_nodes (dom, i, j, zlev, offs, dims)
    ! Nearest neighbours of active nodes at same scale
    ! If at least one of six neighbouring nodes at same scale is active add node to adjacent mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idE, idN, idNE

    ! test node
    id = idx (i, j, offs, dims)

    ! six nearest neighbour nodes at same scale
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    ! Add test node to adjacent mask if at least one of its six neighbours is active
    ! (also needed for div, gradi_e and gradv_e operators)
    if (dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ) &
         call set_at_least(dom%mask_n%elts(id+1), ADJZONE)
  end subroutine mask_adj_same_scale_nodes

  subroutine mask_adj_same_scale (dom, i, j, zlev, offs, dims)
    ! Nearest neighbours at same scale
    ! Add nodes and edges adjacent to adjacent mask if at least one of their nearest neighbours at same scale is active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idE, idN, idNE

    ! node/edge to be checked (may not be active or in adjacent zone)
    id  = idx (i, j, offs, dims)

    ! six neighbour nodes at same scale
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    ! add node to adjacent zone if at least one of six neighbouring nodes is active
    if (dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ) &
         call set_at_least (dom%mask_n%elts(id+1), ADJSPACE)

    ! add edges to adjacent zone if at least one of four neighouring edges is active
    if (dom%mask_e%elts(EDGE*id+DG+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), ADJSPACE)

    if (dom%mask_e%elts(EDGE*id+UP+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id+RT+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), ADJSPACE)

    if (dom%mask_e%elts(EDGE*id+DG+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), ADJSPACE)
  end subroutine mask_adj_same_scale

  subroutine mask_adj_children (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add nearest child neighbours at finer scale to adjacent mask if parent is active and label nodes that can be obtained by restriction
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, idS_par, idW_par, idSW_par, idE_par, idN_par, idNE_par
    integer :: id_chd, idE_chd, idNE_chd, idN2E_chd, id2NE_chd, idN_chd, idW_chd, idNW_chd
    integer :: idS2W_chd, idSW_chd, idS_chd, id2SW_chd, idSE_chd

    ! Includes some edge->node and node->edge cross masking
    id_chd   = idx (i_chd,     j_chd,     offs_chd, dims_chd)

    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW_chd  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSW_chd  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS_chd   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idSE_chd  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)
    idN2E_chd = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    id2SW_chd = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idS2W_chd = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)

    id_par   = idx (i_par,   j_par,   offs_par, dims_par)

    idS_par  = idx (i_par,   j_par-1, offs_par, dims_par)
    idW_par  = idx (i_par-1, j_par,   offs_par, dims_par)
    idSW_par = idx (i_par-1, j_par-1, offs_par, dims_par)
    idE_par  = idx (i_par+1, j_par,   offs_par, dims_par)
    idN_par  = idx (i_par,   j_par+1, offs_par, dims_par)
    idNE_par = idx (i_par+1, j_par+1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= TOLRNZ) then
       ! Nearest neighbour nodes in same cell
       call set_at_least (dom%mask_n%elts(id_chd+1),    ADJZONE)
       call set_at_least (dom%mask_n%elts(idE_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idNE_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idN_chd+1),   ADJZONE)

       ! Needed for prolongation of scalars
       call set_at_least (dom%mask_n%elts(idW_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idNW_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idSW_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idS_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idSE_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idN2E_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id2NE_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id2SW_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idS2W_chd+1), ADJZONE)

       ! Nearest neighbour edges at same scale (also necessary for gradi_e operator)
       ! (at same position as neighbour nodes at finer scale, therefore needed for restriction to coarse node)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    else
       ! If parent node is not active, check three neighbours of parent node
       ! If active, add associated child node and parent edge joining neighbour and node we are checking
       if (dom%mask_n%elts(idN_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idN_chd+1),        ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idNE_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idNE_chd+1),       ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idE_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idE_chd+1),        ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
       end if
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idN_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idS_par+UP+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE) ! Add child node if S neighbour active
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idNE_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idSW_par+DG+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  ! Add child node if SW neighbour active
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idE_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idW_par+RT+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  ! Add child node if W neighbour active
    end if

    ! Add two child edges if at least one neighbouring parent edge is active
    if (dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+UP+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
    end if

    ! Check all five UPLT and LORT triangle edges
    if (dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+DG+1),   ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
    end if

    ! Check all five triangle edges of LORT triangle of current cell and UPLT triangle of southern neighbour
    if (dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_chd+RT+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+RT+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(DG+EDGE*idN_chd+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+UP+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+DG+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+RT+1), ADJZONE)
    end if

    ! Add node to restrict mask at coarse scale if at least one of six neighbouring edges at coarse scale is active
    ! (at positions corresponding to neighbouring nodes at finer scale)
    if (dom%mask_e%elts(EDGE*id_par+RT+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1)  >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idSW_par+DG+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1)  >= ADJSPACE) &
         call set_at_least (dom%mask_n%elts(id_par+1), RESTRCT)
  end subroutine mask_adj_children

  subroutine mask_adj_nodes_edges
    ! Add neighbour edges of active nodes and vice versa (ensure consistency of node and edge adjacent zones)
    implicit none
    integer :: l

    do l = level_end-1, level_start, -1
       call apply_onescale (mask_e_if_both_n, l+1, z_null, 0, 1)
    end do

    call comm_masks_mpi (NONE)

    do l = level_end-1, level_start+1, -1
       call apply_interscale (mask_e_consist, l, z_null, 0, 1)
       call comm_masks_mpi (l+1)
       call apply_interscale (mask_e_consist2, l, z_null, 0, 0)
       call comm_masks_mpi (l)
    end do

    if (level_start < level_end) then
       call apply_interscale (mask_e_consist, level_start, z_null, 0, 1)
       call comm_masks_mpi (level_start+1)
    end if

    do l = level_end-1, level_start+1, -1
       call apply_onescale (mask_n_if_all_e, l+1, z_null, 0, 1)
       call apply_interscale (prolong_node_adjzone, l, z_null, 0, 1)
    end do

    if (level_start+1 <= level_end) call apply_onescale (mask_n_if_all_e, level_start+1, z_null, 0, 1)
  end subroutine mask_adj_nodes_edges

  subroutine prolong_node_adjzone (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent to adjacent mask if child is in adjacent zone for prolongation of nodal quantities in inverse wavelet transform
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par

    id_par = idx (i_par, j_par, offs_par, dims_par) ! id of parent
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) ! id of child

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) call set_at_least(dom%mask_n%elts(id_par+1), ADJZONE)
  end subroutine prolong_node_adjzone

  subroutine mask_restrict_flux (dom, i_par, j_par, zlev, offs_par, dims_par)
    ! Add edges required for flux restriction
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, zlev
    integer, dimension(N_BDRY+1)   :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer                        :: id_par, idE_par, idN_par, idNE_par

    id_par   = idx (i_par,   j_par,   offs_par, dims_par)
    idE_par  = idx (i_par+1, j_par,   offs_par, dims_par)
    idN_par  = idx (i_par,   j_par+1, offs_par, dims_par)
    idNE_par = idx (i_par+1, j_par+1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= ADJSPACE) then
       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    else
       if (dom%mask_n%elts(idN_par+1)  >= ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       if (dom%mask_n%elts(idNE_par+1) >= ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       if (dom%mask_n%elts(idE_par+1)  >= ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    end if
  end subroutine mask_restrict_flux

  

  subroutine complete_masks
    implicit none
    integer :: l

    call apply_onescale (mask_adj_space, level_end, z_null, 0, 1)

    do l = level_end-1, level_start, -1
       call apply_interscale (mask_adj_scale, l, z_null, 0, 0)
       call apply_onescale (mask_e_if_both_n, l+1, z_null, 0, 0)
    end do

    call comm_masks_mpi (NONE)

    do l = level_end-1, level_start+1, -1
       call apply_interscale (mask_e_consist, l, z_null, 0, 1)
       call comm_masks_mpi (l+1)
       call apply_interscale (mask_e_consist2, l, z_null, 0, 0)
       call comm_masks_mpi (l)
    end do

    if (level_start < level_end) then
       call apply_interscale (mask_e_consist, level_start, z_null, 0, 1)
       call comm_masks_mpi (level_start+1)
    end if

    do l = level_end-1, level_start+1, -1
       call apply_onescale (mask_n_if_all_e, l+1, z_null, 0, 1)
       call apply_interscale (prolong_node_adjzone, l,   z_null, 0, 1)
    end do

    if (level_start+1 <= level_end) call apply_onescale (mask_n_if_all_e, level_start+1, z_null, 0, 1)
  end subroutine complete_masks

  subroutine mask_adj_space (dom, i, j, zlev, offs, dims)
    ! Nearest neighbours of active nodes at same scale
    ! If at least one of six neighbouring nodes at same scale is active add node to adjacent mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idE, idN, idNE

    ! test node
    id = idx (i, j, offs, dims)

    ! six nearest neighbour nodes at same scale
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    ! Add test node to adjacent mask if at least one of its six neighbours is active
    ! (also needed for div, gradi_e and gradv_e operators)
    if (dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ) &
         call set_at_least(dom%mask_n%elts(id+1), ADJZONE)
  end subroutine mask_adj_space

  subroutine mask_adj_scale (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add nearest child neighbours at finer scale to adjacent mask if parent is active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, idS_par, idW_par, idSW_par, idE_par, idN_par, idNE_par
    integer :: id_chd, idE_chd, idNE_chd, idN2E_chd, id2NE_chd, idN_chd, idW_chd, idNW_chd
    integer :: idS2W_chd, idSW_chd, idS_chd, id2SW_chd, idSE_chd

    ! Includes some edge->node and node->edge cross masking
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW_chd  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSW_chd  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS_chd   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idSE_chd  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)
    idN2E_chd = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    id2SW_chd = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idS2W_chd = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)

    id_par   = idx (i_par,   j_par,   offs_par, dims_par)

    idS_par  = idx (i_par,   j_par-1, offs_par, dims_par)
    idW_par  = idx (i_par-1, j_par,   offs_par, dims_par)
    idSW_par = idx (i_par-1, j_par-1, offs_par, dims_par)
    idE_par  = idx (i_par+1, j_par,   offs_par, dims_par)
    idN_par  = idx (i_par,   j_par+1, offs_par, dims_par)
    idNE_par = idx (i_par+1, j_par+1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= TOLRNZ) then
       ! Nearest neighbour nodes in same cell at finer scale
       call set_at_least (dom%mask_n%elts(id_chd+1),    ADJZONE)
       call set_at_least (dom%mask_n%elts(idE_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idNE_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idN_chd+1),   ADJZONE)

       ! Needed for prolongation of scalars
       call set_at_least (dom%mask_n%elts(idW_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idNW_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idSW_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idS_chd+1),   ADJZONE)
       call set_at_least (dom%mask_n%elts(idSE_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idN2E_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id2NE_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id2SW_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idS2W_chd+1), ADJZONE)

       ! Nearest neighbour edges at same scale (also necessary for gradi_e operator)
       ! (at same position as neighbour nodes at finer scale, therefore needed for restriction to coarse node)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    else
       ! If parent node is not active, check three neighbours of parent node
       ! If active, add associated child node and parent edge joining neighbour and node we are checking
       if (dom%mask_n%elts(idN_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idN_chd+1),        ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idNE_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idNE_chd+1),       ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idE_par+1) >= TOLRNZ) then
          call set_at_least (dom%mask_n%elts(idE_chd+1),        ADJZONE)
          call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
       end if
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idN_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idS_par+UP+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE) ! Add child node if S neighbour active
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idNE_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idSW_par+DG+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  ! Add child node if SW neighbour active
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idE_chd+1), ADJZONE) ! Add associated chid node
       if (dom%mask_e%elts(EDGE*idW_par+RT+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  ! Add child node if W neighbour active
    end if

    ! Add two child edges if at least one neighbouring parent edge is active
    if (dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+UP+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
    end if

    ! Check all five UPLT and LORT triangle edges
    if (dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+DG+1),   ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
    end if

    ! Check all five triangle edges of LORT triangle of current cell and UPLT triangle of southern neighbour
    if (dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_chd+RT+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+RT+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(DG+EDGE*idN_chd+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ) then

       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+DG+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+RT+1), ADJZONE)
    end if

    ! Add node to restrict mask at coarse scale if at least one of six neighbouring edges at coarse scale is active
    ! (at positions corresponding to neighbouring nodes at finer scale)
    if (dom%mask_e%elts(EDGE*id_par+RT+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)   >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1)  >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idSW_par+DG+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1)  >= ADJSPACE) &
         call set_at_least (dom%mask_n%elts(id_par+1), RESTRCT)
  end subroutine mask_adj_scale

  subroutine mask_e_if_both_n (dom, i, j, zlev, offs, dims)
    ! Add edge to adjacent zone if both neighbouring nodes are active
    ! (also needed for div and gradv_e operator)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN, idNE

    id   = idx (i,   j,     offs, dims)
    idE  = idx (i+1, j,     offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    if ( dom%mask_n%elts(id+1)  >= ADJZONE .and. &
         dom%mask_n%elts(idE+1) >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), ADJZONE)
    end if

    if ( dom%mask_n%elts(idNE+1) >= ADJZONE .and. &
         dom%mask_n%elts(id+1)   >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), ADJZONE)
    end if

    if ( dom%mask_n%elts(id+1)  >= ADJZONE .and. &
         dom%mask_n%elts(idN+1) >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), ADJZONE)
    end if
  end subroutine mask_e_if_both_n

  subroutine mask_e_consist (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id, id_par, idN, idS, idE, idW, idNE, idNW, idSE

    id     = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    idN  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idS  = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idW  = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNE = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idNW = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSE = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id+UP+1)   >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+UP+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1) >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  >= CONSIST) then

       call set_at_least (dom%mask_e%elts(EDGE*idN+UP+1),    ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),     ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id+DG+1)   >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+DG+1) >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1) >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1) >= CONSIST) then

       call set_at_least(dom%mask_e%elts(EDGE*idNE+DG+1),   ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id+DG+1),     ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id+RT+1)   >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+RT+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1) >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)  >= CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)  >= CONSIST) then

       call set_at_least(dom%mask_e%elts(EDGE*idE+RT+1),    ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id+RT+1),     ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), ADJZONE)
    end if
  end subroutine mask_e_consist

  subroutine mask_e_consist2 (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Make parent edges active if neighouring child edges are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par, idN, idE, idS, idW, idNE, idSW, idNW, idSE

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    ! Eight nearest neighbours of child node
    idN  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    idE  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    idNE = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)
    idS  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    idW  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    idSW = idx (i_chd-2, j_chd-2, offs_chd, dims_chd)
    idNW = idx (i_chd-2, j_chd+2, offs_chd, dims_chd)
    idSE = idx (i_chd+2, j_chd-2, offs_chd, dims_chd)

    ! Check if edges of child nodes are in adjacent zone
    if ( dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd+DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+DG+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+UP+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+UP+1)    >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*idW+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+UP+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1)   >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*idSW+RT+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+DG+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd+DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+UP+1)    >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1)   >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)    >= ADJZONE) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), ADJZONE)
    end if
  end subroutine mask_e_consist2

  subroutine mask_adj_space2 (dom, i, j, zlev, offs, dims)
    ! Nearest neighbours at same scale
    ! Add nodes and edges adjacent to adjacent mask if at least one of their nearest neighbours at same scale is active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idS, idW, idSW, idE, idN, idNE

    ! node/edge to be checked (may not be active or in adjacent zone)
    id  = idx (i, j, offs, dims)

    ! six neighbour nodes at same scale
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    ! add node to adjacent zone if at least one of six neighbouring nodes is active
    if (dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ) &
         call set_at_least (dom%mask_n%elts(id+1), ADJSPACE)

    ! add edges to adjacent zone if at least one of four neighouring edges is active
    if (dom%mask_e%elts(EDGE*id+DG+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), ADJSPACE)

    if (dom%mask_e%elts(EDGE*id+UP+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id+RT+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), ADJSPACE)

    if (dom%mask_e%elts(EDGE*id+DG+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), ADJSPACE)
  end subroutine mask_adj_space2

  subroutine mask_n_if_all_e (dom, i, j, zlev, offs, dims)
    ! Add node to adjacent zone if all neighbouring edges are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idW, idS, idSW

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idSW = idx (i-1, j-1, offs, dims)

    if ( dom%mask_e%elts(EDGE*id+UP+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+DG+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+RT+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idS+UP+1)  >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idSW+DG+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idW+RT+1)  >= ADJZONE) then

       call set_at_least (dom%mask_n%elts(id+1), ADJZONE)
    end if
  end subroutine mask_n_if_all_e

  subroutine init_mask_mod
    implicit none
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_domain_mod
    call init_comm_mod
    initialized = .True.
  end subroutine init_mask_mod

  subroutine init_masks
    implicit none
    integer :: d, i, num

    do d = 1, size(grid)
       call init (grid(d)%mask_n, 1)
       call init (grid(d)%mask_e, EDGE)
       call init (grid(d)%level, 1)
    end do

    do d = 1, size(grid)
       num = grid(d)%node%length-1
       call extend (grid(d)%mask_n, num,      TOLRNZ)
       call extend (grid(d)%mask_e, EDGE*num, TOLRNZ)
       call extend (grid(d)%level,  num,      min_level-1)
    end do
  end subroutine init_masks

  subroutine set_at_least (mask, typ)
    implicit none
    integer :: mask, typ

    if (mask < typ) mask = typ
  end subroutine set_at_least
end module mask_mod
