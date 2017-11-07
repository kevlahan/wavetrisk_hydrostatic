module mask_mod
  use param_mod
  use shared_mod
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use wavelet_mod
  implicit none
  real(8) tol_velo
  real(8) tol_mass
  real(8) tol_temp

contains
  subroutine init_mask_mod()
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_domain_mod()
    call init_comm_mod()
    tol_velo = 0.0_8
    tol_mass = 0.0_8
    tol_temp = 0.0_8
    initialized = .True.
  end subroutine init_mask_mod

  subroutine mask_e_consist(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    integer id
    integer idN, idS
    integer idE, idW
    integer idNE, idNW, idSE
    integer id_par

    id     = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN    = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE    = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idS    = idx(i_chd,     j_chd - 1, offs_chd, dims_chd)
    idW    = idx(i_chd - 1, j_chd,     offs_chd, dims_chd)
    idNE   = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idNW   = idx(i_chd - 1, j_chd + 1, offs_chd, dims_chd)
    idSE   = idx(i_chd + 1, j_chd - 1, offs_chd, dims_chd)
    id_par = idx(i_par,     j_par,     offs_par, dims_par)

    if ( dom%mask_e%elts(EDGE*id+UP+1)   .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+UP+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1) .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  .ge. CONSIST) then

       call set_at_least(dom%mask_e%elts(EDGE*idN+UP+1),    ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id+UP+1),     ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id+DG+1)   .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+DG+1) .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1) .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1) .ge. CONSIST) then

       call set_at_least(dom%mask_e%elts(DG+EDGE*idNE+1),   ADJZONE)
       call set_at_least(dom%mask_e%elts(DG+EDGE*id+1),     ADJZONE)
       call set_at_least(dom%mask_e%elts(DG+EDGE*id_par+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id+RT+1)   .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+RT+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1) .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)  .ge. CONSIST .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)  .ge. CONSIST) then

       call set_at_least(dom%mask_e%elts(EDGE*idE+RT+1),    ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id+RT+1),     ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), ADJZONE)
    end if

  end subroutine mask_e_consist

  subroutine mask_adj_space(dom, i, j, zlev, offs, dims)
    ! Label nodes adjacent to active nodes
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer idS
    integer idW
    integer idSW
    integer idE
    integer idN
    integer idNE

    ! id of active node
    id   = idx(i,     j,     offs, dims)

    ! id of six nearest neighbours at same scale
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if ( dom%mask_n%elts(idN+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  .eq. TOLRNZ) then

       call set_at_least(dom%mask_n%elts(id+1), ADJZONE)
    end if
  end subroutine mask_adj_space

  subroutine mask_adj_space2(dom, i, j, zlev, offs, dims)
    ! Label nodes and edges adjacent to active nodes and edges
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer idS
    integer idW
    integer idSW
    integer idE
    integer idN
    integer idNE

    ! id of active node
    id  = idx(i, j, offs, dims)

    ! id of six neighbours at same scale
    idS  = idx(i,     j - 1, offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if ( dom%mask_n%elts(idN+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idE+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  .eq. TOLRNZ) then

       call set_at_least(dom%mask_n%elts(id+1), ADJSPACE)
    end if

    if ( dom%mask_e%elts(EDGE*id+DG+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+RT+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id+UP+1), ADJSPACE)
    end if

    if ( dom%mask_e%elts(EDGE*id+UP+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id+RT+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id+DG+1), ADJSPACE)
    end if

    if ( dom%mask_e%elts(EDGE*id+DG+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+UP+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id+RT+1), ADJSPACE)
    end if

  end subroutine mask_adj_space2

  subroutine mask_active_edges(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    !label active edge data (e.g. velocity) of active children
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_chd
    integer id_par
    integer id, idE, idN, idNE
    integer e

    id_par = idx(i_par, j_par, offs_par, dims_par) ! parent node

    id     = idx(i_chd, j_chd, offs_chd, dims_chd) ! child node

    ! three neighbours of child node
    idE  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idN  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idNE = idx(i_chd + 1, j_chd+ 1,  offs_chd, dims_chd)

    ! check four child nodes of each type of parent velocity node to see of any are active
    if ( dom%mask_e%elts(EDGE*id+RT+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+RT+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+DG+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), TOLRNZ)
    end if

    if ( dom%mask_e%elts(EDGE*id+DG+1)   .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+DG+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), TOLRNZ)
    end if

    if ( dom%mask_e%elts(EDGE*id+UP+1)  .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+UP+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) .eq. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+DG+1) .eq. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), TOLRNZ)
    end if

  end subroutine mask_active_edges

  subroutine init_masks()
    integer d
    integer num
    integer i

    do d = 1, size(grid)
       call init(grid(d)%mask_n, 1)
       call init(grid(d)%mask_e, EDGE)
       call init(grid(d)%level, 1)
    end do

    do d = 1, size(grid)
       num = grid(d)%node%length-1
       call extend(grid(d)%mask_n, num, TOLRNZ)
       call extend(grid(d)%mask_e, EDGE*num, TOLRNZ)
       call extend(grid(d)%level, num, min_level-1)
    end do
  end subroutine init_masks

  subroutine inj_p_adjzone(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! add parent to adjacent zone if child is in adjacent zone
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    integer id_par
    integer id

    id_par = idx(i_par, j_par, offs_par, dims_par) ! id of parent
    id     = idx(i_chd, j_chd, offs_chd, dims_chd) ! id of child

    if (dom%mask_n%elts(id+1) .ge. ADJZONE) call set_at_least(dom%mask_n%elts(id_par+1), ADJZONE)

  end subroutine inj_p_adjzone

  subroutine set_masks (dom, p, i, j, zlev, offs, dims, mask)
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, mask

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .eq. FROZEN) return

    dom%mask_n%elts(id+1) = mask
    do e = 1, EDGE
       dom%mask_e%elts(EDGE*id+e) = mask
    end do
  end subroutine set_masks

  subroutine mask_adjacent (wav)
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: wav
    
    integer :: d, j, k, l

    call update_array_bdry1 (wav, level_start, level_end)
    
    do k = 1, zlevels
       do d = 1, size(grid)
          wc_m => wav(S_MASS,k)%data(d)%elts
          wc_t => wav(S_TEMP,k)%data(d)%elts
          wc_u => wav(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (mask_adj, grid(d), grid(d)%lev(level_end)%elts(j), k, -1, 2)
          end do
          nullify(wc_m, wc_t, wc_u)
       end do
    
       do l = level_end-1, level_start, -1
          do d = 1, size(grid)
             wc_m => wav(S_MASS,k)%data(d)%elts
             wc_t => wav(S_TEMP,k)%data(d)%elts
             wc_u => wav(S_VELO,k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (mask_adj, grid(d), grid(d)%lev(l)%elts(j), k, -1, 2)
             end do
             nullify(wc_m, wc_t, wc_u)
          end do
         
       end do
    end do
  end subroutine mask_adjacent

  subroutine mask_adj (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer e
    integer mask_mass, mask_temp, mask_velo

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .eq. FROZEN) return

    call set_adj_mask(dom%mask_n%elts(id+1), wc_m(id+1), tol_mass)
    call set_adj_mask(dom%mask_n%elts(id+1), wc_t(id+1), tol_temp)
    do e = 1, EDGE
       call set_adj_mask(dom%mask_e%elts(EDGE*id+e), wc_u(EDGE*id+e), tol_velo)
    end do
  end subroutine mask_adj

  subroutine set_adj_mask (mask, wc, tol)
    ! add adjacent zone points to mask
    integer, intent(inout) :: mask
    real(8), intent(in)    :: wc, tol

    if (mask .gt. ADJZONE) mask = ADJZONE
  end subroutine set_adj_mask
  
  subroutine mask_active (wav)
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: wav
    
    integer :: d, j, k, l

    call update_array_bdry1 (wav, level_start, level_end)
    
    do k = 1, zlevels
       do d = 1, size(grid)
          wc_m => wav(S_MASS,k)%data(d)%elts
          wc_t => wav(S_TEMP,k)%data(d)%elts
          wc_u => wav(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(mask_tol, grid(d), grid(d)%lev(level_end)%elts(j), k, -1, 2)
          end do
          nullify(wc_m, wc_t, wc_u)
       end do
    
       do l = level_end-1, level_start, -1
          do d = 1, size(grid)
             wc_m => wav(S_MASS,k)%data(d)%elts
             wc_t => wav(S_TEMP,k)%data(d)%elts
             wc_u => wav(S_VELO,k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch(mask_tol, grid(d), grid(d)%lev(l)%elts(j), k, -1, 2)
             end do
             nullify(wc_m, wc_t, wc_u)
          end do
         
       end do
    end do

    do l = level_end-1, level_start, -1
       call apply_interscale(mask_active_nodes, l, z_null,  0, 1)
       call apply_interscale(mask_active_edges, l, z_null, -1, 1)
       call comm_masks_mpi(l)
    end do
  end subroutine mask_active
  
  subroutine mask_tol (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer e
    integer mask_mass, mask_temp, mask_velo

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .eq. FROZEN) return

    call set_active_mask(dom%mask_n%elts(id+1), wc_m(id+1), tol_mass)
    call set_active_mask(dom%mask_n%elts(id+1), wc_t(id+1), tol_temp)
    do e = 1, EDGE
       call set_active_mask(dom%mask_e%elts(EDGE*id+e), wc_u(EDGE*id+e), tol_velo)
    end do
  end subroutine mask_tol

  subroutine set_active_mask (mask, wc, tol)
    ! add active points to mask
    integer, intent(inout) :: mask
    real(8), intent(in)    :: wc, tol

    if (abs(wc) .ge. tol) mask = TOLRNZ
  end subroutine set_active_mask

  subroutine mask_restrict_flux(dom, i_par, j_par, zlev, offs_par, dims_par)
    type(Domain) dom
    integer i_par, j_par, zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer id_par, idE_par, idN_par, idNE_par

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)
    idE_par  = idx(i_par + 1, j_par,     offs_par, dims_par)
    idN_par  = idx(i_par,     j_par + 1, offs_par, dims_par)
    idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) .ge. ADJSPACE) then
       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    else
       if (dom%mask_n%elts(idN_par+1)  .ge. ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       if (dom%mask_n%elts(idNE_par+1) .ge. ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       if (dom%mask_n%elts(idE_par+1)  .ge. ADJSPACE) call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    end if

  end subroutine mask_restrict_flux

  subroutine mask_adj_scale(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_chd
    integer id_chd
    integer idN_chd
    integer idE_chd
    integer idNE_chd
    integer id_par
    integer idS_par
    integer idW_par
    integer idSW_par
    integer idE_par
    integer idN_par
    integer idNE_par

    ! Includes some u->p and p->u cross masking
    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)
    idS_par  = idx(i_par,     j_par - 1, offs_par, dims_par)
    idW_par  = idx(i_par - 1, j_par,     offs_par, dims_par)
    idSW_par = idx(i_par - 1, j_par - 1, offs_par, dims_par)
    idE_par  = idx(i_par + 1, j_par,     offs_par, dims_par)
    idN_par  = idx(i_par,     j_par + 1, offs_par, dims_par)
    idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) .ge. TOLRNZ) then
       call set_at_least(dom%mask_n%elts(id_chd+1),   ADJZONE)
       call set_at_least(dom%mask_n%elts(idN_chd+1),  ADJZONE)
       call set_at_least(dom%mask_n%elts(idNE_chd+1), ADJZONE)
       call set_at_least(dom%mask_n%elts(idE_chd+1),  ADJZONE)

       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
    else
       if (dom%mask_n%elts(idN_par+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(idN_chd+1),        ADJZONE)
          call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idNE_par+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(idNE_chd+1),       ADJZONE)
          call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       end if
       if (dom%mask_n%elts(idE_par+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(idE_chd+1),        ADJZONE)
          call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
       end if
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1) .ge. TOLRNZ) then
       call set_at_least(dom%mask_n%elts(idN_chd+1), ADJZONE)
       if (dom%mask_e%elts(EDGE*idS_par+UP+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(id_chd+1), ADJZONE)
       end if
    end if

    if ( dom%mask_e%elts(EDGE*id_par+UP+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1) .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+DG+1) .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) .ge. TOLRNZ) then
       call set_at_least(dom%mask_e%elts(EDGE*id_chd+UP+1),  ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) .ge. TOLRNZ) then
       call set_at_least(dom%mask_n%elts(idNE_chd+1), ADJZONE)
       if (dom%mask_e%elts(EDGE*idSW_par+DG+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(id_chd+1), ADJZONE)
       end if
    end if

    if ( dom%mask_e%elts(EDGE*id_par+DG+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+RT+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) .ge. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(DG+EDGE*id_chd+1),   ADJZONE)
       call set_at_least(dom%mask_e%elts(DG+EDGE*idNE_chd+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) .ge. TOLRNZ) then
       call set_at_least(dom%mask_n%elts(idE_chd+1), ADJZONE)
       if (dom%mask_e%elts(EDGE*idW_par+RT+1) .ge. TOLRNZ) then
          call set_at_least(dom%mask_n%elts(id_chd+1), ADJZONE)
       end if
    end if

    if ( dom%mask_e%elts(EDGE*id_par+RT+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1) .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+DG+1) .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) .ge. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_chd+RT+1),  ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id_par+UP+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) .ge. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*idN_chd+RT+1),  ADJZONE)
       call set_at_least(dom%mask_e%elts(DG+EDGE*idN_chd+1),  ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idNE_chd+UP+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id_par+RT+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)  .ge. TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) .ge. TOLRNZ) then

       call set_at_least(dom%mask_e%elts(EDGE*idE_chd+UP+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idE_chd+DG+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idNE_chd+RT+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id_par+RT+1)   .ge. ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+DG+1)   .ge. ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par+UP+1)   .ge. ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1)  .ge. ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idSW_par+DG+1) .ge. ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1)  .ge. ADJSPACE) then

       call set_at_least(dom%mask_n%elts(id_par+1), RESTRCT)
    end if

  end subroutine mask_adj_scale

  subroutine set_at_least(mask, typ)
    integer mask
    integer typ

    if (mask .lt. typ) mask = typ
  end subroutine set_at_least

  subroutine mask_e_consist2(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    integer id
    integer idN, idE
    integer idS, idW
    integer idNE, idSW
    integer idNW, idSE
    integer id_par

    id   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN  = idx(i_chd,     j_chd + 2, offs_chd, dims_chd)
    idE  = idx(i_chd + 2, j_chd,     offs_chd, dims_chd)
    idNE = idx(i_chd + 2, j_chd + 2, offs_chd, dims_chd)
    idS  = idx(i_chd,     j_chd - 2, offs_chd, dims_chd)
    idW  = idx(i_chd - 2, j_chd,     offs_chd, dims_chd)
    idSW = idx(i_chd - 2, j_chd - 2, offs_chd, dims_chd)
    idNW = idx(i_chd - 2, j_chd + 2, offs_chd, dims_chd)
    idSE = idx(i_chd + 2, j_chd - 2, offs_chd, dims_chd)

    id_par = idx(i_par, j_par, offs_par, dims_par)

    if ( dom%mask_e%elts(EDGE*id+RT+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*id+DG+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+DG+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+UP+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+UP+1)  .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+UP+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*idW+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*id+RT+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+UP+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*id+UP+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1) .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*idSW+RT+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW+DG+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*id+DG+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS+UP+1)  .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1) .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*id+UP+1)   .ge. ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)  .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id_par+RT+1), ADJZONE)
    end if

  end subroutine mask_e_consist2

  subroutine mask_active_nodes(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    !mask active nodal data (e.g. mass and potential temperature)
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_chd
    integer id_par, id_chd
    integer idN
    integer idE
    integer idNE
    integer idSW
    integer idS
    integer idW

    id_par = idx(i_par, j_par, offs_par, dims_par)

    id_chd = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN    = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE    = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE   = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idSW   = idx(i_chd - 1, j_chd - 1, offs_chd, dims_chd)
    idS    = idx(i_chd,     j_chd - 1, offs_chd, dims_chd)
    idW    = idx(i_chd - 1, j_chd,     offs_chd, dims_chd)

    if ( dom%mask_n%elts(idE+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idN+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  .eq. TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) .eq. TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  .eq. TOLRNZ) then

       call set_at_least(dom%mask_n%elts(id_par+1), TOLRNZ)
       call set_at_least(dom%mask_n%elts(id_chd+1), TOLRNZ)
    end if
  end subroutine mask_active_nodes

  subroutine mask_n_if_all_e(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer idW
    integer idS
    integer idSW

    id   = idx(i,     j,     offs, dims)
    idW  = idx(i - 1, j,     offs, dims)
    idS  = idx(i,     j - 1, offs, dims)
    idSW = idx(i - 1, j - 1, offs, dims)

    if ( dom%mask_e%elts(EDGE*id+UP+1)   .ge. ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+DG+1)   .ge. ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+RT+1)   .ge. ADJZONE .and. &
         dom%mask_e%elts(EDGE*idS+UP+1)  .ge. ADJZONE .and. &
         dom%mask_e%elts(EDGE*idSW+DG+1) .ge. ADJZONE .and. &
         dom%mask_e%elts(EDGE*idW+RT+1)  .ge. ADJZONE) then

       call set_at_least(dom%mask_n%elts(id+1), ADJZONE)
    end if
  end subroutine mask_n_if_all_e

  subroutine mask_e_if_both_n(dom, i, j, zlev, offs, dims)
    type(Domain) dom
    integer i
    integer j
    integer zlev
    integer, dimension(N_BDRY + 1) :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer id
    integer idE
    integer idN
    integer idNE

    id   = idx(i,     j,     offs, dims)
    idE  = idx(i + 1, j,     offs, dims)
    idN  = idx(i,     j + 1, offs, dims)
    idNE = idx(i + 1, j + 1, offs, dims)

    if ( dom%mask_n%elts(id+1)  .ge. ADJZONE .and. &
         dom%mask_n%elts(idE+1) .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id+RT+1), ADJZONE)
    end if

    if ( dom%mask_n%elts(idNE+1) .ge. ADJZONE .and. &
         dom%mask_n%elts(id+1)   .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id+DG+1), ADJZONE)
    end if

    if ( dom%mask_n%elts(id+1)  .ge. ADJZONE .and. &
         dom%mask_n%elts(idN+1) .ge. ADJZONE) then

       call set_at_least(dom%mask_e%elts(EDGE*id+UP+1), ADJZONE)
    end if

  end subroutine mask_e_if_both_n

  subroutine complete_masks()
    integer l

    call apply_onescale(mask_adj_space, level_end, z_null, 0, 1)

    do l = level_end-1, level_start, -1
       call apply_interscale(mask_adj_scale, l, z_null, 0, 1)
       call apply_onescale(mask_e_if_both_n, l+1, z_null, 0, 0)
    end do

    call comm_masks_mpi(NONE)

    do l = level_end-1, level_start+1, -1
       call apply_interscale(mask_e_consist, l, z_null, 0, 1)
       call comm_masks_mpi(l+1)
       call apply_interscale(mask_e_consist2, l, z_null, 0, 0)
       call comm_masks_mpi(l)
    end do

    if (level_start .lt. level_end) then
       call apply_interscale(mask_e_consist, level_start, z_null, 0, 1)
       call comm_masks_mpi(level_start+1)
    end if

    do l = level_end-1, level_start+1, -1
       call apply_onescale(mask_n_if_all_e, l+1, z_null, 0, 1)
       call apply_interscale(inj_p_adjzone, l,   z_null, 0, 1)
    end do

    if (level_start+1 .le. level_end) call apply_onescale(mask_n_if_all_e, level_start+1, z_null, 0, 1)
  end subroutine complete_masks

end module mask_mod
