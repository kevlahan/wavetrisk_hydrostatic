module mask_mod
  use param_mod
  use shared_mod
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use wavelet_mod

  implicit none
  real(8) :: tol_velo, tol_mass, tol_temp
contains
  subroutine set_masks (dom, p, i, j, zlev, offs, dims, mask)
    ! Sets all nodes and edges to value mask unless node is FROZEN
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, mask
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) .eq. FROZEN) return

    dom%mask_n%elts(id+1) = mask
    do e = 1, EDGE
       dom%mask_e%elts(EDGE*id+e) = mask
    end do
  end subroutine set_masks

  subroutine check_masks (dom, p, i, j, zlev, offs, dims, mask)
    ! Checks if some nodes or edges have value mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev, mask
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1)==mask) write(6,*) 'node ', id, ' has value ', mask
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e)==mask) write(6,*) 'edge ', id, e, ' has value ', mask
    end do
  end subroutine check_masks

  subroutine mask_adjacent_initial
    implicit none
    integer :: k, l

    do k = 1, zlevels
       call apply_onescale (mask_adj, level_end, k, -1, 2)
       do l = level_end-1, level_start, -1
          call apply_onescale (mask_adj, l, k, -1, 2)
       end do
    end do
  end subroutine mask_adjacent_initial

  subroutine mask_adj (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) == FROZEN) return

    if(dom%mask_n%elts(id+1) > ADJZONE) dom%mask_n%elts(id+1) = ADJZONE
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) > ADJZONE) dom%mask_e%elts(EDGE*id+e) = ADJZONE
    end do
  end subroutine mask_adj

  subroutine mask_active (wavelet)
    implicit none
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: wavelet

    integer :: d, j, k, l

    call update_array_bdry1 (wavelet, level_start, level_end)

    do k = 1, zlevels
       do d = 1, size(grid)
          wc_m => wavelet(S_MASS,k)%data(d)%elts
          wc_t => wavelet(S_TEMP,k)%data(d)%elts
          wc_u => wavelet(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch (mask_tol, grid(d), grid(d)%lev(level_end)%elts(j), k, -1, 2)
          end do
          nullify (wc_m, wc_t, wc_u)
       end do
       call comm_masks_mpi (level_end)

       do l = level_end-1, level_start, -1
          do d = 1, size(grid)
             wc_m => wavelet(S_MASS,k)%data(d)%elts
             wc_t => wavelet(S_TEMP,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (mask_tol, grid(d), grid(d)%lev(l)%elts(j), k, -1, 2)
             end do
             nullify (wc_m, wc_t, wc_u)
          end do
          call comm_masks_mpi (l)
       end do
    end do
  end subroutine mask_active

  subroutine mask_tol (dom, i, j, zlev, offs, dims)
    ! Set active wavelets (determines which grid points are active at adjacent finer scale)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id

    id = idx(i, j, offs, dims)

    if (dom%mask_n%elts(id+1) == FROZEN) return

    call set_active_mask (dom%mask_n%elts(id+1), wc_m(id+1), tol_mass)
    call set_active_mask (dom%mask_n%elts(id+1), wc_t(id+1), tol_temp)
    do e = 1, EDGE
       call set_active_mask (dom%mask_e%elts(EDGE*id+e), wc_u(EDGE*id+e), tol_velo)
    end do
  end subroutine mask_tol

  subroutine set_active_mask (mask, wc, tol)
    ! add active points to mask
    implicit none
    integer, intent(inout) :: mask
    real(8), intent(in)    :: wc, tol

    if (abs(wc) >= tol) mask = TOLRNZ
  end subroutine set_active_mask

  subroutine mask_perfect_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add nodes required for scalar wavelet transform based on which scalar wavelets are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd    = idx(i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd   = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx(i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx(i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx(i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx(i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx(i_chd+2, j_chd+2, offs_chd, dims_chd)

    if (dom%mask_n%elts(idE_chd+1) >= ADJZONE) then
       call set_at_least (dom%mask_n%elts(id_chd+1),    ADJZONE)
       call set_at_least (dom%mask_n%elts(id2E_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(id2NE_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id2S_chd+1),  ADJZONE)
    end if

    if (dom%mask_n%elts(idNE_chd+1) >= ADJZONE) then
       call set_at_least (dom%mask_n%elts(id2NE_chd+1), ADJZONE)
       call set_at_least (dom%mask_n%elts(id_chd+1),    ADJZONE)
       call set_at_least (dom%mask_n%elts(id2E_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(id2N_chd+1),  ADJZONE)
    end if

    if (dom%mask_n%elts(idN_chd+1) >= ADJZONE) then
       call set_at_least (dom%mask_n%elts(id_chd+1),    ADJZONE)
       call set_at_least (dom%mask_n%elts(id2N_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(id2W_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(id2NE_chd+1), ADJZONE)
    end if
  end subroutine mask_perfect_scalar

  subroutine mask_perfect_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add nodes at finer scale required for velocity wavelet transform based on which velocity wavelets are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: e, id1, id2, id_chd, idN_chd, idE_chd, idNE_chd, idN2E_chd, id2NE_chd, idN2_chd, idE2_chd

    id_chd    = idx(i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd   = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN2E_chd = idx(i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE_chd = idx(i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN2_chd  = idx(i_chd,   j_chd+2, offs_chd, dims_chd)
    idE2_chd  = idx(i_chd+2, j_chd,   offs_chd, dims_chd)

    ! Outer
    do e = 1, EDGE
       id1 = idx(i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2 = idx(i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       if (dom%mask_e%elts(EDGE*id2+e) >= ADJZONE .or. dom%mask_e%elts(EDGE*id1+e) >= ADJZONE) then
          call set_at_least (dom%mask_e%elts(EDGE*id2+e), ADJZONE)
          call mask_outer_velo (dom, i_par, j_par, e - 1, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)
       end if
    end do

    ! Inner
    if (dom%mask_e%elts(EDGE*idE_chd+UP+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+RT+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+DG+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*idE_chd+DG+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE2_chd+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN2E_chd+UP+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*idN_chd+RT+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id_chd+DG+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*idN_chd+DG+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN2_chd+RT+1), ADJZONE)
    end if

    if (dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE) then
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id2NE_chd+RT+1), ADJZONE)
    end if
  end subroutine mask_perfect_velo

  subroutine mask_outer_velo (dom, i, j, e, offs, dims, i_chd, j_chd, offs_chd, dims_chd)
    ! Interpolate outer velocities to fine edges
    implicit none
    real(8)                        :: interp_outer_velo
    type(Domain)                   :: dom
    integer                        :: i, j, e, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims, dims_chd

    integer :: ide

    ide = idx(i_chd+end_pt(1,2,e+1),j_chd+end_pt(2,2,e+1),offs_chd,dims_chd)

    call set_at_least (dom%mask_e%elts(idx(i,j,offs,dims)*EDGE+e+1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + opp_no(1,1,e+1), j + opp_no(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + opp_no(1,1,e+1), j + opp_no(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + opp_no(1,2,e+1), j + opp_no(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims) + 1), ADJZONE)

    call set_at_least (dom%mask_e%elts(ed_idx(i + opp_no(1,2,e+1), j + opp_no(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims) + 1), ADJZONE)
  end subroutine mask_outer_velo

  subroutine mask_perfect_velo_penta (dom, p, c, offs, dims, zlev)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                        :: id_chd, idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then
       id_chd  = idx(0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx(0, LAST,   offs_chd, dims_chd)

       if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
          call set_at_least (dom%mask_e%elts(idx( 0, PATCH_SIZE, offs, dims)*EDGE + UP + 1), ADJZONE)
          call set_at_least (dom%mask_e%elts(idx(-1, PATCH_SIZE, offs, dims)*EDGE + RT + 1), ADJZONE)
       end if
    else
       if (c == IPLUSJMINUS) then
          id_chd = idx(LAST - 1, 0, offs_chd, dims_chd)
          idE_chd = idx(LAST, 0, offs_chd, dims_chd)

          if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
             call set_at_least (dom%mask_e%elts(idx(PATCH_SIZE, 0, offs, dims)*EDGE + RT + 1), ADJZONE)
             call set_at_least (dom%mask_e%elts(idx(PATCH_SIZE,-1, offs, dims)*EDGE + UP + 1), ADJZONE)
          end if
       end if
    end if

    if (.not. c == IJMINUS) return

    id_chd  = idx(0, 0, offs_chd, dims_chd)
    idN_chd = idx(0, 1, offs_chd, dims_chd)
    idE_chd = idx(1, 0, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
       call mask_penta_corr (dom, offs, dims, offs_chd, dims_chd, 1)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
    end if
    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
       call mask_penta_corr (dom, offs, dims, offs_chd, dims_chd, 2)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
    end if
  end subroutine mask_perfect_velo_penta

  subroutine mask_penta_corr (dom, offs, dims, offs_chd, dims_chd, itype)
    implicit none
    real(8), dimension(2)          :: velo_interp_penta_corr
    type(Domain)                   :: dom
    integer, dimension(N_BDRY+1)   :: offs, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims, dims_chd
    integer                        :: itype

    if (itype==1) then
       call set_at_least (dom%mask_e%elts(idx(0, -1, offs, dims)*EDGE + UP + 1), ADJZONE) 
       call set_at_least (dom%mask_e%elts(idx(-1, -1, offs, dims)*EDGE + 1), ADJZONE)
       call set_at_least (dom%mask_e%elts(ed_idx(end_pt(1,1,UP+1), end_pt(2,1,UP+1), &
            hex_sides(:,hex_s_offs(UP+1) + 0 + 1), offs, dims) + 1), ADJZONE)
       call set_at_least (dom%mask_e%elts(ed_idx(opp_no(1,2,UP+1), opp_no(2,2,UP+1), &
            hex_sides(:,hex_s_offs(UP+1) + 1 + 1), offs, dims) + 1), ADJZONE)
    elseif (itype==2) then
       call set_at_least (dom%mask_e%elts(idx(-1, -1, offs, dims)*EDGE + 1), ADJZONE) 
       call set_at_least (dom%mask_e%elts(idx(-1, 0, offs, dims)*EDGE + RT + 1), ADJZONE)
       call set_at_least (dom%mask_e%elts(ed_idx(opp_no(1,1,RT+1), opp_no(2,1,RT+1), &
            hex_sides(:,hex_s_offs(RT+1) + 1 + 1), offs, dims) + 1), ADJZONE) 
       call set_at_least (dom%mask_e%elts(ed_idx(end_pt(1,1,RT+1), end_pt(2,1,RT+1), &
            hex_sides(:,hex_s_offs(RT+1) + 2 + 1), offs, dims) + 1), ADJZONE)
    end if
  end subroutine mask_penta_corr

  subroutine mask_remap (dom, i, j, zlev, offs, dims)
    ! Add nodes required for remapping velocities (need mass fluxes at edges)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN

    id   = idx(i,   j,   offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)
    idN  = idx(i,   j+1, offs, dims)

    if (dom%mask_e%elts(EDGE*id+RT+1) >= TRSK) then
       call set_at_least (dom%mask_n%elts(id+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idE+1), TRSK)
    end if
    if (dom%mask_e%elts(EDGE*id+DG+1) >= TRSK) then
       call set_at_least (dom%mask_n%elts(id+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    end if
    if (dom%mask_e%elts(EDGE*id+UP+1) >= TRSK) then
       call set_at_least (dom%mask_n%elts(id+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idN+1), TRSK)
    end if
  end subroutine mask_remap

  subroutine mask_trsk (dom, i, j, zlev, offs, dims)
    ! Add additional nodes and edgs required for TRISK operator stencils
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idNE, idN, idNW, idW, idSW, idS, idSE, idN2E, id2E, id2NE, id2N2E, id2N

    id     = idx(i,   j,   offs, dims)
    idE    = idx(i+1, j,   offs, dims)
    idNE   = idx(i+1, j+1, offs, dims)
    idNW   = idx(i-1, j+1, offs, dims)
    idN    = idx(i,   j+1, offs, dims)
    idW    = idx(i-1, j,   offs, dims)
    idSW   = idx(i-1, j-1, offs, dims)
    idS    = idx(i,   j-1, offs, dims)
    idSE   = idx(i+1, j-1, offs, dims)
    idN2E  = idx(i+2, j+1, offs, dims)
    id2E   = idx(i+2, j,   offs, dims)
    id2N   = idx(i,   j+2, offs, dims)
    id2NE  = idx(i+1, j+2, offs, dims)
    id2N2E = idx(i+2, j+2, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       ! Flux divergence stencil
       call set_at_least (dom%mask_n%elts(idE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
       call set_at_least (dom%mask_n%elts(idN+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), TRSK) 

       ! Kinetic energy and circulation stencil
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)

       ! Circulation stencil
       call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+DG+1), TRSK)

       ! Potential vorticity stencil
       call set_at_least (dom%mask_n%elts(idW+1), TRSK)
       call set_at_least (dom%mask_n%elts(idS+1), TRSK)

       ! Stencil for gradients of Bernoulli function
       call set_at_least (dom%mask_n%elts(id+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idE+1), TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idE+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idE+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSE+UP+1), TRSK)

       call set_at_least (dom%mask_n%elts(id+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1), TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1), TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idNE+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNE+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNE+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1), TRSK)

       call set_at_least (dom%mask_n%elts(id+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idN+1), TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1), TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNW+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+DG+1), TRSK)

       ! Qperp stencil (Gassmann)
       call set_at_least (dom%mask_n%elts(id2E+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idN2E+1), TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idN+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idW+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idSW+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idS+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idSE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idE+1),   TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idE+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idE+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSE+UP+1), TRSK)

       call set_at_least (dom%mask_n%elts(idE+1),    TRSK)
       call set_at_least (dom%mask_n%elts(idN2E+1),  TRSK)
       call set_at_least (dom%mask_n%elts(id2N2E+1), TRSK)
       call set_at_least (dom%mask_n%elts(id2NE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idN+1),    TRSK)
       call set_at_least (dom%mask_n%elts(idW+1),    TRSK)
       call set_at_least (dom%mask_n%elts(idSW+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idS+1),    TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1),   TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNE+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNE+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNE+UP+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+UP+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)

       call set_at_least (dom%mask_n%elts(idE+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idNE+1),  TRSK)
       call set_at_least (dom%mask_n%elts(id2NE+1), TRSK)
       call set_at_least (dom%mask_n%elts(id2N+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idNW+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idW+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idSW+1),  TRSK)
       call set_at_least (dom%mask_n%elts(idS+1),   TRSK)
       call set_at_least (dom%mask_n%elts(idN+1),   TRSK)

       call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idN+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idNW+RT+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+DG+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1),  TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+RT+1),   TRSK)
       call set_at_least (dom%mask_e%elts(EDGE*id+DG+1),   TRSK)
    end if
  end subroutine mask_trsk

  subroutine mask_adj_parent_nodes (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent node to adjacent zone if any child is active, also make child active if any child neighbours are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idN, idE, idNE, idSW, idS, idW

    id_par = idx(i_par, j_par, offs_par, dims_par)
    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)

    idN    = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idE    = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE   = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    idSW   = idx(i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS    = idx(i_chd,   j_chd-1, offs_chd, dims_chd)
    idW    = idx(i_chd-1, j_chd,   offs_chd, dims_chd)

    if (dom%mask_n%elts(idE+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idN+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idW+1)  == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idS+1)  == TOLRNZ) then

       call set_at_least (dom%mask_n%elts(id_par+1), TOLRNZ)
       call set_at_least (dom%mask_n%elts(id_chd+1), TOLRNZ)
    end if
  end subroutine mask_adj_parent_nodes

  subroutine mask_adj_parent_edges (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Nearest neighbours of active edges at coarser scale
    !
    ! Add parent edge to active mask if at least one of the four neighbouring child edges is active
    ! (check two additional child neighbour edges for each parent edge compared to wavetrisk)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: e, id_chd,  id_par, idE, idN, idNE, idNW, idS, idSE, idW

    id_par = idx(i_par, j_par, offs_par, dims_par) ! parent node
    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd) ! child node

    ! Three neighbours of child node
    idE  = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idW  = idx(i_chd-1, j_chd,   offs_chd, dims_chd)
    idN  = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idNE = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    idNW = idx(i_chd-1, j_chd+1, offs_chd, dims_chd)
    idS  = idx(i_chd,   j_chd-1, offs_chd, dims_chd)
    idSE = idx(i_chd+1, j_chd-1, offs_chd, dims_chd)

    ! Check six child edges neighbouring each parent edge to see if at least one is active
    ! if true make parent edge active
    if (dom%mask_e%elts(EDGE*id_chd+RT+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idSE+UP+1)    == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), TOLRNZ)

    if (dom%mask_e%elts(EDGE*id_chd+DG+1)   == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)      == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1)      == TOLRNZ) &
         call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), TOLRNZ)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1)  == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+UP+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+DG+1)     == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNW+RT+1)    == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1)     == TOLRNZ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), TOLRNZ)
  end subroutine mask_adj_parent_edges

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
    id = idx(i, j, offs, dims)

    ! six nearest neighbour nodes at same scale
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

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
    id  = idx(i, j, offs, dims)

    ! six neighbour nodes at same scale
    idS  = idx(i,   j-1, offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idSW = idx(i-1, j-1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

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
    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)

    idE_chd   = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd   = idx(i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW_chd  = idx(i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSW_chd  = idx(i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS_chd   = idx(i_chd,   j_chd-1, offs_chd, dims_chd)
    idSE_chd  = idx(i_chd+1, j_chd-1, offs_chd, dims_chd)
    idN2E_chd = idx(i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE_chd = idx(i_chd+1, j_chd+2, offs_chd, dims_chd)
    id2SW_chd = idx(i_chd-1, j_chd-2, offs_chd, dims_chd)
    idS2W_chd = idx(i_chd-2, j_chd-1, offs_chd, dims_chd)

    id_par   = idx(i_par,   j_par,   offs_par, dims_par)

    idS_par  = idx(i_par,   j_par-1, offs_par, dims_par)
    idW_par  = idx(i_par-1, j_par,   offs_par, dims_par)
    idSW_par = idx(i_par-1, j_par-1, offs_par, dims_par)
    idE_par  = idx(i_par+1, j_par,   offs_par, dims_par)
    idN_par  = idx(i_par,   j_par+1, offs_par, dims_par)
    idNE_par = idx(i_par+1, j_par+1, offs_par, dims_par)

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
    ! Add neigbhour edges of active nodes and vice versa (ensure consistency of node and edge adjacent zones)
    implicit none
    integer :: l

    do l = level_end-1, level_start, -1
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
  end subroutine mask_adj_nodes_edges

  subroutine prolong_node_adjzone (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent to adjacent mask if child is in adjacent zone for prolongation of nodal quantities in inverse wavelet transform
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par

    id_par = idx(i_par, j_par, offs_par, dims_par) ! id of parent
    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd) ! id of child

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

    id_par   = idx(i_par,   j_par,   offs_par, dims_par)
    idE_par  = idx(i_par+1, j_par,   offs_par, dims_par)
    idN_par  = idx(i_par,   j_par+1, offs_par, dims_par)
    idNE_par = idx(i_par+1, j_par+1, offs_par, dims_par)

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

  subroutine mask_e_consist (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id, id_par, idN, idS, idE, idW, idNE, idNW, idSE

    id     = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    idN  = idx(i_chd,   j_chd+1, offs_chd, dims_chd)
    idE  = idx(i_chd+1, j_chd,   offs_chd, dims_chd)
    idS  = idx(i_chd,   j_chd-1, offs_chd, dims_chd)
    idW  = idx(i_chd-1, j_chd,   offs_chd, dims_chd)
    idNE = idx(i_chd+1, j_chd+1, offs_chd, dims_chd)
    idNW = idx(i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSE = idx(i_chd+1, j_chd-1, offs_chd, dims_chd)

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

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    ! Eight nearest neighbours of child node
    idN  = idx(i_chd,   j_chd+2, offs_chd, dims_chd)
    idE  = idx(i_chd+2, j_chd,   offs_chd, dims_chd)
    idNE = idx(i_chd+2, j_chd+2, offs_chd, dims_chd)
    idS  = idx(i_chd,   j_chd-2, offs_chd, dims_chd)
    idW  = idx(i_chd-2, j_chd,   offs_chd, dims_chd)
    idSW = idx(i_chd-2, j_chd-2, offs_chd, dims_chd)
    idNW = idx(i_chd-2, j_chd+2, offs_chd, dims_chd)
    idSE = idx(i_chd+2, j_chd-2, offs_chd, dims_chd)

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

  subroutine mask_n_if_all_e (dom, i, j, zlev, offs, dims)
    ! Add node to adjacent zone if all neighbouring edges are active
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idW, idS, idSW

    id   = idx(i,   j,   offs, dims)
    idW  = idx(i-1, j,   offs, dims)
    idS  = idx(i,   j-1, offs, dims)
    idSW = idx(i-1, j-1, offs, dims)

    if ( dom%mask_e%elts(EDGE*id+UP+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+DG+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id+RT+1)   >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idS+UP+1)  >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idSW+DG+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idW+RT+1)  >= ADJZONE) then

       call set_at_least (dom%mask_n%elts(id+1), ADJZONE)
    end if
  end subroutine mask_n_if_all_e

  subroutine mask_e_if_both_n (dom, i, j, zlev, offs, dims)
    ! Add edge to adjacent zone if both neighbouring nodes are active
    ! (also needed for div and gradv_e operator)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, idE, idN, idNE

    id   = idx(i,   j,     offs, dims)
    idE  = idx(i+1, j,     offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

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

  subroutine init_mask_mod
    implicit none
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_domain_mod
    call init_comm_mod
    tol_velo = 0.0_8
    tol_mass = 0.0_8
    tol_temp = 0.0_8
    initialized = .True.
  end subroutine init_mask_mod

  subroutine init_masks
    implicit none
    integer :: d, i, num

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

  subroutine set_at_least (mask, typ)
    implicit none
    integer :: mask, typ

    if (mask < typ) mask = typ
  end subroutine set_at_least
end module mask_mod
