module mask_mod
  use param_mod
  use shared_mod
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use wavelet_mod
  implicit none
  real(8) toll_velo
  real(8) toll_height

contains
  subroutine init_mask_mod()
      logical :: initialized = .False.
      if (initialized) return ! initialize only once
      call init_domain_mod()
      call init_comm_mod()
      toll_velo = 0.0_8
      toll_height = 0.0_8
      initialized = .True.
  end subroutine init_mask_mod

  subroutine mask_u_consist(dom, i_par, j_par, i_chd, j_chd, offs_par, &
          dims_par, offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,9) :: dims_par
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,9) :: dims_chd
      integer id
      integer idN, idS
      integer idE, idW
      integer idNE, idNW, idSE
      integer id_par
      id = idx(i_chd, j_chd, offs_chd, dims_chd)
      idN = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
      idE = idx(i_chd + 1, j_chd, offs_chd, dims_chd)
      idS = idx(i_chd, j_chd - 1, offs_chd, dims_chd)
      idW = idx(i_chd - 1, j_chd, offs_chd, dims_chd)
      idNE = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
      idNW = idx(i_chd - 1, j_chd + 1, offs_chd, dims_chd)
      idSE = idx(i_chd + 1, j_chd - 1, offs_chd, dims_chd)
      id_par = idx(i_par, j_par, offs_par, dims_par)
      if (dom%mask_u%elts(EDGE*id+UP+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idN+UP+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idW+DG+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idN+DG+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idNW+RT+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .ge. CONSIST) then
          call set_at_least(dom%mask_u%elts(EDGE*idN+UP+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*id+UP+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), ADJZONE)
      end if
      if (dom%mask_u%elts(DG+EDGE*id+1) .ge. CONSIST .or. &
              dom%mask_u%elts(DG+EDGE*idNE+1) .ge. CONSIST .or. &
              dom%mask_u%elts(UP+EDGE*idE+1) .ge. CONSIST .or. &
              dom%mask_u%elts(UP+EDGE*idNE+1) .ge. CONSIST .or. &
              dom%mask_u%elts(RT+EDGE*idN+1) .ge. CONSIST .or. &
              dom%mask_u%elts(RT+EDGE*idNE+1) .ge. CONSIST) then
          call set_at_least(dom%mask_u%elts(DG+EDGE*idNE+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(DG+EDGE*id+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(DG+EDGE*id_par+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*id+RT+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idE+RT+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idSE+UP+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idS+DG+1) .ge. CONSIST .or. &
              dom%mask_u%elts(EDGE*idE+DG+1) .ge. CONSIST) then
          call set_at_least(dom%mask_u%elts(EDGE*idE+RT+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*id+RT+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), ADJZONE)
      end if
  end subroutine

  subroutine mask_adj_space(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id
      integer idS
      integer idW
      integer idSW
      integer idE
      integer idN
      integer idNE
      id = idx(i, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      if (dom%mask_p%elts(idN+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idNE+1) .eq. &
              TOLLRNZ .or. dom%mask_p%elts(idE+1) .eq. TOLLRNZ .or. &
              dom%mask_p%elts(idS+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idSW+1) &
              .eq. TOLLRNZ .or. dom%mask_p%elts(idW+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(id+1), ADJZONE)
      end if
  end subroutine

  subroutine mask_adj_space2(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id
      integer idS
      integer idW
      integer idSW
      integer idE
      integer idN
      integer idNE
      id = idx(i, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      if (dom%mask_p%elts(idN+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idNE+1) .eq. &
              TOLLRNZ .or. dom%mask_p%elts(idE+1) .eq. TOLLRNZ .or. &
              dom%mask_p%elts(idS+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idSW+1) &
              .eq. TOLLRNZ .or. dom%mask_p%elts(idW+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(id+1), ADJSPACE)
      end if
      if (dom%mask_u%elts(DG+EDGE*id+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idW+RT+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*idW+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*id+UP+1), ADJSPACE)
      end if
      if (dom%mask_u%elts(EDGE*id+UP+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*id+RT+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*id+DG+1), ADJSPACE)
      end if
      if (dom%mask_u%elts(DG+EDGE*id+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idS+UP+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*idS+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*id+RT+1), ADJSPACE)
      end if
  end subroutine

  subroutine mask_active_velo(dom, i_par, j_par, i_chd, j_chd, offs_par, &
          dims_par, offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer id_par
      integer id, idE, idN, idNE
      integer e
      id_par = idx(i_par, j_par, offs_par, dims_par)
      id = idx(i_chd, j_chd, offs_chd, dims_chd)
      idE = idx(i_chd + 1, j_chd, offs_chd, dims_chd)
      idN = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
      idNE = idx(i_chd + 1, j_chd+ 1, offs_chd, dims_chd)

      if (dom%mask_u%elts(EDGE*id+RT+1) .eq. TOLLRNZ .or. dom%mask_u%elts(EDGE*idE+RT+1) .eq. TOLLRNZ) &
              call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), TOLLRNZ)
      if (dom%mask_u%elts(EDGE*id+DG+1) .eq. TOLLRNZ .or. dom%mask_u%elts(EDGE*idNE+DG+1) .eq. TOLLRNZ) &
              call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), TOLLRNZ)
      if (dom%mask_u%elts(EDGE*id+UP+1) .eq. TOLLRNZ .or. dom%mask_u%elts(EDGE*idN+UP+1) .eq. TOLLRNZ) &
              call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), TOLLRNZ)

      if (dom%mask_u%elts(EDGE*idNE+RT+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idNE+UP+1) .eq. TOLLRNZ) &
          call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), TOLLRNZ)
      if (dom%mask_u%elts(EDGE*idE+DG+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .eq. TOLLRNZ) &
          call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), TOLLRNZ)
      if (dom%mask_u%elts(EDGE*idN+RT+1) .eq. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN+DG+1) .eq. TOLLRNZ) &
          call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), TOLLRNZ)
  end subroutine

  subroutine init_masks()
      integer d
      integer num
      integer i
      do d = 1, size(grid)
          call init(grid(d)%mask_p, 1)
          call init(grid(d)%mask_u, EDGE)
          call init(grid(d)%level, 1)
      end do
      do d = 1, size(grid)
          num = grid(d)%node%length-1
          call extend(grid(d)%mask_p, num, TOLLRNZ)
          call extend(grid(d)%mask_u, EDGE*num, TOLLRNZ)
          call extend(grid(d)%level, num, min_level-1)
      end do
  end subroutine

  subroutine inj_p_adjzone(dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par, &
          offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,9) :: dims_par
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,9) :: dims_chd
      integer id_par
      integer id
      id_par = idx(i_par, j_par, offs_par, dims_par)
      id = idx(i_chd, j_chd, offs_chd, dims_chd)
      if (dom%mask_p%elts(id+1) .ge. ADJZONE) then
          call set_at_least(dom%mask_p%elts(id_par+1), ADJZONE)
      end if
  end subroutine

  subroutine set_masks(dom, p, i, j, offs, dims, mask)
      type(Domain) dom
      integer p
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer mask
      integer id
      integer e
      id = idx(i, j, offs, dims)
      if (dom%mask_p%elts(id+1) .eq. FROZEN) return
      dom%mask_p%elts(id+1) = mask
      do e = 1, EDGE
          dom%mask_u%elts(EDGE*id+e) = mask
      end do
  end subroutine

  subroutine mask_active()
      integer l
      call update_bdry1(wav_coeff(S_HEIGHT), level_start, level_end)
      call update_bdry1(wav_coeff(S_VELO), level_start, level_end)
      call apply_onescale(mask_toll, level_end, -1, 2)
      do l = level_end-1, level_start, -1
          call apply_onescale(mask_toll, l, -1, 2)
          call apply_interscale(mask_active_height, l, 0, 1)
          call apply_interscale(mask_active_velo, l, -1, 1)
          call comm_masks_mpi(l)
      end do
  end subroutine

  subroutine mask_toll(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      integer id
      integer e
      id = idx(i, j, offs, dims)
      if (dom%mask_p%elts(id+1) .eq. FROZEN) return
      call set_active_mask(dom%mask_p%elts(id+1), wav_coeff(S_HEIGHT)%data(dom%id+1)%elts(id+1), &
                           toll_height)
      do e = 1, EDGE
          call set_active_mask(dom%mask_u%elts(EDGE*id+e), &
                               wav_coeff(S_VELO)%data(dom%id+1)%elts(EDGE*id+e), &
                               toll_velo)
      end do
  end subroutine

  subroutine set_active_mask(mask, wc, toll)
      integer, intent(inout) :: mask
      real(8), intent(in) :: wc, toll
      if (abs(wc) .ge. toll) then
          mask = TOLLRNZ
      else
          if (mask .gt. ADJZONE) &
              mask = ADJZONE
      end if
  end subroutine

  subroutine mask_restrict_flux(dom, i_par, j_par, offs_par, dims_par)
      type(Domain) dom
      integer i_par, j_par
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer id_par, idE_par, idN_par, idNE_par
      id_par = idx(i_par, j_par, offs_par, dims_par)
      idE_par = idx(i_par + 1, j_par, offs_par, dims_par)
      idN_par = idx(i_par, j_par + 1, offs_par, dims_par)
      idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)
      if (dom%mask_p%elts(id_par+1) .ge. ADJSPACE) then
          call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), RESTRCT)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), RESTRCT)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), RESTRCT)
      else
          if (dom%mask_p%elts(idN_par+1) .ge. ADJSPACE) then
              call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), RESTRCT)
          end if
          if (dom%mask_p%elts(idNE_par+1) .ge. ADJSPACE) then
              call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), RESTRCT)
          end if
          if (dom%mask_p%elts(idE_par+1) .ge. ADJSPACE) then
              call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), RESTRCT)
          end if
      end if
  end subroutine

  subroutine mask_adj_scale(dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par, &
          offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
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
      id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
      idN_chd = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
      idE_chd = idx(i_chd + 1, j_chd, offs_chd, dims_chd)
      idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
      id_par = idx(i_par, j_par, offs_par, dims_par)
      idS_par = idx(i_par, j_par - 1, offs_par, dims_par)
      idW_par = idx(i_par - 1, j_par, offs_par, dims_par)
      idSW_par = idx(i_par - 1, j_par - 1, offs_par, dims_par)
      idE_par = idx(i_par + 1, j_par, offs_par, dims_par)
      idN_par = idx(i_par, j_par + 1, offs_par, dims_par)
      idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)
      if (dom%mask_p%elts(id_par+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(id_chd+1), ADJZONE)
          call set_at_least(dom%mask_p%elts(idN_chd+1), ADJZONE)
          call set_at_least(dom%mask_p%elts(idNE_chd+1), ADJZONE)
          call set_at_least(dom%mask_p%elts(idE_chd+1), ADJZONE)

          call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), RESTRCT)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), RESTRCT)
          call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), RESTRCT)
      else
          if (dom%mask_p%elts(idN_par+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(idN_chd+1), ADJZONE)
              call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), RESTRCT)
          end if
          if (dom%mask_p%elts(idNE_par+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(idNE_chd+1), ADJZONE)
              call set_at_least(dom%mask_u%elts(EDGE*id_par+DG+1), RESTRCT)
          end if
          if (dom%mask_p%elts(idE_par+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(idE_chd+1), ADJZONE)
              call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), RESTRCT)
          end if
      end if
      if (dom%mask_u%elts(EDGE*id_par+UP+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(idN_chd+1), ADJZONE)
          if (dom%mask_u%elts(EDGE*idS_par+UP+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(id_chd+1), ADJZONE)
          end if
      end if
      if (dom%mask_u%elts(EDGE*id_par+UP+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idW_par+RT+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*idW_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN_par+RT+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*id_chd+UP+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*idN_chd+UP+1), ADJZONE)
      end if
      if (dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(idNE_chd+1), ADJZONE)
          if (dom%mask_u%elts(DG+EDGE*idSW_par+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(id_chd+1), ADJZONE)
          end if
      end if
      if (dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*id_par+UP+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*id_par+RT+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE_par+UP+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN_par+RT+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(DG+EDGE*id_chd+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(DG+EDGE*idNE_chd+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*id_par+RT+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(idE_chd+1), ADJZONE)
          if (dom%mask_u%elts(EDGE*idW_par+RT+1) .ge. TOLLRNZ) then
              call set_at_least(dom%mask_p%elts(id_chd+1), ADJZONE)
          end if
      end if
      if (dom%mask_u%elts(EDGE*id_par+RT+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idS_par+UP+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*idS_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE_par+UP+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*id_chd+RT+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*idE_chd+RT+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*id_par+UP+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idN_par+RT+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*idN_chd+RT+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(DG+EDGE*idN_chd+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*idNE_chd+UP+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*id_par+RT+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(DG+EDGE*id_par+1) .ge. TOLLRNZ .or. &
              dom%mask_u%elts(EDGE*idE_par+UP+1) .ge. TOLLRNZ) then
          call set_at_least(dom%mask_u%elts(EDGE*idE_chd+UP+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(DG+EDGE*idE_chd+1), ADJZONE)
          call set_at_least(dom%mask_u%elts(EDGE*idNE_chd+RT+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*id_par+RT+1) .ge. ADJSPACE .or. &
              dom%mask_u%elts(EDGE*id_par+DG+1) .ge. ADJSPACE .or. &
              dom%mask_u%elts(EDGE*id_par+UP+1) .ge. ADJSPACE .or. &
              dom%mask_u%elts(EDGE*idW_par+RT+1) .ge. ADJSPACE .or. &
              dom%mask_u%elts(EDGE*idSW_par+DG+1) .ge. ADJSPACE .or. &
              dom%mask_u%elts(EDGE*idS_par+UP+1) .ge. ADJSPACE) then
          call set_at_least(dom%mask_p%elts(id_par+1), RESTRCT)
      end if
  end subroutine

  subroutine set_at_least(mask, typ)
      integer mask
      integer typ
      if (mask .lt. typ) then
          mask = typ
      end if
  end subroutine

  subroutine mask_u_consist2(dom, i_par, j_par, i_chd, j_chd, offs_par, &
          dims_par, offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,9) :: dims_par
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,9) :: dims_chd
      integer id
      integer idN, idE
      integer idS, idW
      integer idNE, idSW
      integer idNW, idSE
      integer id_par
      id = idx(i_chd, j_chd, offs_chd, dims_chd)
      idN = idx(i_chd, j_chd + 2, offs_chd, dims_chd)
      idE = idx(i_chd + 2, j_chd, offs_chd, dims_chd)
      idNE = idx(i_chd + 2, j_chd + 2, offs_chd, dims_chd)
      idS = idx(i_chd, j_chd - 2, offs_chd, dims_chd)
      idW = idx(i_chd - 2, j_chd, offs_chd, dims_chd)
      idSW = idx(i_chd - 2, j_chd - 2, offs_chd, dims_chd)
      idNW = idx(i_chd - 2, j_chd + 2, offs_chd, dims_chd)
      idSE = idx(i_chd + 2, j_chd - 2, offs_chd, dims_chd)
      id_par = idx(i_par, j_par, offs_par, dims_par)
      if (dom%mask_u%elts(EDGE*id+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*id+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idN+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idNE+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idW+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idNW+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idSW+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idW+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idSW+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idW+UP+1) .ge. ADJZONE) then
          call set_at_least(dom%mask_u%elts(EDGE*id_par+UP+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*idW+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*id+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idNE+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idS+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idW+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idE+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idN+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idS+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*id+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idNE+UP+1) .ge. ADJZONE) then
          call set_at_least(dom%mask_u%elts(DG+EDGE*id_par+1), ADJZONE)
      end if
      if (dom%mask_u%elts(EDGE*idSW+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idS+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idN+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idNE+RT+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idSW+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idS+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*id+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(DG+EDGE*idE+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idS+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idSE+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*id+UP+1) .ge. ADJZONE .or. &
              dom%mask_u%elts(EDGE*idE+UP+1) .ge. ADJZONE) then
          call set_at_least(dom%mask_u%elts(EDGE*id_par+RT+1), ADJZONE)
      end if
  end subroutine

  subroutine mask_active_height(dom, i_par, j_par, i_chd, j_chd, offs_par, &
          dims_par, offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
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
      id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
      idN = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
      idE = idx(i_chd + 1, j_chd, offs_chd, dims_chd)
      idNE = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
      idSW = idx(i_chd - 1, j_chd - 1, offs_chd, dims_chd)
      idS = idx(i_chd, j_chd - 1, offs_chd, dims_chd)
      idW = idx(i_chd - 1, j_chd, offs_chd, dims_chd)
      if (dom%mask_p%elts(idE+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idNE+1) .eq. &
              TOLLRNZ .or. dom%mask_p%elts(idN+1) .eq. TOLLRNZ .or. &
              dom%mask_p%elts(idW+1) .eq. TOLLRNZ .or. dom%mask_p%elts(idSW+1) &
              .eq. TOLLRNZ .or. dom%mask_p%elts(idS+1) .eq. TOLLRNZ) then
          call set_at_least(dom%mask_p%elts(id_par+1), TOLLRNZ)
          call set_at_least(dom%mask_p%elts(id_chd+1), TOLLRNZ)
      end if
  end subroutine

  subroutine mask_p_if_all_u(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id
      integer idW
      integer idS
      integer idSW
      id = idx(i, j, offs, dims)
      idW = idx(i - 1, j, offs, dims)
      idS = idx(i, j - 1, offs, dims)
      idSW = idx(i - 1, j - 1, offs, dims)
      if (dom%mask_u%elts(EDGE*id+UP+1) .ge. ADJZONE .and. &
              dom%mask_u%elts(EDGE*id+DG+1) .ge. ADJZONE .and. &
              dom%mask_u%elts(EDGE*id+RT+1) .ge. ADJZONE .and. &
              dom%mask_u%elts(EDGE*idS+UP+1) .ge. ADJZONE .and. &
              dom%mask_u%elts(EDGE*idSW+DG+1) .ge. ADJZONE .and. &
              dom%mask_u%elts(EDGE*idW+RT+1) .ge. ADJZONE) &
          call set_at_least(dom%mask_p%elts(id+1), ADJZONE)
  end subroutine

  subroutine mask_u_if_both_p(dom, i, j, offs, dims)
      type(Domain) dom
      integer i
      integer j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,9) :: dims
      integer id
      integer idE
      integer idN
      integer idNE
      id = idx(i, j, offs, dims)
      idE = idx(i + 1, j, offs, dims)
      idN = idx(i, j + 1, offs, dims)
      idNE = idx(i + 1, j + 1, offs, dims)
      if (dom%mask_p%elts(id+1) .ge. ADJZONE .and. dom%mask_p%elts(idE+1) .ge. ADJZONE) &
         call set_at_least(dom%mask_u%elts(id*EDGE+RT+1), ADJZONE)
      if (dom%mask_p%elts(idNE+1) .ge. ADJZONE .and. dom%mask_p%elts(id+1) .ge. ADJZONE) &
         call set_at_least(dom%mask_u%elts(id*EDGE+DG+1), ADJZONE)
      if (dom%mask_p%elts(id+1) .ge. ADJZONE .and. dom%mask_p%elts(idN+1) .ge. ADJZONE) &
         call set_at_least(dom%mask_u%elts(id*EDGE+UP+1), ADJZONE)
  end subroutine

  subroutine complete_masks()
      integer l
      call apply_onescale(mask_adj_space, level_end, 0, 1)
      do l = level_end-1, level_start, -1
          call apply_interscale(mask_adj_scale, l, 0, 1)
          call apply_onescale(mask_u_if_both_p, l+1, 0, 0)
      end do
      call comm_masks_mpi(NONE)
      do l = level_end-1, level_start+1, -1
          call apply_interscale(mask_u_consist, l, 0, 1)
          call comm_masks_mpi(l+1)
          call apply_interscale(mask_u_consist2, l, 0, 0)
          call comm_masks_mpi(l)
      end do
      if (level_start .lt. level_end) then
          call apply_interscale(mask_u_consist, level_start, 0, 1)
          call comm_masks_mpi(level_start+1)
      end if
      do l = level_end-1, level_start+1, -1
          call apply_onescale(mask_p_if_all_u, l+1, 0, 1)
          call apply_interscale(inj_p_adjzone, l, 0, 1)
      end do
      if (level_start+1 .le. level_end) call apply_onescale(mask_p_if_all_u, level_start+1, 0, 1)
  end subroutine
end module mask_mod
