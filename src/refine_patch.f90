module refine_patch_mod
  use shared_mod
  use domain_mod
  use init_mod
  use wavelet_mod
  use mask_mod
  implicit none

contains
  subroutine init_refine_patch_mod()
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_shared_mod()
    call init_domain_mod()
    call init_init_mod()
    call init_wavelet_mod()
    call init_mask_mod()
    initialized = .True.
  end subroutine init_refine_patch_mod

  subroutine attach_bdry(dom, p_par, c, s, side)
    type(Domain) dom
    integer p_par
    integer c
    integer s
    integer side
    integer n_chd
    integer p_chd

    n_chd = find_neigh_patch_Domain(dom, p_par, c, s)
    if (n_chd .eq. 0) then
       n_chd = -add_bdry_patch_Domain(dom, side)
    end if

    p_chd = dom%patch%elts(p_par+1)%children(c+1)
    dom%patch%elts(p_chd+1)%neigh(s+1) = n_chd
  end subroutine attach_bdry

  subroutine refine_patch(dom, p, c0)
    type(Domain) dom
    integer p, c0

    call refine_patch1(dom, p, c0)
    call refine_patch2(dom, p, c0)
  end subroutine refine_patch

  subroutine refine_patch1(dom, p, c0)
    type(Domain) dom
    integer p, c0
    integer s
    integer lev
    integer c
    integer p_chd
    integer k
    integer j
    integer j_chd
    integer j_par
    integer i
    integer i_chd
    integer i_par
    integer id_par
    integer num, d
    type(Coord), dimension(6) :: tmp
    !  Main difficulty: in order to precompute geometry, weights etc
    !         nodes outside the patch are needed that might not be part of the grid yet
    !         Better compute them temporarily that adding additional patches
    !         that are empty in terms of evaluating operators and slow down the inner loops
    !         Stategy A: use temporary bdry_patch

    lev = dom%patch%elts(p+1)%level
    if (lev .eq. max_level) then
       return
    end if

    if (level_end .eq. lev) then
       level_end = level_end + 1
    end if

    num = dom%node%length 
    p_chd = add_patch_Domain(dom, lev + 1)
    c = c0 + 1

    dom%patch%elts(p+1)%children(c) = p_chd
    num = dom%node%length - num

    call extend(dom%level, num, dom%patch%elts(p_chd+1)%level)
  end subroutine refine_patch1

  subroutine refine_patch2(dom, p, c0)
    type(Domain) dom
    integer p, c0
    integer s
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer lev
    integer c
    integer p_chd
    integer k
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_chd
    integer j
    integer j_chd
    integer j_par
    integer i
    integer i_chd
    integer i_par
    integer id_par
    integer num, d
    type(Coord), dimension(6) :: tmp

    call get_offs_Domain(dom, p, offs_par, dims_par)

    lev = dom%patch%elts(p+1)%level
    c = c0 + 1
    p_chd = dom%patch%elts(p+1)%children(c)

    do k = 1, 2
       call attach_bdry(dom, p, c - 1, modulo(c + k - 2, N_CHDRN), side(dom, p, modulo(c + k - 2, N_CHDRN)))
       call attach_bdry(dom, p, c - 1, modulo(c + k, N_CHDRN), -(modulo(c + k, N_CHDRN) + 1))
       call connect_cousin(dom, p, p_chd, modulo(c + k - 2, N_CHDRN), modulo(c + k - 2, N_CHDRN), modulo(c - 2*k + 2, N_CHDRN))
    end do

    call attach_bdry(dom, p, c - 1, c + 3, side(dom, p, c + 3))

    call connect_cousin(dom, p, p_chd, c + 3, c + 3, modulo(c + 1, N_CHDRN))

    call attach_bdry(dom, p, c - 1, modulo(c + 1, N_CHDRN) + 4, -(modulo(c + 1, N_CHDRN) + 4 + 1))
    call attach_bdry(dom, p, c - 1, modulo(c,     N_CHDRN) + 4, side(dom, p, modulo(c, N_CHDRN)))
    call attach_bdry(dom, p, c - 1, modulo(c + 2, N_CHDRN) + 4, side(dom, p, modulo(c - 1, N_CHDRN)))

    call get_offs_Domain(dom, p_chd, offs_chd, dims_chd)

    do j = 0, PATCH_SIZE/2 + 1
       j_chd = (j - 1)*2
       j_par = j - 1 + chd_offs(2,c)
       do i = 0, PATCH_SIZE/2 + 1
          i_chd = (i - 1)*2
          i_par = i - 1 + chd_offs(1,c)
          id_par = idx(i_par, j_par, offs_par, dims_par)

          dom%node%elts(idx(i_chd, j_chd, offs_chd, dims_chd) + 1) = dom%node%elts(id_par+1)

          dom%node%elts(idx(i_chd + 1, j_chd,     offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+RT+1)
          dom%node%elts(idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+DG+1)
          dom%node%elts(idx(i_chd,     j_chd + 1, offs_chd, dims_chd) + 1) = dom%midpt%elts(EDGE*id_par+UP+1)
       end do
    end do

    if (is_penta(dom, p_chd, IPLUSJMINUS-1)) then
       dom%node%elts(idx(PATCH_SIZE, -1, offs_chd, dims_chd) + 1) = mid_pt( &
            dom%node%elts(idx(PATCH_SIZE+1, 0, offs_par, dims_par)+1), &
            dom%node%elts(idx(PATCH_SIZE,   0, offs_par, dims_par)+1))
    end if

    if (is_penta(dom, p_chd, IMINUSJPLUS-1)) then
       dom%node%elts(idx(-1, PATCH_SIZE, offs_chd, dims_chd) + 1) = mid_pt( &
            dom%node%elts(idx(0, PATCH_SIZE+1, offs_par, dims_par)+1), &
            dom%node%elts(idx(0, PATCH_SIZE,   offs_par, dims_par)+1))
    end if

    num = dom%node%length - dom%areas%length
    d = dom%id + 1

    call extend(dom%ccentre, TRIAG*num, ORIGIN)
    call apply_onescale_to_patch2(ccentre, dom, p_chd, z_null, -2, 1)
    call ccentre_penta(dom, p_chd)
    call extend(dom%midpt, EDGE*num, ORIGIN)
    call apply_onescale_to_patch2(midpt, dom, p_chd, z_null, -1, 2)
    call extend(dom%pedlen, EDGE*num, 0.0_8)
    call extend(dom%len, EDGE*num, 0.0_8)
    call apply_onescale_to_patch2(lengths, dom, p_chd, z_null, -1, 2)

    tmp = ORIGIN
    call extend(dom%areas, num, Areas(0.0_8, 0.0_8))
    call apply_onescale_to_patch2(cpt_areas, dom, p_chd, z_null, -1, 2)
    call extend(dom%triarea, EDGE*num, 1.0_8)
    call apply_onescale_to_patch(cpt_triarea, dom, p_chd, z_null, -1, 1)
    call extend(dom%coriolis, TRIAG*num, 0.0_8)
    call apply_onescale_to_patch(coriolis, dom, p_chd, z_null, -1, 1)
    call extend(dom%windstress, EDGE*num, 0.0_8)
    call extend(dom%bernoulli, num, 0.0_8)
    call extend(dom%exner, num, 0.0_8)
    call extend(dom%surf_press, num, 0.0_8)
    call extend(dom%press, num, 0.0_8)
    call extend(dom%surf_geopot, num, 0.0_8)
    call extend(dom%geopot, num, 0.0_8)
    call extend(dom%spec_vol, num, 0.0_8)
    call extend(dom%adj_mass, num, 0.0_8)
    call extend(dom%adj_temp, num, 0.0_8)
    call extend(dom%adj_spec_vol, num, 0.0_8)
    call extend(dom%kin_energy, num, 0.0_8)

    do k = 1, zlevels
       call extend(trend(S_MASS,k)%data(d), num, 0.0_8)
       call extend(trend(S_VELO,k)%data(d), num*EDGE, 0.0_8)
       call extend(trend(S_TEMP,k)%data(d), num, 0.0_8)
       call extend(horiz_massflux(k)%data(d), num*EDGE, 0.0_8)
       call extend(horiz_tempflux(k)%data(d), num*EDGE, 0.0_8)
       call extend(wav_coeff(S_MASS,k)%data(d), num, 0.0_8)
       call extend(wav_coeff(S_VELO,k)%data(d), num*EDGE, 0.0_8)
       call extend(wav_coeff(S_TEMP,k)%data(d), num, 0.0_8)
    end do

    if (penalize) call extend(penal%data(d), num, 1.0_8)
    call extend(dom%qe, EDGE*num, 0.0_8)
    call extend(dom%vort, TRIAG*num, 0.0_8)
    call extend(dom%divu, num, 0.0_8)
    call extend(dom%overl_areas, EDGE*num, Overl_Area(0.0_8, 0.0_8))
    call extend(dom%I_u_wgt, EDGE*num, Iu_Wgt(0.0_8))
    call extend(dom%R_F_wgt, num, RF_Wgt(0.0_8))
    call extend(dom%mask_n, num, 0)
    call extend(dom%mask_e, EDGE*num, 0)

    call apply_interscale_to_patch3(set_WT_wgts, dom, p, c, z_null, 0, 0)
    call apply_interscale_to_patch3(set_RF_wgts, dom, p, c, z_null, 0, 0)

    num = dom%node%length - dom%level%length
    call extend(dom%level, num, dom%patch%elts(p_chd+1)%level)
  end subroutine refine_patch2

  integer function side(dom, p, s)
    type(Domain) dom
    integer p
    integer s
    integer n

    n = dom%patch%elts(p+1)%neigh(s+1)
    if (n .ge. 0) then
       side = -s - 1
       return
    else
       side = dom%bdry_patch%elts(-n+1)%side
       return
    end if
  end function side

  subroutine connect_cousin(dom, p_par, p_chd, s_par, s_chd, c)
    type(Domain) dom
    integer p_par
    integer p_chd
    integer s_par
    integer s_chd
    integer c
    integer n
    integer typ
    !  c: which child on neighbour

    n = dom%patch%elts(p_par+1)%neigh(s_par+1)
    if (n .lt. 0) then
       n = -n
       typ = dom%bdry_patch%elts(n+1)%side
       if (.not.s_chd .eq. typ - 1) then
          return
       end if
       if (dom%neigh(typ) .eq. POLE) then
          call connect_pole(dom, n, p_chd, s_par)
          return
       end if
       call append(dom%send_pa_all, n)
       call append(dom%send_pa_all, c)
       call append(dom%send_pa_all, p_chd)
       call append(dom%send_pa_all, s_par)
    end if
  end subroutine connect_cousin

  subroutine connect_pole(dom, n, p_chd, s_par)
    type(Domain) dom
    integer n
    integer p_chd
    integer s_par
    integer i
    do i = 1, 2
       call append(dom%send_pa_all, n)
       call append(dom%send_pa_all, i - 1)
       call append(dom%send_pa_all, p_chd)
       call append(dom%send_pa_all, s_par)
    end do
  end subroutine connect_pole

  subroutine connect_children(dom, p_par)
    type(Domain) dom
    integer p_par
    integer, dimension(N_CHDRN) :: children
    integer c
    integer p_chd
    integer s
    integer n_chd
    integer n_tmp
    ! children of patch `p_par` are connected to neighbours on same level if they exist
    !        and temporary boundaries are removed
    !        considers the case that not all four children are present
    ! \update: still used to connect old patches to new patches (new patches already connected now)

    children = dom%patch%elts(p_par+1)%children
    do c = 1, N_CHDRN
       p_chd = children(c)
       if (p_chd .eq. 0) then
          cycle
       end if
       do s = 1, N_BDRY
          n_tmp = dom%patch%elts(p_chd+1)%neigh(s)
          if (n_tmp .ge. 1) cycle ! already connected
          n_chd = find_neigh_patch_Domain(dom, p_par, c - 1, s - 1)
          if (n_chd .eq. 0) then
             cycle
          else
             if (.not.n_chd .eq. n_tmp) then
                dom%bdry_patch%elts(-n_tmp+1)%side = 0
                dom%patch%elts(p_chd+1)%neigh(s) = n_chd
             end if
          end if
       end do
    end do
  end subroutine connect_children
end module refine_patch_mod
