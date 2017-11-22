module wavelet_mod
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  real(8), dimension(9) :: Iu_Base_Wgt
contains
  subroutine forward_wavelet_transform (scaling, wavelet)
    ! Forward wavelet transform
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: scaling, wavelet
    
    integer :: k, l, d

    do k = 1, zlevels
       do l = level_end - 1, level_start - 1, -1
          call update_vector_bdry (scaling(S_MASS:S_TEMP,k), l+1)

          ! Compute scalar wavelet coefficients
          do d = 1, n_domain(rank+1)
             mass => scaling(S_MASS,k)%data(d)%elts
             temp => scaling(S_TEMP,k)%data(d)%elts
             wc_m => wavelet(S_MASS,k)%data(d)%elts
             wc_t => wavelet(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             nullify (wc_m, wc_t, mass, temp)
          end do

          call update_vector_bdry (wavelet(S_MASS:S_TEMP,k), l+1)

          ! Restrict scalars (sub-sample and lift) and velocity (average) to coarser grid
          do d = 1, size(grid)
             mass => scaling(S_MASS,k)%data(d)%elts
             temp => scaling(S_TEMP,k)%data(d)%elts
             wc_m => wavelet(S_MASS,k)%data(d)%elts
             wc_t => wavelet(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (restrict_scalar, grid(d), l, k, 0, 1) ! +1 to include poles
             
             velo => scaling(S_VELO,k)%data(d)%elts
             call apply_interscale_d (restrict_velo, grid(d), l, k, 0, 0)
             nullify (mass, temp, velo, wc_m, wc_t)
          end do
       end do

       scaling(:,k)%bdry_uptodate = .False.
       wavelet(S_MASS:S_TEMP,k)%bdry_uptodate = .False.

       call update_bdry (scaling(S_VELO,k), NONE)

       ! Compute velocity wavelet coefficients
       do l = level_end - 1, level_start - 1, -1
          do d = 1, n_domain(rank+1)
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             velo => scaling(S_VELO,k)%data(d)%elts
             call apply_interscale_d (compute_velo_wavelets, grid(d), l, z_null, 0, 0)
             call apply_to_penta_d (compute_velo_wavelets_penta, grid(d), l, z_null)
             nullify (wc_u, velo)
          end do
       end do
       wavelet(S_VELO,k)%bdry_uptodate = .False.
    end do
  end subroutine forward_wavelet_transform

  subroutine inverse_wavelet_transform (wavelet, scaling, l_start0)
    ! Inverse wavelet transform
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: scaling, wavelet
    integer, optional                                             :: l_start0
    
    integer :: l, d, k, v, l_start

    if (present(l_start0)) then
       l_start = l_start0
    else
       l_start = level_start
    end if

    call update_array_bdry1 (wavelet, level_start, level_end)
    call update_array_bdry1 (scaling, l_start,     level_end)
    
    scaling%bdry_uptodate = .False.

    do k = 1, zlevels
       do l = l_start, level_end-1
          ! Prolong scalars to finer nodes existing at coarser grid (undo lifting)
          do d = 1, n_domain(rank+1)
             mass => scaling(S_MASS,k)%data(d)%elts
             temp => scaling(S_TEMP,k)%data(d)%elts
             wc_m => wavelet(S_MASS,k)%data(d)%elts
             wc_t => wavelet(S_TEMP,k)%data(d)%elts
             if (present(l_start0)) then
                call apply_interscale_d2 (IWT_prolong_scalar,       grid(d), l, z_null, 0, 1) ! needs wc
             else
                call apply_interscale_d2 (IWT_prolong_scalar__fast, grid(d), l, z_null, 0, 1) ! needs wc
             end if
             nullify (mass, temp, wc_m, wc_t)
          end do

          if (l .gt. l_start) call update_bdry__finish (scaling(S_VELO,k), l) ! for next outer velocity
       
          call update_vector_bdry__start (scaling(S_MASS:S_TEMP,k), l+1)

          ! Reconstruct outer velocities at finer edges (interpolate and add wavelet coefficients)
          do d = 1, n_domain(rank+1)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             if (present(l_start0)) then
                call apply_interscale_d2 (IWT_reconstruct_outer_velo, grid(d), l, z_null, 0, 1) ! needs val
             else
                call apply_interscale_d2 (IWT_reconstruct_outer_velo, grid(d), l, z_null, 0, 1) ! needs val
             end if
             call apply_to_penta_d (IWT_reconstruct_velo_penta, grid(d), l, z_null)
             nullify (velo, wc_u)
          end do

          call update_vector_bdry__finish (scaling(S_MASS:S_TEMP,k), l+1)
          call update_bdry__start (scaling(S_VELO,k), l+1)

          ! Reconstruct scalars at finer nodes not existing at coarser grid (interpolate and add wavelet coefficients)
          do d = 1, n_domain(rank+1)
             mass => scaling(S_MASS,k)%data(d)%elts
             temp => scaling(S_TEMP,k)%data(d)%elts
             wc_m => wavelet(S_MASS,k)%data(d)%elts
             wc_t => wavelet(S_TEMP,k)%data(d)%elts
             call apply_interscale_d (IWT_reconstruct_scalar, grid(d), l, z_null, 0, 0)
             nullify (mass, temp, wc_m, wc_t)
          end do
       
          call update_bdry__finish (scaling(S_VELO,k), l+1)

          ! Reconstruct inner velocities at finer edges (interpolate and add wavelet coefficients)
          do d = 1, n_domain(rank+1)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d (IWT_reconstruct_velo_inner, grid(d), l, z_null, 0, 0)
             nullify (velo, wc_u)
          end do
          
          if (l .lt. level_end-1) call update_bdry__start (scaling(S_VELO,k), l+1) ! for next outer velocity
       
          scaling(:,k)%bdry_uptodate = .False.
       end do
    end do
  end subroutine inverse_wavelet_transform

  function velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)
    real(8), dimension(2)          :: velo_interp_penta_corr
    type(Domain)                   :: dom
    integer, dimension(N_BDRY+1)   :: offs, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims, dims_chd
    
    integer :: i, j, i_chd, j_chd

    i = 0
    j = 0
    i_chd = 0
    j_chd = 0

    velo_interp_penta_corr = (/ &
         (Iu_Base_Wgt(9) + dble(dom%I_u_wgt%elts(idx__fast(i_chd+end_pt(1,2,UP+1), &
         j_chd+end_pt(2,2,UP+1), offs_chd(1)) + 1)%enc(9))) * &
         ((-velo(idx(0, -1, offs, dims)*EDGE + UP + 1) - &
         (-velo(idx(-1, -1, offs, dims)*EDGE + 1))) - &
         (velo(ed_idx(i + end_pt(1,1,UP+1), j + end_pt(2,1,UP+1), &
         hex_sides(:,hex_s_offs(UP+1) + 0 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + opp_no(1,2,UP+1), j + opp_no(2,2,UP+1), &
         hex_sides(:,hex_s_offs(UP+1) + 1 + 1), offs, dims) + 1))), &
         (Iu_Base_Wgt(6) + dble(dom%I_u_wgt%elts(idx__fast(i_chd + end_pt(1,2,RT+1), j_chd + &
         end_pt(2,2,RT+1), offs_chd(1)) + 1)%enc(6)))* &
         (velo(idx(-1, -1, offs, dims)*EDGE + 1) + &
         velo(idx(-1, 0, offs, dims)*EDGE + RT + 1) - (velo(ed_idx(i + &
         opp_no(1,1,RT+1), j + opp_no(2,1,RT+1), &
         hex_sides(:,hex_s_offs(RT+1) + 1 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + end_pt(1,1,RT+1), j + end_pt(2,1,RT+1), &
         hex_sides(:,hex_s_offs(RT+1) + 2 + 1), offs, dims) + 1)))/)
  end function velo_interp_penta_corr

  subroutine basic_F_restr_wgt (dom, i_par, j_par, e, offs_par, dims_par, i0, j0, offs, dims, typ)
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, e, i0, j0
    integer, dimension(N_BDRY+1)   :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer, dimension(2)          :: typ
    
    integer                        :: i, j, k

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer, dimension(3,16)       :: ije

    integer, dimension(2)          :: ij_nbp_mp, ij_nbp_pp, ij_nbp_pm, ij_nbp_mm
    integer, dimension(3)          :: ije_lcsd
    integer, dimension(4)          :: id_enc
    real(8), dimension(6)          :: wgt

    if (e .eq. UP) then
       id_enc = (/idx(i0-2,j0, offs, dims), idx(i0-2,j0+1, offs, dims), &
            idx(i0+1,j0, offs, dims), idx(i0+1,j0+1, offs, dims)/)
       i = i0
       j = j0+1
    elseif (e .eq. RT) then
       id_enc = (/idx(i0  ,j0, offs, dims), idx(i0  ,j0+1, offs, dims), &
            idx(i0+1,j0-2, offs, dims), idx(i0+1,j0-1, offs, dims)/)
       i = i0+1
       j = j0
    else
       write(0,*) 'Error 447: R_F_wgts'
       stop
    end if

    ije(:,UMZ+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 1 + 1)
    ije(:,UPZ+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 4 + 1)
    ije(:,WMP+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 0 + 1)
    ije(:,VPP+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 5 + 1)
    ije(:,WPM+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 3 + 1)
    ije(:,VMM+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1) + 2 + 1)

    ij_nbp_mp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 0 + 1)
    ij_nbp_pp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 5 + 1)
    ij_nbp_pm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 3 + 1)
    ij_nbp_mm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 2 + 1)

    ije(:,VMP +1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 4 - 2 + 1)
    ije(:,VMPP+1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1)
    ije(:,UZP +1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 0 + 4 + 1)
    ije(:,WPPP+1) = (/ij_nbp_pp(1), ij_nbp_pp(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 4 - 4 + 1)
    ije(:,WPP +1) = (/ij_nbp_pp(1), ij_nbp_pp(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1)
    ije(:,VPM +1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1)
    ije(:,VPMM+1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 4 - 2 + 1)
    ije(:,UZM +1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 3 - 2 + 1)
    ije(:,WMMM+1) = (/ij_nbp_mm(1), ij_nbp_mm(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1)
    ije(:,WMM +1) = (/ij_nbp_mm(1), ij_nbp_mm(2), 0/) + hex_sides(:,hex_s_offs(e+1) + 4 - 4 + 1)

    k = 1
    if (dist(dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(ije(1,UZP+1),ije(2,UZP+1), &
         adj_tri(:,-k+2,ije(3,UZP+1)+1),offs,dims)+1)) .lt. eps()) then ! COINSIDE

       dom%R_F_wgt%elts(id_enc(1)+1)%enc = 0.0_8
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = 0.0_8
    else

       if (typ(k+1) .eq. OUTER1) then
          ije_lcsd = ije(:,VPP+1)
       else if (typ(k+1) .eq. OUTER2) then
          ije_lcsd = ije(:,WMP+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZP+1)
       end if

       wgt = interp_F_wgts(e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, &
            j_par, adj_tri(:,k+1,e+1), offs_par, dims_par) + 1), ije, &
            (/VPP, WMP, UZP, UPZ, WPP, VMP, UMZ, WPPP, VMPP/))

       if (e .eq. RT) then
          wgt = (/wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)/)
       else if (e .eq. UP) then
          wgt = (/wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)/)
       end if

       dom%R_F_wgt%elts(id_enc(1)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = wgt(4:6)
    end if

    k = 0
    if (dist(dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         
         dom%ccentre%elts(tri_idx(ije(1,UZM+1),ije(2,UZM+1), &
         
         adj_tri(:,-k+2,ije(3,UZM+1)+1),offs,dims)+1)) .lt. eps()) then ! COINSIDE

       dom%R_F_wgt%elts(id_enc(3)+1)%enc = 0.0_8
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = 0.0_8
    else

       if (typ(k+1) .eq. OUTER1) then
          ije_lcsd = ije(:,VMM+1)
       else if (typ(k+1) .eq. OUTER2) then
          ije_lcsd = ije(:,WPM+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZM+1)
       end if

       wgt = interp_F_wgts(e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, &
            j_par, adj_tri(:,k+1,e+1), offs_par, dims_par) + 1), ije, &
            (/VMM, WPM, UZM, UMZ, WMM, VPM, UPZ, WMMM, VPMM/))

       if (e .eq. UP) then
          wgt = (/wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)/)
       else if (e .eq. RT) then
          wgt = (/wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)/)
       end if

       dom%R_F_wgt%elts(id_enc(3)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = wgt(4:6)
    end if

  contains

    function interp_F_wgts (e, k1, ije_lcsd, endpt_o, ije, stencil)
      real(8), dimension(6)    :: interp_F_wgts
      integer                  :: e, k1
      integer, dimension(3,16) :: ije
      integer, dimension(9)    :: stencil
      type(Coord)              :: endpt_o
      
      integer                 :: id_tri, info
      integer, dimension(3)   :: ije_lcsd
      integer, dimension(6)   :: ipiv
      real(8), dimension(6,6) :: G
      real(8), dimension(6)   :: b
      type(Coord)             :: endpt, x, y
      
      id_tri = tri_idx(ije_lcsd(1), ije_lcsd(2), adj_tri(:,-k1+2,ije_lcsd(3) + 1), offs, dims)

      call local_coord(dom%ccentre%elts(id_tri+1), &
           dom%ccentre%elts(id_tri+1), dom%midpt%elts(ed_idx(0, 0, &
           ije_lcsd, offs, dims) + 1), x, y)

      G(:,1) = coords_to_row(ije(:,stencil(1) + 1), x, y)
      G(:,2) = coords_to_row(ije(:,stencil(2) + 1), x, y)
      G(:,3) = coords_to_row(ije(:,stencil(3) + 1), x, y)
      G(:,4) = coords_to_row(ije(:,stencil(4) + 1), x, y) - coords_to_row(ije(:,stencil(5) + 1), x, y)
      G(:,5) = coords_to_row(ije(:,stencil(6) + 1), x, y) - coords_to_row(ije(:,stencil(7) + 1), x, y)
      G(:,6) = coords_to_row(ije(:,stencil(8) + 1), x, y) - coords_to_row(ije(:,stencil(9) + 1), x, y)

      endpt = endpt_o
      b = coords_to_row_perp((/dom%ccentre%elts(id_tri+1), endpt/), x, y)
      ipiv = 0
      info = 0
      call dgesv(6, 1, G, 6, ipiv, b, 6, info)
      interp_F_wgts = b
      return
    end function interp_F_wgts

    function coords_to_row_perp (coords, x, y)
      real(8), dimension(6)     :: coords_to_row_perp
      type(Coord), dimension(2) :: coords
      type(Coord)               :: x, y
      
      type(Coord) :: dirvec, midpt

      midpt = mid_pt(coords(1), coords(2))

      dirvec = cross(vector(coords(1), coords(2)), midpt)

      coords_to_row_perp = coords_to_rowd(midpt, dirvec, x, y)*dist(coords(1), coords(2))
    end function coords_to_row_perp

    function coords_to_row (ije0, x, y)
      real(8), dimension(6) :: coords_to_row
      integer, dimension(3) :: ije0
      type(Coord)           :: x, y
      
      integer     :: i, j, e
      type(Coord) :: endpt1, endpt2, midpt
      real(8)     :: pedlen

      i = ije0(1)
      j = ije0(2)
      e = ije0(3)

      endpt1 = dom%node%elts(idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), offs, dims) + 1)
      endpt2 = dom%node%elts(idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), offs, dims) + 1)

      pedlen = dist(dom%ccentre%elts(tri_idx(i, j, adj_tri(:,1,e+1), offs, &
           dims) + 1), dom%ccentre%elts(tri_idx(i, j, adj_tri(:,2,e+1), &
           offs, dims) + 1))

      midpt = mid_pt(endpt1, endpt2)

      coords_to_row = coords_to_rowd(midpt, vector(endpt1, endpt2), x, &
           y)*pedlen
    end function coords_to_row
  end subroutine basic_F_restr_wgt

  subroutine compute_velo_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute velocity wavelet coefficients (except at pentagons)
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: e, id1, id2, id_chd, idN_chd, idE_chd, idNE_chd
    real(8)               :: u
    real(8), dimension(6) :: u_inner

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)

    do e = 1, EDGE
       id1 = idx(i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2 = idx(i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       if (dom%mask_e%elts(EDGE*id2+e) .lt. ADJZONE) cycle

       u = interp_outer_velo (dom, i_par, j_par, e - 1, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

       wc_u(EDGE*id2+e) = velo(EDGE*id2+e) - u
       wc_u(EDGE*id1+e) = u - velo(EDGE*id2+e)
    end do

    u_inner = Interp_velo_inner(dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, &
         offs_chd, dims_chd, EDGE*idE_chd + UP, DG + EDGE*idE_chd, &
         EDGE*idNE_chd + RT, EDGE*idN_chd + RT, DG + EDGE*idN_chd, &
         EDGE*idNE_chd + UP)

    if (dom%mask_e%elts(EDGE*idE_chd+UP+1) .ge. ADJZONE) &
         wc_u(EDGE*idE_chd+UP+1) = velo(EDGE*idE_chd+UP+1) - u_inner(1)
    if (dom%mask_e%elts(DG+EDGE*idE_chd+1) .ge. ADJZONE) &
         wc_u(DG+EDGE*idE_chd+1) = velo(DG+EDGE*idE_chd+1) - u_inner(2)
    if (dom%mask_e%elts(EDGE*idNE_chd+RT+1) .ge. ADJZONE) &
         wc_u(EDGE*idNE_chd+RT+1) = velo(EDGE*idNE_chd+RT+1) - u_inner(3)
    if (dom%mask_e%elts(EDGE*idN_chd+RT+1) .ge. ADJZONE) &
         wc_u(EDGE*idN_chd+RT+1) = velo(EDGE*idN_chd+RT+1) - u_inner(4)
    if (dom%mask_e%elts(DG+EDGE*idN_chd+1) .ge. ADJZONE) &
         wc_u(DG+EDGE*idN_chd+1) = velo(DG+EDGE*idN_chd+1) - u_inner(5)
    if (dom%mask_e%elts(EDGE*idNE_chd+UP+1) .ge. ADJZONE) &
         wc_u(EDGE*idNE_chd+UP+1) = velo(EDGE*idNE_chd+UP+1) - u_inner(6)
  end subroutine compute_velo_wavelets

  subroutine compute_velo_wavelets_penta (dom, p, c, offs, dims, zlev)
    ! Compute velocity wavelet coefficients at pentagons
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer                        :: id_chd, idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(8)                        :: v
    real(8), dimension(2)          :: u

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd .eq. 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c .eq. IMINUSJPLUS) then
       id_chd  = idx(0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx(0, LAST,   offs_chd, dims_chd)

       v = -(Iu_Base_Wgt(8) + dble(dom%I_u_wgt%elts(idN_chd+1)%enc(8)))*( &
            velo(idx(0, PATCH_SIZE, offs, dims)*EDGE + UP + 1) &
            +velo(idx(-1, PATCH_SIZE, offs, dims)*EDGE + RT + 1))

       if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) then
          wc_u(EDGE*id_chd+UP+1)  = wc_u(EDGE*id_chd +UP+1) - v
          wc_u(EDGE*idN_chd+UP+1) = wc_u(EDGE*idN_chd+UP+1) + v
       end if
    else
       if (c .eq. IPLUSJMINUS) then
          id_chd = idx(LAST - 1, 0, offs_chd, dims_chd)
          idE_chd = idx(LAST, 0, offs_chd, dims_chd)

          v = (Iu_Base_Wgt(7) + dble(dom%I_u_wgt%elts(idE_chd+1)%enc(7)))*( &
               velo(idx(PATCH_SIZE, 0, offs, dims)*EDGE + RT + 1) &
               +velo(idx(PATCH_SIZE,-1, offs, dims)*EDGE + UP + 1))

          if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) then
             wc_u(EDGE*id_chd+RT+1)  = wc_u(EDGE*id_chd +RT+1) - v
             wc_u(EDGE*idE_chd+RT+1) = wc_u(EDGE*idE_chd+RT+1) + v
          end if
       end if
    end if

    if (.not. c .eq. IJMINUS) return

    id_chd  = idx(0, 0, offs_chd, dims_chd)
    idN_chd = idx(0, 1, offs_chd, dims_chd)
    idE_chd = idx(1, 0, offs_chd, dims_chd)

    u = velo_interp_penta_corr(dom, offs, dims, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) then
       wc_u(EDGE*id_chd+UP+1)  = wc_u(EDGE*id_chd +UP+1) + u(1)
       wc_u(EDGE*idN_chd+UP+1) = wc_u(EDGE*idN_chd+UP+1) - u(1)
    end if
    if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) then
       wc_u(EDGE*id_chd+RT+1)  = wc_u(EDGE*id_chd +RT+1) + u(2)
       wc_u(EDGE*idE_chd+RT+1) = wc_u(EDGE*idE_chd+RT+1) - u(2)
    end if
  end subroutine compute_velo_wavelets_penta

  subroutine init_wavelets()
    integer :: d, i, k, num, v

    do k = 1, zlevels
       do v = S_MASS, S_VELO
          call init_Float_Field(wav_coeff(v,k), POSIT(v))
          call init_Float_Field(trend_wav_coeff(v,k), POSIT(v))
       end do
    end do

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init(grid(d)%overl_areas, num)
       call init(grid(d)%I_u_wgt, num)

       do i = 1, num
          call init_Iu_Wgt(grid(d)%I_u_wgt%elts(i), Iu_Base_Wgt)
       end do

       call init(grid(d)%R_F_wgt, num)

       do i = 1, num
          call init_RF_Wgt(grid(d)%R_F_wgt%elts(i), (/0.0_4, 0.0_4, 0.0_4/))
       end do

       do k = 1, zlevels
          do v = S_MASS, S_TEMP
             call init(wav_coeff(v,k)%data(d), num)
             call init(trend_wav_coeff(v,k)%data(d), num)
          end do
          call init(wav_coeff(S_VELO,k)%data(d), EDGE*num)
          call init(trend_wav_coeff(S_VELO,k)%data(d), EDGE*num)
       end do
    end do
  end subroutine init_wavelets

  subroutine get_overl_areas (dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par, offs_chd, dims_chd, e, area, typ)
    type(Domain)                   :: dom
    integer                        :: e, i_par, j_par, i_chd, j_chd
    integer, dimension(2)          :: typ
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    real(8), dimension(8)          :: area
    
    type(Coord), dimension(6)   :: hex
    type(Coord), dimension(3,2) :: tri
    type(Coord)                 :: inters_pt0, inters_pt1, pt
    integer                     :: i
    logical                     :: does_inters0, does_inters1, troubles

    area = 0.0_8
    typ = 0

    hex = (/ (dom%ccentre%elts(tri_idx(i_chd, j_chd, no_adj_tri(:,i + &
         hex_s_offs(-e+3) + 1), offs_chd, dims_chd) + 1), i = 0, 6-1) /)

    tri = reshape((/dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,4,e+1), &
         offs_par, dims_par) + 1), dom%ccentre%elts(tri_idx(i_par, j_par, &
         adj_tri(:,2,e+1), offs_par, dims_par) + 1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,3,e+1), &
         offs_par, dims_par) + 1), dom%ccentre%elts(tri_idx(i_par, j_par, &
         bfly_tri(:,2,e+1), offs_par, dims_par) + 1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,1,e+1), &
         offs_par, dims_par) + 1), dom%ccentre%elts(tri_idx(i_par, j_par, &
         bfly_tri(:,1,e+1), offs_par, dims_par) + 1)/), (/3, 2/))

    pt = dom%node%elts(idx(i_chd, j_chd, offs_chd, dims_chd) + 1)

    area(1) = triarea(hex(6), hex(1), pt)
    area(2) = triarea(hex(3), hex(4), pt)

    do i = 1, 2
       call arc_inters(tri(1,i), tri(2,i), hex(3*i-2), hex(3*i-1), &
            inters_pt0, does_inters0, troubles)
       call arc_inters(tri(3,i), tri(2,i), hex(3*i), hex(3*i-1), inters_pt1, &
            does_inters1, troubles)
       if (does_inters0 .and. does_inters1) then
          area(i+4) = triarea(inters_pt0, tri(2,i), hex(3*i-1))
          area(i+6) = triarea(tri(2,i), hex(3*i-1), inters_pt1)
          area(i+2) = area(i+4) + area(i+6)
          area(i) = area(i) + triarea(hex(3*i-2), inters_pt0, pt) + &
               triarea(inters_pt0, pt, tri(2,i))
          area(-i+3) = area(-i+3) + triarea(inters_pt1, hex(3*i), pt) + &
               triarea(tri(2,i), pt, inters_pt1)
          typ(-i+3) = INSIDE
       else
          if (.not. does_inters0 .and. .not. does_inters1) then
             area(i+2) = 0.0_8
             call arc_inters(tri(2,1), tri(2,2), hex(3*i-2), hex(3*i-1), &
                  inters_pt0, does_inters0, troubles)
             call arc_inters(tri(2,2), tri(2,1), hex(3*i-1), hex(3*i), &
                  inters_pt1, does_inters1, troubles)
             if (.not. does_inters0 .and. does_inters1) then
                area(i) = area(i) + triarea(hex(3*i-2), hex(3*i-1), pt) + &
                     triarea(hex(3*i-1), inters_pt1, pt)
                area(-i+3) = area(-i+3) + triarea(inters_pt1, hex(3*i), &
                     pt)
                typ(-i+3) = OUTER2
             else
                if (does_inters0 .and. .not. does_inters1) then
                   area(i) = area(i) + triarea(hex(3*i-2), inters_pt0, &
                        pt)
                   area(-i+3) = area(-i+3) + triarea(hex(3*i-1), &
                        hex(3*i), pt) + triarea(inters_pt0, &
                        hex(3*i-1), pt)
                   typ(-i+3) = OUTER1
                else
                   write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'A', does_inters0, does_inters1
                end if
             end if
          else
             write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'B', does_inters0, does_inters1
          end if
       end if
    end do
    return
  end subroutine get_overl_areas

  subroutine normalize2 (q, u, v)
    real(8), dimension(2) :: q
    real(8)               :: nrm, u, v

    nrm = sqrt(q(1)**2 + q(2)**2)
    u = q(1)/nrm
    v = q(2)/nrm
  end subroutine normalize2

  subroutine IWT_prolong_scalar__fast (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars at fine points existing at coarse scale by undoing lifting
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd
   
    ! Locally filled, IWT reproduces previous value

    id_par = idx(i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) .ge. TOLRNZ) return

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)

    mass(id_chd+1) = prolong(mass(id_par+1), wc_m, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
    temp(id_chd+1) = prolong(temp(id_par+1), wc_t, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
  end subroutine IWT_prolong_scalar__fast

  subroutine IWT_prolong_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars at fine points existing at coarse scale by undoing lifting
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_par = idx(i_par, j_par, offs_par, dims_par)
    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    
    if (dom%mask_n%elts(id_chd+1) .eq. FROZEN) return ! FROZEN mask -> do not overide with wrong value

    mass(id_chd+1) = prolong(mass(id_par+1), wc_m, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
    temp(id_chd+1) = prolong(temp(id_par+1), wc_t, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
  end subroutine IWT_prolong_scalar

  function prolong (scalar, wavelet, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
    ! Prolongation at fine points existing at coarse scale by undoing lifting
    real(8)                        :: prolong
    real(8)                        :: scalar
    real(8), dimension(:), pointer :: wavelet
    type(Domain)                   :: dom
    integer                        :: id_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd

    integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE
    
    idE   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN2E = idx(i_chd + 2, j_chd + 1, offs_chd, dims_chd)
    id2NE = idx(i_chd + 1, j_chd + 2, offs_chd, dims_chd)
    idN   = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idW   = idx(i_chd - 1, j_chd,     offs_chd, dims_chd)
    idNW  = idx(i_chd - 1, j_chd + 1, offs_chd, dims_chd)
    idS2W = idx(i_chd - 2, j_chd - 1, offs_chd, dims_chd)
    idSW  = idx(i_chd - 1, j_chd - 1, offs_chd, dims_chd)
    idS   = idx(i_chd,     j_chd - 1, offs_chd, dims_chd)
    id2SW = idx(i_chd - 1, j_chd - 2, offs_chd, dims_chd)
    idSE  = idx(i_chd + 1, j_chd - 1, offs_chd, dims_chd)

    if (id2NE+1 .le.  ubound(wavelet,dim=1)) then
       prolong = scalar - &
            (wavelet(idE+1)*dom%overl_areas%elts(idE+1)%a(1) + &
            wavelet(idNE+1)*dom%overl_areas%elts(idNE+1)%a(2) + &
            wavelet(idN2E+1)*dom%overl_areas%elts(idN2E+1)%a(3) + &
            wavelet(id2NE+1)*dom%overl_areas%elts(id2NE+1)%a(4) + &
            wavelet(idN+1)*dom%overl_areas%elts(idN+1)%a(1) + &
            wavelet(idW+1)*dom%overl_areas%elts(idW+1)%a(2) + &
            wavelet(idNW+1)*dom%overl_areas%elts(idNW+1)%a(3) + &
            wavelet(idS2W+1)*dom%overl_areas%elts(idS2W+1)%a(4) + &
            wavelet(idSW+1)*dom%overl_areas%elts(idSW+1)%a(1) + &
            wavelet(idS+1)*dom%overl_areas%elts(idS+1)%a(2) + &
            wavelet(id2SW+1)*dom%overl_areas%elts(id2SW+1)%a(3) + &
            wavelet(idSE+1)*dom%overl_areas%elts(idSE+1)%a(4))*dom%areas%elts(id_par+1)%hex_inv
    else
       prolong = scalar - &
            (wavelet(idE+1)*dom%overl_areas%elts(idE+1)%a(1) + &
            wavelet(idNE+1)*dom%overl_areas%elts(idNE+1)%a(2) + &
            wavelet(idN2E+1)*dom%overl_areas%elts(idN2E+1)%a(3) + &
            wavelet(idN+1)*dom%overl_areas%elts(idN+1)%a(1) + &
            wavelet(idW+1)*dom%overl_areas%elts(idW+1)%a(2) + &
            wavelet(idNW+1)*dom%overl_areas%elts(idNW+1)%a(3) + &
            wavelet(idS2W+1)*dom%overl_areas%elts(idS2W+1)%a(4) + &
            wavelet(idSW+1)*dom%overl_areas%elts(idSW+1)%a(1) + &
            wavelet(idS+1)*dom%overl_areas%elts(idS+1)%a(2) + &
            wavelet(id2SW+1)*dom%overl_areas%elts(id2SW+1)%a(3) + &
            wavelet(idSE+1)*dom%overl_areas%elts(idSE+1)%a(4))*dom%areas%elts(id_par+1)%hex_inv
    end if
  end function prolong

  subroutine compute_scalar_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute wavelet coefficients for mass and potential temperature
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd    = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd   = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    id2N_chd  = idx(i_chd,     j_chd + 2, offs_chd, dims_chd)
    id2E_chd  = idx(i_chd + 2, j_chd,     offs_chd, dims_chd)
    id2S_chd  = idx(i_chd,     j_chd - 2, offs_chd, dims_chd)
    id2W_chd  = idx(i_chd - 2, j_chd,     offs_chd, dims_chd)
    id2NE_chd = idx(i_chd + 2, j_chd + 2, offs_chd, dims_chd)

    if (dom%mask_n%elts(idNE_chd+1) .ge. ADJZONE) then
       wc_m(idNE_chd+1) = mass(idNE_chd+1) - Interp_node(dom, mass, idNE_chd, id2NE_chd, id_chd, id2E_chd, id2N_chd)
       wc_t(idNE_chd+1) = temp(idNE_chd+1) - Interp_node(dom, temp, idNE_chd, id2NE_chd, id_chd, id2E_chd, id2N_chd)
    end if

    if (dom%mask_n%elts(idN_chd+1) .ge. ADJZONE) then
       wc_m(idN_chd+1) = mass(idN_chd+1) - Interp_node(dom, mass, idN_chd, id_chd, id2N_chd, id2W_chd, id2NE_chd)
       wc_t(idN_chd+1) = temp(idN_chd+1) - Interp_node(dom, temp, idN_chd, id_chd, id2N_chd, id2W_chd, id2NE_chd)
    end if

    if (dom%mask_n%elts(idE_chd+1) .ge. ADJZONE) then
       wc_m(idE_chd+1) = mass(idE_chd+1) - Interp_node(dom, mass, idE_chd, id_chd, id2E_chd, id2NE_chd, id2S_chd)
       wc_t(idE_chd+1) = temp(idE_chd+1) - Interp_node(dom, temp, idE_chd, id_chd, id2E_chd, id2NE_chd, id2S_chd)
    end if
  end subroutine compute_scalar_wavelets

  subroutine IWT_reconstruct_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct scalars at fine nodes not existing at coarse scale
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    if (dom%mask_n%elts(id_chd+1) .eq. FROZEN) return ! FROZEN mask -> do not overide with wrong value

    idN_chd   = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    id2N_chd  = idx(i_chd,     j_chd + 2, offs_chd, dims_chd)
    id2E_chd  = idx(i_chd + 2, j_chd,     offs_chd, dims_chd)
    id2S_chd  = idx(i_chd,     j_chd - 2, offs_chd, dims_chd)
    id2W_chd  = idx(i_chd - 2, j_chd,     offs_chd, dims_chd)
    id2NE_chd = idx(i_chd + 2, j_chd + 2, offs_chd, dims_chd)

    ! Interpolate scalars and add wavelets to reconstruct values at fine scale
    mass(idNE_chd+1) = Interp_node (dom, mass, idNE_chd, id2NE_chd, id_chd, id2E_chd, id2N_chd) + wc_m(idNE_chd+1)
    mass(idN_chd+1)  = Interp_node (dom, mass, idN_chd, id_chd, id2N_chd, id2W_chd, id2NE_chd)  + wc_m(idN_chd+1)
    mass(idE_chd+1)  = Interp_node (dom, mass, idE_chd, id_chd, id2E_chd, id2NE_chd, id2S_chd)  + wc_m(idE_chd+1)

    temp(idNE_chd+1) = Interp_node (dom, temp, idNE_chd, id2NE_chd, id_chd, id2E_chd, id2N_chd) + wc_t(idNE_chd+1)
    temp(idN_chd+1)  = Interp_node (dom, temp, idN_chd, id_chd, id2N_chd, id2W_chd, id2NE_chd)  + wc_t(idN_chd+1)
    temp(idE_chd+1)  = Interp_node (dom, temp, idE_chd, id_chd, id2E_chd, id2NE_chd, id2S_chd)  + wc_t(idE_chd+1)
  end subroutine IWT_reconstruct_scalar

  function Interp_node (dom, var, id, id1, id2, id3, id4)
    ! Interpolation at nodes
    real(8)                        :: Interp_node
    type(Domain)                   :: dom
    real(8), dimension(:), pointer :: var
    integer                        :: id, id1, id2, id3, id4

    Interp_node = (dom%overl_areas%elts(id+1)%a(1)*var(id1+1) + &
         dom%overl_areas%elts(id+1)%a(2)*var(id2+1) + &
         dom%overl_areas%elts(id+1)%a(3)*var(id3+1) + &
         dom%overl_areas%elts(id+1)%a(4)*var(id4+1))*dom%areas%elts(id+1)%hex_inv
  end function Interp_node

  subroutine local_coord (midpt, endpt1, endpt2, x, y)
    type(Coord) :: midpt, endpt1, endpt2, x, y

    type(Coord) :: y0

    x = direction(endpt1, endpt2)
    y0 = cross(x, midpt)
    y = normalize_Coord(y0)
    return
  end subroutine local_coord

  subroutine IWT_reconstruct_outer_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at outer fine edges not existing on coarse grid
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: e, id_chd, id_par, id1, id2
    real(8) :: u

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)

    do e = 1, EDGE
       id1 = idx(i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2 = idx(i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       id_par = idx(i_par, j_par, offs_par, dims_par)

       u = Interp_outer_velo (dom, i_par, j_par, e - 1, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

       velo(EDGE*id2+e) = u + wc_u(EDGE*id2+e)
       velo(EDGE*id1+e) = 2.0_8*velo(EDGE*id_par+e) - u + wc_u(EDGE*id1+e)
    end do
  end subroutine IWT_reconstruct_outer_velo

  function Interp_outer_velo (dom, i, j, e, offs, dims, i_chd, j_chd, offs_chd, dims_chd)
    ! Interpolate outer velocities to fine edges
    real(8)                        :: interp_outer_velo
    type(Domain)                   :: dom
    integer                        :: i, j, e, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims, dims_chd

    integer :: ide
    real(8) :: wgt(9)
    
    ide = idx(i_chd+end_pt(1,2,e+1),j_chd+end_pt(2,2,e+1),offs_chd,dims_chd)
    
    wgt = Iu_Base_Wgt + dble(dom%I_u_wgt%elts(ide+1)%enc)
    
    interp_outer_velo = sum (wgt* &
         (/velo(idx(i, j, offs, dims)*EDGE + e + 1), &
         velo(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims) + 1), &
         velo(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims) + 1), &
         velo(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims) + 1), &
         velo(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims) + 1), &
         velo(ed_idx(i + opp_no(1,1,e+1), j + opp_no(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims) + 1), &
         velo(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + opp_no(1,1,e+1), j + opp_no(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims) + 1), &
         velo(ed_idx(i + opp_no(1,2,e+1), j + opp_no(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims) + 1), &
         velo(ed_idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims) + 1) - &
         velo(ed_idx(i + opp_no(1,2,e+1), j + opp_no(2,2,e+1), &
         hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims) + 1)/))
  end function Interp_outer_velo

  subroutine IWT_reconstruct_velo_inner (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at inner fine edges not existing on coarse grid
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_chd, idN_chd, idE_chd, idNE_chd
    real(8), dimension(6) :: u_inner

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)

    u_inner = Interp_velo_inner (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, &
         offs_chd, dims_chd, EDGE*idE_chd + UP, DG + EDGE*idE_chd, &
         EDGE*idNE_chd + RT, EDGE*idN_chd + RT, DG + EDGE*idN_chd, &
         EDGE*idNE_chd + UP)

    velo(EDGE*idE_chd +UP+1) = u_inner(1) + wc_u(EDGE*idE_chd +UP+1)
    velo(EDGE*idE_chd +DG+1) = u_inner(2) + wc_u(EDGE*idE_chd +DG+1)
    velo(EDGE*idNE_chd+RT+1) = u_inner(3) + wc_u(EDGE*idNE_chd+RT+1)
    velo(EDGE*idN_chd +RT+1) = u_inner(4) + wc_u(EDGE*idN_chd +RT+1)
    velo(EDGE*idN_chd +DG+1) = u_inner(5) + wc_u(EDGE*idN_chd +DG+1)
    velo(EDGE*idNE_chd+UP+1) = u_inner(6) + wc_u(EDGE*idNE_chd+UP+1)
  end subroutine IWT_reconstruct_velo_inner

  function Interp_velo_inner (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, &
       idE_UP, idE_DG, idNE_RT, idN_RT, idN_DG, idNE_UP)
    ! Interpolate inner velocities to fine edges
    real(8), dimension(6)          :: Interp_velo_inner
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, idE_UP, idE_DG, idNE_RT, idN_RT, idN_DG,  idNE_UP
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    
    integer               :: id, id_par, id1_par, id2_par, k, idN, idUP, idDG, idRT, idE, idNE, idN2E, id2NE, idN2, idE2
    real(8), dimension(2) :: curl_u
    real(8), dimension(6) :: u
   
    id_par = idx(i_par, j_par, offs_par, dims_par)

    curl_u = 0.0_8
    do k = 1, TRIAG
       id1_par = idx(i_par - k + 2, j_par,         offs_par, dims_par)
       id2_par = idx(i_par,         j_par + k - 1, offs_par, dims_par)
       curl_u(k) = (velo(DG+EDGE*id_par+1)*dom%len%elts(DG+EDGE*id_par+1) &
                  + velo(EDGE*id1_par+UP+1)*dom%len%elts(EDGE*id1_par+UP+1) &
                  + velo(EDGE*id2_par+RT+1)*dom%len%elts(EDGE*id2_par+RT+1)) &
                  / dom%triarea%elts(TRIAG*id_par+k)
    end do

    id    = idx(i_chd, j_chd, offs_chd, dims_chd)
    idN   = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
    idUP  = EDGE*id + UP
    idDG  = EDGE*id + DG
    idRT  = EDGE*id + RT
    idE   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN2E = idx(i_chd + 2, j_chd + 1, offs_chd, dims_chd)
    id2NE = idx(i_chd + 1, j_chd + 2, offs_chd, dims_chd)
    idN2  = idx(i_chd,     j_chd + 2, offs_chd, dims_chd)
    idE2  = idx(i_chd + 2, j_chd,     offs_chd, dims_chd)

    u = 0.0_8
    
    u(1) = (dom%triarea%elts(LORT+TRIAG*id+1)*curl_u(1) - &
         velo(idRT+1)*dom%len%elts(idRT+1) - &
         velo(idDG+1)*dom%len%elts(idDG+1))/dom%len%elts(idE_UP+1)

    u(2) = (dom%triarea%elts(LORT+TRIAG*idE+1)*curl_u(1) - &
         velo(EDGE*idE+RT+1)*dom%len%elts(EDGE*idE+RT+1) - &
         velo(EDGE*idE2+UP+1)*dom%len%elts(EDGE*idE2+UP+1))/dom%len%elts(idE_DG+1)

    u(3) = (dom%triarea%elts(LORT+TRIAG*idNE+1)*curl_u(1) - &
         velo(DG+EDGE*idNE+1)*dom%len%elts(DG+EDGE*idNE+1) - &
         velo(EDGE*idN2E+UP+1)*dom%len%elts(EDGE*idN2E+UP+1))/dom%len%elts(idNE_RT+1)

    u(4) = (dom%triarea%elts(TRIAG*id+UPLT+1)*curl_u(2) - &
         velo(idUP+1)*dom%len%elts(idUP+1) - &
         velo(idDG+1)*dom%len%elts(idDG+1))/dom%len%elts(idN_RT+1)

    u(5) = (dom%triarea%elts(TRIAG*idN+UPLT+1)*curl_u(2) - &
         velo(EDGE*idN+UP+1)*dom%len%elts(EDGE*idN+UP+1) - &
         velo(EDGE*idN2+RT+1)*dom%len%elts(EDGE*idN2+RT+1))/dom%len%elts(idN_DG+1)

    u(6) = (dom%triarea%elts(TRIAG*idNE+UPLT+1)*curl_u(2) - &
         velo(DG+EDGE*idNE+1)*dom%len%elts(DG+EDGE*idNE+1) - &
         velo(EDGE*id2NE+RT+1)*dom%len%elts(EDGE*id2NE+RT+1))/dom%len%elts(idNE_UP+1)

    Interp_velo_inner = u
  end function Interp_velo_inner

  subroutine IWT_reconstruct_velo_penta (dom, p, c, offs, dims, z_null)
    ! Interpolate velocity to fine edges at pentagons
    type(Domain)                   :: dom
    integer                        :: p, c, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer                        :: id_chd,  idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(8)                        :: v
    real(8), dimension(2)          :: u

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd .eq. 0) return

    call get_offs_Domain(dom, p_chd, offs_chd, dims_chd)

    if (c .eq. IMINUSJPLUS) then
       id_chd  = idx(0, LAST - 1, offs_chd, dims_chd)
       idN_chd = idx(0, LAST,     offs_chd, dims_chd)

       v = (Iu_Base_Wgt(8) + dble(dom%I_u_wgt%elts(idN_chd+1)%enc(8)))*( &
            velo(EDGE*idx( 0, PATCH_SIZE, offs, dims) + UP + 1) &
            +velo(EDGE*idx(-1, PATCH_SIZE, offs, dims) + RT + 1))

       velo(EDGE*id_chd+UP+1)  = velo(EDGE*id_chd +UP+1) - v
       velo(EDGE*idN_chd+UP+1) = velo(EDGE*idN_chd+UP+1) + v
    else
       if (c .eq. IPLUSJMINUS) then
          id_chd  = idx(LAST - 1, 0, offs_chd, dims_chd)
          idE_chd = idx(LAST,     0, offs_chd, dims_chd)

          v = -(Iu_Base_Wgt(7) + dble(dom%I_u_wgt%elts(idE_chd+1)%enc(7)))*( &
               velo(EDGE*idx(PATCH_SIZE,  0, offs, dims) + RT + 1) &
               +velo(EDGE*idx(PATCH_SIZE, -1, offs, dims) + UP + 1))

          velo(EDGE*id_chd +RT+1) = velo(EDGE*id_chd +RT+1) - v
          velo(EDGE*idE_chd+RT+1) = velo(EDGE*idE_chd+RT+1) + v
       end if
    end if

    if (.not. c .eq. IJMINUS) return

    id_chd  = idx(0, 0, offs_chd, dims_chd)
    idN_chd = idx(0, 1, offs_chd, dims_chd)
    idE_chd = idx(1, 0, offs_chd, dims_chd)

    u = velo_interp_penta_corr(dom, offs, dims, offs_chd, dims_chd)

    velo(EDGE*id_chd +UP+1) = velo(EDGE*id_chd +UP+1) - u(1)
    velo(EDGE*idN_chd+UP+1) = velo(EDGE*idN_chd+UP+1) + u(1)
    velo(EDGE*id_chd +RT+1) = velo(EDGE*id_chd +RT+1) - u(2)
    velo(EDGE*idE_chd+RT+1) = velo(EDGE*idE_chd+RT+1) + u(2)
  end subroutine IWT_reconstruct_velo_penta

  subroutine set_RF_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd,  p_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    
    integer               :: id_chd, idN_chd, idE_chd, idNE_chd
    integer, dimension(2) :: typ
    real(8), dimension(8) :: area
    
    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)

    call get_overl_areas(dom, i_par, j_par, i_chd + 1, j_chd, offs_par, dims_par, offs_chd, dims_chd, RT, area, typ)
    call init_Overl_Area(dom%overl_areas%elts(idE_chd+1), area)
    call basic_F_restr_wgt(dom, i_par, j_par, RT, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)

    call get_overl_areas(dom, i_par, j_par, i_chd + 1, j_chd + 1, offs_par, dims_par, offs_chd, dims_chd, DG, area, typ)
    call init_Overl_Area(dom%overl_areas%elts(idNE_chd+1), area)

    call get_overl_areas(dom, i_par, j_par, i_chd, j_chd + 1, offs_par, dims_par, offs_chd, dims_chd, UP, area, typ)
    call init_Overl_Area(dom%overl_areas%elts(idN_chd+1), area)
    call basic_F_restr_wgt(dom, i_par, j_par, UP, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)

    call set_coarse_overlay()
  end subroutine set_RF_wgts

  subroutine set_coarse_overlay ()
    ! Set overlay quantities on coarsest level
    integer :: d, p

    p = 2
    do d = 1, size(grid)
       call apply_onescale_to_patch(zero_overlay, grid(d), p - 1, z_null, 0, 1)
    end do
  end subroutine set_coarse_overlay

 subroutine zero_overlay (dom, i, j, zlev, offs, dims)
    ! Sets overlay values to zero
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id

    id = idx(i, j, offs, dims)
    dom%overl_areas%elts(id+1)%a     = 0.0_8
    dom%overl_areas%elts(id+1)%split = 0.0_8
  end subroutine zero_overlay

  subroutine set_WT_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd

    id_chd  = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)

    dom%I_u_wgt%elts(idE_chd+1) = outer_velo_weights(dom, p_chd, i_par, j_par, RT, offs_par, dims_par)
    dom%I_u_wgt%elts(id_chd+1)  = outer_velo_weights(dom, p_chd, i_par, j_par, DG, offs_par, dims_par)
    dom%I_u_wgt%elts(idN_chd+1) = outer_velo_weights(dom, p_chd, i_par, j_par, UP, offs_par, dims_par)
  end subroutine set_WT_wgts

  subroutine restrict_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Restrict velocity by averaging
    type(Domain)                   :: dom
    integer                        :: i_par, j_par,  i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id_par

    id_par   = idx(i_par, j_par, offs_par, dims_par)

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) .gt. 0) &
         velo(EDGE*id_par+RT+1) = 0.5_8*(velo(EDGE*id_chd+RT+1) + velo(EDGE*idE_chd+RT+1))

    if (dom%mask_e%elts(EDGE*id_chd+DG+1) .gt. 0) &
         velo(EDGE*id_par+DG+1) = 0.5_8*(velo(EDGE*idNE_chd+DG+1) + velo(EDGE*id_chd+DG+1))

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) .gt. 0) &
         velo(EDGE*id_par+UP+1) = 0.5_8*(velo(EDGE*id_chd+UP+1) + velo(EDGE*idN_chd+UP+1))
  end subroutine restrict_velo

  subroutine check_m (dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par,  offs_chd, dims_chd)
    ! Check_m is an unused subroutine
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par, idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE
    real(8) :: ratio

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx(i_par, j_par, offs_par, dims_par)

    idE   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN2E = idx(i_chd + 2, j_chd + 1, offs_chd, dims_chd)
    id2NE = idx(i_chd + 1, j_chd + 2, offs_chd, dims_chd)
    idN   = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idW   = idx(i_chd - 1, j_chd,     offs_chd, dims_chd)
    idNW  = idx(i_chd - 1, j_chd + 1, offs_chd, dims_chd)
    idS2W = idx(i_chd - 2, j_chd - 1, offs_chd, dims_chd)
    idSW  = idx(i_chd - 1, j_chd - 1, offs_chd, dims_chd)
    idS   = idx(i_chd,     j_chd - 1, offs_chd, dims_chd)
    id2SW = idx(i_chd - 1, j_chd - 2, offs_chd, dims_chd)
    idSE  = idx(i_chd + 1, j_chd - 1, offs_chd, dims_chd)

    ratio = (1.0_8/dom%areas%elts(id_chd+1)%hex_inv + &
         dom%overl_areas%elts(idE+1)%a(1) + &
         dom%overl_areas%elts(idNE+1)%a(2) + &
         dom%overl_areas%elts(idN2E+1)%a(3) + &
         dom%overl_areas%elts(id2NE+1)%a(4) + &
         dom%overl_areas%elts(idN+1)%a(1) + &
         dom%overl_areas%elts(idW+1)%a(2) + &
         dom%overl_areas%elts(idNW+1)%a(3) + &
         dom%overl_areas%elts(idS2W+1)%a(4) + &
         dom%overl_areas%elts(idSW+1)%a(1) + &
         dom%overl_areas%elts(idS+1)%a(2) + &
         dom%overl_areas%elts(id2SW+1)%a(3) + &
         dom%overl_areas%elts(idSE+1)%a(4))*dom%areas%elts(id_par+1)%hex_inv
  end subroutine check_m

  subroutine restrict_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Restrict both mass and potential temperature
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    
    integer :: id_chd, id_par
   
    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
    
    if (dom%mask_n%elts(id_chd+1) .eq. 0) return

    id_par = idx(i_par, j_par, offs_par, dims_par)

    mass(id_par+1) = restrict_s (mass(id_chd+1), wc_m, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
    temp(id_par+1) = restrict_s (temp(id_chd+1), wc_t, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
  end subroutine restrict_scalar

  function restrict_s (scalar, wavelet, dom, id_par, i_chd, j_chd, offs_chd, dims_chd)
    ! Restriction operator at nodes: sub-sample and lift
    real(8)                        :: restrict_s
    real(8)                        :: scalar
    real(8), dimension(:), pointer :: wavelet
    type(Domain)                   :: dom
    integer                        :: id_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    
    integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE, d
   
    idE   = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE  = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN2E = idx(i_chd + 2, j_chd + 1, offs_chd, dims_chd)
    id2NE = idx(i_chd + 1, j_chd + 2, offs_chd, dims_chd)
    idN   = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)
    idW   = idx(i_chd - 1, j_chd,     offs_chd, dims_chd)
    idNW  = idx(i_chd - 1, j_chd + 1, offs_chd, dims_chd)
    idS2W = idx(i_chd - 2, j_chd - 1, offs_chd, dims_chd)
    idSW  = idx(i_chd - 1, j_chd - 1, offs_chd, dims_chd)
    idS   = idx(i_chd,     j_chd - 1, offs_chd, dims_chd)
    id2SW = idx(i_chd - 1, j_chd - 2, offs_chd, dims_chd)
    idSE  = idx(i_chd + 1, j_chd - 1, offs_chd, dims_chd)
    
    restrict_s = scalar + &
         (wavelet(idE+1)*dom%overl_areas%elts(idE+1)%a(1) + &
         wavelet(idNE+1)*dom%overl_areas%elts(idNE+1)%a(2) + &
         wavelet(idN2E+1)*dom%overl_areas%elts(idN2E+1)%a(3) + &
         wavelet(id2NE+1)*dom%overl_areas%elts(id2NE+1)%a(4) + &
         wavelet(idN+1)*dom%overl_areas%elts(idN+1)%a(1) + &
         wavelet(idW+1)*dom%overl_areas%elts(idW+1)%a(2) + &
         wavelet(idNW+1)*dom%overl_areas%elts(idNW+1)%a(3) + &
         wavelet(idS2W+1)*dom%overl_areas%elts(idS2W+1)%a(4) + &
         wavelet(idSW+1)*dom%overl_areas%elts(idSW+1)%a(1) + &
         wavelet(idS+1)*dom%overl_areas%elts(idS+1)%a(2) + &
         wavelet(id2SW+1)*dom%overl_areas%elts(id2SW+1)%a(3) + &
         wavelet(idSE+1)*dom%overl_areas%elts(idSE+1)%a(4))* &
         dom%areas%elts(id_par+1)%hex_inv
  end function restrict_s

  type(Iu_Wgt) function outer_velo_weights (dom, p, i0, j0, e0, offs, dims)
    type(Domain)                   :: dom
    integer                        :: p, i0, j0, e0
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    type(Coord)             :: x, y
    integer                 :: k, id, info
    real(8), dimension(9)   :: weights
    real(8), dimension(6)   :: b, ipiv
    real(8), dimension(6,6) :: G

    call local_coord (dom%midpt%elts(idx(i0, j0, offs, dims)*EDGE + e0 + 1), &
         dom%node%elts(idx(i0 + end_pt(1,1,e0+1), j0 + end_pt(2,1,e0+1), &
         offs, dims) + 1), dom%node%elts(idx(i0 + end_pt(1,2,e0+1), j0 + &
         end_pt(2,2,e0+1), offs, dims) + 1), x, y)

    weights = 0.0_8

    do k = 1, 2
       id = idx(i0, j0, offs, dims)

       G = reshape((/0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
            0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
            0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
            0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, 0.0_8, &
            0.0_8, 0.0_8, 0.0_8, 0.0_8/), (/6, 6/))

       G(:,1) = coords_to_row(i0, j0, (/0, 0/), (/0, 0, e0/), e0)

       G(:,2) = coords_to_row(i0, j0, end_pt(:,-k+3,e0+1), &
            hex_sides(:,hex_s_offs(e0+1) + 2 + 3*k - 3 + 1), e0)

       G(:,3) = coords_to_row(i0, j0, end_pt(:,k,e0+1), &
            hex_sides(:,(hex_s_offs(e0+1) + 3) - (3*k - 3) + 1), e0)

       G(:,4) = coords_to_row(i0, j0, opp_no(:,k,e0+1), &
            hex_sides(:,hex_s_offs(e0+1) + 1 + 3*k - 3 + 1), e0) - &
            coords_to_row(i0, j0, end_pt(:,k,e0+1), &
            hex_sides(:,hex_s_offs(e0+1) + 2 + 3*k - 3 + 1), e0)

       G(:,5) = coords_to_row(i0, j0, end_pt(:,-k+3,e0+1), &
            hex_sides(:,(hex_s_offs(e0+1) + 3) - (3*k - 3) + 1), e0) - &
            coords_to_row(i0, j0, opp_no(:,k,e0+1), &
            hex_sides(:,(hex_s_offs(e0+1) + 4) - (3*k - 3) + 1), e0)

       G(:,6) = coords_to_row(i0, j0, end_pt(:,k,e0+1), &
            hex_sides(:,(hex_s_offs(e0+1) + 5) - (3*k - 3) + 1), e0) - &
            coords_to_row(i0, j0, end_pt(:,-k+3,e0+1), &
            hex_sides(:,hex_s_offs(e0+1) + 0 + 3*k - 3 + 1), e0)

       b = coords_to_rowd(mid_pt(dom%midpt%elts(EDGE*id+e0+1), &
            dom%node%elts(idx2(i0, j0, end_pt(:,2,e0+1), offs, dims) + &
            1)), vector(dom%node%elts(idx2(i0, j0, end_pt(:,1,e0+1), &
            offs, dims) + 1), dom%node%elts(idx2(i0, j0, &
            end_pt(:,2,e0+1), offs, dims) + 1)), x, y)

       ipiv = 0
       info = 0

       call dgesv(6, 1, G, 6, ipiv, b, 6, info)

       weights(1) = weights(1) + b(1)
       weights(2*k:2*k + 1) = weights(2*k:2*k + 1) + b(2:3)
       weights(2*k + 4:2*k + 5) = b(4:5)
       weights(-2*k+6) = weights(-2*k+6) + b(6)
       weights(-2*k+7) = weights(-2*k+7) - b(6)
    end do

    outer_velo_weights = Iu_Wgt(0.5_8*weights - Iu_Base_Wgt)
  contains
    function coords_to_row (i00, j00, n_offs1, n_offs2, e00)
      real(8), dimension(6) :: coords_to_row
      integer               :: e00, i00, j00
      integer, dimension(2) :: n_offs1
      integer, dimension(3) :: n_offs2
      
      type(Coord) :: endpt1, endpt2
      integer     :: i, j, e
      
      i = i00 + n_offs1(1) + n_offs2(1)
      j = j00 + n_offs1(2) + n_offs2(2)

      e = n_offs2(3)

      endpt1 = get_coord(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), e00)
      endpt2 = get_coord(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), e00)

      coords_to_row = coords_to_rowd(mid_pt(endpt1, endpt2), vector(endpt1, endpt2), x, y)
    end function coords_to_row

    function get_coord (i, j, e)
      type(Coord) ::  get_coord
      integer     :: i, j, e

      if (i .eq. -1) then
         if (j .eq. -1 .and. is_penta(dom, p, IJMINUS - 1)) then
            if (e .eq. RT) then
               get_coord = dom%node%elts(nidx(LAST_BDRY, 0, IMINUS, &
                    offs, dims) + 1)
               return
            else
               if (e .eq. UP) then
                  get_coord = dom%node%elts(nidx(0, LAST_BDRY, JMINUS, &
                       offs, dims) + 1)
                  return
               end if
            end if
         else
            if (j .eq. PATCH_SIZE .and. is_penta(dom, p, IMINUSJPLUS - &
                 1)) then
               get_coord = dom%node%elts(nidx(0, 1, JPLUS, offs, dims) + &
                    1)
               return
            else
               get_coord = dom%node%elts(idx(i, j, offs, dims) + 1)
               return
            end if
         end if
      else
         if (i .eq. PATCH_SIZE .and. j .eq. -1 .and. is_penta(dom, p, &
              IPLUSJMINUS - 1)) then
            get_coord = dom%node%elts(nidx(1, 0, IPLUS, offs, dims) + 1)
            return
         else
            if (i .eq. PATCH_SIZE + 1 .and. j .eq. PATCH_SIZE + 1 .and. &
                 is_penta(dom, p, IJPLUS - 1)) then
               get_coord = dom%node%elts(nidx(1, 0, IJPLUS, offs, dims) &
                    + 1)
               return
            else
               get_coord = dom%node%elts(idx(i, j, offs, dims) + 1)
               return
            end if
         end if
      end if
    end function get_coord
  end function outer_velo_weights

  function coord2local(c, x, y)
    real(8), dimension(2) :: coord2local
    type(Coord)           :: c, x, y

    coord2local = (/inner(c, x), inner(c, y)/)
  end function coord2local

  subroutine init_wavelet_mod()
    logical :: initialized = .False.
    
    if (initialized) return ! initialize only once
    call init_shared_mod()
    call init_domain_mod()

    Iu_Base_Wgt = (/16.0_8, -1.0_8, 1.0_8, 1.0_8, -1.0_8, -1.0_8, -1.0_8, 1.0_8, 1.0_8/)/16.0_8
    initialized = .True.
  end subroutine init_wavelet_mod

  function coords_to_rowd (midpt, dirvec, x, y)
    real(8), dimension(6) :: coords_to_rowd
    type(Coord)           :: dirvec, x, y

    type(Coord)           :: midpt
    real(8)               :: u, v
    real(8), dimension(2) :: xy

    call normalize2(coord2local(dirvec, x, y), u, v)
    xy = coord2local(midpt, x, y)
    coords_to_rowd = (/u, u*xy(1), u*xy(2), v, v*xy(1), v*xy(2)/)
  end function coords_to_rowd
end module wavelet_mod
