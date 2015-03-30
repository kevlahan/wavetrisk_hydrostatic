module multi_level_mod
  use comm_mod
  use ops_mod
  use viscous_mod
  use wavelet_mod
  use refine_patch_mod
  use comm_mpi_mod
  implicit none

contains
  subroutine init_multi_level_mod()
      logical :: initialized = .False.
      if (initialized) return ! initialize only once
      call init_comm_mod()
      call init_ops_mod()
      call init_wavelet_mod()
      call init_refine_patch_mod()
      initialized = .True.
  end subroutine init_multi_level_mod

  subroutine add_second_level()
      integer d, c
      do d = 1, size(grid)
          do c = 1, N_CHDRN
              call refine_patch1(grid(d), 1, c-1)
          end do
          do c = 1, N_CHDRN
              call refine_patch2(grid(d), 1, c-1)
          end do
          call connect_children(grid(d), 1)
      end do
      call comm_patch_conn_mpi()
      do d = 1, size(grid)
          call update_comm(grid(d))
      end do
      call comm_communication_mpi()
      call comm_nodes9_mpi(get_areas, set_areas, NONE)
      call apply_to_penta(area_post_comm, NONE)
  end subroutine

  subroutine fill_up_level()
  ! fills up level `level_start + 1` and increases `level_start`
      integer d, k, p_par, c, p_chd
      do d = 1, size(grid)
          do k = 1, grid(d)%lev(level_start)%length
              p_par = grid(d)%lev(level_start)%elts(k)
              do c = 1, N_CHDRN
                  p_chd = grid(d)%patch%elts(p_par+1)%children(c)
                  if (p_chd .eq. 0) then
                      call refine_patch(grid(d), p_par, c - 1)
                  end if
              end do
          end do
      end do
      call post_refine()
      level_start = level_start+1
  end subroutine

  subroutine post_refine()
      integer d, p
      level_end = sync_max(level_end)
      do d = 1, n_domain(rank+1)
          do p = 2, grid(d)%patch%length
              call connect_children(grid(d), p - 1)
          end do
      end do
      call comm_patch_conn_mpi()
      do d = 1, size(grid)
          call update_comm(grid(d))
      end do
      call comm_communication_mpi()
      call comm_nodes9_mpi(get_areas, set_areas, NONE)
      call apply_to_penta(area_post_comm, NONE)
  end subroutine

  subroutine flux_cpt_restr(dom, p_chd, i_par, j_par, i_chd, j_chd, offs_par, &
          dims_par, offs_chd, dims_chd)
      type(Domain) dom
      integer i_par
      integer j_par
      integer i_chd
      integer j_chd
      integer p_chd
      integer, dimension(N_BDRY + 1) :: offs_par
      integer, dimension(2,N_BDRY + 1) :: dims_par
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd
      integer id_par
      integer id_chd
      integer idN_par
      integer idE_par
      integer idNE_par
      real(8) sm_flux(4)
      id_par = idx(i_par, j_par, offs_par, dims_par)
      id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
      id_par = idx(i_par, j_par, offs_par, dims_par)
      idN_par = idx(i_par, j_par + 1, offs_par, dims_par)
      idE_par = idx(i_par + 1, j_par, offs_par, dims_par)
      idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)
      if (dom%mask_p%elts(id_par+1) .ge. RESTRCT) then
          dom%bernoulli%elts(id_par+1) = dom%bernoulli%elts(id_chd+1)
      end if
      if (i_chd .ge. PATCH_SIZE .or. j_chd .ge. PATCH_SIZE) then
          return
      end if
      if (maxval(dom%mask_u%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)) .ge. RESTRCT) then
          call interp_small_fluxes(dom, i_chd, j_chd, offs_chd, dims_chd, sm_flux)
      end if
      if (dom%mask_u%elts(EDGE*id_par+RT+1) .ge. RESTRCT) then
          tflux(EDGE*id_par+RT+1) = part_coarse_fluxes(dom, i_chd+1, j_chd, RT, offs_chd, dims_chd) &
                                  + corr_coarse_fluxes(dom, i_par, j_par, i_chd+1, j_chd, RT) &
                                  + sm_flux(1) + sm_flux(2)
      end if
      if (dom%mask_u%elts(DG+EDGE*id_par+1) .ge. RESTRCT) then
          tflux(DG+EDGE*id_par+1) = part_coarse_fluxes(dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd) &
                                  + corr_coarse_fluxes(dom, i_par, j_par, i_chd + 1, j_chd + 1, DG) &
                                  + sm_flux(2) + sm_flux(3)
      end if
      if (dom%mask_u%elts(EDGE*id_par+UP+1) .ge. RESTRCT) then
          tflux(EDGE*id_par+UP+1) = part_coarse_fluxes(dom, i_chd, j_chd + 1, UP, offs_chd, dims_chd) &
                                  + corr_coarse_fluxes(dom, i_par, j_par, i_chd, j_chd + 1, UP) &
                                  + sm_flux(3) + sm_flux(4)
      end if
  contains
      real(8) function corr_coarse_fluxes(dom, i_par, j_par, i_chd, j_chd, e)
          type(Domain) dom
          integer i_par
          integer j_par
          integer i_chd
          integer j_chd
          integer e
          integer id_mz
          integer id_pz
          integer id_mp
          integer id_pp
          integer id_pm
          integer id_mm
          integer id_mm2
          integer id_pm2
          integer id_pp2
          integer id_mp2
          integer id
          id_mz = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 1 + 1), &
                  offs_chd, dims_chd)
          id_pz = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 4 + 1), &
                  offs_chd, dims_chd)
          id_mp = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 0 + 1), &
                  offs_chd, dims_chd)
          id_pp = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 5 + 1), &
                  offs_chd, dims_chd)
          id_pm = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 3 + 1), &
                  offs_chd, dims_chd)
          id_mm = idx2(i_chd, j_chd, nghb_pt(:,hex_s_offs(e+1) + 2 + 1), &
                  offs_chd, dims_chd)
          id_mm2 = idx2(i_chd, j_chd, bfly_no2(:,1,e+1), offs_chd, dims_chd)
          id_pm2 = idx2(i_chd, j_chd, bfly_no2(:,2,e+1), offs_chd, dims_chd)
          id_pp2 = idx2(i_chd, j_chd, bfly_no2(:,3,e+1), offs_chd, dims_chd)
          id_mp2 = idx2(i_chd, j_chd, bfly_no2(:,4,e+1), offs_chd, dims_chd)
          id = idx(i_chd, j_chd, offs_chd, dims_chd)
          corr_coarse_fluxes = ( &
                  dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
                + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
                + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
                + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
                + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
                *(dheight(id_pz+1) - dheight(id_mz+1)) + &
                  dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
                  *0.5_8*(dheight(id_pp2+1) - dheight(id_mz+1)) + &
                  dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
                  *0.5_8*(dheight(id_pm2+1) - dheight(id_mz+1)) + &
                  dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
                  *0.5_8*(dheight(id_pz+1) - dheight(id_mp2+1)) + &
                  dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
                  *0.5_8*(dheight(id_pz+1) - dheight(id_mm2+1))
      end function
      subroutine get_indices(dom, i, j, e, offs, dims, id)
          type(Domain) dom
          integer i, j, e
          integer, dimension(N_BDRY + 1) :: offs
          integer, dimension(2,N_BDRY + 1) :: dims
          integer, dimension(2) :: ij_mp, ij_pp, ij_pm, ij_mm
          integer id(20)
          id(UMZ+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 1 + 1), offs, dims)
          id(UPZ+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 4 + 1), offs, dims)
          id(WMP+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 0 + 1), offs, dims)
          id(VPP+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 5 + 1), offs, dims)
          id(WPM+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 3 + 1), offs, dims)
          id(VMM+1) = ed_idx(i, j, hex_sides(:,hex_s_offs(e+1) + 2 + 1), offs, dims)
          ij_mp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 0 + 1)
          id(MP+1) = idx(ij_mp(1), ij_mp(2), offs, dims)
          ij_pp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 5 + 1)
          id(PP+1) = idx(ij_pp(1), ij_pp(2), offs, dims)
          ij_pm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 3 + 1)
          id(PM+1) = idx(ij_pm(1), ij_pm(2), offs, dims)
          ij_mm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1) + 2 + 1)
          id(MM+1) = idx(ij_mm(1), ij_mm(2), offs, dims)
          id(VMP+1) = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
          id(VMPP+1) = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1), offs, dims)
          id(UZP+1) = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 0 + 4 + 1), offs, dims)
          id(WPPP+1) = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
          id(WPP+1) = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1), offs, dims)
          id(VPM+1) = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1), offs, dims)
          id(VPMM+1) = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
          id(UZM+1) = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 3) - 2 + 1), offs, dims)
          id(WMMM+1) = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1), offs, dims)
          id(WMM+1) = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
      end subroutine
      subroutine interp_small_fluxes(dom, i, j, offs, dims, flux)
          type(Domain) dom
          integer i, j
          integer, dimension(N_BDRY + 1) :: offs
          integer, dimension(2,N_BDRY + 1) :: dims
          real(8) flux(4)
          real(4) wgt
          integer id(20)
          call get_indices(dom, i+1, j, RT, offs, dims, id)
          flux(1) = - sum(thickflux%data(dom%id+1)%elts(id((/WPM,UZM,VMM/)+1)+1)* &
                      dom%R_F_wgt%elts(idx(i+1,j-2, offs, dims)+1)%enc) &
                    - sum((thickflux%data(dom%id+1)%elts(id((/VPM,WMMM,UMZ/)+1)+1) &
                      -thickflux%data(dom%id+1)%elts(id((/UPZ,VPMM,WMM/)+1)+1))* &
                      dom%R_F_wgt%elts(idx(i+1,j-1, offs, dims)+1)%enc) ! UPLT S
          flux(2) = sum(thickflux%data(dom%id+1)%elts(id((/WMP,UZP,VPP/)+1)+1)* &
                    dom%R_F_wgt%elts(idx(i  ,j, offs, dims)+1)%enc) &
                  + sum((thickflux%data(dom%id+1)%elts(id((/VMP,WPPP,UPZ/)+1)+1) &
                    -thickflux%data(dom%id+1)%elts(id((/UMZ,VMPP,WPP/)+1)+1))* &
                    dom%R_F_wgt%elts(idx(i  ,j+1, offs, dims)+1)%enc) ! LORT
          call get_indices(dom, i, j+1, UP, offs, dims, id)
          flux(3) = - sum(thickflux%data(dom%id+1)%elts(id((/UZM,VMM,WPM/)+1)+1)* &
                      dom%R_F_wgt%elts(idx(i+1,j, offs, dims)+1)%enc) &
                    - sum((thickflux%data(dom%id+1)%elts(id((/WMMM,UMZ,VPM/)+1)+1) &
                      -thickflux%data(dom%id+1)%elts(id((/VPMM,WMM,UPZ/)+1)+1))* &
                      dom%R_F_wgt%elts(idx(i+1,j+1, offs, dims)+1)%enc) ! UPLT
          flux(4) = sum(thickflux%data(dom%id+1)%elts(id((/UZP,VPP,WMP/)+1)+1)* &
                    dom%R_F_wgt%elts(idx(i-2,j, offs, dims)+1)%enc) &
                  + sum((thickflux%data(dom%id+1)%elts(id((/WPPP,UPZ,VMP/)+1)+1) &
                    -thickflux%data(dom%id+1)%elts(id((/VMPP,WPP,UMZ/)+1)+1))* &
                    dom%R_F_wgt%elts(idx(i-2,j+1, offs, dims)+1)%enc) ! LORT W
      end subroutine
      real(8) function part_coarse_fluxes(dom, i, j, e, offs, dims)
          type(Domain) dom
          integer i, j, e
          integer, dimension(N_BDRY + 1) :: offs
          integer, dimension(2,N_BDRY + 1) :: dims
          real(8) area(2), ol_area(4)
          integer id(20)
          call get_indices(dom, i, j, e, offs, dims, id)
          area = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(1:2)
          ol_area(1:2) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%split
          ol_area(3:4) = dom%overl_areas%elts(idx(i, j, offs, dims) + 1)%a(3:4) - ol_area(1:2)
          area(1) = area(1) + ol_area(1) + ol_area(4)
          area(2) = area(2) + ol_area(2) + ol_area(3)
          area = area/sum(area)
          
          ol_area(1) = dom%overl_areas%elts(id(PP+1)+1)%split(1)
          ol_area(2) = dom%overl_areas%elts(id(MM+1)+1)%split(2)
          ol_area(3) = dom%overl_areas%elts(id(MP+1)+1)%a(3) - dom%overl_areas%elts(id(MP+1)+1)%split(1)
          ol_area(4) = dom%overl_areas%elts(id(PM+1)+1)%a(4) - dom%overl_areas%elts(id(PM+1)+1)%split(2)

          part_coarse_fluxes = &
                  sum(thickflux%data(dom%id+1)%elts(id((/UPZ,UMZ/)+1)+1)*area) &
                - sum(thickflux%data(dom%id+1)%elts(id((/VMM,WMP/)+1)+1))*area(2) &
                - sum(thickflux%data(dom%id+1)%elts(id((/WPM,VPP/)+1)+1))*area(1) &
                + ol_area(3)*dheight(id(MP+1)+1) - ol_area(4)*dheight(id(PM+1)+1) - &
                  ol_area(1)*dheight(id(PP+1)+1) + ol_area(2)*dheight(id(MM+1)+1)
      end function
  end subroutine

  subroutine cpt_or_restr_Qperp(dom, l)
      type(Domain) dom
      integer l
      integer k
      integer p_par
      integer c
      integer p_chd
      do k = 1, dom%lev(l)%length
          p_par = dom%lev(l)%elts(k)
          do c = 1, N_CHDRN
              p_chd = dom%patch%elts(p_par+1)%children(c)
              if (p_chd .eq. 0) then
                  call apply_onescale_to_patch(du_source, dom, p_par, 0, 0)
              end if
          end do
          call apply_interscale_to_patch(Qperp_cpt_restr, dom, &
                  dom%lev(l)%elts(k), 0, 0)
      end do
  end subroutine

  subroutine trend_ml(q, dq)
      integer d
      integer k
      integer l
      integer p
      type(Float_Field), target :: q(S_HEIGHT:S_VELO), dq(S_HEIGHT:S_VELO)
      call update_bdry(q(S_HEIGHT), NONE)
      call update_bdry(q(S_VELO), NONE)
      do d = 1, size(grid)
          height => q(S_HEIGHT)%data(d)%elts
          velo => q(S_VELO)%data(d)%elts
          tflux => thickflux%data(d)%elts
          do k = 1, grid(d)%lev(level_end)%length
              call step1(grid(d), grid(d)%lev(level_end)%elts(k))
          end do
          call apply_to_penta_d(post_step1, grid(d), level_end)
! uncomment to diffuse height
!         do k = 1, grid(d)%lev(level_end)%length
!             if (viscosity .ne. 0) call apply_onescale_to_patch(flux_gradP, grid(d), &
!                     grid(d)%lev(level_end)%elts(k), -1, 1)
!         end do
          nullify(height,velo,tflux)
      end do
      if (level_start .lt. level_end) then
          call update_bdry__start(thickflux, level_end) ! <= comm flux (Jmax)
      end if
      do d = 1, size(grid)
          height => q(S_HEIGHT)%data(d)%elts
          velo => q(S_VELO)%data(d)%elts
          tflux => thickflux%data(d)%elts
          dheight => dq(S_HEIGHT)%data(d)%elts
          do k = 1, grid(d)%lev(level_end)%length
              call apply_onescale_to_patch(height_trend, grid(d), &  ! <= cpt dp
                      grid(d)%lev(level_end)%elts(k), 0, 1)
          end do
          nullify(height,velo)
          nullify(dheight, tflux)
      end do
      dq(:)%bdry_uptodate = .False.
      thickflux%bdry_uptodate = .False.
      if (level_start .lt. level_end) then
          call update_bdry__finish(thickflux, level_end) ! <= comm flux (Jmax)
          call update_bdry__start(dq(S_HEIGHT), level_end)  ! <= comm dp (l+1)
      end if
      do d = 1, size(grid)
          if (advect_only) cycle
          dvelo => dq(S_VELO)%data(d)%elts
          velo => q(S_VELO)%data(d)%elts
          height => q(S_HEIGHT)%data(d)%elts
          tflux => thickflux%data(d)%elts
          do k = 1, grid(d)%lev(level_end)%length
              call apply_onescale_to_patch(du_source, grid(d), &
                      grid(d)%lev(level_end)%elts(k), 0, 0)
          end do
          nullify(velo, height, dvelo, tflux)
      end do
      do l = level_end-1, level_start, -1
          call update_bdry__finish(dq(S_HEIGHT), l+1)  ! <= comm dp (l+1)
          do d = 1, size(grid)
              dheight => dq(S_HEIGHT)%data(d)%elts
              height => q(S_HEIGHT)%data(d)%elts
              velo => q(S_VELO)%data(d)%elts
              tflux => thickflux%data(d)%elts
              do k = 1, grid(d)%lev(l)%length
                  p = grid(d)%lev(l)%elts(k)
                  call step1(grid(d), p)
              end do
              call apply_to_penta_d(post_step1, grid(d), l)
! uncomment to diffuse height
!             do k = 1, grid(d)%lev(l)%length
!                 p = grid(d)%lev(l)%elts(k)
!                 if (viscosity .ne. 0) call apply_onescale_to_patch(flux_gradP, grid(d), p, -1, 1)
!             end do
              call cpt_or_restr_flux(grid(d), l)  ! <= cpt flux(l) & use dp (l+1)
              nullify(height,velo,dheight,tflux)
          end do
          call update_bdry__start(thickflux, l)  ! <= comm flux (l)
          do d = 1, size(grid)
              velo => q(S_VELO)%data(d)%elts
              height => q(S_HEIGHT)%data(d)%elts
              do k = 1, grid(d)%lev(l)%length
!                 call apply_onescale_to_patch2(qv, grid(d), grid(d)%lev(l)%elts(k), -1, 1)
                  if (viscosity .ne. 0) call apply_onescale_to_patch(divu, grid(d), grid(d)%lev(l)%elts(k), 0, 1)
              end do
              nullify(velo, height)
          end do
          call update_bdry__finish(thickflux, l)  ! <= comm flux (l)
          do d = 1, size(grid)
              dheight => dq(S_HEIGHT)%data(d)%elts
              tflux => thickflux%data(d)%elts
              do k = 1, grid(d)%lev(l)%length
                  call apply_onescale_to_patch(height_trend, grid(d), grid(d)%lev(l)%elts(k), 0, 1) !!!!! use flux (l) & cpt dp (l)
              end do
              nullify(dheight, tflux)
          end do
          dq(S_HEIGHT)%bdry_uptodate = .False.
          if (l .gt. level_start) call update_bdry__start(dq(S_HEIGHT), l)  ! <= comm dp (l+1)
          if (advect_only) cycle
          do d = 1, size(grid)
              dvelo => dq(S_VELO)%data(d)%elts
              dheight => dq(S_HEIGHT)%data(d)%elts
              velo => q(S_VELO)%data(d)%elts
              height => q(S_HEIGHT)%data(d)%elts
              tflux => thickflux%data(d)%elts
              call cpt_or_restr_Qperp(grid(d), l)
              nullify(dvelo, dheight, tflux)
              nullify(velo, height)
          end do
          dq(S_VELO)%bdry_uptodate = .False.
      end do
      if (advect_only) return
      do d = 1, size(grid)
          dvelo => dq(S_VELO)%data(d)%elts
          do p = 2, grid(d)%patch%length
              call apply_onescale_to_patch(du_gradB, grid(d), p - 1, 0, 0)
          end do
          nullify(dvelo)
      end do
  end subroutine

  subroutine Qperp_cpt_restr(dom, i_par, j_par, i_chd, j_chd, offs_par, &
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
      integer id_chd
      integer idE_chd
      integer idNE_chd
      integer idN_chd
      id_par = idx(i_par, j_par, offs_par, dims_par)
      id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)
      idE_chd = idx(i_chd + 1, j_chd, offs_chd, dims_chd)
      idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
      idN_chd = idx(i_chd, j_chd + 1, offs_chd, dims_chd)
      if (minval(dom%mask_u%elts(EDGE*id_chd + RT + 1:EDGE*id_chd + UP + 1)) .lt. &
              ADJZONE) then
          call du_source(dom, i_par, j_par, offs_par, dims_par)
      end if
      if (dom%mask_u%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) then
          dvelo(EDGE*id_par+RT+1) = dvelo(EDGE*id_chd+RT+1) + dvelo(EDGE*idE_chd+RT+1)
      end if
      if (dom%mask_u%elts(DG+EDGE*id_chd+1) .ge. ADJZONE) then
          dvelo(DG+EDGE*id_par+1) = dvelo(DG+EDGE*idNE_chd+1) + dvelo(DG+EDGE*id_chd+1)
      end if
      if (dom%mask_u%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) then
          dvelo(EDGE*id_par+UP+1) = dvelo(EDGE*id_chd+UP+1) + dvelo(EDGE*idN_chd+UP+1)
      end if
  end subroutine

  subroutine cpt_or_restr_flux(dom, l)
      type(Domain) dom
      integer l
      integer k
      integer p_par
      integer c
      integer p_chd
      logical restrict(N_CHDRN)
      do k = 1, dom%lev(l)%length
          p_par = dom%lev(l)%elts(k)
          restrict = .False.
          do c = 1, N_CHDRN
              p_chd = dom%patch%elts(p_par+1)%children(c)
              if (p_chd .gt. 0) restrict(c) = .True.
          end do
          do c = 1, N_CHDRN
              if (restrict(c)) then
                  call apply_interscale_to_patch3(flux_cpt_restr, dom, p_par, c, 0, 1)
              end if
          end do
      end do
  end subroutine
end module multi_level_mod
