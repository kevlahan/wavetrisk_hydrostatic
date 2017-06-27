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
    call apply_to_penta(area_post_comm, NONE, z_null)
  end subroutine add_second_level

  subroutine fill_up_level()
    ! fills up level `level_start + 1` and increases `level_start`
    integer d, j, p_par, c, p_chd

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(level_start)%length
          p_par = grid(d)%lev(level_start)%elts(j)
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
  end subroutine fill_up_level

  subroutine post_refine()
    integer d, p

    level_end = sync_max(level_end)

    do d = 1, n_domain(rank+1)
       do p = 3, grid(d)%patch%length
          call connect_children(grid(d), p - 1)
       end do
    end do

    call comm_patch_conn_mpi()

    do d = 1, size(grid)
       call update_comm(grid(d))
    end do

    call comm_communication_mpi()
    call comm_nodes9_mpi(get_areas, set_areas, NONE)
    call apply_to_penta(area_post_comm, NONE, z_null)
  end subroutine post_refine

  subroutine flux_cpt_restr(dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    !compute flux restriction of mass and potential temperature by summing coarse, corrective and small fluxes
    type(Domain) dom
    integer i_par
    integer j_par
    integer i_chd
    integer j_chd
    integer zlev
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
    real(8) sm_flux_m(4), sm_flux_t(4), c_c_flux_m, c_c_flux_t, p_c_flux_m, p_c_flux_t

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    id_par   = idx(i_par,     j_par,     offs_par, dims_par)

    idN_par  = idx(i_par,     j_par + 1, offs_par, dims_par)
    idE_par  = idx(i_par + 1, j_par,     offs_par, dims_par)
    idNE_par = idx(i_par + 1, j_par + 1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) .ge. RESTRCT) then
       dom%bernoulli%elts(id_par+1) = dom%bernoulli%elts(id_chd+1)
       dom%exner%elts(id_par+1) = dom%exner%elts(id_chd+1)
    end if

    if (i_chd .ge. PATCH_SIZE .or. j_chd .ge. PATCH_SIZE) return

    if (maxval(dom%mask_e%elts(EDGE*id_par+RT+1:EDGE*id_par+UP+1)) .ge. RESTRCT) then
       call interp_small_fluxes(dom, i_chd, j_chd, offs_chd, dims_chd, sm_flux_m, sm_flux_t)
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) .ge. RESTRCT) then
       call part_coarse_fluxes(dom, i_chd+1, j_chd, RT, offs_chd, dims_chd, p_c_flux_m, p_c_flux_t)
       call corr_coarse_fluxes(dom, i_par, j_par, i_chd+1, j_chd, RT, c_c_flux_m, c_c_flux_t)

       h_mflux(EDGE*id_par+RT+1) = p_c_flux_m + c_c_flux_m + sm_flux_m(1) + sm_flux_m(2)
       h_tflux(EDGE*id_par+RT+1) = p_c_flux_t + c_c_flux_t + sm_flux_t(1) + sm_flux_t(2)
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) .ge. RESTRCT) then
       call part_coarse_fluxes(dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd, p_c_flux_m, p_c_flux_t)
       call corr_coarse_fluxes(dom, i_par, j_par, i_chd + 1, j_chd + 1, DG, c_c_flux_m, c_c_flux_t)

       h_mflux(DG+EDGE*id_par+1) = p_c_flux_m + c_c_flux_m + sm_flux_m(2) + sm_flux_m(3)
       h_tflux(DG+EDGE*id_par+1) = p_c_flux_t + c_c_flux_t + sm_flux_t(2) + sm_flux_t(3)
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1) .ge. RESTRCT) then
       call part_coarse_fluxes(dom, i_chd, j_chd + 1, UP, offs_chd, dims_chd, p_c_flux_m, p_c_flux_t)
       call corr_coarse_fluxes(dom, i_par, j_par, i_chd, j_chd + 1, UP, c_c_flux_m, c_c_flux_t)

       h_mflux(EDGE*id_par+UP+1) = p_c_flux_m + c_c_flux_m + sm_flux_m(3) + sm_flux_m(4)
       h_tflux(EDGE*id_par+UP+1) = p_c_flux_t + c_c_flux_t + sm_flux_t(3) + sm_flux_t(4)
    end if

  contains

    subroutine corr_coarse_fluxes(dom, i_par, j_par, i_chd, j_chd, e, flux_m, flux_t)
      !compute correction of coarse fluxes for mass and potential temperature
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
      real(8) flux_m, flux_t

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

      !mass
      flux_m = ( &
           dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
           + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
           + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
           + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
           *(dmass(id_pz+1) - dmass(id_mz+1)) + &
           dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
           *0.5_8*(dmass(id_pp2+1) - dmass(id_mz+1)) + &
           dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           *0.5_8*(dmass(id_pm2+1) - dmass(id_mz+1)) + &
           dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
           *0.5_8*(dmass(id_pz+1) - dmass(id_mp2+1)) + &
           dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
           *0.5_8*(dmass(id_pz+1) - dmass(id_mm2+1))

      !potential temperature
      flux_t = ( &
           dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
           + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
           + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
           + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
           *(dtemp(id_pz+1) - dtemp(id_mz+1)) + &
           dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
           *0.5_8*(dtemp(id_pp2+1) - dtemp(id_mz+1)) + &
           dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           *0.5_8*(dtemp(id_pm2+1) - dtemp(id_mz+1)) + &
           dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
           *0.5_8*(dtemp(id_pz+1) - dtemp(id_mp2+1)) + &
           dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
           *0.5_8*(dtemp(id_pz+1) - dtemp(id_mm2+1))

    end subroutine corr_coarse_fluxes

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

      id(MM+1)   = idx(ij_mm(1), ij_mm(2), offs, dims)

      id(VMP+1)  = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
      id(VMPP+1) = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1),   offs, dims)
      id(UZP+1)  = ed_idx(ij_mp(1), ij_mp(2), hex_sides(:,hex_s_offs(e+1) + 0 + 4 + 1),   offs, dims)
      id(WPPP+1) = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)
      id(WPP+1)  = ed_idx(ij_pp(1), ij_pp(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1),   offs, dims)
      id(VPM+1)  = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 4 + 1),   offs, dims)
      id(VPMM+1) = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 2 + 1), offs, dims)
      id(UZM+1)  = ed_idx(ij_pm(1), ij_pm(2), hex_sides(:,(hex_s_offs(e+1) + 3) - 2 + 1), offs, dims)
      id(WMMM+1) = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,hex_s_offs(e+1) + 1 + 2 + 1),   offs, dims)
      id(WMM+1)  = ed_idx(ij_mm(1), ij_mm(2), hex_sides(:,(hex_s_offs(e+1) + 4) - 4 + 1), offs, dims)

    end subroutine get_indices

    subroutine interp_small_fluxes(dom, i, j, offs, dims, flux_m, flux_t)
      type(Domain) dom
      integer i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      real(8) flux_m(4), flux_t(4)
      real(4) wgt
      integer id(20)

      call get_indices(dom, i+1, j, RT, offs, dims, id)

      flux_m(1) = - sum(h_mflux(id((/WPM,UZM,VMM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j-2, offs, dims)+1)%enc) &
           - sum((h_mflux(id((/VPM,WMMM,UMZ/)+1)+1) -h_mflux(id((/UPZ,VPMM,WMM/)+1)+1)) * &
           dom%R_F_wgt%elts(idx(i+1,j-1, offs, dims)+1)%enc) ! UPLT S

      flux_t(1) = - sum(h_tflux(id((/WPM,UZM,VMM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j-2, offs, dims)+1)%enc) &
           - sum((h_tflux(id((/VPM,WMMM,UMZ/)+1)+1) -h_tflux(id((/UPZ,VPMM,WMM/)+1)+1)) * &
           dom%R_F_wgt%elts(idx(i+1,j-1, offs, dims)+1)%enc) ! UPLT S

      flux_m(2) = sum(h_mflux(id((/WMP,UZP,VPP/)+1)+1)* dom%R_F_wgt%elts(idx(i  ,j, offs, dims)+1)%enc) &
           + sum((h_mflux(id((/VMP,WPPP,UPZ/)+1)+1) - h_mflux(id((/UMZ,VMPP,WPP/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i  ,j+1, offs, dims)+1)%enc) ! LORT

      flux_t(2) = sum(h_tflux(id((/WMP,UZP,VPP/)+1)+1)* dom%R_F_wgt%elts(idx(i  ,j, offs, dims)+1)%enc) &
           + sum((h_tflux(id((/VMP,WPPP,UPZ/)+1)+1) - h_tflux(id((/UMZ,VMPP,WPP/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i  ,j+1, offs, dims)+1)%enc) ! LORT

      call get_indices(dom, i, j+1, UP, offs, dims, id)

      flux_m(3) = - sum(h_mflux(id((/UZM,VMM,WPM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j, offs, dims)+1)%enc) &
           - sum((h_mflux(id((/WMMM,UMZ,VPM/)+1)+1) - h_mflux(id((/VPMM,WMM,UPZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i+1,j+1, offs, dims)+1)%enc) ! UPLT

      flux_t(3) = - sum(h_tflux(id((/UZM,VMM,WPM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j, offs, dims)+1)%enc) &
           - sum((h_tflux(id((/WMMM,UMZ,VPM/)+1)+1) - h_tflux(id((/VPMM,WMM,UPZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i+1,j+1, offs, dims)+1)%enc) ! UPLT

      flux_m(4) = sum(h_mflux(id((/UZP,VPP,WMP/)+1)+1) * dom%R_F_wgt%elts(idx(i-2,j, offs, dims)+1)%enc) &
           + sum((h_mflux(id((/WPPP,UPZ,VMP/)+1)+1) - h_mflux(id((/VMPP,WPP,UMZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i-2,j+1, offs, dims)+1)%enc) ! LORT W

      flux_t(4) = sum(h_tflux(id((/UZP,VPP,WMP/)+1)+1) * dom%R_F_wgt%elts(idx(i-2,j, offs, dims)+1)%enc) &
           + sum((h_tflux(id((/WPPP,UPZ,VMP/)+1)+1) - h_tflux(id((/VMPP,WPP,UMZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i-2,j+1, offs, dims)+1)%enc) ! LORT W
    end subroutine interp_small_fluxes

    subroutine part_coarse_fluxes(dom, i, j, e, offs, dims, flux_m, flux_t)
      type(Domain) dom
      integer i, j, e
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      real(8) area(2), ol_area(4)
      integer id(20)
      real(8) flux_m, flux_t

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

      !mass
      flux_m = &
           sum(h_mflux(id((/UPZ,UMZ/)+1)+1)*area) &
           - sum(h_mflux(id((/VMM,WMP/)+1)+1))*area(2) &
           - sum(h_mflux(id((/WPM,VPP/)+1)+1))*area(1) &
           + ol_area(3)*dmass(id(MP+1)+1) - ol_area(4)*dmass(id(PM+1)+1) - &
           ol_area(1)*dmass(id(PP+1)+1) + ol_area(2)*dmass(id(MM+1)+1)

      !potential temperature
      flux_t = &
           sum(h_tflux(id((/UPZ,UMZ/)+1)+1)*area) &
           - sum(h_tflux(id((/VMM,WMP/)+1)+1))*area(2) &
           - sum(h_tflux(id((/WPM,VPP/)+1)+1))*area(1) &
           + ol_area(3)*dtemp(id(MP+1)+1) - ol_area(4)*dtemp(id(PM+1)+1) - &
           ol_area(1)*dtemp(id(PP+1)+1) + ol_area(2)*dtemp(id(MM+1)+1)
    end subroutine part_coarse_fluxes

  end subroutine flux_cpt_restr

  subroutine cpt_or_restr_Qperp(dom, l, zlev)
    type(Domain) dom
    integer l
    integer j
    integer zlev
    integer p_par
    integer c
    integer p_chd

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .eq. 0) then
             call apply_onescale_to_patch(du_source, dom, p_par, zlev, 0, 0)
          end if
       end do
       call apply_interscale_to_patch(Qperp_cpt_restr, dom, dom%lev(l)%elts(j), z_null, 0, 0)
    end do

  end subroutine cpt_or_restr_Qperp

  subroutine trend_ml(q, dq)
    !compute mass, velocity and potential temperature trend
    integer d
    integer j
    integer k, l
    integer p
    type(Float_Field), target :: q(S_MASS:S_VELO,1:zlevels), dq(S_MASS:S_VELO,1:zlevels)

    call update_array_bdry(q, NONE)

    ! First integrate pressure down across all grid points in order to compute surface pressure
    do k = zlevels, 1, -1
       do d = 1, size(grid)
          mass => q(S_MASS,k)%data(d)%elts
          temp => q(S_TEMP,k)%data(d)%elts

          do p = 3, grid(d)%patch%length 
             call apply_onescale_to_patch(integrate_pressure_down, grid(d), p - 1, k, 0, 1)
          end do

          nullify(mass, temp)
       end do
    end do

    ! Calculate trend on finest scale
    do d = 1, size(grid)
       do k = 1, zlevels
          mass    => q(S_MASS,k)%data(d)%elts
          velo    => q(S_VELO,k)%data(d)%elts
          temp    => q(S_TEMP,k)%data(d)%elts
          h_mflux => horiz_massflux(k)%data(d)%elts
          h_tflux => horiz_tempflux(k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(integrate_pressure_up, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 1)
          end do

          do j = 1, grid(d)%lev(level_end)%length
             call step1(grid(d), grid(d)%lev(level_end)%elts(j), k)
          end do
          call apply_to_penta_d(post_step1, grid(d), level_end, k)

          nullify(mass, velo, temp, h_mflux, h_tflux)
       end do
    end do

    if (level_start .lt. level_end) then
       call update_vector_bdry__start(horiz_massflux, level_end) ! <= comm flux (Jmax)
       call update_vector_bdry__start(horiz_tempflux, level_end) ! <= comm flux (Jmax)
    end if

    ! compute mass and potential temperature trend
    do d = 1, size(grid)
       do k = 1, zlevels
          mass    =>  q(S_MASS,k)%data(d)%elts
          velo    =>  q(S_VELO,k)%data(d)%elts
          temp    =>  q(S_TEMP,k)%data(d)%elts
          dmass   => dq(S_MASS,k)%data(d)%elts
          dtemp   => dq(S_TEMP,k)%data(d)%elts
          h_mflux => horiz_massflux(k)%data(d)%elts
          h_tflux => horiz_tempflux(k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(masstemp_trend, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
          end do

          nullify(mass, velo, temp, dmass, dtemp, h_mflux, h_tflux)
       end do
    end do

    dq%bdry_uptodate             = .False.
    horiz_massflux%bdry_uptodate = .False.
    horiz_tempflux%bdry_uptodate = .False.

    if (level_start .lt. level_end) then
       call update_vector_bdry__finish(horiz_massflux, level_end) ! <= finish non-blocking communicate mass flux (Jmax)
       call update_vector_bdry__finish(horiz_tempflux, level_end) ! <= finish non-blocking communicate temp flux (Jmax)
       call update_array_bdry__start(dq(S_MASS:S_TEMP,:), level_end) ! <= start non-blocking communicate dmass (l+1)
    end if

       ! Compute velocity trend
    do d = 1, size(grid)
       do k = 1, zlevels
          if (advect_only) cycle
          
          mass    =>  q(S_MASS,k)%data(d)%elts
          velo    =>  q(S_VELO,k)%data(d)%elts
          temp    =>  q(S_TEMP,k)%data(d)%elts
          dvelo   => dq(S_VELO,k)%data(d)%elts
          h_mflux => horiz_massflux(k)%data(d)%elts
          h_tflux => horiz_tempflux(k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(du_Qperp, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 0)
          end do

          nullify(mass, velo, temp, dvelo, h_mflux, h_tflux)
       end do
    end do

       ! Calculate trend on coarser scales
    do l = level_end-1, level_start, -1
       call update_array_bdry__finish(dq(S_MASS:S_TEMP,:), l+1) ! <= finish non-blocking communicate dmass (l+1)

       do d = 1, size(grid)
          do k = 1, zlevels
             mass    =>  q(S_MASS,k)%data(d)%elts
             velo    =>  q(S_VELO,k)%data(d)%elts
             temp    =>  q(S_TEMP,k)%data(d)%elts
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_massflux(k)%data(d)%elts
             h_tflux => horiz_tempflux(k)%data(d)%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch(integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do

             do j = 1, grid(d)%lev(l)%length
                p = grid(d)%lev(l)%elts(j)
                call step1(grid(d), p, k)
             end do

             call apply_to_penta_d(post_step1, grid(d), l, k)
             call cpt_or_restr_flux(grid(d), l)  ! <= compute flux(l) & use dmass (l+1)

             nullify(mass, velo, temp, dmass, dtemp, h_mflux, h_tflux)
          end do
       end do

       call update_vector_bdry__start(horiz_massflux, l)  ! <= start non-blocking communicate flux (l)
       call update_vector_bdry__start(horiz_tempflux, l)  ! <= start non-blocking communicate flux (l)
          
       if (viscosity .ne. 0) then ! Calculate viscous term
          do d = 1, size(grid)
             do k = 1, zlevels
                velo => q(S_VELO,k)%data(d)%elts
                mass => q(S_MASS,k)%data(d)%elts
                
                do j = 1, grid(d)%lev(l)%length 
                   call apply_onescale_to_patch(divu, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
                end do
                
                nullify (velo, mass)
             end do
          end do
       end if
       
       call update_vector_bdry__finish(horiz_massflux, l)  ! <= finish non-blocking communicate flux (l)
       call update_vector_bdry__finish(horiz_tempflux, l)  ! <= finish non-blocking communicate flux (l)

       do d = 1, size(grid)
          do k = 1, zlevels
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_massflux(k)%data(d)%elts
             h_tflux => horiz_tempflux(k)%data(d)%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch(masstemp_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1) ! use flux (l) & cpt dmass (l)
             end do
             nullify(dmass, dtemp, h_mflux, h_tflux)
          end do
       end do

       dq(S_MASS:S_TEMP,:)%bdry_uptodate = .False.
       
       if (l .gt. level_start) then
          call update_array_bdry__start(dq(S_MASS:S_TEMP,:), l)  ! <= start non-blocking communicate dmass (l+1)
       end if
       
       if (advect_only) cycle

       do d = 1, size(grid)
          do k = 1, zlevels
             mass    => q(S_MASS,k)%data(d)%elts
             velo    => q(S_VELO,k)%data(d)%elts
             temp    => q(S_TEMP,k)%data(d)%elts
             dmass   => dq(S_MASS,k)%data(d)%elts
             dvelo   => dq(S_VELO,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_massflux(k)%data(d)%elts
             h_tflux => horiz_tempflux(k)%data(d)%elts

             call cpt_or_restr_Qperp(grid(d), l, k)

             nullify(mass, velo, temp, dmass, dvelo, dtemp, h_mflux, h_tflux)
          end do
       end do
          
       dq(S_VELO,:)%bdry_uptodate = .False.
    end do

    if (advect_only) return

    ! Add gradient terms at all scales
    do d = 1, size(grid)
       do k = 1, zlevels
          mass    =>  q(S_MASS,k)%data(d)%elts
          temp    =>  q(S_TEMP,k)%data(d)%elts
          dvelo   => dq(S_VELO,k)%data(d)%elts

          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch(du_gradB_gradExn, grid(d), p - 1, k, 0, 0)
          end do
          nullify(mass, temp, dvelo)
       end do
    end do
  end subroutine trend_ml

  subroutine Qperp_cpt_restr(dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
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
    integer id_chd
    integer idE_chd
    integer idNE_chd
    integer idN_chd

    id_par   = idx(i_par,     j_par,     offs_par, dims_par)

    id_chd   = idx(i_chd,     j_chd,     offs_chd, dims_chd)
    idE_chd  = idx(i_chd + 1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx(i_chd + 1, j_chd + 1, offs_chd, dims_chd)
    idN_chd  = idx(i_chd,     j_chd + 1, offs_chd, dims_chd)

    if (minval(dom%mask_e%elts(EDGE*id_chd + RT + 1:EDGE*id_chd + UP + 1)) .lt. &
         ADJZONE) then
       call du_source(dom, i_par, j_par, zlev, offs_par, dims_par)
    end if

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) .ge. ADJZONE) then
       dvelo(EDGE*id_par+RT+1) = dvelo(EDGE*id_chd+RT+1) + dvelo(EDGE*idE_chd+RT+1)
    end if

    if (dom%mask_e%elts(DG+EDGE*id_chd+1) .ge. ADJZONE) then
       dvelo(DG+EDGE*id_par+1) = dvelo(DG+EDGE*idNE_chd+1) + dvelo(DG+EDGE*id_chd+1)
    end if

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) .ge. ADJZONE) then
       dvelo(EDGE*id_par+UP+1) = dvelo(EDGE*id_chd+UP+1) + dvelo(EDGE*idN_chd+UP+1)
    end if

  end subroutine Qperp_cpt_restr

  subroutine cpt_or_restr_flux(dom, l)
    type(Domain) dom
    integer l
    integer j
    integer p_par
    integer c
    integer p_chd
    logical restrict(N_CHDRN)

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .False.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd .gt. 0) restrict(c) = .True.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3(flux_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do

  end subroutine cpt_or_restr_flux

end module multi_level_mod
