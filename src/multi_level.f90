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
    ! Compute flux restriction of mass and potential temperature by summing coarse, corrective and small fluxes
    type(Domain) :: dom
    integer :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY + 1) :: offs_par
    integer, dimension(2,N_BDRY + 1) :: dims_par
    integer, dimension(N_BDRY + 1) :: offs_chd
    integer, dimension(2,N_BDRY + 1) :: dims_chd
    integer :: id_par, id_chd, idN_par, idE_par, idNE_par
    real(8), dimension(4) :: sm_flux_m, sm_flux_t

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
       sm_flux_m = interp_flux(h_mflux, dom, i_chd, j_chd, offs_chd, dims_chd)
       sm_flux_t = interp_flux(h_tflux, dom, i_chd, j_chd, offs_chd, dims_chd)
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) .ge. RESTRCT) then
       h_mflux(EDGE*id_par+RT+1) = &
            complete_coarse_flux(dmass, h_mflux, sm_flux_m, dom, i_par, j_par, i_chd, j_chd, RT, offs_chd, dims_chd)
       h_tflux(EDGE*id_par+RT+1) = &
            complete_coarse_flux(dtemp, h_tflux, sm_flux_t, dom, i_par, j_par, i_chd, j_chd, RT, offs_chd, dims_chd)
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) .ge. RESTRCT) then
       h_mflux(DG+EDGE*id_par+1) = &
            complete_coarse_flux(dmass, h_mflux, sm_flux_m, dom, i_par, j_par, i_chd, j_chd, DG, offs_chd, dims_chd)
       h_tflux(DG+EDGE*id_par+1) = &
            complete_coarse_flux(dtemp, h_tflux, sm_flux_t, dom, i_par, j_par, i_chd, j_chd, DG, offs_chd, dims_chd)
    end if

    if (dom%mask_e%elts(EDGE*id_par+UP+1) .ge. RESTRCT) then
       h_mflux(EDGE*id_par+UP+1) = &
            complete_coarse_flux(dmass, h_mflux, sm_flux_m, dom, i_par, j_par, i_chd, j_chd, UP, offs_chd, dims_chd)
       h_tflux(EDGE*id_par+UP+1) = &
            complete_coarse_flux(dtemp, h_tflux, sm_flux_t, dom, i_par, j_par, i_chd, j_chd, UP, offs_chd, dims_chd)
    end if
  contains
    function complete_coarse_flux(scalar, flux, sm_flux, dom, i_par, j_par, i_chd, j_chd, e, offs_chd, dims_chd)
      real(8) :: complete_coarse_flux
      real(8), dimension(:), pointer :: scalar, flux
      real(8), dimension(4) :: sm_flux
      type(Domain) :: dom
      integer :: i_par, j_par, i_chd, j_chd, e
      integer, dimension(N_BDRY + 1) :: offs_chd
      integer, dimension(2,N_BDRY + 1) :: dims_chd

      real(8) :: p_flux, c_flux

      if (e .eq. RT) then
         p_flux = part_flux(scalar, flux, dom, i_chd+1, j_chd, RT, offs_chd, dims_chd)
         c_flux = coarse_flux(scalar, dom, i_par, j_par, i_chd+1, j_chd, RT)
         complete_coarse_flux = p_flux + c_flux + sm_flux(1) + sm_flux(2)
      elseif (e .eq. DG) then
          p_flux = part_flux(scalar, flux, dom, i_chd+1, j_chd+1, DG, offs_chd, dims_chd)
          c_flux = coarse_flux(scalar, dom, i_par, j_par, i_chd+1, j_chd+1, DG)
          complete_coarse_flux = p_flux + c_flux + sm_flux(2) + sm_flux(3)
      elseif (e .eq. UP) then
         p_flux = part_flux(scalar, flux, dom, i_chd, j_chd+1, UP, offs_chd, dims_chd)
         c_flux = coarse_flux(scalar, dom, i_par, j_par, i_chd, j_chd+1, UP)
         complete_coarse_flux = p_flux + c_flux + sm_flux(3) + sm_flux(4)
      end if
    end function complete_coarse_flux

    function coarse_flux (scalar, dom, i_par, j_par, i_chd, j_chd, e)
      real(8) :: coarse_flux
      real(8), dimension(:), pointer :: scalar
      type(Domain) :: dom
      integer :: i_par, j_par, i_chd, j_chd, e
      
      integer :: id_mz, id_pz, id_mp,  id_pp, id_pm, id_mm, id_mm2, id_pm2, id_pp2,  id_mp2, id
          
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
      
      coarse_flux = (dom%overl_areas%elts(id+1)%a(1)*dom%overl_areas%elts(id+1)%a(2)*dom%areas%elts(id+1)%hex_inv &
           + dom%overl_areas%elts(id_mp+1)%a(2)*dom%overl_areas%elts(id_mp+1)%a(3)*dom%areas%elts(id_mp+1)%hex_inv &
           + dom%overl_areas%elts(id_pp+1)%a(1)*dom%overl_areas%elts(id_pp+1)%a(3)*dom%areas%elts(id_pp+1)%hex_inv &
           + dom%overl_areas%elts(id_pm+1)%a(1)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           + dom%overl_areas%elts(id_mm+1)%a(2)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv) &
           *(scalar(id_pz+1) - scalar(id_mz+1)) + &
           dom%overl_areas%elts(id_pp+1)%a(3)*dom%overl_areas%elts(id_pp+1)%a(4)*dom%areas%elts(id_pp+1)%hex_inv &
           *0.5_8*(scalar(id_pp2+1) - scalar(id_mz+1)) + &
           dom%overl_areas%elts(id_pm+1)%a(3)*dom%overl_areas%elts(id_pm+1)%a(4)*dom%areas%elts(id_pm+1)%hex_inv &
           *0.5_8*(scalar(id_pm2+1) - scalar(id_mz+1)) + &
           dom%overl_areas%elts(id_mp+1)%a(3)*dom%overl_areas%elts(id_mp+1)%a(4)*dom%areas%elts(id_mp+1)%hex_inv &
           *0.5_8*(scalar(id_pz+1) - scalar(id_mp2+1)) + &
           dom%overl_areas%elts(id_mm+1)%a(3)*dom%overl_areas%elts(id_mm+1)%a(4)*dom%areas%elts(id_mm+1)%hex_inv &
           *0.5_8*(scalar(id_pz+1) - scalar(id_mm2+1))
    end function coarse_flux
    
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

    function interp_flux (flux, dom, i, j, offs, dims)
      real(8), dimension(4) :: interp_flux
      real(8), dimension(:), pointer :: flux
      type(Domain) :: dom
      integer ::  i, j
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims

      real(4) :: wgt
      integer, dimension(20) :: id

      call get_indices(dom, i+1, j, RT, offs, dims, id)

      interp_flux(1) = - sum(flux(id((/WPM,UZM,VMM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j-2, offs, dims)+1)%enc) &
           - sum((flux(id((/VPM,WMMM,UMZ/)+1)+1) -flux(id((/UPZ,VPMM,WMM/)+1)+1)) * &
           dom%R_F_wgt%elts(idx(i+1,j-1, offs, dims)+1)%enc) ! UPLT S

      interp_flux(2) = sum(flux(id((/WMP,UZP,VPP/)+1)+1)* dom%R_F_wgt%elts(idx(i  ,j, offs, dims)+1)%enc) &
           + sum((flux(id((/VMP,WPPP,UPZ/)+1)+1) - flux(id((/UMZ,VMPP,WPP/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i  ,j+1, offs, dims)+1)%enc) ! LORT

      call get_indices(dom, i, j+1, UP, offs, dims, id)

      interp_flux(3) = - sum(flux(id((/UZM,VMM,WPM/)+1)+1) * dom%R_F_wgt%elts(idx(i+1,j, offs, dims)+1)%enc) &
           - sum((flux(id((/WMMM,UMZ,VPM/)+1)+1) - flux(id((/VPMM,WMM,UPZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i+1,j+1, offs, dims)+1)%enc) ! UPLT
     
      interp_flux(4) = sum(flux(id((/UZP,VPP,WMP/)+1)+1) * dom%R_F_wgt%elts(idx(i-2,j, offs, dims)+1)%enc) &
           + sum((flux(id((/WPPP,UPZ,VMP/)+1)+1) - flux(id((/VMPP,WPP,UMZ/)+1)+1))* &
           dom%R_F_wgt%elts(idx(i-2,j+1, offs, dims)+1)%enc) ! LORT W
    end function interp_flux

    function part_flux(scalar, flux, dom, i, j, e, offs, dims)
      real(8) :: part_flux
      real(8), dimension(:), pointer :: scalar
      real(8), dimension(:), pointer :: flux
      type(Domain) :: dom
      integer :: i, j, e
      integer, dimension(N_BDRY + 1) :: offs
      integer, dimension(2,N_BDRY + 1) :: dims
      
      real(8), dimension(2) :: area
      real(8), dimension(4) :: ol_area
      integer, dimension(20) :: id

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
      
      part_flux = sum(flux(id((/UPZ,UMZ/)+1)+1)*area) - sum(flux(id((/VMM,WMP/)+1)+1))*area(2) &
           - sum(flux(id((/WPM,VPP/)+1)+1))*area(1) &
           + ol_area(3)*scalar(id(MP+1)+1) - ol_area(4)*scalar(id(PM+1)+1) &
           - ol_area(1)*scalar(id(PP+1)+1) + ol_area(2)*scalar(id(MM+1)+1)
    end function part_flux
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
    do k = 1, zlevels
       do d = 1, size(grid)
          mass    => q(S_MASS,k)%data(d)%elts
          velo    => q(S_VELO,k)%data(d)%elts
          temp    => q(S_TEMP,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(integrate_pressure_up, grid(d), grid(d)%lev(level_end)%elts(j), k, 0, 1)
          end do

          do j = 1, grid(d)%lev(level_end)%length
             call step1(grid(d), grid(d)%lev(level_end)%elts(j), k)
          end do
          call apply_to_penta_d(post_step1, grid(d), level_end, k)

          if (viscosity .ne. 0.0_8) then
             do j = 1, grid(d)%lev(level_end)%length
                call apply_onescale_to_patch(flux_grad_scalar, grid(d), grid(d)%lev(level_end)%elts(j), z_null, -1, 1)
             end do
          end if

          nullify(mass, velo, temp, h_mflux, h_tflux)
       end do
    end do

    if (level_start .lt. level_end) then
       call update_array_bdry__start(horiz_flux, level_end) ! <= comm flux (Jmax)
    end if

    ! compute mass and potential temperature trend
    do k = 1, zlevels
       do d = 1, size(grid)
          mass    =>  q(S_MASS,k)%data(d)%elts
          velo    =>  q(S_VELO,k)%data(d)%elts
          temp    =>  q(S_TEMP,k)%data(d)%elts
          dmass   => dq(S_MASS,k)%data(d)%elts
          dtemp   => dq(S_TEMP,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(scalar_trend, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 1)
          end do

          nullify(mass, velo, temp, dmass, dtemp, h_mflux, h_tflux)
       end do
    end do

    dq%bdry_uptodate         = .False.
    horiz_flux%bdry_uptodate = .False.

    if (level_start .lt. level_end) then
       call update_array_bdry__finish(horiz_flux, level_end) ! <= finish non-blocking communicate mass flux (Jmax)
       call update_array_bdry__start(dq(S_MASS:S_TEMP,:), level_end) ! <= start non-blocking communicate dmass (l+1)
    end if

    ! Compute velocity trend
    do k = 1, zlevels
       do d = 1, size(grid)
          if (advect_only) cycle
          
          mass    =>  q(S_MASS,k)%data(d)%elts
          velo    =>  q(S_VELO,k)%data(d)%elts
          temp    =>  q(S_TEMP,k)%data(d)%elts
          dvelo   => dq(S_VELO,k)%data(d)%elts
          h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
          h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

          do j = 1, grid(d)%lev(level_end)%length
             call apply_onescale_to_patch(du_Qperp, grid(d), grid(d)%lev(level_end)%elts(j), z_null, 0, 0)
          end do

          nullify(mass, velo, temp, dvelo, h_mflux, h_tflux)
       end do
    end do

    ! Calculate trend on coarser scales
    do l = level_end-1, level_start, -1
       call update_array_bdry__finish(dq(S_MASS:S_TEMP,:), l+1) ! <= finish non-blocking communicate dmass (l+1)
       do k = 1, zlevels
          do d = 1, size(grid)
             mass    =>  q(S_MASS,k)%data(d)%elts
             velo    =>  q(S_VELO,k)%data(d)%elts
             temp    =>  q(S_TEMP,k)%data(d)%elts
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch(integrate_pressure_up, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do

             do j = 1, grid(d)%lev(l)%length
                call step1(grid(d), grid(d)%lev(l)%elts(j), k)
             end do

             call apply_to_penta_d(post_step1, grid(d), l, k)

             ! Diffuse scalars
             if (viscosity .ne. 0.0_8) then
                do  j = 1, grid(d)%lev(l)%length
                   call apply_onescale_to_patch(flux_grad_scalar, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
                end do
             end if
             
             call cpt_or_restr_flux(grid(d), l)  ! <= compute flux(l) & use dmass (l+1)

             nullify(mass, velo, temp, dmass, dtemp, h_mflux, h_tflux)
          end do
       end do

       call update_array_bdry(horiz_flux, l)
       
       do k = 1, zlevels
          do d = 1, size(grid)
             dmass   => dq(S_MASS,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch(scalar_trend, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1) ! use flux (l) & cpt dmass (l)
             end do
             nullify(dmass, dtemp, h_mflux, h_tflux)
          end do
       end do

       dq(S_MASS:S_TEMP,:)%bdry_uptodate = .False.
       
       if (l .gt. level_start) then
          call update_array_bdry__start(dq(S_MASS:S_TEMP,:), l)  ! <= start non-blocking communicate dmass (l+1)
       end if
       
       if (advect_only) cycle

       do k = 1, zlevels
          do d = 1, size(grid)
             mass    => q(S_MASS,k)%data(d)%elts
             velo    => q(S_VELO,k)%data(d)%elts
             temp    => q(S_TEMP,k)%data(d)%elts
             dmass   => dq(S_MASS,k)%data(d)%elts
             dvelo   => dq(S_VELO,k)%data(d)%elts
             dtemp   => dq(S_TEMP,k)%data(d)%elts
             h_mflux => horiz_flux(S_MASS,k)%data(d)%elts
             h_tflux => horiz_flux(S_TEMP,k)%data(d)%elts

             call cpt_or_restr_Qperp(grid(d), l, k)

             nullify(mass, velo, temp, dmass, dvelo, dtemp, h_mflux, h_tflux)
          end do
       end do
          
       dq(S_VELO,:)%bdry_uptodate = .False.
    end do

    if (advect_only) return

    ! Add gradient terms at all scales
    do k = 1, zlevels
       do d = 1, size(grid)
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
