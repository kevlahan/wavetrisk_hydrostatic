module adapt_mod
  use refine_patch_mod
  use multi_level_mod
  implicit none
  real(8), parameter :: FILLUP_THRESHOLD = 0.9 ! Fillup patch to remove lower level if at least FILLUP_THRESHOLD*100% of nodes are active
  integer, parameter :: DOF_PER_PATCH = PATCH_SIZE*PATCH_SIZE*(EDGE+1)
  integer, parameter :: FILLED_AND_FROZEN = DOF_PER_PATCH + 1
  logical            :: max_level_exceeded
contains
  subroutine init_adapt_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_comm_mod
    call init_refine_patch_mod
    max_level_exceeded = .false.
    initialized = .true.
  end subroutine init_adapt_mod

  subroutine adapt_grid (set_thresholds, type)
    ! Grid adaptation during time stepping
    implicit none
    external          :: set_thresholds
    logical, optional :: type

    logical :: local_type

    if (present(type)) then
       local_type = type
    else
       local_type = .true.
    end if

    if (adapt_trend) then
       call trend_ml (sol, trend)
       call forward_wavelet_transform (trend, trend_wav_coeff)
    end if

    ! Find significant wavelets, adaptive grid and all masks
    call adapt (set_thresholds, local_type)
    
    call inverse_wavelet_transform (wav_coeff, sol)
  end subroutine adapt_grid 

  subroutine adapt (set_thresholds, type)
    ! Determines significant wavelets, adaptive grid and all masks associated with adaptive grid
    ! Assumes that trend wavelets have been calculated
    ! Does NOT transform solution and trend wavelets back onto adaptive grid
    implicit none
    external           :: set_thresholds
    logical, optional  :: type ! Recalculate thresholds
    
    integer :: k, l, d
    logical :: local_type

    ! Recalculate thresholds?
    if (present(type)) then
       local_type = type
    else
       local_type = .true.
    end if

    if (local_type) call set_thresholds
    
    call apply_onescale__int (set_masks, level_start, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, TOLRNZ)
    ! Initialize all other nodes and edges to ZERO
    do l = level_start+1, level_end
       call apply_onescale__int (set_masks, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, ZERO)
    end do

    ! Set all current masks > ADJZONE to at least ADJZONE for next time step
    call mask_adjacent_initial
    call comm_masks_mpi (NONE)

    ! Make nodes and edges with significant wavelet coefficients active
    if (adapt_trend) then
       call mask_active (trend_wav_coeff)
    else
       call mask_active (wav_coeff)
    end if
    call comm_masks_mpi (NONE)
    
    ! Add neighbouring parent wavelet nodes/edges
    do l = level_end-1, level_start, -1
       call apply_interscale (mask_adj_parent_nodes, l, z_null,  0, 1)
       call apply_interscale (mask_adj_parent_edges, l, z_null, -1, 1)
       call comm_masks_mpi (l)
    end do

    ! Add nearest neighbour wavelets of active nodes and edges at same scale
    do l = level_start, level_end
       call apply_onescale (mask_adj_same_scale, l, z_null, 0, 1)
    end do
    call comm_masks_mpi (NONE)

    ! needed if bdry is only 2 layers for scenario:
    ! mass > threshold @ PATCH_SIZE + 2 => flux restr @ PATCH_SIZE + 1
    ! => patch needed (contains flux for corrective part of R_F)
    do l = level_start+1, min(level_end, max_level-1)
       call apply_onescale (mask_restrict_flux, l, z_null, 0, 0)
    end do
    call comm_masks_mpi (NONE)
    
    ! Determine whether any new patches are required
    if (refine()) call post_refine

    if (perfect) then ! Apply perfect reconstruction check and exact TRiSK operators
       ! Add neighbouring wavelets at finer scale                                                                                         
       do l = level_end-1, level_start, -1
          call apply_interscale (mask_adj_children, l, z_null, 0, 0)
       end do
       call comm_masks_mpi (NONE)

       ! Ensure consistency of adjacent zones for nodes and edges                                                                         
       call mask_adj_nodes_edges
       call comm_masks_mpi (NONE)

       ! Ensure that perfect reconstruction criteria for active wavelets are satisfied                                                    
       do l = level_end-1, level_start-1, -1
          call apply_interscale (mask_perfect_scalar, l, z_null, 0, 0)
          call apply_interscale (mask_perfect_velo,   l, z_null, 0, 0)
          do d = 1, size(grid)
             call apply_to_penta_d (mask_perfect_velo_penta, grid(d), l, z_null)
          end do
          call comm_masks_mpi (l)
       end do
    else ! Apply old wavetrisk routines
       call complete_masks
    end if

    ! Add nodes and edges required for TRISK operators
    do l = level_start, level_end
       call apply_onescale (mask_trsk, l, z_null, 0, 0)
    end do
    call comm_masks_mpi (NONE)

    ! Label points required for remap as TRSK
    do l = level_start, level_end
       call apply_onescale (mask_remap, l, z_null, -1, 1)
    end do
    call comm_masks_mpi (NONE)
    
    ! Set insignificant wavelet coefficients to zero
    if (local_type) then
       do k = 1, zlevels
          do l = level_start+1, level_end
             do d = 1, size(grid)
                wc_m => wav_coeff(S_MASS,k)%data(d)%elts
                wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
                wc_u => wav_coeff(S_VELO,k)%data(d)%elts
                call apply_onescale_d (compress, grid(d), l, z_null, 0, 1)
                nullify (wc_m, wc_t, wc_u)
                wc_m => trend_wav_coeff(S_MASS,k)%data(d)%elts
                wc_t => trend_wav_coeff(S_TEMP,k)%data(d)%elts
                wc_u => trend_wav_coeff(S_VELO,k)%data(d)%elts
                call apply_onescale_d (compress, grid(d), l, z_null, 0, 1)
                nullify (wc_m, wc_t, wc_u)
             end do
          end do
       end do
    end if
    wav_coeff%bdry_uptodate = .false.
    trend_wav_coeff%bdry_uptodate = .false.
  end subroutine adapt

   subroutine WT_after_step (q, wav, l_start0)
    !  Everything needed in terms of forward and backward wavelet transform
    !  after one time step (e.g. RK sub-step)
    !    A) compute wavelets and perform backwards transform to conserve mass
    !    B) interpolate values onto adapted grid for next step
    implicit none
    type(Float_Field), dimension(S_MASS:S_VELO,1:zlevels), target :: q, wav
    integer, optional                                             :: l_start0
    
    integer :: d, j, k, l, l_start

    if (present(l_start0)) then
       l_start = l_start0
       if (max_level > min_level) then
          do k = 1, zlevels
             do d = 1, size(grid)
                velo => q(S_VELO,k)%data(d)%elts
                call apply_interscale_d (restrict_velo, grid(d), level_start-1, k, 0, 0)
                nullify (velo)
             end do
          end do
       end if
    else
       l_start = level_start
    end if

    q%bdry_uptodate = .false.
    call update_array_bdry (q, NONE)

    do k = 1, zlevels
       do l = l_start, level_end-1
          do d = 1, size(grid)
             mass => q(S_MASS,k)%data(d)%elts
             temp => q(S_TEMP,k)%data(d)%elts
             velo => q(S_VELO,k)%data(d)%elts
             
             wc_m => wav(S_MASS,k)%data(d)%elts
             wc_t => wav(S_TEMP,k)%data(d)%elts
             wc_u => wav(S_VELO,k)%data(d)%elts
             
             call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             call apply_interscale_d (compute_velo_wavelets,   grid(d), l, z_null, 0, 0)
             call apply_to_penta_d (compute_velo_wavelets_penta, grid(d), l, z_null)
             nullify (mass, temp, velo, wc_m, wc_t, wc_u)
          end do
          wav(:,k)%bdry_uptodate = .false.
       end do

       do l = level_start+1, level_end
          do d = 1, size(grid)
             wc_m => wav(S_MASS,k)%data(d)%elts
             wc_t => wav(S_TEMP,k)%data(d)%elts
             wc_u => wav(S_VELO,k)%data(d)%elts
             call apply_onescale_d (compress, grid(d), l, k, 0, 1)
             nullify (wc_m, wc_t, wc_u)
          end do
          wav(:,k)%bdry_uptodate = .false.
       end do
    end do

    call inverse_wavelet_transform (wav, q)
  end subroutine WT_after_step

  subroutine compress (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: e, id, k, v

    id = idx(i, j, offs, dims)
    
    if (dom%mask_n%elts(id+1) < ADJZONE) then
       wc_m(id+1) = 0.0_8
       wc_t(id+1) = 0.0_8
    end if
    
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) < ADJZONE) wc_u(EDGE*id+e) = 0.0_8
    end do
  end subroutine compress

  logical function refine ()
    implicit none
    ! Determines where new patches are needed
    integer :: c, d, did_refine, old_n_patch, p_chd, p_par
    logical :: required
    
    ! Use threshold masks call refine patch where necessary
    did_refine = FALSE
    do d = 1, size(grid)
       old_n_patch = grid(d)%patch%length
       do p_par = 2, grid(d)%patch%length
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             required = check_child_required (grid(d), p_par - 1, c - 1)
             if (required .and. p_chd <= 0) then ! New patch required and does not yet exist
                if (grid(d)%patch%elts(p_par)%level == max_level) then ! Cannot refine further
                   max_level_exceeded = .true.
                else
                   call refine_patch1 (grid(d), p_par - 1, c - 1)
                   did_refine = TRUE
                end if
             end if
          end do
       end do
       do p_par = 2, old_n_patch
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             if (p_chd+1 > old_n_patch) call refine_patch2 (grid(d), p_par - 1, c - 1)
          end do
       end do
    end do
    refine = sync_max(did_refine) == TRUE
    return
  end function refine

  subroutine patch_count_active (dom, p)
    implicit none
    type(Domain) :: dom
    integer      :: p

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: e, i, id, j

    ! TODO set FILLED_AND_FROZEN in `remove_inside_patches` to save work here
    if (dom%patch%elts(p+1)%active == FILLED_AND_FROZEN) return

    dom%patch%elts(p+1)%active = 0
    call get_offs_Domain (dom, p, offs, dims)
    do j = 1, PATCH_SIZE
       do i = 1, PATCH_SIZE
          id = idx(i-1, j-1, offs, dims)
          if (dom%mask_n%elts(id+1) >= ADJZONE) dom%patch%elts(p+1)%active = dom%patch%elts(p+1)%active + 1
          do e = 1, EDGE
             if (dom%mask_e%elts(EDGE*id+e) >= ADJZONE) dom%patch%elts(p+1)%active = dom%patch%elts(p+1)%active + 1
          end do
       end do
    end do
  end subroutine patch_count_active

  function get_child_and_neigh_patches (dom, p_par, c)
    implicit none
    integer, dimension(4) :: get_child_and_neigh_patches
    type(Domain)          :: dom
    integer               :: p_par, c
    
    integer :: n

    get_child_and_neigh_patches = 0
    get_child_and_neigh_patches(1) = dom%patch%elts(p_par+1)%children(c)
    n = dom%patch%elts(p_par+1)%neigh(c) ! side
    if (n > 0) then
       get_child_and_neigh_patches(2) = dom%patch%elts(n+1)%children(modulo((c+1)-1,4)+1) 
       get_child_and_neigh_patches(3) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
    endif

    n = dom%patch%elts(p_par+1)%neigh(c+4) ! corner
    if (n > 0) get_child_and_neigh_patches(4) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
  end function get_child_and_neigh_patches

  real(8) function check_children_fillup (dom, p_par)
    implicit none
    type(Domain) :: dom
    integer      :: p_par
    
    integer               :: active, p_chd, c, i
    integer, dimension(4) :: cn

    ! do not fill up children + remove if any child or neighbour patch is missing
    if (product(dom%patch%elts(p_par+1)%neigh)*product(dom%patch%elts(p_par+1)%children) == 0) then
       check_children_fillup = -1.0d8 ! for now make removing imposible for this case
       return
    end if

    active = 0
    do c = 1, N_CHDRN
       cn = get_child_and_neigh_patches(dom, p_par, c)
       do i = 2, 4
          active = active + dom%patch%elts(cn(i)+1)%active
       end do
       ! count real child double since neighbours are counted double during accumulation
       active = active + 2*dom%patch%elts(cn(1)+1)%active
    end do
    check_children_fillup = dble(active)/dble(N_CHDRN*5*DOF_PER_PATCH)
  end function check_children_fillup

  logical function remove_inside_patches ()
    ! Removes patches that are not required because they are far enough away from the locally finest level
    
    integer               :: d, k, p, l, c, c1
    integer, dimension(4) :: chdrn
    real(8)               :: children_fullness
    logical               :: changes

    changes = .false.
    do d = 1, size(grid)
       do p = 2, grid(d)%patch%length
          call patch_count_active (grid(d), p-1)
       end do
    end do
    l = level_start - 1
    do k = 1, grid(d)%lev(l)%length
       do d = 1, size(grid)
          p = grid(d)%lev(l)%elts(k)
          children_fullness = 0
          ! Patch has 4 children including first neighbours 16
          do c = 1, N_CHDRN
             ! Always do one child with 3 adjacent neighbours (2 on side 1 on corner)
             chdrn = get_child_and_neigh_patches(grid(d), p, c)
             do c1 = 1, N_CHDRN
                children_fullness = children_fullness + check_children_fillup(grid(d), chdrn(c1))
             end do
          end do
          children_fullness = children_fullness/16
          if (children_fullness > FILLUP_THRESHOLD) then
             do c = 1, N_CHDRN
                chdrn = get_child_and_neigh_patches(grid(d), p, c)
                do c1 = 1, N_CHDRN
                   call apply_onescale_to_patch__int (set_masks, grid(d), chdrn(c1), z_null, 0, 0, FROZEN)
                end do
             end do
             grid(d)%patch%elts(p+1)%active = NONE
             changes = .true.
          else
             grid(d)%patch%elts(p+1)%active = NONE
          end if
       end do
    end do
    remove_inside_patches = changes
  end function remove_inside_patches

  logical function check_child_required (dom, p, c)
    implicit none
    type(Domain) :: dom
    integer      :: p, c

    integer                        :: e, j0, j, i0, i, id, st, en
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    logical                        :: required

    st = -BDRY_THICKNESS
    en =  BDRY_THICKNESS

    call get_offs_Domain (dom, p, offs, dims)

    do j0 = st + 1, PATCH_SIZE/2 + en
       j = j0 - 1 + chd_offs(2,c+1)
       do i0 = st + 1, PATCH_SIZE/2 + en
          i = i0 - 1 + chd_offs(1,c+1)
          id = idx(i, j, offs, dims)
          required = dom%mask_n%elts(id+1) >= ADJSPACE .or. dom%mask_n%elts(id+1) == TRSK
          do e = 1, EDGE
             required = required .or. dom%mask_e%elts(EDGE*id+e) >= RESTRCT
          end do
          if (required) then
             check_child_required = .true.
             return
          end if
       end do
    end do
    check_child_required = .false.
  end function check_child_required

  subroutine init_multi_level_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_comm_mod
    call init_ops_mod
    call init_wavelet_mod
    call init_refine_patch_mod
    initialized = .true.
  end subroutine init_multi_level_mod

  subroutine add_second_level
    implicit none
    integer :: d, c

    do d = 1, size(grid)
       do c = 1, N_CHDRN
          call refine_patch1 (grid(d), 1, c-1)
       end do
       do c = 1, N_CHDRN
          call refine_patch2 (grid(d), 1, c-1)
       end do
       call connect_children (grid(d), 1)
    end do

    call comm_patch_conn_mpi

    do d = 1, size(grid)
       call update_comm (grid(d))
    end do

    call comm_communication_mpi
    call comm_nodes9_mpi (get_areas, set_areas, NONE)
    call apply_to_penta (area_post_comm, NONE, z_null)
  end subroutine add_second_level
  
  subroutine fill_up_level
    ! Fills up level level_start+1 and increases level_start
    implicit none
    integer :: d, j, p_par, c, p_chd

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(level_start)%length
          p_par = grid(d)%lev(level_start)%elts(j)
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd == 0) call refine_patch (grid(d), p_par, c - 1)
          end do
       end do
    end do
    call post_refine
    level_start = level_start+1
  end subroutine fill_up_level

  subroutine fill_up_grid_and_IWT (l)
    ! Fills grid up to level l and does inverse wavelet transform of solution onto grid
    implicit none
    integer :: l
    
    integer :: old_level_start
    
    old_level_start = level_start
    do while (level_start < l)
       call fill_up_level
    end do
    call inverse_wavelet_transform (wav_coeff, sol, old_level_start)
    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE)
    level_start = old_level_start
  end subroutine fill_up_grid_and_IWT
end module adapt_mod
