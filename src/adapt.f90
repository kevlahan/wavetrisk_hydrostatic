module adapt_mod
  use comm_mpi_mod
  use refine_patch_mod
  use multi_level_mod
  implicit none
  logical :: max_level_exceeded
  ! fillup patch to remove lower level if at least FILLUP_THRESHOLD*100% of nodes are active
  real(8), parameter :: FILLUP_THRESHOLD = 0.9
  integer, parameter :: DOF_PER_PATCH = PATCH_SIZE*PATCH_SIZE*(EDGE+1)
  integer, parameter :: FILLED_AND_FROZEN = DOF_PER_PATCH + 1

contains
  subroutine init_adapt_mod
    implicit none
    logical :: initialized = .False.
    if (initialized) return ! initialize only once
    call init_comm_mod
    call init_refine_patch_mod
    max_level_exceeded = .False.
    initialized = .True.
  end subroutine init_adapt_mod

  subroutine adapt_grid (set_thresholds)
    implicit none
    external :: set_thresholds

    if (adapt_trend) then
       call trend_ml (sol, trend)
       call forward_wavelet_transform (trend, trend_wav_coeff)
       call adapt (set_thresholds)
    else
       call adapt (set_thresholds)
    end if
    if (level_end .gt. level_start) call inverse_wavelet_transform (wav_coeff, sol, level_start)
  end subroutine adapt_grid

  subroutine adapt (set_thresholds)
    implicit none
    external :: set_thresholds
    integer :: k, l, d

    ! Ensure all nodes and edges at coarsest scales are active
    do l = level_start-1, level_start
       call apply_onescale__int (set_masks, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, TOLRNZ)
    end do

    ! Initialize all other nodes and edges to ZERO
    do l = level_start+1, level_end
       call apply_onescale__int (set_masks, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, ZERO)
    end do

    ! Set all current masks > ADJZONE to at least ADJZONE for next time step
    call mask_adjacent
    call comm_masks_mpi (NONE)
    
    if (adapt_trend) then 
       if (istep.eq.0) then ! Also adapt on variables when initializing
          call set_thresholds (1)
          call mask_active (wav_coeff)
       end if
       call set_thresholds (0)
       call mask_active (trend_wav_coeff)
    else
       call set_thresholds (1)
       call mask_active (wav_coeff)
    end if
    call comm_masks_mpi (NONE)
    
    do l = level_end-1, level_start, -1
       call apply_interscale (mask_active_nodes, l, z_null,  0, 1)
       call apply_interscale (mask_active_edges, l, z_null, -1, 1)
       call comm_masks_mpi (l)
    end do
    call comm_masks_mpi (NONE)

    do l = level_start, level_end
       call apply_onescale (mask_adj_space2, l, z_null, 0, 1)
    end do
    call comm_masks_mpi (NONE)

    ! needed if bdry is only 2 layers for scenario:
    ! mass > tol @ PATCH_SIZE + 2 => flux restr @ PATCH_SIZE + 1
    ! => patch needed (contains flux for corrective part of R_F)
    do l = level_start, min(level_end, max_level-1)
       call apply_onescale (mask_restrict_flux, l, z_null, 0, 0)
    end do
    call comm_masks_mpi (NONE)
    
    if (refine()) call post_refine
    call complete_masks

    do k = 1, zlevels
       do l = level_start+1, level_end
          call apply_onescale (compress, l, k, 0, 1)
       end do
    end do

    wav_coeff%bdry_uptodate = .False.
    if (adapt_trend) trend_wav_coeff%bdry_uptodate = .False.
  end subroutine adapt

  subroutine compress (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: e, id, k, v

    id = idx(i, j, offs, dims)
    
    if (dom%mask_n%elts(id+1) .lt. ADJZONE) then
       do v = S_MASS, S_TEMP
          wav_coeff(v,zlev)%data(dom%id+1)%elts(id+1) = 0.0_8
          if (adapt_trend) trend_wav_coeff(v,zlev)%data(dom%id+1)%elts(id+1) = 0.0_8
       end do
    end if
    
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) .lt. ADJZONE) then
          wav_coeff(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+e) = 0.0_8
          if (adapt_trend) trend_wav_coeff(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+e) = 0.0_8
       end if
    end do
  end subroutine compress

  logical function refine()
    logical :: required
    integer :: c, d, did_refine, old_n_patch, p_chd, p_par

    !  using tol masks call refine patch where necessary
    did_refine = FALSE
    do d = 1, size(grid)
       old_n_patch = grid(d)%patch%length
       do p_par = 2, grid(d)%patch%length
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             required = check_child_required(grid(d), p_par - 1, c - 1)
             if (p_chd .gt. 0) then
             else if (required) then ! new patch required
                if (grid(d)%patch%elts(p_par)%level .eq. max_level) then
                   max_level_exceeded = .True.
                else
                   call refine_patch1(grid(d), p_par - 1, c - 1)
                   did_refine = TRUE
                end if
             end if
          end do
       end do
       do p_par = 2, old_n_patch
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par)%children(c)
             if (p_chd+1 .gt. old_n_patch) then
                call refine_patch2(grid(d), p_par - 1, c - 1)
             end if
          end do
       end do
    end do
    refine = sync_max(did_refine) .eq. TRUE
    return
  end function refine

  subroutine patch_count_active (dom, p)
    type(Domain) :: dom
    integer      :: p

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer                        :: j, i, id, e

    ! TODO set FILLED_AND_FROZEN in `remove_inside_patches` to save work here
    if (dom%patch%elts(p+1)%active .eq. FILLED_AND_FROZEN) return

    dom%patch%elts(p+1)%active = 0
    call get_offs_Domain(dom, p, offs, dims)
    do j = 1, PATCH_SIZE
       do i = 1, PATCH_SIZE
          id = idx(i-1, j-1, offs, dims)
          if (dom%mask_n%elts(id+1) .ge. ADJZONE) dom%patch%elts(p+1)%active = dom%patch%elts(p+1)%active + 1
          do e = 1, EDGE
             if (dom%mask_e%elts(EDGE*id+e) .ge. ADJZONE) &
                  dom%patch%elts(p+1)%active = dom%patch%elts(p+1)%active + 1
          end do
       end do
    end do
  end subroutine patch_count_active

  function get_child_and_neigh_patches (dom, p_par, c)
    integer, dimension(4) :: get_child_and_neigh_patches
    type(Domain)          :: dom
    integer               :: p_par, c
    
    integer :: n

    get_child_and_neigh_patches = 0
    get_child_and_neigh_patches(1) = dom%patch%elts(p_par+1)%children(c)
    n = dom%patch%elts(p_par+1)%neigh(c) ! side
    if (n .gt. 0) then
       get_child_and_neigh_patches(2) = dom%patch%elts(n+1)%children(modulo((c+1)-1,4)+1) 
       get_child_and_neigh_patches(3) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
    endif

    n = dom%patch%elts(p_par+1)%neigh(c+4) ! corner
    if (n .gt. 0) then
       get_child_and_neigh_patches(4) = dom%patch%elts(n+1)%children(modulo((c+2)-1,4)+1) 
    endif
  end function get_child_and_neigh_patches

  function check_children_fillup (dom, p_par)
    real(8)      :: check_children_fillup
    type(Domain) :: dom
    integer      :: p_par
    
    integer               :: active, p_chd, c, i
    integer, dimension(4) :: cn

    ! do not fill up children + remove if any child or neighbour patch is missing
    if (product(dom%patch%elts(p_par+1)%neigh)*product(dom%patch%elts(p_par+1)%children) .eq. 0) then
       check_children_fillup = -1e8 ! for now make removing imposible for this case
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


  function remove_inside_patches()
    ! removes patches that are not required because they are far enough away from the locally finest level
    logical :: remove_inside_patches
    
    integer               :: d, k, p, l, c, c1
    integer, dimension(4) :: chdrn
    real(8)               :: children_fullness
    logical               :: changes

    changes = .false.
    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call patch_count_active(grid(d), p-1)
       end do
    end do
    l = level_start - 1
    do k = 1, grid(d)%lev(l)%length
       do d = 1, size(grid)
          p = grid(d)%lev(l)%elts(k)
          children_fullness = 0
          ! the patch has 4 children including first neighbours 16
          do c = 1, N_CHDRN
             ! allways do one child with 3 adjacent neighbours (2 on side 1 on corner)
             chdrn = get_child_and_neigh_patches(grid(d), p, c)
             do c1 = 1, N_CHDRN
                children_fullness = children_fullness + check_children_fillup(grid(d), chdrn(c1))
             end do
          end do
          children_fullness = children_fullness/16.0_8
          if (children_fullness .gt. FILLUP_THRESHOLD) then
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

  function check_child_required (dom, p, c)
    logical :: check_child_required
    type(Domain) :: dom
    integer      :: p, c

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer :: e, j0, j, i0, i, id, st, en
    logical :: required

    st = -BDRY_THICKNESS
    en =  BDRY_THICKNESS

    call get_offs_Domain(dom, p, offs, dims)

    do j0 = st + 1, PATCH_SIZE/2 + en
       j = j0 - 1 + chd_offs(2,c+1)
       do i0 = st + 1, PATCH_SIZE/2 + en
          i = i0 - 1 + chd_offs(1,c+1)
          id = idx(i, j, offs, dims)
          required = dom%mask_n%elts(id+1) .ge. ADJSPACE
          do e = 1, EDGE
             required = required .or. dom%mask_e%elts(EDGE*id+e) .ge. RESTRCT
          end do
          if (required) then
             check_child_required = .True.
             return
          end if
       end do
    end do
    check_child_required = .False.
  end function check_child_required

  subroutine init_multi_level_mod
    logical :: initialized = .False.

    if (initialized) return ! initialize only once
    call init_comm_mod
    call init_ops_mod
    call init_wavelet_mod
    call init_refine_patch_mod
    initialized = .True.
  end subroutine init_multi_level_mod

  subroutine add_second_level
    integer :: d, c

    do d = 1, size(grid)
       do c = 1, N_CHDRN
          call refine_patch1(grid(d), 1, c-1)
       end do
       do c = 1, N_CHDRN
          call refine_patch2(grid(d), 1, c-1)
       end do
       call connect_children(grid(d), 1)
    end do

    call comm_patch_conn_mpi

    do d = 1, size(grid)
       call update_comm(grid(d))
    end do

    call comm_communication_mpi
    call comm_nodes9_mpi(get_areas, set_areas, NONE)
    call apply_to_penta(area_post_comm, NONE, z_null)
  end subroutine add_second_level

  subroutine fill_up_level
    ! Fills up level level_start+1 and increases level_start
    integer :: d, j, p_par, c, p_chd

    do d = 1, size(grid)
       do j = 1, grid(d)%lev(level_start)%length
          p_par = grid(d)%lev(level_start)%elts(j)
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd .eq. 0) then
                call refine_patch (grid(d), p_par, c - 1)
             end if
          end do
       end do
    end do
    call post_refine
    level_start = level_start+1
  end subroutine fill_up_level

  subroutine fill_up_grid_and_IWT (l)
    ! Fills grid up to level l and does inverse wavelet transform of solution onto grid
    integer :: l
    
    integer :: old_level_start
    
    old_level_start = level_start
    do while (level_start .lt. l)
       call fill_up_level
    end do
    call inverse_wavelet_transform (wav_coeff, sol, old_level_start)
    call update_array_bdry (sol, NONE)
    level_start = old_level_start
  end subroutine fill_up_grid_and_IWT
end module adapt_mod
