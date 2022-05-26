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

  subroutine adapt (set_thresholds, type)
    ! Determines significant wavelets, adaptive grid and all masks associated with adaptive grid
    implicit none
    external           :: set_thresholds
    logical, optional  :: type ! recalculate thresholds
    
    integer :: k, l, d
    logical :: local_type

    ! Recalculate thresholds?
    if (present(type)) then
       local_type = type
    else
       local_type = .true.
    end if

    if (local_type) call set_thresholds

    ! Initialize all nodes and edges to ZERO at finer scales
    do l = level_start+1, level_end
       call apply_onescale__int (set_masks, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, ZERO)
    end do
    
    ! Make nodes and edges with significant wavelet coefficients active
    call update_array_bdry1 (wav_coeff, level_start, level_end, 15)
    call mask_active
    call comm_masks_mpi (NONE)
    
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

    call complete_masks

    ! Set insignificant wavelet coefficients to zero
    if (local_type) call compress_wavelets (wav_coeff)
  end subroutine adapt

  subroutine compress_wavelets (wav)
    ! Sets wavelets associated with inactive grid points to zero
    implicit none
    type(Float_Field), dimension(:,:), target :: wav

    integer :: d, k, l, v
    
    do k = 1, size(wav,2)
       do l = level_start+1, level_end
          do d = 1, size (grid)
             do v = scalars(1), scalars(2)
                wc_s => wav(v,k)%data(d)%elts
                call apply_onescale_d (compress_scalar, grid(d), l, z_null, 0, 1)
                nullify (wc_s)
             end do
             wc_u => wav(S_VELO,k)%data(d)%elts
             call apply_onescale_d (compress_vector, grid(d), l, z_null, 0, 0)
             nullify (wc_u)
          end do
       end do
    end do
    wav%bdry_uptodate = .false.
  end subroutine compress_wavelets

  subroutine compress_wavelets_scalar (wav)
    ! Sets scalar wavelets associated with inactive grid points to zero
    implicit none
    type(Float_Field), dimension(:), target :: wav

    integer :: d, k, l

    do k = 1, size(wav)
       do d = 1, size (grid)
          do l = level_start+1, level_end
             wc_s => wav(k)%data(d)%elts
             call apply_onescale_d (compress_scalar, grid(d), l, z_null, 0, 1)
             nullify (wc_s)
          end do
       end do
    end do
    wav%bdry_uptodate = .false.
  end subroutine compress_wavelets_scalar

  subroutine compress_wavelets_velo (wav)
    ! Sets wavelets associated with inactive grid points to zero
    implicit none
    type(Float_Field), dimension(:), target :: wav

    integer :: d, k, l
    
    do k = 1, size(wav)
       do d = 1, size (grid)
          do l = level_start+1, level_end
             wc_u => wav(k)%data(d)%elts
             call apply_onescale_d (compress_vector, grid(d), l, z_null, 0, 0)
             nullify (wc_u)
          end do
       end do
    end do
    wav%bdry_uptodate = .false.
  end subroutine compress_wavelets_velo

  subroutine compress_scalar (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    
    if (dom%mask_n%elts(id_i) < ADJZONE) wc_s(id_i) = 0d0
  end subroutine compress_scalar

  subroutine compress_vector (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: e, id

    id = idx (i, j, offs, dims)
    
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) < ADJZONE) wc_u(EDGE*id+e) = 0d0
    end do
  end subroutine compress_vector

  subroutine WT_after_step (scaling, wavelet, l_start0)
    !  Everything needed in terms of forward and backward wavelet transform
    !  after one time step (e.g. RK sub-step)
    !    A) compute wavelets and perform backwards transform to conserve mass
    !    B) interpolate values onto adapted grid for next step
    implicit none
    type(Float_Field), dimension(:,:), target :: scaling, wavelet
    integer, optional                         :: l_start0
    
    integer :: d, j, k, l, l_start, v

    if (present(l_start0)) then
       l_start = l_start0
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             velo => scaling(S_VELO,k)%data(d)%elts
             call apply_interscale_d (restrict_velo, grid(d), level_start-1, k, 0, 0)
             nullify (velo)
          end do
       end do
    else
       l_start = level_start
    end if

    call update_array_bdry (scaling, NONE, 16)

    do k = 1, size(scaling,2)
       do l = l_start, level_end-1
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                scalar => scaling(v,k)%data(d)%elts
                wc_s   => wavelet(v,k)%data(d)%elts
                call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
                nullify (scalar, wc_s)
             end do
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d (compute_velo_wavelets, grid(d), l, z_null, 0, 0)
             call apply_to_penta_d (compute_velo_wavelets_penta, grid(d), l, z_null)
             nullify (velo, wc_u)
          end do
          wavelet(:,k)%bdry_uptodate = .false.
       end do
    end do
    call compress_wavelets (wavelet)
    call inverse_wavelet_transform (wavelet, scaling)
  end subroutine WT_after_step

  subroutine WT_after_scalar (scaling, wavelet, l_start0)
    !  Everything needed in terms of forward and backward scalar wavelet transform
    !  after one time step for a vector of scalars
    !    A) compute wavelets and perform backwards transform to conserve mass
    !    B) interpolate values onto adapted grid for next step
    implicit none
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0
    
    integer :: d, j, k, l, l_start

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_vector_bdry (scaling, NONE, 16)

    do k = 1, size (scaling)
       do l = l_start, level_end-1
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             nullify (scalar, wc_s)
          end do
          wavelet(k)%bdry_uptodate = .false.
       end do
    end do
    call compress_wavelets_scalar (wavelet)
    call inverse_scalar_transform (wavelet, scaling)
  end subroutine WT_after_scalar

  subroutine WT_after_velo (scaling, wavelet, l_start0)
    !  Everything needed in terms of forward and backward velocity wavelet transform
    !  after one time step (e.g. RK sub-step)
    !    A) compute wavelets and perform backwards transform to conserve mass
    !    B) interpolate values onto adapted grid for next step
    implicit none
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0
    
    integer :: d, j, k, l, l_start

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
       do k = 1, size(scaling)
          do d = 1, size(grid)
             velo => scaling(k)%data(d)%elts
             call apply_interscale_d (restrict_velo, grid(d), level_start-1, k, 0, 0)
             nullify (velo)
          end do
       end do
    end if

    call update_vector_bdry (scaling, NONE, 16)

    do k = 1, size(scaling)
       do l = l_start, level_end-1
          do d = 1, size(grid)
             velo => scaling(k)%data(d)%elts
             wc_u => wavelet(k)%data(d)%elts
             call apply_interscale_d (compute_velo_wavelets, grid(d), l, z_null, 0, 0)
             call apply_to_penta_d (compute_velo_wavelets_penta, grid(d), l, z_null)
             nullify (velo, wc_u)
          end do
          wavelet(k)%bdry_uptodate = .false.
       end do
    end do
    call compress_wavelets_velo (wavelet)
    call inverse_velo_transform (wavelet, scaling)
  end subroutine WT_after_velo
  
  logical function refine ()
    ! Determines where new patches are needed
    implicit none
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
    refine = sync_max_int (did_refine) == TRUE
    return
  end function refine

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
  
  subroutine fill_up_grid_and_IWT (l)
    ! Fills grid up to level l and does inverse wavelet transform of solution onto grid
    implicit none
    integer :: l
    
    integer :: old_level_start
    
    old_level_start = level_start
    do while (level_start < l)
       if (rank == 0) write(6,'(a,i2)') 'Filling up level ', level_start+1
       call fill_up_level
    end do
    call inverse_wavelet_transform (wav_coeff, sol, old_level_start)
    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE, 17)
  end subroutine fill_up_grid_and_IWT

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

  !!! Remaining routines were added by Matthias to prepare for allowing domain splitting for better rebalancing !!!
  !!! they are not yet incorporated in the code !!!
  logical function remove_inside_patches ()
    ! Removes patches that are not required because they are far enough away from the locally finest level
    implicit none
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

  real(8) function check_children_fillup (dom, p_par)
    implicit none
    type(Domain) :: dom
    integer      :: p_par
    
    integer               :: active, p_chd, c, i
    integer, dimension(4) :: cn

    ! Do not fill up children + remove if any child or neighbour patch is missing
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
       ! Count real child double since neighbours are counted double during accumulation
       active = active + 2*dom%patch%elts(cn(1)+1)%active
    end do
    check_children_fillup = dble(active)/dble(N_CHDRN*5*DOF_PER_PATCH)
  end function check_children_fillup

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
end module adapt_mod
