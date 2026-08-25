module mask_mod
  ! Module containing routines that define masks on adaptive grid.
  ! (required by adapt_mod and refine_patch_mod)

  use iso_fortran_env, only : error_unit, int64

  use kind_mod,       only : dp
  use arch_mod,       only : rank
  use comm_mpi_mod,   only : comm_masks_mpi, update_bdry1
  use dyn_arrays,     only : extend, init, Float_Array, Int_Array
  use domain_mod,     only : chd_offs, Domain, Float_Field, &
       get_offs_Domain, grid, &
       wav_coeff, idx, id_edge
  use domain_ops_mod, only : apply_bdry, apply_interscale, &
       apply_onescale, apply_onescale_d, apply_onescale__int
  use patch_mod,      only : PATCH_SIZE
  
  use shared_mod, only : ADJSPACE, ADJZONE, BDRY_THICKNESS, EDGE, &
       N_BDRY, N_CHDRN, TOLRNZ, RT, DG, UP, TRSK, RESTRCT, ZERO, z_null, &
       min_level, Laplace_rotu, Laplace_sclr, level_fill, S_VELO, FROZEN, NONE, level_start, level_end, max_level, &
       zlevels, scalars, threshold

  
  private
  public :: init_masks_zero, mask_active, mask_adj_child, mask_adj_same_scale, mask_restrict_same_scale, mask_adj_finer_scale
  public :: init_masks, mask_second_neighbours, complete_masks, mask_trsk
  public :: begin_pre_refinement_mask_shadow
  public :: advance_pre_refinement_mask_shadow
  public :: begin_post_refinement_mask_shadow
  public :: advance_post_refinement_mask_shadow
  public :: begin_wavelet_compression_shadow
  public :: compare_wavelet_compression_shadow

  type(Int_Array), allocatable, save :: pre_refinement_mask_n(:)
  type(Int_Array), allocatable, save :: pre_refinement_mask_e(:)
  logical, allocatable, save :: pre_refinement_request(:)
  logical, save :: pre_refinement_shadow_ready = .false.
  integer, save :: pre_refinement_shadow_stage = -1
  integer(int64), save :: pre_refinement_shadow_allocations = 0_int64
  integer(int64), save :: pre_refinement_allocation_checkpoint = 0_int64

  type(Float_Array), allocatable, target, save :: &
       wavelet_compression_expected(:,:,:)
  real(dp), pointer, save :: wavelet_compression_values(:) => null()
  logical, save :: wavelet_compression_shadow_ready = .false.
  integer(int64), save :: wavelet_compression_shadow_allocations = 0_int64
  integer(int64), save :: wavelet_compression_allocation_checkpoint = &
       0_int64

  
contains


  subroutine begin_pre_refinement_mask_shadow
    ! Retain the exact pre-adaptation mask state in persistent storage.  The
    ! staged masks are advanced independently while the authoritative Domain
    ! masks follow the normal adapt path.

    implicit none

    integer :: d
    integer :: request_extent

    call ensure_pre_refinement_mask_shadow_storage
    do d = 1,size(grid)
       pre_refinement_mask_n(d)%length = grid(d)%mask_n%length
       pre_refinement_mask_e(d)%length = grid(d)%mask_e%length
       pre_refinement_mask_n(d)%elts = grid(d)%mask_n%elts
       pre_refinement_mask_e(d)%elts = grid(d)%mask_e%elts
    end do
    request_extent = pre_refinement_request_extent()
    if (.not. allocated(pre_refinement_request)) then
       allocate(pre_refinement_request(max(1,request_extent)))
       pre_refinement_shadow_allocations = &
            pre_refinement_shadow_allocations+1_int64
    else if (size(pre_refinement_request) < request_extent) then
       deallocate(pre_refinement_request)
       allocate(pre_refinement_request(max(1,request_extent)))
       pre_refinement_shadow_allocations = &
            pre_refinement_shadow_allocations+1_int64
    end if
    pre_refinement_request = .false.
    pre_refinement_shadow_stage = 0
    pre_refinement_shadow_ready = .true.
    pre_refinement_allocation_checkpoint = &
         pre_refinement_shadow_allocations

  end subroutine begin_pre_refinement_mask_shadow


  subroutine advance_pre_refinement_mask_shadow(stage)
    ! Replay one legacy pre-refinement phase on the staged masks, restore the
    ! authoritative masks, then require exact equality.  Stage four also
    ! compares the complete parent/child refinement-requirement manifest.

    implicit none

    integer, intent(in) :: stage

    if (.not. pre_refinement_shadow_ready .or. &
         stage /= pre_refinement_shadow_stage+1 .or. &
         stage < 1 .or. stage > 4) then
       error stop "advance_pre_refinement_mask_shadow: invalid stage"
    end if
    call swap_pre_refinement_masks
    select case (stage)
    case (1)
       call init_masks_zero
    case (2)
       call mask_active
    case (3)
       call mask_adj_same_scale
    case (4)
       call mask_restrict_same_scale
    end select
    call swap_pre_refinement_masks
    call compare_pre_refinement_masks(stage)
    pre_refinement_shadow_stage = stage

    if (pre_refinement_shadow_allocations /= &
         pre_refinement_allocation_checkpoint) then
       error stop &
            "advance_pre_refinement_mask_shadow: persistent storage changed"
    end if
    if (stage == 4) then
       call compare_pre_refinement_requests
       pre_refinement_shadow_ready = .false.
       pre_refinement_shadow_stage = -1
    end if

  end subroutine advance_pre_refinement_mask_shadow


  subroutine begin_post_refinement_mask_shadow
    ! Restart the persistent mask shadow after refine/post_refine has made any
    ! required topology changes.  The remaining mask phases do not alter
    ! topology and can therefore be replayed independently and exactly.

    implicit none

    integer :: d

    if (pre_refinement_shadow_ready) then
       error stop "begin_post_refinement_mask_shadow: shadow is still active"
    end if
    call ensure_pre_refinement_mask_shadow_storage
    do d = 1,size(grid)
       pre_refinement_mask_n(d)%length = grid(d)%mask_n%length
       pre_refinement_mask_e(d)%length = grid(d)%mask_e%length
       pre_refinement_mask_n(d)%elts = grid(d)%mask_n%elts
       pre_refinement_mask_e(d)%elts = grid(d)%mask_e%elts
    end do
    pre_refinement_shadow_stage = 5
    pre_refinement_shadow_ready = .true.
    pre_refinement_allocation_checkpoint = &
         pre_refinement_shadow_allocations

  end subroutine begin_post_refinement_mask_shadow


  subroutine advance_post_refinement_mask_shadow(stage)
    ! Replay one post-refinement mask phase on the staged masks, restore the
    ! authoritative masks, and require exact equality.

    implicit none

    integer, intent(in) :: stage

    if (.not. pre_refinement_shadow_ready .or. &
         stage /= pre_refinement_shadow_stage+1 .or. &
         stage < 6 .or. stage > 8) then
       error stop "advance_post_refinement_mask_shadow: invalid stage"
    end if
    call swap_pre_refinement_masks
    select case (stage)
    case (6)
       call mask_adj_finer_scale
    case (7)
       call complete_masks
    case (8)
       call mask_second_neighbours
    end select
    call swap_pre_refinement_masks
    call compare_pre_refinement_masks(stage)
    pre_refinement_shadow_stage = stage

    if (pre_refinement_shadow_allocations /= &
         pre_refinement_allocation_checkpoint) then
       error stop &
            "advance_post_refinement_mask_shadow: persistent storage changed"
    end if
    if (stage == 8) then
       pre_refinement_shadow_ready = .false.
       pre_refinement_shadow_stage = -1
    end if

  end subroutine advance_post_refinement_mask_shadow


  subroutine begin_wavelet_compression_shadow(wavelet)
    ! Retain the boundary-refreshed, uncompressed wavelet family and apply
    ! the completed masks independently to persistent expected storage.

    implicit none

    type(Float_Field), intent(in), target :: wavelet(:,:)

    integer :: d
    integer :: k
    integer :: l
    integer :: v

    if (wavelet_compression_shadow_ready) then
       error stop "begin_wavelet_compression_shadow: shadow is still active"
    end if
    if (pre_refinement_shadow_ready .or. &
         pre_refinement_shadow_stage /= -1) then
       error stop &
            "begin_wavelet_compression_shadow: mask shadow is incomplete"
    end if
    if (associated(wavelet_compression_values)) then
       error stop &
            "begin_wavelet_compression_shadow: work pointer is associated"
    end if
    if (level_start < level_end .and. &
         .not. all(wavelet%bdry_uptodate)) then
       error stop &
            "begin_wavelet_compression_shadow: boundary state is incomplete"
    end if
    call ensure_wavelet_compression_shadow_storage(wavelet)
    do k = 1,size(wavelet,2)
       do v = 1,size(wavelet,1)
          do d = 1,size(grid)
             wavelet_compression_expected(v,k,d)%length = &
                  wavelet(v,k)%data(d)%length
             wavelet_compression_expected(v,k,d)%elts = &
                  wavelet(v,k)%data(d)%elts
          end do
       end do
    end do

    do k = 1,size(wavelet,2)
       do l = level_start+1,level_end
          do d = 1,size(grid)
             do v = scalars(1),scalars(2)
                wavelet_compression_values => &
                     wavelet_compression_expected(v,k,d)%elts
                call apply_onescale_d( &
                     compress_shadow_scalar,grid(d),l,z_null,0,1)
                nullify(wavelet_compression_values)
             end do
             wavelet_compression_values => &
                  wavelet_compression_expected(S_VELO,k,d)%elts
             call apply_onescale_d( &
                  compress_shadow_vector,grid(d),l,z_null,0,0)
             nullify(wavelet_compression_values)
          end do
       end do
    end do

    wavelet_compression_shadow_ready = .true.
    wavelet_compression_allocation_checkpoint = &
         wavelet_compression_shadow_allocations

  end subroutine begin_wavelet_compression_shadow


  subroutine compare_wavelet_compression_shadow(wavelet)
    ! Require the production compression result to match the independently
    ! masked expected family bit-for-bit, including retained coefficients.

    implicit none

    type(Float_Field), intent(in) :: wavelet(:,:)

    integer :: d
    integer :: i
    integer :: k
    integer :: v
    integer(int64) :: actual_bits
    integer(int64) :: expected_bits

    if (.not. wavelet_compression_shadow_ready) then
       error stop "compare_wavelet_compression_shadow: shadow is unavailable"
    end if
    if (wavelet_compression_shadow_allocations /= &
         wavelet_compression_allocation_checkpoint) then
       error stop &
            "compare_wavelet_compression_shadow: persistent storage changed"
    end if
    if (any(wavelet%bdry_uptodate)) then
       error stop &
            "compare_wavelet_compression_shadow: boundary state is current"
    end if
    do k = 1,size(wavelet,2)
       do v = 1,size(wavelet,1)
          do d = 1,size(grid)
             if (wavelet(v,k)%data(d)%length /= &
                  wavelet_compression_expected(v,k,d)%length) then
                error stop &
                     "compare_wavelet_compression_shadow: extent differs"
             end if
             do i = 1,wavelet(v,k)%data(d)%length
                actual_bits = transfer( &
                     wavelet(v,k)%data(d)%elts(i),0_int64)
                expected_bits = transfer( &
                     wavelet_compression_expected(v,k,d)%elts(i), &
                     0_int64)
                if (actual_bits == expected_bits) cycle
                write(error_unit,'(/,a,i0,a)') &
                     "Rank ",rank,": compressed wavelet differs"
                write(error_unit,'(a,i0,a,i0,a,i0,a,i0)') &
                     "  variable = ",v,", level slot = ",k, &
                     ", Domain = ",d,", index = ",i
                write(error_unit,'(a,2(es24.16,1x))') &
                     "  actual, expected = ", &
                     wavelet(v,k)%data(d)%elts(i), &
                     wavelet_compression_expected(v,k,d)%elts(i)
                flush(error_unit)
                error stop "wavelet-compression comparison failed"
             end do
          end do
       end do
    end do
    wavelet_compression_shadow_ready = .false.

  end subroutine compare_wavelet_compression_shadow


  subroutine ensure_wavelet_compression_shadow_storage(wavelet)

    implicit none

    type(Float_Field), intent(in) :: wavelet(:,:)

    integer :: d
    integer :: k
    integer :: v
    logical :: shape_differs

    shape_differs = .false.
    if (allocated(wavelet_compression_expected)) then
       shape_differs = &
            size(wavelet_compression_expected,1) /= size(wavelet,1) .or. &
            size(wavelet_compression_expected,2) /= size(wavelet,2) .or. &
            size(wavelet_compression_expected,3) /= size(grid)
       if (shape_differs) then
          do d = 1,size(wavelet_compression_expected,3)
             do k = 1,size(wavelet_compression_expected,2)
                do v = 1,size(wavelet_compression_expected,1)
                   if (allocated( &
                        wavelet_compression_expected(v,k,d)%elts)) &
                        deallocate( &
                        wavelet_compression_expected(v,k,d)%elts)
                end do
             end do
          end do
          deallocate(wavelet_compression_expected)
       end if
    end if
    if (.not. allocated(wavelet_compression_expected)) then
       allocate(wavelet_compression_expected( &
            size(wavelet,1),size(wavelet,2),size(grid)))
       wavelet_compression_shadow_allocations = &
            wavelet_compression_shadow_allocations+1_int64
    end if

    do k = 1,size(wavelet,2)
       do v = 1,size(wavelet,1)
          if (.not. allocated(wavelet(v,k)%data)) then
             error stop &
                  "ensure_wavelet_compression_shadow_storage: field absent"
          end if
          if (size(wavelet(v,k)%data) /= size(grid)) then
             error stop &
                  "ensure_wavelet_compression_shadow_storage: layout differs"
          end if
          do d = 1,size(grid)
             if (.not. allocated(wavelet(v,k)%data(d)%elts)) then
                error stop &
                     "ensure_wavelet_compression_shadow_storage: data absent"
             end if
             if (allocated( &
                  wavelet_compression_expected(v,k,d)%elts)) then
                if (size( &
                     wavelet_compression_expected(v,k,d)%elts) /= &
                     size(wavelet(v,k)%data(d)%elts)) then
                   deallocate( &
                        wavelet_compression_expected(v,k,d)%elts)
                end if
             end if
             if (.not. allocated( &
                  wavelet_compression_expected(v,k,d)%elts)) then
                allocate(wavelet_compression_expected(v,k,d)%elts( &
                     size(wavelet(v,k)%data(d)%elts)))
                wavelet_compression_shadow_allocations = &
                     wavelet_compression_shadow_allocations+1_int64
             end if
          end do
       end do
    end do

  end subroutine ensure_wavelet_compression_shadow_storage


  subroutine compress_shadow_scalar(dom,i,j,zlev,offs,dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: zlev
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: id_i

    id_i = idx(i,j,offs,dims)+1
    if (dom%mask_n%elts(id_i) < ADJZONE) &
         wavelet_compression_values(id_i) = 0.0_dp

  end subroutine compress_shadow_scalar


  subroutine compress_shadow_vector(dom,i,j,zlev,offs,dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i
    integer, intent(in) :: j
    integer, intent(in) :: zlev
    integer, intent(in) :: offs(N_BDRY+1)
    integer, intent(in) :: dims(2,N_BDRY+1)

    integer :: e
    integer :: id
    integer :: id_e

    id = idx(i,j,offs,dims)
    do e = 1,EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) < ADJZONE) &
            wavelet_compression_values(id_e) = 0.0_dp
    end do

  end subroutine compress_shadow_vector


  subroutine ensure_pre_refinement_mask_shadow_storage

    implicit none

    integer :: d

    if (allocated(pre_refinement_mask_n) .neqv. &
         allocated(pre_refinement_mask_e)) then
       error stop &
            "ensure_pre_refinement_mask_shadow_storage: inconsistent storage"
    end if
    if (allocated(pre_refinement_mask_n)) then
       if (size(pre_refinement_mask_n) /= size(grid) .or. &
            size(pre_refinement_mask_e) /= size(grid)) then
          do d = 1,size(pre_refinement_mask_n)
             if (allocated(pre_refinement_mask_n(d)%elts)) &
                  deallocate(pre_refinement_mask_n(d)%elts)
             if (allocated(pre_refinement_mask_e(d)%elts)) &
                  deallocate(pre_refinement_mask_e(d)%elts)
          end do
          deallocate(pre_refinement_mask_n,pre_refinement_mask_e)
       end if
    end if
    if (.not. allocated(pre_refinement_mask_n)) then
       allocate(pre_refinement_mask_n(size(grid)))
       allocate(pre_refinement_mask_e(size(grid)))
       pre_refinement_shadow_allocations = &
            pre_refinement_shadow_allocations+2_int64
    end if

    do d = 1,size(grid)
       if (allocated(pre_refinement_mask_n(d)%elts)) then
          if (size(pre_refinement_mask_n(d)%elts) /= &
               size(grid(d)%mask_n%elts)) then
             deallocate(pre_refinement_mask_n(d)%elts)
          end if
       end if
       if (.not. allocated(pre_refinement_mask_n(d)%elts)) then
          allocate(pre_refinement_mask_n(d)%elts( &
               size(grid(d)%mask_n%elts)))
          pre_refinement_shadow_allocations = &
               pre_refinement_shadow_allocations+1_int64
       end if

       if (allocated(pre_refinement_mask_e(d)%elts)) then
          if (size(pre_refinement_mask_e(d)%elts) /= &
               size(grid(d)%mask_e%elts)) then
             deallocate(pre_refinement_mask_e(d)%elts)
          end if
       end if
       if (.not. allocated(pre_refinement_mask_e(d)%elts)) then
          allocate(pre_refinement_mask_e(d)%elts( &
               size(grid(d)%mask_e%elts)))
          pre_refinement_shadow_allocations = &
               pre_refinement_shadow_allocations+1_int64
       end if
    end do

  end subroutine ensure_pre_refinement_mask_shadow_storage


  subroutine swap_pre_refinement_masks

    implicit none

    integer :: d
    integer, allocatable :: temporary(:)

    if (.not. allocated(pre_refinement_mask_n)) then
       error stop "swap_pre_refinement_masks: storage is unavailable"
    end if
    if (.not. allocated(pre_refinement_mask_e)) then
       error stop "swap_pre_refinement_masks: storage is unavailable"
    end if
    if (size(pre_refinement_mask_n) /= size(grid) .or. &
         size(pre_refinement_mask_e) /= size(grid)) then
       error stop "swap_pre_refinement_masks: storage is unavailable"
    end if
    do d = 1,size(grid)
       call move_alloc(grid(d)%mask_n%elts,temporary)
       call move_alloc(pre_refinement_mask_n(d)%elts,grid(d)%mask_n%elts)
       call move_alloc(temporary,pre_refinement_mask_n(d)%elts)

       call move_alloc(grid(d)%mask_e%elts,temporary)
       call move_alloc(pre_refinement_mask_e(d)%elts,grid(d)%mask_e%elts)
       call move_alloc(temporary,pre_refinement_mask_e(d)%elts)
    end do

  end subroutine swap_pre_refinement_masks


  subroutine compare_pre_refinement_masks(stage)

    implicit none

    integer, intent(in) :: stage

    integer :: d
    integer :: mismatch

    do d = 1,size(grid)
       mismatch = find_first_mask_mismatch( &
            grid(d)%mask_n%elts,pre_refinement_mask_n(d)%elts, &
            grid(d)%mask_n%length)
       if (mismatch > 0) then
          write(error_unit,'(/,a,i0,a)') &
               "Rank ",rank,": adaptation node mask differs"
          write(error_unit,'(a,i0,a,i0,a,i0,a,i0)') &
               "  stage = ",stage,", Domain = ",d,", index = ", &
               mismatch,", authoritative = ", &
               grid(d)%mask_n%elts(mismatch)
          write(error_unit,'(a,i0)') "  staged = ", &
               pre_refinement_mask_n(d)%elts(mismatch)
          flush(error_unit)
          error stop "adaptation node-mask comparison failed"
       end if
       mismatch = find_first_mask_mismatch( &
            grid(d)%mask_e%elts,pre_refinement_mask_e(d)%elts, &
            grid(d)%mask_e%length)
       if (mismatch > 0) then
          write(error_unit,'(/,a,i0,a)') &
               "Rank ",rank,": adaptation edge mask differs"
          write(error_unit,'(a,i0,a,i0,a,i0,a,i0)') &
               "  stage = ",stage,", Domain = ",d,", index = ", &
               mismatch,", authoritative = ", &
               grid(d)%mask_e%elts(mismatch)
          write(error_unit,'(a,i0)') "  staged = ", &
               pre_refinement_mask_e(d)%elts(mismatch)
          flush(error_unit)
          error stop "adaptation edge-mask comparison failed"
       end if
    end do

  end subroutine compare_pre_refinement_masks


  integer function find_first_mask_mismatch ( &
       authoritative,staged,extent) result(first)

    implicit none

    integer, intent(in) :: authoritative(:)
    integer, intent(in) :: staged(:)
    integer, intent(in) :: extent
    integer :: i

    first = 0
    if (extent < 0 .or. extent > size(authoritative) .or. &
         extent > size(staged)) then
       error stop "find_first_mask_mismatch: invalid extent"
    end if
    do i = 1,extent
       if (authoritative(i) == staged(i)) cycle
       first = i
       return
    end do

  end function find_first_mask_mismatch


  subroutine compare_pre_refinement_requests

    implicit none

    integer :: c
    integer :: d
    integer :: p_storage
    integer :: position
    integer :: request_extent
    logical :: staged_request

    request_extent = pre_refinement_request_extent()
    if (.not. allocated(pre_refinement_request)) then
       error stop "compare_pre_refinement_requests: storage is unavailable"
    end if
    if (size(pre_refinement_request) < request_extent) then
       error stop "compare_pre_refinement_requests: storage is unavailable"
    end if

    position = 0
    do d = 1,size(grid)
       do p_storage = 3,grid(d)%patch%length
          do c = 0,N_CHDRN-1
             position = position+1
             pre_refinement_request(position) = &
                  child_request_required(grid(d),p_storage-1,c)
          end do
       end do
    end do
    if (position /= request_extent) then
       error stop "compare_pre_refinement_requests: extent differs"
    end if

    call swap_pre_refinement_masks
    position = 0
    do d = 1,size(grid)
       do p_storage = 3,grid(d)%patch%length
          do c = 0,N_CHDRN-1
             position = position+1
             staged_request = &
                  child_request_required(grid(d),p_storage-1,c)
             if (staged_request .neqv. &
                  pre_refinement_request(position)) then
                write(error_unit,'(/,a,i0,a)') &
                     "Rank ",rank, &
                     ": pre-refinement child request differs"
                write(error_unit,'(a,i0,a,i0,a,i0)') &
                     "  Domain = ",d,", parent patch = ", &
                     p_storage-1,", child = ",c
                flush(error_unit)
                error stop &
                     "pre-refinement request-manifest comparison failed"
             end if
          end do
       end do
    end do
    call swap_pre_refinement_masks

  end subroutine compare_pre_refinement_requests


  integer function pre_refinement_request_extent() result(extent)

    implicit none

    integer :: d

    extent = 0
    do d = 1,size(grid)
       extent = extent + N_CHDRN*max(0,grid(d)%patch%length-2)
    end do

  end function pre_refinement_request_extent


  logical function child_request_required(dom,p_parent,child) &
       result(requested)

    implicit none

    type(Domain), intent(in) :: dom
    integer, intent(in) :: p_parent
    integer, intent(in) :: child

    integer :: dims_parent(2,N_BDRY+1)
    integer :: i0
    integer :: i_parent
    integer :: id_parent
    integer :: j0
    integer :: j_parent
    integer :: offs_parent(N_BDRY+1)
    logical :: is_required

    requested = .false.
    if (p_parent < 0 .or. p_parent >= dom%patch%length .or. &
         child < 0 .or. child >= N_CHDRN) then
       error stop "child_request_required: invalid address"
    end if
    call get_offs_Domain(dom,p_parent,offs_parent,dims_parent)
    do j0 = -BDRY_THICKNESS+1, &
         PATCH_SIZE/2+BDRY_THICKNESS
       j_parent = j0-1+chd_offs(2,child+1)
       do i0 = -BDRY_THICKNESS+1, &
            PATCH_SIZE/2+BDRY_THICKNESS
          i_parent = i0-1+chd_offs(1,child+1)
          id_parent = idx( &
               i_parent,j_parent,offs_parent,dims_parent)
          is_required = dom%mask_n%elts(id_parent+1) >= ADJSPACE .or. &
               maxval(dom%mask_e%elts(id_edge(id_parent))) >= RESTRCT
          if (is_required) then
             requested = .true.
             return
          end if
       end do
    end do

  end function child_request_required


  subroutine init_masks_zero
    ! Initialize all node and edge masks to ZERO at finer scales
    implicit none
    integer :: l

    do l = level_start+1, level_end
       call apply_onescale__int (set_masks, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS, ZERO)
    end do
  end subroutine init_masks_zero

  
  subroutine mask_active
    ! Determine active mask
    implicit none
    integer :: l

    call update_bdry1 (wav_coeff, level_start, level_end, 910)

    ! Set active grid at finest scale
    call apply_onescale (mask_tol_vars, level_end, z_null, 0, 1)
    call comm_masks_mpi (level_end)

    ! Set active grid at coarser scales
    do l = level_end-1, level_start, -1
       call apply_onescale (mask_tol_vars, l, z_null, -1, 2)

       ! Add  parents to active mask at scale l if any of 6 child neighbours at scale l+1 are in active mask
       call apply_interscale (mask_parent_nodes, l, z_null, 0, 1) ! also modifies child mask
       call apply_interscale (mask_parent_edges, l, z_null, 0, 0)
       call comm_masks_mpi (l)
    end do
    call comm_masks_mpi (NONE)
  end subroutine mask_active

  
  subroutine mask_adj_same_scale
    ! Add nearest neighbour wavelets of active nodes and edges at same scale
    implicit none
    integer :: l

    do l = level_start, level_end
       call apply_onescale (mask_adj_same_scale_nodes_edges, l, z_null, 0, 1)
    end do
    call comm_masks_mpi (NONE)
  end subroutine mask_adj_same_scale

  
  subroutine mask_restrict_same_scale
    ! Needed if bdry is only 2 layers for scenario:
    ! scalar wavelet coefficient > threshold @ PATCH_SIZE + 2 => flux restr @ PATCH_SIZE + 1
    ! => patch needed (contains flux for corrective part of R_F)
    implicit none
    integer :: l

    do l = level_start, min (level_end, max_level-1)
       call apply_onescale (mask_restrict_flux, l, z_null, 0, 0)
    end do
    call comm_masks_mpi (NONE)
  end subroutine mask_restrict_same_scale

  
  subroutine mask_adj_finer_scale
    ! Add adjacent mask at finer scale
    implicit none
    integer :: l

    do l = level_end-1, level_start, -1
       call apply_interscale (mask_adj_child,           l,   z_null, 0, 1)
       call apply_onescale   (mask_edges_if_both_nodes, l+1, z_null, 0, 0)
    end do
    call comm_masks_mpi (NONE)
  end subroutine mask_adj_finer_scale

  
  subroutine complete_masks
    ! Ensure consistency between node and edge masks
    implicit none
    integer :: l

    do l = level_end-1, level_start+1, -1
       call apply_interscale (mask_edges_consist,  l, z_null, 0, 1) ! adds child/parent edge to adjacent mask if any of neighbour edges in adjacent mask
       call comm_masks_mpi (l+1)
       call apply_interscale (mask_edges_consist2, l, z_null, 0, 0) ! adds child/parent edge to adjacent mask if any of second neighbour edges in adjacent mask
       call comm_masks_mpi (l)
    end do

    if (level_start < level_end) then
       call apply_interscale (mask_edges_consist, level_start, z_null, 0, 1) 
       call comm_masks_mpi (level_start)  
       call comm_masks_mpi (level_start+1)
    end if

    do l = level_end-1, level_start+1, -1
       call apply_onescale   (mask_nodes_if_all_edges, l+1, z_null, 0, 1) ! add node to adjacent mask if all neighbour edges in adjacent mask
       call apply_interscale (prolong_node_adjzone,    l,   z_null, 0, 1) ! add parent to adjacent mask if child is in adjacent mask
    end do
    if (level_start+1 <= level_end) call apply_onescale (mask_nodes_if_all_edges, level_start+1, z_null, 0, 1)
    call comm_masks_mpi (NONE) 
  end subroutine complete_masks

  
  subroutine mask_tol_vars (dom, i, j, zlev, offs, dims)
    ! Add nodes/edges to active mask
    ! (do not adapt on soil layers zmin<=k<=0)
    use utils_mod
    implicit none
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i, j, zlev
    integer, dimension(N_BDRY+1), intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims

    integer  :: d, e, id, id_e, id_i, k, l, v
    real(dp) :: wc
    logical  :: active

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    l    = dom%level%elts(id_i)

    if (dom%mask_n%elts(id_i) == FROZEN) return

    ! Scalars
    active = .false.
    do k = 1, zlevels
       do v = scalars(1), scalars(2)
          wc = wav_coeff(v,k)%data(d)%elts(id_i)
          
          if (abs (wc) >= threshold(v,k) .or. l < level_fill) active = .true.
       end do
    end do

    if (active) then
       dom%mask_n%elts(id_i) = TOLRNZ
    else
       if (dom%mask_n%elts(id_i) > ADJZONE) dom%mask_n%elts(id_i) = ADJZONE
    end if

    ! Vectors
    do e = 1, EDGE
       id_e = EDGE*id + e

       active = .false.
       do k = 1, zlevels
          wc = wav_coeff(S_VELO,k)%data(d)%elts(id_e)
          
          if (abs (wc) >= threshold(S_VELO,k) .or. l < level_fill) active = .true.
       end do

       if (active) then
          dom%mask_e%elts(id_e) = TOLRNZ
       else
          if (dom%mask_e%elts(id_e) > ADJZONE) dom%mask_e%elts(id_e) = ADJZONE
       end if
    end do
  end subroutine mask_tol_vars

  
  subroutine mask_parent_nodes (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent node to active mask if any of its 6 child neighbours is in active mas
    ! (also ensures child is added to active mask if its direct parent is made active)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_par, id_chd, idN_chd, idE_chd, idNE_chd, idSW_chd, idS_chd, idW_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd  = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idSW_chd = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS_chd  = idx (i_chd,   j_chd-1, offs_chd, dims_chd)

    if ( dom%mask_n%elts(idE_chd +1) == TOLRNZ .or. &
         dom%mask_n%elts(idNE_chd+1) == TOLRNZ .or. &
         dom%mask_n%elts(idN_chd +1) == TOLRNZ .or. &
         dom%mask_n%elts(idW_chd +1) == TOLRNZ .or. &
         dom%mask_n%elts(idSW_chd+1) == TOLRNZ .or. &
         dom%mask_n%elts(idS_chd +1) == TOLRNZ ) then

       call set_at_least (dom%mask_n%elts(id_par+1), TOLRNZ)
       call set_at_least (dom%mask_n%elts(id_chd+1), TOLRNZ) ! ensure child is in active zone if parent is
    end if
  end subroutine mask_parent_nodes

  
  subroutine mask_parent_edges (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent edge to active mask if any of its 6 child neighbours is in active mask
     
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_chd, id_par, idE_chd, idN_chd, idNE_chd, idNW_chd, idS_chd, idSE_chd, idW_chd

    id_par = idx (i_par, j_par, offs_par, dims_par) 
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) 

    ! Three neighbours of child node
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd  = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idS_chd  = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idNW_chd = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSE_chd = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    ! Check six child edges neighbouring each parent edge to see if at least one is active
    if ( dom%mask_e%elts(EDGE*id_chd  +RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_chd +RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_chd +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_chd +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idSE_chd+UP+1) == TOLRNZ ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), TOLRNZ)

    if ( dom%mask_e%elts(EDGE*idN_chd +RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE_chd+RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_chd  +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE_chd+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNE_chd+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) == TOLRNZ ) &
         call set_at_least(dom%mask_e%elts(EDGE*id_par+DG+1), TOLRNZ)

    if ( dom%mask_e%elts(EDGE*idN_chd +RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idNW_chd+RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_chd +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_chd +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_chd  +UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_chd +UP+1) == TOLRNZ ) &
         call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), TOLRNZ)
  end subroutine mask_parent_edges

  
  subroutine mask_adj_same_scale_nodes_edges (dom, i, j, zlev, offs, dims)
    ! Add nearest neighbours of active nodes/edges at same scale to adjacent ask
    ! (at least one of 6 neighbouring nodes/edges must be active)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idE, idN, idNE, idS, idW, idSW

    id  = idx (i, j, offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    ! Add node to adjacent mask if at least one of 6 neighbouring nodes is in active mask
    if ( dom%mask_n%elts(idE +1) == TOLRNZ .or. &
         dom%mask_n%elts(idNE+1) == TOLRNZ .or. &
         dom%mask_n%elts(idN +1) == TOLRNZ .or. &
         dom%mask_n%elts(idW +1) == TOLRNZ .or. &
         dom%mask_n%elts(idSW+1) == TOLRNZ .or. &
         dom%mask_n%elts(idS +1) == TOLRNZ ) &
         call set_at_least (dom%mask_n%elts(id+1), ADJSPACE)
    
    ! Add edges to adjacent mask if at least one of 4 neighouring edges is in active mask
    if ( dom%mask_e%elts(EDGE*id +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), ADJSPACE)

    if ( dom%mask_e%elts(EDGE*id +UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id +RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE+UP+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), ADJSPACE)

    if ( dom%mask_e%elts(EDGE*id +DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+RT+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW+DG+1) == TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN+RT+1) == TOLRNZ ) &
         call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), ADJSPACE)
  end subroutine mask_adj_same_scale_nodes_edges

  
  subroutine mask_restrict_flux (dom, i_par, j_par, zlev, offs_par, dims_par)
    ! Add edges required for flux restriction
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1)
   
    integer :: id_par, idE_par, idN_par, idNE_par

    id_par = idx (i_par, j_par, offs_par, dims_par)
    
    idE_par  = idx (i_par+1, j_par,   offs_par, dims_par)
    idNE_par = idx (i_par+1, j_par+1, offs_par, dims_par)
    idN_par  = idx (i_par,   j_par+1, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) >= ADJSPACE) then
       call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
    else
       if (dom%mask_n%elts(idE_par +1) >= ADJSPACE) call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), RESTRCT)
       if (dom%mask_n%elts(idNE_par+1) >= ADJSPACE) call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), RESTRCT)
       if (dom%mask_n%elts(idN_par +1) >= ADJSPACE) call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), RESTRCT)
    end if
  end subroutine mask_restrict_flux

  
  subroutine mask_adj_child (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add nearest node/edge child neighbours at finer scale to adjacent mask if parent is active
    ! (includes some edge->node and node->edge cross masking)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_par, idS_par, idW_par, idSW_par, idE_par, idN_par, idNE_par
    integer :: id_chd, idE_chd, idNE_chd, idN2E_chd, id2NE_chd, idN_chd, idW_chd, idNW_chd
    integer :: idS2W_chd, idSW_chd, idS_chd, id2SW_chd, idSE_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)
    
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW_chd   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idSW_chd  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS_chd   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)

    idNW_chd  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSE_chd  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)
    
    idN2E_chd = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    id2SW_chd = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idS2W_chd = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)

    idE_par  = idx (i_par+1, j_par,   offs_par, dims_par)
    idNE_par = idx (i_par+1, j_par+1, offs_par, dims_par)
    idN_par  = idx (i_par,   j_par+1, offs_par, dims_par)
    idW_par  = idx (i_par-1, j_par,   offs_par, dims_par)
    idSW_par = idx (i_par-1, j_par-1, offs_par, dims_par)
    idS_par  = idx (i_par,   j_par-1, offs_par, dims_par)
  
    if (dom%mask_n%elts(id_par+1) >= TOLRNZ) then
       ! Nearest neighbour nodes in same cell at finer scale
       call set_at_least (dom%mask_n%elts(id_chd  +1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idE_chd+ 1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idNE_chd+1),  ADJZONE)
       call set_at_least (dom%mask_n%elts(idN_chd +1),  ADJZONE)

       ! Needed for prolongation of scalars
       call set_at_least (dom%mask_n%elts(idW_chd  +1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idSW_chd +1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idS_chd  +1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idNW_chd +1), ADJZONE)
       call set_at_least (dom%mask_n%elts(idSE_chd +1), ADJZONE)

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
       call set_at_least (dom%mask_n%elts(idN_chd +1), ADJZONE) ! add associated chid node
       if (dom%mask_e%elts(EDGE*idS_par +UP+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE) 
    end if

    if (dom%mask_e%elts(EDGE*id_par+DG+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idNE_chd+1), ADJZONE) ! add associated chid node
       if (dom%mask_e%elts(EDGE*idSW_par+DG+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  
    end if

    if (dom%mask_e%elts(EDGE*id_par+RT+1) >= TOLRNZ) then
       call set_at_least (dom%mask_n%elts(idE_chd +1), ADJZONE) ! add associated chid node
       if (dom%mask_e%elts(EDGE*idW_par +RT+1) >= TOLRNZ) call set_at_least (dom%mask_n%elts(id_chd+1), ADJZONE)  
    end if

    ! Add two child edges if at least one neighbouring parent edge is active
    if ( dom%mask_e%elts(EDGE*id_par +UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idW_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd +UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
    end if

    ! Check all five UPLT and LORT triangle edges
    if ( dom%mask_e%elts(EDGE*id_par +DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ ) then
       call set_at_least (dom%mask_e%elts(EDGE*id_chd +DG+1),  ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
    end if

    ! Check all five triangle edges of LORT triangle of current cell and UPLT triangle of southern neighbour
    if  (dom%mask_e%elts(EDGE*id_par +RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idS_par+DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_chd+ RT+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id_par +UP+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idN_par+RT+1) >= TOLRNZ ) then

       call set_at_least (dom%mask_e%elts(EDGE*idN_chd +RT+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd +DG+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+UP+1), ADJZONE)
    end if

    if ( dom%mask_e%elts(EDGE*id_par +RT+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*id_par +DG+1) >= TOLRNZ .or. &
         dom%mask_e%elts(EDGE*idE_par+UP+1) >= TOLRNZ ) then

       call set_at_least (dom%mask_e%elts(EDGE*idE_chd +UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idE_chd +DG+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idNE_chd+RT+1), ADJZONE)
    end if

    ! Add node to restrict mask at coarse scale if at least one of six neighbouring edges at coarse scale is active
    ! (at positions corresponding to neighbouring nodes at finer scale)
    if ( dom%mask_e%elts(EDGE*id_par  +RT+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par  +DG+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*id_par  +UP+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idW_par +RT+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idSW_par+DG+1) >= ADJSPACE .or. &
         dom%mask_e%elts(EDGE*idS_par +UP+1) >= ADJSPACE ) &
         call set_at_least (dom%mask_n%elts(id_par+1), RESTRCT)
  end subroutine mask_adj_child

  
  subroutine mask_edges_consist (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Ensures consistency of adjacent zone edge mask for nearest neighbours.
    ! Adds child/parent edge to adjacent mask if any of 6 nearest neighbour edges are in adjacent mask.
    ! (modifies child E, NE, N neighbour edges and parent edges)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_chd, id_par, idN_chd, idS_chd, idE_chd, idW_chd, idNE_chd, idNW_chd, idSE_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idS_chd  = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    idW_chd  = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW_chd = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idSE_chd = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    ! RT child/parent edges
    if ( dom%mask_e%elts(EDGE*id_chd  +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSE_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) >= ADJZONE ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_chd +RT+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idE_chd+RT+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par +RT+1), ADJZONE)
    end if

    ! DG child/parent edges
    if ( dom%mask_e%elts(EDGE*idN_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE ) then

       call set_at_least(dom%mask_e%elts(EDGE*id_chd  +DG+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*idNE_chd+DG+1), ADJZONE)
       call set_at_least(dom%mask_e%elts(EDGE*id_par  +DG+1), ADJZONE)
    end if

    ! UP child/parent edges
    if ( dom%mask_e%elts(EDGE*idNW_chd+RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +UP+1) >= ADJZONE ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_chd +UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*idN_chd+UP+1), ADJZONE)
       call set_at_least (dom%mask_e%elts(EDGE*id_par +UP+1), ADJZONE)
    end if
  end subroutine mask_edges_consist

  
  subroutine mask_edges_consist2 (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Make parent edge active if any of 12 second neighouring child edges is active
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_chd, id_par, idN_chd, idE_chd, idS_chd, idW_chd, idNE_chd, idSW_chd, idNW_chd, idSE_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    idW_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    idSW_chd = idx (i_chd-2, j_chd-2, offs_chd, dims_chd)
    idS_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    idNW_chd = idx (i_chd-2, j_chd+2, offs_chd, dims_chd)
    idSE_chd = idx (i_chd+2, j_chd-2, offs_chd, dims_chd)
    
    ! RT parent edge
    if ( dom%mask_e%elts(EDGE*idSW_chd+RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idSW_chd+DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idS_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +DG+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idS_chd +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSE_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) >= ADJZONE ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+RT+1), ADJZONE)
    end if

    ! DG parent edge
    if ( dom%mask_e%elts(EDGE*idW_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idS_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +DG+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idS_chd +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idE_chd +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+DG+1), ADJZONE)
    end if

    ! UP parent edge
    if ( dom%mask_e%elts(EDGE*id_chd  +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idN_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW_chd +RT+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNW_chd+RT+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idN_chd +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*id_chd  +DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW_chd+DG+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW_chd +DG+1) >= ADJZONE .or. &
         
         dom%mask_e%elts(EDGE*idE_chd +UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idSW_chd+UP+1) >= ADJZONE .or. &
         dom%mask_e%elts(EDGE*idW_chd +UP+1) >= ADJZONE ) then

       call set_at_least (dom%mask_e%elts(EDGE*id_par+UP+1), ADJZONE)
    end if
  end subroutine mask_edges_consist2

  
  subroutine mask_edges_if_both_nodes (dom, i, j, zlev, offs, dims)
    ! Add edge to adjacent mask if both neighbouring nodes are active
    ! (also needed for div and gradv_e operator)
    ! Add edges required for flux restriction

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idE, idN, idNE

    id = idx (i, j, offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (dom%mask_n%elts(id +1) >= ADJZONE) then
       if (dom%mask_n%elts(idE +1) >= ADJZONE) call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), ADJZONE)
       if (dom%mask_n%elts(idNE+1) >= ADJZONE) call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), ADJZONE)
       if (dom%mask_n%elts(idN +1) >= ADJZONE) call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), ADJZONE)
    end if
  end subroutine mask_edges_if_both_nodes

  
  subroutine mask_nodes_if_all_edges (dom, i, j, zlev, offs, dims)
    ! Add node to adjacent mask if all 6 neighbouring edges are active
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idW, idS, idSW

    id = idx (i, j, offs, dims)
    
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    if ( dom%mask_e%elts(EDGE*id  +RT+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id  +DG+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*id  +UP+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idW +RT+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idSW+DG+1) >= ADJZONE .and. &
         dom%mask_e%elts(EDGE*idS +UP+1) >= ADJZONE ) then

       call set_at_least (dom%mask_n%elts(id+1), ADJZONE)
    end if
  end subroutine mask_nodes_if_all_edges

  
  subroutine prolong_node_adjzone (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Add parent to adjacent mask if child is in adjacent mask for prolongation of nodal quantities in inverse wavelet transform
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer,      intent(in)    :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)    :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)

    integer :: id_chd, id_par

    id_par = idx (i_par, j_par, offs_par, dims_par) 
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) 

    if (dom%mask_n%elts(id_chd+1) >= ADJZONE) call set_at_least (dom%mask_n%elts(id_par+1), ADJZONE)
  end subroutine prolong_node_adjzone

  
  subroutine set_masks (dom, p, i, j, zlev, offs, dims, mask)
    ! Sets all nodes and edges to value mask
     
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, mask, p, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id

    id = idx (i, j, offs, dims)

    dom%mask_n%elts(id+1)        = mask
    dom%mask_e%elts(id_edge(id)) = mask
  end subroutine set_masks

  
  subroutine init_masks
    implicit none
    integer :: d, num

    do d = 1, size(grid)
       call init (grid(d)%mask_n, 1)
       call init (grid(d)%mask_e, EDGE)
       call init (grid(d)%level,  1)
    end do

    do d = 1, size(grid)
       num = grid(d)%node%length-1
       call extend (grid(d)%mask_n, num,        TOLRNZ)
       call extend (grid(d)%mask_e, EDGE * num, TOLRNZ)
       call extend (grid(d)%level,  num,        min_level - 1)
    end do
  end subroutine init_masks

  
  subroutine set_at_least (mask, typ)
    implicit none
    integer, intent(inout) :: mask
    integer, intent(in)    :: typ

    if (mask < typ) mask = typ
  end subroutine set_at_least 

  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Second neighbour operator stencils for trend computation 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mask_second_neighbours
    ! Label nodes/edges that are up to second nearest neighbours of adajcent zone nodes as TRSK
    implicit none

    call apply_bdry (second_neigh, z_null, 0, 0) 
  end subroutine mask_second_neighbours

  
  subroutine second_neigh (dom, i, j, zlev,  offs, dims)
    ! Second neighbours (uses qe stencil)
     
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
   
    integer :: id, ii, jj
  
    id = idx (i, j, offs, dims)

    if (maxval (dom%mask_e%elts(id_edge(id))) >= ADJZONE) then
       do ii = -1, 1
          do jj = -1, 1
             call qe_trsk (dom, i+ii, j+jj, zlev, offs, dims)
          end do
       end do
    end if
  end subroutine second_neigh

  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    TRiSK operator stencils 
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mask_trsk
    ! Label nodes/edges required to compute TRiSK operators for active nodes and edges
    ! (used for testing)
    implicit none

    call apply_bdry (nodes_trsk, z_null, 0, 1)
    call apply_bdry (edges_trsk, z_null, 0, 0) 
  end subroutine mask_trsk

  
  subroutine nodes_trsk (dom, i, j, zlev, offs, dims)
    ! TRISK operator stencils needed for acive nodes
     
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id

    id = idx (i, j, offs, dims)
    
    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       call remap_trsk (dom, i, j, zlev,  offs, dims)

       if (Laplace_sclr == 1) call Laplacian_sclr_trsk      (dom, i, j, zlev, offs, dims) 
       if (Laplace_sclr == 2) call hyperLaplacian_sclr_trsk (dom, i, j, zlev, offs, dims)
    end if
  end subroutine nodes_trsk

  
  subroutine edges_trsk (dom, i, j, zlev, offs, dims)
    ! TRISK operator stencils needed for active edges
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    integer :: id

    id = idx (i, j, offs, dims)

    if (maxval (dom%mask_e%elts(id_edge(id))) >= ADJZONE) then
       call gradB_trsk (dom, i, j, zlev,  offs, dims) ! gradient of Bernoulli function
       call gradK_trsk (dom, i, j, zlev,  offs, dims) ! gradient of kinetic energy
       call Qperp_trsk (dom, i, j, zlev,  offs, dims) ! Qperp

       if (Laplace_rotu == 1) call      Laplacian_u_trsk (dom, i, j, zlev,  offs, dims)
       if (Laplace_rotu == 2) call hyperLaplacian_u_trsk (dom, i, j, zlev,  offs, dims)
    end if
  end subroutine edges_trsk

  
  subroutine gradB_trsk (dom, i, j, zlev,  offs, dims)
    ! Gradient of Bernoulli function
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    integer :: id

    id = idx (i, j, offs, dims)
    
    call set_at_least (dom%mask_e%elts(EDGE*id+RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id+UP+1), TRSK)
  end subroutine gradB_trsk

  
  subroutine gradK_trsk (dom, i, j, zlev,  offs, dims)
    ! Gradient of kinetic energy
     
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    call ke_trsk (dom, i,   j,   zlev, offs, dims)
    call ke_trsk (dom, i+1, j,   zlev, offs, dims)
    call ke_trsk (dom, i+1, j+1, zlev, offs, dims)
    call ke_trsk (dom, i,   j+1, zlev, offs, dims)
  end subroutine gradK_trsk

  
  subroutine Qperp_trsk (dom, i, j, zlev,  offs, dims)
    ! Gradient of kinetic energy
      
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
   
    integer :: ii, jj

    call div_grad_trsk (dom, i+1, j,   zlev, offs, dims)
    call div_grad_trsk (dom, i+1, j+1, zlev, offs, dims)
    call div_grad_trsk (dom, i,   j+1, zlev, offs, dims)

    do ii = -1, 1
       do jj = -1, 1
          call qe_trsk (dom, i+ii, j+jj, zlev, offs, dims)
       end do
    end do
  end subroutine Qperp_trsk

  
  subroutine remap_trsk (dom, i, j, zlev,  offs, dims)
    ! Remap stencil
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
   
    integer :: idE, idNE, idN

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    call set_at_least (dom%mask_n%elts(idE +1), TRSK)
    call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    call set_at_least (dom%mask_n%elts(idN +1), TRSK)
  end subroutine remap_trsk

  
  subroutine ke_trsk (dom, i, j, zlev,  offs, dims)
    ! Kinetic energy stencil
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
   
    integer :: id, idS, idSW, idW

    id   = idx (i,   j,   offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_e%elts(EDGE*id  +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +UP+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS +UP+1), TRSK)
  end subroutine ke_trsk

  
  subroutine qe_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil for qe
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
   
    integer :: id, idE, idN, idNE, idS, idSE, idW

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idS    = idx (i,   j-1, offs, dims)
    idSE   = idx (i+1, j-1, offs, dims)

    ! Circulation stencil
    call set_at_least (dom%mask_e%elts(EDGE*id +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id +DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id +UP+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idE+UP+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idN+RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idW+RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS+UP+1), TRSK)

    ! Potential vorticity stencil
    call set_at_least (dom%mask_n%elts(id  +1), TRSK)
    call set_at_least (dom%mask_n%elts(idE +1), TRSK)
    call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    call set_at_least (dom%mask_n%elts(idN +1), TRSK)
    call set_at_least (dom%mask_n%elts(idW +1), TRSK)
    call set_at_least (dom%mask_n%elts(idSE+1), TRSK)
  end subroutine qe_trsk

  
  subroutine Laplacian_sclr_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil of Laplacian hyperdiffusion of scalars
      
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    call div_grad_trsk (dom, i, j, zlev,  offs, dims)
  end subroutine Laplacian_sclr_trsk

  
  subroutine hyperLaplacian_sclr_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil of Laplacian hyperdiffusion of scalars
       
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    integer :: ii, jj
    
    do ii = -1, 1
       do jj = -1, 1
          call div_grad_trsk (dom, i+ii, j+jj, zlev, offs, dims)
       end do
    end do
  end subroutine hyperLaplacian_sclr_trsk

  
  subroutine hyperLaplacian_u_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil of Laplacian hyperdiffusion of velocity
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    integer :: ii, jj

    do ii = -1, 1
       do jj = -1, 1
          call Laplacian_u_trsk (dom, i+ii, j+jj, zlev, offs, dims)
       end do
    end do
  end subroutine hyperLaplacian_u_trsk

  
  subroutine Laplacian_u_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil for Laplacian(u) operators
     
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
  
    call divu_trsk (dom, i,   j,   zlev, offs, dims)
    call divu_trsk (dom, i+1, j,   zlev, offs, dims)
    call divu_trsk (dom, i+1, j+1, zlev, offs, dims)
    call divu_trsk (dom, i,   j+1, zlev, offs, dims)
  end subroutine Laplacian_u_trsk

  
  subroutine divu_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil for divu operator
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idE, idNE, idN, idW, idS, idSW

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idSW   = idx (i-1, j-1, offs, dims)
    idS    = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_e%elts(EDGE*id  +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +UP+1), TRSK) 
    call set_at_least (dom%mask_e%elts(EDGE*idW +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS +UP+1), TRSK)
  end subroutine divu_trsk

  
  subroutine div_grad_trsk (dom, i, j, zlev,  offs, dims)
    ! Stencil for flux-divergence operator
    
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id, idE, idNE, idN, idW, idSW, idS

    id     = idx (i,   j,   offs, dims)
    idE    = idx (i+1, j,   offs, dims)
    idNE   = idx (i+1, j+1, offs, dims)
    idN    = idx (i,   j+1, offs, dims)
    idW    = idx (i-1, j,   offs, dims)
    idSW   = idx (i-1, j-1, offs, dims)
    idS    = idx (i,   j-1, offs, dims)

    call set_at_least (dom%mask_n%elts(id  +1), TRSK)
    call set_at_least (dom%mask_n%elts(idE +1), TRSK)
    call set_at_least (dom%mask_n%elts(idNE+1), TRSK)
    call set_at_least (dom%mask_n%elts(idN +1), TRSK)
    call set_at_least (dom%mask_n%elts(idW +1), TRSK)
    call set_at_least (dom%mask_n%elts(idSW+1), TRSK)
    call set_at_least (dom%mask_n%elts(idS +1), TRSK)

    call set_at_least (dom%mask_e%elts(EDGE*id  +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*id  +UP+1), TRSK) 
    call set_at_least (dom%mask_e%elts(EDGE*idW +RT+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idSW+DG+1), TRSK)
    call set_at_least (dom%mask_e%elts(EDGE*idS +UP+1), TRSK)
  end subroutine div_grad_trsk

  
end module mask_mod
