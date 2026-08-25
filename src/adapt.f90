module adapt_mod
  use iso_fortran_env,   only : error_unit, int64

  use kind_mod,         only : dp
  
  use shared_mod,       only : ADJZONE, EDGE, level_start, level_end, &
       N_BDRY, N_VARIABLE, n_node_old, n_patch_old, NONE, &
       scalars, S_VELO, z_null, zlevels
  
  use arch_mod,         only : rank
  use comm_mpi_mod,     only : update_bdry
  use domain_ops_mod,   only : apply_interscale_d, apply_onescale_d, apply_to_penta_d
  use init_mod,         only : noarg_sub, set_thresholds, update
  use refine_patch_mod, only : fill_up_level, max_level_exceeded, post_refine, refine
  use utils_mod,        only : zero_float


  use domain_mod,  only :  Domain, Float_Field, grid, idx, scalar, sol, velo, &
       wav_coeff, wc_s, wc_u

  use mask_mod,    only : complete_masks, init_masks_zero, mask_active, mask_adj_finer_scale, mask_adj_same_scale, &
       mask_restrict_same_scale, mask_second_neighbours

  use wavelet_mod, only : compute_scalar_wavelets, compute_velo_wavelets, compute_velo_wavelets_penta, &
       inverse_scalar_transform, inverse_velo_transform, inverse_wavelet_transform, restrict_velo

  implicit none

  private
  public :: adapt, compress_wavelets_scalar, fill_up_grid_and_IWT
  public :: WT_after_scalar, WT_after_step, WT_after_velo 
  public :: begin_inverse_wavelet_transform_shadow
  public :: compare_inverse_wavelet_transform_shadow
  
  interface compress_wavelets_scalar
     procedure :: compress_wavelets_scalar_0, compress_wavelets_scalar_1
  end interface compress_wavelets_scalar
  
  interface WT_after_scalar
     procedure :: WT_after_scalar_0, WT_after_scalar_1
  end interface WT_after_scalar

  abstract interface
     subroutine Scalar_Wavelet_Validator (scaling,first_level)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(in) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       integer, intent(in) :: first_level
     end subroutine Scalar_Wavelet_Validator

     subroutine Outer_Vector_Wavelet_Validator (scaling,wavelet_level)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(in) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       integer, intent(in) :: wavelet_level
     end subroutine Outer_Vector_Wavelet_Validator

     subroutine Velocity_Restriction_Validator (scaling,coarse_level)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(in) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       integer, intent(in) :: coarse_level
     end subroutine Velocity_Restriction_Validator

     subroutine Adaptation_Mask_Validator(stage)
       integer, intent(in) :: stage
     end subroutine Adaptation_Mask_Validator

     subroutine Wavelet_Compression_Validator(wavelet,stage)
       import :: Float_Field

       type(Float_Field), intent(in) :: wavelet(:,:)
       integer, intent(in) :: stage
     end subroutine Wavelet_Compression_Validator
  end interface

  type(Float_Field), allocatable, target, save :: &
       inverse_wavelet_input(:,:)
  type(Float_Field), allocatable, target, save :: &
       inverse_scaling_expected(:,:)
  logical, save :: inverse_wavelet_shadow_ready = .false.
  integer(int64), save :: inverse_wavelet_shadow_allocations = 0_int64
  integer(int64), save :: inverse_wavelet_allocation_checkpoint = 0_int64

  
contains

  subroutine begin_inverse_wavelet_transform_shadow(wavelet,scaling)
    ! Reconstruct an expected solution through the independent scalar and
    ! velocity inverse-transform entry points.  This covers scalar lifting,
    ! outer and inner velocity prolongation, and pentagon corrections.

    implicit none

    type(Float_Field), intent(in) :: wavelet(:,:)
    type(Float_Field), intent(in) :: scaling(:,:)

    integer :: v

    if (inverse_wavelet_shadow_ready) then
       error stop &
            "begin_inverse_wavelet_transform_shadow: shadow is still active"
    end if
    if (size(wavelet,1) /= size(scaling,1) .or. &
         size(wavelet,2) /= size(scaling,2)) then
       error stop &
            "begin_inverse_wavelet_transform_shadow: field layout differs"
    end if
    call ensure_inverse_field_storage(inverse_wavelet_input,wavelet)
    call ensure_inverse_field_storage(inverse_scaling_expected,scaling)
    call copy_inverse_field_family(inverse_wavelet_input,wavelet)
    call copy_inverse_field_family(inverse_scaling_expected,scaling)

    do v = scalars(1),scalars(2)
       call inverse_scalar_transform( &
            inverse_wavelet_input(v,:),inverse_scaling_expected(v,:))
    end do
    call inverse_velo_transform( &
         inverse_wavelet_input(S_VELO,:), &
         inverse_scaling_expected(S_VELO,:))

    inverse_wavelet_shadow_ready = .true.
    inverse_wavelet_allocation_checkpoint = &
         inverse_wavelet_shadow_allocations

  end subroutine begin_inverse_wavelet_transform_shadow


  subroutine compare_inverse_wavelet_transform_shadow(scaling)
    ! Require the combined production inverse transform to match the
    ! independently orchestrated scalar/vector reconstruction bit-for-bit.

    implicit none

    type(Float_Field), intent(in) :: scaling(:,:)

    integer :: d
    integer :: i
    integer :: k
    integer :: v
    integer(int64) :: actual_bits
    integer(int64) :: expected_bits

    if (.not. inverse_wavelet_shadow_ready) then
       error stop &
            "compare_inverse_wavelet_transform_shadow: shadow unavailable"
    end if
    if (inverse_wavelet_shadow_allocations /= &
         inverse_wavelet_allocation_checkpoint) then
       error stop &
            "compare_inverse_wavelet_transform_shadow: storage changed"
    end if
    if (size(scaling,1) /= size(inverse_scaling_expected,1) .or. &
         size(scaling,2) /= size(inverse_scaling_expected,2)) then
       error stop &
            "compare_inverse_wavelet_transform_shadow: layout differs"
    end if
    do k = 1,size(scaling,2)
       do v = 1,size(scaling,1)
          if (size(scaling(v,k)%data) /= size(grid)) then
             error stop &
                  "compare_inverse_wavelet_transform_shadow: Domain differs"
          end if
          do d = 1,size(grid)
             if (scaling(v,k)%data(d)%length /= &
                  inverse_scaling_expected(v,k)%data(d)%length) then
                error stop &
                     "compare_inverse_wavelet_transform_shadow: extent differs"
             end if
             do i = 1,scaling(v,k)%data(d)%length
                actual_bits = transfer( &
                     scaling(v,k)%data(d)%elts(i),0_int64)
                expected_bits = transfer( &
                     inverse_scaling_expected(v,k)%data(d)%elts(i), &
                     0_int64)
                if (actual_bits == expected_bits) cycle
                write(error_unit,'(/,a,i0,a)') &
                     "Rank ",rank,": inverse wavelet reconstruction differs"
                write(error_unit,'(a,i0,a,i0,a,i0,a,i0)') &
                     "  variable = ",v,", level slot = ",k, &
                     ", Domain = ",d,", index = ",i
                write(error_unit,'(a,2(es24.16,1x))') &
                     "  production, expected = ", &
                     scaling(v,k)%data(d)%elts(i), &
                     inverse_scaling_expected(v,k)%data(d)%elts(i)
                flush(error_unit)
                error stop "inverse-wavelet reconstruction comparison failed"
             end do
          end do
       end do
    end do
    inverse_wavelet_shadow_ready = .false.

  end subroutine compare_inverse_wavelet_transform_shadow


  subroutine ensure_inverse_field_storage(storage,source)

    implicit none

    type(Float_Field), allocatable, target, intent(inout) :: storage(:,:)
    type(Float_Field), intent(in) :: source(:,:)

    integer :: d
    integer :: k
    integer :: v
    logical :: shape_differs

    shape_differs = .false.
    if (allocated(storage)) then
       shape_differs = size(storage,1) /= size(source,1) .or. &
            size(storage,2) /= size(source,2)
       if (shape_differs) deallocate(storage)
    end if
    if (.not. allocated(storage)) then
       allocate(storage(size(source,1),size(source,2)))
       inverse_wavelet_shadow_allocations = &
            inverse_wavelet_shadow_allocations+1_int64
    end if

    do k = 1,size(source,2)
       do v = 1,size(source,1)
          if (.not. allocated(source(v,k)%data)) then
             error stop "ensure_inverse_field_storage: source is absent"
          end if
          if (allocated(storage(v,k)%data)) then
             if (size(storage(v,k)%data) /= &
                  size(source(v,k)%data)) then
                deallocate(storage(v,k)%data)
             end if
          end if
          if (.not. allocated(storage(v,k)%data)) then
             allocate(storage(v,k)%data(size(source(v,k)%data)))
             inverse_wavelet_shadow_allocations = &
                  inverse_wavelet_shadow_allocations+1_int64
          end if
          do d = 1,size(source(v,k)%data)
             if (.not. allocated(source(v,k)%data(d)%elts)) then
                error stop &
                     "ensure_inverse_field_storage: source data is absent"
             end if
             if (allocated(storage(v,k)%data(d)%elts)) then
                if (size(storage(v,k)%data(d)%elts) /= &
                     size(source(v,k)%data(d)%elts)) then
                   deallocate(storage(v,k)%data(d)%elts)
                end if
             end if
             if (.not. allocated(storage(v,k)%data(d)%elts)) then
                allocate(storage(v,k)%data(d)%elts( &
                     size(source(v,k)%data(d)%elts)))
                inverse_wavelet_shadow_allocations = &
                     inverse_wavelet_shadow_allocations+1_int64
             end if
          end do
       end do
    end do

  end subroutine ensure_inverse_field_storage


  subroutine copy_inverse_field_family(destination,source)

    implicit none

    type(Float_Field), intent(inout) :: destination(:,:)
    type(Float_Field), intent(in) :: source(:,:)

    integer :: d
    integer :: k
    integer :: v

    do k = 1,size(source,2)
       do v = 1,size(source,1)
          destination(v,k)%pos = source(v,k)%pos
          destination(v,k)%bdry_tag = source(v,k)%bdry_tag
          destination(v,k)%bdry_uptodate = &
               source(v,k)%bdry_uptodate
          do d = 1,size(source(v,k)%data)
             destination(v,k)%data(d)%length = &
                  source(v,k)%data(d)%length
             destination(v,k)%data(d)%elts = &
                  source(v,k)%data(d)%elts
          end do
       end do
    end do

  end subroutine copy_inverse_field_family

  
  subroutine adapt ( &
       set_thresholds,type,validate_adaptation_masks, &
       validate_wavelet_compression)
    ! Determines significant wavelets, adaptive grid and all masks associated with adaptive grid
    
    implicit none
    
    procedure (noarg_sub)         :: set_thresholds
    logical, intent(in), optional :: type ! recalculate thresholds
    procedure(Adaptation_Mask_Validator), optional :: &
         validate_adaptation_masks
    procedure(Wavelet_Compression_Validator), optional :: &
         validate_wavelet_compression
    
    logical :: local_type

    n_patch_old = grid%patch%length
    n_node_old  = grid%node%length

    ! Recalculate thresholds
    if (present(type)) then
       local_type = type
    else
       local_type = .true.
    end if

    if (local_type) call set_thresholds

    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(0)

    ! Initialize all masks to ZERO at scales > level_start
    call init_masks_zero
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(1)
    
    ! Active zone at all scales
    call mask_active
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(2)

    ! Adjacent zone at same scale
    call mask_adj_same_scale
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(3)

    ! Mask for restriction at same scale
    call mask_restrict_same_scale
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(4)

    ! Determine whether any new patches are required
    if (refine ()) call post_refine

    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(5)

    ! Adjacent zone at finer scale
    call mask_adj_finer_scale
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(6)

    ! Ensure consistency between node and edge masks
    call complete_masks
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(7)

    ! Label ZERO mask nodes/edges that are second nearest neighbours of nodes in adjacent zone mask where trends are be computed
    call mask_second_neighbours
    if (present(validate_adaptation_masks)) &
         call validate_adaptation_masks(8)

    ! Set insignificant wavelet coefficients to zero
    if (local_type) then
       if (present(validate_wavelet_compression)) then
          call compress_wavelets( &
               wav_coeff,validate_wavelet_compression)
       else
          call compress_wavelets(wav_coeff)
       end if
    end if
    
    ! Evaluate sol_mean, topography and penalization (as defined in test case) on new grid
    call update 
  end subroutine adapt

  
  subroutine compress_wavelets (wav,validate_compression)
    ! Sets wavelets associated with inactive grid points to zero
    
    implicit none
    
    type(Float_Field), intent(inout), target :: wav(:,:)
    procedure(Wavelet_Compression_Validator), optional :: &
         validate_compression

    integer :: d, k, l, v

    call update_bdry (wav, NONE, 901)
    if (present(validate_compression)) &
         call validate_compression(wav,0)
    
    do k = 1, size (wav, 2)
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
    if (present(validate_compression)) &
         call validate_compression(wav,1)
  end subroutine compress_wavelets

  
  subroutine compress_wavelets_scalar_0 (wav)
    ! Sets scalar wavelets associated with inactive grid points to zero
    
    implicit none
    
    type(Float_Field), intent(inout), target :: wav

    integer :: d, l

    call update_bdry (wav, NONE, 902)
    do d = 1, size (grid)
       do l = level_start+1, level_end
          wc_s => wav%data(d)%elts
          call apply_onescale_d (compress_scalar, grid(d), l, z_null, 0, 1)
          nullify (wc_s)
       end do
    end do
    wav%bdry_uptodate = .false.
  end subroutine compress_wavelets_scalar_0

  
  subroutine compress_wavelets_scalar_1 (wav)
    ! Sets scalar wavelets associated with inactive grid points to zero
    
    implicit none
    
    type(Float_Field), intent(inout), target :: wav(:)

    integer :: d, k, l

    call update_bdry (wav, NONE, 903)
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
  end subroutine compress_wavelets_scalar_1

  
  subroutine compress_wavelets_velo (wav)
    ! Sets wavelets associated with inactive grid points to zero
    
    implicit none
    
    type(Float_Field), intent(inout), target :: wav(:)

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
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    
    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    
    if (dom%mask_n%elts(id_i) < ADJZONE) wc_s(id_i) = 0.0_dp
  end subroutine compress_scalar

  
  subroutine compress_vector (dom, i, j, zlev, offs, dims)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    
    integer :: e, id

    id = idx (i, j, offs, dims)
    
    do e = 1, EDGE
       if (dom%mask_e%elts(EDGE*id+e) < ADJZONE) wc_u(EDGE*id+e) = 0.0_dp
    end do
  end subroutine compress_vector

  
  subroutine WT_after_step ( &
       scaling,wavelet,l_start0,validate_scalar_wavelets, &
       validate_outer_vector_wavelets,validate_velocity_restriction)
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    
    implicit none
    
    type(Float_Field), target,   intent(inout) :: scaling(:,:), wavelet(:,:)
    integer,           optional, intent(in)    :: l_start0
    procedure(Scalar_Wavelet_Validator), optional :: &
         validate_scalar_wavelets
    procedure(Outer_Vector_Wavelet_Validator), optional :: &
         validate_outer_vector_wavelets
    procedure(Velocity_Restriction_Validator), optional :: &
         validate_velocity_restriction

    integer :: d, k, l, l_start, v

    call zero_float (wavelet)

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

    call update_bdry (scaling, NONE, 904)
    if (present(l_start0) .and. &
         present(validate_velocity_restriction)) then
       call validate_velocity_restriction(scaling,level_start-1)
    end if

    do l = l_start, level_end-1
       do k = 1, size(scaling,2)
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
             nullify (velo, wc_u)
          end do
       end do
       if (present(validate_outer_vector_wavelets)) then
          call validate_outer_vector_wavelets(scaling,l)
       end if
       do k = 1, size(scaling,2)
          do d = 1,size(grid)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_to_penta_d( &
                  compute_velo_wavelets_penta,grid(d),l,z_null)
             nullify(velo,wc_u)
          end do
          wavelet(:,k)%bdry_uptodate = .false.
       end do
    end do
    if (present(validate_scalar_wavelets)) then
       call validate_scalar_wavelets(scaling,l_start)
    end if
    call compress_wavelets (wavelet)
    call inverse_wavelet_transform (wavelet, scaling)
  end subroutine WT_after_step

  
  subroutine WT_after_scalar_0 (scaling, wavelet, l_start0)
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    
    implicit none
    
    type(Float_Field), target,   intent(inout) :: scaling, wavelet
    integer,           optional, intent(in)    :: l_start0

    integer :: d, l, l_start

    call zero_float (wavelet)

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_bdry (scaling, NONE, 905)

    do l = l_start, level_end-1
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          wc_s   => wavelet%data(d)%elts
          call apply_interscale_d (compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
          nullify (scalar, wc_s)
       end do
       wavelet%bdry_uptodate = .false.
    end do
    call compress_wavelets_scalar (wavelet)
    call inverse_scalar_transform (wavelet, scaling)
  end subroutine WT_after_scalar_0

  
  subroutine WT_after_scalar_1 (scaling, wavelet, l_start0)
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    
    implicit none

    type(Float_Field), target,   intent(inout) :: scaling(:), wavelet(:)
    integer,           optional, intent(in)    :: l_start0
    
    integer :: d, k, l, l_start

    call zero_float (wavelet)

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_bdry (scaling, NONE, 906)

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
  end subroutine WT_after_scalar_1

  
  subroutine WT_after_velo (scaling, wavelet, l_start0)
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    
    implicit none

    type(Float_Field), target,   intent(inout) :: scaling(:), wavelet(:)
    integer,           optional, intent(in)    :: l_start0
    
    integer :: d, k, l, l_start

    call zero_float (wavelet)

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

    call update_bdry (scaling, NONE, 907)

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

  
  subroutine fill_up_grid_and_IWT (l)
    ! Fills grid up to level l and does inverse wavelet transform of solution onto grid
    
    implicit none
    
    integer, intent(in) :: l
    
    integer :: old_level_start
    
    old_level_start = level_start
    do while (level_start < l)
       if (rank == 0) write(6,'(a,i2)') 'Filling up level ', level_start+1
       call fill_up_level
    end do
    call inverse_wavelet_transform (wav_coeff, sol, jmin_in=old_level_start)
    
    sol%bdry_uptodate = .false.
    call update_bdry (sol, NONE, 908)
  end subroutine fill_up_grid_and_IWT

  
end module adapt_mod
