module adapt_mod
  use kind_mod,         only : dp
  
  use shared_mod,       only : ADJZONE, EDGE, level_start, level_end, &
       N_BDRY, N_VARIABLE, n_node_old, n_patch_old, NONE, &
       scalars, S_VELO, z_null, zlevels
  
  use arch_mod,         only : rank
  use comm_mpi_mod,     only : update_bdry, update_bdry1
  use domain_ops_mod,   only : apply_interscale_d, apply_onescale_d, apply_to_penta_d
  use init_mod,         only : noarg_sub, set_thresholds, update
  use refine_patch_mod, only : fill_up_level, max_level_exceeded, post_refine, refine
  use utils_mod,        only : zero_float


  use domain_mod,  only :  Domain, Float_Field, grid, idx, scalar, sol, velo, &
       wav_coeff, wc_s, wc_u

  use mask_mod,    only : complete_masks, init_masks_zero, mask_active, mask_adj_finer_scale, mask_adj_same_scale, &
       mask_restrict_same_scale, mask_second_neighbours

  use wavelet_mod, only : compute_scalar_wavelets, compute_velo_wavelets, compute_velo_wavelets_penta, &
       inverse_scalar_transform, inverse_velo_transform, &
       inverse_wavelet_transform, restrict_velo

  implicit none

  private
  public :: adapt, compress_wavelets_scalar, fill_up_grid_and_IWT
  public :: WT_after_scalar, WT_after_step, WT_after_velo 
  
  interface compress_wavelets_scalar
     procedure :: compress_wavelets_scalar_0, compress_wavelets_scalar_1
  end interface compress_wavelets_scalar
  
  interface WT_after_scalar
     procedure :: WT_after_scalar_0, WT_after_scalar_1
  end interface WT_after_scalar

  abstract interface
     subroutine Scalar_Wavelet_Validator (scaling,wavelet,first_level)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(in) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       type(Float_Field), intent(inout) :: &
            wavelet(1:N_VARIABLE,1:zlevels)
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

     subroutine Wavelet_Compression_Handler (wavelet)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(inout) :: &
            wavelet(1:N_VARIABLE,1:zlevels)
     end subroutine Wavelet_Compression_Handler

     subroutine Inverse_Scalar_Boundary_Handler (scaling,level)
       import :: Float_Field, N_VARIABLE, zlevels

       type(Float_Field), intent(inout) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       integer, intent(in) :: level
     end subroutine Inverse_Scalar_Boundary_Handler

  end interface

  abstract interface
     subroutine Inverse_Wavelet_Handler ( &
          wavelet,scaling,jmin,jmax, &
          refresh_scalar_boundary,refresh_vector_boundary)
       import :: Float_Field, Inverse_Scalar_Boundary_Handler, &
            N_VARIABLE, zlevels

       type(Float_Field), intent(inout) :: &
            wavelet(1:N_VARIABLE,1:zlevels)
       type(Float_Field), intent(inout) :: &
            scaling(1:N_VARIABLE,1:zlevels)
       integer, intent(in) :: jmin
       integer, intent(in) :: jmax
       procedure(Inverse_Scalar_Boundary_Handler) :: &
            refresh_scalar_boundary
       procedure(Inverse_Scalar_Boundary_Handler) :: &
            refresh_vector_boundary
     end subroutine Inverse_Wavelet_Handler
  end interface

  
contains

  subroutine refresh_native_inverse_scalar_boundary (scaling,level)
    ! Keep the established Domain boundary operator above the compact block
    ! layer so parallel_block_mpi_mod retains its original module ordering.

    implicit none

    type(Float_Field), intent(inout) :: &
         scaling(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: level

    call update_bdry1( &
         scaling(scalars(1):scalars(2),:),level,level,834)

  end subroutine refresh_native_inverse_scalar_boundary


  subroutine refresh_native_inverse_vector_boundary (scaling,level)
    ! Refresh the legacy geometric boundary after a compact vector phase.

    implicit none

    type(Float_Field), intent(inout) :: &
         scaling(1:N_VARIABLE,1:zlevels)
    integer, intent(in) :: level

    call update_bdry1(scaling(S_VELO,:),level,level,835)

  end subroutine refresh_native_inverse_vector_boundary


  subroutine adapt (set_thresholds, type,validate_significance)
    ! Determines significant wavelets, adaptive grid and all masks associated with adaptive grid
    
    implicit none
    
    procedure (noarg_sub)         :: set_thresholds
    logical, intent(in), optional :: type ! recalculate thresholds
    procedure(noarg_sub), optional :: validate_significance
    
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

    ! Initialize all masks to ZERO at scales > level_start
    call init_masks_zero

    ! A retained compact wavelet image installs the authoritative direct
    ! threshold seed after mask initialization.  The callback also validates
    ! the seed and releases topology-dependent block state before refinement.
    if (present(validate_significance)) call validate_significance
    
    ! Active zone at all scales
    if (present(validate_significance)) then
       call mask_active(.true.)
    else
       call mask_active
    end if

    ! Adjacent zone at same scale
    call mask_adj_same_scale

    ! Mask for restriction at same scale
    call mask_restrict_same_scale

    ! Determine whether any new patches are required
    if (refine ()) call post_refine

    ! Adjacent zone at finer scale
    call mask_adj_finer_scale

    ! Ensure consistency between node and edge masks
    call complete_masks

    ! Label ZERO mask nodes/edges that are second nearest neighbours of nodes in adjacent zone mask where trends are be computed
    call mask_second_neighbours

    ! Set insignificant wavelet coefficients to zero
    if (local_type) call compress_wavelets (wav_coeff)
    
    ! Evaluate sol_mean, topography and penalization (as defined in test case) on new grid
    call update 
  end subroutine adapt

  
  subroutine compress_wavelets (wav)
    ! Sets wavelets associated with inactive grid points to zero
    
    implicit none
    
    type(Float_Field), intent(inout), target :: wav(:,:)

    integer :: d, k, l, v

    call update_bdry (wav, NONE, 901)
    
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
       validate_outer_vector_wavelets,validate_velocity_restriction, &
       native_wavelet_output,prepare_native_compression, &
       activate_native_compression,activate_native_inverse)
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
    logical, optional, intent(in) :: native_wavelet_output
    procedure(Wavelet_Compression_Handler), optional :: &
         prepare_native_compression
    procedure(Wavelet_Compression_Handler), optional :: &
         activate_native_compression
    procedure(Inverse_Wavelet_Handler), optional :: &
         activate_native_inverse

    integer :: d, k, l, l_start, v
    logical :: native_output

    native_output = .false.
    if (present(native_wavelet_output)) then
       native_output = native_wavelet_output
    end if
    if (native_output .and. &
         (.not. present(l_start0) .or. &
          .not. present(validate_scalar_wavelets) .or. &
          .not. present(validate_outer_vector_wavelets) .or. &
          .not. present(prepare_native_compression) .or. &
          .not. present(activate_native_compression) .or. &
          .not. present(activate_native_inverse))) then
       error stop "native wavelet output requires complete production callbacks"
    end if

    call zero_float (wavelet)

    if (present(l_start0)) then
       l_start = l_start0
       if (l_start == level_start-1) then
          do k = 1, size(scaling,2)
             do d = 1, size(grid)
                velo => scaling(S_VELO,k)%data(d)%elts
                call apply_interscale_d ( &
                     restrict_velo,grid(d),level_start-1,k,0,0)
                nullify (velo)
             end do
          end do
       end if
    else
       l_start = level_start
    end if
    if (native_output .and. l_start /= level_start-1 .and. &
         l_start /= level_start) then
       error stop "native wavelet output has an invalid first level"
    end if
    if (native_output .and. l_start == level_start-1 .and. &
         .not. present(validate_velocity_restriction)) then
       error stop "complete native wavelet output requires restriction"
    end if

    call update_bdry (scaling, NONE, 904)
    if (l_start == level_start-1 .and. &
         present(validate_velocity_restriction)) then
       call validate_velocity_restriction(scaling,level_start-1)
    end if

    do l = l_start, level_end-1
       if (.not. native_output) then
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
       end if
       if (present(validate_outer_vector_wavelets)) then
          call validate_outer_vector_wavelets(scaling,l)
       end if
       if (.not. native_output) then
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
       end if
    end do
    if (present(validate_scalar_wavelets)) then
       call validate_scalar_wavelets(scaling,wavelet,l_start)
    end if
    if (native_output) then
       call prepare_native_compression(wavelet)
       call activate_native_compression(wavelet)
       ! Establish the geometric aliases consumed by the native inverse.
       ! No legacy reconstruction is executed on the production path.
       call update_bdry1( &
            wavelet,max(l_start,level_start),level_end,802)
       call update_bdry1(scaling,l_start,level_end,803)
       scaling%bdry_uptodate = .false.
       call activate_native_inverse( &
            wavelet,scaling,l_start,level_end, &
            refresh_native_inverse_scalar_boundary, &
            refresh_native_inverse_vector_boundary)
    else
       call compress_wavelets (wavelet)
       call inverse_wavelet_transform(wavelet,scaling)
    end if
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
