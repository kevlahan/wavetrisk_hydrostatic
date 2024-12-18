module adapt_mod
  use refine_patch_mod
  use multi_level_mod
  implicit none

  interface compress_wavelets_scalar
     procedure :: compress_wavelets_scalar_0, compress_wavelets_scalar_1
  end interface compress_wavelets_scalar
  
  interface WT_after_scalar
     procedure :: WT_after_scalar_0, WT_after_scalar_1
  end interface WT_after_scalar 
contains
  subroutine adapt (set_thresholds, type)
    ! Determines significant wavelets, adaptive grid and all masks associated with adaptive grid
    implicit none
    external           :: set_thresholds
    logical, optional  :: type ! recalculate thresholds
    
    integer :: k, l, d
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
    
    ! Active zone at all scales
    call mask_active

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

    ! Add nodes/edges required for TRiSK operators
    ! (only affects ZERO mask nodes/edges that are second nearest neighbours of nodes in adjacent mask)
    call mask_trsk

    ! Set insignificant wavelet coefficients to zero
    if (local_type) call compress_wavelets (wav_coeff)
    
    ! Evaluate sol_mean, topography and penalization (as defined in test case) on new grid
    call update 
  end subroutine adapt

 subroutine init_adapt_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_comm_mod
    call init_refine_patch_mod
    
    max_level_exceeded = .false.
    initialized        = .true.
  end subroutine init_adapt_mod

  subroutine compress_wavelets (wav)
    ! Sets wavelets associated with inactive grid points to zero
    implicit none
    type(Float_Field), dimension(:,:), target :: wav

    integer :: d, k, l, v

    call update_bdry (wav, NONE)
    
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
    type(Float_Field), target :: wav

    integer :: d, k, l

    call update_bdry (wav, NONE)
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
    type(Float_Field), dimension(:), target :: wav

    integer :: d, k, l

    call update_bdry (wav, NONE)
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
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    implicit none
    type(Float_Field), dimension(:,:), target :: scaling, wavelet
    integer, optional                         :: l_start0

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

    call update_bdry (scaling, NONE, 16)

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

  subroutine WT_after_scalar_0 (scaling, wavelet, l_start0)
    ! Compute wavelets and interpolate solution onto adaptive grid (including ZERO mask cells)
    implicit none
    type(Float_Field), target :: scaling, wavelet
    integer, optional         :: l_start0

    integer :: d, j, k, l, l_start

    call zero_float (wavelet)

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_bdry (scaling, NONE)

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
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0
    
    integer :: d, j, k, l, l_start

    call zero_float (wavelet)

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_bdry (scaling, NONE)

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
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0
    
    integer :: d, j, k, l, l_start

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

    call update_bdry (scaling, NONE)

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
    call inverse_wavelet_transform (wav_coeff, sol, jmin_in=old_level_start)
    sol%bdry_uptodate = .false.
    call update_bdry (sol, NONE, 17)
  end subroutine fill_up_grid_and_IWT    
end module adapt_mod
