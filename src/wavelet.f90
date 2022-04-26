module wavelet_mod
  use domain_mod
  use comm_mpi_mod
  use utils_mod
  implicit none
  real(8), dimension(9) :: Iu_Base_Wgt
contains
  subroutine forward_wavelet_transform (scaling, wavelet)
    ! Forward wavelet transform
    implicit none
    type(Float_Field), dimension(:,:), target :: scaling, wavelet

    integer :: k, l, d, v

    do l = level_end-1, level_start-1, -1
       ! Compute scalar wavelet coefficients
       call update_array_bdry (scaling(scalars(1):scalars(2),:), l+1, 1)

       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                scalar => scaling(v,k)%data(d)%elts
                wc_s   => wavelet(v,k)%data(d)%elts
                call apply_interscale_d (Compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
                nullify (scalar, wc_s)
             end do
          end do
       end do
       call update_array_bdry (wavelet(scalars(1):scalars(2),:), l+1, 2)

       ! Restrict scalars (sub-sample and lift) and velocity (average) to coarser grid
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                scalar => scaling(v,k)%data(d)%elts
                wc_s   => wavelet(v,k)%data(d)%elts
                call apply_interscale_d (Restrict_scalar, grid(d), l, z_null, 0, 1) ! +1 to include poles
                nullify (scalar, wc_s)
             end do
             velo => scaling(S_VELO,k)%data(d)%elts
             call apply_interscale_d (Restrict_velo, grid(d), l, z_null, 0, 0)
             nullify (velo)
          end do
       end do
    end do

    scaling%bdry_uptodate                          = .false.
    wavelet(scalars(1):scalars(2),:)%bdry_uptodate = .false.

    call update_vector_bdry (scaling(S_VELO,:), NONE, 3)

    ! Compute vector wavelet coefficients
    do l = level_end-1, level_start-1, -1
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             velo => scaling(S_VELO,k)%data(d)%elts
             call apply_interscale_d (Compute_velo_wavelets, grid(d), l, z_null, 0, 0)
             call apply_to_penta_d (Compute_velo_wavelets_penta, grid(d), l, z_null)
             nullify (wc_u, velo)
          end do
       end do
    end do
    wavelet(S_VELO,:)%bdry_uptodate = .false.
  end subroutine forward_wavelet_transform

  subroutine forward_scalar_transform (scaling, wavelet)
    ! Forward scalar wavelet transform
    implicit none
    type(Float_Field), dimension(:), target :: scaling, wavelet

    integer :: k, l, d

    do l = level_end-1, level_start-1, -1
       call update_vector_bdry (scaling, l+1, 61)

       ! Compute scalar wavelet coefficients
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d (Compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             nullify (scalar, wc_s)
          end do
       end do
       call update_vector_bdry (wavelet, l+1, 62)

       ! Restrict scalars (sub-sample and lift) to coarser grid
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d (Restrict_scalar, grid(d), l, z_null, 0, 1) ! +1 to include poles
             nullify (scalar, wc_s)
          end do
       end do
       scaling%bdry_uptodate = .false.
       wavelet%bdry_uptodate = .false.
    end do
  end subroutine forward_scalar_transform

  subroutine inverse_wavelet_transform (wavelet, scaling, l_start0)
    ! Inverse wavelet transform
    implicit none
    type(Float_Field), dimension(:,:), target :: scaling, wavelet
    integer, optional                         :: l_start0

    integer :: l, d, k, l_start, v

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_array_bdry1 (wavelet, level_start, level_end, 4)
    call update_array_bdry1 (scaling, l_start,     level_end, 5)

    scaling%bdry_uptodate = .false.

    do l = l_start, level_end-1
       ! Prolong scalars to finer nodes existing at coarser grid (undo lifting)
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                scalar => scaling(v,k)%data(d)%elts
                wc_s   => wavelet(v,k)%data(d)%elts
                call apply_interscale_d2 (Prolong_scalar, grid(d), l, z_null, 0, 1) 
                nullify (scalar, wc_s)
             end do
          end do
       end do

       if (l > l_start) call update_vector_bdry__finish (scaling(S_VELO,:), l) ! for next outer velocity

       call update_array_bdry__start (scaling(scalars(1):scalars(2),:), l+1)

       ! Prolong outer velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d2 (Reconstruct_outer_velo, grid(d), l, z_null, 0, 1) 
             call apply_to_penta_d (Reconstruct_velo_penta, grid(d), l, z_null)
             nullify (velo, wc_u)
          end do
       end do

       call update_array_bdry__finish (scaling(scalars(1):scalars(2),:), l+1)
       call update_vector_bdry__start (scaling(S_VELO,:), l+1)

       ! Prolong scalars at finer nodes not existing at coarser grid (interpolate and add wavelet coefficients)
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                scalar => scaling(v,k)%data(d)%elts
                wc_s   => wavelet(v,k)%data(d)%elts
                call apply_interscale_d (Reconstruct_scalar, grid(d), l, z_null, 0, 0)
                nullify (scalar, wc_s)
             end do
          end do
       end do

       call update_vector_bdry__finish (scaling(S_VELO,:), l+1)

       ! Prolong inner velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d (Reconstruct_inner_velo, grid(d), l, z_null, 0, 0)
             nullify (velo, wc_u)
          end do
       end do

       if (l < level_end-1) call update_vector_bdry__start (scaling(S_VELO,:), l+1) ! for next outer velocity

       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_wavelet_transform

  subroutine inverse_scalar_transform (wavelet, scaling, l_start0)
    ! Inverse scalar wavelet transform
    implicit none
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0

    integer :: l, d, k, l_start

    if (present(l_start0)) then
       l_start = l_start0
    else
       l_start = level_start
    end if

    call update_vector_bdry1 (wavelet, level_start, level_end, 64)
    call update_vector_bdry1 (scaling, l_start,     level_end, 65)

    scaling%bdry_uptodate = .false.
    
    do l = l_start, level_end-1
       ! Prolong scalar to finer nodes existing at coarser grid (undo lifting)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d2 (Prolong_scalar, grid(d), l, z_null, 0, 1) ! needs wc
             nullify (scalar, wc_s)
          end do
       end do
       call update_vector_bdry (scaling, l+1, 66)

       ! Prolong scalars at finer nodes not existing at coarser grid (interpolate and add wavelet coefficients)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d (Reconstruct_scalar, grid(d), l, z_null, 0, 0)
             nullify (scalar, wc_s)
          end do
       end do
       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_scalar_transform

  subroutine inverse_velo_transform (wavelet, scaling, l_start0)
    ! Inverse velocity wavelet transform
    implicit none
    type(Float_Field), dimension(:), target :: scaling, wavelet
    integer, optional                       :: l_start0

    integer :: l, d, k, l_start

    if (.not. present(l_start0)) then
       l_start = level_start
    else
       l_start = l_start0
    end if

    call update_vector_bdry1 (wavelet, level_start, level_end, 4)
    call update_vector_bdry1 (scaling, l_start,     level_end, 5)

    scaling%bdry_uptodate = .false.

    do l = l_start, level_end-1
       if (l > l_start) call update_vector_bdry__finish (scaling, l) ! for next outer velocity

       ! Prolong outer velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             velo => scaling(k)%data(d)%elts
             wc_u => wavelet(k)%data(d)%elts
             call apply_interscale_d2 (Reconstruct_outer_velo, grid(d), l, z_null, 0, 1) ! needs val
             call apply_to_penta_d (Reconstruct_velo_penta, grid(d), l, z_null)
             nullify (velo, wc_u)
          end do
       end do

       call update_vector_bdry (scaling, l+1, 33)

       ! Prolong inner velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             velo => scaling(k)%data(d)%elts
             wc_u => wavelet(k)%data(d)%elts
             call apply_interscale_d (Reconstruct_inner_velo, grid(d), l, z_null, 0, 0)
             nullify (velo, wc_u)
          end do
       end do

       if (l < level_end-1) call update_vector_bdry__start (scaling, l+1) ! for next outer velocity

       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_velo_transform

  subroutine Restrict_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Restrict both scalar and potential temperature
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) == 0) return

    id_par = idx (i_par, j_par, offs_par, dims_par)

    scalar(id_par+1) = restrict_s ()
  contains
    real(8) function restrict_s ()
      ! Restriction operator at nodes: sub-sample and lift
      integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE

      idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
      idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
      idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
      id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
      idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
      idW   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
      idNW  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
      idS2W = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)
      idSW  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
      idS   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
      id2SW = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
      idSE  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

      restrict_s = scalar(id_chd+1) + &
           (wc_s(idE+1)   * dom%overl_areas%elts(idE+1)%a(1)   + &
           wc_s(idNE+1)  * dom%overl_areas%elts(idNE+1)%a(2)  + &
           wc_s(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
           wc_s(id2NE+1) * dom%overl_areas%elts(id2NE+1)%a(4) + &
           wc_s(idN+1)   * dom%overl_areas%elts(idN+1)%a(1)   + &
           wc_s(idW+1)   * dom%overl_areas%elts(idW+1)%a(2)   + &
           wc_s(idNW+1)  * dom%overl_areas%elts(idNW+1)%a(3)  + &
           wc_s(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
           wc_s(idSW+1)  * dom%overl_areas%elts(idSW+1)%a(1)  + &
           wc_s(idS+1)   * dom%overl_areas%elts(idS+1)%a(2)   + &
           wc_s(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
           wc_s(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4)) * dom%areas%elts(id_par+1)%hex_inv
    end function restrict_s
  end subroutine Restrict_scalar

  subroutine Prolong_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars to fine points existing at coarse scale by undoing lifting
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd, idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) == FROZEN) return ! FROZEN mask -> do not overide with wrong value

    idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idW   = idx (i_chd-1, j_chd,   offs_chd, dims_chd)
    idNW  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idS2W = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)
    idSW  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS   = idx (i_chd,   j_chd-1, offs_chd, dims_chd)
    id2SW = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idSE  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    if (id2NE+1 <= ubound (wc_s, dim=1)) then
       scalar(id_chd+1) = scalar(id_par+1) - ( &
            wc_s(idE+1)   * dom%overl_areas%elts(idE+1)%a(1) + &
            wc_s(idNE+1)  * dom%overl_areas%elts(idNE+1)%a(2) + &
            wc_s(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
            wc_s(id2NE+1) * dom%overl_areas%elts(id2NE+1)%a(4) + &
            wc_s(idN+1)   * dom%overl_areas%elts(idN+1)%a(1) + &
            wc_s(idW+1)   * dom%overl_areas%elts(idW+1)%a(2) + &
            wc_s(idNW+1)  * dom%overl_areas%elts(idNW+1)%a(3) + &
            wc_s(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
            wc_s(idSW+1)  * dom%overl_areas%elts(idSW+1)%a(1) + &
            wc_s(idS+1)   * dom%overl_areas%elts(idS+1)%a(2) + &
            wc_s(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
            wc_s(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4) &
            ) * dom%areas%elts(id_par+1)%hex_inv
    else
       scalar(id_chd+1) = scalar(id_par+1) - ( &
            wc_s(idE+1)   * dom%overl_areas%elts(idE+1)%a(1) + &
            wc_s(idNE+1)  * dom%overl_areas%elts(idNE+1)%a(2) + &
            wc_s(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
            wc_s(idN+1)   * dom%overl_areas%elts(idN+1)%a(1) + &
            wc_s(idW+1)   * dom%overl_areas%elts(idW+1)%a(2) + &
            wc_s(idNW+1)  * dom%overl_areas%elts(idNW+1)%a(3) + &
            wc_s(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
            wc_s(idSW+1)  * dom%overl_areas%elts(idSW+1)%a(1) + &
            wc_s(idS+1)   * dom%overl_areas%elts(idS+1)%a(2) + &
            wc_s(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
            wc_s(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4) &
            ) * dom%areas%elts(id_par+1)%hex_inv
    end if
  end subroutine Prolong_scalar

  subroutine Reconstruct_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars to fine nodes not existing at coarse scale by interpolating and adding wavelet coefficients
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)

    ! Interpolate scalars and add wavelets to reconstruct values at fine scale
    scalar(idE_chd+1)  = Interp_node (dom, idE_chd,  id_chd,    id2E_chd, id2NE_chd, id2S_chd)  + wc_s(idE_chd+1)
    scalar(idNE_chd+1) = Interp_node (dom, idNE_chd, id2NE_chd, id_chd,   id2E_chd,  id2N_chd)  + wc_s(idNE_chd+1)
    scalar(idN_chd+1)  = Interp_node (dom, idN_chd,  id_chd,    id2N_chd, id2W_chd,  id2NE_chd) + wc_s(idN_chd+1)
  end subroutine Reconstruct_scalar

  subroutine Compute_scalar_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute wavelet coefficients for scalars
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd    = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)

    if (dom%mask_n%elts(idE_chd+1) >= ADJZONE) &
         wc_s(idE_chd+1) = scalar(idE_chd+1) - Interp_node (dom, idE_chd, id_chd, id2E_chd, id2NE_chd, id2S_chd)

    if (dom%mask_n%elts(idNE_chd+1) >= ADJZONE) &
         wc_s(idNE_chd+1) = scalar(idNE_chd+1) - Interp_node (dom, idNE_chd, id2NE_chd, id_chd, id2E_chd, id2N_chd)

    if (dom%mask_n%elts(idN_chd+1) >= ADJZONE) &
         wc_s(idN_chd+1) = scalar(idN_chd+1) - Interp_node (dom, idN_chd, id_chd, id2N_chd, id2W_chd, id2NE_chd)
  end subroutine Compute_scalar_wavelets

  subroutine Compute_velo_wavelets (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute velocity wavelet coefficients (except at pentagons)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: e, idN_chd, idE_chd, idNE_chd, id1, id2
    real(8)               :: u
    real(8), dimension(6) :: u_inner

    ! Velocity wavelets at 6 outer child edges
    do e = 1, EDGE
       id1 = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2 = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       if (dom%mask_e%elts(EDGE*id2+e) < ADJZONE) cycle

       u = Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)

       wc_u(EDGE*id2+e) = velo(EDGE*id2+e) - u
       wc_u(EDGE*id1+e) = - wc_u(EDGE*id2+e) ! to ensure Restriction(Prolongation) = Identity
    end do

    ! Velocity wavelets at 6 inner child edges
    u_inner = Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*idE_chd+UP+1)  >= ADJZONE) wc_u(EDGE*idE_chd+UP+1)  = velo(EDGE*idE_chd+UP+1)  - u_inner(1)
    if (dom%mask_e%elts(EDGE*idE_chd+DG+1)  >= ADJZONE) wc_u(EDGE*idE_chd+DG+1)  = velo(EDGE*idE_chd+DG+1)  - u_inner(2)
    if (dom%mask_e%elts(EDGE*idNE_chd+RT+1) >= ADJZONE) wc_u(EDGE*idNE_chd+RT+1) = velo(EDGE*idNE_chd+RT+1) - u_inner(3)
    if (dom%mask_e%elts(EDGE*idN_chd+RT+1)  >= ADJZONE) wc_u(EDGE*idN_chd+RT+1)  = velo(EDGE*idN_chd+RT+1)  - u_inner(4)
    if (dom%mask_e%elts(EDGE*idN_chd+DG+1)  >= ADJZONE) wc_u(EDGE*idN_chd+DG+1)  = velo(EDGE*idN_chd+DG+1)  - u_inner(5)
    if (dom%mask_e%elts(EDGE*idNE_chd+UP+1) >= ADJZONE) wc_u(EDGE*idNE_chd+UP+1) = velo(EDGE*idNE_chd+UP+1) - u_inner(6)
  end subroutine Compute_velo_wavelets
  
  subroutine Compute_velo_wavelets_penta (dom, p, c, offs, dims, zlev)
    ! Compute velocity wavelet coefficients at pentagons
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                        :: id_chd, idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(8)                        :: v
    real(8), dimension(2)          :: u

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then 
       ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
       ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
       id_chd  = idx (0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx (0, LAST,   offs_chd, dims_chd)

       v = -(Iu_Base_Wgt(8) + dble(dom%I_u_wgt%elts(idN_chd+1)%enc(8)))*( &
             velo(idx(0, PATCH_SIZE, offs, dims)*EDGE + UP+1) &
            + velo(idx(-1, PATCH_SIZE, offs, dims)*EDGE + RT+1))

       if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
          wc_u(EDGE*id_chd+UP+1)  = wc_u(EDGE*id_chd +UP+1) - v
          wc_u(EDGE*idN_chd+UP+1) = wc_u(EDGE*idN_chd+UP+1) + v
       end if
    else 
       if (c == IPLUSJMINUS) then
          ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
          ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
          id_chd = idx (LAST-1, 0, offs_chd, dims_chd)
          idE_chd = idx (LAST, 0, offs_chd, dims_chd)

          v = (Iu_Base_Wgt(7) + dble(dom%I_u_wgt%elts(idE_chd+1)%enc(7)))*( &
               velo(idx(PATCH_SIZE, 0, offs, dims)*EDGE+RT+1) +velo(idx(PATCH_SIZE,-1, offs, dims)*EDGE+UP+1))

          if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
             wc_u(EDGE*id_chd+RT+1)  = wc_u(EDGE*id_chd +RT+1) - v
             wc_u(EDGE*idE_chd+RT+1) = wc_u(EDGE*idE_chd+RT+1) + v
          end if
       end if
    end if

    if (.not. c == IJMINUS) return
    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case

    id_chd  = idx (0, 0, offs_chd, dims_chd)
    idN_chd = idx (0, 1, offs_chd, dims_chd)
    idE_chd = idx (1, 0, offs_chd, dims_chd)

    u = velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
       wc_u(EDGE*id_chd+UP+1)  = wc_u(EDGE*id_chd +UP+1) + u(1)
       wc_u(EDGE*idN_chd+UP+1) = wc_u(EDGE*idN_chd+UP+1) - u(1)
    end if
    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
       wc_u(EDGE*id_chd+RT+1)  = wc_u(EDGE*id_chd +RT+1) + u(2)
       wc_u(EDGE*idE_chd+RT+1) = wc_u(EDGE*idE_chd+RT+1) - u(2)
    end if
  end subroutine Compute_velo_wavelets_penta

  function velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)
    implicit none
    real(8), dimension(2)          :: velo_interp_penta_corr
    type(Domain)                   :: dom
    integer, dimension(N_BDRY+1)   :: offs, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims, dims_chd

    integer :: i, j, i_chd, j_chd

    i = 0
    j = 0
    i_chd = 0
    j_chd = 0

    velo_interp_penta_corr = (/  &
         (Iu_Base_Wgt(9) &
         + dble (dom%I_u_wgt%elts(idx__fast(i_chd+end_pt(1,2,UP+1), j_chd+end_pt(2,2,UP+1), offs_chd(1))+1)%enc(9))) * &
         ((-velo(idx(0, -1, offs, dims)*EDGE+UP+1) - (-velo(idx(-1, -1, offs, dims)*EDGE+1))) - &
          (velo(ed_idx(i+end_pt(1,1,UP+1), j+end_pt(2,1,UP+1), hex_sides(:,hex_s_offs(UP+1)+0+1), offs, dims)+1) - &
           velo(ed_idx(i+opp_no(1,2,UP+1), j+opp_no(2,2,UP+1), hex_sides(:,hex_s_offs(UP+1)+1+1), offs, dims)+1))), &
         (Iu_Base_Wgt(6) &
         + dble (dom%I_u_wgt%elts(idx__fast(i_chd + end_pt(1,2,RT+1), j_chd + end_pt(2,2,RT+1), offs_chd(1))+1)%enc(6))) * &
         (velo(idx(-1, -1, offs, dims)*EDGE+1) + velo(idx(-1, 0, offs, dims)*EDGE + RT+1) &
         - (velo(ed_idx(i+opp_no(1,1,RT+1), j+opp_no(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+1+1), offs, dims)+1) &
         -  velo(ed_idx(i+end_pt(1,1,RT+1), j+end_pt(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+2+1), offs, dims)+1))) /)
  end function velo_interp_penta_corr

  subroutine Reconstruct_outer_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at 6 outer fine edges by interpolating and adding wavelet coefficients
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: e, id_chd, id_par, id1, id2
    real(8) :: u
    
    id_par = idx (i_par, j_par, offs_par, dims_par)

    do e = 1, EDGE
       id1    = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2    = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       u = Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)
       velo(EDGE*id2+e) = u + wc_u(EDGE*id2+e)
       
       velo(EDGE*id1+e) = 2d0 * velo(EDGE*id_par+e) - velo(EDGE*id2+e)  ! to ensure Restriction(Prolongation) = Identity
    end do
  end subroutine Reconstruct_outer_velo

  subroutine Reconstruct_inner_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at 6 inner fine edges by interpolating and adding wavelet coefficients
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_chd, idN_chd, idE_chd, idNE_chd
    real(8), dimension(6) :: u_inner

    id_chd   = idx (i_chd, j_chd, offs_chd, dims_chd)

    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    
    u_inner = Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

    velo(EDGE*idE_chd +UP+1) = u_inner(1) + wc_u(EDGE*idE_chd +UP+1)
    velo(EDGE*idE_chd +DG+1) = u_inner(2) + wc_u(EDGE*idE_chd +DG+1)
    velo(EDGE*idNE_chd+RT+1) = u_inner(3) + wc_u(EDGE*idNE_chd+RT+1)
    velo(EDGE*idN_chd +RT+1) = u_inner(4) + wc_u(EDGE*idN_chd +RT+1)
    velo(EDGE*idN_chd +DG+1) = u_inner(5) + wc_u(EDGE*idN_chd +DG+1)
    velo(EDGE*idNE_chd+UP+1) = u_inner(6) + wc_u(EDGE*idNE_chd+UP+1)
  end subroutine Reconstruct_inner_velo

  function Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)
    ! Interpolate inner velocities to fine edges
    implicit none
    real(8), dimension(6)          :: Interp_inner_velo
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer                       :: id, id_par, id1_par, id2_par, t, idN, idE, idNE, idN2E, id2NE, idN2, idE2
    real(8), dimension(LORT:UPLT) :: curl_u

    id_par = idx(i_par, j_par, offs_par, dims_par)

    curl_u = 0d0
    do t = LORT, UPLT
       id1_par = idx (i_par-t+1, j_par,   offs_par, dims_par)
       id2_par = idx (i_par,     j_par+t, offs_par, dims_par)
       curl_u(t) = (velo(EDGE*id_par +DG+1)*dom%len%elts(EDGE*id_par+DG+1) &
                  + velo(EDGE*id1_par+UP+1)*dom%len%elts(EDGE*id1_par+UP+1) &
                  + velo(EDGE*id2_par+RT+1)*dom%len%elts(EDGE*id2_par+RT+1)) / dom%triarea%elts(TRIAG*id_par+t+1)
    end do

    id = idx (i_chd, j_chd, offs_chd, dims_chd)
    
    idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN2  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    idE2  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)

    Interp_inner_velo(1) = (dom%triarea%elts(TRIAG*id+LORT+1) * curl_u(LORT) &
         - velo(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1) - velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)) &
         / dom%len%elts(EDGE*idE+UP+1)

    Interp_inner_velo(2) = (dom%triarea%elts(TRIAG*idE+LORT+1) * curl_u(LORT) &
         - velo(EDGE*idE+RT+1) * dom%len%elts(EDGE*idE+RT+1) - velo(EDGE*idE2+UP+1) * dom%len%elts(EDGE*idE2+UP+1)) &
         / dom%len%elts(EDGE*idE+DG+1)

    Interp_inner_velo(3) = (dom%triarea%elts(TRIAG*idNE+LORT+1) * curl_u(LORT) &
         - velo(EDGE*idNE+DG+1) * dom%len%elts(EDGE*idNE+DG+1) - velo(EDGE*idN2E+UP+1) * dom%len%elts(EDGE*idN2E+UP+1)) &
         / dom%len%elts(EDGE*idNE+RT+1)

    Interp_inner_velo(4) = (dom%triarea%elts(TRIAG*id+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1) - velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)) &
         / dom%len%elts(EDGE*idN+RT+1)

    Interp_inner_velo(5) = (dom%triarea%elts(TRIAG*idN+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*idN+UP+1) * dom%len%elts(EDGE*idN+UP+1) - velo(EDGE*idN2+RT+1) * dom%len%elts(EDGE*idN2+RT+1)) &
         / dom%len%elts(EDGE*idN+DG+1)

    Interp_inner_velo(6) = (dom%triarea%elts(TRIAG*idNE+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*idNE+DG+1) * dom%len%elts(EDGE*idNE+DG+1) - velo(EDGE*id2NE+RT+1) * dom%len%elts(EDGE*id2NE+RT+1)) &
         / dom%len%elts(EDGE*idNE+UP+1)
  end function Interp_inner_velo

  real(8) function Interp_outer_velo (dom, i_par, j_par, ide, e, offs, dims)
    ! Interpolate outer velocities to fine edge EDGE*ide + e
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, ide, e
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    real(8), dimension(9) :: wgt

    wgt = Iu_Base_Wgt + dble (dom%I_u_wgt%elts(ide+1)%enc)

    Interp_outer_velo = sum (wgt * (/ &
         velo(EDGE * idx(i_par, j_par, offs, dims) + e), &
         velo(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+2+1), offs, dims) + 1), &
         velo(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+3+1), offs, dims) + 1), &
         velo(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+5+1), offs, dims) + 1), &
         velo(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+0+1), offs, dims) + 1), &
         
         velo(ed_idx(i_par + opp_no(1,1,e), j_par + opp_no(2,1,e), hex_sides(:,hex_s_offs(e)+1+1), offs, dims) + 1) - &
         velo(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+2+1), offs, dims) + 1), &
         
         velo(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+3+1), offs, dims) + 1) - &
         velo(ed_idx(i_par + opp_no(1,1,e), j_par + opp_no(2,1,e), hex_sides(:,hex_s_offs(e)+4+1), offs, dims) + 1), &
         
         velo(ed_idx(i_par + opp_no(1,2,e), j_par + opp_no(2,2,e), hex_sides(:,hex_s_offs(e)+4+1), offs, dims) + 1) - &
         velo(ed_idx(i_par + end_pt(1,2,e), j_par + end_pt(2,2,e), hex_sides(:,hex_s_offs(e)+5+1), offs, dims) + 1), &
         
         velo(ed_idx(i_par + end_pt(1,1,e), j_par + end_pt(2,1,e), hex_sides(:,hex_s_offs(e)+0+1), offs, dims) + 1) - &
         velo(ed_idx(i_par + opp_no(1,2,e), j_par + opp_no(2,2,e), hex_sides(:,hex_s_offs(e)+1+1), offs, dims) + 1) /))
  end function Interp_outer_velo

  subroutine Reconstruct_velo_penta (dom, p, c, offs, dims, z_null)
    ! Interpolate velocity to fine edges at pentagons
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, c, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                        :: id_chd,  idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(8)                        :: v
    real(8), dimension(2)          :: u

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then
       id_chd  = idx (0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx (0, LAST,   offs_chd, dims_chd)

       v = (Iu_Base_Wgt(8) + dble (dom%I_u_wgt%elts(idN_chd+1)%enc(8)))*( &
              velo(EDGE*idx( 0, PATCH_SIZE, offs, dims)+UP+1) &
            + velo(EDGE*idx(-1, PATCH_SIZE, offs, dims)+RT+1))

       velo(EDGE*id_chd+UP+1)  = velo(EDGE*id_chd +UP+1) - v
       velo(EDGE*idN_chd+UP+1) = velo(EDGE*idN_chd+UP+1) + v
    else
       if (c == IPLUSJMINUS) then
          id_chd  = idx (LAST-1, 0, offs_chd, dims_chd)
          idE_chd = idx (LAST,   0, offs_chd, dims_chd)

          v = -(Iu_Base_Wgt(7) + dble(dom%I_u_wgt%elts(idE_chd+1)%enc(7)))*( &
                 velo(EDGE*idx(PATCH_SIZE,  0, offs, dims)+RT+1) &
               + velo(EDGE*idx(PATCH_SIZE, -1, offs, dims)+UP+1))

          velo(EDGE*id_chd +RT+1) = velo(EDGE*id_chd +RT+1) - v
          velo(EDGE*idE_chd+RT+1) = velo(EDGE*idE_chd+RT+1) + v
       end if
    end if

    if (.not. c == IJMINUS) return

    id_chd  = idx (0, 0, offs_chd, dims_chd)
    idN_chd = idx (0, 1, offs_chd, dims_chd)
    idE_chd = idx (1, 0, offs_chd, dims_chd)

    u = velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)

    velo(EDGE*id_chd +UP+1) = velo(EDGE*id_chd +UP+1) - u(1)
    velo(EDGE*idN_chd+UP+1) = velo(EDGE*idN_chd+UP+1) + u(1)
    velo(EDGE*id_chd +RT+1) = velo(EDGE*id_chd +RT+1) - u(2)
    velo(EDGE*idE_chd+RT+1) = velo(EDGE*idE_chd+RT+1) + u(2)
  end subroutine Reconstruct_velo_penta

  subroutine Restrict_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Circulation conserving velocity  restriction (note that fine edges are exact bisections of coarse edges)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par,  i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id_par

    id_par   = idx (i_par, j_par, offs_par, dims_par)
    id_chd   = idx (i_chd, j_chd, offs_chd, dims_chd)
    
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) > 0) velo(EDGE*id_par+RT+1) = 0.5d0 * (velo(EDGE*id_chd+RT+1)   + velo(EDGE*idE_chd+RT+1))
    if (dom%mask_e%elts(EDGE*id_chd+DG+1) > 0) velo(EDGE*id_par+DG+1) = 0.5d0 * (velo(EDGE*idNE_chd+DG+1) + velo(EDGE*id_chd+DG+1))
    if (dom%mask_e%elts(EDGE*id_chd+UP+1) > 0) velo(EDGE*id_par+UP+1) = 0.5d0 * (velo(EDGE*id_chd+UP+1)   + velo(EDGE*idN_chd+UP+1))
  end subroutine Restrict_velo

  subroutine basic_F_restr_wgt (dom, i_par, j_par, e, offs_par, dims_par, i0, j0, offs, dims, typ)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, e, i0, j0
    integer, dimension(N_BDRY+1)   :: offs_par
    integer, dimension(2,N_BDRY+1) :: dims_par
    integer, dimension(2)          :: typ

    integer                        :: i, j, k

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer, dimension(3,16)       :: ije

    integer, dimension(2)          :: ij_nbp_mp, ij_nbp_pp, ij_nbp_pm, ij_nbp_mm
    integer, dimension(3)          :: ije_lcsd
    integer, dimension(4)          :: id_enc
    real(8), dimension(6)          :: wgt

    if (e == UP) then
       id_enc = (/idx(i0-2, j0, offs, dims), idx(i0-2,j0+1, offs, dims), idx(i0+1, j0, offs, dims), idx(i0+1, j0+1, offs, dims)/)
       i = i0
       j = j0+1
    elseif (e == RT) then
       id_enc = (/idx(i0, j0, offs, dims), idx(i0, j0+1, offs, dims), idx(i0+1, j0-2, offs, dims), idx(i0+1, j0-1, offs, dims)/)
       i = i0+1
       j = j0
    else
       write(0,'(A)') 'Error 447: R_F_wgts'
       stop
    end if

    ije(:,UMZ+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+1+1)
    ije(:,UPZ+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+4+1)
    ije(:,WMP+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+0+1)
    ije(:,VPP+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+5+1)
    ije(:,WPM+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+3+1)
    ije(:,VMM+1) = (/i, j, 0/) + hex_sides(:,hex_s_offs(e+1)+2+1)

    ij_nbp_mp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1)+0+1)
    ij_nbp_pp = (/i, j/) + nghb_pt(:,hex_s_offs(e+1)+5+1)
    ij_nbp_pm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1)+3+1)
    ij_nbp_mm = (/i, j/) + nghb_pt(:,hex_s_offs(e+1)+2+1)

    ije(:,VMP +1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides (:,hex_s_offs(e+1)+4-2+1)
    ije(:,VMPP+1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides (:,hex_s_offs(e+1)+1+4+1)
    ije(:,UZP +1) = (/ij_nbp_mp(1), ij_nbp_mp(2), 0/) + hex_sides (:,hex_s_offs(e+1)+0+4+1)
    ije(:,WPPP+1) = (/ij_nbp_pp(1), ij_nbp_pp(2), 0/) + hex_sides (:,hex_s_offs(e+1)+4-4+1)
    ije(:,WPP +1) = (/ij_nbp_pp(1), ij_nbp_pp(2), 0/) + hex_sides (:,hex_s_offs(e+1)+1+2+1)
    ije(:,VPM +1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides (:,hex_s_offs(e+1)+1+4+1)
    ije(:,VPMM+1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides (:,hex_s_offs(e+1)+4-2+1)
    ije(:,UZM +1) = (/ij_nbp_pm(1), ij_nbp_pm(2), 0/) + hex_sides (:,hex_s_offs(e+1)+3-2+1)
    ije(:,WMMM+1) = (/ij_nbp_mm(1), ij_nbp_mm(2), 0/) + hex_sides (:,hex_s_offs(e+1)+1+2+1)
    ije(:,WMM +1) = (/ij_nbp_mm(1), ij_nbp_mm(2), 0/) + hex_sides (:,hex_s_offs(e+1)+4-4+1)

    k = 1
    if (dist(dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(ije(1,UZP+1),ije(2,UZP+1), &
         adj_tri(:,-k+2,ije(3,UZP+1)+1),offs,dims)+1)) < eps()) then ! COINCIDE
       dom%R_F_wgt%elts(id_enc(1)+1)%enc = 0d0
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = 0d0
    else

       if (typ(k+1) == OUTER1) then
          ije_lcsd = ije(:,VPP+1)
       else if (typ(k+1) == OUTER2) then
          ije_lcsd = ije(:,WMP+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZP+1)
       end if

       wgt = interp_F_wgts (e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,k+1,e+1), offs_par, dims_par)+1), &
            ije, (/VPP, WMP, UZP, UPZ, WPP, VMP, UMZ, WPPP, VMPP/))

       if (e == RT) then
          wgt = (/wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)/)
       else if (e == UP) then
          wgt = (/wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)/)
       end if

       dom%R_F_wgt%elts(id_enc(1)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = wgt(4:6)
    end if

    k = 0
    if (dist(dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(ije(1,UZM+1),ije(2,UZM+1), adj_tri(:,-k+2,ije(3,UZM+1)+1),offs,dims)+1)) < eps()) then ! Coincide
       dom%R_F_wgt%elts(id_enc(3)+1)%enc = 0d0
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = 0d0
    else

       if (typ(k+1) == OUTER1) then
          ije_lcsd = ije(:,VMM+1)
       else if (typ(k+1) == OUTER2) then
          ije_lcsd = ije(:,WPM+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZM+1)
       end if

       wgt = interp_F_wgts (e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,k+1,e+1), offs_par, dims_par)+1), &
            ije, (/VMM, WPM, UZM, UMZ, WMM, VPM, UPZ, WMMM, VPMM/))

       if (e == UP) then
          wgt = (/wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)/)
       else if (e == RT) then
          wgt = (/wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)/)
       end if

       dom%R_F_wgt%elts(id_enc(3)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = wgt(4:6)
    end if
  contains
    function interp_F_wgts (e, k1, ije_lcsd, endpt_o, ije, stencil)
      implicit none
      real(8), dimension(6)    :: interp_F_wgts
      integer                  :: e, k1
      integer, dimension(3,16) :: ije
      integer, dimension(9)    :: stencil
      type(Coord)              :: endpt_o

      integer                 :: id_tri, info
      integer, dimension(3)   :: ije_lcsd
      integer, dimension(6)   :: ipiv
      real(8), dimension(6,6) :: G
      real(8), dimension(6)   :: b
      type(Coord)             :: endpt, x, y

      id_tri = tri_idx (ije_lcsd(1), ije_lcsd(2), adj_tri(:,-k1+2,ije_lcsd(3)+1), offs, dims)

      call local_coord (dom%ccentre%elts(id_tri+1), dom%ccentre%elts(id_tri+1), &
           dom%midpt%elts(ed_idx(0, 0, ije_lcsd, offs, dims)+1), x, y)

      G(:,1) = coords_to_row (ije(:,stencil(1)+1), x, y)
      G(:,2) = coords_to_row (ije(:,stencil(2)+1), x, y)
      G(:,3) = coords_to_row (ije(:,stencil(3)+1), x, y)
      G(:,4) = coords_to_row (ije(:,stencil(4)+1), x, y) - coords_to_row (ije(:,stencil(5)+1), x, y)
      G(:,5) = coords_to_row (ije(:,stencil(6)+1), x, y) - coords_to_row (ije(:,stencil(7)+1), x, y)
      G(:,6) = coords_to_row (ije(:,stencil(8)+1), x, y) - coords_to_row (ije(:,stencil(9)+1), x, y)

      endpt = endpt_o
      b = coords_to_row_perp ((/dom%ccentre%elts(id_tri+1), endpt/), x, y)
      ipiv = 0
      info = 0
      call dgesv (6, 1, G, 6, ipiv, b, 6, info)
      interp_F_wgts = b
    end function interp_F_wgts

    function coords_to_row_perp (coords, x, y)
      implicit none
      real(8), dimension(6)     :: coords_to_row_perp
      type(Coord), dimension(2) :: coords
      type(Coord)               :: x, y

      type(Coord) :: dirvec, midpt

      midpt = mid_pt (coords(1), coords(2))
      dirvec = cross (vector(coords(1), coords(2)), midpt)
      coords_to_row_perp = coords_to_rowd (midpt, dirvec, x, y)*dist(coords(1), coords(2))
    end function coords_to_row_perp

    function coords_to_row (ije0, x, y)
      implicit none
      real(8), dimension(6) :: coords_to_row
      integer, dimension(3) :: ije0
      type(Coord)           :: x, y

      integer     :: i, j, e
      type(Coord) :: endpt1, endpt2, midpt
      real(8)     :: pedlen

      i = ije0(1)
      j = ije0(2)
      e = ije0(3)

      endpt1 = dom%node%elts(idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), offs, dims)+1)
      endpt2 = dom%node%elts(idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), offs, dims)+1)

      pedlen = dist(dom%ccentre%elts(tri_idx(i, j, adj_tri(:,1,e+1), offs, dims)+1), &
          dom%ccentre%elts(tri_idx(i, j, adj_tri(:,2,e+1), offs, dims)+1))

      midpt = mid_pt (endpt1, endpt2)

      coords_to_row = coords_to_rowd (midpt, vector(endpt1, endpt2), x, y)*pedlen
    end function coords_to_row
  end subroutine basic_F_restr_wgt

  subroutine init_wavelets
    implicit none
    integer :: d, i, k, num, v

    call init_Float_Field (wav_topography, AT_NODE)
    
    do k = 1, zmax
       do v = 1, N_VARIABLE
          call init_Float_Field (wav_coeff(v,k), POSIT(v))
       end do
    end do
    
    if (vert_diffuse) then
       do k = 1, zlevels
          call init_Float_Field (wav_tke(k), AT_NODE)
       end do
    end if

    do d = 1, size(grid)
       num = grid(d)%node%length
       call init (grid(d)%overl_areas, num)
       call init (grid(d)%I_u_wgt, num)

       do i = 1, num
          call init_Iu_Wgt (grid(d)%I_u_wgt%elts(i), Iu_Base_Wgt)
       end do

       call init (grid(d)%R_F_wgt, num)

       do i = 1, num
          call init_RF_Wgt (grid(d)%R_F_wgt%elts(i), (/0d0, 0d0, 0d0/))
       end do

       call init (wav_topography%data(d), num)

       do k = 1, zmax
          do v = scalars(1), scalars(2)
             call init (wav_coeff(v,k)%data(d), num)
          end do
          call init (wav_coeff(S_VELO,k)%data(d), EDGE*num)
       end do
       
       if (vert_diffuse) then
          do k = 1, zlevels
             call init (wav_tke(k)%data(d), num)
          end do
       end if
    end do
  end subroutine init_wavelets

  subroutine get_overl_areas (dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par, offs_chd, dims_chd, e, area, typ)
    implicit none
    type(Domain)                   :: dom
    integer                        :: e, i_par, j_par, i_chd, j_chd
    integer, dimension(2)          :: typ
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd
    real(8), dimension(8)          :: area

    type(Coord), dimension(6)   :: hex
    type(Coord), dimension(3,2) :: tri
    type(Coord)                 :: inters_pt0, inters_pt1, pt
    integer                     :: i
    logical                     :: does_inters0, does_inters1, troubles

    area = 0d0
    typ = 0

    hex = (/ (dom%ccentre%elts(tri_idx(i_chd, j_chd, no_adj_tri(:,i + &
         hex_s_offs(-e+3)+1), offs_chd, dims_chd)+1), i = 0, 6-1) /)

    tri = reshape((/dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,4,e+1), &
         offs_par, dims_par)+1), dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,2,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,3,e+1), &
         offs_par, dims_par)+1), dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,2,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,1,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,1,e+1), offs_par, dims_par)+1)/), (/3, 2/))

    pt = dom%node%elts(idx(i_chd, j_chd, offs_chd, dims_chd)+1)

    area(1) = triarea(hex(6), hex(1), pt)
    area(2) = triarea(hex(3), hex(4), pt)

    do i = 1, 2
       call arc_inters (tri(1,i), tri(2,i), hex(3*i-2), hex(3*i-1), inters_pt0, does_inters0, troubles)
       call arc_inters (tri(3,i), tri(2,i), hex(3*i),   hex(3*i-1), inters_pt1, does_inters1, troubles)
       if (does_inters0 .and. does_inters1) then
          area(i+4) = triarea(inters_pt0, tri(2,i), hex(3*i-1))
          area(i+6) = triarea(tri(2,i), hex(3*i-1), inters_pt1)
          area(i+2) = area(i+4) + area(i+6)
          area(i) = area(i) + triarea(hex(3*i-2), inters_pt0, pt)     + triarea(inters_pt0, pt, tri(2,i))
          area(-i+3) = area(-i+3) + triarea(inters_pt1, hex(3*i), pt) + triarea(tri(2,i), pt, inters_pt1)
          typ(-i+3) = INSIDE
       else
          if (.not. does_inters0 .and. .not. does_inters1) then
             area(i+2) = 0d0
             call arc_inters (tri(2,1), tri(2,2), hex(3*i-2), hex(3*i-1), inters_pt0, does_inters0, troubles)
             call arc_inters (tri(2,2), tri(2,1), hex(3*i-1), hex(3*i),   inters_pt1, does_inters1, troubles)
             if (.not. does_inters0 .and. does_inters1) then
                area(i) = area(i) + triarea(hex(3*i-2), hex(3*i-1), pt) + triarea(hex(3*i-1), inters_pt1, pt)
                area(-i+3) = area(-i+3) + triarea(inters_pt1, hex(3*i), pt)
                typ(-i+3) = OUTER2
             else
                if (does_inters0 .and. .not. does_inters1) then
                   area(i) = area(i) + triarea(hex(3*i-2), inters_pt0, pt)
                   area(-i+3) = area(-i+3) + triarea(hex(3*i-1), hex(3*i), pt) + triarea(inters_pt0, hex(3*i-1), pt)
                   typ(-i+3) = OUTER1
                else
                   write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'A', does_inters0, does_inters1
                end if
             end if
          else
             write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'B', does_inters0, does_inters1
          end if
       end if
    end do
    return
  end subroutine get_overl_areas

  subroutine normalize2 (q, u, v)
    implicit none
    real(8), dimension(2) :: q
    real(8)               :: nrm, u, v

    nrm = sqrt (q(1)**2 + q(2)**2)
    u = q(1)/nrm
    v = q(2)/nrm
  end subroutine normalize2

  real(8) function Interp_node (dom, id, id1, id2, id3, id4)
    ! Interpolation at nodes
    type(Domain) :: dom
    integer      :: id, id1, id2, id3, id4

    Interp_node = ( &
         dom%overl_areas%elts(id+1)%a(1) * scalar(id1+1) + &
         dom%overl_areas%elts(id+1)%a(2) * scalar(id2+1) + &
         dom%overl_areas%elts(id+1)%a(3) * scalar(id3+1) + &
         dom%overl_areas%elts(id+1)%a(4) * scalar(id4+1) &
         ) * dom%areas%elts(id+1)%hex_inv
  end function Interp_node

  subroutine local_coord (midpt, endpt1, endpt2, x, y)
    ! Local coordinate (x,y) on a tangent plane
    implicit none
    type(Coord) :: midpt, endpt1, endpt2, x, y

    type(Coord) :: y0

    x = direction (endpt1, endpt2)
    y0 = cross (x, midpt)
    y = normalize_Coord (y0)
  end subroutine local_coord

  subroutine set_RF_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd,  p_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_chd, idN_chd, idE_chd, idNE_chd
    integer, dimension(2) :: typ
    real(8), dimension(8) :: area

    id_chd   = idx (i_chd,     j_chd,     offs_chd, dims_chd)
    idN_chd  = idx (i_chd,     j_chd+1, offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,     offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)

    call get_overl_areas (dom, i_par, j_par, i_chd+1, j_chd, offs_par, dims_par, offs_chd, dims_chd, RT, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idE_chd+1), area)
    call basic_F_restr_wgt (dom, i_par, j_par, RT, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)

    call get_overl_areas (dom, i_par, j_par, i_chd+1, j_chd+1, offs_par, dims_par, offs_chd, dims_chd, DG, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idNE_chd+1), area)

    call get_overl_areas (dom, i_par, j_par, i_chd, j_chd+1, offs_par, dims_par, offs_chd, dims_chd, UP, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idN_chd+1), area)
    call basic_F_restr_wgt (dom, i_par, j_par, UP, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)

    !call set_coarse_overlay ! May be needed for ifort compilers to avoid floating point error on coarsest level
  end subroutine set_RF_wgts

  subroutine set_coarse_overlay 
    ! Set overlay quantities on coarsest level
    integer :: d, p

    p = 2
    do d = 1, size(grid)
       call apply_onescale_to_patch (zero_overlay, grid(d), p-1, z_null, 0, 1)
    end do
  end subroutine set_coarse_overlay

  subroutine zero_overlay (dom, i, j, zlev, offs, dims)
    ! Sets overlay values to zero
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx(i, j, offs, dims)
    dom%overl_areas%elts(id+1)%a     = 0d0
    dom%overl_areas%elts(id+1)%split = 0d0
  end subroutine zero_overlay

  subroutine set_WT_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Sets local weights used in outer velocity interpolation when refining the grid
    implicit none
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd

    id_chd  = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd = idx (i_chd+1, j_chd,   offs_chd, dims_chd)

    dom%I_u_wgt%elts(idE_chd+1) = outer_velo_weights (dom, p_chd, i_par, j_par, RT, offs_par, dims_par)
    dom%I_u_wgt%elts(id_chd+1)  = outer_velo_weights (dom, p_chd, i_par, j_par, DG, offs_par, dims_par)
    dom%I_u_wgt%elts(idN_chd+1) = outer_velo_weights (dom, p_chd, i_par, j_par, UP, offs_par, dims_par)
  end subroutine set_WT_wgts

  type(Iu_Wgt) function outer_velo_weights (dom, p, i0, j0, e0, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i0, j0, e0
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    type(Coord)             :: x, y
    integer                 :: k, id, info
    integer, dimension(6)   :: ipiv
    real(8), dimension(9)   :: weights
    real(8), dimension(6)   :: b
    real(8), dimension(6,6) :: G

    call local_coord (&
         dom%midpt%elts(idx(i0, j0, offs, dims)*EDGE + e0 + 1), &
         dom%node%elts(idx(i0 + end_pt(1,1,e0+1), j0 + end_pt(2,1,e0+1), offs, dims) + 1), &
         dom%node%elts(idx(i0 + end_pt(1,2,e0+1), j0 + end_pt(2,2,e0+1), offs, dims) + 1), x, y)

    weights = 0d0

    do k = 1, 2
       id = idx (i0, j0, offs, dims)

       G = 0d0

       G(:,1) = coords_to_row (i0, j0, (/0, 0/), (/0, 0, e0/), e0)

       G(:,2) = coords_to_row (i0, j0, end_pt(:,-k+3,e0+1), hex_sides(:,hex_s_offs(e0+1)+2+3*k-3+1), e0)

       G(:,3) = coords_to_row (i0, j0, end_pt(:,k,e0+1),    hex_sides(:,(hex_s_offs(e0+1)+3)-(3*k-3)+1), e0)

       G(:,4) = coords_to_row (i0, j0, opp_no(:,k,e0+1), hex_sides(:,hex_s_offs(e0+1)+1+3*k-3+1), e0) - &
                coords_to_row (i0, j0, end_pt(:,k,e0+1), hex_sides(:,hex_s_offs(e0+1)+2+3*k-3+1), e0)

       G(:,5) = coords_to_row (i0, j0, end_pt(:,-k+3,e0+1), hex_sides(:,(hex_s_offs(e0+1)+3)-(3*k-3)+1), e0) - &
                coords_to_row (i0, j0, opp_no(:, k,  e0+1), hex_sides(:,(hex_s_offs(e0+1)+4)-(3*k-3)+1), e0)

       G(:,6) = coords_to_row (i0, j0, end_pt(:, k,  e0+1), hex_sides(:,(hex_s_offs(e0+1)+5)-(3*k-3)+1), e0) - &
                coords_to_row (i0, j0, end_pt(:,-k+3,e0+1), hex_sides(:, hex_s_offs(e0+1)+0+3*k-3+1),    e0)

       b = coords_to_rowd (mid_pt(dom%midpt%elts(EDGE*id+e0+1), dom%node%elts(idx2(i0, j0, end_pt(:,2,e0+1), offs, dims) + 1)), &
            vector (dom%node%elts(idx2(i0, j0, end_pt(:,1,e0+1), offs, dims)+1), &
            dom%node%elts(idx2(i0, j0, end_pt(:,2,e0+1), offs, dims)+1)), x, y)

       ipiv = 0
       info = 0

       call dgesv (6, 1, G, 6, ipiv, b, 6, info)

       weights(1)           = weights(1) + b(1)
       weights(2*k:2*k+1)   = weights(2*k:2*k+1) + b(2:3)
       weights(2*k+4:2*k+5) = b(4:5)
       weights(-2*k+6)      = weights(-2*k+6) + b(6)
       weights(-2*k+7)      = weights(-2*k+7) - b(6)
    end do

    outer_velo_weights = Iu_wgt (0.5d0 * weights - Iu_Base_Wgt)
  contains
    function coords_to_row (i00, j00, n_offs1, n_offs2, e00)
      implicit none
      real(8), dimension(6) :: coords_to_row
      integer               :: e00, i00, j00
      integer, dimension(2) :: n_offs1
      integer, dimension(3) :: n_offs2

      type(Coord) :: endpt1, endpt2
      integer     :: i, j, e

      i = i00 + n_offs1(1) + n_offs2(1)
      j = j00 + n_offs1(2) + n_offs2(2)

      e = n_offs2(3)

      endpt1 = get_coord (i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), e00)
      endpt2 = get_coord (i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), e00)

      coords_to_row = coords_to_rowd (mid_pt(endpt1, endpt2), vector(endpt1, endpt2), x, y)
    end function coords_to_row

    function get_coord (i, j, e)
      type(Coord) ::  get_coord
      integer     :: i, j, e

      if (i == -1) then
         if (j == -1 .and. is_penta (dom, p, IJMINUS-1)) then
            if (e == RT) then
               get_coord = dom%node%elts(nidx(LAST_BDRY, 0, IMINUS, offs, dims)+1)
               return
            else
               if (e == UP) then
                  get_coord = dom%node%elts(nidx(0, LAST_BDRY, JMINUS, offs, dims)+1)
                  return
               end if
            end if
         else
            if (j == PATCH_SIZE .and. is_penta(dom, p, IMINUSJPLUS-1)) then
               get_coord = dom%node%elts(nidx(0, 1, JPLUS, offs, dims)+1)
               return
            else
               get_coord = dom%node%elts(idx(i, j, offs, dims)+1)
               return
            end if
         end if
      else
         if (i == PATCH_SIZE .and. j == -1 .and. is_penta(dom, p, IPLUSJMINUS-1)) then
            get_coord = dom%node%elts(nidx(1, 0, IPLUS, offs, dims)+1)
            return
         else
            if (i == PATCH_SIZE+1 .and. j == PATCH_SIZE+1 .and. is_penta(dom, p, IJPLUS-1)) then
               get_coord = dom%node%elts(nidx(1, 0, IJPLUS, offs, dims)+1)
               return
            else
               get_coord = dom%node%elts(idx(i, j, offs, dims)+1)
               return
            end if
         end if
      end if
    end function get_coord
  end function outer_velo_weights

  subroutine check_m (dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par,  offs_chd, dims_chd)
    ! Check that scalar is indeed conserved by restriction
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par, idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE
    real(8) :: ratio

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    idE   = idx (i_chd+1, j_chd,     offs_chd, dims_chd)
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd)
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd)
    idN   = idx (i_chd,     j_chd+1, offs_chd, dims_chd)
    idW   = idx (i_chd-1, j_chd,     offs_chd, dims_chd)
    idNW  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd)
    idS2W = idx (i_chd-2, j_chd-1, offs_chd, dims_chd)
    idSW  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd)
    idS   = idx (i_chd,     j_chd-1, offs_chd, dims_chd)
    id2SW = idx (i_chd-1, j_chd-2, offs_chd, dims_chd)
    idSE  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd)

    ratio = (1/dom%areas%elts(id_chd+1)%hex_inv + &
         dom%overl_areas%elts(idE+1)%a(1) + &
         dom%overl_areas%elts(idNE+1)%a(2) + &
         dom%overl_areas%elts(idN2E+1)%a(3) + &
         dom%overl_areas%elts(id2NE+1)%a(4) + &
         dom%overl_areas%elts(idN+1)%a(1) + &
         dom%overl_areas%elts(idW+1)%a(2) + &
         dom%overl_areas%elts(idNW+1)%a(3) + &
         dom%overl_areas%elts(idS2W+1)%a(4) + &
         dom%overl_areas%elts(idSW+1)%a(1) + &
         dom%overl_areas%elts(idS+1)%a(2) + &
         dom%overl_areas%elts(id2SW+1)%a(3) + &
         dom%overl_areas%elts(idSE+1)%a(4))*dom%areas%elts(id_par+1)%hex_inv
  end subroutine check_m

  function coord2local(c, x, y)
    implicit none
    real(8), dimension(2) :: coord2local
    type(Coord)           :: c, x, y

    coord2local = (/inner(c, x), inner(c, y)/)
  end function coord2local

  subroutine init_wavelet_mod
    implicit none
    logical :: initialized = .false.

    if (initialized) return ! initialize only once
    call init_shared_mod
    call init_domain_mod

    Iu_Base_Wgt = (/16d0, -1d0, 1d0, 1d0, -1d0, -1d0, -1d0, 1d0, 1d0/) / 16d0
    initialized = .true.
  end subroutine init_wavelet_mod

  function coords_to_rowd (midpt, dirvec, x, y)
    implicit none
    real(8), dimension(6) :: coords_to_rowd
    type(Coord)           :: dirvec, x, y

    type(Coord)           :: midpt
    real(8)               :: u, v
    real(8), dimension(2) :: xy

    call normalize2 (coord2local(dirvec, x, y), u, v)
    xy = coord2local (midpt, x, y)
    coords_to_rowd = (/u, u*xy(1), u*xy(2), v, v*xy(1), v*xy(2)/)
  end function coords_to_rowd
end module wavelet_mod
