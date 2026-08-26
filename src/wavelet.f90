module wavelet_mod
  
  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, EDGE, Coord, scalars, hex_sides, hex_s_offs, nghb_pt, adj_tri, no_adj_tri, bfly_tri, &
       end_pt, opp_no, N_BDRY, ORIGIN, IJMINUS, IJPLUS, IMINUS, IMINUSJPLUS, IPLUS, IPLUSJMINUS, JMINUS, JPLUS, RT, DG, UP, &
       Z_NULL, INSIDE, OUTER1, OUTER2, N_VARIABLE, vert_diffuse, zmin, zmax, zlevels, AT_NODE, radius, &
       UPZ, UZM, UZP, VMM, VMP, VMPP, VPP, UMZ, VPM, VPMM, WMM, WMMM, WMP, WPM, WPP, WPPP, &
       S_VELO, TRIAG, FROZEN, LORT, UPLT, NONE, POSIT, level_start, level_end, eps


  use arch_mod,       only : abort_run
  use comm_mpi_mod,   only : update_bdry, update_bdry1, update_bdry__finish, update_bdry__start
  use domain_mod,     only : get_offs_Domain, init_Field
  use domain_ops_mod, only : apply_interscale_d, apply_interscale_d2, apply_onescale_to_patch, apply_to_penta_d
  use dyn_arrays,     only : init
  use geom_mod,       only : arc_intersect_test, cross, direction, dist, inner, mid_pt, normalize_Coord, vector, triarea
  use patch_mod,      only : init_Overl_Area, init_Iu_Wgt, init_RF_Wgt, Iu_wgt, LAST, LAST_BDRY, PATCH_SIZE
  use utils_mod,      only : zero_float

  use domain_mod, only : Domain, Float_Field, grid, velo, wav_coeff, wav_tke, wc_s, wc_u, scalar, idx, is_penta, &
       nidx, tri_idx, idx2, idx__fast, ed_idx

  implicit none

  
  private
  public :: forward_wavelet_transform, forward_scalar_transform
  public :: inverse_wavelet_transform, inverse_scalar_transform, inverse_velo_transform
  public :: inverse_velo_outer_level
  public :: Compute_scalar_wavelets, Compute_velo_wavelets, Compute_velo_wavelets_penta, Restrict_velo
  public :: Restrict_scalar, scalar_restriction, init_wavelets, set_RF_wgts, set_WT_wgts
  public :: check_m, Prolong_full_weighting

  
  real(dp), parameter :: Iu_Base_Wgt(9) = [16.0_dp, -1.0_dp, 1.0_dp, 1.0_dp, -1.0_dp, -1.0_dp, -1.0_dp, 1.0_dp, 1.0_dp] / 16.0_dp
  
  logical,  parameter :: lapack = .true. ! use lapack or local LU routine

  
  interface
     
     subroutine dgesv(n, nrhs, a, lda, ipiv, b, ldb, info)
       use kind_mod, only : dp
       implicit none

       integer, intent(in)    :: n
       integer, intent(in)    :: nrhs
       integer, intent(in)    :: lda
       integer, intent(in)    :: ldb

       real(dp), intent(inout) :: a(lda,*)
       integer,  intent(out)   :: ipiv(*)
       real(dp), intent(inout) :: b(ldb,*)

       integer, intent(out) :: info
     end subroutine dgesv
     
  end interface

  interface forward_scalar_transform
     procedure :: forward_scalar_transform_0, forward_scalar_transform_1
  end interface forward_scalar_transform

  interface inverse_scalar_transform
     procedure :: inverse_scalar_transform_0, inverse_scalar_transform_1
  end interface inverse_scalar_transform

  interface inverse_velo_transform
     procedure inverse_velo_transform_0, inverse_velo_transform_1
  end interface inverse_velo_transform

  
contains
  

  subroutine forward_wavelet_transform (scaling, wavelet, jmin_in, jmax_in)
    ! Forward wavelet transform
    implicit none

    type(Float_Field), intent(inout), target :: scaling(:,:), wavelet(:,:)
    integer, optional, intent(in)            :: jmin_in, jmax_in
    
    integer :: d, jmin, jmax, k, l, v
    
    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if
    
    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine .. abort"
       call abort_run
    end if

    call zero_float (wavelet)
    wavelet%bdry_uptodate = .false.
    
    do l = jmax-1, jmin-1, -1
       ! Compute scalar wavelet coefficients
       call update_bdry (scaling(scalars(1):scalars(2),:), l+1, 918)

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
       call update_bdry (wavelet(scalars(1):scalars(2),:), l+1, 919)

       ! Restrict scalars (sub-sample and lift) and velocity (average) from fine grid l+1 to coarse grid l
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

    call update_bdry (scaling(S_VELO,:), NONE, 920)

    ! Compute vector wavelet coefficients
    do l = jmax-1, jmin-1, -1
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
  

  subroutine forward_scalar_transform_0 (scaling, wavelet, jmin_in, jmax_in)
    ! Forward scalar wavelet transform
    implicit none

    type(Float_Field), intent(inout), target :: scaling, wavelet
    integer, optional, intent(in)            :: jmin_in, jmax_in

    integer :: d, jmin, jmax, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine .. abort"
       call abort_run
    end if

    call zero_float (wavelet)
    wavelet%bdry_uptodate = .false.

    do l = jmax-1, jmin-1, -1
       call update_bdry (scaling, l+1, 921)

       ! Compute scalar wavelet coefficients
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          wc_s   => wavelet%data(d)%elts
          call apply_interscale_d (Compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
          nullify (scalar, wc_s)
       end do
       call update_bdry (wavelet, l+1, 922)

       ! Restrict scalars (sub-sample and lift) from fine grid l+1 to coarse grid l
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          wc_s   => wavelet%data(d)%elts
          call apply_interscale_d (Restrict_scalar, grid(d), l, z_null, 0, 1) ! +1 to include poles
          nullify (scalar, wc_s)
       end do
       scaling%bdry_uptodate = .false.
       wavelet%bdry_uptodate = .false.
    end do
  end subroutine forward_scalar_transform_0
  

  subroutine forward_scalar_transform_1 (scaling, wavelet, jmin_in, jmax_in)
    ! Forward scalar wavelet transform
    implicit none

    type(Float_Field), intent(inout), target :: scaling(:), wavelet(:)
    integer, optional, intent(in)            :: jmin_in, jmax_in

    integer :: d, jmin, jmax, k, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine .. abort"
       call abort_run
    end if

    call zero_float (wavelet)
    wavelet%bdry_uptodate = .false.
    
    do l = jmax-1, jmin-1, -1
       call update_bdry (scaling, l+1, 923)

       ! Compute scalar wavelet coefficients
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d (Compute_scalar_wavelets, grid(d), l, z_null, 0, 0)
             nullify (scalar, wc_s)
          end do
       end do
       call update_bdry (wavelet, l+1, 924)

       ! Restrict scalars (sub-sample and lift) from fine grid l+1 to coarse grid l
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
  end subroutine forward_scalar_transform_1
  

  subroutine inverse_wavelet_transform (wavelet, scaling, jmin_in, jmax_in)
    ! Inverse wavelet transform
    implicit none
    integer, optional                         :: jmin_in, jmax_in
    type(Float_Field), dimension(:,:), target :: scaling, wavelet

    integer :: d, jmin, jmax, k, l, v

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine ... abort"
       call abort_run
    end if

    call update_bdry1 (wavelet, max (jmin, level_start), jmax, 802)
    call update_bdry1 (scaling, jmin,                    jmax, 803)

    scaling%bdry_uptodate = .false.
    do l = jmin, jmax-1
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

       if (l > jmin) call update_bdry__finish (scaling(S_VELO,:), l) ! for next outer velocity

       call update_bdry__start (scaling(scalars(1):scalars(2),:), l+1)

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

       call update_bdry__finish (scaling(scalars(1):scalars(2),:), l+1)
       call update_bdry__start (scaling(S_VELO,:),                 l+1)

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

       call update_bdry__finish (scaling(S_VELO,:), l+1)

       ! Prolong inner velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling,2)
          do d = 1, size(grid)
             velo => scaling(S_VELO,k)%data(d)%elts
             wc_u => wavelet(S_VELO,k)%data(d)%elts
             call apply_interscale_d (Reconstruct_inner_velo, grid(d), l, z_null, 0, 0)
             nullify (velo, wc_u)
          end do
       end do

       if (l < jmax-1) call update_bdry__start (scaling(S_VELO,:), l+1) ! for next outer velocity

       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_wavelet_transform
  

  subroutine inverse_velo_outer_level (wavelet,scaling,level)
    ! Reconstruct regular outer edges and their pentagon correction together.
    ! Their Domain topology/orientation transport remains indivisible here;
    ! Stage 136 takes ownership of the subsequent inner-edge reconstruction.

    implicit none

    type(Float_Field), intent(inout), target :: wavelet(:,:), scaling(:,:)
    integer, intent(in) :: level

    integer :: d
    integer :: k

    do k = 1,size(scaling,2)
       do d = 1,size(grid)
          velo => scaling(S_VELO,k)%data(d)%elts
          wc_u => wavelet(S_VELO,k)%data(d)%elts
          call apply_interscale_d2( &
               Reconstruct_outer_velo,grid(d),level,z_null,0,1)
          call apply_to_penta_d( &
               Reconstruct_velo_penta,grid(d),level,z_null)
          nullify(velo,wc_u)
       end do
       scaling(S_VELO,k)%bdry_uptodate = .false.
    end do

  end subroutine inverse_velo_outer_level


  subroutine inverse_scalar_transform_0 (wavelet, scaling, jmin_in, jmax_in)
    ! Inverse scalar wavelet transform
    implicit none

    integer, optional, intent(in)            :: jmin_in, jmax_in
    type(Float_Field), intent(inout), target :: scaling, wavelet

    integer :: d, jmin, jmax, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine ... abort"
       call abort_run
    end if

    call update_bdry1 (wavelet, max (jmin, level_start), jmax, 809)
    call update_bdry1 (scaling, jmin,                    jmax, 810)

    scaling%bdry_uptodate = .false.
    do l = jmin, jmax-1
       ! Prolong scalar to finer nodes existing at coarser grid (undo lifting)
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          wc_s   => wavelet%data(d)%elts
          call apply_interscale_d2 (Prolong_scalar, grid(d), l, z_null, 0, 1) ! needs wc
          nullify (scalar, wc_s)
       end do
       call update_bdry (scaling, l+1, 925)

       ! Prolong scalars at finer nodes not existing at coarser grid (interpolate and add wavelet coefficients)
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          wc_s   => wavelet%data(d)%elts
          call apply_interscale_d (Reconstruct_scalar, grid(d), l, z_null, 0, 0)
          nullify (scalar, wc_s)
       end do
       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_scalar_transform_0
  

  subroutine inverse_scalar_transform_1 (wavelet, scaling, jmin_in, jmax_in)
    ! Inverse scalar wavelet transform
    implicit none

    integer, optional, intent(in)            :: jmin_in, jmax_in
    type(Float_Field), intent(inout), target :: scaling(:), wavelet(:)

    integer :: d, jmin, jmax, k, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine ... abort"
       call abort_run
    end if

    call update_bdry1 (wavelet, max (jmin, level_start), jmax, 811)
    call update_bdry1 (scaling, jmin,                    jmax, 812)

    scaling%bdry_uptodate = .false.
    do l = jmin, jmax-1
       ! Prolong scalar to finer nodes existing at coarser grid (undo lifting)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             scalar => scaling(k)%data(d)%elts
             wc_s   => wavelet(k)%data(d)%elts
             call apply_interscale_d2 (Prolong_scalar, grid(d), l, z_null, 0, 1) ! needs wc
             nullify (scalar, wc_s)
          end do
       end do
       call update_bdry (scaling, l+1, 926)

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
  end subroutine inverse_scalar_transform_1
  
  
  subroutine inverse_velo_transform_0 (wavelet, scaling, jmin_in, jmax_in)
    ! Inverse velocity wavelet transform
    implicit none
    integer, optional         :: jmin_in, jmax_in
    type(Float_Field), target :: scaling, wavelet
    
    integer :: d, jmin, jmax, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine ... abort"
       call abort_run
    end if

    call update_bdry1 (wavelet, max (jmin, level_start), jmax, 813)
    call update_bdry1 (scaling, jmin,                    jmax, 814)

    scaling%bdry_uptodate = .false.
    do l = jmin, jmax-1
       if (l > jmin) call update_bdry__finish (scaling, l) ! for next outer velocity

       ! Prolong outer velocities at finer edges (interpolate and add wavelet coefficients)
       do d = 1, size(grid)
          velo => scaling%data(d)%elts
          wc_u => wavelet%data(d)%elts
          call apply_interscale_d2 (Reconstruct_outer_velo, grid(d), l, z_null, 0, 1) ! needs val
          call apply_to_penta_d (Reconstruct_velo_penta, grid(d), l, z_null)
          nullify (velo, wc_u)
       end do

       call update_bdry (scaling, l+1, 927)

       ! Prolong inner velocities at finer edges (interpolate and add wavelet coefficients)
       do d = 1, size(grid)
          velo => scaling%data(d)%elts
          wc_u => wavelet%data(d)%elts
          call apply_interscale_d (Reconstruct_inner_velo, grid(d), l, z_null, 0, 0)
          nullify (velo, wc_u)
       end do

       if (l < jmax-1) call update_bdry__start (scaling, l+1) ! for next outer velocity

       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_velo_transform_0
  
  
  subroutine inverse_velo_transform_1 (wavelet, scaling, jmin_in, jmax_in)
    ! Inverse velocity wavelet transform
    implicit none

    integer, optional, intent(in)            :: jmin_in, jmax_in
    type(Float_Field), intent(inout), target :: scaling(:), wavelet(:)

    integer :: d, jmin, jmax, k, l

    if (present(jmin_in)) then
       jmin = jmin_in
    else
       jmin = level_start
    end if

    if (present(jmax_in)) then
       jmax = jmax_in
    else
       jmax = level_end
    end if

    if (jmax < jmin) then
       write (6,'(a)') "ERROR: jmax < jmin in wavelet routine ... abort"
       call abort_run
    end if

    call update_bdry1 (wavelet, max (jmin, level_start), jmax, 817)
    call update_bdry1 (scaling, jmin,                    jmax, 818)

    scaling%bdry_uptodate = .false.
    do l = jmin, jmax-1
       if (l > jmin) call update_bdry__finish (scaling, l) ! for next outer velocity

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

       call update_bdry (scaling, l+1, 928)

       ! Prolong inner velocities at finer edges (interpolate and add wavelet coefficients)
       do k = 1, size(scaling)
          do d = 1, size(grid)
             velo => scaling(k)%data(d)%elts
             wc_u => wavelet(k)%data(d)%elts
             call apply_interscale_d (Reconstruct_inner_velo, grid(d), l, z_null, 0, 0)
             nullify (velo, wc_u)
          end do
       end do

       if (l < jmax-1) call update_bdry__start (scaling, l+1) ! for next outer velocity
       
       scaling%bdry_uptodate = .false.
    end do
  end subroutine inverse_velo_transform_1
  

  subroutine Restrict_scalar ( &
       dom, i_par, j_par, i_chd, j_chd, zlev, &
       offs_par, dims_par, offs_chd, dims_chd)

    ! Restrict scalar variables by subsampling and lifting.
    implicit none

    type(Domain), intent(inout) :: dom

    integer, intent(in) :: i_par, j_par
    integer, intent(in) :: i_chd, j_chd
    integer, intent(in) :: zlev

    integer, intent(in) :: offs_par(N_BDRY+1)
    integer, intent(in) :: offs_chd(N_BDRY+1)

    integer, intent(in) :: dims_par(2,N_BDRY+1)
    integer, intent(in) :: dims_chd(2,N_BDRY+1)

    integer :: id_chd, id_par

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) == 0) return

    id_par = idx(i_par, j_par, offs_par, dims_par)

    scalar(id_par+1) = restrict_s()

  contains

    function restrict_s() result(val)
      ! Restriction operator at nodes: subsample and lift.
      implicit none

      real(dp) :: val

      integer :: idE, idNE, idN2E, id2NE
      integer :: idN, idW, idNW, idS2W
      integer :: idSW, idS, id2SW, idSE

      idE = idx( &
           i_chd+1, j_chd, offs_chd, dims_chd)

      idNE = idx( &
           i_chd+1, j_chd+1, offs_chd, dims_chd)

      idN2E = idx( &
           i_chd+2, j_chd+1, offs_chd, dims_chd)

      id2NE = idx( &
           i_chd+1, j_chd+2, offs_chd, dims_chd)

      idN = idx( &
           i_chd, j_chd+1, offs_chd, dims_chd)

      idW = idx( &
           i_chd-1, j_chd, offs_chd, dims_chd)

      idNW = idx( &
           i_chd-1, j_chd+1, offs_chd, dims_chd)

      idS2W = idx( &
           i_chd-2, j_chd-1, offs_chd, dims_chd)

      idSW = idx( &
           i_chd-1, j_chd-1, offs_chd, dims_chd)

      idS = idx( &
           i_chd, j_chd-1, offs_chd, dims_chd)

      id2SW = idx( &
           i_chd-1, j_chd-2, offs_chd, dims_chd)

      idSE = idx( &
           i_chd+1, j_chd-1, offs_chd, dims_chd)

      val = scalar(id_chd+1) + ( &
           wc_s(idE+1)   * dom%overl_areas%elts(idE+1)%a(1)   + &
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
           wc_s(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4)  &
           ) * dom%areas%elts(id_par+1)%hex_inv
    end function restrict_s

  end subroutine Restrict_scalar
  

  subroutine Prolong_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars to fine points existing at coarse scale by undoing lifting
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

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

    scalar(id_chd+1) = scalar(id_par+1) - ( &
         wc_s(idE+1)   * dom%overl_areas%elts(idE  +1)%a(1) + &
         wc_s(idNE+1)  * dom%overl_areas%elts(idNE +1)%a(2) + &
         wc_s(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
         wc_s(idN+1)   * dom%overl_areas%elts(idN  +1)%a(1) + &
         wc_s(idW+1)   * dom%overl_areas%elts(idW  +1)%a(2) + &
         wc_s(idNW+1)  * dom%overl_areas%elts(idNW +1)%a(3) + &
         wc_s(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
         wc_s(idSW+1)  * dom%overl_areas%elts(idSW+ 1)%a(1) + &
         wc_s(idS+1)   * dom%overl_areas%elts(idS  +1)%a(2) + &
         wc_s(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
         wc_s(idSE+1)  * dom%overl_areas%elts(idSE +1)%a(4) &
         ) * dom%areas%elts(id_par+1)%hex_inv
  end subroutine Prolong_scalar
  

  subroutine Reconstruct_scalar (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars to fine nodes not existing at coarse scale by interpolating and adding wavelet coefficients
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

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

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
    
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

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer                :: e, idN_chd, idE_chd, idNE_chd, id1, id2
    real(dp)               :: u
    real(dp), dimension(6) :: u_inner

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
    ! Compute velocity wavelet coefficients at pentagons.
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, c, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id_chd, idE_chd, idN_chd, p_chd
    integer :: offs_chd(N_BDRY+1)
    integer :: dims_chd(2,N_BDRY+1)

    real(dp) :: v
    real(dp) :: u(2)

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then
       ! Parts 3 and 4 of hexagon IMINUSJPLUS, at the upper-left
       ! corner of the lozenge, combine to form a pentagon.
       !
       ! pedlen(EDGE*idW + RT + 1) = 0 in this case.
       id_chd  = idx(0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx(0, LAST,   offs_chd, dims_chd)

       v = -( &
            Iu_Base_Wgt(8) + &
            real(dom%I_u_wgt%elts(idN_chd+1)%enc(8), kind=dp) &
            ) * ( &
            velo(idx(0,  PATCH_SIZE, offs, dims)*EDGE + UP + 1) + &
            velo(idx(-1, PATCH_SIZE, offs, dims)*EDGE + RT + 1)   &
            )

       if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
          wc_u(EDGE*id_chd +UP+1) = &
               wc_u(EDGE*id_chd +UP+1) - v

          wc_u(EDGE*idN_chd+UP+1) = &
               wc_u(EDGE*idN_chd+UP+1) + v
       end if

    else if (c == IPLUSJMINUS) then
       ! Parts 5 and 6 of hexagon IPLUSJMINUS, at the lower-right
       ! corner of the lozenge, combine to form a pentagon.
       !
       ! pedlen(EDGE*idS + UP + 1) = 0 in this case.
       id_chd  = idx(LAST-1, 0, offs_chd, dims_chd)
       idE_chd = idx(LAST,   0, offs_chd, dims_chd)

       v = ( &
            Iu_Base_Wgt(7) + &
            real(dom%I_u_wgt%elts(idE_chd+1)%enc(7), kind=dp) &
            ) * ( &
            velo(idx(PATCH_SIZE,  0, offs, dims)*EDGE + RT + 1) + &
            velo(idx(PATCH_SIZE, -1, offs, dims)*EDGE + UP + 1)   &
            )

       if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
          wc_u(EDGE*id_chd +RT+1) = &
               wc_u(EDGE*id_chd +RT+1) - v

          wc_u(EDGE*idE_chd+RT+1) = &
               wc_u(EDGE*idE_chd+RT+1) + v
       end if
    end if

    if (c /= IJMINUS) return

    ! Parts 4 and 5 of hexagon IJMINUS, at the lower-left corner
    ! of the lozenge, combine to form a pentagon.
    !
    ! pedlen(EDGE*idSW + DG + 1) = 0 in this case.
    id_chd  = idx(0, 0, offs_chd, dims_chd)
    idN_chd = idx(0, 1, offs_chd, dims_chd)
    idE_chd = idx(1, 0, offs_chd, dims_chd)

    u = velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+UP+1) >= ADJZONE) then
       wc_u(EDGE*id_chd +UP+1) = &
            wc_u(EDGE*id_chd +UP+1) + u(1)

       wc_u(EDGE*idN_chd+UP+1) = &
            wc_u(EDGE*idN_chd+UP+1) - u(1)
    end if

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) >= ADJZONE) then
       wc_u(EDGE*id_chd +RT+1) = &
            wc_u(EDGE*id_chd +RT+1) + u(2)

       wc_u(EDGE*idE_chd+RT+1) = &
            wc_u(EDGE*idE_chd+RT+1) - u(2)
    end if
  end subroutine Compute_velo_wavelets_penta
  

  subroutine scalar_restriction (q, itype, fine, coarse)
    ! Restrict scalar data from level fine to level coarse.
    ! Defaults to level_end and level_start when limits are omitted.
    !
    ! itype = "ss" : subsamplings
    ! itype = "fw" : full-weighting restriction
    implicit none

    type(Float_Field), intent(inout), target :: q
    character(len=*),  intent(in)            :: itype
    integer, optional, intent(in)            :: fine, coarse

    integer :: d, jmax, jmin, l

    jmax = level_end
    if (present(fine)) jmax = fine

    jmin = level_start
    if (present(coarse)) jmin = coarse

    if (jmax < jmin) then
       write(6,'(A)') &
            "ERROR: fine level is less than coarse level in scalar_restriction"
       call abort_run()
    end if

    do l = jmax-1, jmin, -1
       call update_bdry(q, l+1, 929)

       do d = 1, size(grid)
          scalar => q%data(d)%elts

          select case (itype)
          case ("fw")
             ! Full-weighting restriction; en=1 includes poles.
             call apply_interscale_d( &
                  Restrict_fwr, grid(d), l, z_null, 0, 1)

          case ("ss")
             ! Subsampling restriction; en=1 includes poles.
             call apply_interscale_d( &
                  Restrict_ss, grid(d), l, z_null, 0, 1)

          case default
             nullify(scalar)

             write(6,'(A,A)') &
                  "ERROR: invalid restriction type in scalar_restriction: ", &
                  trim(itype)

             call abort_run()
          end select

          nullify(scalar)
       end do

       q%bdry_uptodate = .false.
    end do
  end subroutine scalar_restriction
  

  subroutine Restrict_ss (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Sub sampling restriction, as for Bernoulli and Exner functions
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

    scalar(id_par+1) = scalar(id_chd+1)
  end subroutine Restrict_ss
  
  
  subroutine Restrict_fwr (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Full weighting (mass conserving) restriction as alternative to sub-sampling for topography
    
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_chd, id_par
    integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    if (dom%mask_n%elts(id_chd+1) == 0) return

    id_par = idx (i_par, j_par, offs_par, dims_par)

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

    scalar(id_par+1) = (scalar(id_chd+1)/dom%areas%elts(id_chd+1)%hex_inv + &
         scalar(idE+1)   * dom%overl_areas%elts(idE+1)%a(1)   + &
         scalar(idNE+1)  * dom%overl_areas%elts(idNE+1)%a(2)  + &
         scalar(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
         scalar(id2NE+1) * dom%overl_areas%elts(id2NE+1)%a(4) + &
         scalar(idN+1)   * dom%overl_areas%elts(idN+1)%a(1)   + &
         scalar(idW+1)   * dom%overl_areas%elts(idW+1)%a(2)   + &
         scalar(idNW+1)  * dom%overl_areas%elts(idNW+1)%a(3)  + &
         scalar(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
         scalar(idSW+1)  * dom%overl_areas%elts(idSW+1)%a(1)  + &
         scalar(idS+1)   * dom%overl_areas%elts(idS+1)%a(2)   + &
         scalar(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
         scalar(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4) ) * dom%areas%elts(id_par+1)%hex_inv
  end subroutine Restrict_fwr
  

  subroutine Prolong_full_weighting (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Prolong scalars to fine points existing at coarse scale by undoing lifting
    
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd
  
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

    scalar(id_chd+1) = (scalar(id_par+1)/dom%areas%elts(id_par+1)%hex_inv - &
         (scalar(idE+1)  * dom%overl_areas%elts(idE+1)%a(1)   + &
         scalar(idNE+1)  * dom%overl_areas%elts(idNE+1)%a(2)  + &
         scalar(idN2E+1) * dom%overl_areas%elts(idN2E+1)%a(3) + &
         scalar(id2NE+1) * dom%overl_areas%elts(id2NE+1)%a(4) + &
         scalar(idN+1)   * dom%overl_areas%elts(idN+1)%a(1)   + &
         scalar(idW+1)   * dom%overl_areas%elts(idW+1)%a(2)   + &
         scalar(idNW+1)  * dom%overl_areas%elts(idNW+1)%a(3)  + &
         scalar(idS2W+1) * dom%overl_areas%elts(idS2W+1)%a(4) + &
         scalar(idSW+1)  * dom%overl_areas%elts(idSW+1)%a(1)  + &
         scalar(idS+1)   * dom%overl_areas%elts(idS+1)%a(2)   + &
         scalar(id2SW+1) * dom%overl_areas%elts(id2SW+1)%a(3) + &
         scalar(idSE+1)  * dom%overl_areas%elts(idSE+1)%a(4))) * dom%areas%elts(id_chd+1)%hex_inv
  end subroutine Prolong_full_weighting

  
  function velo_interp_penta_corr (dom, offs, dims, offs_chd, dims_chd) result(val)
    
    implicit none

    type(Domain),                   intent(inout) :: dom
    integer, dimension(N_BDRY+1),   intent(in)    :: offs, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims, dims_chd
    real(dp)                                      :: val(2)
    
    integer :: i, j, i_chd, j_chd

    i     = 0
    j     = 0
    i_chd = 0
    j_chd = 0

    val = [  &
         (Iu_Base_Wgt(9) &
         + real (dom%I_u_wgt%elts(idx__fast(i_chd+end_pt(1,2,UP+1), j_chd+end_pt(2,2,UP+1), offs_chd(1))+1)%enc(9),kind=dp)) * &
         ((-velo(idx(0, -1, offs, dims)*EDGE+UP+1) - (-velo(idx(-1, -1, offs, dims)*EDGE+1))) - &
          (velo(ed_idx(i+end_pt(1,1,UP+1), j+end_pt(2,1,UP+1), hex_sides(:,hex_s_offs(UP+1)+0+1), offs, dims)+1) - &
           velo(ed_idx(i+opp_no(1,2,UP+1), j+opp_no(2,2,UP+1), hex_sides(:,hex_s_offs(UP+1)+1+1), offs, dims)+1))), &
         (Iu_Base_Wgt(6) &
         + real (dom%I_u_wgt%elts(idx__fast(i_chd + end_pt(1,2,RT+1), j_chd + end_pt(2,2,RT+1), offs_chd(1))+1)%enc(6),kind=dp)) * &
         (velo(idx(-1, -1, offs, dims)*EDGE+1) + velo(idx(-1, 0, offs, dims)*EDGE + RT+1) &
         - (velo(ed_idx(i+opp_no(1,1,RT+1), j+opp_no(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+1+1), offs, dims)+1) &
         -  velo(ed_idx(i+end_pt(1,1,RT+1), j+end_pt(2,1,RT+1), hex_sides(:,hex_s_offs(RT+1)+2+1), offs, dims)+1))) ]
  end function velo_interp_penta_corr

  
  subroutine Reconstruct_outer_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at 6 outer fine edges by interpolating and adding wavelet coefficients
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer  :: e, id_par, id1, id2
    real(dp) :: u
    
    id_par = idx (i_par, j_par, offs_par, dims_par)

    do e = 1, EDGE
       id1    = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2    = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)

       u = Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)
       velo(EDGE*id2+e) = u + wc_u(EDGE*id2+e)
       
       velo(EDGE*id1+e) = 2 * velo(EDGE*id_par+e) - velo(EDGE*id2+e)  ! to ensure Restriction(Prolongation) = Identity
    end do
  end subroutine Reconstruct_outer_velo

  
  subroutine Reconstruct_inner_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at 6 inner fine edges by interpolating and adding wavelet coefficients
    implicit none
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer                :: id_chd, idN_chd, idE_chd, idNE_chd
    real(dp), dimension(6) :: u_inner

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

  
  function Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd) result(val)
    ! Interpolate inner velocities to fine edges
    
    implicit none
    
    type(Domain), intent (in) :: dom
    integer,      intent(in)  :: i_par, j_par, i_chd, j_chd
    integer,      intent(in)  :: offs_par(N_BDRY+1), offs_chd(N_BDRY+1)
    integer,      intent(in)  :: dims_par(2,N_BDRY+1), dims_chd(2,N_BDRY+1)
    real(dp)                  :: val(6)
    
    integer                        :: id, id_par, id1_par, id2_par, t, idN, idE, idNE, idN2E, id2NE, idN2, idE2
    real(dp), dimension(LORT:UPLT) :: curl_u

    id_par = idx(i_par, j_par, offs_par, dims_par)

    curl_u = 0.0_dp
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

    val(1) = (dom%triarea%elts(TRIAG*id+LORT+1) * curl_u(LORT) &
         - velo(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1) - velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)) &
         / dom%len%elts(EDGE*idE+UP+1)

    val(2) = (dom%triarea%elts(TRIAG*idE+LORT+1) * curl_u(LORT) &
         - velo(EDGE*idE+RT+1) * dom%len%elts(EDGE*idE+RT+1) - velo(EDGE*idE2+UP+1) * dom%len%elts(EDGE*idE2+UP+1)) &
         / dom%len%elts(EDGE*idE+DG+1)

    val(3) = (dom%triarea%elts(TRIAG*idNE+LORT+1) * curl_u(LORT) &
         - velo(EDGE*idNE+DG+1) * dom%len%elts(EDGE*idNE+DG+1) - velo(EDGE*idN2E+UP+1) * dom%len%elts(EDGE*idN2E+UP+1)) &
         / dom%len%elts(EDGE*idNE+RT+1)

    val(4) = (dom%triarea%elts(TRIAG*id+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1) - velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)) &
         / dom%len%elts(EDGE*idN+RT+1)

    val(5) = (dom%triarea%elts(TRIAG*idN+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*idN+UP+1) * dom%len%elts(EDGE*idN+UP+1) - velo(EDGE*idN2+RT+1) * dom%len%elts(EDGE*idN2+RT+1)) &
         / dom%len%elts(EDGE*idN+DG+1)

    val(6) = (dom%triarea%elts(TRIAG*idNE+UPLT+1) * curl_u(UPLT) &
         - velo(EDGE*idNE+DG+1) * dom%len%elts(EDGE*idNE+DG+1) - velo(EDGE*id2NE+RT+1) * dom%len%elts(EDGE*id2NE+RT+1)) &
         / dom%len%elts(EDGE*idNE+UP+1)
  end function Interp_inner_velo

  
  function Interp_outer_velo (dom, i_par, j_par, ide, e, offs, dims) result(val)
    ! Interpolate outer velocities to fine edge EDGE*ide + e
    implicit none
    
    type(Domain),                   intent(in) :: dom
    integer,                        intent(in) :: i_par, j_par, ide, e
    integer, dimension(N_BDRY+1),   intent(in) :: offs
    integer, dimension(2,N_BDRY+1), intent(in) :: dims
    real(dp)                                   :: val

    real(dp), dimension(9) :: wgt

    wgt = Iu_Base_Wgt + real (dom%I_u_wgt%elts(ide+1)%enc, kind=dp)

   val = sum (wgt * [ &
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
         velo(ed_idx(i_par + opp_no(1,2,e), j_par + opp_no(2,2,e), hex_sides(:,hex_s_offs(e)+1+1), offs, dims) + 1) ])
  end function Interp_outer_velo

  
  subroutine Reconstruct_velo_penta (dom, p, c, offs, dims, z_null)
    ! Interpolate velocity to fine edges at pentagons
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: p, c, z_null
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims

    integer                        :: id_chd,  idE_chd, idN_chd, p_chd
    integer, dimension(N_BDRY+1)   :: offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_chd
    real(dp)                       :: v
    real(dp), dimension(2)         :: u

    p_chd = dom%patch%elts(p+1)%children(c-4)

    if (p_chd == 0) return

    call get_offs_Domain (dom, p_chd, offs_chd, dims_chd)

    if (c == IMINUSJPLUS) then
       id_chd  = idx (0, LAST-1, offs_chd, dims_chd)
       idN_chd = idx (0, LAST,   offs_chd, dims_chd)

       v = (Iu_Base_Wgt(8) + real (dom%I_u_wgt%elts(idN_chd+1)%enc(8),kind=dp))*( &
              velo(EDGE*idx( 0, PATCH_SIZE, offs, dims)+UP+1) &
            + velo(EDGE*idx(-1, PATCH_SIZE, offs, dims)+RT+1))

       velo(EDGE*id_chd+UP+1)  = velo(EDGE*id_chd +UP+1) - v
       velo(EDGE*idN_chd+UP+1) = velo(EDGE*idN_chd+UP+1) + v
    else
       if (c == IPLUSJMINUS) then
          id_chd  = idx (LAST-1, 0, offs_chd, dims_chd)
          idE_chd = idx (LAST,   0, offs_chd, dims_chd)

          v = -(Iu_Base_Wgt(7) + real (dom%I_u_wgt%elts(idE_chd+1)%enc(7),kind=dp))*( &
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

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id_par

    id_par   = idx (i_par, j_par, offs_par, dims_par)
    id_chd   = idx (i_chd, j_chd, offs_chd, dims_chd)
    
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)

    if (dom%mask_e%elts(EDGE*id_chd+RT+1) > 0) velo(EDGE*id_par+RT+1) = 0.5 * (velo(EDGE*id_chd+RT+1)   + velo(EDGE*idE_chd+RT+1))
    if (dom%mask_e%elts(EDGE*id_chd+DG+1) > 0) velo(EDGE*id_par+DG+1) = 0.5 * (velo(EDGE*idNE_chd+DG+1) + velo(EDGE*id_chd+DG+1))
    if (dom%mask_e%elts(EDGE*id_chd+UP+1) > 0) velo(EDGE*id_par+UP+1) = 0.5 * (velo(EDGE*id_chd+UP+1)   + velo(EDGE*idN_chd+UP+1))
  end subroutine Restrict_velo

  
  subroutine basic_F_restr_wgt (dom, i_par, j_par, e, offs_par, dims_par, i0, j0, offs, dims, typ)
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, e, i0, j0
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par
    integer, dimension(2),          intent(in)    :: typ

    integer                        :: i, j, k

    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    integer, dimension(3,16)       :: ije

    integer, dimension(2)          :: ij_nbp_mp, ij_nbp_pp, ij_nbp_pm, ij_nbp_mm
    integer, dimension(3)          :: ije_lcsd
    integer, dimension(4)          :: id_enc
    real(dp), dimension(6)         :: wgt

    if (e == UP) then
       id_enc = [idx(i0-2, j0, offs, dims), idx(i0-2,j0+1, offs, dims), idx(i0+1, j0, offs, dims), idx(i0+1, j0+1, offs, dims)]
       i = i0
       j = j0+1
    elseif (e == RT) then
       id_enc = [idx(i0, j0, offs, dims), idx(i0, j0+1, offs, dims), idx(i0+1, j0-2, offs, dims), idx(i0+1, j0-1, offs, dims)]
       i = i0+1
       j = j0
    else
       write(0,'(a)') 'Error 447: R_F_wgts'
       stop
    end if

    ije(:,UMZ+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+1+1)
    ije(:,UPZ+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+4+1)
    ije(:,WMP+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+0+1)
    ije(:,VPP+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+5+1)
    ije(:,WPM+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+3+1)
    ije(:,VMM+1) = [i, j, 0] + hex_sides(:,hex_s_offs(e+1)+2+1)

    ij_nbp_mp = [i, j] + nghb_pt(:,hex_s_offs(e+1)+0+1)
    ij_nbp_pp = [i, j] + nghb_pt(:,hex_s_offs(e+1)+5+1)
    ij_nbp_pm = [i, j] + nghb_pt(:,hex_s_offs(e+1)+3+1)
    ij_nbp_mm = [i, j] + nghb_pt(:,hex_s_offs(e+1)+2+1)

    ije(:,VMP +1) = [ij_nbp_mp(1), ij_nbp_mp(2), 0] + hex_sides (:,hex_s_offs(e+1)+4-2+1)
    ije(:,VMPP+1) = [ij_nbp_mp(1), ij_nbp_mp(2), 0] + hex_sides (:,hex_s_offs(e+1)+1+4+1)
    ije(:,UZP +1) = [ij_nbp_mp(1), ij_nbp_mp(2), 0] + hex_sides (:,hex_s_offs(e+1)+0+4+1)
    ije(:,WPPP+1) = [ij_nbp_pp(1), ij_nbp_pp(2), 0] + hex_sides (:,hex_s_offs(e+1)+4-4+1)
    ije(:,WPP +1) = [ij_nbp_pp(1), ij_nbp_pp(2), 0] + hex_sides (:,hex_s_offs(e+1)+1+2+1)
    ije(:,VPM +1) = [ij_nbp_pm(1), ij_nbp_pm(2), 0] + hex_sides (:,hex_s_offs(e+1)+1+4+1)
    ije(:,VPMM+1) = [ij_nbp_pm(1), ij_nbp_pm(2), 0] + hex_sides (:,hex_s_offs(e+1)+4-2+1)
    ije(:,UZM +1) = [ij_nbp_pm(1), ij_nbp_pm(2), 0] + hex_sides (:,hex_s_offs(e+1)+3-2+1)
    ije(:,WMMM+1) = [ij_nbp_mm(1), ij_nbp_mm(2), 0] + hex_sides (:,hex_s_offs(e+1)+1+2+1)
    ije(:,WMM +1) = [ij_nbp_mm(1), ij_nbp_mm(2), 0] + hex_sides (:,hex_s_offs(e+1)+4-4+1)

    k = 1
    if (dist(dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(ije(1,UZP+1),ije(2,UZP+1), &
         adj_tri(:,-k+2,ije(3,UZP+1)+1),offs,dims)+1)) < eps(radius)) then ! COINCIDE
       dom%R_F_wgt%elts(id_enc(1)+1)%enc = 0.0_dp
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = 0.0_dp
    else

       if (typ(k+1) == OUTER1) then
          ije_lcsd = ije(:,VPP+1)
       else if (typ(k+1) == OUTER2) then
          ije_lcsd = ije(:,WMP+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZP+1)
       end if

       wgt = interp_F_wgts (e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,k+1,e+1), offs_par, dims_par)+1), &
            ije, [VPP, WMP, UZP, UPZ, WPP, VMP, UMZ, WPPP, VMPP])

       if (e == RT) then
          wgt = [wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)]
       else if (e == UP) then
          wgt = [wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)]
       end if

       dom%R_F_wgt%elts(id_enc(1)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(2)+1)%enc = wgt(4:6)
    end if

    k = 0
    if (dist (&
         dom%ccentre%elts(tri_idx(i_par,j_par,adj_tri(:,k+1,e+1),offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(ije(1,UZM+1),ije(2,UZM+1), adj_tri(:,-k+2,ije(3,UZM+1)+1),offs,dims)+1)) &
         < eps(radius)) then ! Coincide
       dom%R_F_wgt%elts(id_enc(3)+1)%enc = 0.0_dp
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = 0.0_dp
    else

       if (typ(k+1) == OUTER1) then
          ije_lcsd = ije(:,VMM+1)
       else if (typ(k+1) == OUTER2) then
          ije_lcsd = ije(:,WPM+1)
       else ! INSIDE
          ije_lcsd = ije(:,UZM+1)
       end if

       wgt = interp_F_wgts (e, k, ije_lcsd, dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,k+1,e+1), offs_par, dims_par)+1), &
            ije, [VMM, WPM, UZM, UMZ, WMM, VPM, UPZ, WMMM, VPMM])

       if (e == UP) then
          wgt = [wgt(3), wgt(1), wgt(2), wgt(6), wgt(4), wgt(5)]
       else if (e == RT) then
          wgt = [wgt(2), wgt(3), wgt(1), wgt(5), wgt(6), wgt(4)]
       end if

       dom%R_F_wgt%elts(id_enc(3)+1)%enc = wgt(1:3)
       dom%R_F_wgt%elts(id_enc(4)+1)%enc = wgt(4:6)
    end if
    
  contains
    
    function interp_F_wgts (e, k1, ije_lcsd, endpt_o, ije, stencil) result(val)

      implicit none

      integer,     intent(in) :: e, k1
      integer,     intent(in) :: ije_lcsd(3)
      integer,     intent(in) :: ije(3,16)
      integer,     intent(in) :: stencil(9)
      type(Coord), intent(in) :: endpt_o
      real(dp)                :: val(6)
      
      integer                  :: id_tri, info
      integer,  dimension(6)   :: ipiv

      real(dp), dimension(6,6) :: G
      real(dp), dimension(6)   :: b
      type(Coord)              :: endpt, x, y

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
      b = coords_to_row_perp ([dom%ccentre%elts(id_tri+1), endpt], x, y)

      if (lapack) then
         call dgesv (6, 1, G, 6, ipiv, b, 6, info)
      else
         call LU (6, G, b, info)
      end if
     
      val = b
    end function interp_F_wgts

    function coords_to_row_perp (coords, x, y) result(val)
      
      implicit none

      type(Coord), intent(in) :: coords(2)
      type(Coord), intent(in) :: x, y
      real(dp)                :: val(6)
      
      type(Coord) :: dirvec, midpt

      midpt = mid_pt (coords(1), coords(2))
      dirvec = cross (vector(coords(1), coords(2)), midpt)
      
      val = coords_to_rowd (midpt, dirvec, x, y) * dist (coords(1), coords(2))
    end function coords_to_row_perp

    function coords_to_row (ije0, x, y) result(val)

      implicit none

      integer,     intent(in) :: ije0(3)
      type(Coord), intent(in) :: x, y
      real(dp)                :: val(6)
      
      integer     :: i, j, e
      type(Coord) :: endpt1, endpt2, midpt
      real(dp)    :: pedlen

      i = ije0(1)
      j = ije0(2)
      e = ije0(3)

      endpt1 = dom%node%elts(idx(i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), offs, dims)+1)
      endpt2 = dom%node%elts(idx(i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), offs, dims)+1)

      pedlen = dist(dom%ccentre%elts(tri_idx(i, j, adj_tri(:,1,e+1), offs, dims)+1), &
          dom%ccentre%elts(tri_idx(i, j, adj_tri(:,2,e+1), offs, dims)+1))

      midpt = mid_pt (endpt1, endpt2)

      val = coords_to_rowd (midpt, vector(endpt1, endpt2), x, y)*pedlen
    end function coords_to_row
  end subroutine basic_F_restr_wgt
  

  subroutine init_wavelets
    implicit none
    integer :: d, i, k, num, v

    do k = zmin, zmax
       do v = 1, N_VARIABLE
          call init_Field (wav_coeff(v,k), POSIT(v))
       end do
    end do
    
    if (vert_diffuse) then
       do k = 1, zlevels
          call init_Field (wav_tke(k), AT_NODE)
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
          call init_RF_Wgt (grid(d)%R_F_wgt%elts(i), [0.0_dp, 0.0_dp, 0.0_dp])
       end do

       do k = zmin, zmax
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
    
    type(Domain),                   intent(in)  :: dom
    integer,                        intent(in)  :: e, i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1),   intent(in)  :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)  :: dims_par, dims_chd
    integer, dimension(2),          intent(out) :: typ
    real(dp), dimension(8),         intent(out) :: area

    type(Coord), dimension(6)   :: hex
    type(Coord), dimension(3,2) :: tri
    type(Coord)                 :: intersection_pt0, intersection_pt1, pt
    integer                     :: i, id_chd
    logical                     :: intersects0, intersects1, degenerate

    area = 0.0_dp
    typ = 0

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    hex = [ (dom%ccentre%elts(tri_idx(i_chd, j_chd, no_adj_tri(:,i + &
         hex_s_offs(-e+3)+1), offs_chd, dims_chd)+1), i = 0, 6-1) ]

    tri = reshape([dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,4,e+1), &
         offs_par, dims_par)+1), dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,2,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,3,e+1), &
         offs_par, dims_par)+1), dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,2,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, adj_tri(:,1,e+1), offs_par, dims_par)+1), &
         dom%ccentre%elts(tri_idx(i_par, j_par, bfly_tri(:,1,e+1), offs_par, dims_par)+1)], [3, 2])

    pt = dom%node%elts(id_chd+1)

    area(1) = triarea(hex(6), hex(1), pt)
    area(2) = triarea(hex(3), hex(4), pt)

    do i = 1, 2
       call arc_intersect_test (tri(1,i), tri(2,i), hex(3*i-2), hex(3*i-1), intersection_pt0, intersects0, degenerate)
       call arc_intersect_test (tri(3,i), tri(2,i), hex(3*i),   hex(3*i-1), intersection_pt1, intersects1, degenerate)
       if (intersects0 .and. intersects1) then
          area(i+4) = triarea(intersection_pt0, tri(2,i), hex(3*i-1))
          area(i+6) = triarea(tri(2,i), hex(3*i-1), intersection_pt1)
          area(i+2) = area(i+4) + area(i+6)
          area(i) = area(i) + triarea(hex(3*i-2), intersection_pt0, pt)     + triarea(intersection_pt0, pt, tri(2,i))
          area(-i+3) = area(-i+3) + triarea(intersection_pt1, hex(3*i), pt) + triarea(tri(2,i), pt, intersection_pt1)
          typ(-i+3) = INSIDE
       else
          if (.not. intersects0 .and. .not. intersects1) then
             area(i+2) = 0.0_dp
             call arc_intersect_test (tri(2,1), tri(2,2), hex(3*i-2), hex(3*i-1), intersection_pt0, intersects0, degenerate)
             call arc_intersect_test (tri(2,2), tri(2,1), hex(3*i-1), hex(3*i),   intersection_pt1, intersects1, degenerate)
             if (.not. intersects0 .and. intersects1) then
                area(i) = area(i) + triarea(hex(3*i-2), hex(3*i-1), pt) + triarea(hex(3*i-1), intersection_pt1, pt)
                area(-i+3) = area(-i+3) + triarea(intersection_pt1, hex(3*i), pt)
                typ(-i+3) = OUTER2
             else
                if (intersects0 .and. .not. intersects1) then
                   area(i) = area(i) + triarea(hex(3*i-2), intersection_pt0, pt)
                   area(-i+3) = area(-i+3) + triarea(hex(3*i-1), hex(3*i), pt) + triarea(intersection_pt0, hex(3*i-1), pt)
                   typ(-i+3) = OUTER1
                else
                   write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'A', intersects0, intersects1
                end if
             end if
          else
             write(0,*) 'ERROR: overlap area', dom%id, offs_chd(1), i_chd, j_chd, 'B', intersects0, intersects1
          end if
       end if
    end do
    return
  end subroutine get_overl_areas
  

  subroutine normalize2 (q, u, v)
    implicit none
    real(dp), dimension(2), intent(in)  :: q
    real(dp),               intent(out) :: u, v

    real(dp) :: nrm

    nrm = sqrt (q(1)**2 + q(2)**2)
    u = q(1)/nrm
    v = q(2)/nrm
  end subroutine normalize2
  

  real(dp) function Interp_node (dom, id, id1, id2, id3, id4)
    ! Interpolation at nodes
    type(Domain), intent(in) :: dom
    integer,      intent(in) :: id, id1, id2, id3, id4

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
    
    type(Coord), intent(in)  :: midpt, endpt1, endpt2
    type(Coord), intent(out) :: x, y

    type(Coord) :: y0

    x = direction (endpt1, endpt2)
    y0 = cross (x, midpt)
    y = normalize_Coord (y0)
  end subroutine local_coord
  

  subroutine set_RF_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i_par, j_par, i_chd, j_chd,  p_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer                :: id_chd, idN_chd, idE_chd, idNE_chd
    integer,  dimension(2) :: typ
    real(dp), dimension(8) :: area

    id_chd   = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)

    call get_overl_areas (dom, i_par, j_par, i_chd+1, j_chd, offs_par, dims_par, offs_chd, dims_chd, RT, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idE_chd+1), area)
    call basic_F_restr_wgt (dom, i_par, j_par, RT, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)

    call get_overl_areas (dom, i_par, j_par, i_chd+1, j_chd+1, offs_par, dims_par, offs_chd, dims_chd, DG, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idNE_chd+1), area)

    call get_overl_areas (dom, i_par, j_par, i_chd, j_chd+1, offs_par, dims_par, offs_chd, dims_chd, UP, area, typ)
    call init_Overl_Area (dom%overl_areas%elts(idN_chd+1), area)
    call basic_F_restr_wgt (dom, i_par, j_par, UP, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd, typ)
  end subroutine set_RF_wgts
  

  subroutine set_WT_wgts (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Sets local weights used in outer velocity interpolation when refining the grid
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd

    id_chd  = idx (i_chd,   j_chd,   offs_chd, dims_chd)
    idN_chd = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd = idx (i_chd+1, j_chd,   offs_chd, dims_chd)

    dom%I_u_wgt%elts(idE_chd+1) = outer_velo_weights (dom, p_chd, i_par, j_par, RT, offs_par, dims_par)
    dom%I_u_wgt%elts(id_chd+1)  = outer_velo_weights (dom, p_chd, i_par, j_par, DG, offs_par, dims_par)
    dom%I_u_wgt%elts(idN_chd+1) = outer_velo_weights (dom, p_chd, i_par, j_par, UP, offs_par, dims_par)
  end subroutine set_WT_wgts
  

  function outer_velo_weights (dom, p, i0, j0, e0, offs, dims) result(val)
    
    implicit none
    
    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: p, i0, j0, e0
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    type(Iu_Wgt)                                  :: val

    integer                  :: k, id, info
    integer,  dimension(6)   :: ipiv
    real(dp), dimension(9)   :: weights
    real(dp), dimension(6)   :: b
    real(dp), dimension(6,6) :: G
    type(Coord)              :: x, y
    
    call local_coord (&
         dom%midpt%elts(idx(i0, j0, offs, dims)*EDGE + e0 + 1), &
         dom%node%elts(idx(i0 + end_pt(1,1,e0+1), j0 + end_pt(2,1,e0+1), offs, dims) + 1), &
         dom%node%elts(idx(i0 + end_pt(1,2,e0+1), j0 + end_pt(2,2,e0+1), offs, dims) + 1), x, y)

    weights = 0.0_dp

    do k = 1, 2
       id = idx (i0, j0, offs, dims)

       G = 0.0_dp

       G(:,1) = coords_to_row (i0, j0, [0, 0], [0, 0, e0], e0)

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

       if (lapack) then
          call dgesv (6, 1, G, 6, ipiv, b, 6, info)
       else
          call LU (6, G, b, info)
       end if
       
       weights(1)           = weights(1) + b(1)
       weights(2*k:2*k+1)   = weights(2*k:2*k+1) + b(2:3)
       weights(2*k+4:2*k+5) = b(4:5)
       weights(-2*k+6)      = weights(-2*k+6) + b(6)
       weights(-2*k+7)      = weights(-2*k+7) - b(6)
    end do

    val = Iu_wgt (0.5 * weights - Iu_Base_Wgt)
    
  contains
    
    function coords_to_row (i00, j00, n_offs1, n_offs2, e00) result(val)

      implicit none
      
      integer, intent(in) :: e00, i00, j00
      integer, intent(in) :: n_offs1(2)
      integer, intent(in) :: n_offs2(3)
      real(dp)            :: val(6)
      
      type(Coord) :: endpt1, endpt2
      integer     :: i, j, e

      i = i00 + n_offs1(1) + n_offs2(1)
      j = j00 + n_offs1(2) + n_offs2(2)

      e = n_offs2(3)
      
      endpt1 = get_coord (i + end_pt(1,1,e+1), j + end_pt(2,1,e+1), e00)
      endpt2 = get_coord (i + end_pt(1,2,e+1), j + end_pt(2,2,e+1), e00)

      val = coords_to_rowd (mid_pt(endpt1, endpt2), vector(endpt1, endpt2), x, y)
    end function coords_to_row

    function get_coord (i, j, e) result(val)
      integer, intent(in) :: i, j, e
      type(Coord)         :: val

      integer :: id
      
      id = idx (i, j, offs, dims)

      val = ORIGIN

      if (i == -1) then
         if (j == -1 .and. is_penta (dom, p, IJMINUS-1)) then
            if (e == RT) then
               val = dom%node%elts(nidx(LAST_BDRY, 0, IMINUS, offs, dims)+1)
               return
            else
               if (e == UP) then
                  val = dom%node%elts(nidx(0, LAST_BDRY, JMINUS, offs, dims)+1)
                  return
               end if
            end if
         else
            if (j == PATCH_SIZE .and. is_penta(dom, p, IMINUSJPLUS-1)) then
               val = dom%node%elts(nidx(0, 1, JPLUS, offs, dims)+1)
               return
            else
               val = dom%node%elts(id+1)
               return
            end if
         end if
      else
         if (i == PATCH_SIZE .and. j == -1 .and. is_penta(dom, p, IPLUSJMINUS-1)) then
            val = dom%node%elts(nidx(1, 0, IPLUS, offs, dims)+1)
            return
         else
            if (i == PATCH_SIZE+1 .and. j == PATCH_SIZE+1 .and. is_penta(dom, p, IJPLUS-1)) then
               val = dom%node%elts(nidx(1, 0, IJPLUS, offs, dims)+1)
               return
            else
               val = dom%node%elts(id+1)
               return
            end if
         end if
      end if
    end function get_coord
    
  end function outer_velo_weights
  

  subroutine check_m (dom, i_par, j_par, i_chd, j_chd, offs_par, dims_par,  offs_chd, dims_chd)
    ! Check that scalar is indeed conserved by restriction
    implicit none
    
    type(Domain),                   intent(in) :: dom
    integer,                        intent(in) :: i_par, j_par, i_chd, j_chd
    integer, dimension(N_BDRY+1),   intent(in) :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1), intent(in) :: dims_par, dims_chd

    integer  :: id_chd, id_par, idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE
    real(dp) :: ratio

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)
    id_par = idx (i_par, j_par, offs_par, dims_par)

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

    ratio = (1.0_dp/dom%areas%elts(id_chd+1)%hex_inv + &
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
  

  function coord2local (c, x, y) result(val)

    implicit none
    
    type(Coord), intent(in) :: c, x, y
    real(dp)                :: val(2)
    
    val = [inner(c, x), inner(c, y)]
  end function coord2local
  

  function coords_to_rowd (midpt, dirvec, x, y) result(val)

    implicit none
    
    type(Coord), intent(in) :: dirvec, x, y
    real(dp)                :: val(6)
    
    type(Coord)            :: midpt
    real(dp)               :: u, v
    real(dp), dimension(2) :: xy

    call normalize2 (coord2local (dirvec, x, y), u, v)
    
    xy = coord2local (midpt, x, y)
    
    val = [u, u*xy(1), u*xy(2), v, v*xy(1), v*xy(2)]
  end function coords_to_rowd
  

  subroutine LU (n, A, b, info)
    ! Solve linear system using LU decomposition with partial pivoting
    ! result is returned in b and A is replaced by the factors L and U from the factorization
    ! A = P L U (the unit diagonal elements of L are not stored).
    implicit none
    integer, intent(in)    :: n
    integer, intent(out)   :: info
    real(8), intent(inout) :: A(n,n), b(n)

    integer  :: i, j, k, p
    integer  :: ipiv(n)
    real(dp) :: tmp, pivot_abs, cand_abs

    info = 0
    ipiv = 0

    ! Argument checks 
    if (n <  0) then; info = -1;  return; end if
    if (n == 0) return

    ! LU factorization (A = P * L * U) with partial pivoting
    do k = 1, n
       ! pivot index p = argmax_{i=k..n} |A(i,k)|
        p  = k
        pivot_abs = abs (A(k,k))
        do i = k+1, n
          cand_abs = abs (A(i,k))
          if (cand_abs > pivot_abs) then
            pivot_abs = cand_abs
            p = i
          end if
        end do
        ipiv(k) = p

        ! Singular if pivot is zero
        if (abs (A(p,k)) <= eps (1.0_dp)) then
          info = k
          return
        end if

        ! Swap rows p and k in A
        if (p /= k) then
          do j = 1, n
            tmp = A(p,j);  A(p,j) = A(k,j);  A(k,j) = tmp
          end do
        end if

        ! Compute multipliers below the pivot
        A(k+1:n,k) = A(k+1:n,k) / A(k,k)

        ! Rank-1 update of the trailing submatrix
        do i = k+1, n
          tmp = A(i,k)
          if (abs (tmp) > eps (1.0_dp)) A(i,k+1:n) = A(i,k+1:n) - tmp * A(k,k+1:n)
        end do
    end do

    ! Apply row permutations to B  (B := P * B)
    do k = 1, n
       p = ipiv(k)
       if (p /= k) then
          tmp = b(p);  b(p) = b(k);  b(k) = tmp
       end if
    end do

    ! Solve L * Y = B  (unit-lower L in A)
    do k = 1, n-1
      do i = k+1, n
         tmp = A(i,k)
         if (abs(tmp) > eps (1.0_dp)) b(i) = b(i) - tmp * b(k)
      end do
    end do
    
    ! Solve U * X = Y  (upper U in A), result X overwrites B
    do i = n, 1, -1
       ! divide row i
       b(i) = b(i) / A(i,i)
       ! eliminate above
       do k = 1, i-1
          tmp = A(k,i)
          if (abs(tmp) > eps (1.0_dp)) b(k) = b(k) - tmp * b(i)
       end do
    end do
  end subroutine LU

  
end module wavelet_mod
