module lin_solve_mod
  use ops_mod
  use comm_mpi_mod
  implicit none
  integer                            :: ii, m
  real(8)                            :: dp_loc, linf_loc, l2_loc
  real(8), pointer                   :: mu1, mu2
  real(8), dimension(:), pointer     :: scalar2, scalar3
  real(8), dimension(:), allocatable :: w
  character(4)                       :: var_type
contains
  real(8) function linf (s, l)
    ! Returns l_inf norm of scalar s at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: s

    integer :: d, j

    linf_loc = -1d16
    do d = 1, size(grid)
       scalar => s%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_linf_scalar, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar)
    end do
    linf = sync_max_real (linf_loc)
  end function linf

  subroutine cal_linf_scalar (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id) >= ADJZONE) linf_loc = max (linf_loc, abs (scalar(id)))
  end subroutine cal_linf_scalar

  real(8) function l2 (s, l)
    ! Returns l_2 norm of scalar s at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: s
    
    integer :: d, j

    l2_loc = 0.0_8
    do d = 1, size(grid)
       scalar => s%data(d)%elts
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_l2_scalar, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_l2_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar)
    end do
    l2 = sqrt (sum_real (l2_loc))
  end function l2

  subroutine cal_l2_scalar (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) l2_loc = l2_loc + scalar(id_i)**2
  end subroutine cal_l2_scalar

  subroutine cal_l2_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e

    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) l2_loc = l2_loc + scalar(id_e)**2
    end do
  end subroutine cal_l2_velo

  real(8) function dp (s1, s2, l)
    ! Calculates dot product of s1 and s2 at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: s1, s2

    integer :: d, j

    dp_loc = 0.0_8
    do d = 1, size(grid)
       scalar => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_dotproduct, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_dotproduct_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar, scalar2)
    end do
    dp = sum_real (dp_loc)
  end function dp

  subroutine cal_dotproduct (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id) >= ADJZONE) dp_loc = dp_loc + scalar(id) * scalar2(id)
  end subroutine cal_dotproduct

   subroutine cal_dotproduct_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e
    
    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) dp_loc = dp_loc +  scalar(id_e) * scalar2(id_e)
    end do
  end subroutine cal_dotproduct_velo

  subroutine lc (s1, a, s2, s3, l)
    ! Calculates linear combination of scalars s3 = (s1 + a*s2) at scale l
    implicit none
    integer                   :: l
    real(8), target           :: a
    type(Float_Field), target :: s1, s2, s3

    integer :: d, j

    do d = 1, size(grid)
       scalar  => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       scalar3 => s3%data(d)%elts
       mu1     => a
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_lc, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar, scalar2, scalar3, mu1)
    end do
    s1%bdry_uptodate = .false.
  end subroutine lc

  function lcf (s1, a, s2, l)
    ! Calculates linear combination of scalars s3 = (s1 + b*s2) at scale l
    implicit none
    integer                   :: l
    real(8), target           :: a
    type(Float_Field), target :: lcf, s1, s2

    integer :: d, j

    lcf = s1
    do d = 1, size(grid)
       scalar  => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       scalar3 => lcf%data(d)%elts
       mu1     => a
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_lc, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_lc_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar, scalar2, scalar3, mu1)
    end do
    lcf%bdry_uptodate = .false.
  end function lcf

  subroutine cal_lc (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar3(id_i) = scalar(id_i) + mu1 * scalar2(id_i)
  end subroutine cal_lc

  subroutine cal_lc_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e

    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) scalar3(id_e) = scalar(id_e) +  mu1 * scalar2(id_e)
    end do
  end subroutine cal_lc_velo

  subroutine lc2 (u, alpha, p, omega, s, l)
    ! Calculates linear combination of scalars u = u + alpha*p + omega*s at scale l
    implicit none
    integer                   :: l
    real(8), target           :: alpha, omega
    type(Float_Field), target :: u, p, s

    integer :: d, j

    do d = 1, size(grid)
       scalar  => u%data(d)%elts
       scalar2 => p%data(d)%elts
       scalar3 => s%data(d)%elts
       mu1     => alpha
       mu2     => omega
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_lc2, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_lc2_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar, scalar2, scalar3, mu1, mu2)
    end do
    u%bdry_uptodate = .false.
  end subroutine lc2

  subroutine cal_lc2 (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar(id_i) = scalar(id_i) + mu1 * scalar2(id_i) + mu2*scalar3(id_i)
  end subroutine cal_lc2
  
  subroutine cal_lc2_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e

    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) scalar(id_e) = scalar(id_e) + mu1 * scalar2(id_e) + mu2*scalar3(id_e)
    end do
  end subroutine cal_lc2_velo

  function residual (f, u, Lu, l)
    ! Calculates f - Lu(u) at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: f, residual, u

    integer :: d, j

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    residual = Lu (u, l)

    do d = 1, size(grid)
       scalar  => f%data(d)%elts
       scalar2 => residual%data(d)%elts
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_res, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_res_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar, scalar2)
    end do
    residual%bdry_uptodate = .false.
  end function residual

  subroutine cal_res (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar2(id_i) = scalar(id_i) - scalar2(id_i)
  end subroutine cal_res

  subroutine cal_res_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e

    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) scalar2(id_e) = scalar(id_e) - scalar2(id_e)
    end do
  end subroutine cal_res_velo

  subroutine set_zero (s, l)
    ! Sets scalar s to zero
    implicit none
    integer                   :: l
    type(Float_Field), target :: s

    integer :: d, j

    do d = 1, size(grid)
       scalar => s%data(d)%elts
       select case (var_type)
       case ("sclr")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_set_zero, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
       case ("velo")
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_set_zero_velo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 0)
          end do
       end select
       nullify (scalar)
    end do
    s%bdry_uptodate = .false.
  end subroutine set_zero

  subroutine cal_set_zero (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id_i

    id_i = idx (i, j, offs, dims) + 1
    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar(id_i) = 0.0_8
  end subroutine cal_set_zero

  subroutine cal_set_zero_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: e, id, id_e
    
    id = idx (i, j, offs, dims)
    do e = 1, EDGE
       id_e = EDGE*id+e
       if (dom%mask_e%elts(id_e) >= ADJZONE) scalar(id_e) = 0.0_8
    end do
  end subroutine cal_set_zero_velo

  subroutine prolongation (scaling, fine)
    ! Prolong from coarse scale fine-1 to scale fine
    use wavelet_mod
    implicit none
    integer                   :: fine
    type(Float_Field), target :: scaling

    integer :: coarse, d, l
    
    coarse = fine - 1

    select case (var_type)
    case ("sclr")
       ! Prolong scalar to finer nodes existing at coarser grid (subsample) 
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          call apply_interscale_d2 (subsample, grid(d), coarse, z_null, 0, 1)
          nullify (scalar)
       end do

       ! Reconstruct scalar at finer nodes not existing at coarser grid by interpolation
       do d = 1, size(grid)
          scalar => scaling%data(d)%elts
          call apply_interscale_d (interpolate, grid(d), coarse, z_null, 0, 0)
          nullify (scalar)
       end do
       scaling%bdry_uptodate = .false.
       call update_bdry (scaling, fine, 66)
    case ("velo")
       ! Reconstruct outer velocities at finer edges (interpolate and add wavelet coefficients)
       do d = 1, size(grid)
          velo => scaling%data(d)%elts
          call apply_interscale_d2 (interpolate_outer_velo, grid(d), coarse, z_null, 0, 1) ! needs val
          call apply_to_penta_d (reconstruct_velo_penta, grid(d), coarse, z_null)
          nullify (velo)
       end do
       scaling%bdry_uptodate = .false.
       call update_bdry (scaling, fine, 200)

       ! Reconstruct inner velocities at finer edges (interpolate and add wavelet coefficients)
       do d = 1, size(grid)
          velo => scaling%data(d)%elts
          call apply_interscale_d (interpolate_inner_velo, grid(d), coarse, z_null, 0, 0)
          nullify (velo)
       end do
       scaling%bdry_uptodate = .false.
       call update_bdry (scaling, fine, 201) ! for next outer velocity
    end select
  end subroutine prolongation

  subroutine subsample (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Sub-sample to prolong coarse points to fine grid
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_par = idx (i_par, j_par, offs_par, dims_par)
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    scalar(id_chd+1) = scalar(id_par+1) 
  end subroutine subsample

  subroutine interpolate (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct scalars at fine nodes not existing at coarse scale by interpolation
    use wavelet_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)

    ! Interpolate scalars and add wavelets to reconstruct values at fine scale
    scalar(idE_chd+1)  = Interp_node (dom,  idE_chd,    id_chd, id2E_chd, id2NE_chd,  id2S_chd)
    scalar(idNE_chd+1) = Interp_node (dom, idNE_chd, id2NE_chd,   id_chd,  id2E_chd,  id2N_chd) 
    scalar(idN_chd+1)  = Interp_node (dom,  idN_chd,    id_chd, id2N_chd,  id2W_chd, id2NE_chd)  
  end subroutine interpolate

   subroutine interpolate_outer_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at outer fine edges not existing on coarse grid
    use wavelet_mod
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: e, id_chd, id_par, id1, id2
    real(8) :: u

    id_chd = idx(i_chd, j_chd, offs_chd, dims_chd)

    do e = 1, EDGE
       id1    = idx (i_chd + end_pt(1,1,e), j_chd + end_pt(2,1,e), offs_chd, dims_chd)
       id2    = idx (i_chd + end_pt(1,2,e), j_chd + end_pt(2,2,e), offs_chd, dims_chd)
       id_par = idx (i_par, j_par, offs_par, dims_par)

       u = Interp_outer_velo (dom, i_par, j_par, id2, e, offs_par, dims_par)

       velo(EDGE*id2+e) = u 
       velo(EDGE*id1+e) = 2d0 * velo(EDGE*id_par+e) - u 
    end do
  end subroutine interpolate_outer_velo

  subroutine interpolate_inner_velo (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct velocity at inner fine edges not existing on coarse grid
    use wavelet_mod
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer               :: id_chd, idN_chd, idE_chd, idNE_chd
    real(8), dimension(6) :: u_inner

    id_chd   = idx (i_chd,     j_chd,     offs_chd, dims_chd)

    idN_chd  = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd  = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)

    u_inner = Interp_inner_velo (dom, i_par, j_par, offs_par, dims_par, i_chd, j_chd, offs_chd, dims_chd)

    velo(EDGE*idE_chd +UP+1) = u_inner(1) 
    velo(EDGE*idE_chd +DG+1) = u_inner(2) 
    velo(EDGE*idNE_chd+RT+1) = u_inner(3) 
    velo(EDGE*idN_chd +RT+1) = u_inner(4) 
    velo(EDGE*idN_chd +DG+1) = u_inner(5) 
    velo(EDGE*idNE_chd+UP+1) = u_inner(6) 
  end subroutine interpolate_inner_velo

  subroutine multiscale (u, f, Lu)
    ! Solves linear equation L(u) = f using a simple multiscale algorithm with Scheduled Rexation Jacobi iterations as the smoother
    implicit none
    type(Float_Field), target :: f, u

    integer                                       :: l, n
    real(8)                                       :: k_max, k_min
    real(8), dimension(level_start:level_end)     :: nrm_f
    real(8), dimension(1:2,level_start:level_end) :: r_error

    integer, dimension(level_start:level_end) :: iter

    interface
       function Lu (u, l)
         ! Returns result of linear operator applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface
    
    ! Optimal Scheduled Relaxation Jacobi parameters (Adsuara, et al J Comput Phys v 332, 2016)
    ! (k_min and k_max are determined empirically to give optimal convergence on fine non uniform grids)
    k_max = 1.8d0

    ! Values for k_min optimized for J5J9
    if (tol_jacobi <= 1d-5) then
       k_min = 1d-2
    else
       k_min = 3d-2
    end if
   
    m = 20; allocate (w(1:m))
    if (fine_iter < m) then ! do not use SJR
       w = 1d0
    else 
       fine_iter = m * (fine_iter / m) ! ensure that fine_iter is an integer multiple of m
       do n = 1, m
          w(n) = 2d0 / (k_min + k_max - (k_max - k_min) * cos (MATH_PI * (2d0*dble(n) - 1d0)/(2d0*dble(m))))
       end do
    end if

    var_type = "sclr"

    call update_bdry (f, NONE, 55)
    
    do l = level_start, level_end
       nrm_f(l) = l2 (f, l)
       if (nrm_f(l) < tol_elliptic) nrm_f(l) = 1d0
       if (log_iter) r_error(1,l) = l2 (residual (f, u, Lu, l), l) / nrm_f(l)
    end do

    l = level_start
    call bicgstab (u, f, nrm_f(l), Lu, l, coarse_iter, r_error(2,l), iter(l))
    do l = level_start+1, level_end
       call prolongation (u, l)
       call jacobi (u, f, nrm_f(l), Lu, l, fine_iter, r_error(2,l), iter(l))
    end do

    if (log_iter) then
       do l = level_start, level_end
          if (rank == 0) write (6, '("residual at scale ", i2, " = ", 2(es10.4,1x)," after ", i6, " iterations")') &
               l, r_error(:,l), iter(l)
       end do
    end if

    deallocate (w)
  end subroutine multiscale

 
  subroutine jacobi (u, f, nrm_f, Lu, l, max_iter, nrm_res, iter)
    ! Max_iter Jacobi iterations for smoothing multigrid iterations
    ! uses Scheduled Relaxation Jacobi (SJR) iterations (Yang and Mittal JCP 274, 2014)
    implicit none
    integer                   :: iter, l, max_iter
    real(8)                   :: nrm_f, nrm_res
    type(Float_Field), target :: f, u

    integer                   :: i
    type(Float_Field), target :: res
    
    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    ii = 0
    iter = 0
    res = residual (f, u, Lu, l)
    nrm_res = l2 (res, l) / nrm_f

    do while (iter < max_iter)
       if (nrm_res <= tol_jacobi) exit
       
       ii = ii + 1
       if (ii > m) then
          if (nrm_res < 2d0 * tol_jacobi) then ! avoid starting a new SRJ cycle if error is small enough
             exit
          else
             ii = 1
          end if
       end if
       iter = iter + 1
       call lc_jacobi (u, res, l)
       res = residual (f, u, Lu, l)
       nrm_res = l2 (res, l) / nrm_f
    end do
    u%bdry_uptodate = .false.
  end subroutine jacobi

  subroutine lc_jacobi (s1, s2, l)
    ! Jacobi iteration
    implicit none
    integer                   :: l
    type(Float_Field), target :: s1, s2

    integer :: d, j

    do d = 1, size(grid)
       scalar  => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_jacobi, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar, scalar2, scalar3)
    end do
    s1%bdry_uptodate = .false.
  end subroutine lc_jacobi

  subroutine cal_jacobi (dom, i, j, zlev, offs, dims)
    ! Jacobi iteration using local approximation of diagonal of elliptic operator
    ! (Laplacian_scalar(S_TEMP) is the old free surface perturbation)
    ! ** using exact value of diagonal of Laplacian typically does NOT improve results **
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer            :: d, id, id_i, idE, idNE, idN, idW, idSW, idS
    real(8)            :: depth, depth_e, Laplace_diag, wgt
    logical, parameter :: exact = .false.

    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       d = dom%id + 1
       depth = abs (dom%topo%elts(id_i)) + Laplacian_scalar(S_TEMP)%data(d)%elts(id_i) / phi_node (d, id_i, zlevels)

       if (.not. exact) then ! average value 
          wgt = 2d0 * sqrt (3d0) * depth
       else ! true local value
          idE  = idx (i+1, j,   offs, dims) 
          idNE = idx (i+1, j+1, offs, dims) 
          idN  = idx (i,   j+1, offs, dims) 
          idW  = idx (i-1, j,   offs, dims) 
          idSW = idx (i-1, j-1, offs, dims) 
          idS  = idx (i,   j-1, offs, dims)
       
          depth_e = abs (dom%topo%elts(idE+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idE+1) / phi_node (d, idE+1, zlevels)
          wgt = dom%pedlen%elts(EDGE*id+RT+1) / dom%len%elts(EDGE*id+RT+1) * interp (depth_e, depth)

          depth_e = abs (dom%topo%elts(idNE+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idNE+1) / phi_node (d, idNE+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*id+DG+1) / dom%len%elts(EDGE*id+DG+1) * interp (depth_e, depth)

          depth_e = abs (dom%topo%elts(idN+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idN+1) / phi_node (d, idN+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*id+UP+1) / dom%len%elts(EDGE*id+UP+1) * interp (depth_e, depth)

          depth_e = abs (dom%topo%elts(idW+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idW+1) / phi_node (d, idW+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idW+RT+1) / dom%len%elts(EDGE*idW+RT+1) * interp (depth_e, depth)

          depth_e = abs (dom%topo%elts(idSW+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idSW+1) / phi_node (d, idSW+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idSW+DG+1) / dom%len%elts(EDGE*idSW+DG+1) * interp (depth_e, depth)

          depth_e = abs (dom%topo%elts(idS+1)) + Laplacian_scalar(S_TEMP)%data(d)%elts(idS+1) / phi_node (d, idS+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idS+UP+1) / dom%len%elts(EDGE*idS+UP+1) * interp (depth_e, depth)
       end if
       Laplace_diag = - grav_accel * wgt * dt**2 * dom%areas%elts(id_i)%hex_inv

       scalar(id_i) = scalar(id_i) + w(ii) * scalar2(id_i) / (theta1*theta2 * Laplace_diag - 1d0)
    end if
  end subroutine cal_jacobi

  subroutine bicgstab (u, f, nrm_f, Lu, l, max_iter, nrm_res, iter)
    ! Solves the linear system Lu(u) = f at scale l using bi-cgstab algorithm (van der Vorst 1992).
    ! This is a conjugate gradient type algorithm.
    implicit none
    integer                   :: iter, l, max_iter
    real(8)                   :: nrm_f, nrm_res
    type(Float_Field), target :: f, u
    
    integer                   :: i
    real(8)                   :: alph, b, nrm_res0, omga, rho, rho_old
    type(Float_Field), target :: Ap, As, res, res0, p, s

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    ! Initialize
    res  = residual (f, u, Lu, l)
    res0 = res
    rho  = dp (res0, res, l)
    
    p    = res0
    Ap   = Lu (p, l) 
    
    iter = 0
    do while (iter < max_iter)
       alph = rho / dp (Ap, res0, l)

       s = lcf (res, -alph, Ap, l)
       As = Lu (s, l)

       omga = dp (As, s, l) / dp (As, As, l)
       
       call lc2 (u, alph, p, omga, s, l)
       
       res = lcf (s, -omga, As, l)
       nrm_res = l2 (res, l) / nrm_f
       if (nrm_res <= tol_elliptic) exit
       
       rho_old = rho
       rho = dp (res0, res, l)
       
       b = (alph/omga) * (rho/rho_old)
       p = lcf (res, b, lcf (p, -omga, Ap, l), l) 
       Ap = Lu (p, l)
       
       iter = iter + 1
    end do
    u%bdry_uptodate = .false.
  end subroutine bicgstab
end module lin_solve_mod
