module lin_solve_mod
  use ops_mod
  use comm_mpi_mod
  use wavelet_mod
  implicit none
  integer                            :: ii, m
  real(8)                            :: dp_loc, linf_loc, l2_loc
  real(8), dimension(:), allocatable :: w
contains
  real(8) function linf (s, l)
    ! Returns l_inf norm of scalar s at scale l
    implicit none
    integer           :: l
    type(Float_Field) :: s

    linf_loc = -1d16

    call apply_onescale (cal_linf_scalar, l, z_null, 0, 1)
   
    linf = sync_max_real (linf_loc)
  contains
    subroutine cal_linf_scalar (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1

      id = idx (i, j, offs, dims) + 1

      if (dom%mask_n%elts(id) >= ADJZONE) linf_loc = max (linf_loc, abs (s%data(d)%elts(id)))
    end subroutine cal_linf_scalar
  end function linf

  real(8) function l2 (s, l)
    ! Returns l_2 norm of scalar s at scale l
    implicit none
    integer           :: l
    type(Float_Field) :: s

    integer :: d, j

    l2_loc = 0d0
  
    call apply_onescale (cal_l2_scalar, l, z_null, 0, 1)
    
    l2 = sqrt (sum_real (l2_loc))
  contains
    subroutine cal_l2_scalar (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i

      d = dom%id + 1

      id_i = idx (i, j, offs, dims) + 1

      if (dom%mask_n%elts(id_i) >= ADJZONE) l2_loc = l2_loc +  s%data(d)%elts(id_i)**2
    end subroutine cal_l2_scalar
  end function l2
 
  real(8) function dp (s1, s2, l)
    ! Calculates dot product of s1 and s2 at scale l
    implicit none
    integer          :: l
    type(Float_Field):: s1, s2

    call update_bdry (s1, l, 50)
    call update_bdry (s2, l, 50)

    dp_loc = 0d0

    call apply_onescale (cal_dotproduct, l, z_null, 0, 1)
    
    dp = sum_real (dp_loc)
  contains
    subroutine cal_dotproduct (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1

      id = idx (i, j, offs, dims) + 1
      dp_loc = dp_loc + s1%data(d)%elts(id) * s2%data(d)%elts(id)
    end subroutine cal_dotproduct
  end function dp
 
  function lcf (a1, s1, a2, s2, l)
    ! Calculates linear combination of scalars lcf = a1*s1 + a2*s2 at scale l
    implicit none
    integer           :: l
    real(8)           :: a1, a2
    type(Float_Field) :: lcf, s1, s2
    
    lcf = s1

    call apply_onescale (cal_lc, l, z_null, 0, 1)

    lcf%bdry_uptodate = .false.
  contains
    subroutine cal_lc (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1

      lcf%data(d)%elts(id_i) = a1 * s1%data(d)%elts(id_i) + a2 * s2%data(d)%elts(id_i)
    end subroutine cal_lc
  end function lcf

  function residual (f, u, Lu, l)
    ! Calculates f - Lu(u) at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: f, residual, u

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    residual = lcf (1d0, f, -1d0, Lu (u, l), l)
    residual%bdry_uptodate = .false.
  end function residual

  subroutine restrict (scaling, coarse)
    ! Restriction operator for scalars
    implicit none
    integer                   :: coarse
    type(Float_Field), target :: scaling

    integer :: d, l

    do d = 1, size(grid)
       scalar => scaling%data(d)%elts
       call apply_interscale_d (Restrict_full_weighting, grid(d), coarse, z_null, 0, 1) ! +1 to include poles
       nullify (scalar)
    end do
    scaling%bdry_uptodate = .false.
    call update_bdry (scaling, coarse, 66)
  end subroutine restrict

  function restrict_fun (scaling, coarse)
    ! Restriction operator for scalars
    implicit none
    integer                   :: coarse
    type(Float_Field), target :: scaling
    type(Float_Field), target :: restrict_fun

    integer :: d, l

    call restrict (scaling, coarse)
    restrict_fun = scaling
  end function restrict_fun

  subroutine prolong (scaling, fine)
    ! Prolong from coarse scale fine-1 to scale fine
    use wavelet_mod
    implicit none
    integer                   :: fine
    type(Float_Field), target :: scaling, wavelet

    integer :: d, l

    call update_bdry1 (scaling, fine-1, fine, 5)

    ! Scaling scalar to finer nodes existing at coarser grid (extend) 
    do d = 1, size(grid)
       scalar => scaling%data(d)%elts
       call apply_interscale_d2 (copy_coarse, grid(d), fine-1, z_null, 0, 1)
       nullify (scalar)
    end do
    scaling%bdry_uptodate = .false.
    call update_bdry (scaling, fine, 66)

    ! Reconstruct scalar at finer nodes not existing at coarser grid by interpolation
    do d = 1, size(grid)
       scalar => scaling%data(d)%elts
       call apply_interscale_d (interpolate, grid(d), fine-1, z_null, 0, 0)
       nullify (scalar)
    end do
    scaling%bdry_uptodate = .false.
  end subroutine prolong

  function prolong_fun (scaling, fine)
    ! Prolong from coarse scale fine-1 to scale fine
    use wavelet_mod
    implicit none
    integer                   :: fine
    type(Float_Field), target :: prolong_fun, scaling

    integer :: d, l

    prolong_fun = scaling

    call prolong (prolong_fun, fine)
  end function prolong_fun

  subroutine copy_coarse (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Sub-sample to prolong coarse points to fine grid
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd
    integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE

    id_par = idx (i_par, j_par, offs_par, dims_par) + 1
    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) + 1

    if (dom%mask_n%elts(id_chd) == FROZEN) return ! FROZEN mask -> do not overide with wrong value
    
    scalar(id_chd) = scalar(id_par)
  end subroutine copy_coarse

  subroutine interpolate (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Reconstruct scalars at fine nodes not existing at coarse scale by interpolation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, idN_chd, idE_chd, idNE_chd, id2N_chd, id2E_chd, id2S_chd, id2W_chd, id2NE_chd
    integer :: id_par, idE_par, idNE_par, idN_par

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd)

    idN_chd   = idx (i_chd,   j_chd+1, offs_chd, dims_chd)
    idE_chd   = idx (i_chd+1, j_chd,   offs_chd, dims_chd)
    idNE_chd  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd)
    id2N_chd  = idx (i_chd,   j_chd+2, offs_chd, dims_chd)
    id2E_chd  = idx (i_chd+2, j_chd,   offs_chd, dims_chd)
    id2S_chd  = idx (i_chd,   j_chd-2, offs_chd, dims_chd)
    id2W_chd  = idx (i_chd-2, j_chd,   offs_chd, dims_chd)
    id2NE_chd = idx (i_chd+2, j_chd+2, offs_chd, dims_chd)


    ! Interpolate scalars to reconstruct values at fine scale
    scalar(idE_chd +1) = Interp_node (dom,  idE_chd,    id_chd, id2E_chd, id2NE_chd,  id2S_chd)
    scalar(idNE_chd+1) = Interp_node (dom, idNE_chd, id2NE_chd,   id_chd,  id2E_chd,  id2N_chd) 
    scalar(idN_chd +1) = Interp_node (dom,  idN_chd,    id_chd, id2N_chd,  id2W_chd, id2NE_chd)
  end subroutine interpolate

  subroutine wlt_interpolate (scaling, wavelet, scale_in, scale_out)
    ! Use forward/inverse wavelet transformation for restriction and prolongation
    ! note: must have previously computed the appropriate wavelet coeffcients to prolong
    implicit none
    integer,                   intent (in)    :: scale_in, scale_out
    type(Float_Field), target, intent (inout) :: scaling, wavelet

    if (scale_in > scale_out) then     ! restriction
       call forward_scalar_transform (scaling, wavelet, scale_out, scale_in)
    elseif (scale_in < scale_out) then ! prolongation
       call inverse_scalar_transform (wavelet, scaling, scale_in, scale_out)
    else                               ! do nothing
    end if
  end subroutine wlt_interpolate

  function wlt_int (scaling, wavelet, scale_in, scale_out)
    ! Use forward/inverse wavelet transformation for restriction and prolongation
    ! note: must have previously computed the appropriate wavelet coeffcients to prolong
    implicit none
    type(Float_Field), target                 :: wlt_int
    integer,                   intent (in)    :: scale_in, scale_out
    type(Float_Field), target, intent (inout) :: scaling, wavelet

    wlt_int = scaling
    if (scale_in > scale_out) then     ! restriction
       call forward_scalar_transform (wlt_int, wavelet, scale_out, scale_in)
    elseif (scale_in < scale_out) then ! prolongation
       call inverse_scalar_transform (wavelet, wlt_int, scale_in, scale_out)
    else                               ! do nothing
    end if
  end function wlt_int

  subroutine FMG (u, f, Lu, Lu_diag)
    ! Solves linear equation L(u) = f using the full multigrid algorithm with V-cycles
    use adapt_mod
    implicit none
    type(Float_Field), target :: f, u

    integer                                       :: iter, j, l
    integer, dimension(level_start:level_end)     :: iterations
    integer                                       :: max_vcycle
    real(8)                                       :: vcycle_tol
    real(8), dimension(level_start:level_end)     :: nrm_f
    real(8), dimension(level_start:level_end,1:2) :: nrm_res

    interface
       function Lu (u, l)
         ! Linear operator applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         ! Diagonal of linear operator applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface

    vcycle_tol  = 1d-4
    fine_tol    = vcycle_tol
    max_vcycle  = 10
    coarse_iter = 1000
    coarse_tol  = 1d-9

    call update_bdry (f, NONE, 55)
    call update_bdry (u, NONE, 55)
 
    if (log_iter) then
       nrm_res = 0d0; iterations = 0
       do j = level_start, level_end
          nrm_f(j) = l2 (f, j); if (nrm_f(j) == 0d0) nrm_f(j) = 1d0
          nrm_res(j,1) = l2 (residual (f, u, Lu, j), j) / nrm_f(j)
       end do
    end if
    if (maxval (nrm_res(:,1)) < coarse_tol / 1d1) return

    call bicgstab (u, f, Lu, level_start, coarse_tol, coarse_iter, nrm_res(level_start,2), iterations(level_start))
    do j = level_start+1, level_end
       if (nrm_res(j,1) == 0d0) exit
       call prolong (u, j)
       do iter = 1, max_vcycle
          call v_cycle (u, f, Lu, Lu_diag, level_start, j, nrm_res(j,2))
          iterations(j) = iter
          if (nrm_res(j,2) < vcycle_tol) exit
       end do
    end do

    if (log_iter) then
       if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual    Iterations"
       do j = level_start, level_end
          if (rank == 0) write (6,'(i2,12x,2(es8.2,10x),i4)') j, nrm_res(j,:), iterations(j)
       end do
    end if
  end subroutine FMG

  subroutine v_cycle (u, f, Lu, Lu_diag, jmin, jmax, err)
    ! Solves linear equation L(u) = f using the standard multigrid algorithm with V-cycles
    implicit none
    integer                   :: jmin, jmax
    real(8)                   :: err
    type(Float_Field), target :: u, f

    integer :: j
    integer :: down_iter = 2, up_iter = 2, pre_iter = 2

    real(8) :: w1

    type(Float_Field), target :: corr, res

    interface
       function Lu (u, l)
         ! Returns result of linear operator applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         ! Returns diagonal of linear operator applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface

    w1 = 1d0

    corr = u

    ! Down V-cycle
    res = residual (f, u, Lu, jmax)
    do j = jmax, jmin+1, -1
       call zero_float_field (corr, S_MASS, j)
       call Jacobi (corr, res, Lu, Lu_diag, j, up_iter)
       res = residual (res, corr, Lu, j)
       call restrict (res, j-1)
    end do

    ! Exact solution on coarsest grid
    call zero_float_field (corr, S_MASS, jmin)
    call bicgstab (corr, res, Lu, jmin, coarse_tol, coarse_iter)

    ! Up V-cycle
    do j = jmin+1, jmax
       corr = lcf (1d0, corr, w1, prolong_fun (corr, j), j)
       call Jacobi (corr, res, Lu, Lu_diag, j, down_iter)
    end do

    ! V-cycle correction to solution
    u = lcf (1d0, u, 1d0, corr, jmax)

    ! Post-smooth to reduce zero eigenvalue error mode
    call Jacobi (u, f, Lu, Lu_diag, jmax, pre_iter, err)

    err = l2 (residual (f, u, Lu, jmax),jmax) / l2 (f, jmax)
  end subroutine v_cycle

  subroutine SJR (u, f, Lu, Lu_diag)
    ! Solves linear equation L(u) = f using a simple multiscale algorithm with Scheduled Rexation Jacobi iterations as the smoother
    implicit none
    type(Float_Field), target :: f, u

    integer                                       :: j, n
    integer, dimension(level_start:level_end)     :: iterations
    real(8)                                       :: k_max, k_min
    real(8), dimension(level_start:level_end)     :: nrm_f
    real(8), dimension(level_start:level_end,1:2) :: nrm_res

    integer, dimension(level_start:level_end) :: iter

    interface
       function Lu (u, l)
         ! Returns result of linear operator Lu applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         ! Returns diagonal of linear operator Lu applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface
    
    if (log_iter) then
       nrm_res = 0d0; iterations = 0
       do j = level_start, level_end
          nrm_f(j) = l2 (f, j); if (nrm_f(j) == 0d0) nrm_f(j) = 1d0
          nrm_res(j,1) = l2 (residual (f, u, Lu, j), j) / nrm_f(j)
       end do
    end if
    if (maxval (nrm_res(:,1)) < coarse_tol / 1d1) return

    ! Optimal Scheduled Relaxation Jacobi parameters (Adsuara, et al J Comput Phys v 332, 2016)
    ! (k_min and k_max are determined empirically to give optimal convergence on fine non uniform grids)
    k_max = 1.8d0

    ! Values for k_min optimized for J5J9
    if (fine_tol <= 1d-5) then
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

    call update_bdry (f, NONE, 55)

    j = level_start
    call bicgstab (u, f, Lu, j, coarse_tol, coarse_iter, nrm_res(j,2), iterations(j))
    do j = level_start+1, level_end
       call prolong (u, j)
       call SJR_iter (u, f, nrm_f(j), Lu, Lu_diag, j, fine_iter, nrm_res(j,2), iterations(j))
    end do
    deallocate (w)

    if (log_iter) then
       if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual    Iterations"
       do j = level_start, level_end
          if (rank == 0) write (6,'(i2,12x,2(es8.2,10x),i4)') j, nrm_res(j,:), iterations(j)
       end do
    end if
  end subroutine SJR

  subroutine SJR_iter (u, f, nrm_f, Lu, Lu_diag, l, max_iter, nrm_res, iter)
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
         ! Returns result of linear operator Lu applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         ! Returns diagonal of linear operator Lu applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface

    ! Initialize
    ii = 0
    iter = 0
    nrm_res = l2 (residual (f, u, Lu, l), l) / nrm_f

    do while (iter < max_iter)
       if (nrm_res <= fine_tol) exit
       ii = ii + 1
       if (ii > m) then
          if (nrm_res < 2d0 * fine_tol) then ! avoid starting a new SJR cycle if error is small enough
             exit
          else
             ii = 1
          end if
       end if
       iter = iter + 1
       call Jacobi_iteration (u, f, w(ii), Lu, Lu_diag, l)
       nrm_res = l2 (residual (f, u, Lu, l), l) / nrm_f
    end do
    u%bdry_uptodate = .false.
  end subroutine SJR_iter

  subroutine Jacobi (u, f, Lu, Lu_diag, l, iter_max, err_out)
    ! Damped Jacobi iterations
    implicit none
    integer                   :: l, iter_max
    real(8), optional         :: err_out

    type(Float_Field), target :: f, u

    integer            :: iter
    real(8)            :: nrm_f
    real(8), parameter :: omega = 1d0

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         ! Returns diagonal of linear operator Lu applied to u at scale l
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface

    call update_bdry (f, l, 50)
    call update_bdry (u, l, 50)

    do iter = 1, iter_max
       call Jacobi_iteration (u, f, omega, Lu, Lu_diag, l)
    end do

    nrm_f = l2 (f, l); if (nrm_f == 0d0) nrm_f = 1d0
    if (present(err_out))  err_out = l2 (residual (f, u, Lu, l), l) / nrm_f
  end subroutine Jacobi

  subroutine Jacobi_iteration (u, f, omega, Lu, Lu_diag, l)
    ! Performs a single weighted Jacobi iteration for equation Lu(u) = f
    implicit none
    integer            :: l
    real(8)            :: omega ! weight
    type(Float_Field)  :: u, f

    type(Float_Field) :: Au, diag

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
       function Lu_diag (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu_diag, u
       end function Lu_diag
    end interface

    Au = Lu (u, l)
    diag = Lu_diag (u, l)

    call apply_onescale (cal_jacobi, l, z_null, 0, 1)

    u%bdry_uptodate = .false.
    call update_bdry (u, l, 50)
  contains
    subroutine cal_jacobi (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id
 
      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1
      
      if (dom%mask_n%elts(id) >= ADJZONE) &
           u%data(d)%elts(id) =  u%data(d)%elts(id) + omega * (f%data(d)%elts(id) - Au%data(d)%elts(id)) / diag%data(d)%elts(id)
    end subroutine cal_jacobi
  end subroutine Jacobi_iteration
  
  subroutine bicgstab (u, f,  Lu, l, tol, iter_max, err_out, iter_out)
    ! Solves the linear system Lu(u) = f at scale l using bi-cgstab algorithm (van der Vorst 1992).
    ! This is a conjugate gradient type algorithm.
    implicit none
    integer                   :: l, iter_max
    integer, optional         :: iter_out
    real(8)                   :: tol
    real(8), optional         :: err_out
    type(Float_Field), target :: f, u

    integer                   :: iter
    real(8)                   :: alph, b, err, nrm_f, nrm_res0, omga, rho, rho_old
    type(Float_Field), target :: Ap, As, res, res0, p, s

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    call update_bdry (f, l, 88)
    call update_bdry (u, l, 88)

    nrm_f = l2 (f, l); if (nrm_f == 0d0) nrm_f = 1d0

    ! Initialize float fields
    res  = residual (f, u, Lu, l)
    res0 = res
    p    = res0
    s    = res0
    Ap   = Lu (p, l)
    As   = Ap

    rho  = dp (res0, res, l)
    
    do iter = 1, iter_max
       alph = rho / dp (Ap, res0, l)

       s = lcf (1d0, res, -alph, Ap, l)
       As = Lu (s, l)

       omga = dp (As, s, l) / dp (As, As, l)

       u = lcf (1d0, u, 1d0, lcf (alph, p, omga, s, l), l)

       res = lcf (1d0, s, -omga, As, l)

       err = l2 (res, l) / nrm_f
       if (err <= tol) exit

       rho_old = rho
       rho = dp (res0, res, l)

       b = (alph/omga) * (rho/rho_old)
       p = lcf (1d0, res, b, lcf (1d0, p, -omga, Ap, l), l)
       Ap = Lu (p, l)
    end do
    u%bdry_uptodate = .false.
    call update_bdry (u, l, 88)

    if (present(err_out))  err_out  = err
    if (present(iter_out)) iter_out = iter
  end subroutine bicgstab
  
  function Lu2 (Lu, u, l)
    use domain_mod
    implicit none
    integer                   :: l
    type(Float_Field), target :: Lu2, u

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    Lu2 = restrict_fun (Lu (prolong_fun (u, l+1), l+1), l)
  end function Lu2

  function elliptic_fun (u, l)
    ! Test elliptic equation
    ! Computes Laplacian(u) + u at scale l
    use ops_mod
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_fun, u

    integer :: d, j
    real(8) :: sigma

    call update_bdry (u, l, 50)

    elliptic_fun = u
    call zero_float_field (elliptic_fun, S_MASS, l)

    ! Compute scalar fluxes
    do d = 1, size(grid)
       scalar => u%data(d)%elts
       h_flux => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
       end do
       nullify (scalar, h_flux)
    end do
    horiz_flux%bdry_uptodate = .false.
    call update_bdry (horiz_flux(S_MASS), l, 11)

    ! Compute divergence of fluxes
    do d = 1, size(grid)
       dscalar => elliptic_fun%data(d)%elts
       h_flux  => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, h_flux)
    end do

    ! Add constant term
    sigma = radius/4d0
    elliptic_fun = lcf (1d0, elliptic_fun, 4d0/sigma**2, u, l)

    elliptic_fun%bdry_uptodate = .false.
    call update_bdry (elliptic_fun, l, 12)
  end function elliptic_fun
  
  function elliptic_fun_diag (q, l)
    ! Local approximation of diagonal of elliptic operator
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_fun_diag, q

    elliptic_fun_diag = q
    call apply_onescale (cal_elliptic_fun_diag, l, z_null, 0, 1)
  contains
    subroutine cal_elliptic_fun_diag  (dom, i, j, zlev, offs, dims)
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer            :: d, id, id_i, idE, idNE, idN, idW, idSW, idS
      real(8)            :: wgt
      logical, parameter :: exact = .true.

      real(8) :: sigma

      d = dom%id + 1
      id = idx (i, j, offs, dims)
      id_i = id + 1

      sigma = radius/4d0

      if (exact) then ! true local value
         idE  = idx (i+1, j,   offs, dims) 
         idNE = idx (i+1, j+1, offs, dims) 
         idN  = idx (i,   j+1, offs, dims) 
         idW  = idx (i-1, j,   offs, dims) 
         idSW = idx (i-1, j-1, offs, dims) 
         idS  = idx (i,   j-1, offs, dims)
         wgt = &
              dom%pedlen%elts(EDGE*id+RT+1)   / dom%len%elts(EDGE*id+RT+1)   + &
              dom%pedlen%elts(EDGE*id+DG+1)   / dom%len%elts(EDGE*id+DG+1)   + &
              dom%pedlen%elts(EDGE*id+UP+1)   / dom%len%elts(EDGE*id+UP+1)   + &
              dom%pedlen%elts(EDGE*idW+RT+1)  / dom%len%elts(EDGE*idW+RT+1)  + &
              dom%pedlen%elts(EDGE*idSW+DG+1) / dom%len%elts(EDGE*idSW+DG+1) + &
              dom%pedlen%elts(EDGE*idS+UP+1)  / dom%len%elts(EDGE*idS+UP+1)
      else ! average value (error less than about 5%)
         wgt = 2d0 * sqrt (3d0)
      end if
      elliptic_fun_diag%data(d)%elts(id_i) = - wgt * dom%areas%elts(id_i)%hex_inv + 4d0/sigma**2
    end subroutine cal_elliptic_fun_diag
  end function elliptic_fun_diag

  real(8) function relative_error (u, l)
    implicit none
    integer                   :: l
    type(Float_Field), target :: u

    real(8)                   :: nrm_err, nrm_sol
    type(Float_Field), target :: err

    err = u
    call apply_onescale (cal_err, l, z_null, 0, 1)

    ! Compute relative error
    nrm_sol = l2 (u,   l)
    nrm_err = l2 (err, l)
    relative_error = nrm_err / nrm_sol
  contains
    subroutine cal_err (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1

      id = idx (i, j, offs, dims) 
      if (dom%mask_n%elts(id+1) >= ADJZONE)  err%data(d)%elts(id+1) = abs (u%data(d)%elts(id+1) -  exact_sol (dom%node%elts(id+1)))
    end subroutine cal_err
  end function relative_error

  real(8) function exact_sol (p)
    implicit none
    type (Coord) :: p

    real(8) :: r, sigma

    sigma = radius/4d0

    r = dist (p, sph2cart (0d0, 0d0))

    exact_sol = exp (-(r/sigma)**2) 
  end function exact_sol
end module lin_solve_mod
