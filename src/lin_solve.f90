module lin_solve_mod
  ! Module providing two adaptive multi-grid linear equation solvers: FMG (full multi-grid) and SRJ (Scheduled Relaxation Jacobi)
  use ops_mod
  use comm_mpi_mod
  use wavelet_mod
  implicit none
  real(8)                            :: dp_loc, l2_loc
  real(8)                            :: s_test 
  logical                            :: test_elliptic = .false.
contains
  subroutine FMG (u, f, Lu, Lu_diag)
    ! Solves linear equation L(u) = f using the full multi-grid (FMG) algorithm with V-cycles
    implicit none
    type(Float_Field) :: f, u

    integer                                       :: max_vcycle = 10
    integer                                       :: down_iter = 1, up_iter = 2, post_iter = 4, pre_iter = 2
    
    integer                                       :: iter, l
    integer, dimension(level_start:level_end)     :: iterations
    real(8)                                       :: cpu, err, rel_err, t0
    real(8), dimension(level_start:level_end)     :: nrm_f
    real(8), dimension(level_start:level_end,1:2) :: nrm_res

    type(Float_Field) :: res
    
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

    call update_bdry (f, l)
    call update_bdry (u, l)

    res = u; call zero_float_field (res, AT_NODE)
    call zero_float_field (wav_coeff(S_MASS,1), AT_NODE)
    
    iterations = 0
    do l = level_end, level_start, -1
       nrm_f(l) = l2 (f, l); if (nrm_f(l) == 0d0) nrm_f(l) = 1d0
       call res_err (f, u, Lu, nrm_f(l), l, nrm_res(l,1))
    end do
    nrm_res(:,2) = nrm_res(:,1)

    if (log_iter) t0 = MPI_Wtime ()  

    ! Full multigrid iterations
    l = level_start
    call bicgstab (u, f, Lu, l, coarse_tol, coarse_iter, nrm_res(l,2), iterations(l))
    do l = l+1, level_end
       if (nrm_res(l,1) < fine_tol) cycle

       call prolong (u, l)
       
       call Jacobi (u, f, Lu, Lu_diag, l, pre_iter) ! pre-smooth to reduce zero eigenvalue error mode
       
       do iter = 1, max_vcycle
          call v_cycle
          
          iterations(l) = iter; if (nrm_res(l,2) <= 2d0 * fine_tol) exit
       end do
    end do

    if (log_iter) then
       cpu = MPI_Wtime () - t0
       if (rank==0) write (6,'(/,a,es10.4,a,/)') "FMG solver CPU time = ", cpu, "s"
       if (test_elliptic) then
          if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual    Relative error      Iterations"
          do l = level_start, level_end
             rel_err = relative_error (u, l)
             if (rank == 0) write (6,'(i2,12x,3(es8.2,10x),i4)') l, nrm_res(l,:), rel_err, iterations(l)
          end do
       else
          if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual     Iterations"
          do l = level_start, level_end
             if (rank == 0) write (6,'(i2,12x,2(es8.2,10x),i4)') l, nrm_res(l,:), iterations(l)
          end do
       end if
    end if
  contains
    subroutine v_cycle 
      ! Standard V-cycle iterations
      implicit none
      integer           :: j
      type(Float_Field) :: corr

      ! Initialize
      corr = u
      do j = level_start, l
         call zero_float_field (corr, AT_NODE, level_start, l)
      end do

      ! Down V-cycle
      call residual (f, u, Lu, l, res)
      do j = l, level_start+1, -1
         call Jacobi (corr, res, Lu, Lu_diag, j, down_iter)
         call residual (res, corr, Lu, j, res)
         call restrict (res, j-1)
      end do

      ! Coarsest scale
      call bicgstab (corr, res, Lu, level_start, coarse_tol, coarse_iter) ! exact solution on coarsest grid

      ! Up V-cycle
      do j = level_start+1, l
         call lc (corr, 1d0, corr, 1d0, prol_fun(corr,j), j)
         call Jacobi (corr, res, Lu, Lu_diag, j, up_iter)                 ! smoother
      end do

      call lc (u, 1d0, u, 1d0, corr, l)                                   ! V-cycle correction to solution
      call Jacobi (u, f, Lu, Lu_diag, l, post_iter)                       ! post-smooth to reduce zero eigenvalue error mode
      call res_err (f, u, Lu, nrm_f(l), l, nrm_res(l,2))                  ! normalized residual error
    end subroutine v_cycle
  end subroutine FMG

  subroutine SRJ (u, f, Lu, Lu_diag)
    ! Solves linear equation L(u) = f using a simple multiscale algorithm with scheduled relaxation Jacobi (SRJ) iterations (Adsuara et al J Comput Phys v 332, 2017)
    ! parameters m, k_min and k_max were chosen empirically to give optimal convergence on fine non uniform grids
    implicit none
    type(Float_Field) :: f, u

    integer, parameter                            :: m     = 10
    real(8)                                       :: k_min = 5d-2   ! empirically optimized, 0 < k_min <= k_max
    real(8)                                       :: k_max = 2d0
    
    integer                                       :: l, n
    integer, dimension(level_start:level_end)     :: iterations
    real(8)                                       :: cpu, t0
    real(8), dimension(1:m)                       :: w
    real(8), dimension(level_start:level_end)     :: nrm_f
    real(8), dimension(level_start:level_end,1:2) :: nrm_res

    interface
       function Lu (u, l)
         ! Linear operator Lu applied to u at scale l
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

    call update_bdry (f, l)
    call update_bdry (u, l)

    call zero_float_field (wav_coeff(S_MASS,1), AT_NODE)

    iterations = 0
    do l = level_start, level_end
       nrm_f(l) = l2 (f, l); if (nrm_f(l) == 0d0) nrm_f(l) = 1d0
       call res_err (f, u, Lu, nrm_f(l), l, nrm_res(l,1))
    end do

    if (log_iter) t0 = MPI_Wtime ()  

    if (fine_iter < m) then ! do not use SRJ
       w = 1d0
    else 
       fine_iter = m * (fine_iter / m) ! ensure that fine_iter is an integer multiple of m
       do n = 1, m
          w(n) = 2d0 / (k_min + k_max - (k_max - k_min) * cos (MATH_PI * (2d0*dble(n) - 1d0)/(2d0*dble(m))))
       end do
    end if

    l = level_start
    call bicgstab (u, f, Lu, l, coarse_tol, coarse_iter, nrm_res(l,2), iterations(l))
    do l = l+1, level_end
       call prolong (u, l)
       call SRJ_iter
    end do

    if (log_iter) then
       cpu = MPI_Wtime () - t0; if (rank==0) write (6,'(/,a,es10.4,a,/)') "SJR solver CPU time = ", cpu, "s"
       if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual    Iterations"
       do l = level_start, level_end
          if (rank == 0) write (6,'(i2,12x,2(es8.2,10x),i4)') l, nrm_res(l,:), iterations(l)
       end do
    end if
  contains
    subroutine SRJ_iter
      ! Jacobi iterations for smoothing multigrid iterations
      ! uses Scheduled Relaxation Jacobi (SRJ) iterations (Yang and Mittal JCP 274, 2014)
      implicit none
      integer :: ii, iter

      ii = 0
      do iter = 1, fine_iter
         ii = ii + 1
         if (ii > m) then
            if (nrm_res(l,2) < fine_tol) then ! avoid starting a new SJR cycle if error is small enough
               exit
            else
               ii = 1
            end if
         end if

         call Jacobi_iteration (u, f, w(ii), Lu, Lu_diag, l)
         
         call res_err (f, u, Lu, nrm_f(l), l, nrm_res(l,2))
         iterations(l) = iter; if (nrm_res(l,2) <= 2d0 * fine_tol) exit
      end do
      u%bdry_uptodate = .false.
    end subroutine SRJ_iter
  end subroutine SRJ

  subroutine restrict (scaling, coarse)
    ! Wavelet restriction
    implicit none
    integer           :: coarse
    type(Float_Field) :: scaling

    call forward_scalar_transform (scaling, wav_coeff(S_TEMP,1), coarse, coarse+1)
    scaling%bdry_uptodate = .false.
  end subroutine restrict
  
  subroutine prolong (scaling, fine)
    ! Prolong from coarse scale fine-1 to scale fine
    use wavelet_mod
    implicit none
    integer           :: fine
    type(Float_Field) :: scaling

    call inverse_scalar_transform (wav_coeff(S_MASS,1), scaling, fine-1, fine)

    scaling%bdry_uptodate = .false.
  end subroutine prolong

  function prol_fun (scaling, fine)
    ! Prolong from coarse scale fine-1 to scale fine
    use wavelet_mod
    implicit none
    integer           :: fine
    type(Float_Field) :: scaling
    type(Float_Field) :: prol_fun

    call zero_float_field (wav_coeff(S_TEMP,1), AT_NODE, fine)

    prol_fun = scaling
    call inverse_scalar_transform (wav_coeff(S_TEMP,1), prol_fun, fine-1, fine)

    prol_fun%bdry_uptodate = .false.
  end function prol_fun

  subroutine Jacobi (u, f, Lu, Lu_diag, l, iter_max)
    ! Damped Jacobi iterations
    implicit none
    integer  :: l, iter_max

    type(Float_Field) :: f, u

    integer            :: iter
    real(8), parameter :: omega = 1d0

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

    call update_bdry (f, l)

    do iter = 1, iter_max
       call Jacobi_iteration (u, f, omega, Lu, Lu_diag, l)
    end do
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

    Au   = Lu      (u, l)
    diag = Lu_diag (u, l)

    call apply_onescale (cal_jacobi, l, z_null, 0, 1)

    u%bdry_uptodate = .false.
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
    integer           :: l, iter_max
    integer, optional :: iter_out
    real(8)           :: tol
    real(8), optional :: err_out
    type(Float_Field) :: f, u

    integer           :: iter
    real(8)           :: alph, b, err, nrm_f, nrm_res0, omga, rho, rho_old
    type(Float_Field) :: Ap, As, corr, res, res0, p, s

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    call update_bdry (f, l)
    call update_bdry (u, l)

    nrm_f = l2 (f, l); if (nrm_f == 0d0) nrm_f = 1d0

    res  = u; call zero_float_field (res,  AT_NODE, lmin=level_start)
    corr = u; call zero_float_field (corr, AT_NODE, lmin=level_start)

    ! Initialize float fields
    call residual (f, u, Lu, l, res)
    
    res0 = res
    p    = res0
    s    = res0
    Ap   = Lu (p, l)
    As   = Ap

    rho  = dp (res0, res, l)
    
    do iter = 1, iter_max
       alph = rho / dp (Ap, res0, l)

       call lc (s, 1d0, res, -alph, Ap, l)
       As = Lu (s, l)

       omga = dp (As, s, l) / dp (As, As, l)

       call lc (corr, alph, p, omga, s, l)
       call lc (u, 1d0, u, 1d0, corr, l)

       call lc (res, 1d0, s, -omga, As, l)

       err = l2 (res, l) / nrm_f
       if (err <= tol) exit

       rho_old = rho
       rho = dp (res0, res, l)
       
       b = (alph/omga) * (rho/rho_old)

       call lc (corr, 1d0, p, -omga, Ap, l)
       call lc (p, 1d0, res, b, corr, l)
       Ap = Lu (p, l)
    end do
    u%bdry_uptodate = .false.

    if (present(err_out))  err_out  = err
    if (present(iter_out)) iter_out = iter
  end subroutine bicgstab

  real(8) function l2 (s, l)
    ! Returns l_2 norm of scalar s at scale l
    implicit none
    integer,           intent(in) :: l
    type(Float_Field), intent(in) :: s

    call update_bdry (s, l)

    l2_loc = 0d0

    call apply_onescale (cal_l2, l, z_null, 0, 1)

    l2 = sqrt (sum_real (l2_loc))
  contains
    subroutine cal_l2 (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      l2_loc = l2_loc +  s%data(d)%elts(id)**2
    end subroutine cal_l2
  end function l2

  real(8) function dp (s1, s2, l)
    ! Calculates dot product of s1 and s2 at scale l
    implicit none
    integer,           intent(in) :: l
    type(Float_Field), intent(in) :: s1, s2

    call update_bdry (s1, l)
    call update_bdry (s2, l)

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

  subroutine lc (s, a1, s1, a2, s2, l)
    ! Calculates linear combination of scalars s = a1*s1 + a2*s2 at scale l
    implicit none
    integer,           intent(in)    :: l
    real(8),           intent(in)    :: a1, a2
    type(Float_Field), intent(inout) :: s
    type(Float_Field), intent(in)    :: s1, s2

    call update_bdry (s1, l)
    call update_bdry (s2, l)
    
    call apply_onescale (cal_lc, l, z_null, 0, 1)

    s%bdry_uptodate = .false.
  contains
    subroutine cal_lc (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      s%data(d)%elts(id) = a1 * s1%data(d)%elts(id) + a2 * s2%data(d)%elts(id)
    end subroutine cal_lc
  end subroutine lc

  subroutine residual (f, u, Lu, l, res)
    ! Residual f - Lu(u) at scale l
    implicit none
    integer,           intent(in)    :: l
    type(Float_Field), intent(in)    :: f, u
    type(Float_Field), intent(inout) :: res
    
    type(Float_Field) :: Au

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    call update_bdry (f, l)

    Au = Lu (u, l)

    call apply_onescale (cal_res, l, z_null, 0, 1)

    res%bdry_uptodate = .false.
  contains
    subroutine cal_res (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      if (dom%mask_n%elts(id) >= ADJZONE) res%data(d)%elts(id) = f%data(d)%elts(id) - Au%data(d)%elts(id)
    end subroutine cal_res
  end subroutine residual

  subroutine res_err (f, u, Lu, nrm_f, l, err)
    ! Normalized residual error ||f - Lu(u)||/||f|| at scale l
    implicit none
    integer,           intent(in)    :: l
    real(8),           intent(in)    :: nrm_f
    real(8),           intent(out)   :: err
    type(Float_Field), intent(in)    :: f, u
    
    type(Float_Field) :: res

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface
    
    res = u; call zero_float_field (res, AT_NODE, l)

    call residual (f, u, Lu, l, res)

    err = l2 (res, l) / nrm_f
  end subroutine res_err
  
  function elliptic_fun (u, l)
    ! Test elliptic equation
    !
    ! L(u) = Laplacian (u) - 10/s_test^2 u
    use ops_mod
    implicit none
    integer                   :: l
    type(Float_Field), target ::  u
    type(Float_Field), target :: elliptic_fun

    integer :: d, j

    call update_bdry (u, l)

    elliptic_fun = u; call zero_float_field (elliptic_fun, AT_NODE, lmin=l)

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
    call update_bdry (horiz_flux(S_MASS), l)

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
    call lc (elliptic_fun, 1d0, elliptic_fun, -10d0/s_test**2, u, l)

    elliptic_fun%bdry_uptodate = .false.
    call update_bdry (elliptic_fun, l)
  end function elliptic_fun
  
  function elliptic_fun_diag (q, l)
    ! Local approximation of diagonal of elliptic operator
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_fun_diag, q

    elliptic_fun_diag = q
    call apply_onescale (cal_elliptic_fun_diag, l, z_null, 0, 1)

    elliptic_fun_diag%bdry_uptodate = .false.
    call update_bdry (elliptic_fun_diag, l)
  contains
    subroutine cal_elliptic_fun_diag  (dom, i, j, zlev, offs, dims)
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer            :: d, id, id_i, idE, idNE, idN, idW, idSW, idS
      real(8)            :: wgt
      logical, parameter :: exact = .true.

      d = dom%id + 1
      id = idx (i, j, offs, dims)
      id_i = id + 1

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
      elliptic_fun_diag%data(d)%elts(id_i) = - wgt * dom%areas%elts(id_i)%hex_inv - 10d0/s_test**2
    end subroutine cal_elliptic_fun_diag
  end function elliptic_fun_diag

  real(8) function relative_error (u, l)
    ! Relative error of test elliptic problem
    implicit none
    integer           :: l
    type(Float_Field) :: u

    real(8)           :: nrm_err, nrm_sol
    type(Float_Field) :: err

    call update_bdry (u, l)

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
      id = idx (i, j, offs, dims) + 1

      err%data(d)%elts(id) = u%data(d)%elts(id) -  exact_sol (dom%node%elts(id))
    end subroutine cal_err
  end function relative_error

  real(8) function exact_sol (p)
    implicit none
    type (Coord) :: p

    real(8) :: r

    r = geodesic (p, sph2cart (0d0, 0d0))

    exact_sol = s_test**2 * exp (-(r/s_test)**2) 
  end function exact_sol
end module lin_solve_mod
