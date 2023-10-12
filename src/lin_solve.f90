module lin_solve_mod
  use ops_mod
  use comm_mpi_mod
  use wavelet_mod
  implicit none
  integer                            :: ii, m
  real(8)                            :: dp_loc, linf_loc, l2_loc
  real(8), pointer                   :: mu1, mu2
  real(8), dimension(:), pointer     :: scalar2, scalar3
  real(8), dimension(:), allocatable :: w
  type(Float_Field), target          :: float_var
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

    if (dom%mask_n%elts(id) >= ZERO) linf_loc = max (linf_loc, abs (scalar(id)))
  end subroutine cal_linf_scalar

  real(8) function l2 (s, l)
    ! Returns l_2 norm of scalar s at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: s
    
    integer :: d, j

    l2_loc = 0d0
    do d = 1, size(grid)
       scalar => s%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_l2_scalar, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
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

    call update_bdry (s1, l, 50)
    call update_bdry (s2, l, 50)

    dp_loc = 0d0
    do d = 1, size(grid)
       scalar => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_dotproduct, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
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
  
  subroutine lc (a1, s1, a2, s2, s3, l)
    ! Calculates linear combination of scalars s3 = a1*s1 + a*s2 at scale l
    implicit none
    integer                   :: l
    real(8), target           :: a1, a2
    type(Float_Field), target :: s1, s2, s3

    integer :: d, j

    do d = 1, size(grid)
       scalar  => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       scalar3 => s3%data(d)%elts
       mu1     => a1
       mu2     => a2
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_lc, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar, scalar2, scalar3, mu1, mu2)
    end do
    s1%bdry_uptodate = .false.
    call update_bdry (s1, l, 80)
  end subroutine lc

  function lcf (a1, s1, a2, s2, l)
    ! Calculates linear combination of scalars lcf = a1*s1 + a2*s2 at scale l
    implicit none
    integer                   :: l
    real(8), target           :: a1, a2
    type(Float_Field), target :: lcf, s1, s2

    integer :: d, j
    
    lcf = s1
    do d = 1, size(grid)
       scalar  => s1%data(d)%elts
       scalar2 => s2%data(d)%elts
       scalar3 => lcf%data(d)%elts
       mu1     => a1
       mu2     => a2

       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_lc, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
    end do
    nullify (scalar, scalar2, scalar3, mu1, mu2)
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
    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar3(id_i) = mu1*scalar(id_i) + mu2 * scalar2(id_i)
  end subroutine cal_lc

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
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_lc2, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar, scalar2, scalar3, mu1, mu2)
    end do
    u%bdry_uptodate = .false.
    call update_bdry (u, l, 84)
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

  subroutine restrict (scaling, coarse)
    ! Restriction operator for scalars
    implicit none
    integer                   :: coarse
    type(Float_Field), target :: scaling

    integer :: d, l

    do d = 1, size(grid)
       scalar => scaling%data(d)%elts
       !call apply_interscale_d (cal_restriction, grid(d), coarse, z_null, 0, 1) ! +1 to include poles
       call apply_interscale_d (Restrict_full_weighting, grid(d), coarse, z_null, 0, 1) ! +1 to include poles
       nullify (scalar)
    end do
    scaling%bdry_uptodate = .false.
    call update_bdry (scaling, coarse, 66)
  end subroutine restrict
  
  subroutine cal_restriction (dom, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Restriction on hexagons
    ! (accounts for non-nested hexagons at coarse and fine grids)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_chd, id_par
    integer :: idE, idNE, idN2E, id2NE, idN, idW, idNW, idS2W, idSW, idS, id2SW, idSE

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) + 1

    if (dom%mask_n%elts(id_chd) == 0) return

    id_par = idx (i_par, j_par, offs_par, dims_par) + 1

    idE   = idx (i_chd+1, j_chd,   offs_chd, dims_chd) + 1
    idNE  = idx (i_chd+1, j_chd+1, offs_chd, dims_chd) + 1
    idN2E = idx (i_chd+2, j_chd+1, offs_chd, dims_chd) + 1
    id2NE = idx (i_chd+1, j_chd+2, offs_chd, dims_chd) + 1
    idN   = idx (i_chd,   j_chd+1, offs_chd, dims_chd) + 1
    idW   = idx (i_chd-1, j_chd,   offs_chd, dims_chd) + 1
    idNW  = idx (i_chd-1, j_chd+1, offs_chd, dims_chd) + 1
    idS2W = idx (i_chd-2, j_chd-1, offs_chd, dims_chd) + 1
    idSW  = idx (i_chd-1, j_chd-1, offs_chd, dims_chd) + 1
    idS   = idx (i_chd,   j_chd-1, offs_chd, dims_chd) + 1
    id2SW = idx (i_chd-1, j_chd-2, offs_chd, dims_chd) + 1
    idSE  = idx (i_chd+1, j_chd-1, offs_chd, dims_chd) + 1

    scalar(id_par) = ( &
         scalar(id_chd) / dom%areas%elts(id_chd)%hex_inv    + &
         scalar(idE)    * dom%overl_areas%elts(idE  )%a(1)  + &
         scalar(idNE)   * dom%overl_areas%elts(idNE )%a(2)  + &
         scalar(idN2E)  * dom%overl_areas%elts(idN2E)%a(3)  + &
         scalar(id2NE)  * dom%overl_areas%elts(id2NE)%a(4)  + &
         scalar(idN)    * dom%overl_areas%elts(idN  )%a(1)  + &
         scalar(idW)    * dom%overl_areas%elts(idW  )%a(2)  + &
         scalar(idNW)   * dom%overl_areas%elts(idNW )%a(3)  + &
         scalar(idS2W)  * dom%overl_areas%elts(idS2W)%a(4)  + &
         scalar(idSW)   * dom%overl_areas%elts(idSW )%a(1)  + &
         scalar(idS)    * dom%overl_areas%elts(idS  )%a(2)  + &
         scalar(id2SW)  * dom%overl_areas%elts(id2SW)%a(3)  + &
         scalar(idSE)   * dom%overl_areas%elts(idSE )%a(4) ) * dom%areas%elts(id_par)%hex_inv
  end subroutine cal_restriction

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

  subroutine cpt_residual (f, u, Lu, res, l)
    ! Restrict residual from scale l+1 if possible, otherwise compute at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: f, res, u

    integer                     :: d, j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface
    
    float_var = residual (f, u, Lu, l)

    ! Determine which residual values can be obtained by scalar restriction
    do d = 1, size(grid)
       scalar  => res%data(d)%elts       ! residual from level l+1 (restricted values only)
       scalar2 => float_var%data(d)%elts ! residual computed at level l
       do j = 1, grid(d)%lev(l)%length
          p_par = grid(d)%lev(l)%elts(j)
          restrict = .false.
          do c = 1, N_CHDRN
             p_chd = grid(d)%patch%elts(p_par+1)%children(c)
             if (p_chd > 0) restrict(c) = .true.
          end do
          do c = 1, N_CHDRN
             if (restrict(c)) then
                call apply_interscale_to_patch3 (residual_cpt, grid(d), p_par, c, z_null, 0, 1)
             end if
          end do
       end do
       nullify (scalar, scalar2)
    end do
  end subroutine cpt_residual

  subroutine residual_cpt (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer ::  id_par

    id_par = idx (i_par, j_par, offs_par, dims_par)

    if (dom%mask_n%elts(id_par+1) < RESTRCT) scalar(id_par) = scalar2(id_par)
  end subroutine residual_cpt

  subroutine multigrid (u, f, Lu, Lu_diag)
    ! Solves linear equation L(u) = f using the standard multigrid algorithm with V-cycles
    implicit none
    type(Float_Field), target :: f, u

    integer,  parameter                  :: vtype = 1, vcycle_iter = 10
    integer                              :: iter, j
    real(8), allocatable, dimension(:)   :: nrm_f
    real(8), allocatable, dimension(:,:) :: nrm_res

    character(*), parameter  :: type = "FMG" ! "FMG" or "V-cycle"
    !character(*), parameter  :: type = "V-cycle" ! "FMG" or "V-cycle"
    
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

    log_iter = .true.
   
    call update_bdry (f, NONE, 55)
    call update_bdry (u, NONE, 55)

    if (log_iter) then
       allocate (nrm_f(level_start:level_end), nrm_res(level_start:level_end,1:2))
       do j = level_start, level_end
          nrm_f(j) = l2 (f, j); if (nrm_f(j) == 0d0) nrm_f(j) = 1d0
          nrm_res(j,1) = l2 (residual (f, u, Lu, j), j) / nrm_f(j)
       end do
    end if

    select case (type)
    case ("FMG") ! full multigrid
       call bicgstab (u, f, Lu, level_start, coarse_iter)
       do j = level_start+1, level_end
          call prolong (u, j)
          do iter = 1, vcycle_iter
             call v_cycle (u, f, Lu, Lu_diag, level_start, j, vtype)
          end do
       end do
    case ("V-cycle") ! V-cycle
       do iter = 1, vcycle_iter
          call v_cycle (u, f, Lu, Lu_diag, level_start, level_end, vtype)
       end do
    end select
    
    if (log_iter) then
       if (rank == 0) write (6,'(a)') "Scale     Initial residual   Final residual"
       do j = level_start, level_end
          nrm_res(j,2) = l2 (residual (f, u, Lu, j), j) / nrm_f(j)
          if (rank == 0) write (6,'(i2,12x,2(es8.2,10x))') j, nrm_res(j,:)
       end do
       deallocate (nrm_res, nrm_f)
    end if
  end subroutine multigrid

  subroutine v_cycle (u, f, Lu, Lu_diag, jmin, jmax, vtype)
    ! Solves linear equation L(u) = f using the standard multigrid algorithm with V-cycles
    implicit none
    integer                   :: vtype, jmin, jmax
    type(Float_Field), target :: u, f

    integer                       :: j
    integer, parameter            :: pre_iter = 2, down_iter = 2, up_iter = 2
    real(8), parameter            :: omega = 1d0
    
    type(Float_Field), target     :: corr, res
    
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

    ! Pre-smooth
    call Jacobi (u, f, 1d0, Lu, Lu_diag, jmax, pre_iter)
    
    corr = u

    ! Down V-cycle
    res = residual (f, u, Lu, jmax)
    do j = jmax, jmin+1, -1
       call zero_float_field (corr, S_MASS, j)
       call Jacobi (corr, res, omega, Lu, Lu_diag, j, down_iter)

       res = residual (res, corr, Lu, j)
       call restrict (res, j-1)
    end do

    ! Exact solution on coarsest grid
    call zero_float_field (corr, S_MASS, jmin)
    call bicgstab (corr, res, Lu, jmin, coarse_iter)

    ! Up V-cycle
    do j = jmin+1, jmax
       corr = lcf (1d0, corr, 1d0, prolong_fun(corr,j), j)
       call Jacobi (corr, res, omega, Lu, Lu_diag, j, up_iter)
    end do

    ! V-cycle correction to solution
    u = lcf (1d0, u, 1d0, corr, jmax)
    
    ! Post-smooth to reduce zero eigenvalue error mode
    call Jacobi (u, f, omega, Lu, Lu_diag, jmax, pre_iter)
  end subroutine v_cycle

  subroutine multiscale (u, f, Lu, Lu_diag)
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

    log_iter = .true.
    
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

    call update_bdry (f, NONE, 55)
    
    do l = level_start, level_end
       nrm_f(l) = l2 (f, l)
       if (nrm_f(l) < tol_elliptic) nrm_f(l) = 1d0
       if (log_iter) r_error(1,l) = l2 (residual (f, u, Lu, l), l) / nrm_f(l)
    end do

    l = level_start
    call bicgstab (u, f, Lu, l, coarse_iter, r_error(2,l), iter(l))
    do l = level_start+1, level_end
       call prolong (u, l)
       call SJR (u, f, nrm_f(l), Lu, Lu_diag, l, fine_iter, r_error(2,l), iter(l))
    end do

    if (log_iter) then
       do l = level_start, level_end
          if (rank == 0) write (6, '("residual at scale ", i2, " = ", 2(es10.4,1x)," after ", i6, " iterations")') &
               l, r_error(:,l), iter(l)
       end do
    end if

    deallocate (w)
  end subroutine multiscale

  subroutine SJR (u, f, nrm_f, Lu, Lu_diag, l, max_iter, nrm_res, iter)
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
       call Jacobi_iteration (u, f, w(ii), Lu, Lu_diag, res, l)
       res = residual (f, u, Lu, l)
       nrm_res = l2 (res, l) / nrm_f
    end do
    u%bdry_uptodate = .false.
  end subroutine SJR

  subroutine Jacobi (u, f, omega, Lu, Lu_diag, l, max_iter)
    ! Damped Jacobi iterations
    implicit none
    integer                   :: l, max_iter
    real(8)                   :: omega
    type(Float_Field), target :: f, u

    integer                   :: iter
    type(Float_Field), target :: res

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

    res = residual (f, u, Lu, l)
    do iter = 1, max_iter
       call Jacobi_iteration (u, f, omega, Lu, Lu_diag, res, l)
    end do
    u%bdry_uptodate = .false.
  end subroutine Jacobi

  subroutine Jacobi_iteration (u, f, omega, Lu, Lu_diag, res, l)
    ! Performs a single weighted Jacobi iteration for equation Lu(u) = f
    implicit none
    integer                   :: l
    real(8), target           :: omega ! weight
    type(Float_Field), target :: u, f, res

    integer                   :: d, j
    type(Float_Field), target :: diag

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

    diag = Lu_diag (u, l)

    do d = 1, size(grid)
       mu1     => omega
       scalar  => u%data(d)%elts
       scalar2 => res%data(d)%elts
       scalar3 => diag%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_jacobi, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (mu1, scalar, scalar2, scalar3)
    end do

    res = residual (f, u, Lu, l)
  end subroutine Jacobi_iteration
  
  subroutine cal_jacobi (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
    scalar(id) = scalar(id) + mu1 * scalar2(id) / scalar3(id)
  end subroutine cal_jacobi

  subroutine bicgstab (u, f,  Lu, l, max_iter, nrm_res, total_iter)
    ! Solves the linear system Lu(u) = f at scale l using bi-cgstab algorithm (van der Vorst 1992).
    ! This is a conjugate gradient type algorithm.
    implicit none
    integer                   :: iter, l, max_iter
    integer, optional         :: total_iter
    real(8), optional         :: nrm_res
    type(Float_Field), target :: f, u
    
    real(8)                   :: alph, b, nrm_f, nrm_res0, omga, rho, rho_old
    type(Float_Field), target :: Ap, As, res, res0, p, s

    interface
       function Lu (u, l)
         use domain_mod
         implicit none
         integer                   :: l
         type(Float_Field), target :: Lu, u
       end function Lu
    end interface

    nrm_f = l2 (f, l); if (nrm_f == 0d0) nrm_f = 1d0

    ! Initialize float fields
    res  = residual (f, u, Lu, l)
    res0 = res
    p    = res0
    s    = res0
    Ap   = Lu (p, l)
    As   = Ap
    
    rho  = dp (res0, res, l)
    
    iter = 0
    do while (iter < max_iter)
       alph = rho / dp (Ap, res0, l)

       call equals_float_field (s, lcf (1d0, res, -alph, Ap, l), S_MASS, l)
       call equals_float_field (As, Lu (s, l), S_MASS, l)

       omga = dp (As, s, l) / dp (As, As, l)
       
       call lc2 (u, alph, p, omga, s, l)
       
       call equals_float_field (res, lcf (1d0, s, -omga, As, l), S_MASS, l)
       
       if (l2 (res,l)/nrm_f <= tol_elliptic) exit
       
       rho_old = rho
       rho = dp (res0, res, l)
       
       b = (alph/omga) * (rho/rho_old)
       call equals_float_field (p, lcf (1d0, res, b, lcf (1d0, p, -omga, Ap, l), l), S_MASS, l)
       call equals_float_field (Ap, Lu (p, l), S_MASS, l)
       
       iter = iter + 1
    end do
    u%bdry_uptodate = .false.
    call update_bdry (u, l, 88)

    if (present(total_iter)) total_iter = iter
    if (present(nrm_res))    nrm_res = l2 (res,l)/nrm_f 
  end subroutine bicgstab

  function elliptic_fun (u, l)
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

    integer :: d, j

    elliptic_fun_diag = q
    call zero_float_field (elliptic_fun_diag, S_MASS, l)
    
    do d = 1, size(grid)
       scalar => elliptic_fun_diag%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_elliptic_fun_diag, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar)
    end do
  end function elliptic_fun_diag

  subroutine cal_elliptic_fun_diag  (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i
    real(8) :: sigma
  
    id = idx (i, j, offs, dims)
    id_i = id + 1

    sigma = radius/4d0

    if (dom%mask_n%elts(id_i) >= ADJZONE) scalar(id_i) = 2d0 * sqrt (3d0) * dom%areas%elts(id_i)%hex_inv - 4d0/sigma**2
  end subroutine cal_elliptic_fun_diag

  real(8) function relative_error (u, l)
    implicit none
    integer                   :: l
    type(Float_Field), target :: u

    integer                   :: d, j
    real(8)                   :: nrm_err, nrm_sol
    type(Float_Field), target :: err

    err = u
    call zero_float_field (err, S_MASS, l)

    ! Compute absolute error
    do d = 1, size(grid)
       scalar  => err%data(d)%elts
       scalar2 => u%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_err, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (scalar, scalar2)
    end do

    ! Compute relative error
    nrm_sol = l2 (u,   l)
    nrm_err = l2 (err, l)
    relative_error = nrm_err / nrm_sol
  end function relative_error

  subroutine cal_err (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) 
    if (dom%mask_n%elts(id+1) >= ADJZONE) scalar(id+1) = abs (scalar2(id+1) -  exact_sol (dom%node%elts(id+1)))
  end subroutine cal_err

  real(8) function exact_sol (p)
    implicit none
    type (Coord) :: p

    real(8) :: r, sigma

    sigma = radius/4d0

    r = dist (p, sph2cart (0d0, 0d0))
    
    exact_sol = exp (-(r/sigma)**2) 
  end function exact_sol
end module lin_solve_mod
