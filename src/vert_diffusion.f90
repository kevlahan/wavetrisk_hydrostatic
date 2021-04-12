module vert_diffusion_mod
  ! Provides a backwards Euler step for vertical diffusion of buoyancy (temp variable) and
  ! velocity for ocean models. Backwards Euler is unconditionally stable.
  !
  ! Eddy diffusivity and eddy viscosity may be computed either analytically (tke_closure=.false.) or
  ! using a TKE closure model (tke_closure=.true.).
  !
  ! User must supply the following functions in test_case_mod.f90:
  !        source_b, source_t   (bottom and top buoyancy sources)
  !        wind_tau, wind_drag  (magnitude of wind stress tau and wind drag)
  !        r                    (bottom friction)
  use utils_mod
  use test_case_mod
  implicit none
  real(8) :: friction

  ! Parameters for analytic eddy viscosity/diffusion scheme
  real(8), parameter :: Kv_bottom = 2d-3
  real(8), parameter :: Kt_const  = 1d-6

  ! Parameters for TKE closure eddy viscosity/diffusion scheme
  real(8), parameter :: c_e      = 1.0d0
  real(8), parameter :: c_eps    = 1/sqrt(2.0_8)
  real(8), parameter :: C_l      = 2d5
  real(8), parameter :: c_m      = 0.1d0
  real(8), parameter :: C_sfc    = 67.83d0
  real(8), parameter :: e_min    = 1d-6 
  real(8), parameter :: e_0      = 1d-6/sqrt(2.0d0) 
  real(8), parameter :: e_sfc_0  = 1d-4    
  real(8), parameter :: eps_s    = 1d-20   
  real(8), parameter :: kappa_VK = 0.4     ! von Karman constant
  real(8), parameter :: l_0      = 0.04    
  real(8), parameter :: Neps_sq  = 1d-20   
  real(8), parameter :: Ri_c     = 2 / (2 + c_eps/c_m) ! 0.22
  real(8), parameter :: Kv_0     = 1.2d-4 ! NEMO value
  real(8), parameter :: Kt_0     = 1.2d-5 ! NEMO value

  logical, parameter :: tke_closure = .false.
  
  abstract interface
     real(8) function fun1 (eta, ri, z)
       implicit none
       real(8) :: eta, ri, z
     end function fun1
     function fun2 (eta, ri, z)
       use shared_mod
       implicit none
       real(8), dimension(1:EDGE) :: fun2, eta, z
       real(8)                    :: ri
     end function fun2
     real(8) function fun3 (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function fun3
     function fun4 (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
       real(8), dimension(1:EDGE)     :: fun4
     end function fun4
     real(8) function fun5 (p)
       use geom_mod
       implicit none
       type(Coord) :: p
     end function fun5
  end interface
  procedure (fun3), pointer :: bottom_temp_source => null ()
  procedure (fun3), pointer :: top_temp_source    => null ()
  procedure (fun4), pointer :: wind_drag          => null ()
  procedure (fun5), pointer :: tau_mag            => null ()
contains
  subroutine implicit_vertical_diffusion (r, wind_tau, wind_d, source_b, source_t)
    ! Backwards euler step for vertical diffusion
    use adapt_mod
    implicit none
    integer :: l
    real(8) :: r

    interface
       real(8) function source_b (dom, i, j, z_lev, offs, dims)
         use domain_mod
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, z_lev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function source_b
       real(8) function source_t (dom, i, j, z_lev, offs, dims)
         use domain_mod
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, z_lev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
       end function source_t
       real(8) function wind_tau (p)
         use geom_mod
         implicit none
         type(Coord) :: p
       end function wind_tau
       function wind_d (dom, i, j, z_lev, offs, dims)
         use domain_mod
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, z_lev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims
         real(8), dimension(1:EDGE)     :: wind_d
       end function wind_d
    end interface

    bottom_temp_source => source_b
    top_temp_source    => source_t
    wind_drag          => wind_d
    tau_mag            => wind_tau

    friction = r
    
    call update_array_bdry (sol, NONE, 27)
    
    ! Compute eddy diffusivity and eddy viscosity at nodes and layer interfaces at all grid points
    do l = level_end, level_start, -1
       call apply_onescale (turbulence_model, l, z_null, 0, 1)
    end do

    if (tke_closure) then
       tke%bdry_uptodate = .false.
       call WT_after_scalar (tke, wav_tke, level_start-1)
    end if

    Kv%bdry_uptodate = .false.
    Kt%bdry_uptodate = .false.
    call update_vector_bdry (Kv, NONE, 28)
    call update_vector_bdry (Kt, NONE, 29)
    
    ! Apply vertical diffusion to each vertical column using a backwards Euler step
    do l = level_end, level_start, -1
       call apply_onescale (backwards_euler_temp, l, z_null, 0, 1)
       call apply_onescale (backwards_euler_velo, l, z_null, 0, 0)
    end do
    sol%bdry_uptodate = .false.
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine implicit_vertical_diffusion

  subroutine turbulence_model (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for TKE closure equation (or analytic)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, info, k, l
    real(8)                         :: dudz_2, eta, filt, Ri, z
    real(8), dimension(0:zlevels)   :: e, Kv_i, Kt_i, l_eps, l_m, Nsq,  dudzsq
    real(8), dimension(1:zlevels)   :: dz
    real(8), dimension(1:zlevels-1) :: dz_l, diag, rhs, S1, S2
    real(8), dimension(1:zlevels-2) :: diag_l, diag_u
    type(Coord)                     :: p

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1
    p = dom%node%elts(id)
       
    if (tke_closure) then
       do k = 1, zlevels
          dz(k) = dz_i (dom, i, j, k, offs, dims)
       end do

       do l = 1, zlevels-1
          dz_l(l) = interp (dz(l), dz(l+1))
       end do

       call update_Kv_Kt

       ! RHS terms
       do l = 1, zlevels-1
          S1(l) = Kv_i(l) * dudzsq(l) - Kt_i(l) * Nsq(l)
          if (S1(l) < 0.0_8) then ! Patankar "trick"
             S1(l) = Kv_i(l) * dudz_2
             S2(l) = 1 + dt * (c_eps/l_eps(l) * sqrt (e(l)) + Kt_i(l) * Nsq(l) / e(l))
          else
             S2(l) = 1 + dt * c_eps/l_eps(l) * sqrt (e(l))
          end if
       end do

       ! Tridiagonal matrix and rhs entries for linear system
       l = 1
       diag_u(l) = - coeff (dz(l+1), interp(Kv_i(l), Kv_i(l+1))) ! super-diagonal
       diag(l)   = 1 - diag_u(l) 
       rhs(l)    = (e(l) + dt * S1(l)) / S2(l)

       do l = 2, zlevels-2
          diag_u(l)   = - coeff (dz(l+1), interp(Kv_i(l), Kv_i(l+1))) ! super-diagonal
          diag_l(l-1) = - coeff (dz(l),   interp(Kv_i(l), Kv_i(l-1))) ! sub-diagonal
          diag(l)     = 1 - (diag_u(l) + diag_l(l-1))
          rhs(l)      = (e(l) + dt * S1(l)) / S2(l)
       end do

       l = zlevels-1
       diag_l(l-1) = - coeff (dz(l), interp(Kv_i(l), Kv_i(l-1))) ! sub-diagonal
       diag(l)     = 1 - diag_l(l-1)
       rhs(l)      = (e(l) + dt * S1(l)) / S2(l)

       ! Solve tridiagonal linear system
       call dgtsv (zlevels-1, 1, diag_l, diag, diag_u, rhs, zlevels-1, info)
       if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed in TKE computation with info = ", info

       ! Backwards Euler step
       ! e(0) = e_0  at bathymetry
       do l = 1, zlevels-1
          tke(l)%data(d)%elts(id) = max (rhs(l), e_min)
       end do
       filt = 1 - penal_node(zlevels)%data(d)%elts(id)
       tke(zlevels)%data(d)%elts(id) = max (C_sfc * tau_mag (p)*filt / ref_density, e_sfc_0) ! at free surface

       call update_Kv_Kt
    else ! Analytic eddy diffusivity and eddy viscosity
       Kt_i = Kt_analytic ()

       if (mode_split) then
          eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
       else
          eta = free_surface (dom, i, j, z_null, offs, dims)
       end if
       z = dom%topo%elts(id) 

       Kv_i(0) = Kv_analytic (z, eta)
       do l = 1, zlevels
          z = z + dz_i (dom, i, j, l, offs, dims)
          Kv_i(l) = Kv_analytic (z, eta)
       end do
    end if

    ! Assign diffusivity and eddy viscosity 
    do l = 0, zlevels
       Kt(l)%data(d)%elts(id) = Kt_i(l)
       Kv(l)%data(d)%elts(id) = Kv_i(l)
    end do
  contains
    subroutine update_Kv_Kt
      ! TKE and N^2 at interfaces
      implicit none

      e(0) = e_0
      do l = 1, zlevels
         e(l) = tke(l)%data(d)%elts(id) ! impose minimum TKE value
      end do
      
      do l = 0, zlevels
         Nsq(l) = N_sq (dom, i, j, l, offs, dims, dz)
         dudzsq(l) = dudz_sq (dom, i, j, l, offs, dims)
      end do

      ! Length scales
      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_m)

      ! Eddy viscosity and eddy diffusivity at interfaces
      do l = 0, zlevels
         Ri = Richardson (Nsq(l), dudzsq(l))
         Kv_i(l) = Kv_tke (l_m(l), e(l))
         Kt_i(l) = Kt_tke (Kv_i(l), Ri)
      end do
    end subroutine update_Kv_Kt
    
    real(8) function coeff (dz, Kv)
      ! Computes entries of vertical Laplacian matrix
      implicit none
      real(8) :: dz, Kv
      
      coeff = dt * c_e * Kv / (dz_l(l) * dz) / S2(l)
    end function coeff
  end subroutine turbulence_model

  subroutine backwards_euler_temp (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for temp variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, info, k, l
    real(8)                         :: eta, full_mass, theta
    
    real(8), dimension(1:zlevels)   :: diag, dz, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u, dz_l

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1

    ! Layer thicknesses
    do k = 1, zlevels
       dz(k) = dz_i (dom, i, j, k, offs, dims)
    end do
    
    do l = 1, zlevels-1
       dz_l(l) = interp (dz(l), dz(l+1))
    end do

    ! Bottom layer
    k = 1
    diag_u(k) = - coeff(1) ! super-diagonal
    diag(k)   = 1 - diag_u(k)
    rhs(k)    = b() + dt * flux_mean(1)/dz(k) + dt * bottom_temp_source (dom, i, j, z_null, offs, dims)
    
    do k = 2, zlevels-1
       diag_u(k)   = - coeff( 1) ! super-diagonal
       diag_l(k-1) = - coeff(-1) ! sub-diagonal
       diag(k)     = 1  - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b() + dt * (flux_mean(1) - flux_mean(-1))/dz(k)
    end do

    ! Top layer
    k = zlevels
    diag_l(k-1) = - coeff(-1) ! sub-diagonal
    diag(k)     = 1 - diag_u(k-1)
    rhs(k)      = b() - dt * flux_mean(-1)/dz(k) + dt * top_temp_source (dom, i, j, z_null, offs, dims)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed in vertical diffusion of buoyancy, info = ", info

    ! Backwards Euler step
    do k = 1, zlevels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       sol(S_TEMP,k)%data(d)%elts(id) = full_mass * (b_mean() + rhs(k)) - sol_mean(S_TEMP,k)%data(d)%elts(id)
    end do
  contains
    real(8) function coeff (l)
      ! Computes coefficient above (l = 1) or below (l = -1) for vertical Laplacian matrix
      implicit none
      integer :: l
      integer :: kk

      kk = k + min(0,l)
      coeff = dt / (dz_l(kk) * dz(k)) * Kt(kk)%data(d)%elts(id)
    end function coeff
    
    real(8) function b ()
      ! Fluctuating buoyancy
      implicit none
      real(8) :: full_mass, full_theta
      
      full_mass  = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
      full_theta = sol_mean(S_TEMP,k)%data(d)%elts(id) + sol(S_TEMP,k)%data(d)%elts(id)
      b = full_theta / full_mass - b_mean()
    end function b
    
    real(8) function b_mean ()
      ! Mean buoyancy
      implicit none
      b_mean = sol_mean(S_TEMP,k)%data(d)%elts(id) / sol_mean(S_MASS,k)%data(d)%elts(id) 
    end function b_mean
    
    real(8) function flux_mean (l)
      ! Computes flux of mean buoyancy at interface below (l=-1) or above (l=1) vertical level k
      implicit none
      integer :: l

      integer :: kk
      real(8) :: b_0, b_l, mass_0, mass_l, temp_0, temp_l

      kk = k + min(0,l)

      mass_0 = sol_mean(S_MASS,k)%data(d)%elts(id) 
      temp_0 = sol_mean(S_TEMP,k)%data(d)%elts(id) 
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,k+l)%data(d)%elts(id) 
      temp_l = sol_mean(S_TEMP,k+l)%data(d)%elts(id) 
      b_l = temp_l / mass_l

      flux_mean = l * (b_l - b_0) / dz_l(kk) * Kt(kk)%data(d)%elts(id)
    end function flux_mean
  end subroutine backwards_euler_temp

  subroutine backwards_euler_velo (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for velocity variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                                :: d, e, id, info, k, l
    real(8), dimension(1:EDGE,1:zlevels)   :: diag, dz, rhs
    real(8), dimension(1:EDGE,1:zlevels-1) :: diag_l, diag_u, dz_l
    real(8), dimension(1:zlevels)          :: dd, r
    real(8), dimension(1:zlevels-1)        :: dl, du

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    ! Layer thicknesses
    do k = 1, zlevels
       dz(:,k) = dz_e (dom, i, j, k, offs, dims)
    end do
    
    do l = 1, zlevels-1
       dz_l(:,l) = interp_e (dz(:,l), dz(:,l+1))
    end do
    
    ! Bottom layer
    k = 1
    diag_u(:,k) = - coeff(1) ! super-diagonal
    diag(:,k)   = 1 - diag_u(:,k) + dt * friction / dz(:,k)
    rhs(:,k)    = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) 

    do k = 2, zlevels-1
       diag_u(:,k)   = - coeff( 1) ! super-diagonal
       diag_l(:,k-1) = - coeff(-1) ! sub-diagonal
       diag(:,k)     = 1 - (diag_u(:,k) + diag_l(:,k-1))
       rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    end do
    
    ! Top layer
    k = zlevels
    diag_l(:,k-1) = - coeff(-1) ! sub-diagonal
    diag(:,k)     = 1 - diag_l(:,k-1)
    rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) + dt * wind_drag (dom, i, j, z_null, offs, dims)

    ! Solve tridiagonal linear system
    do e = 1, EDGE
       dl = diag_l(e,:); dd = diag(e,:); du = diag_u(e,:); r = rhs(e,:)
       call dgtsv (zlevels, 1, dl, dd, du, r, zlevels, info)
       if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed in vertical diffusion of velocity, info = ", info
       do k = 1, zlevels
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = r(k)
       end do
    end do
  contains
    function coeff (l)
      ! Computes coefficient above (l = 1) or below (l = -1) for vertical Laplacian matrix
      implicit none
      integer                    :: l
      real(8), dimension(1:EDGE) :: coeff

      integer                    :: idE, idN, idNE, kk
      real(8), dimension(1:EDGE) :: Kv_e

      idE  = idx (i+1, j,   offs, dims) + 1
      idNE = idx (i+1, j+1, offs, dims) + 1
      idN  = idx (i,   j+1, offs, dims) + 1

      kk = k + min(0,l)
      Kv_e(RT+1) = interp (Kv(kk)%data(d)%elts(id+1), Kv(kk)%data(d)%elts(idE))
      Kv_e(DG+1) = interp (Kv(kk)%data(d)%elts(id+1), Kv(kk)%data(d)%elts(idNE))
      Kv_e(UP+1) = interp (Kv(kk)%data(d)%elts(id+1), Kv(kk)%data(d)%elts(idN))
            
      coeff = dt / (dz_l(:,kk) * dz(:,k)) * Kv_e
    end function coeff
  end subroutine backwards_euler_velo

  real(8) function N_sq  (dom, i, j, l, offs, dims, dz)
    ! Brunt-Vaisala number N^2 at interface 0 <= l <= zlevels
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:zlevels)  :: dz

    integer :: d, id, ll
    real(8) :: mass0, mass1, temp0, temp1

    ! Top and bottom values (constant density gradients)
    ll = l
    if (l == 0      ) ll = 1
    if (l == zlevels) ll = zlevels - 1
    
    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1
   
    mass0 = sol_mean(S_MASS,ll)%data(d)%elts(id) + sol(S_MASS,ll)%data(d)%elts(id)
    temp0 = sol_mean(S_TEMP,ll)%data(d)%elts(id) + sol(S_TEMP,ll)%data(d)%elts(id)

    mass1 = sol_mean(S_MASS,ll+1)%data(d)%elts(id) + sol(S_MASS,ll+1)%data(d)%elts(id)
    temp1 = sol_mean(S_TEMP,ll+1)%data(d)%elts(id) + sol(S_TEMP,ll+1)%data(d)%elts(id)

    N_sq = grav_accel * (temp1/mass1 - temp0/mass0) / interp (dz(ll), dz(ll+1)) ! -g drho/dz / rho0
  end function N_sq

  real(8) function dudz_sq  (dom, i, j, l, offs, dims)
    ! Twice the kinetic energy at interface  0 <= l <= zlevels, ||du_h/dz||^2
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                      :: d, id, idS, idSW, idW, ll
    real(8), dimension(1:2*EDGE) :: du, dz, prim_dual

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i,   j-1, offs, dims)
    idS  = idx (i-1, j-1, offs, dims)

    ! Product of primal and dual edges
    prim_dual(1:EDGE) = dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    prim_dual(EDGE+1) = dom%len%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
    prim_dual(EDGE+2) = dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
    prim_dual(EDGE+3) = dom%len%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)

     ! Top and bottom values (constant density gradients)
    ll = l
    if (l == 0      ) ll = 1
    if (l == zlevels) ll = zlevels - 1
    
    ! Layer thicknesses centred on interface
    dz(1:EDGE)        = interp_e (dz_e    (dom, i, j, ll, offs, dims), dz_e    (dom, i, j, ll+1, offs, dims)) 
    dz(EDGE+1:2*EDGE) = interp_e (dz_SW_e (dom, i, j, ll, offs, dims), dz_SW_e (dom, i, j, ll+1, offs, dims)) 

    ! Velocity differences
    du(1:EDGE) = sol(S_VELO,ll+1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - sol(S_VELO,ll)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    du(EDGE+1) = sol(S_VELO,ll+1)%data(d)%elts(EDGE*idW+RT+1)  - sol(S_VELO,ll)%data(d)%elts(EDGE*idW+RT+1)
    du(EDGE+2) = sol(S_VELO,ll+1)%data(d)%elts(EDGE*idSW+DG+1) - sol(S_VELO,ll)%data(d)%elts(EDGE*idSW+DG+1)
    du(EDGE+3) = sol(S_VELO,ll+1)%data(d)%elts(EDGE*idS+UP+1)  - sol(S_VELO,ll)%data(d)%elts(EDGE*idS+UP+1)

    ! Energy term (using TRiSK form of kinetic energy)
    dudz_sq = sum ((du/dz)**2 * prim_dual) * dom%areas%elts(id+1)%hex_inv/2
  end function dudz_sq 

  real(8) function Richardson (Nsq, dudzsq)
    ! Richardson number at interface  0 <= l <= zlevels
    implicit none
    real(8) :: Nsq, dudzsq

    Richardson = Nsq / (dudzsq + eps_s)
  end function Richardson

  real(8) function Prandtl (Ri)
    ! Computes Prandtl number given Richardson number
    implicit none
    real(8) :: Ri

    Prandtl = max (0.1d0, Ri_c / max (Ri_c, Ri))
  end function prandtl

  subroutine l_scales (dz, Nsq, tau, tke, l_eps, l_m)
    ! Computes length scales l_eps and l_m at interfaces 0:zlevels for TKE closure for a single vertical column
    implicit none
    real(8),                       intent (in)  :: tau
    real(8), dimension(1:zlevels), intent (in)  :: dz            ! layer thicknesses
    real(8), dimension(0:zlevels), intent (in)  :: Nsq, tke
    real(8), dimension(0:zlevels), intent (out) :: l_eps, l_m

    real(8), dimension(0:zlevels) :: l_dwn, l_up

    integer :: l

    do l = 0, zlevels
       l_up(l) = sqrt (2 * tke(l) / max (Nsq(l), Neps_sq))
    end do
    l_dwn = l_up

    l_dwn(0) = l_0
    do l = 1, zlevels
       l_dwn(l) = min (l_dwn(l-1) + dz(l), l_dwn(l))
    end do

    l_up(zlevels) = kappa_VK * C_l / (ref_density * grav_accel) * tau
    do l = zlevels-1, 0, -1
       l_up(l) = min (l_up(l+1) + dz(l+1), l_up(l)) 
    end do

    ! Returned length scales
    l_eps = sqrt (l_up * l_dwn)
    l_m   = min (l_up, l_dwn)
  end subroutine l_scales
    
  real(8) function Kt_tke (Kv, Ri)
    ! TKE closure eddy diffusivity
    implicit none
    real(8) :: Kv, Ri

    Kt_tke = max (Kv/Prandtl(Ri), Kt_0)
  end function Kt_tke

  real(8) function Kv_tke (l_m, tke)
    ! TKE closure eddy viscosity
    implicit none
    real(8) :: l_m, tke

    Kv_tke = max (c_m * l_m * sqrt(tke), Kv_0)
  end function Kv_tke

  real(8) function Kt_analytic ()
    ! Analytic eddy diffusivity
    implicit none

    Kt_analytic = Kt_const
  end function Kt_analytic

  real(8) function Kv_analytic (z, eta)
    ! Analytic eddy viscosity
    real(8) :: eta, z

    Kv_analytic = Kv_bottom * (1.0_8 + 4 * exp ((z - eta) / abs(max_depth) ))
  end function Kv_analytic

  real(8) function A_vT (Ri)
    implicit none
    real(8) :: Ri

    A_vT = 1d-4 / (1 + 5*Ri)**2 + Kt_0
  end function A_vT

  real(8) function A_vm (Ri)
    implicit none
    real(8) :: Ri

    A_vm = A_vT (Ri) / (1 + 5*Ri) + Kv_0
  end function A_vm
end module vert_diffusion_mod
