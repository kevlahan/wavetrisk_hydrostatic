module vert_diffusion_mod
  ! Provides backwards Euler and forwards Euler time stepping for vertical diffusion of buoyancy (temp variable) and
  ! velocity for ocean models. Backwards Euler is unconditionally stable and is the default choice.
  !
  ! This implementation assumes that sol_mean(S_TEMP,:) is zero (no buoyancy in the mean component)
  ! and flux boundary conditions at bathymetry and free surface.
  !
  ! Eddy diffusivity and eddy viscosity may be computed either analytically (tke_closure=.false.) or
  ! using a TKE closure model (tke_closure=.true.).
  !
  ! User must supply the following functions in test_case_mod.f90:
  !        flux_b, flux_t   (bottom and top buoyancy sources)
  !        wind_tau, wind_f (magnitude of wind stress tau and wind drag)
  !        r                (bottom friction)
  use init_mod
  use utils_mod
  implicit none

  logical, parameter :: implicit     = .true. ! use backwards Euler (T) or forward Euler (F)
  logical, parameter :: surf_wave    = .true. ! use surface wave parameterization (T) or default (F)
  logical, parameter :: enhance_diff = .false. ! enhanced vertical diffusion in unstable regions (T) or rely on TKE closure values (F)
  
  ! Parameters for analytic eddy viscosity/diffusion scheme (e_min, Kt_0, Kv_0 defaults set in shared.f90)
  real(8), parameter :: Kv_bottom = 2d-3
  real(8), parameter :: Kt_const  = 1d-6

  ! Parameters for TKE closure eddy viscosity/diffusion scheme
  real(8), parameter :: c_e      = 1d0               
  real(8), parameter :: c_eps    = 7d-1                ! factor in Ekman depth equation
  real(8), parameter :: C_l      = 2d5                 ! Charnock constant
  real(8), parameter :: c_m      = 1d-1                ! coefficient for eddy viscosity
  real(8), parameter :: C_sfc    = 67.83d0
  real(8), parameter :: e_0      = 1d-6/sqrt(2.0d0)    ! botom boundary condition for TKE
  real(8), parameter :: e_sfc_0  = 1d-4                ! minimum TKE at free surface
  real(8), parameter :: eps_s    = 1d-20               
  real(8), parameter :: kappa_VK = 4d-1                ! von Karman constant
  real(8), parameter :: Kt_enh   = 1d2                 ! enhanced vertical diffusion
  real(8), parameter :: Kv_enh   = 1d2                 ! enhanced vertical diffusion
  real(8), parameter :: l_0      = 4d-2                ! default turbulent length scale
  real(8), parameter :: Neps_sq  = 1d-20               ! background shear
  real(8), parameter :: Nsq_min  = 1d-12               ! threshold for enhanced diffusion
  real(8), parameter :: Ri_c     = 2 / (2 + c_eps/c_m) ! 0.22
  real(8), parameter :: z_0      = 1d-1                ! roughness parameter of free surface
contains
  subroutine vertical_diffusion
    ! Backwards euler step for vertical diffusion
    use adapt_mod
    use time_integr_mod
    implicit none
    integer :: l
    
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
    
    ! Apply vertical diffusion to each vertical column
    if (implicit) then ! backwards Euler step 
       do l = level_end, level_start, -1
          call apply_onescale (backwards_euler_temp, l, z_null, 0, 1)
          call apply_onescale (backwards_euler_velo, l, z_null, 0, 0)
       end do
    else ! forwards Euler step 
       call Euler (sol, wav_coeff, trend_vertical_diffusion, dt)
    end if
    
    sol%bdry_uptodate = .false.
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine vertical_diffusion

  subroutine turbulence_model (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for TKE closure equation (or analytic)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, info, k, l
    real(8)                         :: eta, filt, turb, z
    real(8), dimension(0:zlevels)   :: e, l_eps, l_m, Nsq,  dudzsq
    real(8), dimension(1:zlevels)   :: dz
    real(8), dimension(1:zlevels-1) :: dz_l, diag, rhs, S1, S2
    real(8), dimension(1:zlevels-2) :: diag_l, diag_u
    type(Coord)                     :: p

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1
    p = dom%node%elts(id)

    if (tke_closure) then
       call init_diffuse
       
       ! RHS terms
       do l = 1, zlevels-1
          if (e(l) == 0.0_8) then
             turb = 0.0_8
          else
             turb = c_eps * sqrt (e(l)) / l_eps(l)
          end if
          
          S1(l) = Kv(l)%data(d)%elts(id) * dudzsq(l) - Kt(l)%data(d)%elts(id) * Nsq(l)
          if (S1(l) <= 0.0_8) then ! Patankar "trick"
             S1(l) = Kv(l)%data(d)%elts(id) * dudzsq(l)
             S2(l) = - (turb + Kt(l)%data(d)%elts(id) * Nsq(l) / e(l))
          else
             S2(l) = - turb
          end if
       end do

       ! Tridiagonal matrix and rhs entries for linear system
       l = 1
       diag_u(l) = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id))) ! super-diagonal
       diag(l)   = 1 - diag_u(l) - dt * S2(l)
       rhs(l)    = e(l) + dt * S1(l)

       do l = 2, zlevels-2
          diag_u(l)   = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id))) ! super-diagonal
          diag_l(l-1) = - coeff (dz(l),   interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id))) ! sub-diagonal
          diag(l)     = 1 - (diag_u(l) + diag_l(l-1)) - dt * S2(l)
          rhs(l)      = e(l) + dt * S1(l)
       end do

       l = zlevels-1
       diag_l(l-1) = - coeff (dz(l), interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id)))  ! sub-diagonal
       diag(l)     = 1 - diag_l(l-1) - dt * S2(l)
       rhs(l)      = e(l) + dt * S1(l)

       ! Solve tridiagonal linear system
       call dgtsv (zlevels-1, 1, diag_l, diag, diag_u, rhs, zlevels-1, info)
       if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in TKE computation with info = ", info

       ! Backwards Euler step
       e(0) = e_0
       do l = 1, zlevels-1
          e(l) = rhs(l)
       end do
       e(zlevels) = max (C_sfc * tau_mag (p) / ref_density, e_sfc_0) ! free surface
       
       call update_Kv_Kt
    else ! Analytic eddy diffusivity and eddy viscosity
       if (mode_split) then
          eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
       else
          eta = free_surface (dom, i, j, z_null, offs, dims)
       end if
       z = dom%topo%elts(id)

       Kt(0)%data(d)%elts(id) = Kt_analytic ()
       Kv(0)%data(d)%elts(id) = Kv_analytic (z, eta)
       do l = 1, zlevels
          z = z + dz_i (dom, i, j, l, offs, dims)
          Kt(l)%data(d)%elts(id) = Kt_analytic ()
          Kv(l)%data(d)%elts(id) = Kv_analytic (z, eta)
       end do
    end if
  contains
    subroutine init_diffuse
      ! Initializations
      implicit none
      integer :: k, l
      real(8) :: Ri

      do k = 1, zlevels
         dz(k) = dz_i (dom, i, j, k, offs, dims)
      end do

      do l = 1, zlevels-1
         dz_l(l) = interp (dz(l), dz(l+1))
      end do

      do l = 0, zlevels
         Nsq(l)    = N_sq    (dom, i, j, l, offs, dims, dz)
         dudzsq(l) = dudz_sq (dom, i, j, l, offs, dims)
      end do

      e(0) = e_0
      do l = 1, zlevels
         e(l) = tke(l)%data(d)%elts(id)
      end do

      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_m)

      ! Recompute eddy viscosity and eddy diffusivity after restart
      if (istep == 1) then
         do l = 0, zlevels
            Ri = Richardson (Nsq(l), dudzsq(l))
            Kv(l)%data(d)%elts(id) = Kv_tke (l_m(l), e(l), Nsq(l))
            Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Ri, Nsq(l))
         end do
      end if
    end subroutine init_diffuse

    subroutine update_Kv_Kt
      ! Update eddy diffusivity and eddy viscosity
      implicit none
      real(8) :: Ri
      
      ! Length scales
      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_m)

      ! Eddy viscosity and eddy diffusivity at interfaces
      do l = 0, zlevels
         Ri = Richardson (Nsq(l), dudzsq(l))
         Kv(l)%data(d)%elts(id) = Kv_tke (l_m(l), e(l), Nsq(l))
         Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Ri, Nsq(l))
      end do

      ! Assign tke
      do l = 1, zlevels
         tke(l)%data(d)%elts(id) = e(l)
      end do
    end subroutine update_Kv_Kt
    
    real(8) function coeff (dz, Kv)
      ! Computes entries of vertical Laplacian matrix
      implicit none
      real(8) :: dz, Kv
      
      coeff = dt * c_e * Kv / (dz_l(l) * dz)
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

    real(8), dimension(0:zlevels)   :: z
    real(8), dimension(1:zlevels)   :: diag, dz, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u, dz_l

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1

    if (mode_split) then
       eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
    else
       eta = free_surface (dom, i, j, z_null, offs, dims)
    end if

    ! Layer thicknesses and interface positions
    z(0) = dom%topo%elts(id)
    do k = 1, zlevels
       dz(k) = dz_i (dom, i, j, k, offs, dims)
       z(k) = z(k-1) + dz(k)
    end do
    
    do l = 1, zlevels-1
       dz_l(l) = interp (dz(l), dz(l+1))
    end do

    ! Bottom layer
    k = 1
    diag_u(k) = - coeff (1) ! super-diagonal
    diag(k)   = 1 - diag_u(k)
    rhs(k)    = b() + dt * ( - bottom_buoy_flux (dom, i, j, z_null, offs, dims) + solar_flux ()) / dz(k)

    do k = 2, zlevels-1
       diag_u(k)   = - coeff ( 1) ! super-diagonal
       diag_l(k-1) = - coeff (-1) ! sub-diagonal
       diag(k)     = 1 - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b() + dt * solar_flux () / dz(k)
    end do

    ! Top layer
    k = zlevels
    diag_l(k-1) = - coeff (-1) ! sub-diagonal
    diag(k)     = 1 - diag_u(k-1)
    rhs(k)      = b() + dt * (top_buoy_flux (dom, i, j, z_null, offs, dims) + Q_sr/(ref_density*c_p)) / dz(k)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in vertical diffusion of buoyancy, info = ", info

    ! Backwards Euler step
    do k = 1, zlevels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       sol(S_TEMP,k)%data(d)%elts(id) = full_mass * rhs(k)
    end do
  contains
    real(8) function coeff (l)
      ! Coefficient at interface above (l = 1) or below (l = -1) vertical level k for vertical Laplacian matrix
      implicit none
      integer :: l
      integer :: kk

      kk = k + min (0,l)
      coeff = dt * Kt(kk)%data(d)%elts(id) / (dz_l(kk) * dz(k))
    end function coeff
    
    real(8) function b ()
      ! Fluctuating buoyancy
      implicit none
      real(8) :: full_mass
      
      full_mass = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
      b = sol(S_TEMP,k)%data(d)%elts(id) / full_mass
    end function b

    real(8) function solar_flux ()
      ! Net solar flux in layer 1 <= k < zlevels
      implicit none

      solar_flux = Q_sr * (irradiance (eta-z(k)) - irradiance (eta-z(k-1))) / (ref_density * c_p)
    end function solar_flux
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
    diag_u(:,k) = - coeff (1) ! super-diagonal
    diag(:,k)   = 1 - diag_u(:,k) + dt * bottom_friction / dz(:,k)
    rhs(:,k)    = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) 

    do k = 2, zlevels-1
       diag_u(:,k)   = - coeff ( 1) ! super-diagonal
       diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
       diag(:,k)     = 1 - (diag_u(:,k) + diag_l(:,k-1))
       rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    end do
    
    ! Top layer
    k = zlevels
    diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
    diag(:,k)     = 1 - diag_l(:,k-1)
    rhs(:,k) = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) + dt * wind_flux (dom, i, j, z_null, offs, dims) / dz(:,k)

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

      integer :: kk

      kk = k + min (0, l)

      coeff = dt  * Kv(kk)%data(d)%elts(id+1) / (dz_l(:,kk) * dz(:,k))
    end function coeff
  end subroutine backwards_euler_velo

  real(8) function irradiance (depth)
    ! Downward irradiance
    implicit none
    real(8) :: depth ! depth below free surface
    
    real(8), parameter :: R = 0.58d0
    real(8), parameter :: xi_0 = 0.35 * METRE
    real(8), parameter :: xi_1 =   23 * METRE

    irradiance = Q_sr * (R * exp (-depth/xi_0) + (1 - R) * exp (-depth/xi_1))
  end function irradiance

  real(8) function N_sq  (dom, i, j, l, offs, dims, dz)
    ! Brunt-Vaisala number N^2 = -g drho/dz / rho0 at interface 0 <= l <= zlevels
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:zlevels)  :: dz

    if (l < zlevels .and. l > 0) then
       N_sq = eval (l)
    elseif (l == 0) then
       N_sq = eval (1)
    elseif (l == zlevels) then
       N_sq = eval(zlevels-1)
    end if
  contains
    real(8) function eval (l)
      implicit none
      integer :: l

      integer :: d, id
      real(8) :: dz_l
      real(8) :: b_above, b_below ! buoyancy above and below interface l

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      b_below = sol(S_TEMP,l)%data(d)%elts(id)   / (sol_mean(S_MASS,l)%data(d)%elts(id)   + sol(S_MASS,l)%data(d)%elts(id))
      b_above = sol(S_TEMP,l+1)%data(d)%elts(id) / (sol_mean(S_MASS,l+1)%data(d)%elts(id) + sol(S_MASS,l+1)%data(d)%elts(id))
      dz_l = interp (dz(l), dz(l+1))

      eval = grav_accel * (b_above - b_below) / dz_l ! - g drho/dz / rho0
    end function eval
  end function N_sq

  real(8) function dudz_sq  (dom, i, j, l, offs, dims)
    ! ||du_h/dz||^2 at interfaces 0 <= l <= zlevels
    ! Computed from twice TRiSK form of kinetic energy to use data from only a single column
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    if (l < zlevels .and. l > 0) then
       dudz_sq = eval (l)
    elseif (l == 0) then
       dudz_sq = eval (1)
    elseif (l == zlevels) then
       dudz_sq = eval (zlevels-1)
    end if
  contains
    real(8) function eval (l)
      implicit none
      integer :: l

      integer :: d, id
      real(8), dimension(1:EDGE) :: dudz, prim_dual

      d = dom%id + 1
      id = idx (i, j, offs, dims)

      prim_dual = dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1)

      dudz = (sol(S_VELO,l+1)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - sol(S_VELO,l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)) &
           / interp_e (dz_e (dom, i, j, l, offs, dims), dz_e (dom, i, j, l+1, offs, dims))

      eval = sum (dudz**2 * prim_dual) * dom%areas%elts(id+1)%hex_inv
    end function eval
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

    ! NEMO
    if (Ri < 0.2_8) then
       Prandtl = 1.0_8
    elseif (Ri >= 0.2_8 .and. Ri <= 2.0_8) then
       Prandtl = 5 * Ri
    else
       Prandtl = 10.0_8
    end if
!!$    Prandtl = max (0.1d0, Ri_c / max (Ri_c, Ri)) ! CROCO
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
       l_dwn(l) = sqrt (2 * tke(l) / max (Nsq(l), Neps_sq))
    end do
    l_up = l_dwn

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
    
  real(8) function Kt_tke (Kv, Ri, Nsq)
    ! TKE closure eddy diffusivity
    implicit none
    real(8) :: Kv, Nsq, Ri

    Kt_tke = max (Kv/Prandtl(Ri), Kt_0)

    if (enhance_diff .and. Nsq <= Nsq_min) Kt_tke = Kt_enh ! enhanced vertical diffusion
  end function Kt_tke

  real(8) function Kv_tke (l_m, tke, Nsq)
    ! TKE closure eddy viscosity
    implicit none
    real(8) :: l_m, Nsq, tke

    Kv_tke = max (c_m * l_m * sqrt(tke), Kv_0)
    
    if (enhance_diff .and. Nsq <= Nsq_min) Kv_tke = Kv_enh ! enhanced vertical diffusion
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

  subroutine trend_vertical_diffusion (q, dq)
    ! Trend for eddy diffusivity and eddy viscosity for forward Euler time step
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (q, NONE, 27)

    ! Scalars
    do d = 1, size(grid)
       do k = 1, zlevels
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          mass   =>  q(S_MASS,k)%data(d)%elts
          temp   =>  q(S_TEMP,k)%data(d)%elts
          velo   =>  q(S_VELO,k)%data(d)%elts
          dmass  => dq(S_MASS,k)%data(d)%elts
          dtemp  => dq(S_TEMP,k)%data(d)%elts
          dvelo  => dq(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (trend_scalars_vert_diffuse, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_velo_vert_diffuse,    grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dtemp, dvelo, mass, mean_m, mean_t, temp, velo)
       end do
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_vertical_diffusion

  subroutine trend_scalars_vert_diffuse (dom, i, j, zlev, offs, dims)
    ! Vertical eddy diffusivity of buoyancy
    ! (layer height is not diffused)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i
    real(8) :: dz_k

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1
      
    dmass(id_i) = 0.0_8

    dz_k = dz_i (dom, i, j, zlev, offs, dims)
    
    if (zlev > 1 .and. zlev < zlevels) then
       dtemp(id_i) = scalar_flux(1) - scalar_flux(-1)
    elseif (zlev == 1) then
       dtemp(id_i) = scalar_flux(1) - Kt(1)%data(d)%elts(id_i) * bottom_buoy_flux (dom, i, j, z_null, offs, dims)
    elseif (zlev == zlevels) then
       dtemp(id_i) = Kt(zlevels)%data(d)%elts(id_i) * top_buoy_flux (dom, i, j, z_null, offs, dims) - scalar_flux(-1)
    end if
    dtemp(id_i) = porous_density (dom, i, j, zlev, offs, dims) * dtemp(id_i)
  contains
    real(8) function scalar_flux (l)
      ! Computes flux at interface below (l=-1) or above (l=1) vertical level zlev
      implicit none
      integer :: l

      real(8) :: b_0, b_l, dz_l, mass_0, mass_l, temp_0, temp_l, visc

      visc = Kt(zlev+min(0,l))%data(d)%elts(id_i)

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dz_l = interp (dz_k,  dz_i(dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      scalar_flux = l * visc * (b_l - b_0) / dz_l
    end function scalar_flux
  end subroutine trend_scalars_vert_diffuse

  subroutine trend_velo_vert_diffuse (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: d, id
    real(8), dimension(1:EDGE) :: dz_k

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    
    dz_k = dz_e (dom, i, j, zlev, offs, dims)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (velo_flux(1) - velo_flux(-1)) / dz_k
    elseif  (zlev == 1) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (velo_flux(1) - bottom_friction * velo(EDGE*id+RT+1:EDGE*id+UP+1)) / dz_k
    elseif (zlev == zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (wind_flux (dom, i, j, z_null, offs, dims) - velo_flux(-1)) / dz_k 
    end if
  contains
    function velo_flux (l)
      ! Flux at upper interface (l=1) or lower interface (l=-1)
      implicit none
      integer               :: l
      real(8), dimension(3) :: velo_flux

      real(8), dimension(3) :: dz_l, visc

      visc = Kv(zlev+min(0,l))%data(d)%elts(id+1)
      dz_l = 0.5 * (dz_k + dz_e (dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      velo_flux = l * visc * (sol(S_VELO,zlev+l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - velo(EDGE*id+RT+1:EDGE*id+UP+1)) / dz_l
    end function velo_flux
  end subroutine trend_velo_vert_diffuse
end module vert_diffusion_mod
