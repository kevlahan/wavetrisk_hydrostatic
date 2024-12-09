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
  !        bottom_buoy_flux, top_buoy_flux (bottom and top buoyancy fluxes)
  !        wind_flux                       (magnitude of wind stress tau)
  !        bottom_friction                 (bottom friction)
  !
  ! Set tke_closure = .false. and Kt_const = Kv_bottom = 0 to turn off vertical diffusion, but include wind stress/bottom friction.
  use init_mod
  use utils_mod
  implicit none

  ! Modifiable parameters for 
  real(8)            :: Kt_max       = 1d2     ! maximum eddy diffusion
  logical            :: enhance_diff = .false. ! enhanced vertical diffusion in unstable regions (T) or rely on TKE closure values (F)
  logical            :: patankar     = .false. ! ensure positivity of TKE using "Patankar trick" if shear is weak and stratification is strong (T)
                                               ! or enforce minimum value e_0 of TKE (F)
                                               ! (Patankar trick can produce noisy solutions if velocity is very small)

  ! Parameters for TKE closure eddy viscosity/diffusion scheme
  real(8) :: c_e       = 1.0d0               
  real(8) :: c_eps     = 7.0d-1                ! factor in Ekman depth equation
  real(8) :: C_l       = 2.0d5                 ! Charnock constant
  real(8) :: c_m       = 1.0d-1                ! coefficient for eddy viscosity
  real(8) :: C_sfc     = 6.783d1
  real(8) :: e_0       = 1.0d-6/sqrt(2d0)      ! bottom boundary condition for TKE
  real(8) :: e_min     = 1.0d-6                ! minimum TKE for vertical diffusion
  real(8) :: e_sfc_0   = 1.0d-4                ! minimum TKE at free surface
  real(8) :: eps_s     = 1.0d-20               
  real(8) :: kappa_VK  = 4.0d-1                ! von Karman constant
  real(8) :: Kt_0      = 1.2d-5                ! minimum/initial eddy diffusion
  real(8) :: Kt_const  = 1.0d-6                ! analytic value for eddy diffusion (tke_closure = .false.)
  real(8) :: Kt_enh    = 1.0d2                 ! enhanced eddy diffusion
  real(8) :: Kv_0      = 1.2d-4                ! minimum/initial eddy viscosity
  real(8) :: Kv_bottom = 2.0d-3                ! analytic value for eddy viscosity (tke_closure = .false.)
  real(8) :: Kv_enh    = 1.0d2                 ! enhanced eddy viscosity
  real(8) :: l_0       = 4.0d-2                ! default turbulent length scale
  real(8) :: Neps_sq   = 1.0d-20               ! background shear
  real(8) :: Nsq_min   = 1.0d-12               ! threshold for enhanced diffusion
  real(8) :: Q_sr      = 0.0d0                 ! penetrative part of solar short wave radiation
  real(8) :: rb_0      = 4.0d-4                ! bottom friction
  real(8) :: z_0       = 1.0d-1                ! roughness parameter of free surface
contains
  subroutine vertical_diffusion
    ! Backwards Euler split step for vertical diffusion
    use adapt_mod
    use time_integr_mod
    implicit none
    integer :: l

    ! Compute eddy diffusivity and eddy viscosity at nodes and layer interfaces at all grid points
    do l = level_end, level_start, -1
       call apply_onescale (turbulent_diffusion, l, z_null, 0, 1)
    end do

    if (tke_closure) tke%bdry_uptodate = .false.

    ! Apply vertical diffusion to each vertical column
    do l = level_end, level_start, -1
       call apply_onescale (backwards_euler_temp, l, z_null, 0, 1)
       call apply_onescale (backwards_euler_velo, l, z_null, 0, 0)
    end do
    sol%bdry_uptodate = .false.
  end subroutine vertical_diffusion

  subroutine turbulent_diffusion (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for TKE closure equation (or analytic)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, info, k, l
    real(8)                         :: eta, filt, turb, z
    real(8), dimension(0:zlevels)   :: e, l_eps, l_m, Nsq,  dUdZ2
    real(8), dimension(1:zlevels)   :: dz, Umag
    real(8), dimension(1:zlevels-1) :: dzl, diag, rhs, S1, S2
    real(8), dimension(1:zlevels-2) :: diag_l, diag_u
    type(Coord)                     :: p

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1
    p = dom%node%elts(id)

    if (tke_closure) then
       call init_diffuse

       ! RHS terms
       do l = 1, zlevels-1
          if (e(l) == 0d0) then
             turb = 0d0
          else
             turb = c_eps * sqrt (e(l)) / l_eps(l)
          end if

          S1(l) = Kv(l)%data(d)%elts(id) * dUdZ2(l) - Kt(l)%data(d)%elts(id) * Nsq(l)
          if (patankar .and. S1(l) <= 0d0) then ! Patankar "trick"
             S1(l) = Kv(l)%data(d)%elts(id) * dUdZ2(l)
             S2(l) = - turb - Kt(l)%data(d)%elts(id) * Nsq(l) / e(l)
          else
             S2(l) = - turb 
          end if
       end do

       ! Tridiagonal matrix and rhs entries for linear system
       l = 1
       diag_u(l) = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id)))      ! super-diagonal
       diag(l)   = 1d0 - diag_u(l) - dt * S2(l)
       rhs(l)    = e(l) + dt * S1(l)

       do l = 2, zlevels-2
          diag_u(l)   = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id))) ! super-diagonal
          diag_l(l-1) = - coeff (dz(l),   interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id))) ! sub-diagonal
          diag(l)     = 1d0 - (diag_u(l) + diag_l(l-1)) - dt * S2(l)
          rhs(l)      = e(l) + dt * S1(l)
       end do

       l = zlevels-1
       diag_l(l-1) = - coeff (dz(l), interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id)))      ! sub-diagonal
       diag(l)     = 1d0 - diag_l(l-1) - dt * S2(l)
       rhs(l)      = e(l) + dt * S1(l)

       ! Solve tridiagonal linear system
       call dgtsv (zlevels-1, 1, diag_l, diag, diag_u, rhs, zlevels-1, info)
       if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in TKE computation with info = ", info

       ! Backwards Euler step
       e(0) = e_0
       do l = 1, zlevels-1
          e(l) = rhs(l)
          if (.not. patankar) e(l) = min (rhs(l), e_min) ! ensure TKE is non-negative
       end do
       e(zlevels) = max (C_sfc * tau_mag (p) / ref_density, e_sfc_0) ! free surface

       call update_Kv_Kt
    else ! Analytic eddy diffusivity and eddy viscosity
       if (mode_split) then
          eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
       else
          eta = free_surface (dom, i, j, z_null, offs, dims, sol)
       end if
       z = topography%data(d)%elts(id)

       Kt(0)%data(d)%elts(id) = Kt_analytic ()
       Kv(0)%data(d)%elts(id) = Kv_analytic (z, eta)
       do l = 1, zlevels
          z = z + dz_i (dom, i, j, l, offs, dims, sol)
          Kt(l)%data(d)%elts(id) = Kt_analytic ()
          Kv(l)%data(d)%elts(id) = Kv_analytic (z, eta)
       end do
    end if
  contains
    subroutine init_diffuse
      ! Initializations
      use io_mod, only : kinetic_energy
      implicit none
      integer :: k, l
      real(8) :: Ri

      do k = 1, zlevels
         dz(k) = dz_i (dom, i, j, k, offs, dims, sol)
         Umag(k) = sqrt (2d0 * kinetic_energy (dom, i, j, k, offs, dims))
      end do

      do l = 1, zlevels-1
         dzl(l) = interp (dz(l), dz(l+1))
      end do

      do l = 0, zlevels
         dUdZ2(l) = dUdZ_sq (l)
      end do

      do l = 0, zlevels
         Nsq(l) = N_sq (dom, i, j, l, offs, dims, dz)
      end do

      e(0) = e_0
      do l = 1, zlevels
         e(l) = max (tke(l)%data(d)%elts(id), e_min)
      end do

      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_m)

      ! Recompute eddy viscosity and eddy diffusivity after restart
      if (istep == 1) then
         do l = 0, zlevels
            Ri = Richardson (Nsq(l), dUdZ2(l))
            Kv(l)%data(d)%elts(id) = Kv_tke (e(l), l_m(l), Nsq(l))
            Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Nsq(l), Ri)
         end do
      end if
    end subroutine init_diffuse

    real(8) function dUdZ_sq (l)
      ! ||du_h/dz||^2 at interfaces 0 <= l <= zlevels
      ! (computed from twice TRiSK form of kinetic energy to use data from only a single colum)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, l
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      real(8) :: dU, dZ

      if (l == 0 .or. l == zlevels) then
         dUdZ_sq = 0d0
      else
         dUdZ_sq = ( (Umag(l+1) - Umag(l)) / dzl(l) )**2
      end if
    end function dUdZ_sq

    subroutine update_Kv_Kt
      ! Update eddy diffusivity and eddy viscosity
      implicit none
      real(8) :: Ri

      ! Length scales
      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_m)

      ! Eddy viscosity and eddy diffusivity at interfaces
      do l = 0, zlevels
         Ri = Richardson (Nsq(l), dUdZ2(l))
         Kv(l)%data(d)%elts(id) = Kv_tke (e(l), l_m(l), Nsq(l))
         Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Nsq(l), Ri)
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

      coeff = dt * c_e * Kv / (dzl(l) * dz)
    end function coeff
  end subroutine turbulent_diffusion

  subroutine backwards_euler_temp (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for temp variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id, info, k, l
    real(8)                         :: eta, rho_dz, theta

    real(8), dimension(0:zlevels)   :: z
    real(8), dimension(1:zlevels)   :: diag, dz, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u, dzl

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1

    if (mode_split) then
       eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
    else
       eta = free_surface (dom, i, j, z_null, offs, dims, sol)
    end if

    ! Layer thicknesses and interface positions
    z(0) = topography%data(d)%elts(id)
    do k = 1, zlevels
       dz(k) = dz_i (dom, i, j, k, offs, dims, sol)
       z(k) = z(k-1) + dz(k)
    end do

    do l = 1, zlevels-1
       dzl(l) = interp (dz(l), dz(l+1))
    end do

    ! Bottom layer
    k = 1
    diag_u(k) = - coeff (1) ! super-diagonal
    diag(k)   = 1d0 - diag_u(k)
    rhs(k)    = b() + dt * ( - bottom_buoy_flux (dom, i, j, z_null, offs, dims) + solar_flux ()) / dz(k)

    do k = 2, zlevels-1
       diag_u(k)   = - coeff ( 1) ! super-diagonal
       diag_l(k-1) = - coeff (-1) ! sub-diagonal
       diag(k)     = 1d0 - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b () + dt * solar_flux () / dz(k)
    end do

    ! Top layer
    k = zlevels
    diag_l(k-1) = - coeff (-1) ! sub-diagonal
    diag(k)     = 1d0 - diag_u(k-1)
    rhs(k)      = b() + dt * (top_buoy_flux (dom, i, j, z_null, offs, dims) + Q_sr/(ref_density*c_p)) / dz(k)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in vertical diffusion of buoyancy, info = ", info

    ! Backwards Euler step
    do k = 1, zlevels
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       sol(S_TEMP,k)%data(d)%elts(id) = rho_dz * rhs(k)
    end do
  contains
    real(8) function coeff (l)
      ! Coefficient at interface above (l = 1) or below (l = -1) vertical level k for vertical Laplacian matrix
      implicit none
      integer :: l
      integer :: kk

      kk = k + min (0,l)
      coeff = dt * Kt(kk)%data(d)%elts(id) / (dzl(kk) * dz(k))
    end function coeff

    real(8) function b ()
      ! Fluctuating buoyancy
      implicit none
      real(8) :: rho_dz

      rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
      b = sol(S_TEMP,k)%data(d)%elts(id) / rho_dz
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
    real(8), dimension(1:EDGE,1:zlevels-1) :: diag_l, diag_u, dzl
    real(8), dimension(1:zlevels)          :: dd, r
    real(8), dimension(1:zlevels-1)        :: dl, du

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    ! Layer thicknesses
    do k = 1, zlevels
       dz(:,k) = dz_e (dom, i, j, k, offs, dims, sol)
    end do

    do l = 1, zlevels-1
       dzl(:,l) = interp_e (dz(:,l), dz(:,l+1))
    end do

    ! Bottom layer
    k = 1
    diag_u(:,k) = - coeff (1) ! super-diagonal
    diag(:,k)   = 1d0 - diag_u(:,k) + dt * bottom_friction / dz(:,k)
    rhs(:,k)    = sol(S_VELO,k)%data(d)%elts(id_edge(id)) 

    do k = 2, zlevels-1
       diag_u(:,k)   = - coeff ( 1) ! super-diagonal
       diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
       diag(:,k)     = 1d0 - (diag_u(:,k) + diag_l(:,k-1))
       rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(id_edge(id))
    end do

    ! Top layer
    k = zlevels
    diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
    diag(:,k)     = 1d0 - diag_l(:,k-1)
    rhs(:,k) = sol(S_VELO,k)%data(d)%elts(id_edge(id)) + dt * wind_flux (dom, i, j, z_null, offs, dims) / dz(:,k)

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

      coeff = dt  * Kv(kk)%data(d)%elts(id+1) / (dzl(:,kk) * dz(:,k))
    end function coeff
  end subroutine backwards_euler_velo

  real(8) function irradiance (depth)
    ! Downward irradiance
    implicit none
    real(8) :: depth ! depth below free surface

    real(8), parameter :: R = 0.58d0
    real(8), parameter :: xi_0 = 0.35d0 * METRE
    real(8), parameter :: xi_1 =   23d0 * METRE

    irradiance = Q_sr * (R * exp (-depth/xi_0) + (1d0 - R) * exp (-depth/xi_1))
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
      real(8) :: dzl
      real(8) :: b_above, b_below ! buoyancy above and below interface l

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      b_below = sol(S_TEMP,l)%data(d)%elts(id)   / (sol_mean(S_MASS,l)%data(d)%elts(id)   + sol(S_MASS,l)%data(d)%elts(id))
      b_above = sol(S_TEMP,l+1)%data(d)%elts(id) / (sol_mean(S_MASS,l+1)%data(d)%elts(id) + sol(S_MASS,l+1)%data(d)%elts(id))
      dzl = interp (dz(l), dz(l+1))

      eval = grav_accel * (b_above - b_below) / dzl ! - g drho/dz / rho0
    end function eval
  end function N_sq

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

    real(8) :: Ri_c

    Ri_c = 2d0 / (2d0 + c_eps/c_m) 

    ! NEMO
    if (Ri < 0.2d0) then
       Prandtl = 1d0
    elseif (Ri >= 0.2d0 .and. Ri <= 2d0) then
       Prandtl = 5d0 * Ri
    else
       Prandtl = 10d0
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
       l_dwn(l) = sqrt (2d0 * tke(l) / max (Nsq(l), Neps_sq))
    end do
    l_up = l_dwn

    l_dwn(0) = l_0
    do l = 1, zlevels
       l_dwn(l) = min (l_dwn(l-1) + dz(l), l_dwn(l))
    end do

    l_up(zlevels) = max (l_0, kappa_VK * C_l / (ref_density * grav_accel) * tau)
    do l = zlevels-1, 0, -1
       l_up(l) = min (l_up(l+1) + dz(l+1), l_up(l)) 
    end do

    ! Returned length scales
    l_eps = sqrt (l_up * l_dwn)
    l_m   = min (l_up, l_dwn)
  end subroutine l_scales

  real(8) function Kt_tke (Kv, Nsq, Ri)
    ! TKE closure eddy diffusivity
    implicit none
    real(8) :: Kv, Nsq, Ri

    Kt_tke = min (Kt_max, max (Kv/Prandtl(Ri), Kt_0))

    if (enhance_diff .and. Nsq <= Nsq_min) Kt_tke = Kt_enh ! enhanced vertical diffusion
  end function Kt_tke

  real(8) function Kv_tke (e, l_m, Nsq)
    ! TKE closure eddy viscosity
    implicit none
    real(8) :: e, l_m, Nsq

    Kv_tke = max (c_m * l_m * sqrt(e), Kv_0)

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

    Kv_analytic = Kv_bottom * (1d0 + 4d0 * exp ((z - eta) / abs(max_depth) ))
  end function Kv_analytic

  subroutine trend_vertical_diffusion (q, dq)
    ! Trend for eddy diffusivity and eddy viscosity for forward Euler time step
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q, dq

    integer :: d, k, p

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

    dmass(id_i) = 0d0

    dz_k = dz_i (dom, i, j, zlev, offs, dims, sol)

    if (zlev > 1 .and. zlev < zlevels) then
       dtemp(id_i) = scalar_flux(1) - scalar_flux(-1)
    elseif (zlev == 1) then
       dtemp(id_i) = scalar_flux(1) - Kt(1)%data(d)%elts(id_i) * bottom_buoy_flux (dom, i, j, z_null, offs, dims)
    elseif (zlev == zlevels) then
       dtemp(id_i) = Kt(zlevels)%data(d)%elts(id_i) * top_buoy_flux (dom, i, j, z_null, offs, dims) - scalar_flux(-1)
    end if
    dtemp(id_i) = porous_density (d, id_i, zlev) * dtemp(id_i)
  contains
    real(8) function scalar_flux (l)
      ! Computes flux at interface below (l=-1) or above (l=1) vertical level zlev
      implicit none
      integer :: l
      integer :: dzl

      real(8) :: b_0, b_l, mass_0, mass_l, temp_0, temp_l, visc

      visc = Kt(zlev+min(0,l))%data(d)%elts(id_i)

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dzl = interp (dz_k,  dz_i(dom, i, j, zlev+l, offs, dims, sol)) ! thickness of layer centred on interface

      scalar_flux = l * visc * (b_l - b_0) / dzl
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

    dz_k = dz_e (dom, i, j, zlev, offs, dims, sol)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(id_edge(id)) = (velo_flux(1) - velo_flux(-1)) / dz_k
    elseif  (zlev == 1) then
       dvelo(id_edge(id)) = (velo_flux(1) - bottom_friction * velo(id_edge(id))) / dz_k
    elseif (zlev == zlevels) then
       dvelo(id_edge(id)) = (wind_flux (dom, i, j, z_null, offs, dims) - velo_flux(-1)) / dz_k 
    end if
  contains
    function velo_flux (l)
      ! Flux at upper interface (l=1) or lower interface (l=-1)
      implicit none
      integer               :: l
      real(8), dimension(3) :: velo_flux

      real(8), dimension(3) :: dzl, visc

      visc = Kv(zlev+min(0,l))%data(d)%elts(id+1)
      dzl = 0.5d0 * (dz_k + dz_e (dom, i, j, zlev+l, offs, dims, sol)) ! thickness of layer centred on interface

      velo_flux = l * visc * (sol(S_VELO,zlev+l)%data(d)%elts(id_edge(id)) - velo(id_edge(id))) / dzl
    end function velo_flux
  end subroutine trend_velo_vert_diffuse
end module vert_diffusion_mod
