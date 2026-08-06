module vert_diffusion_mod
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !   Provides backwards Euler and forwards Euler time stepping for vertical diffusion of buoyancy (temp variable) and
  !   velocity for ocean models. Backwards Euler is unconditionally stable and is the default choice.  Uses either a
  !   one-equation model turbulent kinetic energy closure scheme to compute eddy viscosity Kv and eddy diffusivity Kt if
  !   tke_closure = .true. (as in NEMO 10.1.3), or an analytic profile if tke_closure = .false.
  !
  !   This implementation assumes that sol_mean(S_TEMP,:) is zero (no buoyancy in the mean component)
  !   and flux boundary conditions at bathymetry and free surface.
  !
  !   Eddy diffusivity and eddy viscosity may be computed either analytically (tke_closure=.false.) or
  !   using a TKE closure model (tke_closure=.true.).
  !
  !   User must supply the following functions in test_case_mod.f90:
  !          bottom_buoy_flux, top_buoy_flux (bottom and top buoyancy fluxes)
  !          wind_flux                       (magnitude of wind stress tau)
  !          bottom_friction                 (bottom friction)
  !
  !   Set tke_closure = .false. and Kt_const = Kv_bottom = 0 to turn off vertical diffusion
  !   but still include wind stress/bottom friction.
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use kind_mod,   only : dp
  use shared_mod, only : Coord, EDGE, METRE, N_BDRY, N_VARIABLE, a_0, c_p, dt, eps, grav_accel, istep, level_start, level_end, &
       mode_split,ref_density, tke_closure, S_MASS, S_TEMP, S_VELO,  zlevels, zmax, z_null      

  use arch_mod,       only : rank
  use domain_ops_mod, only : apply_onescale, apply_onescale_to_patch
  use init_mod,       only : bottom_buoy_flux, bottom_friction, tau_mag, top_buoy_flux, wind_flux
  use utils_mod,      only : dz_e, dz_i, free_surface, interp, interp_e, phi_edge, porous_density

   use domain_mod, only : Domain, Float_Field, grid, Kv, Kt, dtemp, dmass, dvelo, mass, temp, velo, mean_m, mean_t, &
       sol, sol_mean, tke, topography, id_edge, idx
  
  implicit none

  private
  public :: trend_vertical_diffusion, vertical_diffusion
  public :: enhance_diff, patankar, C_e, C_eps, C_l, C_k, C_srf, e_0, e_min, e_min_srf, eps_s, kappa_VK
  public :: Kt_enh, Kt_max, Kt_min, Kt_mol, Kv_max, Kv_min, Kv_mol, rb_0, l_0, l_min, Neps_sq, Nsq_min, z_0, z_mixed
  public :: Q_sr, R_lw, xi_lw, xi_sw, alpha_k

  ! Parameters for TKE closure 
  logical :: enhance_diff = .false.              ! enhanced vertical diffusion in unstable regions with very small Nsq < Nsq_min 
  logical :: patankar     = .false.              ! ensure positivity of TKE using "Patankar trick" if shear is weak and stratification is strong (T)
                                                 ! or enforce minimum value e_0 of TKE (F) (can produce noisy solutions if velocity is very small)
  real(dp) :: C_e         = 1.0_dp               
  real(dp) :: C_eps       = 7.0e-1_dp            ! factor in Ekman depth equation
  real(dp) :: C_l         = 2.0e5_dp             ! Charnock constant
  real(dp) :: C_k         = 1.0e-1_dp            ! coefficient for eddy viscosity
  real(dp) :: C_srf       = 6.783e1_dp           ! coefficient of surface value for TKE 

  real(dp) :: e_0         = 1.0e-6/sqrt(2.0_dp)  ! bottom boundary condition for TKE: e_min/sqrt(2)
  real(dp) :: e_min       = 1.0e-6_dp            ! minimum TKE 
  real(dp) :: e_min_srf   = 1.0e-4_dp            ! minimum TKE at free surface

  real(dp) :: eps_s       = 1.0e-20_dp           ! background shear
  real(dp) :: kappa_VK    = 4.0e-1_dp            ! von Karman constant

  real(dp) :: Kt_enh      = 1.0_dp               ! enhanced eddy diffusion for Nsq < Nsq_min
  real(dp) :: Kt_max      = 1.0e-2_dp            ! maximum eddy diffusion
  real(dp) :: Kt_min      = 1.2e-5_dp            ! minimum/initial eddy diffusion 
  real(dp) :: Kt_mol      = 1.0e-7_dp            ! molecular diffusivity of seawater (not used)

  real(dp) :: Kv_max      = 1.0e-2_dp            ! maximum eddy viscosity
  real(dp) :: Kv_min      = 1.2e-4_dp            ! minimum/initial eddy viscosity 
  real(dp) :: Kv_mol      = 1.0e-6_dp            ! molecular viscosity of seawater (not used)

  real(dp) :: z_mixed     =   -200 * METRE       ! lower boundary of mixed layer (used with tke_closure = .false.)
  
  real(dp) :: l_0         = 4.0e-2 * METRE       ! surface buoyancy minimum length scale
  real(dp) :: l_min       = 1.0e-2 * METRE       ! minimum mixing length: Kv_mol/(C_k sqrt(e_min)) 

  real(dp) :: Neps_sq     = 1.0e-20_dp           ! background shear
  real(dp) :: Nsq_min     = 1.0e-12_dp           ! threshold for enhanced diffusion

  real(dp) :: rb_0        = 4.0e-4_dp            ! bottom friction
  real(dp) :: z_0         = 1.0e-1_dp            ! roughness parameter of free surface

  ! NEMO 2-band solar forcing (5.4) modified for buoyancy
  real(dp) :: Q_sr        = 0.0_dp               ! incoming solar heat flux (set to 171.25 W/m^2 average value to include solar forcing)
  real(dp) :: R_lw        = 0.58_dp              ! fraction of Qsr that in almost non-penetrative (longwave) wavebands
  real(dp) :: xi_lw       = 0.35 * METRE         ! longwave penetration depth 
  real(dp) :: xi_sw       =   23 * METRE         ! shortwave penetration depth (400 nm to 700 nm)
  real(dp) :: alpha_k                            ! thermal expansion coefficient

  
  interface
     
   subroutine dgtsv(n, nrhs, dl, d, du, b, ldb, info)
      use kind_mod, only : dp
      implicit none

      integer, intent(in)    :: n
      integer, intent(in)    :: nrhs
      integer, intent(in)    :: ldb

      real(dp), intent(inout) :: dl(*)
      real(dp), intent(inout) :: d(*)
      real(dp), intent(inout) :: du(*)
      real(dp), intent(inout) :: b(ldb,*)

      integer, intent(out) :: info
    end subroutine dgtsv
    
end interface


contains

  
  subroutine vertical_diffusion
    ! Backwards Euler split step for vertical diffusion
    use adapt_mod
    use time_integr_mod
    
    implicit none
    integer :: l

    alpha_k = a_0 / ref_density

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
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, z_null
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                          :: d, id, info, l
    real(dp)                         :: eta, turb, z
    real(dp), dimension(0:zlevels)   :: e, l_eps, l_k, Nsq,  dUdZ2
    real(dp), dimension(1:zlevels)   :: dz, Umag
    real(dp), dimension(1:zlevels-1) :: dzl, diag, rhs, S1, S2
    real(dp), dimension(1:zlevels-2) :: diag_l, diag_u
    type(Coord)                      :: p

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1
    p = dom%node%elts(id)

    if (tke_closure) then ! TKE model for eddy diffusivity and eddy viscosity
       call init_diffuse

       ! RHS terms
       do l = 1, zlevels-1
          if (e(l) < eps (1.0_dp)) then
             turb = 0.0_dp
          else
             turb = C_eps * sqrt (e(l)) / l_eps(l)
          end if

          S1(l) = Kv(l)%data(d)%elts(id) * dUdZ2(l) - Kt(l)%data(d)%elts(id) * Nsq(l)
          if (patankar .and. S1(l) <= 0.0_dp) then ! Patankar "trick"
             S1(l) = Kv(l)%data(d)%elts(id) * dUdZ2(l)
             S2(l) = - turb - Kt(l)%data(d)%elts(id) * Nsq(l) / e(l)
          else
             S2(l) = - turb 
          end if
       end do

       ! Tridiagonal matrix and rhs entries for linear system
       l = 1
       diag_u(l) = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id)))      ! super-diagonal
       diag(l)   = 1.0_dp - diag_u(l) - dt * S2(l)
       rhs(l)    = e(l) + dt * S1(l)

       do l = 2, zlevels-2
          diag_u(l)   = - coeff (dz(l+1), interp (Kv(l)%data(d)%elts(id), Kv(l+1)%data(d)%elts(id))) ! super-diagonal
          diag_l(l-1) = - coeff (dz(l),   interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id))) ! sub-diagonal
          diag(l)     = 1.0_dp - (diag_u(l) + diag_l(l-1)) - dt * S2(l)
          rhs(l)      = e(l) + dt * S1(l)
       end do

       l = zlevels-1
       diag_l(l-1) = - coeff (dz(l), interp (Kv(l)%data(d)%elts(id), Kv(l-1)%data(d)%elts(id)))      ! sub-diagonal
       diag(l)     = 1.0_dp - diag_l(l-1) - dt * S2(l)
       rhs(l)      = e(l) + dt * S1(l)

       ! Solve tridiagonal linear system
       call dgtsv (zlevels-1, 1, diag_l, diag, diag_u, rhs, zlevels-1, info)
       if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in TKE computation with info = ", info

       ! Backwards Euler step
       e(0) = e_0
       do l = 1, zlevels-1
          e(l) = rhs(l)
          if (.not. patankar) e(l) = max (rhs(l), e_min)               ! ensure TKE is non-negative
       end do
       e(zlevels) = max (C_srf * tau_mag (p) / ref_density, e_min_srf) ! free surface boundary conduition

       call update_Kv_Kt
    else ! analytic eddy diffusivity and eddy viscosity
       if (mode_split) then
          eta = sol(S_MASS,zlevels+1)%data(d)%elts(id)
       else
          eta = free_surface (dom, i, j, z_null, offs, dims, sol)
       end if
       z = topography%data(d)%elts(id)

       Kt(0)%data(d)%elts(id) = Kt_analytic (z, eta)
       Kv(0)%data(d)%elts(id) = Kv_analytic (z, eta)
       do l = 1, zlevels
          z = z + dz_i (dom, i, j, l, offs, dims, sol)
          Kt(l)%data(d)%elts(id) = Kt_analytic (z, eta)
          Kv(l)%data(d)%elts(id) = Kv_analytic (z, eta)
       end do
    end if
    
  contains
    
    subroutine init_diffuse
      ! Initializations
      use io_mod, only : kinetic_energy
      implicit none
      integer  :: k, l
      real(dp) :: Ri

      do k = 1, zlevels
         dz(k) = dz_i (dom, i, j, k, offs, dims, sol)
         Umag(k) = sqrt (2 * kinetic_energy (dom, i, j, k, offs, dims))
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

      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_k)

      ! Recompute eddy viscosity and eddy diffusivity after restart
      if (istep == 1) then
         do l = 0, zlevels
            Ri = Richardson (Nsq(l), dUdZ2(l))
            Kv(l)%data(d)%elts(id) = Kv_tke (e(l), l_k(l), Nsq(l))
            Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Nsq(l), Ri)
         end do
      end if
      
    end subroutine init_diffuse
    

    function dUdZ_sq (l) result(val)
      ! ||du_h/dz||^2 at interfaces 0 <= l <= zlevels
      ! (computed from twice TRiSK form of kinetic energy to use data from only a single colum)
      
      implicit none
      
      integer, intent(in) :: l
      real(dp)            :: val

      if (l == 0 .or. l == zlevels) then
         val = 0.0_dp
      else
         val = ( (Umag(l+1) - Umag(l)) / dzl(l) )**2
      end if
    end function dUdZ_sq

    subroutine update_Kv_Kt
      ! Update eddy diffusivity and eddy viscosity
      
      implicit none
      
      real(dp) :: Ri ! Richardson number

      ! Length scales
      call l_scales (dz, Nsq, tau_mag (p), e, l_eps, l_k)

      ! Eddy viscosity and eddy diffusivity at interfaces
      do l = 0, zlevels
         Ri = Richardson (Nsq(l), dUdZ2(l))
         Kv(l)%data(d)%elts(id) = Kv_tke (e(l), l_k(l), Nsq(l))
         Kt(l)%data(d)%elts(id) = Kt_tke (Kv(l)%data(d)%elts(id), Nsq(l), Ri)
      end do

      ! Assign tke
      do l = 1, zlevels
         tke(l)%data(d)%elts(id) = max (e(l), e_min)
      end do
    end subroutine update_Kv_Kt
    

    function coeff (dz, Kv) result(val)
      ! Computes entries of vertical Laplacian matrix
      
      implicit none
      
      real(dp), intent(in) :: dz  ! layer depth
      real(dp), intent(in) :: Kv  ! eddy diffusivity
      real(dp)             :: val

      val = dt * C_e * Kv / (dzl(l) * dz)
    end function coeff
  end subroutine turbulent_diffusion
  

  subroutine backwards_euler_temp (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for buoyancy
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, z_null
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                          :: d, id, info, k, l
    real(dp)                         :: eta, rho_dz

    real(dp), dimension(0:zlevels)   :: z
    real(dp), dimension(1:zlevels)   :: diag, dz, rhs
    real(dp), dimension(1:zlevels-1) :: diag_l, diag_u, dzl

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
    diag(k)   = 1.0_dp - diag_u(k)
    rhs(k)    = b() + dt * ( - bottom_buoy_flux (dom, i, j, z_null, offs, dims) + solar_flux ()) / dz(k)

    do k = 2, zlevels-1
       diag_u(k)   = - coeff ( 1) ! super-diagonal
       diag_l(k-1) = - coeff (-1) ! sub-diagonal
       diag(k)     = 1.0_dp - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b () + dt * solar_flux () / dz(k)
    end do

    ! Top layer
    k = zlevels
    diag_l(k-1) = - coeff (-1) ! sub-diagonal
    diag(k)     = 1.0_dp - diag_u(k-1)
    rhs(k)      = b() + dt * (top_buoy_flux (dom, i, j, z_null, offs, dims) + Q_sr * alpha_k /(ref_density*c_p)) / dz(k)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed in vertical diffusion of buoyancy, info = ", info

    ! Backwards Euler step
    do k = 1, zlevels
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       sol(S_TEMP,k)%data(d)%elts(id) = rho_dz * rhs(k)
    end do
    
  contains
    
    function coeff (l) result(val)
      ! Coefficient at interface above (l = 1) or below (l = -1) vertical level k for vertical Laplacian matrix
      implicit none
      integer, intent(in) :: l
      real(dp)            :: val
      
      integer :: kk

      kk = k + min (0,l)
      
      val = dt * Kt(kk)%data(d)%elts(id) / (dzl(kk) * dz(k))
    end function coeff

    function b () result(val)
      ! Fluctuating buoyancy
      
      implicit none
      
      real(dp) :: rho_dz
      real(dp) :: val

      rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
      
      val = sol(S_TEMP,k)%data(d)%elts(id) / rho_dz
    end function b

    function solar_flux () result(val)
      ! Net solar flux in layer 1 <= k < zlevels from NEMO modified for buoyancy
      
      implicit none
      
      real(dp) :: val
      
      val = (irradiance (eta-z(k)) - irradiance (eta-z(k-1))) * alpha_k / (ref_density * c_p)
    end function solar_flux
  end subroutine backwards_euler_temp
  

  function irradiance (depth) result(val)
    ! Downward irradiance
    
    implicit none
    
    real(dp), intent(in) :: depth ! depth below free surface
    real(dp)             :: val

    val = Q_sr * (R_lw * exp (-depth/xi_lw) + (1.0_dp - R_lw) * exp (-depth/xi_sw))
  end function irradiance
  

  subroutine backwards_euler_velo (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for velocity variable
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, z_null
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                                 :: d, e, id, info, k, l
    real(dp), dimension(1:EDGE,1:zlevels)   :: diag, dz, rhs
    real(dp), dimension(1:EDGE,1:zlevels-1) :: diag_l, diag_u, dzl
    real(dp), dimension(1:zlevels)          :: dd, r
    real(dp), dimension(1:zlevels-1)        :: dl, du

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
    diag(:,k)   = 1.0_dp - diag_u(:,k) + dt * phi_edge (d, id, 1) * bottom_friction / dz(:,k)
    rhs(:,k)    = sol(S_VELO,k)%data(d)%elts(id_edge(id)) 

    do k = 2, zlevels-1
       diag_u(:,k)   = - coeff ( 1) ! super-diagonal
       diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
       diag(:,k)     = 1.0_dp - (diag_u(:,k) + diag_l(:,k-1))
       rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(id_edge(id))
    end do

    ! Top layer
    k = zlevels
    diag_l(:,k-1) = - coeff (-1) ! sub-diagonal
    diag(:,k)     = 1.0_dp - diag_l(:,k-1)
    rhs(:,k) = sol(S_VELO,k)%data(d)%elts(id_edge(id)) + dt * wind_flux (dom, i, j, zlevels, offs, dims) / dz(:,k)

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
    
    function coeff (l) result(val)
      ! Computes coefficient above (l = 1) or below (l = -1) for vertical Laplacian matrix
      
      implicit none
      
      integer, intent(in) :: l
      real(dp)            :: val(EDGE)

      integer :: kk

      kk = k + min (0, l)

      val = dt * Kv(kk)%data(d)%elts(id+1) / (dzl(:,kk) * dz(:,k))
    end function coeff
    
  end subroutine backwards_euler_velo
  

  function N_sq  (dom, i, j, l, offs, dims, dz) result(val)
    ! Brunt-Vaisala number N^2 = -g drho/dz / rho0 at interface 0 <= l <= zlevels
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, l
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp),     intent(in)    :: dz(1:zlevels)
    real(dp)                    :: val

    val = 0.0_dp
    if (l < zlevels .and. l > 0) then
       val = eval (l)
    elseif (l == 0) then
       val = eval (1)
    elseif (l == zlevels) then
       val = eval(zlevels-1)
    end if
  contains
    real(dp) function eval (l)
      implicit none
      integer, intent(in) :: l

      integer  :: d, id
      real(dp) :: dzl
      real(dp) :: b_above, b_below ! buoyancy above and below interface l

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      b_below = sol(S_TEMP,l)%data(d)%elts(id)   / (sol_mean(S_MASS,l)%data(d)%elts(id)   + sol(S_MASS,l)%data(d)%elts(id))
      b_above = sol(S_TEMP,l+1)%data(d)%elts(id) / (sol_mean(S_MASS,l+1)%data(d)%elts(id) + sol(S_MASS,l+1)%data(d)%elts(id))
      dzl = interp (dz(l), dz(l+1))

      eval = grav_accel * (b_above - b_below) / dzl ! - g drho/dz / rho0
    end function eval
  end function N_sq

  
  function Richardson (Nsq, dudzsq) result(val)
    ! Richardson number at interface  0 <= l <= zlevels
    
    implicit none
    
    real(dp), intent(in) :: Nsq, dudzsq
    real(dp)             :: val

   val= Nsq / (dudzsq + eps_s)
  end function Richardson

  
  function Prandtl (Ri) result(val)
    ! Computes Prandtl number given Richardson number
    
    implicit none
    
    real(dp), intent(in) :: Ri
    real(dp)             :: val

    real(dp) :: Ri_c ! critical Richardson number

    Ri_c = 2 / (2.0_dp + C_eps/C_k) 

    ! NEMO
    if (Ri < 0.2e0_dp) then
       val = 1.0_dp
    elseif (Ri >= 0.2_dp .and. Ri <= 2_dp) then
       val = 5 * Ri
    else
       val = 10.0_dp
    end if
!!$    val = max (0.1e0_dp, Ri_c / max (Ri_c, Ri)) ! CROCO
  end function Prandtl
  

  subroutine l_scales (dz, Nsq, tau, tke, l_eps, l_k)
    ! Computes length scales l_eps and l_m at interfaces 0:zlevels for TKE closure for a single vertical column
    
    implicit none
    
    real(dp),                       intent (in)  :: tau   ! wind stress
    real(dp), dimension(1:zlevels), intent (in)  :: dz    ! layer thicknesses
    real(dp), dimension(0:zlevels), intent (in)  :: Nsq   ! Brunt-Vaisala frequency
    real(dp), dimension(0:zlevels), intent (in)  :: tke   ! turbulent kinetic energy
    real(dp), dimension(0:zlevels), intent (out) :: l_k   ! dissipation length scale (velocity)
    real(dp), dimension(0:zlevels), intent (out) :: l_eps ! mixing length scale (buoyancy)

    integer  :: l
    real(dp) :: l_dwn(0:zlevels), l_up(0:zlevels)

    ! First order approximation for mixing length
    do l = 0, zlevels
       l_dwn(l) = sqrt (2 * tke(l) / max (Nsq(l), Neps_sq))
    end do
    l_up = l_dwn

    ! Upwards integration from bottom
    l_dwn(0) = l_0
    do l = 1, zlevels
       l_dwn(l) = min (l_dwn(l), l_dwn(l-1) + dz(l))
    end do

    ! Downward integration from surface
    l_up(zlevels) = max (l_0, kappa_VK * C_l / (ref_density * grav_accel) * tau)
    do l = zlevels-1, 0, -1
       l_up(l) = min (l_up(l), l_up(l+1) + dz(l+1)) 
    end do

    ! Returned mixing length scales
    l_eps = max (l_min, sqrt (l_up * l_dwn)) 
    l_k   = max (l_min, min  (l_up,  l_dwn))
  end subroutine l_scales

  
  function Kt_tke (Kv, Nsq, Ri) result(val)
    ! TKE closure eddy diffusivity

    implicit none

    real(dp), intent(in) :: Kv  ! eddy viscosity
    real(dp), intent(in) :: Nsq ! Brunt-Vaisala frequency squared
    real(dp), intent(in) :: Ri  ! Richardson number
    real(dp)             :: val

   val = min (Kt_max, max (Kv/Prandtl(Ri), Kt_min))

    if (enhance_diff .and. Nsq <= Nsq_min) val = Kt_enh ! enhanced vertical diffusion
  end function Kt_tke
  

  function Kv_tke (e, l_k, Nsq) result(val)
    ! TKE closure eddy viscosity
    implicit none
    real(dp), intent(in) :: e   ! tke
    real(dp), intent(in) :: l_k ! mixing length for eddy viscosity dissipation
    real(dp), intent(in) :: Nsq ! Brunt-Vaisala frequency squared
    real(dp)             :: val

    val = min (Kv_max, max (C_k * l_k * sqrt(e), Kv_min))
  end function Kv_tke

  
  function Kt_analytic (z, eta) result(val)
    ! Analytic eddy diffusivity
    
    implicit none
    
    real(dp), intent(in) :: eta, z
    real(dp)             :: val
    
    val = Kt_min + Kt_max * exp (-20 * (z - eta) / z_mixed)
  end function Kt_analytic


  function Kv_analytic (z, eta) result(val)
    ! Analytic eddy viscosity

    implicit none

    real(dp), intent(in) :: eta, z
    real(dp)             :: val

    val = Kv_min + Kv_max * exp (-20 * (z - eta) / z_mixed)
  end function Kv_analytic

  
  subroutine trend_vertical_diffusion (q, dq)
    ! Trend for eddy diffusivity and eddy viscosity for forward Euler time step
    
    implicit none
    
    type(Float_Field), target, intent(in)    :: q(1:N_VARIABLE,1:zmax)
    type(Float_Field), target, intent(inout) :: dq(1:N_VARIABLE,1:zmax)

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
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: d, id_i
    real(dp) :: dz_k

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    dmass(id_i) = 0.0_dp

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
    
    function scalar_flux (l) result(val)
      ! Computes flux at interface below (l=-1) or above (l=1) vertical level zlev
      
      implicit none
      
      integer, intent(in) :: l
      real(dp)            :: val

      real(dp) :: b_0, b_l, dz_l, mass_0, mass_l, temp_0, temp_l, visc

      visc = Kt(zlev+min(0,l))%data(d)%elts(id_i)

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dz_l = interp (dz_k,  dz_i(dom, i, j, zlev+l, offs, dims, sol)) ! thickness of layer centred on interface

      val = l * visc * (b_l - b_0) / dz_l
    end function scalar_flux
    
  end subroutine trend_scalars_vert_diffuse
  
  
  subroutine trend_velo_vert_diffuse (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: d, id
    real(dp) :: dz_k(EDGE)

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    dz_k = dz_e (dom, i, j, zlev, offs, dims, sol)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(id_edge(id)) = (velo_flux(1) - velo_flux(-1)) / dz_k
    elseif  (zlev == 1) then
       dvelo(id_edge(id)) = (velo_flux(1) - phi_edge (d, id, 1) * bottom_friction * velo(id_edge(id))) / dz_k
    elseif (zlev == zlevels) then
       dvelo(id_edge(id)) = (wind_flux (dom, i, j, z_null, offs, dims) - velo_flux(-1)) / dz_k 
    end if

  contains
    function velo_flux (l) result(val)
      ! Flux at upper interface (l=1) or lower interface (l=-1)

      implicit none

      integer, intent(in) :: l
      real(dp)            :: val(3)

      real(dp) :: dzl(3), visc(3)

      visc = Kv(zlev+min(0,l))%data(d)%elts(id+1)
      dzl = 0.5 * (dz_k + dz_e (dom, i, j, zlev+l, offs, dims, sol)) ! thickness of layer centred on interface

      val= l * visc * (sol(S_VELO,zlev+l)%data(d)%elts(id_edge(id)) - velo(id_edge(id))) / dzl
    end function velo_flux

  end subroutine trend_velo_vert_diffuse


end module vert_diffusion_mod
