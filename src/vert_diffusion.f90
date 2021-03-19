module vert_diffusion_mod
  ! Provides a backwards Euler step for vertical diffusion of buoyancy (temp variable) and
  ! velocity for ocean models. Backwards Euler is unconditionally stable.
  !
  ! implicit_vertical_diffusion = backwards Euler step
  ! trend_vertical_diffusion    = trend routines for forward Euler step
  !
  ! User must supply the following functions in test_case_mod.f90:
  !
  ! vertical boundary condition routines
  !        top_temp_source, bottom_temp_source
  !        wind_drag
  !
  ! scalar variable
  !        bottom_friction
  use utils_mod
  use test_case_mod
  implicit none
  real(8) :: friction
  
  ! Parameters for TKE closure eddy viscosity/diffusion scheme
  real(8), parameter :: c_e      = 1.0d0
  real(8), parameter :: c_eps    = 0.7d0
  real(8), parameter :: C_l      = 2d5
  real(8), parameter :: c_m      = 0.1d0
  real(8), parameter :: C_sfc    = 67.83d0
  real(8), parameter :: e_0      = 1d-6/sqrt(2.0d0) * METRE/SECOND**2
  real(8), parameter :: e_sfc_0  = 1d-4    * METRE/SECOND**2
  real(8), parameter :: eps_s    = 1d-20   * 1/SECOND**2
  real(8), parameter :: kappa_vk = 0.40    ! von Karman constant
  real(8), parameter :: K_m0     = 1.2d-4  * METRE**2/SECOND
  real(8), parameter :: K_t0     = 1.2d-5  * METRE**2/SECOND
  real(8), parameter :: l_0      = 0.04    * METRE
  real(8), parameter :: Neps_sq  = 1d-20   * 1/SECOND**2
  real(8), parameter :: Ri_c     = 2 / (2 + c_eps/c_m)

  ! Parameters for analytic eddy viscosity/diffusion scheme
  real(8), parameter :: D_max     = 150 * METRE
  real(8), parameter :: Km_bottom = 2d-3
  real(8), parameter :: Kt_const  = 1d-6

  real(8), parameter :: tau_mag = 0.1 ! hard code for now

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
  end interface
  procedure (fun3), pointer :: bottom_temp_source => null ()
  procedure (fun3), pointer :: top_temp_source    => null ()
  procedure (fun4), pointer :: wind_drag          => null ()
contains
  subroutine implicit_vertical_diffusion (r, wind_d, source_b, source_t)
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

    friction = r
    
    call update_array_bdry (sol, NONE, 27)

    ! Backwards Euler step for TKE closure to compute eddy viscosity and eddy diffusivity
    if (tke_closure) then
       do l = level_end, level_start, -1
          call apply_onescale (backwards_euler_tke, l, z_null, 0, 1)
       end do
       tke%bdry_uptodate = .false.
       call WT_after_scalar (tke, wav_tke, level_start-1)
    end if

    ! Backwards Euler step for each vertical colum
    do l = level_end, level_start, -1
       call apply_onescale (backwards_euler_temp, l, z_null, 0, 1)
       call apply_onescale (backwards_euler_velo, l, z_null, 0, 0)
    end do
    sol%bdry_uptodate = .false.
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine implicit_vertical_diffusion

  subroutine backwards_euler_tke (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for TKE closure equation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, info, k
    real(8)                         :: dz, eta, full_mass, theta
    real(8), dimension(0:zlevels)   :: K_t
    real(8), dimension(1:zlevels)   :: diag, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u

    id_i = idx (i, j, offs, dims) + 1
    d = dom%id + 1

    ! Bottom layer
    k = 1
    dz = dz_i (dom, i, j, k, offs, dims)
    diag_u(k) = - coeff(1) ! super-diagonal
    diag(k)   = 1.0_8 - diag_u(k)
    rhs(k)    = b() + dt * flux_mean(1)/dz + dt * bottom_temp_source (dom, i, j, z_null, offs, dims)
    
    do k = 2, zlevels-1
       dz = dz_i (dom, i, j, k, offs, dims)
       diag_u(k)   = - coeff( 1) ! super-diagonal
       diag_l(k-1) = - coeff(-1) ! sub-diagonal
       diag(k)     = 1.0_8  - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b() + dt * (flux_mean(1) - flux_mean(-1))/dz
    end do

    ! Top layer
    k = zlevels
    dz  = dz_i (dom, i, j, k, offs, dims)
    diag_l(k-1) = - coeff(-1) ! sub-diagonal
    diag(k)     = 1.0_8 - diag_u(k-1)
    rhs(k)      = b() - dt * flux_mean(-1)/dz + dt * top_temp_source (dom, i, j, z_null, offs, dims)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed with info = ", info

    ! Backwards Euler step
    do k = 1, zlevels-1
       !tke(k)%data(d)%elts(id_i) = 
    end do
    tke(zlevels)%data(d)%elts(id_i) = max (C_sfc * tau_mag / ref_density, e_sfc_0) !!! use tau_0 for now
  contains
    real(8) function coeff (l)
      ! Computes coefficient above (l = 1) or below (l = -1) for vertical Laplacian matrix
      implicit none
      integer :: l
      
      real(8) :: dz_l

      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims))
      coeff = dt / (dz_l * dz) * K_t (k+l)
    end function coeff
    
    real(8) function b ()
      ! Fluctuating buoyancy
      implicit none
      real(8) :: full_mass, full_theta
      full_mass  = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i)
      full_theta = sol_mean(S_TEMP,k)%data(d)%elts(id_i) + sol(S_TEMP,k)%data(d)%elts(id_i)
      b = full_theta / full_mass - b_mean()
    end function b
    
    real(8) function b_mean ()
      ! Mean buoyancy
      implicit none
      b_mean = sol_mean(S_TEMP,k)%data(d)%elts(id_i) / sol_mean(S_MASS,k)%data(d)%elts(id_i) 
    end function b_mean
    
    real(8) function flux_mean (l)
      ! Computes flux of mean buoyancy at interface below (l=-1) or above (l=1) vertical level k
      implicit none
      integer :: l

      real(8) :: b_0, b_l, dz_l, mass_0, mass_l, temp_0, temp_l

      mass_0 = sol_mean(S_MASS,k)%data(d)%elts(id_i) 
      temp_0 = sol_mean(S_TEMP,k)%data(d)%elts(id_i) 
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,k+l)%data(d)%elts(id_i) 
      temp_l = sol_mean(S_TEMP,k+l)%data(d)%elts(id_i) 
      b_l = temp_l / mass_l

      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

      flux_mean = l * K_t(k+l) * (b_l - b_0) / dz_l
    end function flux_mean
  end subroutine backwards_euler_tke

  subroutine backwards_euler_temp (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for temp variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, info, k
    real(8)                         :: dz, eta, full_mass, theta
    real(8), dimension(0:zlevels)   :: K_t
    real(8), dimension(1:zlevels)   :: diag, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u

    id_i = idx (i, j, offs, dims) + 1
    d = dom%id + 1

    if (.not. tke_closure) then
       K_t = Kt_analytic ()
    end if
    
    if (mode_split) then
       eta = sol(S_MASS,zlevels+1)%data(d)%elts(id_i)
    else
       eta = free_surface (dom, i, j, z_null, offs, dims)
    end if

    ! Bottom layer
    k = 1
    dz = dz_i (dom, i, j, k, offs, dims)
    diag_u(k) = - coeff(1) ! super-diagonal
    diag(k)   = 1.0_8 - diag_u(k)
    rhs(k)    = b() + dt * flux_mean(1)/dz + dt * bottom_temp_source (dom, i, j, z_null, offs, dims)
    
    do k = 2, zlevels-1
       dz = dz_i (dom, i, j, k, offs, dims)
       diag_u(k)   = - coeff( 1) ! super-diagonal
       diag_l(k-1) = - coeff(-1) ! sub-diagonal
       diag(k)     = 1.0_8  - (diag_u(k) + diag_l(k-1))
       rhs(k)      = b() + dt * (flux_mean(1) - flux_mean(-1))/dz
    end do

    ! Top layer
    k = zlevels
    dz  = dz_i (dom, i, j, k, offs, dims)
    diag_l(k-1) = - coeff(-1) ! sub-diagonal
    diag(k)     = 1.0_8 - diag_u(k-1)
    rhs(k)      = b() - dt * flux_mean(-1)/dz + dt * top_temp_source (dom, i, j, z_null, offs, dims)

    ! Solve tridiagonal linear system
    call dgtsv (zlevels, 1, diag_l, diag, diag_u, rhs, zlevels, info)
    if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed with info = ", info

    ! Backwards Euler step
    do k = 1, zlevels
       full_mass = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i)
       sol(S_TEMP,k)%data(d)%elts(id_i) = full_mass * (b_mean() + rhs(k)) - sol_mean(S_TEMP,k)%data(d)%elts(id_i)
    end do
  contains
    real(8) function coeff (l)
      ! Computes coefficient above (l = 1) or below (l = -1) for vertical Laplacian matrix
      implicit none
      integer :: l
      
      real(8) :: dz_l

      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims))
      coeff = dt / (dz_l * dz) * K_t (k+l)
    end function coeff
    
    real(8) function b ()
      ! Fluctuating buoyancy
      implicit none
      real(8) :: full_mass, full_theta
      full_mass  = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i)
      full_theta = sol_mean(S_TEMP,k)%data(d)%elts(id_i) + sol(S_TEMP,k)%data(d)%elts(id_i)
      b = full_theta / full_mass - b_mean()
    end function b
    
    real(8) function b_mean ()
      ! Mean buoyancy
      implicit none
      b_mean = sol_mean(S_TEMP,k)%data(d)%elts(id_i) / sol_mean(S_MASS,k)%data(d)%elts(id_i) 
    end function b_mean
    
    real(8) function flux_mean (l)
      ! Computes flux of mean buoyancy at interface below (l=-1) or above (l=1) vertical level k
      implicit none
      integer :: l

      real(8) :: b_0, b_l, dz_l, mass_0, mass_l, temp_0, temp_l

      mass_0 = sol_mean(S_MASS,k)%data(d)%elts(id_i) 
      temp_0 = sol_mean(S_TEMP,k)%data(d)%elts(id_i) 
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,k+l)%data(d)%elts(id_i) 
      temp_l = sol_mean(S_TEMP,k+l)%data(d)%elts(id_i) 
      b_l = temp_l / mass_l

      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

      flux_mean = l * K_t(k+l) * (b_l - b_0) / dz_l
    end function flux_mean
  end subroutine backwards_euler_temp

  subroutine backwards_euler_velo (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for velocity variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                                :: d, e, id, idE, idN, idNE, info, k
    real(8), dimension(1:EDGE)             :: dz_k, eta, z
    real(8), dimension(1:EDGE,0:zlevels)   :: K_m
    real(8), dimension(1:EDGE,1:zlevels)   :: diag, rhs
    real(8), dimension(1:EDGE,1:zlevels-1) :: diag_l, diag_u
    real(8), dimension(1:zlevels)          :: dd, r
    real(8), dimension(1:zlevels-1)        :: dl, du

    d = dom%id + 1

    id   = idx (i, j,     offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    eta = eta_e (dom, i, j, z_null, offs, dims)

    if (.not. tke_closure) then
       z(RT+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idE+1))  
       z(DG+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idNE+1)) 
       z(UP+1) = interp (dom%topo%elts(id+1), dom%topo%elts(idN+1))
       K_m(:,0) = Km_analytic (z, eta)
       do k = 1, zlevels
          z = z + dz_e (dom, i, j, k, offs, dims)
          K_m(:,k) = Km_analytic (z, eta)
       end do
    end if

    ! Bottom layer
    k = 1
    dz_k = dz_e (dom, i, j, k, offs, dims)
    diag_u(:,k) = - coeff(1) ! super-diagonal
    diag(:,k)   = 1.0_8 - diag_u(:,k) + dt * friction/dz_k
    rhs(:,k)    = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) 

    do k = 2, zlevels-1
       dz_k = dz_e (dom, i, j, k, offs, dims)
       diag_u(:,k)   = - coeff( 1) ! super-diagonal
       diag_l(:,k-1) = - coeff(-1) ! sub-diagonal
       diag(:,k)     = 1.0_8 - (diag_u(:,k) + diag_l(:,k-1))
       rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    end do
    
    ! Top layer
    k = zlevels
    dz_k = dz_e (dom, i, j, k, offs, dims)
    diag_l(:,k-1) = - coeff(-1) ! sub-diagonal
    diag(:,k)     = 1.0_8 - diag_l(:,k-1)
    rhs(:,k)      = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) + dt * wind_drag (dom, i, j, z_null, offs, dims)

    ! Solve tridiagonal linear system
    do e = 1, EDGE
       dl = diag_l(e,:); dd = diag(e,:); du = diag_u(e,:); r = rhs(e,:)
       call dgtsv (zlevels, 1, dl, dd, du, r, zlevels, info)
       if (info /= 0 .and. rank == 0) write (6,'(a,i2)') "Warning: dgtsv failed with info = ", info
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
            
      real(8), dimension(1:EDGE) :: dz_l
      
      dz_l = 0.5 * (dz_k + dz_e (dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

      coeff = dt / (dz_l * dz_k) * K_m(:,k+l)
    end function coeff
  end subroutine backwards_euler_velo

  real(8) function N_sq  (dom, i, j, zlev, offs, dims, l)
    ! Brunt-Vaisala number N^2 at interface below (l=-1) or above (l=1) level zlev
    ! (defined at nodes)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i
    real(8) :: drho, mass_0, mass_l, rho_0, rho_l, temp_0, temp_l
    
    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1

    mass_0 = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)
    temp_0 = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)
    rho_0 = (1.0_8 - temp_0 / mass_0) * porous_density (dom, i, j, zlev, offs, dims)

    mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
    temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
    rho_l = (1.0_8 - temp_l / mass_l) * porous_density (dom, i, j, zlev+l, offs, dims)

    drho = rho_l - rho_0

    N_sq = - grav_accel * drho / ref_density
  end function N_sq

  real(8) function dudz_sq  (dom, i, j, zlev, offs, dims, l)
    ! Kinetic at interface below (l=-1) or above (l=1) level zlev
    ! (defined at nodes)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                      :: d, id, idS, idSW, idW
    real(8), dimension(1:2*EDGE) :: du, dz, prim_dual

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i,   j-1, offs, dims)
    idS  = idx (i-1, j-1, offs, dims)

    ! Energy term (using TRiSK form of kinetic energy)
    du(1:EDGE) = sol(S_VELO,zlev+l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) &
         - sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    du(EDGE+1) = sol(S_VELO,zlev+l)%data(d)%elts(EDGE*idW+RT+1)  - sol(S_VELO,zlev)%data(d)%elts(EDGE*idW+RT+1)
    du(EDGE+2) = sol(S_VELO,zlev+l)%data(d)%elts(EDGE*idSW+DG+1) - sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1)
    du(EDGE+3) = sol(S_VELO,zlev+l)%data(d)%elts(EDGE*idS+UP+1)  - sol(S_VELO,zlev)%data(d)%elts(EDGE*idS+UP+1)

    prim_dual(1:EDGE) = dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1)
    prim_dual(EDGE+1) = dom%len%elts(EDGE*idW+RT+1)  * dom%pedlen%elts(EDGE*idW+RT+1)
    prim_dual(EDGE+2) = dom%len%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)
    prim_dual(EDGE+3) = dom%len%elts(EDGE*idS+UP+1)  * dom%pedlen%elts(EDGE*idS+UP+1)

    ! Layer thicknesses centred on interface
    dz(1:EDGE)        = interp_e (dz_e (dom, i, j, zlev, offs, dims), dz_e (dom, i, j, zlev+l, offs, dims)) 
    dz(EDGE+1:2*EDGE) = interp_e (dz_SW_e (dom, i, j, zlev, offs, dims), dz_SW_e (dom, i, j, zlev+l, offs, dims)) 

    dudz_sq = sum ((du/dz)**2 * prim_dual) * dom%areas%elts(id+1)%hex_inv/2
  end function dudz_sq 

  real(8) function richardson (dom, i, j, zlev, offs, dims, l)
    ! Richardson number at interface below (l=-1) or above (l=1) level zlev
    ! (defined at nodes)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer                        :: l
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    real(8) :: KE_dudz, N2

    KE_dudz = dudz_sq (dom, i, j, zlev, offs, dims, l)

    if (KE_dudz == 0.0_8) then
       richardson = 1d10
    else
       N2 = N_sq  (dom, i, j, zlev, offs, dims, l)
       richardson =  l * N2 / (KE_dudz + eps_s)
    end if
  end function richardson

  real(8) function prandtl (Ri)
    ! Computes Prandtl number given the Richardson number
    implicit none
    real(8) :: Ri

    prandtl = max (0.1d0, Ri_c / max (Ri_c, Ri))
  end function prandtl

  subroutine l_scales (dz, Nsq, tau, tke)
    ! Computes length scales for TKE closure for a single vertical column
    implicit none
    real(8)                       :: tau
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: l_dwn, l_up, Nsq, tke

    integer :: k

    do k = 0, zlevels
       l_up(k) = sqrt (2 * tke(k) / max (Nsq(k), Neps_sq))
    end do
    l_dwn = l_up

    l_dwn(0) = l_0
    do k = 1, zlevels
       l_dwn(k) = min (l_dwn(k-1) + dz(k), l_dwn(k))
    end do

    l_up(zlevels) = kappa_vk*C_l / (ref_density * grav_accel) * tau
    do k = zlevels-1, 0, -1
       l_up(k) = min (l_up(k+1) + dz(k+1), l_up(k)) 
    end do
  end subroutine l_scales

  real(8) function Km_tke (l_m, tke)
    ! TKE closure eddy viscosity
    implicit none
    real(8) :: l_m, tke

    Km_tke = max (c_m * l_m * sqrt(tke), K_m0)
  end function Km_tke
  
  function Km_analytic (z, eta)
    ! Analytic eddy viscosity
    real(8), dimension(1:EDGE) :: Km_analytic
    real(8), dimension(1:EDGE) :: eta, z

    Km_analytic = Km_bottom * (1.0_8 + 4 * exp ((z - eta) / D_max ))
  end function Km_analytic

  real(8) function Kt_tke (Km, Pr)
    ! TKE closure eddy diffusivity
    implicit none
    real(8) :: Km, Pr
    
    Kt_tke = max (Km/Pr, K_t0)
  end function Kt_tke

  real(8) function Kt_analytic ()
    ! Analytic eddy diffusivity
    implicit none

    Kt_analytic = Kt_const
  end function Kt_analytic
end module vert_diffusion_mod
