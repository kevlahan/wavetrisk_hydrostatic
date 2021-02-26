module vert_diffusion_mod
  ! Provides a backwards Euler step for vertical diffusion of buoyancy (temp variable) and
  ! velocity for ocean models. Backwards Euler is unconditionally stable.
  !
  ! implicit_vertical_diffusion = backwards Euler step
  ! trend_vertical_diffusion    = trend routines for forward Euler step
  !
  ! User must supply the following functions in test_case_mod.f90:
  !
  ! viscosity routines
  !        eddy_viscosity, eddy_diffusion
  !
  ! vertical boundary condition routines
  !        top_temp_source, bottom_temp_source
  !        wind_drag
  !
  ! scalar variable
  !        bottom_friction
  use utils_mod
  implicit none
  real(8) :: friction
  
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
  procedure (fun1), pointer :: eddy_diffusivity   => null ()
  procedure (fun2), pointer :: eddy_viscosity     => null ()
  procedure (fun3), pointer :: bottom_temp_source => null ()
  procedure (fun3), pointer :: top_temp_source    => null ()
  procedure (fun4), pointer :: wind_drag          => null ()
contains
  subroutine implicit_vertical_diffusion (eddy_d, eddy_v, r, wind_d, source_b, source_t)
    ! Backwards euler step for vertical diffusion
    use adapt_mod
    implicit none
    integer :: l
    real(8) :: r

    interface
       real(8) function eddy_d (eta, ri, z)
         implicit none
         real(8) :: eta, ri, z
       end function eddy_d
       function eddy_v (eta, ri, z)
         use shared_mod
         implicit none
         real(8), dimension(1:EDGE) :: eddy_v, eta, z
         real(8)                    :: ri
       end function eddy_v
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

    eddy_diffusivity   => eddy_d
    eddy_viscosity     => eddy_v
    bottom_temp_source => source_b
    top_temp_source    => source_t
    wind_drag          => wind_d

    friction = r

    call update_array_bdry (sol, NONE, 27)

    ! Backwards Euler step for each vertical colum
    do l = level_end, level_start, -1
       call apply_onescale (backwards_euler_temp, l, z_null, 0, 1)
       call apply_onescale (backwards_euler_velo, l, z_null, 0, 0)
    end do
    sol%bdry_uptodate = .false.
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine implicit_vertical_diffusion

  subroutine backwards_euler_temp (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for temp variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                         :: d, id_i, info, k
    real(8)                         :: dz, eta, full_mass, theta
    real(8), dimension(1:zlevels)   :: diag, rhs
    real(8), dimension(1:zlevels-1) :: diag_l, diag_u

    id_i = idx (i, j, offs, dims) + 1
    d = dom%id + 1
    
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
      
      real(8) :: dz_l, ri, z

      ri  = richardson (dom, i, j, k, offs, dims, l)
      z   = zl_i (dom, i, j, k, offs, dims, l)
      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims))
      coeff = dt / (dz_l * dz) * eddy_diffusivity (eta, ri, z) 
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

      real(8) :: b_0, b_l, dz_l, mass_0, mass_l, ri, temp_0, temp_l, z

      mass_0 = sol_mean(S_MASS,k)%data(d)%elts(id_i) 
      temp_0 = sol_mean(S_TEMP,k)%data(d)%elts(id_i) 
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,k+l)%data(d)%elts(id_i) 
      temp_l = sol_mean(S_TEMP,k+l)%data(d)%elts(id_i) 
      b_l = temp_l / mass_l

      ri  = richardson (dom, i, j, k, offs, dims, l)
      z   = zl_i (dom, i, j, k, offs, dims, l)
      dz_l = interp (dz, dz_i(dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

      flux_mean = l * eddy_diffusivity (eta, ri, z) * (b_l - b_0) / dz_l
    end function flux_mean
  end subroutine backwards_euler_temp

  subroutine backwards_euler_velo (dom, i, j, z_null, offs, dims)
    ! Backwards Euler step for velocity variable
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                                :: d, e, id, info, k
    real(8), dimension(1:EDGE)             :: dz_k, eta
    real(8), dimension(1:EDGE,1:zlevels)   :: diag, rhs
    real(8), dimension(1:EDGE,1:zlevels-1) :: diag_l, diag_u
    real(8), dimension(1:zlevels)          :: dd, r
    real(8), dimension(1:zlevels-1)        :: dl, du

    id = idx (i, j, offs, dims)
    d = dom%id + 1
    eta = eta_e (dom, i, j, z_null, offs, dims)

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
            
      real(8)                    :: ri
      real(8), dimension(1:EDGE) :: dz_l, z
      
      ri  = richardson (dom, i, j, k, offs, dims, l)
      z   = zl_e  (dom, i, j, k, offs, dims, l)
      dz_l = 0.5 * (dz_k + dz_e (dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

      coeff = dt / (dz_l * dz_k) * eddy_viscosity (eta, ri, z)
    end function coeff
  end subroutine backwards_euler_velo
end module vert_diffusion_mod
