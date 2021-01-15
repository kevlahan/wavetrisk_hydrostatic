module vert_diffusion_mod
  ! Provides a backwards Euler step and trend routines for forward Euler for vertical diffusion of buoyancy (temp variable) and
  ! velocity for ocean models. Forward Euler requires nu dt/dz_min^2 <= 1/2 for stability. Backwards Euler is unconditionally stable.
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
  use test_case_mod
  implicit none
contains
  subroutine implicit_vertical_diffusion
    ! Backwards euler step for vertical diffusion
    use adapt_mod
    implicit none
    integer :: d, p

    call update_array_bdry (sol, NONE, 27)
    
    ! Backwards Euler step for each vertical colum
    do d = 1, size(grid)
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (backwards_euler_temp, grid(d), p-1, z_null, 0, 1)
          call apply_onescale_to_patch (backwards_euler_velo, grid(d), p-1, z_null, 0, 0)
       end do
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
    eta = sol(S_MASS,zlevels+1)%data(d)%elts(id_i)

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
    if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed with info = ", info

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
      dz_l = 0.5 * (dz + dz_i(dom, i, j, k+l, offs, dims))
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
      dz_l = 0.5 * (dz + dz_i(dom, i, j, k+l, offs, dims)) ! thickness of layer centred on interface

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
    diag(:,k)   = 1.0_8 - diag_u(:,k) + dt * bottom_friction/dz_k
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
       if (info /= 0) write (6,'(a,i2)') "Warning: dgtsv failed with info = ", info
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
  
  subroutine trend_vertical_diffusion (q, dq)
    ! Trend for eddy diffusivity and eddy viscosity 
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q, dq

    integer :: d, k, p

    call update_array_bdry (q, NONE, 27)

    ! Scalars
    do d = 1, size(grid)
       scalar => q(S_MASS,zlevels+1)%data(d)%elts ! free surface
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
             call apply_onescale_to_patch (trend_scalars, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_velo,    grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dtemp, dvelo, mass, mean_m, mean_t, temp, velo)
       end do
       nullify (scalar)
    end do
    dq%bdry_uptodate = .false.
  end subroutine trend_vertical_diffusion

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Vertical eddy diffusivity of buoyancy
    ! (layer height is not diffused)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer            :: id_i
    real(8)            :: dz_k
    
    id_i = idx (i, j, offs, dims) + 1
    dz_k = dz_i (dom, i, j, zlev, offs, dims)
    
    if (zlev > 1 .and. zlev < zlevels) then
       dtemp(id_i) = flux_temp(1) - flux_temp(-1)
    elseif (zlev == 1) then
       dtemp(id_i) =  flux_temp(1) + bottom_temp_source(dom, i, j, z_null, offs, dims)
    elseif (zlev == zlevels) then
       dtemp(id_i) =  top_temp_source (dom, i, j, z_null, offs, dims) - flux_temp(-1)
    end if
    dtemp(id_i) = porous_density (dom, i, j, zlev, offs, dims) * dtemp(id_i)
  contains
    real(8) function flux_temp (l)
      ! Computes flux at interface below (l=-1) or above (l=1) vertical level zlev
      implicit none
      integer :: l

      integer :: d
      real(8) :: b_0, b_l, dz_l, eta, mass_0, mass_l, ri, temp_0, temp_l, visc, z

      d = dom%id + 1

      eta = scalar(id_i)
      ri  = richardson (dom, i, j, zlev, offs, dims, l)
      z   = zl_i (dom, i, j, zlev, offs, dims, l)
      visc = eddy_diffusivity (eta, ri, z)

      mass_0 = mean_m(id_i) + mass(id_i)
      temp_0 = mean_t(id_i) + temp(id_i)
      b_0 = temp_0 / mass_0

      mass_l = sol_mean(S_MASS,zlev+l)%data(d)%elts(id_i) + sol(S_MASS,zlev+l)%data(d)%elts(id_i)
      temp_l = sol_mean(S_TEMP,zlev+l)%data(d)%elts(id_i) + sol(S_TEMP,zlev+l)%data(d)%elts(id_i)
      b_l = temp_l / mass_l

      dz_l = 0.5 * (dz_k + dz_i(dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      flux_temp = l * visc * (b_l - b_0) / dz_l
    end function flux_temp
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id
    real(8), dimension(1:EDGE) :: dz_k, eta

    id = idx (i, j, offs, dims)
    
    dz_k = dz_e (dom, i, j, zlev, offs, dims)
    eta = eta_e (dom, i, j, zlev, offs, dims)

    if (zlev > 1 .and. zlev < zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = (flux_velo(1) - flux_velo(-1)) / dz_k
    elseif  (zlev == 1) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = flux_velo(1) / dz_k - bottom_friction * velo(EDGE*id+RT+1:EDGE*id+UP+1)
    elseif (zlev == zlevels) then
       dvelo(EDGE*id+RT+1:EDGE*id+UP+1) = wind_drag (dom, i, j, z_null, offs, dims) - flux_velo(-1) / dz_k 
    end if
  contains
    function flux_velo (l)
      ! Flux at upper interface (l=1) or lower interface (l=-1)
      implicit none
      integer                    :: l
      real(8), dimension(1:EDGE) :: flux_velo

      integer                    :: d
      real(8)                    :: ri
      real(8), dimension(1:EDGE) :: dz_l, visc, z
      
      d = dom%id + 1
      
      ri  = richardson (dom, i, j, zlev, offs, dims, l)
      z   = zl_e  (dom, i, j, zlev, offs, dims, l)
      visc = eddy_viscosity (eta, ri, z)

      dz_l = 0.5 * (dz_k + dz_e (dom, i, j, zlev+l, offs, dims)) ! thickness of layer centred on interface

      flux_velo = l * visc * (sol(S_VELO,zlev+l)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) - velo(EDGE*id+RT+1:EDGE*id+UP+1)) / dz_l
    end function flux_velo
  end subroutine trend_velo
end module vert_diffusion_mod
