module barotropic_2d_mod
  ! Files needed to solve barotropic free surface
  use ops_mod
  use multi_level_mod
  use lin_solve_mod
  use utils_mod
  implicit none
  ! Add Laplacian diffusion to free surface perturbation eta
  real(8), parameter :: C_eta = 5d-3
  logical, parameter :: diff_eta = .true.
contains
  subroutine scalar_star (dt, q)
    ! Explicit Euler step for scalars
    implicit none
    real(8)                                   :: dt
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, ibeg, iend, k, v

    do d = 1, size(grid)
       ibeg = (1+2*(POSIT(S_MASS)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend = sol(S_MASS,1)%data(d)%length
       do k = 1, zlevels
          do v = scalars(1), scalars(2)
             q(v,k)%data(d)%elts(ibeg:iend) = sol(v,k)%data(d)%elts(ibeg:iend) + dt * trend(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
    end do
    q(scalars(1):scalars(2),1:zlevels)%bdry_uptodate = .false.
  end subroutine scalar_star

  subroutine u_star (dt, q)
    ! Explicit Euler step for intermediate velocity u_star
    ! remove external pressure gradient
    implicit none
    real(8)                                   :: dt
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, ibeg, iend, k

    ! External pressure gradient
    call grad_eta

    do d = 1, size(grid)
       ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend = q(S_VELO,1)%data(d)%length
       do k = 1, zlevels
          q(S_VELO,k)%data(d)%elts(ibeg:iend) = sol(S_VELO,k)%data(d)%elts(ibeg:iend) &
               + dt * (trend(S_VELO,k)%data(d)%elts(ibeg:iend) + theta1 * horiz_flux(S_TEMP)%data(d)%elts(ibeg:iend))
       end do
    end do
    q(S_VELO,1:zlevels)%bdry_uptodate = .false.
  end subroutine u_star

  subroutine barotropic_correction (q)
    ! Update baroclinic variables mass and mass-weighted buoyancy with new free surface perturbation
    ! uses Bleck and Smith (J. Geophys. Res. 95, 3273–3285 1990) layer dilation method
    ! NOTE: individual layers no longer conserve mass (although total mass is conserved)
    implicit none
    type(Float_Field), dimension(:,:), target :: q
    
    integer :: d, k, p

    ! Sum mass perturbations
    call total_height (q(S_MASS,1:zlevels), exner_fun(1))

    do d = 1, size(grid)
       scalar    => sol(S_MASS,zlevels+1)%data(d)%elts ! free surface perturbation
       scalar_2d => exner_fun(1)%data(d)%elts          ! sum of mass perturbations
       do k = 1, zlevels
          mass   => q(S_MASS,k)%data(d)%elts
          temp   => q(S_TEMP,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (cal_barotropic_correction, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, temp, mean_m, mean_t)
       end do
       nullify (scalar, scalar_2d)
    end do
    q(scalars(1):scalars(2),1:zlevels)%bdry_uptodate = .false.
  end subroutine barotropic_correction

  subroutine cal_barotropic_correction (dom, i, j, zlev, offs, dims)
    ! Correct baroclinic mass and buoyancy based on baroclinic estimate of free surface using layer dilation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: eta, full_mass, mean_theta, theta

    id = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id) >= ADJZONE) then
       d = dom%id + 1
       full_mass = mean_m(id) + mass(id)
       theta = (mean_t(id) + temp(id)) / full_mass ! full buoyancy
       
       ! Correct mass perturbation
       eta = scalar(id) / phi_node (d, id, zlevels) ! free surface perturbation
       mass(id) = (eta - grid(d)%topo%elts(id)) / scalar_2d(id) * full_mass - mean_m(id)

       ! Correct mass-weighted buoyancy
       temp(id) = (mean_m(id) + mass(id)) * theta - mean_t(id)
!!$       temp(id) = mass(id) * mean_t(id) / mean_m(id) ! assume time-independent buoyancy (e.g. no remap, constant density in each layer)
    end if
  end subroutine cal_barotropic_correction

  subroutine eta_update
    ! Theta step for free surface update
    use lin_solve_mod
    implicit none

    ! RHS of elliptic equation
    call rhs_elliptic

    ! Solve elliptic equation 
    call multiscale (sol(S_MASS,zlevels+1), sol(S_TEMP,zlevels+1), elliptic_lo)

    ! Diffuse free surface to increase stability and avoid discontinuities due to wave steepening
    if (diff_eta) then
       call update_bdry (sol(S_MASS,zlevels+1), NONE, 600)
       call diffuse_eta
    end if
  end subroutine eta_update

  subroutine diffuse_eta
    ! Explicit Laplacian diffusion split step for free surface
    implicit none

    integer :: d, ibeg, iend
    real(8) :: dt_nu

    call Laplacian_eta

    dt_nu = C_eta * dx_min**2
    do d = 1, size(grid)
       ibeg = (1+2*(POSIT(S_MASS)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend = sol(S_MASS,zlevels+1)%data(d)%length
       sol(S_MASS,zlevels+1)%data(d)%elts(ibeg:iend) = sol(S_MASS,zlevels+1)%data(d)%elts(ibeg:iend) &
            + dt_nu * Laplacian_scalar(S_MASS)%data(d)%elts(ibeg:iend)
    end do
    sol(S_MASS,zlevels+1)%bdry_uptodate = .false.
  end subroutine diffuse_eta

  subroutine u_update
    ! Explicit Euler velocity update with new external pressure gradient
    ! penalization is advanced using a backwards Euler scheme
    implicit none
    integer :: d, ibeg, iend, k, l
    
    ! External pressure gradient
    call grad_eta
    
    do d = 1, size(grid)
       ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend = sol(S_VELO,1)%data(d)%length
       do k = 1, zlevels
          sol(S_VELO,k)%data(d)%elts(ibeg:iend) = (sol(S_VELO,k)%data(d)%elts(ibeg:iend) &
               - theta1 * dt * horiz_flux(S_TEMP)%data(d)%elts(ibeg:iend)) 
       end do
    end do
    sol(S_VELO,1:zlevels)%bdry_uptodate = .false.
  end subroutine u_update

  subroutine rhs_elliptic 
    ! Forms rhs of elliptic equation for free surface, -eta^* in q(S_TEMP_zlevels+1)
    implicit none
    integer :: d, j, l

    ! Flux divergence of vertically integrated velocity u_star, stored in trend(S_MASS, zlevels+1)
    call flux_divergence (sol, trend(S_MASS,zlevels+1))
    
    ! RHS of elliptic equation, -eta^*
    do l = level_end, level_start, -1
       do d = 1, size(grid)
          dscalar => trend(S_MASS,zlevels+1)%data(d)%elts  ! flux divergence at intermediate time step, div F^*
          dmass   => trend(S_TEMP,zlevels+1)%data(d)%elts  ! flux divergence at previous time step, div F^n
          mass    =>   sol(S_MASS,zlevels+1)%data(d)%elts
          mass1   =>   sol(S_TEMP,zlevels+1)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_rhs_elliptic, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dmass, dscalar, mass, mass1)
       end do
    end do
    sol(S_TEMP,zlevels+1)%bdry_uptodate = .false.
  end subroutine rhs_elliptic

  subroutine cal_rhs_elliptic (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
        
    mass1(id) = - mass(id) + dt * (theta2 * dscalar(id) + (1d0 - theta2) * dmass(id)) / ref_density
  end subroutine cal_rhs_elliptic

  function elliptic_lo (q, l)
    ! Calculates linear operator L(eta) for barotropic elliptic equation for free surface perturbation at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_lo, q

    integer :: d, j

    call update_bdry (q, l, 33)

    elliptic_lo = q
    ! Calculate external pressure gradient flux
    do d = 1, size(grid)
       h_flux => horiz_flux(S_MASS)%data(d)%elts
       scalar => q%data(d)%elts
       mass => sol(S_MASS,zlevels+1)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=2)
       end do
       nullify (h_flux, mass, scalar)
    end do
    horiz_flux(S_MASS)%bdry_uptodate = .false.
    call update_bdry (horiz_flux(S_MASS), l, 213)

    ! Calculate divergence
    do d = 1, size(grid)
       dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
       h_flux  => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, h_flux)
    end do
    Laplacian_scalar(S_MASS)%bdry_uptodate = .false.
    call update_bdry (Laplacian_scalar(S_MASS), l, 101)

    ! Form complete linear operator 
    do d = 1, size(grid)
       dscalar => elliptic_lo%data(d)%elts
       mass    => q%data(d)%elts
       h_flux  => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (complete_elliptic_lo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, mass, h_flux)
    end do
    elliptic_lo%bdry_uptodate = .false.
  end function elliptic_lo

  subroutine complete_elliptic_lo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    dscalar(id) = theta1*theta2 * dt**2 * Laplacian_scalar(S_MASS)%data(dom%id+1)%elts(id) - mass(id)
  end subroutine complete_elliptic_lo

  subroutine flux_divergence (q, div_flux)
    ! Returns flux divergence of vertical integrated velocity in divF using solution q, stored in div_flux
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q
    type(Float_Field),                                 target :: div_flux

    integer :: d, j, l

    call update_vector_bdry (q(S_MASS,1:zlevels), NONE, 50)
    call update_vector_bdry (q(S_VELO,1:zlevels), NONE, 50)

    do l = level_end, level_start, -1
       ! Calculate vertically integrated velocity flux
       do d = 1, size(grid)
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (q=q, dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=4)
          end do
          if (l < level_end) then
             dscalar => div_flux%data(d)%elts
             call cpt_or_restr_flux (grid(d), l) ! restrict flux if possible
             nullify (dscalar)
          end if
          nullify (h_flux)
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 211)

       ! Calculate divergence of vertically integrated velocity flux
       do d = 1, size(grid)
          dscalar => div_flux%data(d)%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       div_flux%bdry_uptodate = .false.
       call update_bdry (div_flux, l, 212)
    end do
  end subroutine flux_divergence

  subroutine Laplacian_eta
    ! Computes Laplacian diffusion of free surface, stored in Laplacian_scalar(S_TEMP)
    implicit none
    
    integer :: d, j, l

    do l = level_end, level_start, -1
       do d = 1, size(grid)
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          scalar => sol(S_MASS,zlevels+1)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
          end do
          nullify (scalar)
          if (l < level_end) then
             dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l) ! restrict flux if possible
             nullify (dscalar)
          end if
          nullify (h_flux)
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       call update_bdry (horiz_flux(S_MASS), l, 211)

       ! Calculate divergence of vertically integrated velocity flux
       do d = 1, size(grid)
          dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       Laplacian_scalar(S_MASS)%bdry_uptodate = .false.
       call update_bdry (Laplacian_scalar(S_MASS), l, 212)
    end do
  end subroutine Laplacian_eta

  subroutine total_height (q, q_2d)
    ! Total height q_2d computed from pseudo-densities q
    implicit none
    type(Float_Field),               target :: q_2d
    type(Float_Field), dimension(:), target :: q

    integer :: d, k, p
    
    do d = 1, size(grid)
       scalar_2d => q_2d%data(d)%elts
       q_2d%data(d)%elts = 0.0_8
       do k = 1, zlevels
          mass   => q(k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (cal_height, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, mean_m)
       end do
       nullify (scalar_2d)
    end do
  end subroutine total_height

  subroutine cal_height (dom, i, j, zlev, offs, dims)
    ! Vertical integration of edge quantity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: dz
    
    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1
    
    if (dom%mask_n%elts(id) >= ADJZONE) then
       dz =  (mean_m(id) + mass(id)) / porous_density (d, id, zlev)
       scalar_2d(id) = scalar_2d(id) + dz
    end if
  end subroutine cal_height

  subroutine cpt_or_restr_eta (dom, l)
    implicit none
    type(Domain) :: dom
    integer      :: l

    integer                     :: j, p_par, c, p_chd
    logical, dimension(N_CHDRN) :: restrict

    do j = 1, dom%lev(l)%length
       p_par = dom%lev(l)%elts(j)
       restrict = .false.
       do c = 1, N_CHDRN
          p_chd = dom%patch%elts(p_par+1)%children(c)
          if (p_chd > 0) restrict(c) = .true.
       end do
       do c = 1, N_CHDRN
          if (restrict(c)) then
             call apply_interscale_to_patch3 (eta_cpt_restr, dom, p_par, c, z_null, 0, 1)
          end if
       end do
    end do
  end subroutine cpt_or_restr_eta

  subroutine eta_cpt_restr (dom, p_chd, i_par, j_par, i_chd, j_chd, zlev, offs_par, dims_par, offs_chd, dims_chd)
    ! Compute or restrict eta for calculation of grad(eta)
    implicit none
    type(Domain)                   :: dom
    integer                        :: p_chd, i_par, j_par, i_chd, j_chd, zlev
    integer, dimension(N_BDRY+1)   :: offs_par, offs_chd
    integer, dimension(2,N_BDRY+1) :: dims_par, dims_chd

    integer :: id_par, id_chd

    id_chd = idx (i_chd, j_chd, offs_chd, dims_chd) + 1
    id_par = idx (i_par, j_par, offs_par, dims_par) + 1

    if (dom%mask_n%elts(id_par) >= RESTRCT) scalar(id_par) = scalar(id_chd)
  end subroutine eta_cpt_restr

  subroutine grad_eta
    ! Calculates grad eta (external pressure gradient due to free surface perturbation)
    implicit none
    integer :: d, j, k, l

    call update_bdry (sol(S_MASS,zlevels+1), NONE, 601)
    
    ! Calculate external pressure gradient
    sol(S_TEMP,zlevels+1) = sol(S_MASS,zlevels+1) ! copy eta to avoid modification by restriction
    do l = level_end, level_start, -1 
       do d = 1, size(grid)
          h_flux => horiz_flux(S_TEMP)%data(d)%elts
          scalar => sol(S_TEMP,zlevels+1)%data(d)%elts
          if (l < level_end) call cpt_or_restr_eta (grid(d), l) ! restrict eta if possible
          do j = 1, grid(d)%lev(l)%length
             call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=3)
          end do
          nullify (h_flux, scalar)
       end do
    end do
  end subroutine grad_eta

  subroutine vertical_velocity
    ! Computes vertical velocity at nodes, stored in  trend(S_MASS,k)
    use adapt_mod
    implicit none
    integer :: d, j, k, l, p

    call update_array_bdry (sol, NONE, 500)

    ! - Divergence of vertically integrated thickness flux, stored in trend(S_MASS,1)
    do l = level_end, level_start, -1
       do d = 1, size(grid)
          h_flux => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call step1 (q=sol, dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=4)
          end do
          if (l < level_end) then
             dscalar => trend(S_MASS,1)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)
             nullify (dscalar)
          end if
          nullify (h_flux)
       end do
       horiz_flux(S_MASS)%bdry_uptodate = .false.
       if (level_start /= level_end) call update_bdry (horiz_flux(S_MASS), l, 211)
       do d = 1, size(grid)
          dscalar =>    trend(S_MASS,1)%data(d)%elts
          h_flux  => horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       trend(S_MASS,1)%bdry_uptodate = .false.
       if (level_start /= level_end) call update_bdry (trend(S_MASS,1), l, 212)
    end do

    ! - Divergence of thickness flux at each vertical level, stored in exner_fun(1:zlevels)
    do k = 1, zlevels
       do l = level_end, level_start, -1
          do d = 1, size(grid)
             mass   => sol(S_MASS,k)%data(d)%elts
             velo   => sol(S_VELO,k)%data(d)%elts
             mean_m => sol_mean(S_MASS,k)%data(d)%elts
             h_flux => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call step1 (q=sol, dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=5)
             end do
             nullify (mass, mean_m, velo)
             if (l < level_end) then
                dscalar => exner_fun(k)%data(d)%elts
                call cpt_or_restr_flux (grid(d), l)
                nullify (dscalar)
             end if
             nullify (h_flux)
          end do
          horiz_flux(S_MASS)%bdry_uptodate = .false.
          if (level_start /= level_end) call update_bdry (horiz_flux(S_MASS), l, 211)
          do d = 1, size(grid)
             dscalar => exner_fun(k)%data(d)%elts
             h_flux  => horiz_flux(S_MASS)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (dscalar, h_flux)
          end do
          exner_fun(k)%bdry_uptodate = .false.
          if (level_start /= level_end) call update_bdry (exner_fun(k), l, 212)
       end do
    end do

    ! Integrate up to find vertical velocity, stored in trend(S_TEMP,1:zlevels)
    do l = level_end, level_start, -1
       do d = 1, size(grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_vertical_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (velo1, velo2)
       end do
    end do
    trend(S_TEMP,:)%bdry_uptodate = .false.
    call update_vector_bdry (trend(S_TEMP,:), NONE, 500)
  end subroutine vertical_velocity

  subroutine cal_vertical_velocity (dom, i, j, zlev, offs, dims)
    ! Vertical velocity
    ! (recall that we compute -div quantities)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                        :: d, id, id_i, k
    real(8)                        :: deta_dt, rho
    real(8), dimension (0:zlevels) :: w 
    
    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       d = dom%id + 1
       deta_dt = trend(S_MASS,1)%data(d)%elts(id_i)
       
       ! Vertical velocity at layer interfaces
       w(0) = 0.0_8; w(zlevels) = 0.0_8 ! impose zero vertical velocity at bottom and top
       do k = 1, zlevels-1
          w(k) = w(k-1) - exner_fun(k)%data(d)%elts(id_i)
       end do

       ! Interpolate to nodes and remove density
       do k = zlevels, 1, -1
          rho = porous_density (d, id_i, k)
          w(k) = interp (w(k-1), w(k)) / rho
       end do
       
       ! Compute vertical velocity relative to z coordinate
       do k = 1, zlevels
          trend(S_TEMP,k)%data(d)%elts(id_i) = w(k) + proj_vel_vertical ()
       end do
    end if
  contains
    real(8) function proj_vel_vertical ()
      ! Computes grad_zonal(z) * u_zonal + grad_merid(z) * u_merid at hexagon centres for vertical velocity computation.
      ! Uses Perot formula as also used for kinetic energy:
      ! u = sum ( u.edge_normal * hexagon_edge_length * (edge_midpoint-hexagon_center) ) / cell_area
      implicit none
      integer     :: idN, idE, idNE, idS, idSW, idW

      velo => sol(S_VELO,k)%data(d)%elts
      
      idE  = idx (i+1, j,   offs, dims)
      idNE = idx (i+1, j+1, offs, dims)
      idN  = idx (i,   j+1, offs, dims)
      idW  = idx (i-1, j,   offs, dims)
      idSW = idx (i-1, j-1, offs, dims)
      idS  = idx (i,   j-1, offs, dims)
    
      proj_vel_vertical =  &
           (vert_vel (i,j,i+1,j,EDGE*id +RT+1) + vert_vel (i+1,j+1,i,j,EDGE*id  +DG+1) + vert_vel (i,j,i,j+1,EDGE*id +UP+1) + &
           (vert_vel (i-1,j,i,j,EDGE*idW+RT+1) + vert_vel (i,j,i-1,j-1,EDGE*idSW+DG+1) + vert_vel (i,j-1,i,j,EDGE*idS+UP+1))) / 6
      
      nullify (velo)
    end function proj_vel_vertical

    real(8) function vert_vel (i1, j1, i2, j2, ide)
      implicit none
      integer :: i1, j1, i2, j2, ide

      real(8) :: dl, dz

      dz =  z_i (dom, i2, j2, k, offs, dims) - z_i (dom, i1, j1, k, offs, dims)
      dl = dom%len%elts(ide)

      vert_vel = dz / sqrt (dl**2 + dz**2) * velo(ide)
    end function vert_vel
  end subroutine cal_vertical_velocity
  
  subroutine cal_omega (dom, i, j, zlev, offs, dims)
    ! Velocity flux across interfaces
    ! (recall that we compute -div quantities)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                        :: d, id_i, k
    real(8)                        :: deta_dt, eta, rho, z_s
    real(8), dimension (0:zlevels) :: omega, z
    
    id_i = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       d = dom%id + 1

       if (sigma_z) then
          eta = sol(S_MASS,zlevels+1)%data(d)%elts(id_i)
          z_s = dom%topo%elts(id_i)
          z = z_coords (eta, z_s) ! set a_vert
       end if

       deta_dt = trend(S_MASS,1)%data(d)%elts(id_i)
       
       omega(0) = 0.0_8; omega(zlevels) = 0.0_8 ! impose zero flux at bottom and top
       do k = 1, zlevels-1
          omega(k) = omega(k-1) + (a_vert(k+1)-a_vert(k)) * deta_dt - exner_fun(k)%data(d)%elts(id_i)
       end do
       do k = zlevels, 1, -1
          rho = porous_density (d, id_i, k)
          omega(k) = interp (omega(k-1), omega(k)) / rho
       end do
       do k = 1, zlevels
          trend(S_TEMP,k)%data(d)%elts(id_i) = omega(k)
       end do
    end if
  end subroutine cal_omega
end module barotropic_2d_mod
