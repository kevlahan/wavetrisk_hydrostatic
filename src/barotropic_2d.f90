module barotropic_2d_mod
  ! Files needed to solve barotropic free surface
  use ops_mod
  use multi_level_mod
  use lin_solve_mod
  use utils_mod
  implicit none
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
       iend = sol(S_VELO,1)%data(d)%length
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
    
    integer :: d, j, k, l

    call update_vector_bdry (q(S_MASS,1:zlevels+1), NONE)
    call update_vector_bdry (q(S_TEMP,1:zlevels),   NONE)

    do l = level_end, level_start, -1
       call total_height (q(S_MASS,1:zlevels), exner_fun(1), l) ! sum mass perturbations
       do d = 1, size(grid)
          scalar    => sol(S_MASS,zlevels+1)%data(d)%elts       ! free surface perturbation
          scalar_2d => exner_fun(1)%data(d)%elts                ! sum of mass perturbations
          do k = 1, zlevels
             mass   =>        q(S_MASS,k)%data(d)%elts
             temp   =>        q(S_TEMP,k)%data(d)%elts
             mean_m => sol_mean(S_MASS,k)%data(d)%elts
             mean_t => sol_mean(S_TEMP,k)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_barotropic_correction, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
             nullify (mass, temp, mean_m, mean_t)
          end do
          nullify (scalar, scalar_2d)
       end do
    end do
    q(S_MASS:S_TEMP,:)%bdry_uptodate = .false.
  end subroutine barotropic_correction

  subroutine cal_barotropic_correction (dom, i, j, zlev, offs, dims)
    ! Correct baroclinic mass and buoyancy based on baroclinic estimate of free surface using layer dilation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: d, id
    real(8) :: eta, full_mass, full_temp, mean_theta, theta, dz
    
    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1
    
    full_mass = mean_m(id) + mass(id)
    full_temp = mean_t(id) + temp(id)

    ! Full buoyancy
    theta = full_temp / full_mass 

    ! Free surface perturbation    
    eta = scalar(id) / phi_node (d, id, zlevels) 

    ! Correct mass perturbation
    mass(id) = (eta - topography%data(d)%elts(id)) / scalar_2d(id) * full_mass - mean_m(id)

    ! Update full mass
    full_mass = mean_m(id) + mass(id)

    ! Correct mass-weighted buoyancy
    temp(id) = full_mass * theta - mean_t(id)
  end subroutine cal_barotropic_correction

  subroutine eta_update
    ! Theta step for free surface update
    use lin_solve_mod
    implicit none

    ! RHS of elliptic equation
    call rhs_elliptic
    call equals_float_field (Laplacian_scalar(S_TEMP), sol(S_MASS,zlevels+1), AT_NODE) ! save old free surface height for elliptic operator
    
    ! Solve elliptic equation
    call elliptic_solver (sol(S_MASS,zlevels+1), sol(S_TEMP,zlevels+1), elliptic_lo, elliptic_lo_diag) 
  end subroutine eta_update

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
          sol(S_VELO,k)%data(d)%elts(ibeg:iend) = sol(S_VELO,k)%data(d)%elts(ibeg:iend) &
               - theta1 * dt * horiz_flux(S_TEMP)%data(d)%elts(ibeg:iend)
       end do
    end do
    sol(S_VELO,1:zlevels)%bdry_uptodate = .false.
  end subroutine u_update

  subroutine rhs_elliptic 
    ! Forms rhs of elliptic equation for free surface, -eta^* in q(S_TEMP_zlevels+1)
    ! trend(S_TEMP,zlevels+1) is flux divergence of vertically integrated velocity at previous time step
    ! (computed in RK routine)
    implicit none
    integer :: l

    call update_bdry (sol(S_MASS,zlevels+1), NONE)
    
    ! Flux divergence of vertically integrated velocity u_star, stored in trend(S_MASS,zlevels+1)
    call flux_divergence (sol, trend(S_MASS,zlevels+1))
    
    ! RHS of elliptic equation, -eta^*
    do l = level_end, level_start, -1
       call apply_onescale (cal_rhs_elliptic, l, z_null, 0, 1)
    end do
    sol(S_TEMP,zlevels+1)%bdry_uptodate = .false.
  contains
    subroutine cal_rhs_elliptic (dom, i, j, zlev, offs, dims)
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      sol(S_TEMP,zlevels+1)%data(d)%elts(id) = - sol(S_MASS,zlevels+1)%data(d)%elts(id) &
           + dt * (theta2 * trend(S_MASS,zlevels+1)%data(d)%elts(id) + (1d0 - theta2) * trend(S_TEMP,zlevels+1)%data(d)%elts(id)) &
           / ref_density
    end subroutine cal_rhs_elliptic
  end subroutine rhs_elliptic

    function elliptic_lo (q, l)
    ! Calculates linear operator L(eta) for barotropic elliptic equation for free surface perturbation at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_lo, q

    integer :: d, j

    call update_bdry (q, l)

    elliptic_lo = q; call zero_float_field (elliptic_lo, AT_NODE)
    
    ! Calculate external pressure gradient flux
    do d = 1, size(grid)
       h_flux =>       horiz_flux(S_MASS)%data(d)%elts
       mass   => Laplacian_scalar(S_TEMP)%data(d)%elts ! old free surface perturbation
       scalar =>                        q%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=2)
       end do
       nullify (h_flux, mass, scalar)
    end do
    horiz_flux(S_MASS)%bdry_uptodate = .false.
    call update_bdry (horiz_flux(S_MASS), l)

    ! Calculate divergence
    do d = 1, size(grid)
       dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
       h_flux  =>       horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, h_flux)
    end do
    Laplacian_scalar(S_MASS)%bdry_uptodate = .false.
    call update_bdry (Laplacian_scalar(S_MASS), l)

    ! Form complete linear operator 
    do d = 1, size(grid)
       dscalar => Laplacian_scalar(S_MASS)%data(d)%elts
       scalar  =>              elliptic_lo%data(d)%elts
       mass    =>                        q%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (complete_elliptic_lo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, mass, scalar)
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

    scalar(id) = theta1 * theta2 * dt**2 * dscalar(id) - mass(id)
  end subroutine complete_elliptic_lo

  function elliptic_lo_diag (q, l)
    ! Local approximation of diagonal of elliptic operator
    ! (Laplacian_scalar(S_TEMP) is the old free surface perturbation)
    ! ** using exact value of diagonal of Laplacian typically does NOT improve results **
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_lo_diag, q

    integer :: d, j

    elliptic_lo_diag = q; call zero_float_field (elliptic_lo_diag, AT_NODE)

    do d = 1, size(grid)
       dscalar => Laplacian_scalar(S_TEMP)%data(d)%elts
       scalar  =>         elliptic_lo_diag%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_elliptic_lo_diag, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, scalar)
    end do
    elliptic_lo_diag%bdry_uptodate = .false.
  end function elliptic_lo_diag

  subroutine cal_elliptic_lo_diag  (dom, i, j, zlev, offs, dims)
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer            :: d, id, id_i, idE, idNE, idN, idW, idSW, idS
    real(8)            :: depth, depth_e, Laplace_diag, wgt
    logical, parameter :: exact = .false.

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       depth = abs (topography%data(d)%elts(id_i)) + dscalar(id_i) / phi_node (d, id_i, zlevels)

       if (.not. exact) then ! average value 
          wgt = 2d0 * sqrt (3d0) * depth
       else ! true local value
          idE  = idx (i+1, j,   offs, dims) 
          idNE = idx (i+1, j+1, offs, dims) 
          idN  = idx (i,   j+1, offs, dims) 
          idW  = idx (i-1, j,   offs, dims) 
          idSW = idx (i-1, j-1, offs, dims) 
          idS  = idx (i,   j-1, offs, dims)

          depth_e = abs (topography%data(d)%elts(idE+1)) + dscalar(idE+1) / phi_node (d, idE+1, zlevels)
          wgt = dom%pedlen%elts(EDGE*id+RT+1) / dom%len%elts(EDGE*id+RT+1) * interp (depth_e, depth)

          depth_e = abs (topography%data(d)%elts(idNE+1)) + dscalar(idNE+1) / phi_node (d, idNE+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*id+DG+1) / dom%len%elts(EDGE*id+DG+1) * interp (depth_e, depth)

          depth_e = abs (topography%data(d)%elts(idN+1)) + dscalar(idN+1) / phi_node (d, idN+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*id+UP+1) / dom%len%elts(EDGE*id+UP+1) * interp (depth_e, depth)

          depth_e = abs (topography%data(d)%elts(idW+1)) + dscalar(idW+1) / phi_node (d, idW+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idW+RT+1) / dom%len%elts(EDGE*idW+RT+1) * interp (depth_e, depth)

          depth_e = abs (topography%data(d)%elts(idSW+1)) + dscalar(idSW+1) / phi_node (d, idSW+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idSW+DG+1) / dom%len%elts(EDGE*idSW+DG+1) * interp (depth_e, depth)

          depth_e = abs (topography%data(d)%elts(idS+1)) + dscalar(idS+1) / phi_node (d, idS+1, zlevels)
          wgt = wgt + dom%pedlen%elts(EDGE*idS+UP+1) / dom%len%elts(EDGE*idS+UP+1) * interp (depth_e, depth)
       end if

       Laplace_diag = - grav_accel * wgt * dt**2 * dom%areas%elts(id_i)%hex_inv
       scalar(id_i) = theta1 * theta2 * Laplace_diag - 1d0
    end if
  end subroutine cal_elliptic_lo_diag

  subroutine flux_divergence (q, div_flux)
    ! Returns flux divergence of vertical integrated velocity in divF using solution q, stored in div_flux
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q
    type(Float_Field),                                 target :: div_flux

    integer :: d, j, l

    call update_vector_bdry (q(S_MASS,1:zlevels), NONE)
    call update_vector_bdry (q(S_VELO,1:zlevels), NONE)

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
       call update_bdry (horiz_flux(S_MASS), l)

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
       call update_bdry (div_flux, l)
    end do
  end subroutine flux_divergence

  subroutine total_height (q, q_2d, l)
    ! Total height q_2d computed from pseudo-densities q
    implicit none
    integer                                 :: l
    type(Float_Field),               target :: q_2d
    type(Float_Field), dimension(:), target :: q

    integer :: d, j, k
    
    do d = 1, size(grid)
       scalar_2d => q_2d%data(d)%elts
       
       q_2d%data(d)%elts = 0d0
       do k = 1, zlevels
          mass   =>               q(k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_height, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
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
    real(8) :: dz, full_mass
    
    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    full_mass = mean_m(id) + mass(id)
    
    dz = full_mass / porous_density (d, id, zlev)
    
    scalar_2d(id) = scalar_2d(id) + dz
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
    integer :: d, j, l

    call update_bdry (sol(S_MASS,zlevels+1), NONE)

    call equals_float_field (sol(S_TEMP,zlevels+1), sol(S_MASS,zlevels+1), AT_NODE)

    ! Calculate external pressure gradient
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
end module barotropic_2d_mod
