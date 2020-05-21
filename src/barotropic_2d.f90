module barotropic_2d_mod
  ! Files needed to solve barotropic free surface 
  use ops_mod
  use multi_level_mod

  implicit none
  real(8) :: dt1

  ! Parameters 
  integer, parameter :: iter = 25     ! number of iterations used in elliptic solver at coarsest scale (minimum suggested = 15)
  real(8), parameter :: w0   = 1.0_8  ! relaxation parameter for linear solver
  logical, parameter :: log  = .false. ! print out residual errors for elliptic solver
contains
  subroutine u_star (q)
    ! Explicit Euler step for intermediate velocity u_star
    ! remove external pressure gradient
    implicit none
    type(Float_Field), dimension(:,:), target :: q
    
    integer :: d, ibeg, iend, k, l

    ! External pressure gradient
    do l = level_end, level_start, -1
       call grad_eta (q(:,zlevels+1), horiz_flux(S_TEMP), l)
    end do
    
    do k = 1, zlevels
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = sol(S_VELO,k)%data(d)%length
          q(S_VELO,k)%data(d)%elts(ibeg:iend) = sol(S_VELO,k)%data(d)%elts(ibeg:iend) &
               + dt1 * (trend(S_VELO,k)%data(d)%elts(ibeg:iend) - horiz_flux(S_TEMP)%data(d)%elts(ibeg:iend))
       end do
    end do
  end subroutine u_star

  subroutine u_update (q)
    ! Explicit Euler velocity update with new external pressure gradient
    ! penalization is advanced using a backwards Euler scheme
    implicit none
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, ibeg, iend, k, l
    
    ! External pressure gradient
    do l = level_end, level_start, -1
       call grad_eta (q(:,zlevels+1), horiz_flux(S_TEMP), l)
    end do
    
    do k = 1, zlevels
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = sol(S_VELO,k)%data(d)%length
          q(S_VELO,k)%data(d)%elts(ibeg:iend) = (q(S_VELO,k)%data(d)%elts(ibeg:iend) &
               + dt1 * horiz_flux(S_TEMP)%data(d)%elts(ibeg:iend)) 
       end do
    end do
  end subroutine u_update

  subroutine scalar_update (q)
    ! Explicit Euler step for density
    implicit none
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, ibeg, iend, j, k, l, v

    do k = 1, zlevels
       do d = 1, size(grid)
          do v = scalars(1), scalars(2)
             ibeg = (1+2*(POSIT(v)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
             iend = sol(v,k)%data(d)%length
             q(v,k)%data(d)%elts(ibeg:iend) = sol(v,k)%data(d)%elts(ibeg:iend) + dt1 * trend(v,k)%data(d)%elts(ibeg:iend)
          end do
       end do
    end do
  end subroutine scalar_update

  subroutine eta_star (q)
    ! Explicit Euler step for intermediate free surface eta_star
    implicit none
    type(Float_Field), dimension(:,:), target :: q
    
    integer :: d, j, k, l

    ! Vertically integrated horizontal flux of intermediate horizontal velocity u_star
    call sum_vertical_flux (q, horiz_flux(S_MASS))
    
    do l = level_end, level_start, -1
       ! Calculate divergence of vertically integrated velocity flux, stored in trend(S_MASS,zlevels+1)
       do d = 1, size(grid)
          dscalar => trend(S_MASS,zlevels+1)%data(d)%elts
          h_flux  =>      horiz_flux(S_MASS)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_div, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, h_flux)
       end do
       trend(S_MASS,zlevels+1)%bdry_uptodate = .false.
       call update_bdry (trend(S_MASS,zlevels+1), l, 303)

       ! Restrict flux F if possible
       if (l < level_end) then
          do d = 1, size(grid)
             dscalar => trend(S_MASS,zlevels+1)%data(d)%elts
             h_flux  =>      horiz_flux(S_MASS)%data(d)%elts
             call cpt_or_restr_flux (grid(d), l)
             nullify (dscalar, h_flux)
          end do
       end if

       ! Euler step for intermediate free surface perturbation eta_star
       do d = 1, size(grid)
          dscalar => trend(S_MASS,zlevels+1)%data(d)%elts
          mass    =>   sol(S_MASS,zlevels+1)%data(d)%elts
          mass1   =>     q(S_MASS,zlevels+1)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (etastar_euler, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, mass, mass1)
       end do
       q(S_MASS,zlevels+1)%bdry_uptodate = .false.
       call update_bdry (q(S_MASS,zlevels+1), l, 304)
    end do
  end subroutine eta_star

  subroutine etastar_euler (dom, i, j, zlev, offs, dims)
    ! eta_star Euler step
    ! results are an edge flux
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1
        
    if (dom%mask_n%elts(id) >= ADJZONE) mass1(id) = mass(id) - dt1 * dscalar(id)
  end subroutine etastar_euler

  subroutine eta_update (q)
    ! Backwards Euler step for eta update
    use lin_solve_mod
    implicit none
    type(Float_Field), dimension(:,:), target :: q

    call rhs_elliptic (q)
    call multiscale (q(S_MASS,zlevels+1), trend(S_MASS,zlevels+1), elliptic_lo, elliptic_diag, iter, w0, log)
  end subroutine eta_update

  subroutine rhs_elliptic (q)
    ! Sets rhs of elliptic equation for free surface as trend(S_MASS,zlevels+1)
    implicit none
    type(Float_Field), dimension(:,:), target :: q
    
    integer :: d, j, k, l, p

    call update_array_bdry (q, NONE, 300)

    do l = level_end, level_start, -1
       ! Complete rhs of elliptic equation, eta^n - dt*div(F)
       do d = 1, size(grid)
          dscalar => trend(S_MASS,zlevels+1)%data(d)%elts
          mass    =>     q(S_MASS,zlevels+1)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_elliptic_rhs, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (dscalar, mass)
       end do
       trend(S_MASS,zlevels+1)%bdry_uptodate = .false.
       call update_bdry (trend(S_MASS,zlevels+1), l, 304)
    end do
  end subroutine rhs_elliptic

  subroutine cal_elliptic_rhs (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id) >= ADJZONE) then
       dscalar(id) = - mass(id)
    else
       dscalar(id) = 0.0_8
    end if
  end subroutine cal_elliptic_rhs

  function elliptic_lo (q, l)
    ! Calculates linear operator L(eta) for barotropic elliptic equation for free surface perturbation at scale l
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_lo, q

    integer :: d, j

    elliptic_lo = q

    ! Calculate external pressure gradient flux
    call external_pressure_gradient_flux (q, horiz_flux(S_MASS), l)

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

    ! Form complete linear operator  div(g(H+eta^n)grad(eta^(n+1)) - eta^(n+1)/dt^2
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
    call update_bdry (elliptic_lo, l, 101)
  end function elliptic_lo

  subroutine complete_elliptic_lo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id) >= ADJZONE) then
       dscalar(id) = dt**2*Laplacian_scalar(S_MASS)%data(dom%id+1)%elts(id) - mass(id)
    else
       dscalar(id) = 0.0_8
    end if
  end subroutine complete_elliptic_lo

  function elliptic_diag (q, l)
    ! Multiplies float array eta by inverse of diagonal part of barotropic elliptic equation linear operator l
    implicit none
    integer                   :: l
    type(Float_Field), target :: elliptic_diag, q

    integer :: d, j

    elliptic_diag = q
    
    do d = 1, size(grid)
       scalar => q%data(d)%elts
       diag   => elliptic_diag%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_elliptic_inv_diag, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (diag, scalar)
    end do
    call update_bdry (elliptic_diag, l, 100)
  end function elliptic_diag

  subroutine cal_elliptic_inv_diag (dom, i, j, zlev, offs, dims)
    ! Approximate inverse of diagonal part of barotropic elliptic equation linear operator
    ! note that eta^n is saved in exner_fun(2)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                 :: d, id, id_i, idE, idN, idNE, idS, idSW, idW
    real(8)                 :: f1, f2, Laplacian_diag
    real(8), dimension(1:6) :: pl
    real(8), dimension(0:6) :: phi
    
    id = idx (i, j, offs, dims)
    id_i = id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       d = dom%id+1
       
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims) 
       idS  = idx (i,   j-1, offs, dims)

       phi(0) = phi_node (d, id+1,   zlevels)
       phi(1) = phi_node (d, idE+1,  zlevels)
       phi(2) = phi_node (d, idNE+1, zlevels)
       phi(3) = phi_node (d, idN+1,  zlevels)
       phi(4) = phi_node (d, idW+1,  zlevels)
       phi(5) = phi_node (d, idSW+1, zlevels)
       phi(6) = phi_node (d, idS+1,  zlevels)

       pl(1) = dom%pedlen%elts(EDGE*id+RT+1)   / dom%len%elts(EDGE*id+RT+1)
       pl(2) = dom%pedlen%elts(EDGE*id+DG+1)   / dom%len%elts(EDGE*id+DG+1)
       pl(3) = dom%pedlen%elts(EDGE*id+UP+1)   / dom%len%elts(EDGE*id+UP+1)
       pl(4) = dom%pedlen%elts(EDGE*idW+RT+1)  / dom%len%elts(EDGE*idW+RT+1)
       pl(5) = dom%pedlen%elts(EDGE*idSW+DG+1) / dom%len%elts(EDGE*idSW+DG+1)
       pl(6) = dom%pedlen%elts(EDGE*idS+UP+1)  / dom%len%elts(EDGE*idS+UP+1)

       f1 = abs (dom%topo%elts(id_i)) * sum (pl)
            
       f2 = abs ( &
            dom%topo%elts(idE+1)*phi(1)*pl(1) + dom%topo%elts(idNE+1)*phi(2)*pl(2) + dom%topo%elts(idN+1)*phi(3)*pl(3) + &
            dom%topo%elts(idW+1)*phi(4)*pl(4) + dom%topo%elts(idSW+1)*phi(5)*pl(5) + dom%topo%elts(idS+1)*phi(6)*pl(6)) / phi(0)

       Laplacian_diag = - dom%areas%elts(id_i)%hex_inv * grav_accel * (f1 + f2)/2

       diag(id_i) = scalar(id_i) / (dt**2*Laplacian_diag - 1.0_8)
    else
       diag(id_i) = 0.0_8
    end if
  end subroutine cal_elliptic_inv_diag

  subroutine barotropic_correction (q)
    ! Update baroclinic variables mass and mass-weighted buoyancy with new free surface perturbation
    ! uses Bleck and Smith (J. Geophys. Res. 95, 3273–3285 1990) layer dilation method
    ! NOTE: individual layers no longer conserve mass (although total mass is conserved)
    implicit none
    type(Float_Field), dimension(:,:), target :: q
    
    integer :: d, k, p

    ! Sum mass perturbations
    call sum_vertical_mass (q(S_MASS,1:zlevels), exner_fun(1))

    do d = 1, size(grid)
       scalar    => q(S_MASS,zlevels+1)%data(d)%elts ! free surface perturbation
       scalar_2d => exner_fun(1)%data(d)%elts        ! sum of mass perturbations
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
    q(S_MASS,1:zlevels)%bdry_uptodate = .false.
    call update_vector_bdry (q(S_MASS,1:zlevels), NONE, 500)
  end subroutine barotropic_correction

  subroutine cal_barotropic_correction (dom, i, j, zlev, offs, dims)
    ! Correct baroclinic mass and buoyancy based on baroclinic estimate of free surface using layer dilation
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: full_mass, full_theta

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1 

    if (dom%mask_n%elts(id) >= ADJZONE) then
       full_mass  = mean_m(id) + mass(id)
       full_theta = (mean_t(id) + temp(id)) / full_mass ! buoyancy

       ! Correct mass
       mass(id) = ref_density*(scalar(id) - phi_node (d, id, zlev)*grid(d)%topo%elts(id))/scalar_2d(id) * full_mass &
            - mean_m(id)

       ! Correct mass-weighted buoyancy
       full_mass = mean_m(id) + mass(id)
       temp(id) = full_mass * full_theta - mean_t(id)
    end if
  end subroutine cal_barotropic_correction

  subroutine sum_vertical_mass (q, q_2d)
    ! Vertical sum of flux of q, returned in q_2d
    ! Assumes linearized free surface (i.e. remove free surface perturbation from sum)
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
             call apply_onescale_to_patch (cal_sum_vertical_mass, grid(d), p-1, z_null, 0, 1)
          end do
          nullify (mass, mean_m)
       end do
       nullify (scalar_2d)
    end do

    q_2d%bdry_uptodate = .false.
    call update_bdry (q_2d, NONE, 301)
  end subroutine sum_vertical_mass

  subroutine cal_sum_vertical_mass (dom, i, j, zlev, offs, dims)
    ! Vertical integration of edge quantity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    
    id = idx (i, j, offs, dims) + 1
    
    if (dom%mask_n%elts(id) >= ADJZONE) scalar_2d(id) = scalar_2d(id) + (mean_m(id) + mass(id))
  end subroutine cal_sum_vertical_mass

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
    ! (free surface eta is stored in pointer bernoulli)
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

  subroutine external_pressure_gradient_flux (q, dq, l) 
    ! Compute external pressure gradient flux term, depth^n * grad(eta^(n+1))
    implicit none
    integer                   :: l
    type(Float_Field), target :: dq, q

    integer :: d, j
    
    do d = 1, size(grid)
       h_flux => dq%data(d)%elts
       scalar => q%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=2)
       end do
       nullify (h_flux, scalar)
    end do

    dq%bdry_uptodate = .false.
    call update_bdry (dq, l, 301)
  end subroutine external_pressure_gradient_flux

  subroutine grad_eta (q, dq, l)
    ! Calculates grad eta (external pressure gradient due to free surface perturbation)
    implicit none
    integer                                 :: l
    type(Float_Field), dimension(:), target :: q
    type(Float_Field),               target :: dq
    
    integer :: d, j, k

    ! Copy eta to avoid modification by restriction
    q(S_TEMP) = q(S_MASS)

    do d = 1, size(grid)
       h_flux => dq%data(d)%elts
       scalar => q(S_TEMP)%data(d)%elts

       if (l < level_end) call cpt_or_restr_eta (grid(d), l) ! restrict grad(eta) if possible (modifies eta)

       ! Calculate external pressure gradient
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=3)
       end do
       nullify (h_flux, scalar)
    end do

    dq%bdry_uptodate = .false.
    call update_bdry (dq, l, 301)
  end subroutine grad_eta

  subroutine sum_vertical_flux (q, q_2d)
    ! Vertical sum of flux of q, returned in q_2d
    ! Assumes linearized free surface (i.e. remove free surface perturbation from sum)
    implicit none
    type(Float_Field),                 target :: q_2d
    type(Float_Field), dimension(:,:), target :: q

    integer :: d, k, p
    
    do d = 1, size(grid)
       scalar  => q(S_MASS,zlevels+1)%data(d)%elts
       velo_2d => q_2d%data(d)%elts

       q_2d%data(d)%elts = 0.0_8
       do k = 1, zlevels
          mass   => q(S_MASS,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          velo   => q(S_VELO,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (cal_sum_vertical_flux, grid(d), p-1, k, 0, 0)
          end do
          nullify (mass, mean_m, velo)
       end do
       nullify (scalar, velo_2d)
    end do

    q_2d%bdry_uptodate = .false.
    call update_bdry (q_2d, NONE, 301)
  end subroutine sum_vertical_flux

  subroutine cal_sum_vertical_flux (dom, i, j, zlev, offs, dims)
    ! Vertical integration of edge quantity
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: e, id, id_e, idE, idN, idNE
    real(8), dimension(1:EDGE) :: dz

    id = idx (i, j, offs, dims) 

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    
    ! Deduce (porous) layer heights from total mass
    dz(RT+1) = interp (mean_m(id+1)+mass(id+1), mean_m(idE+1) +mass(idE+1))  / ref_density 
    dz(DG+1) = interp (mean_m(id+1)+mass(id+1), mean_m(idNE+1)+mass(idNE+1)) / ref_density 
    dz(UP+1) = interp (mean_m(id+1)+mass(id+1), mean_m(idN+1) +mass(idN+1))  / ref_density

    ! Subtract free surface to linearize calculation
    if (zlev == zlevels) then
       dz(1) = dz(1) - interp (scalar(id+1), scalar(idE+1))  
       dz(2) = dz(2) - interp (scalar(id+1), scalar(idNE+1)) 
       dz(3) = dz(3) - interp (scalar(id+1), scalar(idN+1))  
    end if
 
    do e = 1, EDGE
       id_e = EDGE*id + e
       velo_2d(id_e) = velo_2d(id_e) + velo(id_e) * dom%pedlen%elts(id_e) * dz(e)
    end do
  end subroutine cal_sum_vertical_flux

  subroutine cal_div (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id

    id = idx (i, j, offs, dims) + 1

    if (dom%mask_n%elts(id) >= ADJZONE) then
       dscalar(id) = div (h_flux, dom, i, j, offs, dims)
    else
       dscalar(id) = 0.0_8
    end if
  end subroutine cal_div
end module barotropic_2d_mod

