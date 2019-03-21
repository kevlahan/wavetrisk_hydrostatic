module remap_mod
  use time_integr_mod
  use wavelet_mod
  implicit none
  abstract interface
     subroutine interpolation (N, var_new, z_new, var_old, z_old, neumann)
       implicit none
       integer                 :: N
       real(8), dimension(1:N) :: var_new, var_old
       real(8), dimension(0:N) :: z_new, z_old
       logical                 :: neumann
     end subroutine interpolation
  end interface
  procedure (interpolation), pointer :: interp_scalar => null ()
  procedure (interpolation), pointer :: interp_velo   => null ()
contains
  subroutine remap_vertical_coordinates
    ! Remap the Lagrangian layers to initial vertical grid given a_vert and b_vert vertical coordinate parameters 
    ! Conserves mass, potential temperature and velocity divergence
    ! remap0 is too diffusive; remap1, remap2W are very stable and remap2PPM, remap2S, remap4 are less stable.
    implicit none
    integer            :: l
    logical, parameter :: standard = .true. ! .false. uses Lin (2004) scheme which interpolates total energy (otherwise interpolate potential temperature)

    ! Choose interpolation method:
    ! [these methods are modified from routines provided by Alexander Shchepetkin (IGPP, UCLA)]
    !
    ! remap0    = piecewise-constant reconstruction (checks that grid is not changing too fast)
    ! remap1    = piecewise-linear reconstruction with van Leer limiter to ensure monotonicity
    ! remap2PPM = reconstruction by PPM code of Colella and Woodward (1984) 
    ! remap2S   = basic parabolic spline reconstruction
    ! remap2W   = parabolic WENO reconstruction
    ! remap4    = parabolic WENO reconstruction enhanced by quartic power-law reconciliation step
    !                  (ensures continuity of both value and first derivative at each interface)
    interp_scalar => remap2W
    interp_velo   => remap2W
    
    if (standard) then ! Standard (remap potential temperature)
       do l = level_start, level_end
          call apply_onescale (remap_scalars, l, z_null, 0, 1)
          call apply_onescale (remap_velo,    l, z_null, 0, 0)
       end do
    else ! Lin (2004) (remap total energy)
       do l = level_start, level_end
          call apply_onescale (remap_total_energy, l, z_null, 0, 1)
          call apply_onescale (remap_velo,         l, z_null, 0, 0)
          call update_vector_bdry (sol(S_VELO,:), l)
          call apply_onescale (recover_theta, l, z_null, 0, 1)
       end do
    end if

    ! Wavelet transform and interpolate back onto adapted grid
    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine remap_vertical_coordinates

  subroutine remap_scalars (dom, i, j, z_null, offs, dims)
    ! Remap mass-weighted potential temperature
    ! (potential temperature is remapped and then multiplied by new mass)
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                        :: d, id_i, k
    real(8)                        :: p_s
    real(8), dimension (1:zlevels) :: theta_new, theta_old 
    real(8), dimension (0:zlevels) :: p_new, p_old

    d    = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    ! Save old mass
    do k = 1, zlevels
       trend(S_MASS,k)%data(d)%elts(id_i) = sol(S_MASS,k)%data(d)%elts(id_i)
    end do
    
    call find_coordinates (p_new, p_old, d, id_i, p_s)

    do k = 1, zlevels
       theta_old(zlevels-k+1) = sol(S_TEMP,k)%data(d)%elts(id_i) / sol(S_MASS,k)%data(d)%elts(id_i)
       sol(S_MASS,k)%data(d)%elts(id_i) = a_vert_mass(k) + b_vert_mass(k) * p_s/grav_accel ! New mass
    end do
  
    call interp_scalar (zlevels, theta_new, p_new, theta_old, p_old, .false.)

    ! Mass-weighted potential temperature
    do k = 1, zlevels
       sol(S_TEMP,k)%data(d)%elts(id_i) = theta_new(zlevels-k+1) * sol(S_MASS,k)%data(d)%elts(id_i)
    end do
  end subroutine remap_scalars

  subroutine remap_velo (dom, i, j, z_null, offs, dims)
    ! Remap velocity
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                        :: d, e, id, id_i, k
    real(8)                        :: p_s
    real(8), dimension (1:zlevels) :: flux_new, flux_old 
    real(8), dimension (0:zlevels) :: p_new, p_old

    d    = dom%id + 1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    call find_coordinates (p_new, p_old, d, id_i, p_s)

    do e = 1, EDGE
       do k = 1, zlevels
          flux_old(zlevels-k+1) = sol(S_VELO,k)%data(d)%elts(EDGE*id+e)
       end do

       call interp_velo (zlevels, flux_new, p_new, flux_old, p_old, .false.)

       do k = 1, zlevels
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = flux_new(zlevels-k+1)
       end do
    end do
  end subroutine remap_velo

  subroutine remap_momentum (dom, i, j, z_null, offs, dims)
    ! Remap momentum and then recover remapped velocity
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                        :: d, e, id, id_i, k
    integer, dimension(1:EDGE)     :: id_r
    real(8)                        :: mass_e, p_s
    real(8), dimension (1:zlevels) :: flux_new, flux_old 
    real(8), dimension (0:zlevels) :: p_new, p_old

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    id_r(RT+1) = idx (i+1, j,   offs, dims) + 1
    id_r(DG+1) = idx (i+1, j+1, offs, dims) + 1
    id_r(UP+1) = idx (i,   j+1, offs, dims) + 1

    call find_coordinates (p_new, p_old, d, id_i, p_s)

    do e = 1, EDGE
       do k = 1, zlevels
          mass_e = trend(S_MASS,k)%data(d)%elts(id_i) + trend(S_MASS,k)%data(d)%elts(id_r(e))
          flux_old(zlevels-k+1) = sol(S_VELO,k)%data(d)%elts(EDGE*id+e) * mass_e
       end do

       call interp_velo (zlevels, flux_new, p_new, flux_old, p_old, .false.)

       do k = 1, zlevels
          mass_e = sol(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_r(e))
          sol(S_VELO,k)%data(d)%elts(EDGE*id+e) = flux_new(zlevels-k+1) / mass_e
       end do
    end do
  end subroutine remap_momentum

  subroutine remap_total_energy (dom, i, j, z_null, offs, dims)
    ! Remap total energy
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                        :: d, id, id_i, k, kb
    real(8)                        :: p_s, phi_lower, phi_upper, theta, T_mean
    real(8), dimension (1:zlevels) :: energy_new, energy_old
    real(8), dimension (0:zlevels) :: p_new, p_old

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    ! Save old mass
    do k = 1, zlevels
       trend(S_MASS,k)%data(d)%elts(id_i) = sol(S_MASS,k)%data(d)%elts(id_i)
    end do

    call find_coordinates (p_new, p_old, d, id_i, p_s)

    ! Integrate up
    phi_lower = surf_geopot (dom%node%elts(id_i))
    do k = zlevels, 1, -1
       kb = zlevels-k+1
       
       theta = sol(S_TEMP,kb)%data(d)%elts(id_i) / sol(S_MASS,kb)%data(d)%elts(id_i)
       
       T_mean = theta * (p_old(k)**kappa - p_old(k-1)**kappa) / log (p_old(k)/p_old(k-1)) / kappa 

       phi_upper = phi_lower +  R_d * T_mean * log (p_old(k)/p_old(k-1))
       
       velo => sol(S_VELO,kb)%data(d)%elts
       energy_old(k) = c_p*T_mean + (p_old(k)*phi_lower-p_old(k-1)*phi_upper) / (p_old(k)-p_old(k-1)) + ke (dom, i, j, offs, dims)
       nullify (velo)
       phi_lower = phi_upper
    end do
    
    call interp_scalar (zlevels, energy_new, p_new, energy_old, p_old, .false.)
    
    do k = 1, zlevels
       exner_fun(k)%data(d)%elts(id_i) = energy_new(k) ! Store new total energy in exner function
       sol(S_MASS,k)%data(d)%elts(id_i) = a_vert_mass(k) + b_vert_mass(k) * p_s/grav_accel ! New mass
    end do
  end subroutine remap_total_energy
  
  subroutine recover_Theta (dom, i, j, z_null, offs, dims)
    ! Recover mass-weighted potential temperature from remapped total energy
    type (Domain)                   :: dom
    integer                         :: i, j, z_null
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                        :: d, id_i, k, kb
    real(8)                        :: phi_lower, p_s, T_mean, theta
    real(8), dimension (0:zlevels) :: p_new, p_old

    d    = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    call find_coordinates (p_new, p_old, d, id_i, p_s)

    ! Integrate up from surface
    phi_lower = surf_geopot (dom%node%elts(id_i))
    do k = zlevels, 1, -1
       kb = zlevels-k+1
       
       velo => sol(S_VELO,kb)%data(d)%elts
       T_mean = (exner_fun(k)%data(d)%elts(id_i) - ke (dom, i, j, offs, dims) - phi_lower) &
            / (1.0_8 - kappa * p_new(k-1) * log (p_new(k)/p_new(k-1)) / (p_new(k)-p_new(k-1))) / c_p
       nullify (velo)

       theta = kappa * log (p_new(k)/p_new(k-1)) / (p_new(k)**kappa - p_new(k-1)**kappa) * T_mean
       
       sol(S_TEMP,kb)%data(d)%elts(id_i) = sol(S_MASS,kb)%data(d)%elts(id_i) * theta          ! New temp

       phi_lower = phi_lower +  R_d * T_mean * log (p_new(k)/p_new(k-1))
    end do
  end subroutine recover_Theta

  real(8) function ke (dom, i, j, offs, dims)
    type (Domain)                   :: dom
    integer                         :: i, j
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    
    ! Kinetic energy
    integer :: id, idS, idSW, idW
    real(8) :: u_prim_RT,   u_prim_UP,   u_prim_DG,    u_dual_RT,   u_dual_UP,   u_dual_DG
    real(8) :: u_prim_RT_W, u_prim_UP_S, u_prim_DG_SW, u_dual_RT_W, u_dual_UP_S, u_dual_DG_SW

    id    = idx (i,   j,   offs, dims)
    idW   = idx (i-1, j,   offs, dims)
    idS   = idx (i,   j-1, offs, dims)
    idSW  = idx (i-1, j-1, offs, dims) 

    u_prim_RT    = velo(EDGE*id  +RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_dual_RT    = velo(EDGE*id  +RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id+UP+1)
    u_prim_DG    = velo(EDGE*id  +DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_dual_DG    = velo(EDGE*id  +DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP    = velo(EDGE*id  +UP+1) * dom%pedlen%elts(EDGE*id+UP+1)
    u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS+UP+1)
    u_dual_UP_S  = velo(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS+UP+1)
    u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW+RT+1)
    u_dual_RT_W  = velo(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW+RT+1)
    u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)

    ke = (u_prim_UP*u_dual_UP + u_prim_DG*u_dual_DG + u_prim_RT*u_dual_RT + &
         u_prim_UP_S*u_dual_UP_S + u_prim_DG_SW*u_dual_DG_SW + u_prim_RT_W*u_dual_RT_W) * dom%areas%elts(id+1)%hex_inv/4
  end function ke

  subroutine find_coordinates (p_new, p_old, d, id_i, p_s)
    ! Calculates old and new pressure-based z coordinates
    implicit none
    integer                       :: d, id_i
    real(8)                       :: p_s
    real(8), dimension(0:zlevels) :: p_new, p_old

    integer :: k

    p_old(0) = p_top
    do k = 1, zlevels
       p_old(k) = p_old(k-1) + grav_accel * trend(S_MASS,zlevels-k+1)%data(d)%elts(id_i)
    end do
    p_s = p_old(zlevels)
    
    p_new = a_vert(zlevels+1:1:-1) + b_vert(zlevels+1:1:-1) * p_s
  end subroutine find_coordinates

  subroutine remap0 (N, var_new, z_new, var_old, z_old, neumann)
    !
    ! The simplest remapping procedure which assumes that the distribution of var_old(z) is piecewise constant 
    ! in each grid box, (i.e. similar to donor-cell first-order upstream advection).
    !
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new, z_old
    logical                 :: neumann
    
    integer                 :: k
    real(8)                 :: dz
    real(8), parameter      :: Zero=0.0_8
    real(8), dimension(0:N) :: FC
    
    do k = 1, N-1
       dz = z_new(k) - z_old(k)
       FC(k) = min (dz, Zero) * var_old(k) + max (dz, Zero) * var_old(k+1)
    end do
    FC(0) = 0.0_8
    FC(N) = 0.0_8
    
    do k = 1, N
       var_new(k) = ((z_old(k)-z_old(k-1))*var_old(k) + FC(k)-FC(k-1)) / (z_new(k)-z_new(k-1))
    end do
  end subroutine remap0

  subroutine remap1 (N, var_new, z_new, var_old, z_old, neumann)
    !
    ! Remapping procedure using piecewise-linear reconstruction: 
    !---------- --------- ----- --------- ------ ---------------
    ! MINMOD version: The distribution in each grid-box is assumed to be
    !        linear with the slope equal to the smaller of two elementary
    !        slopes computed by linear differencing (set to 0, if the two
    !        slopes are of different signs). This method is commonly known
    !        as that of van Leer. It guarantees monotonicity.
    !                                                              2*X*Y
    ! default version: minmod(X,Y) is replaced with harmonic mean -------
    !                                                              X + Y
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new, z_old
    logical                 :: neumann

    integer                 :: k, iter
    real(8)                 :: cff, cff1, dh, dL, dR, dz
    real(8), parameter      :: Zero=0.0_8, Half=0.5_8
    real(8), dimension(0:N) :: aL, aR, FC
    real(8), dimension(1:N) :: Hz
    logical, parameter      :: ENHANCE = .true.   

    do k = 1, N
       Hz(k) = z_old(k) - z_old(k-1)
    end do

    do k = N-1, 1, -1
       FC(k) = (var_old(k+1) - var_old(k)) / (Hz(k+1) + Hz(k))
    end do

    if (NEUMANN) then
       FC(0) = Zero
       FC(N) = Zero
    else
       FC(0) = FC(1)
       FC(N) = FC(N-1)
    end if

    do k = 1, N
       if (FC(k)*FC(k-1) < Zero) then
          cff = Zero
       else
          cff = 2*FC(k)*FC(k-1) / (FC(k) + FC(k-1))
       end if
       aR(k) = var_old(k) + cff*Hz(k)
       aL(k) = var_old(k) - cff*Hz(k)
    end do !--> discard FC
    
    if (ENHANCE) then ! reconcile side values at each interface, then use minmod limiter again to maintain monotonicity
       do iter = 1, 10                        
          do k = 1, N-1                        
             FC(k) = Half * (aR(k) + aL(k+1))  
          end do
          FC(N) = aR(N)
          FC(0) = aL(1)
          do k = 1, N
             dR = FC(k) - var_old(k)
             dL = var_old(k) - FC(k-1)
             if (dR*dL < Zero) then
                cff = Zero
             else if (abs(dR) > abs(dL)) then
                cff = dL
             else
                cff = dR
             end if
             aR(k) = var_old(k) + cff
             aL(k) = var_old(k) - cff
          end do
       end do
    end if

    ! Remapping: compute
    do k = 1, N-1 ! finite volume fluxes
       dz = z_new(k)-z_old(k)
       if (dz > Zero) then
          cff = aL(k+1)
          cff1 = aR(k+1) - aL(k+1)
          dh = Hz(k+1)
       else
          cff  = aR(k)
          cff1 = aR(k) - aL(k)
          dh = Hz(k)
       end if
       FC(k) = dz * (cff + Half*cff1*dz/dh)
    end do
    FC(0) = Zero
    FC(N) = Zero

    do k = 1, N
       var_new(k) = (Hz(k)*var_old(k) + FC(k)-FC(k-1)) / (z_new(k)-z_new(k-1))
    end do
  end subroutine remap1

  subroutine remap2PPM (N, var_new, z_new, var_old, z_old, neumann)
    !
    ! Reconstruction by PPM code of Colella--Woodward, 1984.
    !
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new, z_old
    logical                 :: neumann

    integer                 :: k, k1, k2
    real(8)                 :: alpha, cff, cffL, cffR, dL, dR, dz
    real(8), parameter      :: Zero=0.0_8, Half=0.5_8, One=1.0_8, ThreeHalfth=1.5_8, Two=2.0_8, Three=3.0_8
    real(8), dimension(1:N) :: Hz
    real(8), dimension(0:N) :: aL, aR, CF, FC, FC1
    logical, parameter      :: LIMIT_INTERIOR = .false.
    logical, parameter      :: LIMIT_SLOPES   = .true. 

    do k = 1, N
       Hz(k) = z_old(k) - z_old(k-1)
    end do

    do k = 1, N-1
       cff = One / (Hz(k) + Hz(k+1))
       CF(k) = cff * (var_old(k+1)*Hz(k) + var_old(k)*Hz(k+1))
       FC(k) = cff * (var_old(k+1) - var_old(k))
    enddo
    
    do k = 2, N-1
       cff = Hz(k) * ((Two*Hz(k-1) + Hz(k))*FC(k) + (Two*Hz(k+1) + Hz(k))*FC(k-1)) / (Hz(k-1) + Hz(k) + Hz(k+1))
       if (LIMIT_SLOPES) then
          cffR = Two * (var_old(k+1) - var_old(k))
          cffL = Two * (var_old(k) - var_old(k-1))
          if (cffR*cffL > Zero) then
             if (abs(cffL) < abs(cff)) cff = cffL
             if (abs(cffR) < abs(cff)) cff = cffR
          else
             cff = Zero
          end if
       end if
       FC1(k) = cff
    end do
    FC1(N) = Hz(N) * FC(N-1)
    FC1(1) = Hz(1) * FC(1)

    do k = 1, N-1
       k1 = max (k-1, 1)
       k2 = min (k+2, N)
       cff = Hz(k1) + Hz(k) + Hz(k+1) + Hz(k2)
       cffL = (Hz(k1)  + Hz(k)) / (cff*(Two*Hz(k) + Hz(k+1)))
       cffR = (Hz(k+1) + Hz(k2)) / (cff*(Hz(k) + Two*Hz(k+1)))

       aR(k) = CF(k) + Two*Hz(k)*Hz(k+1) * (cffL-cffR)*FC(k) - cffL*Hz(k)*FC1(k+1) + cffR*Hz(k+1)*FC1(k)
       aL(k+1) = aR(k)    
    end do
    
    if (NEUMANN) then
       aR(N) = ThreeHalfth*var_old(N) - Half*aL(N)
       aL(1) = ThreeHalfth*var_old(1) - Half*aR(1)
    else
       aR(N) = Two*var_old(N) - aL(N)
       aL(1) = Two*var_old(1) - aR(1)
    end if
    
    if (LIMIT_INTERIOR) then
       do k = 1, N
          dR = aR(k) - var_old(k)
          dL = var_old(k) - aL(k)
          if (dR*dL < Zero) then
             dR = Zero
             dL = Zero
          end if
          if (abs(dR) > Two*abs(dL)) dR = Two * dL
          if (abs(dL) > Two*abs(dR)) dL = Two * dR
          aR(k) = var_old(k) + dR
          aL(k) = var_old(k) - dL
       end do
    end if
    
    !
    ! Remapping step: This operation consists essentially of three
    !---------------- stages: (1) whithin each grid box compute averaged 
    ! slope (stored as CF) and curvature (stored as FC); then (2) compute
    ! interfacial fluxes FC; and (3) apply these fluxes to complete
    ! remapping step.
    !
    do k = 1, N
       CF(k) = Half * (aR(k) - aL(k))
       FC(k) = Half * (aR(k) + aL(k)) - var_old(k)
    end do

    do k = 1, N-1, +1       !<-- irreversible 
       dz = z_new(k) - z_old(k)
       if (dz > Zero) then
          alpha = Hz(k+1)   
          cff   = aL(k+1)
          cffR  = CF(k+1)
          cffL  = FC(k+1)
       else
          alpha = -Hz(k)
          cff   =  aR(k)
          cffR  = -CF(k)
          cffL  =  FC(k)
       end if
       alpha = dz / alpha
       FC(k) = dz * (cff + alpha*(cffR - cffL*(Three - Two*alpha)))
    end do
    FC(0) = Zero
    FC(N) = Zero

    do k = 1, N
       var_new(k) = (Hz(k)*var_old(k) + FC(k)-FC(k-1)) / (z_new(k) - z_new(k-1))
    end do
  end subroutine remap2PPM

  subroutine remap2S (N, var_new, z_new, var_old, z_old, dummy)
    !
    ! Basic parabolic spline reconstruction
    !------ --------- ------ --------------
    !
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new,  z_old
    logical                 :: dummy

    integer                 :: k
    real(8)                 :: alpha, cff, cff1, dz
    real(8), parameter      :: Half = 0.5_8, ThreeHalfth=1.5_8, Zero = 0.0_8, One = 1.0_8, Two = 2.0_8, Three = 3.0_8
    real(8), dimension(1:N) :: Hz
    real(8), dimension(0:N) :: dL, dR, FC, r
    character(255)          :: bc = "PARABOLIC_CONTINUATION" ! options are 'NEUMANN', 'LINEAR_CONTINUATION', 'PARABOLIC_CONTINUATION'

    do k = 1, N
       Hz(k) = z_old(k) - z_old(k-1)
    end do

    select case (bc)
    case ("PARABOLIC_CONTINUATION")
       cff = Hz(1) / Hz(2)
       FC(1) = One + cff
       r(0) = Two*var_old(1) + cff*( var_old(1) + cff*var_old(2) )/FC(1)
    case ('LINEAR_CONTINUATION') 
       FC(1) = One
       r(0) = Two * var_old(1)
    case ('NEUMANN')
       FC(1) = Half
       r(0) = ThreeHalfth * var_old(1)
    end select

    do k = 1, N-1, +1
       cff = One / (Two*Hz(k) + Hz(k+1)*(Two-FC(k)))
       FC(k+1) = cff * Hz(k)
       r(k) = cff * (Three*(var_old(k+1)*Hz(k) + var_old(k)*Hz(k+1)) - Hz(k+1)*r(k-1))
    end do
    
    select case (bc)
    case ('PARABOLIC_CONTINUATION')
       cff = Hz(N) / Hz(N-1)
       cff1 = One + cff
       r(N) = (cff*(var_old(N) + cff*var_old(N-1)) + cff1*(Two*var_old(N) - cff1*r(N-1)))/(cff1*(One - cff1*FC(N)))
    case ('LINEAR_CONTINUATION') 
       r(N) = (Two*var_old(N) - r(N-1)) / (One - FC(N))
    case ('NEUMANN')
       r(N) = (Three*var_old(N) - r(N-1)) / (Two - FC(N))
    end select

    do k = N-1, 0, -1
       r(k) = r(k) - FC(k+1)*r(k+1)
    end do
    !
    ! Remapping step: This operation consists essentially of three
    !---------------- stages: (1) whithin each grid box compute averaged 
    ! slope (stored as dR) and curvature (stored as dL); then (2) compute
    ! interfacial fluxes FC; and (3) apply these fluxes to complete
    ! remapping step.
    !
    do k = 1, N
       dR(k) = Half*(r(k) - r(k-1))
       dL(k) = Half*(r(k) + r(k-1)) - var_old(k)
    end do

    do k = 1, N-1
       dz = z_new(k) - z_old(k)
       if (dz > Zero) then
          alpha = Hz(k+1)   
          cff = dR(k+1)
          cff1 = dL(k+1)
       else
          alpha = Hz(k)
          cff  =  dR(k)
          cff1 = -dL(k)
       endif
       alpha = dz / alpha
       FC(k) = dz * (r(k) + alpha*(cff - cff1*(Three - Two*abs(alpha))))
    end do
    FC(0) = 0.0_8
    FC(N) = 0.0_8

    do k = 1, N
       var_new(k) = (Hz(k)*var_old(k) + FC(k)-FC(k-1)) / (z_new(k) - z_new(k-1))
    end do
  end subroutine remap2S

  subroutine remap2W (N, var_new, z_new, var_old, z_old, neumann)
    !
    ! Parabolic WENO reconstruction: The second and third loops below
    !---------- ---- --------------- compute left and right side limits
    ! aL, aR of the field var_old assuming monotonized parabolic distributions
    ! within each grid box. Also computed are dL, dR, which are then used
    ! as a measure of quadratic variation during sabsequent WENO
    ! reconciliation of side limits.
    !
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new, z_old
    logical                 :: neumann

    integer                 :: k
    real(8)                 :: alpha, cff, cffL, cffR, deltaL, deltaR, dz
    real(8), parameter      :: Zero=0.0_8, Half=0.5_8, One=1.0_8, ThreeHalfth=1.5_8, Two=2.0_8, Three=3.0_8, eps = 1.0d-8
    real(8), dimension(1:N) :: Hz
    real(8), dimension(0:N) :: aL, aR, dL, dR, FC, r
    logical, parameter      :: LIMIT_INTERIOR = .false.

    do k = 1, N
       Hz(k) = z_old(k) - z_old(k-1)
    end do

    do k = N-1, 1, -1
       FC(k) = (var_old(k+1) - var_old(k)) / (Hz(k+1) + Hz(k))
    end do

    do k = 2, N-1
       deltaR = Hz(k) * FC(k)
       deltaL = Hz(k) * FC(k-1)

       if (deltaR*deltaL < Zero) then
          deltaR = Zero
          deltaL = Zero
       end if
       cff = Hz(k-1) + Two*Hz(k) + Hz(k+1)
       cffR = cff * FC(k)
       cffL = cff * FC(k-1)
       if (abs(deltaR) > abs(cffL)) deltaR = cffL
       if (abs(deltaL) > abs(cffR)) deltaL = cffR

       cff = (deltaR - deltaL) / (Hz(k-1) + Hz(k) + Hz(k+1))
       deltaR = deltaR - cff*Hz(k+1)
       deltaL = deltaL + cff*Hz(k-1)

       aR(k) = var_old(k) + deltaR
       aL(k) = var_old(k) - deltaL

       dR(k) = (Two*deltaR - deltaL)**2
       dL(k) = (Two*deltaL - deltaR)**2
    end do
    
    if (LIMIT_INTERIOR) then
       aR(N) = var_old(N)     ! Boundary conditions for strictly monotonic option: The only way to
       aL(N) = var_old(N)     ! avoid extrapolation toward the
       dR(N) = Zero           ! boundary is to assume that field 
       dL(N) = Zero           ! is simply constant within topmost 
       aR(1) = var_old(1)     ! and bottommost grid boxes. Note 
       aL(1) = var_old(1)     ! that even for NEUMANN boundary 
       dR(1) = Zero           ! conditions, the extrapolated 
       dL(1) = Zero           ! values aR(N) and aL(0) exceed corresponding grid box values.
    else 
       aL(N) = aR(N-1)
       if (NEUMANN) then
          aR(N) = ThreeHalfth*var_old(N) - Half*aL(N)
       else
          aR(N) = Two*var_old(N) - aL(N)
       end if
       dR(N) = (Two*aR(N) + aL(N) - Three*var_old(N))**2
       dL(N) = (Three*var_old(N) - Two*aL(N) - aR(N))**2

       aR(1)=aL(2)
       if (NEUMANN) then
          aL(1) = ThreeHalfth*var_old(1) - Half*aR(1)
       else
          aL(1) = Two*var_old(1) - aR(1)
       end if
       dR(1) = (Two*aR(1) + aL(1) - Three*var_old(1))**2
       dL(1) = (Three*var_old(1) - Two*aL(1) - aR(1))**2
    end if

    ! Reconcile interfacial values aR, aL using WENO 
    do k = 1, N-1                        
       deltaL = max (dL(k),   eps)
       deltaR = max (dR(k+1), eps)
       r(k) = (deltaR*aR(k) + deltaL*aL(k+1)) / (deltaR + deltaL)
    end do

    if (NEUMANN) then
       r(N) = ThreeHalfth*var_old(N) - Half*r(N-1)
       r(0) = ThreeHalfth*var_old(1) - Half*r(1)
    else
       r(N) = Two*var_old(N) - r(N-1)
       r(0) = Two*var_old(1) - r(1)
    end if
    !
    ! Remapping step: This operation consists essentially of three
    !---------------- stages: (1) within each grid box compute averaged 
    ! slope (stored as dR) and curvature (stored as dL); then (2) compute
    ! interfacial fluxes FC; and (3) apply these fluxes to complete remapping step.
    !
    do k = 1, N
       if (LIMIT_INTERIOR) then               ! Constrain parabolic
          deltaR = r(k) - var_old(k)          ! segment monotonicity
          deltaL = var_old(k) - r(k-1)        ! like in PPM 
          cffR = Two * deltaR
          cffL = Two * deltaL
          if (deltaR*deltaL < Zero) then
             deltaR = Zero
             deltaL = Zero
          else if (abs(deltaR) > abs(cffL)) then
             deltaR = cffL
          else if (abs(deltaL) > abs(cffR)) then
             deltaL = cffR
          endif
          aR(k) = var_old(k) + deltaR
          aL(k) = var_old(k) - deltaL
       else
          aR(k) = r(k)
          aL(k) = r(k-1)
       end if
       dL(k) = Half * (aR(k) - aL(k))
       dR(k) = Half * (aR(k) + aL(k)) - var_old(k)
    end do

    do k = 1, N-1
       dz = z_new(k) - z_old(k)
       if (dz > Zero) then
          alpha = Hz(k+1)
          cff  = aL(k+1)
          cffL = dL(k+1)
          cffR = dR(k+1)
       else
          alpha = -Hz(k)
          cff = aR(k)
          cffL = -dL(k)
          cffR = dR(k)
       end if
       alpha = dz/alpha
       FC(k) = dz * (cff + alpha*(cffL - cffR*(Three - Two*alpha)))
    end do
    FC(0) = Zero
    FC(N) = Zero

    do k = 1, N
       var_new(k) = (Hz(k)*var_old(k) + FC(k)-FC(k-1)) / (z_new(k) - z_new(k-1))
    end do
  end subroutine remap2W

  subroutine remap4 (N, var_new, z_new, var_old, z_old, neumann)
    ! Parabolic WENO enhanced by quartic power-law reconciliation step. 
    ! Profile in each grid box is then given by a quartic polynomial.
    !
    ! (1) continuity of both value and first derivative at each interface
    ! (2) essentially non-oscillatory
    implicit none
    integer                 :: N
    real(8), dimension(1:N) :: var_new, var_old
    real(8), dimension(0:N) :: z_new, z_old
    logical                 :: neumann

    integer                 :: k
    real(8)                 :: alpha, Ampl, cff, cffL, cffR, deltaL, deltaR, dz, Hdd, rr
    real(8), parameter      :: eps=1.0d-8
    real(8), parameter      :: Half=0.5_8, OneFifth=0.2_8, ThreeHalfth=1.5_8 
    real(8), parameter      :: Zero=0.0_8, One=1.0_8, Two=2.0_8, Three=3.0_8, Four=4.0_8, Six=6.0_8
    real(8), dimension(1:N) :: Hz
    real(8), dimension(0:N) :: aL, aR, dL, dR, d, FC, r, r1
    !
    ! Parabolic WENO reconstruction: The second and third loops below
    !---------- ---- --------------- compute left and right side limits
    ! aL,aR of the field var_old assuming monotonized parabolic distributions
    ! within each grid box. Also computed are dL,dR, which are then used
    ! as a measure of quadratic variation during sabsequent WENO
    ! reconciliation of side limits.    
    !
    do k = 1, N
       Hz(k) = z_old(k) - z_old(k-1)
    end do

    do k = N-1, 1, -1
       FC(k) = One / (Hz(k+1) + Hz(k))
       d(k)  = FC(k) * (var_old(k+1) - var_old(k))
    end do

    do k = 2, N-1
       deltaR = Hz(k)*d(k)
       deltaL = Hz(k)*d(k-1)

       if (deltaR*deltaL < Zero) then
          deltaR = Zero
          deltaL = Zero
       end if
       cff = Hz(k-1) + Two*Hz(k) + Hz(k+1)
       cffR = cff*d(k)
       cffL = cff*d(k-1)
       if (abs(deltaR) > abs(cffL)) deltaR = cffL
       if (abs(deltaL) > abs(cffR)) deltaL = cffR

       cff = (deltaR - deltaL) / (Hz(k-1) + Hz(k) + Hz(k+1))
       deltaR = deltaR - cff*Hz(k+1)
       deltaL = deltaL + cff*Hz(k-1)

       aR(k) = var_old(k) + deltaR
       aL(k) = var_old(k) - deltaL

       dR(k) = (Two*deltaR - deltaL)**2
       dL(k) = (Two*deltaL - deltaR)**2
    end do

    aL(N) = aR(N-1)
    aR(N) = Two*var_old(N) - aL(N)

    dR(N) = (Two*aR(N)+aL(N) - Three*var_old(N))**2
    dL(N) = (Three*var_old(N) - Two*aL(N) - aR(N))**2

    aR(1) = aL(2)
    aL(1) = Two*var_old(1) - aR(1)

    dR(1) = (Two*aR(1) + aL(1) - Three*var_old(1))**2
    dL(1) = (Three*var_old(1) - Two*aL(1) - aR(1))**2

    do k = 1, N-1
       deltaL = max (dL(k),   eps)
       deltaR = max (dR(k+1), eps)
       r1(k) = (deltaR*aR(k) + deltaL*aL(k+1)) / (deltaR + deltaL)
    end do      !--> discard aR,aL,dR,dL

    if (NEUMANN) then
       r1(N) = ThreeHalfth*var_old(N) - Half*r1(N-1)
       r1(0) = ThreeHalfth*var_old(1) - Half*r1(1)
    else
       r1(N) = Two*var_old(N) - r1(N-1)
       r1(0) = Two*var_old(1) - r1(1)
    end if
    !
    ! Power-law reconciliation: Starts with computation of side
    !------ --- --------------- limits dR,dL of the first derivative
    ! assuming parabolic distributions within each grid box. In this
    ! version of the code, before doing so (see "else" branch of 3-way
    ! switch below), in situation when interfacial deviations deltaR
    ! and deltaL are differ by more than a factor of two (hence
    ! monotonic parabolic fit becomes impossible), the parabolic
    ! assumption is switched to power-law function,  such that its
    ! derivative is zero at one end and, consequently, larger than
    ! taht of (would be) limited parabolic on the other end. The basic
    ! parabolic version of the code is commented out, but left here
    ! for reference.
    !
    do k = 1, N
       deltaR = r1(k) - var_old(k)
       deltaL = var_old(k) - r1(k-1)
       cff = deltaR * deltaL
       if (cff > eps) then
          cff = (deltaR + deltaL)/cff
       else
          cff = Zero
       end if
       cffL = cff * deltaL
       cffR = cff * deltaR

       if (cffL > Three) then
          cffL = cffL * deltaL
          cffR = Zero
       elseif (cffR > Three) then
          cffL = Zero
          cffR = cffR * deltaR
       else
          cffL = Four*deltaL - Two*deltaR
          cffR = Four*deltaR - Two*deltaL
       end if
       cff = One / Hz(k)
       dR(k) = cff * cffR
       dL(k) = cff * cffL
    end do
    !
    ! Compute final value of derivative at each interface by reconciling
    ! two side limits dR(k) and dL(k+1) coming from adjacent grid boxes.
    ! the difference between these two also causes change of interfacial
    ! value r(k) by Ampl. The commented code (left here for reference)
    ! computes the exact value of Ampl assuming power law reconciliation
    ! and solving associated quadratic equation. The code segment below  
    ! it corresponds to Pade fit to exact solution, which avoids
    ! computation of sqrt for the sake of computational efficiency.
    !
    do k = N-1, 1, -1
       d(k) = FC(k) * (Hz(k+1)*dL(k+1) + Hz(k)*dR(k))

       cffR = 8*(dR(k)   + Two*dL(k))
       cffL = 8*(dL(k+1) + Two*dR(k+1))
       if (abs(d(k)) > abs(cffR)) d(k) = cffR
       if (abs(d(k)) > abs(cffL)) d(k) = cffL

       if ((dL(k+1)-dR(k))*(var_old(k+1)-var_old(k)) > Zero) then
          Hdd = Hz(k) * (d(k) - dR(k))
          rr = var_old(k) - r1(k-1)
       else
          Hdd = Hz(k+1) * (dL(k+1) - d(k))
          rr = r1(k+1) - var_old(k+1)
       end if
       rr = abs (rr)

       Ampl = OneFifth * Hdd * rr
       Hdd = abs (Hdd)
       cff = rr**2 + 0.0763636363636363636*Hdd * (rr + 0.004329004329004329*Hdd)
       if (cff > eps) then
          Ampl = Ampl * (rr + 0.0363636363636363636*Hdd) / cff
       else
          Ampl = Zero
       end if
       r(k) = r1(k) + Ampl
    end do
    
    if (NEUMANN) then
       r(0) = ThreeHalfth*var_old(1) - Half*r(1)
       r(N) = ThreeHalfth*var_old(N) - Half*r(N-1)
       d(0) = Zero
       d(N) = Zero
    else
       r(0) = Two*var_old(1) - r(1)
       r(N) = Two*var_old(N) - r(N-1)
       d(0) = d(1)
       d(N) = d(N-1)
    end if

    ! Remapping step: This operation consists essentially of three
    !---------- ----- stages: (1) within each grid box compute three
    ! auxiliary fields, which are analogous to approximations for the
    ! second, third and fourth derivatives (stored in scratch arrays
    ! aR, dR and dL respectively); then (2) compute interfacial fluxes
    ! FC; and (3) apply these fluxes to complete remapping step.
    !
    do k = 1, N
       aR(k) = r(k) + r(k-1) - Two*var_old(k)
       dR(k) = Hz(k) * (d(k)+d(k-1)) - Two*(r(k)-r(k-1))
       dL(k) = Hz(k) * (d(k)-d(k-1)) - Six*aR(k)
    end do

    do k = 1, N-1
       dz = z_new(k) - z_old(k)
       if (dz > Zero) then
          alpha = Hz(k+1)   
          cff   = aR(k+1)
          cffR  = dR(k+1)
          cffL  = dL(k+1)
       else
          alpha = -Hz(k)
          cff   =  aR(k)
          cffR  = -dR(k)
          cffL  =  dL(k)
       end if
       alpha = dz / alpha
       FC(k) = dz * (r(k) + Half*dz*d(k) + alpha**2*(cff - cffR*(0.5_8-0.25*alpha) + cffL*(1.0_8-alpha*(1.25_8-0.5*alpha))))
    end do
    FC(0) = 0.0_8
    FC(N) = 0.0_8

    do k = 1, N
       var_new(k) = (Hz(k)*var_old(k) + FC(k)-FC(k-1)) / (z_new(k) - z_new(k-1))
    end do
  end subroutine remap4
end module remap_mod
