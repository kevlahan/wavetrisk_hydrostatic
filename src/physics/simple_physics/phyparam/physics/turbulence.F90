module turbulence
  use phys_const, only : Emin_turb, Stefan
  implicit none
  save
  private
  real    :: Karman    = 4e-1         ! Karman constant
  real    :: dVdZ2_min = 1e-6         ! minimum vertical gradient^2 of speed
  real    :: Ri_c      = 4e-1         ! critical Reynolds number
  real    :: LmixMin   = 1e+2         ! minimum mixing length
  real    :: p_0       = 1e+5         ! reference pressure [Pa]
  real    :: Kv_max    = 1e+0         ! maximum turbulent diffusivity
  public ::  dVdZ2_min, LmixMin, p_0, Ri_c, Vdif
contains
  subroutine Vdif (ngrid, nlay, mask, pTime, pTimestep, pCapCal, pEmis, pFluxSrf, pZ0, pPlay, pPint, pZlay, pZint, pUmag, &
       pU, pV, pW, pTheta, pTsrf)
    use phys_const, only: g, R, Rcp, Cpp
    !============================================================================================
    !
    !   Complete physics time integration using input non-diffusion physics tendencies and computed
    !   vertical turbulent diffusion tendencies.
    !
    !   Physics is advanced using explicit Forward Euler split step for non-diffusion
    !   physics tendencies followed by implicit Backwards Euler for computed
    !   vertical diffusion tendencies:
    !
    !   X(t+dt) = X(t) + dt * T_phys(X(t),t)  +  dt * T_diff(X(t+dt),t+dt)
    !
    !        ->    Y = X(t) + dt * T_phys(X(t),t)
    !
    !              X(t+dt) = Y + dt * T_diff(Y,t+dt)
    !
    !
    !   where  T_phys and T_diff are physics and vertical diffusion tendencies.
    !
    !============================================================================================
    ! Input variables
    integer,                       intent(in) :: ngrid         ! number of columns
    integer,                       intent(in) :: nlay          ! number of vertical layers
    integer,                       intent(in) :: mask          ! adaptive grid mask value
    real,                          intent(in) :: pTime         ! current time t
    real,                          intent(in) :: pTimestep     ! time step dt
    real, dimension(ngrid),        intent(in) :: pCapCal       ! surface conduction flux
    real, dimension(ngrid),        intent(in) :: pEmis         ! emissivity of surface
    real, dimension(ngrid),        intent(in) :: pFluxSrf      ! surface flux
    real, dimension(ngrid),        intent(in) :: pZ0           ! surface roughness length   [m]
    real, dimension(ngrid,nlay),   intent(in) :: pPlay         ! pressure at layers         [Pa]
    real, dimension(ngrid,nlay+1), intent(in) :: pPint         ! pressure at interfaces     [Pa]
    real, dimension(ngrid,nlay),   intent(in) :: pZlay         ! z coordinate at layers     [m]
    real, dimension(ngrid,nlay+1), intent(in) :: pZint         ! z coordinate at interfaces [m]
    real, dimension(ngrid,nlay),   intent(in) :: pUmag         ! speed at layers [m/s]

    ! Output variables: advanced using backwards Euler step
    real, dimension(ngrid,nlay),   intent(inout) :: pU, pV, pW ! velocities [m/s]
    real, dimension(ngrid,nlay),   intent(inout) :: pTheta     ! potential temperature       [K]
    real, dimension(ngrid),        intent(inout) :: pTsrf      ! surface temperature         [K]

    ! Local variables
    integer                       :: ig, il
    integer                       :: it1, it2
    real, dimension(ngrid)        :: rho, T, z2, zdPlanck, zCd_theta, zCdv, zU2
    real, dimension(ngrid,nlay)   :: dZ, zb, zb0, z1, zc, zd

    real, dimension(ngrid,nlay+1) :: zKv, zKtheta              ! turbulent diffusivities

    !  Initialization: computation of rho*dz and dt * rho/dz = dt * rho**2 G/dP with rho = P/(R T) = p/(R Theta*(p/ps)**kappa)
    dZ(:,1:nlay) = (pPint(:,1:nlay) - pPint(:,2:nlay+1)) / G

    zb0(:,1) = pTimestep * pPint(:,1) / (R * pTsrf)
    do il = 2, nlay
       T         = 0.5 * (pTheta(:,il) + pTheta(:,il-1)) * (pPint(:,il)/pPint(:,1))**Rcp ! temperature at interface below level il
       rho       = pPint(:,il) / (T*R)                                            ! density at lower interface below level il
       zb0(:,il) = pTimestep * rho**2 * G / (pPlay(:,il-1) - pPlay(:,il))         ! dt * rho/dz at interface below level il
    end do

    ! Surface drag coefficients Cd
    call Vdif_Cd (ngrid, mask, pZ0, G, pZlay(:,1), pUmag(:,1), pTheta(:,1), pTsrf, zCdv, zCd_theta)

    ! Turbulent diffusion coefficients K
    call Vdif_K (ngrid, nlay, mask, G, pZlay, pZint, pZ0, pUmag, pTheta, zKv, zKtheta)

    ! Vertical integration for turbulent diffusion of velocities
    zb(:,1)      = zCdv          * zb0(:,1)       ! boundary condition
    zb(:,2:nlay) = zKv(:,2:nlay) * zb0(:,2:nlay)  ! dt * Kv * rho/dz

    z1(:,nlay) = 1.0 / (dZ(:,nlay) + zb(:,nlay))
    zd(:,nlay) = zb(:,nlay) * z1(:,nlay)

    do il = nlay-1, 1, -1
       z1(:,il) = 1.0 / (dZ(:,il) + zb(:,il) + zb(:,il+1) * (1.0 - zd(:,il+1)))
       zd(:,il) = zb(:,il) * z1(:,il)
    end do

    ! Vertical integration of velocities
    call velo_integration (ngrid, nlay, dZ, z1, zb, zd, pU)
    call velo_integration (ngrid, nlay, dZ, z1, zb, zd, pV)
    call velo_integration (ngrid, nlay, dZ, z1, zb, zd, pW)

    ! Vertical integration for potential temperature
    zb(:,1)      = zCd_theta * zb0(:,1)               ! boundary condition
    zb(:,2:nlay) = zKtheta(:,2:nlay) * zb0(:,2:nlay)

    z1(:,nlay) = 1.0 / (dZ(:,nlay) + zb(:,nlay))
    zd(:,nlay) = zb(:,nlay) * z1(:,nlay)

    do il = nlay-1, 1, -1
       z1(:,il) = 1.0 / (dZ(:,il) + zb(:,il) + zb(:,il+1) * (1.0 - zd(:,il+1)))
       zd(:,il) = zb(:,il) * z1(:,il)
    end do

    zc(:,nlay) = dZ(:,nlay) * pTheta(:,nlay) * z1(:,nlay)
    do il = nlay-1, 1, -1
       zc(:,il) = (dZ(:,il) * pTheta(:,il)+ zb(:,il+1) * zc(:,il+1)) * z1(:,il)
    end do

    ! Add Planck contribution to implicit scheme
    zdPlanck = pTimestep * 4.0 * Stefan * pEmis * pTsrf**3

    ! Surface temperature at t+dt
    z1(:,1) = pCapCal * pTsrf + Cpp * zb(:,1) * zc(:,1)         + zdPlanck * pTsrf + pFluxSrf * pTimestep
    z2      = pCapCal         + Cpp * zb(:,1) * (1.0 - zd(:,1)) + zdPlanck
    pTsrf   = z1(:,1) / z2

    ! Potential temperature at t+dt
    pTheta(:,1) = zc(:,1) + zd(:,1) * pTsrf
    do il = 2, nlay
       pTheta(:,il) = zc(:,il) + zd(:,il) * pTheta(:,il-1)
    end do
  end subroutine Vdif

  subroutine velo_integration (ngrid, nlay, dZ, z1, zb, zd, vel)
    ! Vertical integration of velocity component vel
    implicit none
    integer                     :: ngrid, nlay
    real, dimension(ngrid,nlay) :: dZ, z1, zb, zd, vel

    integer                     :: il
    real, dimension(ngrid,nlay) :: zc

    zc(:,nlay) = vel(:,nlay) * z1(:,nlay) * dZ(:,nlay)
    do il = nlay-1, 1, -1
       zc(:,il) = (vel(:,il) * dZ(:,il) + zb(:,il+1) * zc(:,il+1)) * z1(:,il)
    end do

    vel(:,1) = zc(:,1)
    do il = 2, nlay
       vel(:,il) = zc(:,il) + zd(:,il) * vel(:,il-1)
    end do
  end subroutine velo_integration
  
  subroutine Vdif_Cd (ngrid, mask, pZ0, pG, pZ, pUmag, pTheta, pTs, pCd_v, pCd_theta)
    ! Simplied formulation  from  J-F Louis (1979) (Boundary Layer Meteorology 17, 187-202)

    ! Input
    integer,                intent(in) :: ngrid  ! number of columns
    integer,                intent(in) :: mask   ! adaptive grid mask
    real,                   intent(in) :: pG     ! gravity                                [m/s^2]
    real, dimension(ngrid), intent(in) :: pZ0    ! surface roughness length scale         [m]
    real, dimension(ngrid), intent(in) :: pZ     ! height of the first atmospheric layer  [m]
    real, dimension(ngrid), intent(in) :: pUmag  ! speed in surface layer                 [m/s]
    real, dimension(ngrid), intent(in) :: pTheta ! potential temperature in surface layer [K]
    real, dimension(ngrid), intent(in) :: pTs    ! temperature at surface                 [K]

    ! Output
    real, dimension(ngrid), intent(out) :: pCd_v     ! drag coefficient for velocity
    real, dimension(ngrid), intent(out) :: pCd_theta ! drag coefficient for potential temperature

    integer         :: ig
    real            :: a, c1_theta, c1_u, F_theta, F_u, Ri, z1
    
    real, parameter :: R = 0.74 ! (velocity drag) / (temperature drag) empirical relation

    ! Parameters fitted to implicit relationship between Monin-Obukhov scale height L and 
    ! layer bulk Richardson number
    real, parameter :: b = 4.7, C_u = 7.4, C_theta = 5.3

    do ig = 1, ngrid
       z1  = pZ(ig) / pZ0(ig)

       a = Karman / log (z1)

       c1_u     = 2.0 * C_u     * a**2 * b * sqrt (z1)
       c1_theta = 2.0 * C_theta * a**2 * b * sqrt (z1)
       
       Ri = pG * pZ(ig) * (pTheta(ig) - pTs(ig)) / (pTheta(ig) * (pUmag(ig)**2 + Emin_turb)) ! bulk Richardson number at surface

       if (Ri < 0.0) then ! unstable regime
          F_u     = 1.0 - 2.0 * b * Ri / (1.0 + c1_u     * sqrt (-Ri))
          F_theta = 1.0 - 2.0 * b * Ri / (1.0 + c1_theta * sqrt (-Ri))
       else               ! stable regime
          F_u     = 1.0 / (1.0 + b * Ri)**2
          F_theta = F_u
       end if

       ! Drag coefficients
       pCd_v(ig)     = a**2 * pUmag(ig) * F_u   
       pCd_theta(ig) = a**2 * pUmag(ig) * F_theta / R
    end do
  end subroutine Vdif_Cd

  subroutine Vdif_K (ngrid, nlay, mask, pG, pZlay, pZint, pZ0, pUmag, pTheta, pKv, pKtheta)
    ! Turbulent diffusion coefficients for velocity (pKU) and potential temperature (pKtheta) using mixing length model
    ! (assumes Pr = 1 so that pKtheta = pKU: same diffusion coefficients for velocity and potential temperature)
    implicit none
    integer,                       intent(in)  :: ngrid        ! number of columns
    integer,                       intent(in)  :: nlay         ! number of layers
    integer,                       intent(in)  :: mask         ! adaptive grid mask
    real,                          intent(in)  :: pG           ! gravitational acceleration
    real, dimension(ngrid),        intent(in)  :: pZ0          ! roughness length
    real, dimension(ngrid,nlay),   intent(in)  :: pZlay        ! height of layers
    real, dimension(ngrid,nlay+1), intent(in)  :: pZint        ! height of interfaces
    real, dimension(ngrid,nlay+1), intent(in)  :: pUmag        ! speed
    real, dimension(ngrid,nlay),   intent(in)  :: pTheta       ! potential temperature at layers

    real, dimension(ngrid,nlay+1), intent(out) :: pKtheta, pKv ! turbulent diffusivities at interfaces

    integer :: ig, il
    real    :: Lmix, Ri
    real    :: dTheta, dU, dUdZ2, Theta, N2
    real    :: dZ, z1

    ! Turbulent diffusivities at interfaces
    pKv(:,1)      = 0d0
    pKv(:,nlay+1) = 0d0
    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pZint(ig,il) + pZ0(ig)

          Lmix = Karman * z1 / (1.0 + Karman * z1 / LmixMin) ! mixing length

          dTheta = pTheta(ig,il) - pTheta(ig,il-1)
          dU     = pUmag (ig,il) - pUmag (ig,il-1)
          dZ     = pZlay (ig,il) - pZlay (ig,il-1)
          dUdZ2  = (dU/dZ)**2

          if (dUdZ2 > dVdZ2_min) then
             Theta = 0.5 * (pTheta(ig,il) + pTheta(ig,il-1)) ! potential temperature at interface below layer il
             N2    = pG/Theta * dTheta/dZ                    ! Brunt-Vaisala frequency squared
             Ri    = N2 / dUdZ2                              ! gradient Richardson number

             ! Turbulent diffusivity
             pKv(ig,il) = Lmix * sqrt (max (Lmix**2 * dUdZ2 * (1.0 - Ri/Ri_c), Emin_turb))
          else
             pKv(ig,il) = Lmix * sqrt (Emin_turb)
          end if
       end do
    end do
    pKtheta = pKv
  end subroutine Vdif_k
end module turbulence
