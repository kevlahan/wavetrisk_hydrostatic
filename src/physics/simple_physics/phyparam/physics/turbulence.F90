module turbulence
  implicit none
  save
  private
  integer :: ADJZONE   = 8    ! adjzone mask
  real    :: Karman    = 4e-1 ! Karman constant
  real    :: Emin_turb = 1e-8 ! minimum turbulent kinetic energy
  real    :: dVdZ2_min = 1e-6 ! minimum vertical gradient^2 of speed
  real    :: Ri_c      = 4e-1 ! critical Reynolds number
  real    :: LmixMin   = 1e+2 ! minimum mixing length
  real    :: p_0       = 1e5  ! reference pressure [Pa]
  public :: ADJZONE, dVdZ2_min, Emin_turb, LmixMin, p_0, Ri_c, Vdif
contains
  subroutine Vdif (ngrid, nlay, mask, pTime, pTimestep, pCapCal, pEmis, pFluxSrf, pZ0, pPlay, pPlev, pZlay, pZlev, &
       pU, pV, pTheta, pTsrf)
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
    integer,                       intent(in) :: ngrid     ! number of columns
    integer,                       intent(in) :: nlay      ! number of vertical layers
    integer,                       intent(in) :: mask      ! adaptive grid mask value
    real,                          intent(in) :: pTime     ! current time t
    real,                          intent(in) :: pTimestep ! time step dt
    real, dimension(ngrid),        intent(in) :: pCapCal   ! surface conduction flux
    real, dimension(ngrid),        intent(in) :: pEmis     ! emissivity of surface
    real, dimension(ngrid),        intent(in) :: pFluxSrf  ! surface flux
    real, dimension(ngrid),        intent(in) :: pZ0       ! surface roughness length   [m]
    real, dimension(ngrid,nlay),   intent(in) :: pPlay     ! pressure at layers         [Pa]
    real, dimension(ngrid,nlay+1), intent(in) :: pPlev     ! pressure at interfaces     [Pa]
    real, dimension(ngrid,nlay),   intent(in) :: pZlay     ! z coordinate at layers     [m]
    real, dimension(ngrid,nlay+1), intent(in) :: pZlev     ! z coordinate at interfaces [m]

    ! Output variables: advanced using backwards Euler step
    real, dimension(ngrid,nlay),   intent(inout) :: pU, pV ! zonal and meridional velocities [m/s]
    real, dimension(ngrid,nlay),   intent(inout) :: pTheta ! potential temperature           [K]
    real, dimension(ngrid),        intent(inout) :: pTsrf  ! surface temperature             [K]

    ! Local variables
    integer                       :: ig, il
    integer                       :: it1, it2
    real                          :: z4st, zcst1
    real, dimension(ngrid)        :: rho, T, z2, zdPlanck, zCdh, zCdv, zU2
    real, dimension(ngrid,nlay)   :: dZ, zb, zb0, z1, zc, zd

    real, dimension(ngrid,nlay+1) :: zKv, zKh              ! turbulent diffusivities

    !  Initialization: computation of rho*dz and dt * rho/dz = dt * rho**2 G/dP with rho = P/(R T)  = p/(R Theta*(p/ps)**kappa)
    dZ(:,1:nlay) = (pPlev(:,1:nlay) - pPlev(:,2:nlay+1)) / G

    zb0(:,1) = pTimestep * pPlev(:,1) / (R * pTsrf)                          ! dt * ho/dz at surface
    do il = 2, nlay
       T   = 0.5 * (pTheta(:,il) + pTheta(:,il-1)) * (pPlev(:,il)/p_0)**Rcp  ! temperature at interface below level il
       rho = pPlev(:,il) / (T*R)                                             ! density at lower interface below level il

       zb0(:,il) = pTimestep * rho**2 * G / (pPlay(:,il-1) - pPlay(:,il))    ! dt * rho/dz at interface below level il
    end do

    ! Surface drag coefficients Cd
    call Vdif_Cd (ngrid, mask, pZ0, G, pZlay(:,1), pU(:,1), pV(:,1), pTheta(:,1), pTsrf, zCdv, zCdh)

    ! Turbulent diffusion coefficients K
    call Vdif_K (ngrid, nlay, mask, pTimestep, G, pZlev, pZlay, pZ0, pU, pV, pTheta, zKv, zKh)

    ! Vertical integration for turbulent diffusion of velocities
    zb(:,1)      = zCdv          * zb0(:,1)       ! boundary condition
    zb(:,2:nlay) = zKv(:,2:nlay) * zb0(:,2:nlay)  ! dt * Kv * rho/dz

    z1(:,nlay) = 1.0 / (dZ(:,nlay) + zb(:,nlay))
    zd(:,nlay) = zb(:,nlay) * z1(:,nlay)

    do il = nlay-1, 1, -1
       z1(:,il) = 1.0 / (dZ(:,il) + zb(:,il) + zb(:,il+1) * (1.0 - zd(:,il+1)))
       zd(:,il) = zb(:,il) * z1(:,il)
    end do

    ! Zonal velocity
    zc(:,nlay) = pU(:,nlay) * z1(:,nlay) * dZ(:,nlay)
    do il = nlay-1, 1, -1
       zc(:,il) = (pU(:,il) * dZ(:,il) + zb(:,il+1) * zc(:,il+1)) * z1(:,il)
    end do

    pU(:,1) = zc(:,1)
    do il = 2, nlay
       pU(:,il) = zc(:,il) + zd(:,il) * pU(:,il-1)
    end do

    ! Meridional velocity
    zc(:,nlay) = pV(:,nlay) * z1(:,nlay) * dZ(:,nlay)
    do il = nlay-1, 1, -1
       zc(:,il) = (pV(:,il) * dZ(:,il) + zb(:,il+1) * zc(:,il+1)) * z1(:,il)
    end do

    pV(:,1) = zc(:,1)
    do il = 2, nlay
       pV(:,il) = zc(:,il) + zd(:,il) * pV(:,il-1)
    end do

    ! Vertical integration for potential temperature
    zb(:,1)      = zCdh * zb0(:,1)               ! boundary condition
    zb(:,2:nlay) = zKh(:,2:nlay) * zb0(:,2:nlay)

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
    z4st = 4.0 * 5.67e-8 * pTimestep
    zdPlanck = z4st * pEmis * pTsrf**3

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

  pure subroutine Vdif_Cd (ngrid, mask, pZ0, pG, pZ, pU, pV, pTheta, pTs, pCdv, pCdh)
    ! Surface drag coefficient using approach developed by Louis for ECMWF.

    ! Input
    integer,                intent(in) :: ngrid  ! number of columns
    integer,                intent(in) :: mask   ! adaptive grid mask
    real,                   intent(in) :: pG     ! gravity                                 [m/s^2]
    real, dimension(ngrid), intent(in) :: pZ0    ! surface roughness length scale          [m]
    real, dimension(ngrid), intent(in) :: pZ     ! height of the first atmospheric layer   [m]
    real, dimension(ngrid), intent(in) :: pU     ! zonal velocity in surface layer         [m/s]
    real, dimension(ngrid), intent(in) :: pV     ! meridional velocity in surface layer    [m/s]
    real, dimension(ngrid), intent(in) :: pTheta ! potential temperature in surface layer  [K]
    real, dimension(ngrid), intent(in) :: pTs    ! temperature at surface                  [K]

    ! Output
    real, dimension(ngrid), intent(out) :: pCdv ! drag coefficient for velocity
    real, dimension(ngrid), intent(out) :: pCdh ! drag coefficient for potential temperature

    integer                :: ig
    real                   :: z1, Ri, Cd0, Z
    real, dimension(ngrid) :: U

    real, parameter :: b = 5.0, c = 5.0, d = 5.0, c2b = 2.0 * b, c3bc = 3.0 * b * c, c3b = 3.0 * b

    do ig = 1, ngrid
       U(ig) = sqrt (pU(ig)**2 + pV(ig)**2 + Emin_turb)

       z1  = 1.0 + pZ(ig) / pZ0(ig)

       Cd0 = (Karman/log(z1))**2 * U(ig)

       ! Gradient Richardson number at surface (T_s = Theta_s)
       Ri = 0.5 * pG/pTs(ig) * (pTheta(ig) - pTs(ig)) / U(ig)**2 * pZ(ig)

       if (Ri < 0.0) then
          z1 = b * Ri / (1.0 + c3bc * Cd0 * sqrt (- z1 * Ri))

          pCdv(ig) = Cd0 * (1.0 - 2.0 * z1)
          pCdh(ig) = Cd0 * (1.0 - 3.0 * z1)
       else
          Z = sqrt (1.0 + d * Ri)

          pCdv(ig) = Cd0 / (1.0 + c2b * Ri / Z)
          pCdh(ig) = Cd0 / (1.0 + c3b * Ri * Z)
       end if
    end do
    pCdv = pCdv * U
    pCdh = pCdh * U
  end subroutine Vdif_Cd

  subroutine Vdif_K (ngrid, nlay, mask, pTimestep, pG, pZlev, pZlay, pZ0, pU, pV, pTheta, pKv, pKh)
    ! Turbulent diffusion coefficients for velocity (pKU) and potential temperature (pKh) using mixing length model
    ! (assumes Pr = 1 so that pKh = pKU: same diffusion coefficients for velocity and potential temperature)
    integer,                       intent(in)  :: ngrid, nlay, mask
    real,                          intent(in)  :: pTimestep, pG
    real, dimension(ngrid),        intent(in)  :: pZ0           ! roughness length
    real, dimension(ngrid,nlay),   intent(in)  :: pZlay
    real, dimension(ngrid,nlay),   intent(in)  :: pU, pV, pTheta
    real, dimension(ngrid,nlay+1), intent(in)  :: pZlev

    real, dimension(ngrid,nlay+1), intent(out) :: pKh, pKv      ! turbulent diffusivities

    integer :: ig, il
    real    :: Lmix, Ri
    real    :: dU, dV, dH, H, N2
    real    :: dVdZ2, dZ, z1

    ! Turbulent diffusivities at lower interfaces corresponding to layer il
    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pZlev(ig,il) + pZ0(ig)

          Lmix = Karman * z1 / (1.0 + Karman * z1 / LmixMin) ! mixing length

          dU = pU(ig,il)     - pU(ig,il-1)
          dV = pV(ig,il)     - pV(ig,il-1)
          dH = pTheta(ig,il) - pTheta(ig,il-1)

          dZ = pZlay(ig,il) - pZlay(ig,il-1)

          dVdZ2 = (dU**2 + dV**2) / dZ**2

          if (dVdZ2 > dVdZ2_min) then
             H = 0.5 * (pTheta(ig,il) + pTheta(ig,il-1)) ! potential temperature at interface below layer il

             N2 = pG/H * dH/dZ                   ! Brunt-Vaisala frequency squared

             Ri = N2 / dVdZ2                     ! gradient Richardson number

             ! Turbulent diffusivity
             pKv(ig,il)= Lmix * sqrt (max (Lmix**2 * dVdZ2 * (1.0 - Ri/Ri_c), Emin_turb))
          else
             pKv(ig,il) = Lmix * sqrt (Emin_turb)
          end if
       end do
    end do
    pKh = pKv
  end subroutine Vdif_k
end module turbulence
