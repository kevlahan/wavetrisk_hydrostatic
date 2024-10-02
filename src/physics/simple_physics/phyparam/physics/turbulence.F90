module turbulence
#include "use_logging.h"
  implicit none
  save
  private
  real :: Karman    = 4e-1 ! Karman constant
  real :: Emin_turb = 1e0!1e-8 ! minimum turbulent kinetic energy
  real :: dVdZ_min  = 1e-3 ! minimum vertical gradient of speed
  real :: Ri_c      = 4e-1 ! critical Reynolds number
  real :: LmixMin   = 1e+2  ! minimum mixing length

  public :: dVdZ_min, Emin_turb, LmixMin, Ri_c, Vdif
contains
  subroutine Vdif (ngrid, nlay, ptime, ptimestep, &
       pCapCal, pZ0, pPlay, pPlev, pZlay, pZlev, &
       pU, pV, pH, &
       pTsrf, pEmis, &
       pdUfi, pdVfi, pdHfi, &
       pFluxSrf, &
       pdUdif, pdVdif, pdHdif, pdTSrf, &
       lwrite)
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
    integer,                       intent(in) :: ngrid                   ! number of columns
    integer,                       intent(in) :: nlay                    ! number of vertical layers

    real,                          intent(in) :: pTime                   ! current time t
    real,                          intent(in) :: pTimestep               ! time step dt

    ! Surface variables
    real, dimension(ngrid),        intent(in) :: pCapCal                 ! surface conduction flux
    real, dimension(ngrid),        intent(in) :: pEmis                   ! emissivity of surface
    real, dimension(ngrid),        intent(in) :: pFluxSrf                ! surface flux
    real, dimension(ngrid),        intent(in) :: pZ0                     ! surface roughness length [m]
    real, dimension(ngrid),        intent(in) :: pTsrf                   ! surface temperature      [K]

    ! Layer centre variables
    real, dimension(ngrid,nlay),   intent(in) :: pU, pV, pH              ! prognostic variables at time t from physics
    real, dimension(ngrid,nlay),   intent(in) :: pdUfi, pdVfi, pdHfi     ! physics tendencies at time t (explicit scheme)
    real, dimension(ngrid,nlay),   intent(in) :: pPlay                   ! pressure at middle of layers     [Pa]
    real, dimension(ngrid,nlay),   intent(in) :: pZlay                   ! z coordinate at middle of layers [m]

    ! Layer interface variables (interface pressure and z coordinate)
    real, dimension(ngrid,nlay+1), intent(in) :: pPlev                   ! pressure at interfaces between layers [Pa]
    real, dimension(ngrid,nlay+1), intent(in) :: pZlev                   ! pressure at middle of layers [Pa]

    logical,                       intent(in) :: lwrite                  ! write diagnostics flag

    ! Output variables: pseudo-tendencies at t+dt
    real, dimension(ngrid),        intent(out) :: pdTsrf                 ! surface temperature
    real, dimension(ngrid,nlay),   intent(out) :: pdUdif, pdVdif         ! zonal and meridional velocities
    real, dimension(ngrid,nlay),   intent(out) :: pdHdif                 ! potential temperature

    ! Local variables
    integer                       :: il, ig, ilay
    integer                       :: unit, ierr, it1, it2
    integer                       :: cluvdb, pUtdat, pUtvdim, setname, setvdim

    real                          :: z4st, zcst1
    real, dimension(ngrid)        :: z2, zdPlanck, zCdH, zCdV, zTsrf2, zU2
    real, dimension(ngrid,nlay)   :: zH, zU, zV                         ! intermediate prognostic variables at t+dt
    real, dimension(ngrid,nlay)   :: dZ, zb, zb0, z1, zc, zd
    real, dimension(ngrid,nlay+1) :: zKV, zKH

    !------------------------------------------------------------------------------------------------------
    !   Initializations:
    !
    !   Computation of rho*dz and dt * rho/dz = dt * rho**2 g/dp with rho = p/RT  = p/(R Theta) (p/ps)**kappa
    !
    !------------------------------------------------------------------------------------------------------

    ! Layer thicknesses
    dZ(:,1:nlay) = (pPlev(:,1:nlay) - pPlev(:,2:nlay+1)) / g

    zcst1 = 4.0 * g * pTimestep / R**2
    do il = 2, nlay
       zb0(:,il) = pPlev(:,il) * (pPlev(:,1) / pPlev(:,il))**Rcp / (pH(:,il-1) + pH(:,il))
       zb0(:,il) = zcst1 * zb0(:,il)**2 / (pPlay(:,il-1) - pPlay(:,il))
    end do
    zb0(:,1) = pTimestep * pPlev(:,1) / (R * pTsrf)

    if (lwrite) then
       ig = ngrid/2 + 1
       WRITELOG (*,*) 'Pression (mbar) altitude (km), u, v, theta, rho dz'
       do ilay = 1, nlay
          WRITELOG (*,*) 0.01 * pPlay(ig,ilay), 0.001 * pZlay(ig,ilay), pU(ig,ilay), pV(ig,ilay), pH(ig,ilay), dZ(ig,ilay)
       end do
       WRITELOG (*,*) 'Pression (mbar) altitude (km),zb'
       do il = 1, nlay
          WRITELOG (*,*) 0.01 * pPlev(ig,il), 0.001 * pZlev(ig,il), zb0(ig,il)
       end do
       LOG_DBG ('Vdif')
    endif

    !-----------------------------------------------------------------------------------
    !   2. Forward Euler split step using non-diffusion physics tendencies from input
    !-----------------------------------------------------------------------------------
    zU = pU + pdUfi * pTimestep
    zV = pV + pdVfi * pTimestep
    zH = pH + pdHfi * pTimestep

    !-----------------------------------------------------------------------
    !   3. Compute drag coefficients Cd
    !-----------------------------------------------------------------------
    call Vdif_Cd (ngrid, pZ0, G, pZlay(:,1), pU(:,1), pV(:,1), pH(:,1), pTsrf, zCdV, zCdH)

    !-----------------------------------------------------------------------
    !   4. Compute turbulent diffusion coefficients K
    !-----------------------------------------------------------------------
    call Vdif_K (ngrid, nlay, pTimestep, G, pZlev, pZlay, pZ0, pU, pV, pH, zKV, zKH)

    if (lwrite) then
       WRITELOG (*,*) 'Drag coefficients for velocity and potential temperature'
       WRITELOG (*,*) zCdV(ngrid/2+1), zCdH(ngrid/2+1)

       WRITELOG (*,*) 'Turbulent diffusivities for velocity and potential temperature'
       do il = 1, nlay
          WRITELOG (*,*) zKV(ngrid/2+1,il), zKH(ngrid/2+1,il)
       end do
       LOG_DBG ('Vdif')
    endif

    !-----------------------------------------------------------------------
    !   Vertical integration for turbulent diffusion of velocities
    !-----------------------------------------------------------------------
    zb(:,1)      = zCdV          * zb0(:,1)       ! boundary condition
    zb(:,2:nlay) = zKV(:,2:nlay) * zb0(:,2:nlay)  

    z1(:,nlay) = 1.0 / (dZ(:,nlay) + zb(:,nlay))
    zd(:,nlay) = zb(:,nlay) * z1(:,nlay)

    do ilay = nlay-1, 1, -1
       z1(:,ilay) = 1.0 / (dZ(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zd(:,ilay) = zb(:,ilay) * z1(:,ilay)
    end do
    
    ! Zonal velocity
    zc(:,nlay) = zU(:,nlay) * dZ(:,nlay) * z1(:,nlay)
    do ilay = nlay-1, 1, -1
       zc(:,ilay) = (zU(:,ilay) * dZ(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1(:,ilay)
    end do

    zU(:,1) = zc(:,1)
    do ilay = 2, nlay
       zU(:,ilay) = zc(:,ilay) + zd(:,ilay) * zU(:,ilay-1)
    end do

    
    ! Meridional velocity
    zc(:,nlay) = zV(:,nlay) * dZ(:,nlay) * z1(:,nlay)
    do ilay = nlay-1, 1, -1
       zc(:,ilay) = (zV(:,ilay) * dZ(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1(:,ilay)
    end do

    zV(:,1) = zc(:,1)
    do ilay = 2, nlay
       zV(:,ilay) = zc(:,ilay) + zd(:,ilay) * zV(:,ilay-1)
    end do

    !------------------------------------------------------------------------
    !   Vertical integration for potential temperature H
    !------------------------------------------------------------------------
    zb(:,1)      = zCdH * zb0(:,1)               ! boundary condition
    zb(:,2:nlay) = zKH(:,2:nlay) * zb0(:,2:nlay)
    
    z1(:,nlay) = 1.0 / (dZ(:,nlay) + zb(:,nlay))
    zd(:,nlay) = zb(:,nlay) * z1(:,nlay)

    do ilay = nlay-1, 1, -1
       z1(:,ilay) = 1.0 / (dZ(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zd(:,ilay) = zb(:,ilay) * z1(:,ilay)
    end do

    zc(:,nlay) = dZ(:,nlay) * zH(:,nlay) * z1(:,nlay)
    do ilay = nlay-1, 1, -1
       zc(:,ilay) = (dZ(:,ilay) * zH(:,ilay)+ zb(:,ilay+1) * zc(:,ilay+1)) * z1(:,ilay)
    end do

    !------------------------------------------------------------------------
    !   Add Planck contribution to implicit scheme
    !------------------------------------------------------------------------
    z4st = 4.0 * 5.67e-8 * pTimestep
    zdPlanck = z4st * pEmis * pTsrf**3

    !-----------------------------------------------------------------------
    !   Evolution of surface temperature
    !-----------------------------------------------------------------------
    z1(:,1) = pCapCal * pTsrf + Cpp * zb(:,1) * zc(:,1)         + zdPlanck * pTsrf + pFluxSrf * pTimestep
    z2      = pCapCal         + Cpp * zb(:,1) * (1.0 - zd(:,1)) + zdPlanck

    zTsrf2 = z1(:,1) / z2

    ! Diffusion tendency for potential temperature
    zH(:,1) = zc(:,1) + zd(:,1) * zTsrf2
    do ilay = 2, nlay
       zH(:,ilay) = zc(:,ilay) + zd(:,ilay) * zH(:,ilay-1)
    end do

    !-----------------------------------------------------------------------
    !   Complete pseudo-tendencies/solution at t+dt
    !-----------------------------------------------------------------------
    pdUdif = ((zU - pdUfi * pTimestep) - pU   ) / pTimestep
    pdVdif = ((zV - pdVfi * pTimestep) - pV   ) / pTimestep
    pdHdif = ((zH - pdHfi * pTimestep) - pH   ) / pTimestep
    pdTsrf = (zTsrf2                   - pTsrf) / pTimestep

    if (lwrite) then
       WRITELOG (*,*)
       WRITELOG (*,*) 'Diagnostique de la diffusion verticale'
       WRITELOG (*,*) 'Potential temperature avant et apres diffusion verticale'
       WRITELOG (*,*) pTsrf(ngrid/2+1), zTsrf2(ngrid/2+1)
       do il = 1, nlay
          WRITELOG (*,*) pH(ngrid/2+1,il), zH(ngrid/2+1,il)
       end do
       LOG_DBG ('Vdif')
    end if
  end subroutine Vdif

  pUre subroutine Vdif_Cd (ngrid, pZ0, pG, pZ, pU, pV, pH, pTs, pCdv, pCdh)
    !   Surface drag coefficient using approach developed by Louis for ECMWF.

    ! Input
    integer,                intent(in) :: ngrid ! number of columns
    real,                   intent(in) :: pG    ! gravity                                 [m/s^2]
    real, dimension(ngrid), intent(in) :: pZ0   ! surface roughness length scale          [m]
    real, dimension(ngrid), intent(in) :: pZ    ! height of the first atmospheric layer   [m]
    real, dimension(ngrid), intent(in) :: pU    ! zonal velocity in surface layer         [m/s]
    real, dimension(ngrid), intent(in) :: pV    ! meridional velocity in surface layer    [m/s]
    real, dimension(ngrid), intent(in) :: pH    ! potential temperature in surface layer  [K]
    real, dimension(ngrid), intent(in) :: pTs   ! temperature at surface                  [K]

    ! Output
    real, dimension(ngrid), intent(out) :: pCdV ! drag coefficient for velocity 
    real, dimension(ngrid), intent(out) :: pCdH ! drag coefficient for potential temperature

    integer         :: ig
    real            :: U, z1, Ri, Cd0, Z
    
    real, parameter :: b = 5.0, c = 5.0, d = 5.0, c2b = 2.0 * b, c3bc = 3.0 * b * c, c3b = 3.0 * b
    
    do ig = 1, ngrid
       U = sqrt (pU(ig)**2 + pV(ig)**2 + Emin_turb)
       
       z1  = 1.0 + pZ(ig) / pZ0(ig)

       Cd0 = (Karman/log(z1))**2 * U
       
       ! Gradient Richardson number at surface (Ts = Hs)
       Ri = 0.5 * pG/pTs(ig) * (pH(ig) - pTs(ig)) / U**2 * pZ(ig)
       
       if (Ri < 0.0) then
          z1 = b * Ri / (1.0 + c3bc * Cd0 * sqrt (- z1 * Ri))

          pCdv(ig) = Cd0 * (1.0 - 2.0 * z1)
          pCdh(ig) = Cd0 * (1.0 - 3.0 * z1)
       else
          Z = sqrt (1.0 + d * Ri)

          pCdV(ig) = Cd0 / (1.0 + c2b * Ri / Z)
          pCdH(ig) = Cd0 / (1.0 + c3b * Ri * Z)
       end if
    end do
    pCdV = pCdV * U
    pCdH = pCdH * U
  end subroutine Vdif_Cd

  pure subroutine Vdif_K (ngrid, nlay, pTimestep, pG, pZlev, pZlay, pZ0, pU, pV, pH, pKV, pKH)
    ! Turbulent diffusion coefficients for velocity (pKU) and potential temperature (pKH) using mixing length model
    ! (assumes Pr = 1 so that pKH = pKU: same diffusion coefficients for velocity and potential temperature)
    integer,                       intent(in)  :: ngrid, nlay
    real,                          intent(in)  :: pTimestep, pG
    real, dimension(ngrid),        intent(in)  :: pZ0           ! roughness length
    real, dimension(ngrid,nlay),   intent(in)  :: pZlay
    real, dimension(ngrid,nlay),   intent(in)  :: pU, pV, pH
    real, dimension(ngrid,nlay+1), intent(in)  :: pZlev
    
    real, dimension(ngrid,nlay+1), intent(out) :: pKH, pKV      ! turbulent diffusivities

    integer :: ig, il
    real    :: Lmix, Ri
    real    :: dU, dV, dH, Hlev, N2
    real    :: dVdZ, dZ, z1

    ! Turbulent diffusivities at lower interfaces corresponding to layer il
    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pZlev(ig,il) + pZ0(ig)

          Lmix = Karman * z1 / (1.0 + Karman * z1 / LmixMin) ! mixing length
       
          dU = pU(ig,il) - pU(ig,il-1)
          dV = pV(ig,il) - pV(ig,il-1)
          dH = pH(ig,il) - pH(ig,il-1)
          
          dZ = pZlay(ig,il) - pZlay(ig,il-1)

          dVdZ = sqrt (dU**2 + dV**2) / dZ 

          if (dVdZ > dVdZ_min) then
             Hlev = 0.5 * (pH(ig,il) + pH(ig,il-1)) ! potential temperature

             N2 = pG/Hlev * dH/dZ                   ! Brunt-Vaisala frequency squared
             
             Ri = N2 / dVdZ**2                      ! gradient Richardson number

             ! Turbulent diffusivity
             pkV(ig,il)= Lmix * sqrt (max (Lmix**2 * dVdZ**2 * (1.0 - Ri/Ri_c), Emin_turb))

          else ! minimum diffusivity for low energy cases
             pKV(ig,il) = Lmix * sqrt (Emin_turb)           
          end if
       end do
    end do
    pKH = pKV
  end subroutine Vdif_k
end module turbulence
