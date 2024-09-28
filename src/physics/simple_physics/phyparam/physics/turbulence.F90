module turbulence
#include "use_logging.h"
  implicit none
  save
  private
  real, parameter :: Karman = 0.4
  real            :: LmixMin = 100.0, Emin_turb = 1e-8
  public          :: Vdif, LmixMin, Emin_turb
contains
  subroutine Vdif (ngrid, nlay, ptime, ptimestep, &
       pcapcal, pz0, pplay, pplev, pZlay, pzlev, &
       pU, pV, pH, &
       ptsrf, pemis, &
       pdUfi, pdVfi, pdHfi, &
       pFluxSrf, &
       pdUdif, pdVdif, pdHdif, pdTSrf, &
       lwrite)
    use phys_const, only: g,r,rcp,Cpp
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
    real, dimension(ngrid),        intent(in) :: pz0                     ! z coordinate at surface [m]
    real, dimension(ngrid),        intent(in) :: pTsrf                   ! surface temperature [K]

    ! Layer centre variables
    real, dimension(ngrid,nlay),   intent(in) :: pU, pV, pH              ! prognostic variables at time t from physics
    real, dimension(ngrid,nlay),   intent(in) :: pdUfi, pdVfi, pdHfi     ! physics tendencies at time t (explicit scheme)
    real, dimension(ngrid,nlay),   intent(in) :: pPlay                   ! pressure at middle of layers [Pa]
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
    integer                       :: ilev, ig, ilay
    integer                       :: unit, ierr, it1, it2
    integer                       :: cluvdb, pUtdat, pUtvdim, setname, setvdim

    real                          :: z4st, zcst1
    real, dimension(ngrid)        :: z1, z2, zdPlanck, zcdh, zcdv, zTsrf2, zU2
    real, dimension(ngrid,nlay)   :: zH, zU, zV                         ! intermediate prognostic variables at t+dt
    real, dimension(ngrid,nlay)   :: za, zb, zb0, zc, zd
    real, dimension(ngrid,nlay+1) :: zkV, zkH

    !------------------------------------------------------------------------------------------------------
    !   Initializations:
    !
    !   Computation of rho*dz and dt*rho/dz=dt*rho**2 g/dp with rho = p/RT  = p/(R Theta) (p/ps)**kappa
    !
    !------------------------------------------------------------------------------------------------------
    za(:,1:nlay) = (pPlev(:,1:nlay) - pPlev(:,2:nlay+1)) / g

    zcst1 = 4.0 * g * pTimestep / r**2
    do ilev = 2, nlay
       zb0(:,ilev) = pPlev(:,ilev) * (pPlev(:,1) / pPlev(:,ilev))**rcp / (pH(:,ilev-1) + pH(:,ilev))
       zb0(:,ilev) = zcst1 * zb0(:,ilev) * zb0(:,ilev) / (pPlay(:,ilev-1) - pPlay(:,ilev))
    end do
    zb0(:,1) = pTimestep * pPlev(:,1) / (r*pTsrf)

    if (lwrite) then
       ig = ngrid/2 + 1
       WRITELOG (*,*) 'Pression (mbar) altitude (km),u,v,theta, rho dz'
       do ilay = 1, nlay
          WRITELOG (*,*) 0.01 * pPlay(ig,ilay), 0.001 * pZlay(ig,ilay), pU(ig,ilay), pV(ig,ilay), pH(ig,ilay), za(ig,ilay)
       end do
       WRITELOG (*,*) 'Pression (mbar) altitude (km),zb'
       do ilev = 1, nlay
          WRITELOG (*,*) 0.01 * pPlev(ig,ilev), 0.001 * pZlev(ig,ilev), zb0(ig,ilev)
       end do
       LOG_DBG('Vdif')
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
    ! Surface drag coefficents zcdv() and zcdh()
    CALL Vdif_Cd (ngrid, pz0, g, pZlay, pU, pV, pTsrf, pH, zcdv, zcdh)

    ! zkV and zkH
    CALL Vdif_k (ngrid, nlay, pTimestep, g, pZlev, pZlay, pz0, pU, pV, pH, zkV, zkH)

    zU2 = pU(:,1)**2 + pV(:,1)**2
    zcdv = zcdv * sqrt (zU2)
    zcdh = zcdh * sqrt (zU2)

    if (lwrite) then
       WRITELOG (*,*) 'Diagnostique diffusion verticale'
       WRITELOG (*,*) 'LmixXmin',LmixMin

       WRITELOG (*,*) 'Coefficients Cd pour V et H'
       WRITELOG (*,*) zcdv(ngrid/2+1), zcdh(ngrid/2+1)

       WRITELOG (*,*) 'Coefficients K pour V et H'
       do ilev = 1, nlay
          WRITELOG (*,*) zkV(ngrid/2+1,ilev), zkH(ngrid/2+1,ilev)
       end do
       LOG_DBG ('Vdif')
    endif

    !-----------------------------------------------------------------------
    !   Vertical integration for U
    !-----------------------------------------------------------------------
    zb(:,2:nlay) = zkV(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1) = zcdv * zb0(:,1)

    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = zU(:,nlay) * za(:,nlay) * z1
    zd(:,nlay) =              zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))

       zc(:,ilay) = (za(:,ilay) * zU(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    ! Diffusion tendency for zonal velocity
    zU(:,1) = zc(:,1)
    do ilay = 2, nlay
       zU(:,ilay) = zc(:,ilay) + zd(:,ilay) * zU(:,ilay-1)
    end do

    !-----------------------------------------------------------------------
    !   Vertical integration for V
    !-----------------------------------------------------------------------
    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = zV(:,nlay) * za(:,nlay) * z1
    zd(:,nlay) =              zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))

       zc(:,ilay) = (za(:,ilay) * zV(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    ! Diffusion tendency for meridional velocity
    zV(:,1) = zc(:,1)
    do ilay = 2, nlay
       zV(:,ilay) = zc(:,ilay) + zd(:,ilay) * zV(:,ilay-1)
    end do

    !------------------------------------------------------------------------
    !   Vertical integration for potential temperature H
    !------------------------------------------------------------------------
    zb(:,2:nlay) = zkH(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1)      = zcdh * zb0(:,1)

    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))

    zc(:,nlay) = za(:,nlay) * zH(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))

       zc(:,ilay) = (za(:,ilay) * zH(:,ilay)+ zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    !------------------------------------------------------------------------
    !   Add Planck contribution to implicit scheme
    !------------------------------------------------------------------------
    z4st = 4.0 * 5.67e-8 * pTimestep
    zdPlanck = z4st * pEmis * pTsrf**3

    !-----------------------------------------------------------------------
    !   Compute evolution of surface temperature
    !-----------------------------------------------------------------------
    z1 = pCapCal * pTsrf + Cpp * zb(:,1) * zc(:,1)         + zdPlanck * pTsrf + pFluxSrf * pTimestep
    z2 = pCapCal         + Cpp * zb(:,1) * (1.0 - zd(:,1)) + zdPlanck

    zTsrf2 = z1 / z2

    ! Diffusion tendency for potential temperature
    zH(:,1) = zc(:,1) + zd(:,1) * zTsrf2
    do ilay = 2, nlay
       zH(:,ilay) = zc(:,ilay) + zd(:,ilay) * zH(:,ilay-1)
    end do

    !-----------------------------------------------------------------------
    !   Complete pseudo-tendencies/solution at t+dt
    !-----------------------------------------------------------------------
    pdUdif = (-pU    + zU - pdUfi * pTimestep) / pTimestep
    pdVdif = (-pV    + zV - pdVfi * pTimestep) / pTimestep
    pdHdif = (-pH    + zH - pdHfi * pTimestep) / pTimestep
    pdTsrf = (-pTsrf + zTsrf2 )                / pTimestep

    if (lwrite) then
       WRITELOG (*,*)
       WRITELOG (*,*) 'Diagnostique de la diffusion verticale'
       WRITELOG (*,*) 'Potential temperature avant et apres diffusion verticale'
       WRITELOG (*,*) pTsrf(ngrid/2+1), zTsrf2(ngrid/2+1)
       do ilev = 1, nlay
          WRITELOG (*,*) pH(ngrid/2+1,ilev), zH(ngrid/2+1,ilev)
       end do
       LOG_DBG ('Vdif')
    end if
  end subroutine Vdif

  pUre subroutine Vdif_Cd (ngrid, pz0, pg, pz, pU, pV, pts, pH, pCdv, pCdh)
    !=======================================================================
    !
    !   Subject: Computation of  surface drag coefficient using
    !   -------  approch developed by Louis for ECMWF.
    !
    !   Author: Frederic Hourdin  15 /10 /93
    !   -------
    !=======================================================================
    integer,                intent(in) :: ngrid ! number of atmospHeric columns

    real,                   intent(in) :: pg    ! gravity [m/s^2]

    real, dimension(ngrid), intent(in) :: pZ0   ! surface roughness length scale [m]
    real, dimension(ngrid), intent(in) :: pZ    ! height of the first atmospHeric layer [m]
    real, dimension(ngrid), intent(in) :: pU    ! zonal velocity (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pV    ! meridional velocity (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pTs   ! surface temperature [K]
    real, dimension(ngrid), intent(in) :: pH    ! potential temperature T * (P/Ps)^kappa

    real, dimension(ngrid), intent(out) :: pCdv ! Cd for the wind
    real, dimension(ngrid), intent(out) :: pCdh ! Cd for potential temperature

    integer         :: ig
    real            :: zU2, z1, zRi, zCd0, zz
    real, parameter :: b = 5.0, c = 5.0, d = 5.0, umin2 = 1e-12, c2b = 2.0 * b, c3bc = 3.0 * b * c, c3b = 3.0 * b

    !-----------------------------------------------------------------------
    !   couche de surface:
    !-----------------------------------------------------------------------
    !        zU2 = pU**2 + pV**2 + umin2
    !        pCdv = pZ0 * (1.0 + sqrt (zU2))
    !        pCdh = pCdv
!!!! WARNING, verifier la formule originale de Louis!
    do ig = 1, ngrid
       zU2 = pU(ig)**2 +  pV(ig)**2 + umin2

       zRi = pg  * pZ(ig) * (pH(ig) - pTs(ig)) / (pH(ig) * zU2)
       z1  = 1.0 + pZ(ig) / pZ0(ig)

       zCd0 = Karman / log (z1)
       zCd0 = zCd0**2 * sqrt (zU2)

       if (zRi < 0.0) then
          z1 = b * zRi / (1.0 + c3bc * zCd0 * sqrt (- z1 * zRi))

          pCdv(ig) = zCd0 * (1.0 - 2.0 * z1)
          pCdh(ig) = zCd0 * (1.0 - 3.0 * z1)
       else
          zz = sqrt (1.0 + d * zRi)

          pCdv(ig) = zCd0 / (1.0 + c2b * zRi / zz)
          pCdh(ig) = zCd0 / (1.0 + c3b * zRi * zz)
       end if
    end do
  end subroutine Vdif_Cd

  pUre subroutine Vdif_k (ngrid, nlay, pTimestep, pg, pZlev, pZlay, pZ0, pU, pV, pH, pkV, pKh)
    ! FIXME : pkH := pkV
    integer,                       intent(in)  :: ngrid, nlay
    real,                          intent(in)  :: pTimestep, pg
    real, dimension(ngrid),        intent(in)  :: pZ0
    real, dimension(ngrid,nlay),   intent(in)  :: pZlay
    real, dimension(ngrid,nlay),   intent(in)  :: pH, pU, pV
    real, dimension(ngrid,nlay+1), intent(in)  :: pZlev
    real, dimension(ngrid,nlay+1), intent(out) :: pkH, pkV

    integer :: ig, il
    real    :: Lmix, zdU, zdV, zRi, zdVodz2, zdz, z1

    pkV = 0.0
    pkH = 0.0

    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pZlev(ig,il) + pZ0(ig)
          Lmix = Karman * z1 / (1.0 + Karman * z1 / LmixMin)
          !              Lmix = LmixMin
          ! WARNING test Lmix = LmixMin
          zdU = pU(ig,il)    -    pU(ig,il-1)
          zdV = pV(ig,il)    -    pV(ig,il-1)
          zdz = pZlay(ig,il) - pZlay(ig,il-1)

          zdVodz2 = (zdU**2 + zdV**2) / zdz**2
          if (zdVodz2 < 1.0e-5) then
             pkV(ig,il) = Lmix * sqrt (Emin_turb)
          else
             zRi = 2.0 * pg * (pH(ig,il) - pH(ig,il-1)) / (zdz * (pH(ig,il) + pH(ig,il-1)) * zdVodz2)

             pkV(ig,il)= Lmix * sqrt (max (Lmix**2 * zdVodz2 * (1.0 - zRi/0.4), Emin_turb))
          end if
          pkH(ig,il) = pkV(ig,il)
       end do
    end do
  end subroutine Vdif_k
end module turbulence
