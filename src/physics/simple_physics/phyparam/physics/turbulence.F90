module turbulence
#include "use_logging.h"
  implicit none
  save
  private
  real, parameter :: karman = 0.4
  real            :: lmixmin = 100.0, emin_turb = 1e-8
  pUblic          :: Vdif, lmixmin, emin_turb
contains
  subroutine Vdif (ngrid, nlay, ptime, ptimestep, &
       pcapcal, pz0, pplay, pplev, pzlay, pzlev, &
       pU, pV, pH, &
       ptsrf, pemis, &
       pdUfi, pdVfi, pdHfi, &
       pfluxsrf, &
       pdUdif, pdVdif, pdhdif, &
       pdtsrf, &
       lwrite, tendency)
    use pHys_const, only: g,r,rcp,cpp
    !=======================================================================
    !
    !   Vertical diffusion for vertical turbulent diffusion
    !
    !
    !   PHysical are advanced using Forward Euler (first-order explicit) and
    !   vertical diffusion is advanced using Backward Euler (first-order implicit)
    !
    !      x(t+dt) = x(t) + dt * (dx/dt) pHys(t)  +  dt * (dx/dt) difv(t+dt)
    !
    !   If flag tendency = true, returns a pseudo "tendency" T(t) :
    !
    !   T(t) = ( -x(t) + dt * (dx/dt) pHys(t)  +  dt * (dx/dt) difv(t+dt) ) / dt
    !
    !   that can be used in Forward Euler:
    !
    !   x(t+dt) = x(t) + dt T(t)
    !
    !   If flag tendency = false, returns the solutions at the new time t+dt.
    !
    !=======================================================================
    ! InpUt variables
    integer,                       intent(in) :: ngrid                   ! number of columns
    integer,                       intent(in) :: nlay                    ! number of vertical layers

    real,                          intent(in) :: ptime                   ! current time t
    real,                          intent(in) :: ptimestep               ! time step dt

    ! Surface variables
    real, dimension(ngrid),        intent(in) :: pcapcal                 ! surface conduction flux
    real, dimension(ngrid),        intent(in) :: pemis                   ! emissivity of surface
    real, dimension(ngrid),        intent(in) :: pfluxsrf                ! surface flux
    real, dimension(ngrid),        intent(in) :: pz0                     ! z coordinate at surface [m]
    real, dimension(ngrid),        intent(in) :: ptsrf                   ! surface temperature [K]

    ! Layer centre variables
    real, dimension(ngrid,nlay),   intent(in) :: pU, pV, pH              ! prognostic variables at time t
    real, dimension(ngrid,nlay),   intent(in) :: pdUfi, pdVfi, pdHfi     ! pHysical tendencies at time t
    real, dimension(ngrid,nlay),   intent(in) :: pplay                   ! pressure at middle of layers [Pa]
    real, dimension(ngrid,nlay),   intent(in) :: pzlay                   ! z coordinate at middle of layers [m]

    ! Layer interface variables (interface pressure and z coordinate)
    real, dimension(ngrid,nlay+1), intent(in) :: pplev                   ! pressure at interfaces between layers [Pa]
    real, dimension(ngrid,nlay+1), intent(in) :: pzlev                   ! pressure at middle of layers [Pa]

    logical,                       intent(in) :: lwrite                  ! write diagnostics flag
    logical,                       intent(in) :: tendency                ! returns pseudo-tendencies if T, solutions if F

    ! OutpUt variables
    real, dimension(ngrid),        intent(out) :: pdTsrf                 ! new surface temperature at time t+dt
    real, dimension(ngrid,nlay),   intent(out) :: pdHdif, pdUdif, pdVdif ! pseudo-tendencies for Forward Euler

    ! Local variables
    integer                       :: ilev, ig, ilay
    integer                       :: unit, ierr, it1, it2
    integer                       :: cluvdb, pUtdat, pUtvdim, setname, setvdim

    real                          :: z4st, zcst1
    real, dimension(ngrid)        :: z1, z2, zdplanck, zCdh, zCdv, zTsrf2, zu2
    real, dimension(ngrid,nlay)   :: za, zb, zb0, zc, zd, zh, zu, zv
    real, dimension(ngrid,nlay+1) :: zKv, zKh

    !------------------------------------------------------------------------------------------------------
    !   Initializations:
    !------------------------------------------------------------------------------------------------------
    !   CompUtation of rho*dz and dt*rho/dz=dt*rho**2 g/dp with rho = p/RT  =p/ (R Theta) (p/ps)**kappa
    !------------------------------------------------------------------------------------------------------
    za(:,1:nlay) = (pplev(:,1:nlay) - pplev(:,2:nlay+1)) / g

    zcst1 = 4.0 * g * ptimestep / r**2
    do ilev = 2, nlay
       zb0(:,ilev) = pplev(:,ilev) * (pplev(:,1) / pplev(:,ilev))**rcp / (pH(:,ilev-1) + pH(:,ilev))
       zb0(:,ilev) = zcst1 * zb0(:,ilev) * zb0(:,ilev) / (pplay(:,ilev-1) - pplay(:,ilev))
    end do
    zb0(:,1) = ptimestep * pplev(:,1) / (r*pTsrf)

    if (lwrite) then
       ig = ngrid/2 + 1
       WRITELOG (*,*) 'Pression (mbar) altitude (km),u,v,theta, rho dz'
       do ilay = 1, nlay
          WRITELOG (*,*) 0.01 * pplay(ig,ilay), 0.001 * pzlay(ig,ilay), pU(ig,ilay), pV(ig,ilay), pH(ig,ilay), za(ig,ilay)
       end do
       WRITELOG (*,*) 'Pression (mbar) altitude (km),zb'
       do ilev = 1, nlay
          WRITELOG (*,*) 0.01 * pplev(ig,ilev), 0.001 * pzlev(ig,ilev), zb0(ig,ilev)
       end do
       LOG_DBG('Vdif')
    endif

    !-----------------------------------------------------------------------
    !   2. Add pHysical tendencies (Forward Euler step)
    !-----------------------------------------------------------------------
    zu = pU + pdUfi * ptimestep
    zv = pV + pdVfi * ptimestep
    zh = pH + pdHfi * ptimestep

    !-----------------------------------------------------------------------
    !   3. CompUte  drag coefficients Cd
    !-----------------------------------------------------------------------
    ! CompUte surface drag coefficents zCdv() and zCdh()
    CALL Vdif_Cd (ngrid, pz0, g, pzlay, pU, pV, pTsrf, pH, zCdv, zCdh)

    ! CompUte zKv and zKh
    CALL Vdif_k (ngrid, nlay, ptimestep, g, pzlev, pzlay, pz0, pU, pV, pH, zKv, zKh)

    zu2 = pU(:,1)**2 + pV(:,1)**2
    zCdv = zCdv * sqrt (zu2)
    zCdh = zCdh * sqrt (zu2)

    if (lwrite) then
       WRITELOG (*,*) 'Diagnostique diffusion verticale'
       WRITELOG (*,*) 'LMIXMIN',lmixmin

       WRITELOG (*,*) 'Coefficients Cd pour v et h'
       WRITELOG (*,*) zCdv(ngrid/2+1), zCdh(ngrid/2+1)

       WRITELOG (*,*) 'Coefficients K pour v et h'
       do ilev = 1, nlay
          WRITELOG (*,*) zKv(ngrid/2+1,ilev), zKh(ngrid/2+1,ilev)
       end do
       LOG_DBG ('Vdif')
    endif

    !-----------------------------------------------------------------------
    !   Vertical integration for u
    !-----------------------------------------------------------------------
    zb(:,2:nlay) = zKv(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1) = zCdv * zb0(:,1)

    z1 = 1.0/ (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zu(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zu(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    zu(:,1) = zc(:,1)
    zu(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zu(:,1:nlay-1)

    !-----------------------------------------------------------------------
    !   Vertical integration for v
    !-----------------------------------------------------------------------
    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zv(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zv(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    zv(:,1)      = zc(:,1)
    zv(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zv(:,1:nlay-1)

    !------------------------------------------------------------------------
    !   Vertical integration for h
    !------------------------------------------------------------------------
    zb(:,2:nlay) = zKh(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1)      = zCdh * zb0(:,1)

    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zh(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zh(:,ilay)+ zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    !------------------------------------------------------------------------
    !   Add Planck contribution to the implicit scheme
    !------------------------------------------------------------------------
    z4st = 4.0 * 5.67e-8 * ptimestep
    zdplanck = z4st * pemis * pTsrf**3

    !-----------------------------------------------------------------------
    !   CompUte the evolution of the surface temperature
    !-----------------------------------------------------------------------
    z1 = pcapcal * pTsrf + cpp * zb(:,1) * zc(:,1)         + zdplanck * pTsrf + pfluxsrf * ptimestep
    z2 = pcapcal         + cpp * zb(:,1) * (1.0 - zd(:,1)) + zdplanck

    zTsrf2 = z1 / z2
    zh(:,1) = zc(:,1) + zd(:,1) * zTsrf2
    pdTsrf = (zTsrf2 - pTsrf) / ptimestep

    !-----------------------------------------------------------------------
    !   Final vertical integration
    !-----------------------------------------------------------------------
    zh(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zh(:,1:nlay-1)

    !-----------------------------------------------------------------------
    !   CompUte the vertical diffusion tendencies:
    !   "tendency" (-u(t) + u(t+dt))/dt to be used as a tendency
    !   in a Forward Euler step
    !-----------------------------------------------------------------------

    if (tendency) then ! pseudo tendencies
       pdUdif = (-pU + zu - pdUfi * ptimestep) / ptimestep
       pdVdif = (-pV + zv - pdVfi * ptimestep) / ptimestep
       pdHdif = (-pH + zh - pdHfi * ptimestep) / ptimestep
    else               ! new solution at t+dt
       pdUdif = zu - pdUfi * ptimestep
       pdVdif = zv - pdVfi * ptimestep
       pdHdif = zh - pdHfi * ptimestep
    end if

    if (lwrite) then
       WRITELOG (*,*)
       WRITELOG (*,*) 'Diagnostique de la diffusion verticale'
       WRITELOG (*,*) 'h avant et apres diffusion verticale'
       WRITELOG (*,*) pTsrf(ngrid/2+1), zTsrf2(ngrid/2+1)
       do ilev = 1, nlay
          WRITELOG (*,*) pH(ngrid/2+1,ilev), zh(ngrid/2+1,ilev)
       end do
       LOG_DBG ('Vdif')
    end if
  end subroutine Vdif

  pUre subroutine Vdif_Cd (ngrid, pz0, pg, pz, pU, pV, pts, pH, pCdv, pCdh)
    !=======================================================================
    !
    !   Subject: compUtation of the surface drag coefficient using the
    !   -------  approch developed by Louis for ECMWF.
    !
    !   Author: Frederic Hourdin  15 /10 /93
    !   -------
    !=======================================================================
    integer,                intent(in) :: ngrid ! number of atmospHeric columns

    real,                   intent(in) :: pg    ! gravity [m/s^2]

    real, dimension(ngrid), intent(in) :: pz0   ! surface roughness length scale [m]
    real, dimension(ngrid), intent(in) :: pz    ! height of the first atmospHeric layer [m]
    real, dimension(ngrid), intent(in) :: pU    ! zonal velocity (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pV    ! meridional velocity (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pTs   ! surface temperature [K]
    real, dimension(ngrid), intent(in) :: pH    ! potential temperature T * (P/Ps)^kappa

    real, dimension(ngrid), intent(out) :: pCdv ! Cd for the wind
    real, dimension(ngrid), intent(out) :: pCdh ! Cd for potential temperature

    integer         :: ig
    real            :: zu2, z1, zri, zCd0, zz
    real, parameter :: b = 5.0, c = 5.0, d = 5.0, umin2 = 1e-12, c2b = 2.0 * b, c3bc = 3.0 * b * c, c3b = 3.0 * b

    !-----------------------------------------------------------------------
    !   couche de surface:
    !-----------------------------------------------------------------------
    !        zu2 = pU**2 + pV**2 + umin2
    !        pCdv = pz0 * (1.0 + sqrt (zu2))
    !        pCdh = pCdv
!!!! WARNING, verifier la formule originale de Louis!
    do ig = 1, ngrid
       zu2 = pU(ig)**2 +  pV(ig)**2 + umin2

       zri = pg  * pz(ig) * (pH(ig) - pTs(ig)) / (pH(ig) * zu2)
       z1  = 1.0 + pz(ig) / pz0(ig)

       zCd0 = karman / log (z1)
       zCd0 = zCd0**2 * sqrt (zu2)

       if (zri < 0.0) then
          z1 = b * zri / (1.0 + c3bc * zCd0 * sqrt (- z1 * zri))

          pCdv(ig) = zCd0 * (1.0 - 2.0 * z1)
          pCdh(ig) = zCd0 * (1.0 - 3.0 * z1)
       else
          zz = sqrt (1.0 + d * zri)

          pCdv(ig) = zCd0 / (1.0 + c2b * zri / zz)
          pCdh(ig) = zCd0 / (1.0 + c3b * zri * zz)
       end if
    end do
  end subroutine Vdif_Cd

  pUre subroutine Vdif_k (ngrid, nlay, ptimestep, pg, pzlev, pzlay, pz0, pU, pV, pH, pKv, pKh)
    ! FIXME : pKh := pKv
    integer,                       intent(in)  :: ngrid, nlay
    real,                          intent(in)  :: ptimestep, pg
    real, dimension(ngrid),        intent(in)  :: pz0
    real, dimension(ngrid,nlay),   intent(in)  :: pH, pU, pV, pzlay
    real, dimension(ngrid,nlay+1), intent(in)  :: pzlev
    real, dimension(ngrid,nlay+1), intent(out) :: pKh, pKv

    integer :: ig, il
    real    :: lmix, zdu, zdv, zri, zdvodz2, zdz, z1

    pKv(:,1)      = 0.0
    pKh(:,1)      = 0.0
    pKv(:,nlay+1) = 0.0
    pKh(:,nlay+1) = 0.0

    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pzlev(ig,il) + pz0(ig)
          lmix = karman * z1 / (1.0 + karman * z1 / lmixmin)
          !           lmix = lmixmin
          ! WARNING test lmix = lmixmin
          zdu = pU(ig,il)    -    pU(ig,il-1)
          zdv = pV(ig,il)    -    pV(ig,il-1)
          zdz = pzlay(ig,il) - pzlay(ig,il-1)

          zdvodz2 = (zdu**2 + zdv**2) / zdz**2
          if (zdvodz2 < 1.0e-5) then
             pKv(ig,il) = lmix * sqrt (emin_turb)
          else
             zri = 2.0 * pg * (pH(ig,il) - pH(ig,il-1)) / (zdz * (pH(ig,il) + pH(ig,il-1)) * zdvodz2)

             pKv(ig,il)= lmix * sqrt (max (lmix * lmix * zdvodz2 * (1.0 - zri/0.4), emin_turb))
          end if
          pKh(ig,il) = pKv(ig,il)
       end do
    end do
  end subroutine Vdif_k
end module turbulence
