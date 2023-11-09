MODULE radiative_mod
#include "use_logging.h"
  USE comgeomfi
  USE callkeys, ONLY : diurnal, lverbose
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL, PARAMETER :: pi=2*ASIN(1.), solarcst=1370., stephan=5.67e-08, &
       &             height_scale=10000., ps_rad=1.e5

  PUBLIC :: radiative_tendencies

CONTAINS

  SUBROUTINE radiative_tendencies(lwrite, ngrid, igout, nlayer, &
       gmtime, zdtime, zday, pplev, pplay, pt, &
       pdt, fluxrad)
    USE planet
    USE phys_const,     ONLY : planet_rad, unjours
    USE astronomy,      ONLY : orbite, solarlong
    USE solar,          ONLY : solang, zenang, mucorr
    USE soil_mod,       ONLY : albedo, emissiv, tsurf
    USE radiative_sw,   ONLY : sw
    USE radiative_lw,   ONLY : lw
    USE writefield_mod, ONLY : writefield

    LOGICAL, INTENT(IN)    :: lwrite            ! true if output wanted
    INTEGER, INTENT(IN)    :: ngrid, igout, nlayer
    REAL,    INTENT(IN)    :: gmtime            ! fraction of the day
    REAL,    INTENT(IN)    :: zdtime            ! time step (s)
    REAL,    INTENT(IN)    :: zday              ! elapsed days (and fraction thereof)
    REAL,    INTENT(IN)    :: pplev(ngrid,nlayer+1), pplay(ngrid, nlayer)
    REAL,    INTENT(IN)    :: pt(ngrid, nlayer) ! temperature (K)
    REAL,    INTENT(INOUT) :: pdt(ngrid,nlayer) ! tendency on temperature (K/s)
    REAL,    INTENT(OUT)   :: fluxrad(ngrid)    ! net surface flux

    ! local variables
    REAL :: fract(ngrid),        & ! day fraction
         &  zfluxsw(ngrid),      & ! short-wave flux at surface
         &  zfluxlw(ngrid),      & ! long-wave flux at surface
         &  zdtsw(ngrid,nlayer), & ! short-wave temperature tendency
         &  zdtlw(ngrid,nlayer), & ! long-wave temperature tendency
         &  mu0(ngrid)             ! cosine of zenithal angle

    REAL :: zls, zinsol, tsurf2, ztim1,ztim2,ztim3, dist_sol, declin
    REAL :: zplanck(ngrid)
    INTEGER :: ig, l

    !$acc data copyin(albedo, emissiv, tsurf, pplev, pplay, pt) &
    !$acc      copyout(fluxrad) &
    !$acc      copy(pdt) &
    !$acc      create(fract, zfluxsw, zfluxlw, zdtsw, zdtlw, mu0, zplanck)

    !    2.1 Insolation
    !    --------------------------------------------------

    ! Compute solar longitude zls (rad)
    CALL solarlong(zday,zls)
    CALL orbite(zls,dist_sol,declin)

    ! IF(diurnal) then fract=0 (night) or 1 (day)
    ! ELSE fract = daylight fraction
    ! => incoming = insol*mu*fract

    IF(diurnal) THEN
       IF ( .TRUE. ) then
          ztim1=SIN(declin)
          ztim2=COS(declin)*COS(2.*pi*(zday-.5))
          ztim3=-COS(declin)*SIN(2.*pi*(zday-.5))
          CALL solang(ngrid,sinlon,coslon,sinlat,coslat, &
               &         ztim1,ztim2,ztim3,                   &
               &         mu0,fract)
       ELSE
          CALL zenang(ngrid,zls,gmtime,zdtime,lati,long,mu0,fract)
       ENDIF

       IF(lverbose) THEN
          WRITELOG(*,*) 'day, declin, sinlon,coslon,sinlat,coslat'
          WRITELOG(*,*) zday, declin,                       &
               &        sinlon(igout),coslon(igout),  &
               &        sinlat(igout),coslat(igout)
          LOG_DBG('radiative_tendencies')
       ENDIF
    ELSE
       WRITELOG(*,*) 'declin,ngrid,planet_rad', declin, ngrid, planet_rad
       LOG_DBG('radiative_tendencies')
       CALL mucorr(ngrid,declin,lati,mu0,fract,height_scale,planet_rad)
    ENDIF

    zinsol=solarcst/(dist_sol*dist_sol)

    !    2.2 Radiative tendencies and fluxes:
    !    --------------------------------------------------

    CALL sw(ngrid,nlayer,diurnal,coefvis,albedo, &
         &              pplev,ps_rad,                 &
         &              mu0,fract,zinsol,             &
         &              zfluxsw,zdtsw,                &
         &              lverbose, lwrite)


    CALL lw(ngrid,nlayer,coefir,emissiv, &
         &             pplev,ps_rad,tsurf,pt, &
         &             zfluxlw,zdtlw,         &
         &             lverbose, lwrite)


    !    2.4 surface fluxes
    !    ------------------------------


    !$acc kernels default(none)

    DO ig=1,ngrid
       fluxrad(ig)=emissiv(ig)*zfluxlw(ig) &
            &         +zfluxsw(ig)*(1.-albedo(ig))
       tsurf2 = tsurf(ig)*tsurf(ig)
       zplanck(ig)=emissiv(ig)*stephan*tsurf2*tsurf2
       fluxrad(ig)=fluxrad(ig)-zplanck(ig)
    ENDDO


    !    2.5 Temperature tendencies
    !    --------------------------

    DO l=1,nlayer
       DO ig=1,ngrid
          pdt(ig,l)=pdt(ig,l)+zdtsw(ig,l)+zdtlw(ig,l)
       ENDDO
    ENDDO

    !$acc end kernels

    IF(lverbose) THEN
       !acc update host(albedo, emissiv, mu0, fract, fluxrad, zplanck, pt, pplay, pplev, zdtsw, zdtlw)
       WRITELOG(*,*) 'Diagnostics for radiation'
       WRITELOG(*,*) 'albedo, emissiv, mu0,fract,Frad,Planck'
       WRITELOG(*,*) albedo(igout),emissiv(igout),mu0(igout), &
            &        fract(igout),                       &
            &        fluxrad(igout),zplanck(igout)
       WRITELOG(*,*) 'Tlay Play Plev dT/dt SW dT/dt LW (K/day)'
       WRITELOG(*,*) 'unjours',unjours
       DO l=1,nlayer
          WRITELOG(*,'(3f15.5,2e15.2)') pt(igout,l),     &
               &         pplay(igout,l),pplev(igout,l),  &
               &         zdtsw(igout,l),zdtlw(igout,l)
       ENDDO
       LOG_DBG('radiative_tendencies')
    ENDIF

    IF(lwrite) THEN
       !acc update host(mu0, zfluxsw, zfluxlw, zdtsw, zdtlw)
       call writefield('mu0','Cosine zenithal angle','',mu0)
       call writefield('swsurf','SW surf','W/m2',zfluxsw)
       call writefield('lwsurf','LW surf','W/m2',zfluxlw)
       call writefield('dtsw','dtsw',' ',zdtsw)
       call writefield('dtlw','dtlw',' ',zdtlw)
    END IF

    !$acc end data

  END SUBROUTINE radiative_tendencies

END MODULE radiative_mod
