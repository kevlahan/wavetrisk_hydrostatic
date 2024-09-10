MODULE turbulence

#include "use_logging.h"

  IMPLICIT NONE
  SAVE
  PRIVATE

  REAL, PARAMETER :: karman=0.4
  REAL :: lmixmin=100., emin_turb=1e-8

  PUBLIC :: vdif, lmixmin, emin_turb

CONTAINS

  PURE SUBROUTINE vdif_cd( ngrid,pz0,pg,pz,pu,pv,pts,ph,pcdv,pcdh)
    !=======================================================================
    !
    !   Subject: computation of the surface drag coefficient using the
    !   -------  approch developed by Loui for ECMWF.
    !
    !   Author: Frederic Hourdin  15 /10 /93
    !   -------
    !=======================================================================
    !   Arguments:
    !   ----------

    INTEGER,INTENT(IN) :: ngrid ! number of atmospheric columns
    REAL,INTENT(IN) :: pg ! gravity (m s -2)
    REAL,INTENT(IN) :: pz0(ngrid) ! surface roughness length scale (m)
    REAL,INTENT(IN) :: pz(ngrid) ! height of the first atmospheric layer (m)
    REAL,INTENT(IN) :: pu(ngrid) ! zonal wind (m/s) in the layer
    REAL,INTENT(IN) :: pv(ngrid) ! meridional wind (m/s) in the layer
    REAL,INTENT(IN) :: pts(ngrid) ! surface temperature (K)
    REAL,INTENT(IN) :: ph(ngrid) ! potential temperature T*(p/ps)^kappa
    REAL,INTENT(OUT) :: pcdv(ngrid) ! Cd for the wind
    REAL,INTENT(OUT) :: pcdh(ngrid) ! Cd for potential temperature

    !   Local:
    !   ------

    REAL, PARAMETER :: b=5., c=5., d=5., umin2=1e-12, &
         c2b=2.*b, c3bc=3.*b*c, c3b=3.*b
    INTEGER ig
    REAL zu2,z1,zri,zcd0,zz

    !-----------------------------------------------------------------------
    !   couche de surface:
    !   ------------------

    !     DO ig=1,ngrid
    !        zu2=pu(ig)*pu(ig)+pv(ig)*pv(ig)+umin2
    !        pcdv(ig)=pz0(ig)*(1.+sqrt(zu2))
    !        pcdh(ig)=pcdv(ig)
    !     ENDDO
    !     RETURN

    !$acc data copyin (pz0, pz, pu, pv, pts, ph)  &
    !$acc &    copyout(pcdv, pcdh)

    !$acc kernels default(none) async
!!!! WARNING, verifier la formule originale de Louis!
    DO ig=1,ngrid
       zu2=pu(ig)*pu(ig)+pv(ig)*pv(ig)+umin2
       zri=pg*pz(ig)*(ph(ig)-pts(ig))/(ph(ig)*zu2)
       z1=1.+pz(ig)/pz0(ig)
       zcd0=karman/log(z1)
       zcd0=zcd0*zcd0*sqrt(zu2)
       IF(zri.LT.0.) THEN
          z1=b*zri/(1.+c3bc*zcd0*sqrt(-z1*zri))
          pcdv(ig)=zcd0*(1.-2.*z1)
          pcdh(ig)=zcd0*(1.-3.*z1)
       ELSE
          zz=sqrt(1.+d*zri)
          pcdv(ig)=zcd0/(1.+c2b*zri/zz)
          pcdh(ig)=zcd0/(1.+c3b*zri*zz)
       ENDIF
    ENDDO
    !$acc end kernels
    ! appelant fait des appels async donc on supprime ce wait. !$acc wait
    !$acc end data

    !-----------------------------------------------------------------------

  END SUBROUTINE vdif_cd

  PURE SUBROUTINE vdif_k(ngrid,nlay,   &
       ptimestep,pg,pzlev,pzlay,pz0,pu,pv,ph,pkv,pkh)
    ! FIXME : pkh := pkv
    INTEGER, INTENT(IN) :: ngrid,nlay
    REAL, INTENT(IN) ::  ptimestep, pg, &
         pzlay(ngrid,nlay),pzlev(ngrid,nlay+1), pz0(ngrid), &
         pu(ngrid,nlay),pv(ngrid,nlay),ph(ngrid,nlay)
    REAL, INTENT(OUT) :: pkv(ngrid,nlay+1),pkh(ngrid,nlay+1)

    INTEGER ig,il
    REAL zdu,zdv,zri,zdvodz2,zdz,z1,lmix

    !$acc data copyin (pzlay, pzlev, pz0, pu, pv, ph)  &
    !$acc &    copyout(pkv, pkh)

    !$acc kernels default(none) async
    DO ig=1,ngrid
       pkv(ig,1)=0.
       pkh(ig,1)=0.
       pkv(ig,nlay+1)=0.
       pkh(ig,nlay+1)=0.
    ENDDO

    DO il=2,nlay
       DO ig=1,ngrid
          z1=pzlev(ig,il)+pz0(ig)
          lmix=karman*z1/(1.+karman*z1/lmixmin)
          !           lmix=lmixmin
          ! WARNING test lmix=lmixmin
          zdu=pu(ig,il)-pu(ig,il-1)
          zdv=pv(ig,il)-pv(ig,il-1)
          zdz=pzlay(ig,il)-pzlay(ig,il-1)
          zdvodz2=(zdu*zdu+zdv*zdv)/(zdz*zdz)
          IF(zdvodz2.LT.1.e-5) THEN
             pkv(ig,il)=lmix*sqrt(emin_turb)
          ELSE
             zri=2.*pg*(ph(ig,il)-ph(ig,il-1)) &
                  / (zdz* (ph(ig,il)+ph(ig,il-1)) *zdvodz2  )
             pkv(ig,il)= &
                  lmix*sqrt(MAX(lmix*lmix*zdvodz2*(1-zri/.4),emin_turb))
          ENDIF
          pkh(ig,il)=pkv(ig,il)
       ENDDO
    ENDDO
    !$acc end kernels
    ! voir vdif_cd !$acc wait
    !$acc end data

  END SUBROUTINE vdif_k

  SUBROUTINE vdif(ngrid,nlay,ptime, &
       ptimestep,pcapcal,pz0, &
       pplay,pplev,pzlay,pzlev, &
       pu,pv,ph,ptsrf,pemis, &
       pdufi,pdvfi,pdhfi,pfluxsrf, &
       pdudif,pdvdif,pdhdif,pdtsrf, &
       lwrite)
    USE phys_const, ONLY: g,r,rcp,cpp

    !=======================================================================
    !
    !   Vertical diffusion
    !   Implicit scheme
    !   We start by adding physical tendencies to variables x
    !   and then solve:
    !      x(t+1) =  x(t) + dt * (dx/dt)phys(t)  +  dt * (dx/dt)difv(t+1)
    !
    !=======================================================================

    INTEGER,INTENT(IN) :: ngrid,nlay
    REAL,INTENT(IN) :: ptime,ptimestep
    REAL,INTENT(IN) :: pcapcal(ngrid),pz0(ngrid)
    REAL,INTENT(IN) :: pplay(ngrid,nlay),pplev(ngrid,nlay+1)
    REAL,INTENT(IN) :: pzlay(ngrid,nlay),pzlev(ngrid,nlay+1)
    REAL,INTENT(IN) :: pu(ngrid,nlay),pv(ngrid,nlay),ph(ngrid,nlay)
    REAL,INTENT(IN) :: ptsrf(ngrid),pemis(ngrid)
    REAL,INTENT(IN) :: pdufi(ngrid,nlay),pdvfi(ngrid,nlay),pdhfi(ngrid,nlay)
    REAL,INTENT(IN) :: pfluxsrf(ngrid)
    REAL,INTENT(OUT) :: pdudif(ngrid,nlay),pdvdif(ngrid,nlay),pdhdif(ngrid,nlay)
    REAL,INTENT(OUT) :: pdtsrf(ngrid)
    LOGICAL,INTENT(IN) :: lwrite
    !
    !   local:
    !   ------

    INTEGER ilev,ig,ilay
    INTEGER unit,ierr,it1,it2
    INTEGER cluvdb,putdat,putvdim,setname,setvdim
    REAL z4st,zdplanck(ngrid),zu2
    REAL zkv(ngrid,nlay+1),zkh(ngrid,nlay+1)
    REAL zcdv(ngrid),zcdh(ngrid)
    REAL zu(ngrid,nlay),zv(ngrid,nlay)
    REAL zh(ngrid,nlay)
    REAL ztsrf2(ngrid)
    REAL z1(ngrid),z2(ngrid)
    REAL za(ngrid,nlay),zb(ngrid,nlay)
    REAL zb0(ngrid,nlay)
    REAL zc(ngrid,nlay),zd(ngrid,nlay)
    REAL zcst1

    !$acc data create (zdplanck, zkv, zkh, zcdv, zcdh, zu, zv, zh, ztsrf2, z1, z2,  za, zb, zb0, zc, zd) &
    !$acc &    copyin (pcapcal, pz0, pplay, pplev, pzlay, pzlev, pu, pv, ph, ptsrf, pemis, pdufi, pdvfi, pdhfi,pfluxsrf) &
    !$acc  &    copyout(pdudif, pdvdif, pdhdif, pdtsrf)
    !
    !-----------------------------------------------------------------------
    !   initializations:
    !   ----------------

    !   computation of rho*dz and dt*rho/dz=dt*rho**2 g/dp:
    !   with rho=p/RT=p/ (R Theta) (p/ps)**kappa
    !   ---------------------------------

    !$acc kernels default(none) async
    DO ilay=1,nlay
       DO ig=1,ngrid
          za(ig,ilay) = (pplev(ig,ilay)-pplev(ig,ilay+1))/g
       ENDDO
    ENDDO

    zcst1=4.*g*ptimestep/(r*r)
    DO ilev=2,nlay
       DO ig=1,ngrid
          zb0(ig,ilev)=pplev(ig,ilev) &
               *(pplev(ig,1)/pplev(ig,ilev))**rcp &
               /(ph(ig,ilev-1)+ph(ig,ilev))
          zb0(ig,ilev)=zcst1*zb0(ig,ilev)*zb0(ig,ilev) &
               / (pplay(ig,ilev-1)-pplay(ig,ilev))
       ENDDO
    ENDDO
    DO ig=1,ngrid
       zb0(ig,1)=ptimestep*pplev(ig,1)/(r*ptsrf(ig))
    ENDDO
    !$acc end kernels

    IF(lwrite) THEN
       !$acc wait
       !$acc update self(pplay, pzlay, pu, pv, ph, za, pplev, pzlev, zb0)
       ig=ngrid/2+1
       WRITELOG(*,*) 'Pression (mbar) altitude (km),u,v,theta, rho dz'
       DO ilay=1,nlay
          WRITELOG(*,*) .01*pplay(ig,ilay),.001*pzlay(ig,ilay), &
               pu(ig,ilay),pv(ig,ilay),ph(ig,ilay),za(ig,ilay)
       ENDDO
       WRITELOG(*,*) 'Pression (mbar) altitude (km),zb'
       DO ilev=1,nlay
          WRITELOG(*,*) .01*pplev(ig,ilev),.001*pzlev(ig,ilev), &
               zb0(ig,ilev)
       ENDDO
       LOG_DBG('vdif')
    ENDIF

    !-----------------------------------------------------------------------
    !   2. Add physical tendencies:
    !   ------------------------------

    !$acc kernels default(none) async
    DO ilev=1,nlay
       DO ig=1,ngrid
          zu(ig,ilev)=pu(ig,ilev)+pdufi(ig,ilev)*ptimestep
          zv(ig,ilev)=pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep
          zh(ig,ilev)=ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep
       ENDDO
    ENDDO
    !$acc end kernels
    ! on sait que vdif_cd et vdif_k sont sur des async, on a donc pas besoin de acc wait !$acc wait

    !-----------------------------------------------------------------------
    !   3. compute  cd :
    !   ----------------
    !
    ! Compute surface drag coefficents zcdv() and zcdh()
    CALL vdif_cd( ngrid,pz0,g,pzlay,pu,pv,ptsrf,ph,zcdv,zcdh)
    ! Compute zkv() and zkh()
    CALL vdif_k(ngrid,nlay,ptimestep,g,pzlev,pzlay,pz0,pu,pv,ph,zkv,zkh)

    !$acc kernels default(none) async
    DO ig=1,ngrid
       zu2=pu(ig,1)*pu(ig,1)+pv(ig,1)*pv(ig,1)
       zcdv(ig)=zcdv(ig)*sqrt(zu2)
       zcdh(ig)=zcdh(ig)*sqrt(zu2)
    ENDDO
    !$acc end kernels


    IF(lwrite) THEN
       !$acc wait
       !$acc update self(zcdv, zcdh,  zkv, zkh)
       WRITELOG(*,*) 'Diagnostique diffusion verticale'
       WRITELOG(*,*) 'LMIXMIN',lmixmin
       WRITELOG(*,*) 'coefficients Cd pour v et h'
       WRITELOG(*,*) zcdv(ngrid/2+1),zcdh(ngrid/2+1)
       WRITELOG(*,*) 'coefficients K pour v et h'
       DO ilev=1,nlay
          WRITELOG(*,*) zkv(ngrid/2+1,ilev),zkh(ngrid/2+1,ilev)
       ENDDO
       LOG_DBG('vdif')
    ENDIF

    !-----------------------------------------------------------------------
    !   vertical integration for u:
    !   -----------------------------

    !$acc kernels default(none) async
    DO ilay=2,nlay
       DO ig=1,ngrid
          zb(ig,ilay)=zkv(ig,ilay)*zb0(ig,ilay)
       ENDDO
    ENDDO

    DO ig=1,ngrid
       zb(ig,1)=zcdv(ig)*zb0(ig,1)
    ENDDO

    DO ig=1,ngrid
       z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
       zc(ig,nlay)=za(ig,nlay)*zu(ig,nlay)*z1(ig)
       zd(ig,nlay)=zb(ig,nlay)*z1(ig)
    ENDDO

    DO ilay=nlay-1,1,-1
       DO ig=1,ngrid
          z1(ig) = 1./(za(ig,ilay)+zb(ig,ilay) &
               + zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
          zc(ig,ilay) = (za(ig,ilay)*zu(ig,ilay) &
               + zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
          zd(ig,ilay) = zb(ig,ilay)*z1(ig)
       ENDDO
    ENDDO

    DO ig=1,ngrid
       zu(ig,1)=zc(ig,1)
    ENDDO
    DO ilay=2,nlay
       DO ig=1,ngrid
          zu(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zu(ig,ilay-1)
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   vertical integration for v:
    !   -----------------------------
    !
    DO ig=1,ngrid
       z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
       zc(ig,nlay)=za(ig,nlay)*zv(ig,nlay)*z1(ig)
       zd(ig,nlay)=zb(ig,nlay)*z1(ig)
    ENDDO

    DO ilay=nlay-1,1,-1
       DO ig=1,ngrid
          z1(ig)=1./(za(ig,ilay)+zb(ig,ilay) &
               + zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
          zc(ig,ilay)=(za(ig,ilay)*zv(ig,ilay) &
               + zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
          zd(ig,ilay)=zb(ig,ilay)*z1(ig)
       ENDDO
    ENDDO

    DO ig=1,ngrid
       zv(ig,1)=zc(ig,1)
    ENDDO
    DO ilay=2,nlay
       DO ig=1,ngrid
          zv(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zv(ig,ilay-1)
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    !   vertical integration for h:
    !   -----------------------------
    !
    DO ilay=2,nlay
       DO ig=1,ngrid
          zb(ig,ilay)=zkh(ig,ilay)*zb0(ig,ilay)
       ENDDO
    ENDDO

    DO ig=1,ngrid
       zb(ig,1)=zcdh(ig)*zb0(ig,1)
    ENDDO

    DO ig=1,ngrid
       z1(ig)=1./(za(ig,nlay)+zb(ig,nlay))
       zc(ig,nlay)=za(ig,nlay)*zh(ig,nlay)*z1(ig)
       zd(ig,nlay)=zb(ig,nlay)*z1(ig)
    ENDDO

    DO ilay=nlay-1,1,-1
       DO ig=1,ngrid
          z1(ig)=1./(za(ig,ilay)+zb(ig,ilay) &
               + zb(ig,ilay+1)*(1.-zd(ig,ilay+1)))
          zc(ig,ilay)=(za(ig,ilay)*zh(ig,ilay) &
               + zb(ig,ilay+1)*zc(ig,ilay+1))*z1(ig)
          zd(ig,ilay)=zb(ig,ilay)*z1(ig)
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   add planck contribution to the implicit scheme:
    !   --------------------------------------------------

    z4st=4.*5.67e-8*ptimestep
    DO ig=1,ngrid
       zdplanck(ig)=z4st*pemis(ig)*ptsrf(ig)*ptsrf(ig)*ptsrf(ig)
    ENDDO
    !$acc end kernels

    !-----------------------------------------------------------------------
    !   compute the evolution of the surface temperature:
    !   -----------------------------------------------

    !$acc kernels default(none) async
    DO ig=1,ngrid
       z1(ig) = pcapcal(ig)*ptsrf(ig)+cpp*zb(ig,1)*zc(ig,1) &
            + zdplanck(ig)*ptsrf(ig)+ pfluxsrf(ig)*ptimestep
       z2(ig) = pcapcal(ig)+cpp*zb(ig,1)*(1.-zd(ig,1))+zdplanck(ig)
       ztsrf2(ig) = z1(ig)/z2(ig)
       zh(ig,1) = zc(ig,1)+zd(ig,1)*ztsrf2(ig)
       pdtsrf(ig) = (ztsrf2(ig)-ptsrf(ig))/ptimestep
    ENDDO
    !$acc end kernels

    !-----------------------------------------------------------------------
    !   final vertical integration:
    !   -----------------------------

    !$acc kernels default(none) async
    DO ilay=2,nlay
       DO ig=1,ngrid
          zh(ig,ilay)=zc(ig,ilay)+zd(ig,ilay)*zh(ig,ilay-1)
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   compute the vertical diffusion tendencies:
    !   -----------------------------------------------------

    DO ilev = 1, nlay
       DO ig=1,ngrid
          pdudif(ig,ilev) = ( zu(ig,ilev) &
               - (pu(ig,ilev)+pdufi(ig,ilev)*ptimestep) )/ptimestep
          pdvdif(ig,ilev) = ( zv(ig,ilev) &
               - (pv(ig,ilev)+pdvfi(ig,ilev)*ptimestep) )/ptimestep
          pdhdif(ig,ilev) = ( zh(ig,ilev) &
               - (ph(ig,ilev)+pdhfi(ig,ilev)*ptimestep) )/ptimestep
       ENDDO
    ENDDO
    !$acc end kernels

    !$acc wait
    IF(lwrite) THEN
       !$acc update self(ptsrf, ztsrf2, ph, zh)
       WRITELOG(*,*)
       WRITELOG(*,*) 'Diagnostique de la diffusion verticale'
       WRITELOG(*,*) 'h avant et apres diffusion verticale'
       WRITELOG(*,*) ptsrf(ngrid/2+1),ztsrf2(ngrid/2+1)
       DO  ilev=1,nlay
          WRITELOG(*,*) ph(ngrid/2+1,ilev),zh(ngrid/2+1,ilev)
       ENDDO
       LOG_DBG('vdif')
    END IF

    !$acc end data
  END SUBROUTINE vdif

END MODULE turbulence
