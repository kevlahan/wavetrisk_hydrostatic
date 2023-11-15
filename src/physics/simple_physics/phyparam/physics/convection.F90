MODULE convection

#include "use_logging.h"

  IMPLICIT NONE

CONTAINS

  PURE SUBROUTINE sigma_levels(ngrid, nlay, i, pplev, ppopsk, sig, dsig, sdsig)
    !$acc routine seq
    INTEGER, INTENT(IN) :: ngrid, nlay, i
    REAL, INTENT(IN)    :: pplev(ngrid,nlay+1), ppopsk(ngrid,nlay)
    REAL, INTENT(OUT)   :: sig(nlay+1), dsig(nlay), sdsig(nlay)
    !   Calcul des niveaux sigma sur cette colonne
    INTEGER :: l
    DO l=1,nlay+1
       sig(l)=pplev(i,l)/pplev(i,1)
    ENDDO
    DO l=1,nlay
       dsig(l)=sig(l)-sig(l+1)
       sdsig(l)=ppopsk(i,l)*dsig(l)
    ENDDO
  END SUBROUTINE sigma_levels

  SUBROUTINE adjust_onecolumn(ngrid, nlay, i, sig, dsig, sdsig, zu2, zv2, zh2)
    !$acc routine seq
    INTEGER, INTENT(IN) :: ngrid, nlay, i
    REAL, INTENT(IN)    :: sig(nlay+1), sdsig(nlay), dsig(nlay)
    REAL, INTENT(INOUT) :: zu2(ngrid,nlay), zv2(ngrid,nlay), zh2(ngrid,nlay)

    INTEGER :: l,l1,l2
    LOGICAL :: extend
    REAL :: zhm,zsm,zum,zvm,zalpha

    l2=1
    DO WHILE( l2<nlay )
       l2 = l2+1
       ! starting from l1=l2, we shall extend l1..l2 downwards and upwards
       ! until zhm = potential temperature mass-weighted over l1..l2 be
       !     higher than pot. temp. at level l1-1
       !     and lower than pot. temp. at level l2+1
       ! for this we incrementally compute zsm = mass of layers l1 .. l2 and zhm
       zsm = sdsig(l2)
       zhm = zh2(i, l2)
       l1 = l2
       extend = .TRUE.
       DO WHILE(extend)
          ! extend downwards if l1>1 and zhm is lower than layer l1-1
          extend = .FALSE.
          IF (l1 > 1) THEN
             IF (zhm < zh2(i, l1-1)) extend = .TRUE.
          END IF
          IF(extend) THEN
             l1 = l1-1
             zsm = zsm + sdsig(l1)
             zhm = zhm + sdsig(l1) * (zh2(i, l1) - zhm) / zsm
             CYCLE
          END IF
          ! extend upwards if l2<nlay and zhm is higher than layer l2+1
          extend=.FALSE.
          IF (l2 < nlay) THEN
             IF (zh2(i, l2+1) < zhm) extend=.TRUE.
          END IF
          IF(extend) THEN
             l2 = l2+1
             zsm = zsm + sdsig(l2)
             zhm = zhm + sdsig(l2) * (zh2(i, l2) - zhm) / zsm
          END IF
       END DO

       ! at this point the profile l1-1 (l1 ...l2) l2+1 is stable

       IF(l1<l2) THEN
#ifndef _OPENACC
          WRITELOG(*,*) 'Mixing between layers ',l1,l2
          LOG_DBG('convadj')
#endif

          ! perform the mixing of layers l1...l2 :
          ! * potential temperature is fully mixed, ie replaced by mass-weighted average
          ! * momentum u,v is mixed with weight zalpha, i.e. u:=zalpha*mean(u)+(1-zalpha)*u
          ! where zalpha = sum( abs(theta-mean(theta))*dsig) / mean(theta)*sum(dsig)
          ! for large deviations of theta from its mean, this formula could in principe yield zalpha>1.
          zalpha=0.
          zum=0.
          zvm=0.
          DO l = l1, l2
             zalpha=zalpha+ABS(zh2(i,l)-zhm)*dsig(l)
             zh2(i, l) = zhm
             ! we must use compute zum, zvm from zu2, zv2 to conserve momentum (see below)
             zum=zum+dsig(l)*zu2(i,l)
             zvm=zvm+dsig(l)*zv2(i,l)
          END DO
          zalpha=zalpha/(zhm*(sig(l1)-sig(l2+1)))
          zum=zum/(sig(l1)-sig(l2+1))
          zvm=zvm/(sig(l1)-sig(l2+1))

#ifndef _OPENACC
          IF(zalpha.GT.1.) THEN
             WRITELOG(*,*) 'zalpha=',zalpha
             CALL log_gridpoint(i)
             LOG_WARN('convadj')
             STOP
          END IF
#endif
          IF(zalpha.LT.1.e-5) zalpha=1.e-5
          ! zum, zvm are mass-weighted averages of zu2, zv2 => mixing preserves total momentum
          DO l=l1,l2
             zu2(i,l)=zu2(i,l)+zalpha*(zum-zu2(i,l))
             zv2(i,l)=zv2(i,l)+zalpha*(zvm-zv2(i,l))
          ENDDO
       END IF

    END DO

  END SUBROUTINE adjust_onecolumn

  SUBROUTINE convadj(ngrid,nlay,ptimestep, &
       &                   pplay,pplev,ppopsk,   &
       &                  pu,pv,ph,              &
       &                  pdufi,pdvfi,pdhfi,     &
       &                  pduadj,pdvadj,pdhadj)

    !=======================================================================
    !
    !   dry convective adjustment
    !   h = pot. temp. , static instability if ph(above)<ph(below)
    !   tendencies pdfX are first added to profiles pX, X=u,v,h, yieldig pX2
    !   adjustment is performed on profiles pX2
    !   then the difference between adjusted and non-adujsted profiles
    !   is converted back as tendencies
    !
    !=======================================================================

    INTEGER, INTENT(IN) :: ngrid,nlay
    REAL, INTENT(IN) ::  ptimestep
    REAL, INTENT(IN) ::  pplay(ngrid,nlay), pplev(ngrid,nlay+1), ppopsk(ngrid,nlay)
    REAL, INTENT(IN) ::  ph(ngrid,nlay), pu(ngrid,nlay), pv(ngrid,nlay), &
         &               pdhfi(ngrid,nlay), pdufi(ngrid,nlay), pdvfi(ngrid,nlay)
    REAL, INTENT(OUT) :: pdhadj(ngrid,nlay), pduadj(ngrid,nlay), pdvadj(ngrid,nlay)

    !   local:
    REAL sig(nlay+1),sdsig(nlay),dsig(nlay)
    REAL zu(ngrid,nlay), zv(ngrid,nlay), zh(ngrid,nlay)
    REAL zu2(ngrid,nlay), zv2(ngrid,nlay), zh2(ngrid,nlay)
    LOGICAL vtest(ngrid)
    INTEGER jadrs(ngrid)
    INTEGER jcnt, ig, l, jj

    !$acc data copyin(ngrid, nlay)                                           &
    !$acc &    copyin(pplay, pplev, ppopsk, ph, pu, pv, pdhfi, pdufi, pdvfi) &
    !$acc &    copyout(pdhadj, pduadj, pdvadj)                               &
    !$acc &    create(sig, sdsig, dsig, zu, zv, zh, zu2, zv2, zh2, vtest)    &
    !$acc &

    !$acc kernels default(none)

    ! Add physics tendencies to u,v,h
    !$acc loop collapse(2)
    DO l=1,nlay
       DO ig=1,ngrid
          zh(ig,l)=ph(ig,l)+pdhfi(ig,l)*ptimestep
          zu(ig,l)=pu(ig,l)+pdufi(ig,l)*ptimestep
          zv(ig,l)=pv(ig,l)+pdvfi(ig,l)*ptimestep
       END DO
    END DO
    zu2(:,:)=zu(:,:)
    zv2(:,:)=zv(:,:)
    zh2(:,:)=zh(:,:)

    !-----------------------------------------------------------------------
    !   detection des profils a modifier:
    !   ---------------------------------
    !   si le profil est a modifier
    !   (i.e. ph(niv_sup) < ph(niv_inf) )
    !   alors le tableau "vtest" est mis a .TRUE. ;
    !   sinon, il reste a sa valeur initiale (.FALSE.)
    !   cette operation est vectorisable

    vtest(:)=.FALSE.
    !$acc loop seq
    DO l=2,nlay
       !$acc loop
       DO ig=1,ngrid
          IF(zh2(ig,l).LT.zh2(ig,l-1)) vtest(ig)=.TRUE.
       END DO
    END DO

    !----------------------------------
    !  Ajustement des profils instables
    !  --------------------------------

    !$acc loop private(sig,dsig,sdsig)
    DO ig=1,ngrid
       IF(vtest(ig)) THEN
          CALL sigma_levels(ngrid, nlay, ig, pplev, ppopsk, sig, dsig, sdsig)
          CALL adjust_onecolumn(ngrid, nlay, ig, sig, dsig, sdsig, zu2, zv2, zh2)
       ENDIF
    END DO

    !$acc loop collapse(2)
    DO l=1,nlay
       DO ig=1,ngrid
          pdhadj(ig,l)=(zh2(ig,l)-zh(ig,l))/ptimestep
          pduadj(ig,l)=(zu2(ig,l)-zu(ig,l))/ptimestep
          pdvadj(ig,l)=(zv2(ig,l)-zv(ig,l))/ptimestep
       END DO
    END DO

    !$acc end kernels
    !$acc end data

  END SUBROUTINE convadj

END MODULE convection
