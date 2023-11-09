MODULE mellor_yamada

  IMPLICIT NONE

CONTAINS

  PURE SUBROUTINE my_25(imax,kmax,dt,gp,zi,z,u,v,teta,cd,q2,long,km,kh)

    !**********************************************************************
    !****** FERMETURE MELLOR-YAMADA DE NIVEAU 2.5 (QUASI-EQUILIBRE) *******
    !* q2 au interfaces entre mailles.
    !**********************************************************************

    INTEGER, INTENT(IN) :: imax,kmax
    REAL, INTENT(IN) :: dt, gp
    REAL, INTENT(IN) :: z(imax,kmax), &
         u(imax,kmax), v(imax,kmax), teta(imax,kmax), cd(imax)
    REAL, INTENT(INOUT) :: zi(imax,kmax+1), q2(imax,kmax+1)
    REAL, INTENT(OUT) :: long(imax,kmax+1), km(imax,kmax+1), kh(imax,kmax+1)

    !************* DECLARATIONS *******************************************

    INTEGER, PARAMETER :: impmax=5
    REAL, PARAMETER :: kappa=0.4, &
         a1=0.92, a2=0.74, b1=16.6, b2=10.1, c1=0.08,           &
         e1=1.8, e2=1.33, &
         khmin=1.0e-5, kmmin=1.0e-5, kqmin=1.0e-5,              &
         q2min=0.001, q2lmin=0.001, &
         ghmax=0.023, ghmin=-0.28

    REAL longc, au, ateta, az, adz, akq, acd, &
         adzi, aq2, al, akm, akh, am2, al0, &
         am, azi, bet, beta, dest, du, dv, gh, &
         prod, q2c, sh, sm, sq, us, vt, vt1, vt2

    REAL unsdz(imax,kmax),unsdzi(imax,kmax+1)
    REAL kq(imax,kmax),                                               &
         &     m2(imax,kmax+1),n2(imax,kmax+1),ri(imax,kmax+1)
    REAL a(imax,kmax),b(imax,kmax),c(imax,kmax),f(imax,kmax),         &
         &     alph(imax,kmax)
    REAL ksdz2inf,ksdz2sup

    INTEGER :: i,k

    !************* INCREMENTS VERTICAUX ***********************************

    DO i=1,imax
       zi(i,kmax+1)=zi(i,kmax)+2.0*(z(i,kmax)-zi(i,kmax))
    END DO
    DO k=1,kmax
       DO  i=1,imax
          unsdz(i,k)=1.0/(zi(i,k+1)-zi(i,k))
       END DO
    END DO

    DO k=2,kmax
       DO i=1,imax
          unsdzi(i,k)=1.0/(z(i,k)-z(i,k-1))
       END do
    END DO

    DO i=1,imax
       unsdzi(i,1)=0.5/(z(i,1)-zi(i,1))
       unsdzi(i,kmax+1)=0.5/(zi(i,kmax+1)-z(i,kmax))
    END do

    !**********************************************************************


    !************* DIFFUSIVITES KH, KM et KQ ******************************
    ! Ci-dessous, une premiere estimation des diffusivites turbulentes km *
    ! et kh est effectuee pour utilisation dans les taux de production    *
    ! et de destruction de q2 et q2l. On calcule aussi kq.                *

    DO k=2,kmax
       DO i=1,imax
          beta=2.0/(teta(i,k)+teta(i,k-1))
          n2(i,k)=beta*gp*unsdzi(i,k)*(teta(i,k)-teta(i,k-1))
          n2(i,k)=amax1(0.0,n2(i,k))
          du=unsdzi(i,k)*(u(i,k)-u(i,k-1))
          dv=unsdzi(i,k)*(v(i,k)-v(i,k-1))
          m2(i,k)=du*du+dv*dv
          ri(i,k)=n2(i,k)/(m2(i,k)+1.0e-10)
          ri(i,k)=amax1(-0.1,min(4.0,ri(i,k)))
       END DO
    END DO

    DO k=2,kmax
       DO i=1,imax
          vt=kappa*(zi(i,k)-zi(i,1))
          long(i,k)=vt/(1.0+vt/160.0)
          IF(n2(i,k).gt.0.0) THEN
             long(i,k)=min(0.53*sqrt(q2(i,k))/sqrt(n2(i,k)),long(i,k))
          END IF
          gh=amax1(ghmin,                                                 &
               &     min(ghmax,-long(i,k)*long(i,k)*n2(i,k)/q2(i,k)))
          sm=a1*(1.0-3.0*c1-6.0*a1/b1-3.0*a2*gh*                          &
               &           ((b2-3.0*a2)*(1.0-6.0*a1/b1)-3.0*c1*(b2+6.0*a1)))/     &
               &        ((1.0-3.0*a2*gh*(6.0*a1+b2))*(1.0-9.0*a1*a2*gh))
          km(i,k)=amax1(kmmin,long(i,k)*sqrt(q2(i,k))*sm)
          sh=a2*(1.0-6.0*a1/b1)/(1.0-3.0*a2*gh*(6.0*a1+b2))
          kh(i,k)=amax1(khmin,long(i,k)*sqrt(q2(i,k))*sh)
       END DO
    END DO

    DO i=1,imax
       us=sqrt(cd(i)*(u(i,1)*u(i,1)+v(i,1)*v(i,1)))
       vt1=(b1**0.666667)*us*us
       vt2=(b1**0.6666667)*kappa*kappa*                                 &
            &     m2(i,2)*(zi(i,2)-zi(i,1))*(zi(i,2)-zi(i,1))
       !       q2(i,1)=amax1(q2min,vt1,vt2)
       q2(i,1)=amax1(q2min,vt1)
       long(i,1)=0.0
       long(i,kmax+1)=long(i,kmax)
       sq=0.2
       kq(i,1)=amax1(kqmin,kappa*(z(i,1)-zi(i,1))*us*sq)
    END DO

    DO k=2,kmax
       DO i=1,imax
          longc=0.5*(long(i,k)+long(i,k+1))
          q2c=0.5*(q2(i,k)+q2(i,k+1))
          sq=0.2
          kq(i,k)=amax1(kqmin,longc*sqrt(q2c)*sq)
       END DO
    END DO

    !**********************************************************************


    !************* CALCUL DE Q2 *******************************************

    DO k=2,kmax
       DO i=1,imax
          prod=2.0*(km(i,k)*m2(i,k)+amax1(0.0,-kh(i,k)*n2(i,k)))
          dest=2.0*(amax1(0.0,kh(i,k)*n2(i,k))                            &
               &           +q2(i,k)*sqrt(q2(i,k))/(b1*long(i,k)))
          IF(k.lt.kmax) THEN
             ksdz2sup=unsdzi(i,k)*unsdz(i,k)*kq(i,k)
          ELSE
             ksdz2sup=0.0
          END IF
          ksdz2inf=unsdzi(i,k)*unsdz(i,k-1)*kq(i,k-1)
          b(i,k)=-ksdz2inf*dt
          a(i,k)=1.0+dt*(dest/q2(i,k)+ksdz2inf+ksdz2sup)
          c(i,k)=-ksdz2sup*dt
          f(i,k)=q2(i,k)+dt*prod
       END DO
    END DO
    DO i=1,imax
       f(i,2)=f(i,2)                                                    &
            &       +dt*unsdzi(i,2)*unsdz(i,1)*kq(i,1)*q2(i,1)
    END DO

    DO i=1,imax
       alph(i,2)=a(i,2)
    END DO

    DO k=3,kmax
       DO i=1,imax
          bet=b(i,k)/alph(i,k-1)
          alph(i,k)=a(i,k)-bet*c(i,k-1)
          f(i,k)=f(i,k)-bet*f(i,k-1)
       END DO
    END DO

    DO i=1,imax
       q2(i,kmax)=amax1(q2min,f(i,kmax)/alph(i,kmax))
       q2(i,kmax+1)=q2(i,kmax)
    END DO

    DO k=kmax-1,2,-1
       DO i=1,imax
          q2(i,k)=amax1(q2min,(f(i,k)-c(i,k)*q2(i,k+1))/alph(i,k))
       END DO
    END DO

    DO i=1,imax
       q2(i,2)=amax1(q2(i,2),q2(i,1))
    END DO

    !**********************************************************************


    !************* EVALUATION FINALE DE KH ET KM **************************

    DO k=2,kmax
       DO i=1,imax
          IF(n2(i,k).gt.0.0) THEN
             long(i,k)=min(0.53*sqrt(q2(i,k))/sqrt(n2(i,k)),long(i,k))
          END IF
          gh=amax1(ghmin,                                                 &
               &         min(ghmax,-long(i,k)*long(i,k)*n2(i,k)/q2(i,k)))
          sm=a1*(1.0-3.0*c1-6.0*a1/b1-3.0*a2*gh*                          &
               &           ((b2-3.0*a2)*(1.0-6.0*a1/b1)-3.0*c1*(b2+6.0*a1)))/     &
               &        ((1.0-3.0*a2*gh*(6.0*a1+b2))*(1.0-9.0*a1*a2*gh))
          km(i,k)=amax1(kmmin,long(i,k)*sqrt(q2(i,k))*sm)
          sh=a2*(1.0-6.0*a1/b1)/(1.0-3.0*a2*gh*(6.0*a1+b2))
          kh(i,k)=amax1(khmin,long(i,k)*sqrt(q2(i,k))*sh)
       END DO
    END DO

    DO i=1,imax
       km(i,1)=kmmin
       km(i,kmax+1)=km(i,kmax)
       kh(i,1)=khmin
       kh(i,kmax+1)=kh(i,kmax)
    END DO

    !**********************************************************************

    am=1.0/float(imax)
    DO k=kmax,1,-1
       au=0.0
       ateta=0.0
       az=0.0
       adz=0.0
       akq=0.0
       acd=0.0
       DO i=1,imax
          au=au+am*sqrt(u(i,k)*u(i,k)+v(i,k)*v(i,k))
          ateta=ateta+am*teta(i,k)
          az=az+am*z(i,k)
          adz=adz+am*(1.0/unsdz(i,k))
          akq=akq+am*kq(i,k)
          acd=acd+am*cd(i)
       END DO
    END DO

    DO k=kmax+1,1,-1
       azi=0.0
       adzi=0.0
       aq2=0.0
       al=0.0
       akm=0.0
       akh=0.0
       am2=0.0
       al0=0.0
       DO i=1,imax
          azi=azi+am*zi(i,k)
          adzi=adzi+am*(1.0/unsdzi(i,k))
          aq2=aq2+am*q2(i,k)
          al=al+am*long(i,k)
          akm=akm+am*km(i,k)
          akh=akh+am*kh(i,k)
          am2=am2+am*m2(i,k)
          !         al0=al0+am*long0d(i)
       END DO
    END DO

  END SUBROUTINE my_25

END MODULE mellor_yamada
