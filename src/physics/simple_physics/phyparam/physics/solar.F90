MODULE solar

#include "use_logging.h"

  IMPLICIT NONE
  PRIVATE

  REAL, PARAMETER :: pi_local = 4.0 * ATAN(1.0), deux_pi_local = 2.0 * pi_local

  PUBLIC :: solang, zenang, mucorr

CONTAINS

  PURE SUBROUTINE solang ( kgrid,psilon,pcolon,psilat,pcolat, &
       ptim1,ptim2,ptim3,pmu0,pfract )

    !-----------------------------------------------------------------------
    !          CALCULATES THE SOLAR ANGLE FOR ALL THE POINTS OF THE GRID
    !
    !     ==== INPUTS  ===
    !
    ! PSILON(KGRID)   : SINUS OF THE LONGITUDE
    ! PCOLON(KGRID)   : COSINUS OF THE LONGITUDE
    ! PSILAT(KGRID)   : SINUS OF THE LATITUDE
    ! PCOLAT(KGRID)   : COSINUS OF THE LATITUDE
    ! PTIM1           : SIN(DECLI)
    ! PTIM2           : COS(DECLI)*COS(TIME)
    ! PTIM3           : SIN(DECLI)*SIN(TIME)
    !
    !     ==== OUTPUTS ===
    !
    ! PMU0 (KGRID)    : SOLAR ANGLE
    ! PFRACT(KGRID)   : DAY FRACTION OF THE TIME INTERVAL
    !
    !     REFERENCE.
    !     ----------
    !
    !         RADIATIVE PROCESSES IN METEOROLOGIE AND CLIMATOLOGIE
    !         PALTRIDGE AND PLATT
    !
    !     AUTHOR.
    !     -------
    !        FREDERIC HOURDIN
    !
    !     MODIFICATIONS.
    !     --------------
    !        ORIGINAL :90-01-14
    !                  92-02-14 CALCULATIONS DONE THE ENTIER GRID (J.Polcher)
    !-----------------------------------------------------------------------

    INTEGER, INTENT(IN)  :: kgrid
    REAL,    INTENT(IN)  :: ptim1,ptim2,ptim3
    REAL,    INTENT(IN)  :: psilon(kgrid), pcolon(kgrid)
    REAL,    INTENT(IN)  :: psilat(kgrid), pcolat(kgrid)
    REAL,    INTENT(OUT) :: pmu0(kgrid), pfract(kgrid)

    INTEGER jl
    REAL ztim1,ztim2,ztim3
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !------------------------------------------------------------------------
    !
    !------------------------------------------------------------------------
    !
    !*     1.     INITIALISATION
    !             --------------
    !
    !
    !$acc kernels present(psilon, pcolon, psilat, pcolat, pmu0, pfract)
    DO jl=1,kgrid
       pmu0(jl)=0.
       pfract(jl)=0.
    ENDDO
    !
    !*     1.1     COMPUTATION OF THE SOLAR ANGLE
    !              ------------------------------
    !
    !$acc loop private(ztim1, ztim2, ztim3)
    DO jl=1,kgrid
       ztim1=psilat(jl)*ptim1
       ztim2=pcolat(jl)*ptim2
       ztim3=pcolat(jl)*ptim3
       pmu0(jl)=ztim1+ztim2*pcolon(jl)+ztim3*psilon(jl)
    ENDDO
    !
    !*     1.2      DISTINCTION BETWEEN DAY AND NIGHT
    !               ---------------------------------
    !
    DO jl=1,kgrid
       IF (pmu0(jl).gt.0.) THEN
          pfract(jl)=1.
       ELSE
          pmu0(jl)=0.
          pfract(jl)=0.
       ENDIF
    ENDDO
    !$acc end kernels

  END SUBROUTINE solang

  SUBROUTINE zenang(klon,longi,gmtime,pdtrad,lat,long,              &
       &                  pmu0,frac)
    USE astronomy

    !=============================================================
    ! Auteur : O. Boucher (LMD/CNRS)
    !          d apres les routines zenith et angle de Z.X. Li
    ! Objet  : calculer les valeurs moyennes du cos de l angle zenithal
    !          et l ensoleillement moyen entre gmtime1 et gmtime2
    !          connaissant la declinaison, la latitude et la longitude.
    ! Rque   : Different de la routine angle en ce sens que zenang
    !          fournit des moyennes de pmu0 et non des valeurs
    !          instantanees, du coup frac prend toutes les valeurs
    !          entre 0 et 1.
    ! Date   : premiere version le 13 decembre 1994
    !          revu pour  GCM  le 30 septembre 1996
    !===============================================================
    ! longi----INPUT : la longitude vraie de la terre dans son plan
    !                  solaire a partir de l equinoxe de printemps (degre)
    ! gmtime---INPUT : temps universel en fraction de jour
    ! pdtrad---INPUT : pas de temps du rayonnement (secondes)
    ! lat------INPUT : latitude en degres
    ! long-----INPUT : longitude en degres
    ! pmu0-----OUTPUT: angle zenithal moyen entre gmtime et gmtime+pdtrad
    ! frac-----OUTPUT: ensoleillement moyen entre gmtime et gmtime+pdtrad
    !================================================================
    integer klon
    !================================================================
    real longi, gmtime, pdtrad
    real lat(klon), long(klon), pmu0(klon), frac(klon)
    !================================================================
    integer i
    real gmtime1, gmtime2
    real incl
    real omega1, omega2, omega
    ! omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
    ! omega : heure en radian du coucher de soleil
    ! -omega est donc l heure en radian de lever du soleil
    real omegadeb, omegafin
    real zfrac1, zfrac2, z1_mu, z2_mu
    ! declinaison en radian
    real lat_sun
    ! longitude solaire en radian
    real lon_sun
    ! latitude du pt de grille en radian
    real latr
    !================================================================
    !
    !     incl=R_incl * pi_local / 180.
    WRITELOG(*,*) 'Obliquite =' ,obliquit
    LOG_INFO('solar')

    incl=obliquit * pi_local / 180.
    !
    !     lon_sun = longi * pi_local / 180.0
    lon_sun = longi
    lat_sun = ASIN (SIN(lon_sun)*SIN(incl) )
    !
    gmtime1=gmtime*86400.
    gmtime2=gmtime*86400.+pdtrad
    !
    DO i = 1, klon
       !
       !     latr = lat(i) * pi_local / 180.
       latr = lat(i)
       !
       !--pose probleme quand lat=+/-90 degres
       !
       !      omega = -TAN(latr)*TAN(lat_sun)
       !      omega = ACOS(omega)
       !      IF (latr.GE.(pi_local/2.+lat_sun)
       !     .    .OR. latr.LE.(-pi_local/2.+lat_sun)) THEN
       !         omega = 0.0       ! nuit polaire
       !      ENDIF
       !      IF (latr.GE.(pi_local/2.-lat_sun)
       !     .          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
       !         omega = pi_local  ! journee polaire
       !      ENDIF
       !
       !--remplace par cela (le cas par defaut est different)
       !
       !--nuit polaire
       omega=0.0
       IF (latr.GE.(pi_local/2.-lat_sun)                                 &
            &          .OR. latr.LE.(-pi_local/2.-lat_sun)) THEN
          ! journee polaire
          omega = pi_local
       ENDIF
       IF (latr.LT.(pi_local/2.+lat_sun).AND.                            &
            &    latr.GT.(-pi_local/2.+lat_sun).AND.                           &
            &    latr.LT.(pi_local/2.-lat_sun).AND.                            &
            &    latr.GT.(-pi_local/2.-lat_sun)) THEN
          omega = -TAN(latr)*TAN(lat_sun)
          omega = ACOS(omega)
       ENDIF
       !
       omega1 = gmtime1 + long(i)*86400.0/360.0
       omega1 = omega1 / 86400.0*deux_pi_local
       omega1 = MOD (omega1+deux_pi_local, deux_pi_local)
       omega1 = omega1 - pi_local
       !
       omega2 = gmtime2 + long(i)*86400.0/360.0
       omega2 = omega2 / 86400.0*deux_pi_local
       omega2 = MOD (omega2+deux_pi_local, deux_pi_local)
       omega2 = omega2 - pi_local
       !
       !--on est dans la meme journee locale
       IF (omega1.LE.omega2) THEN
          !
          IF (omega2.LE.-omega .OR. omega1.GE.omega .OR. omega.LT.1e-5)  &
               THEN
             !--nuit
             frac(i)=0.0
             pmu0(i)=0.0
             !--jour+nuit/jou
          ELSE
             omegadeb=MAX(-omega,omega1)
             omegafin=MIN(omega,omega2)
             frac(i)=(omegafin-omegadeb)/(omega2-omega1)
             pmu0(i)=SIN(latr)*SIN(lat_sun) + COS(latr)*COS(lat_sun)*    &
                  (SIN(omegafin)-SIN(omegadeb))/ (omegafin-omegadeb)
          ENDIF
          !
          !---omega1 GT omega2 -- a cheval sur deux journees
       ELSE
          !
          !-------------------entre omega1 et pi
          !--nuit
          IF (omega1.GE.omega) THEN
             zfrac1=0.0
             z1_mu =0.0
             !--jour+nuit
          ELSE
             omegadeb=MAX(-omega,omega1)
             omegafin=omega
             zfrac1=omegafin-omegadeb
             z1_mu =SIN(latr)*SIN(lat_sun) + COS(latr)*COS(lat_sun)*     &
                  (SIN(omegafin)-SIN(omegadeb))/ (omegafin-omegadeb)
          ENDIF
          !---------------------entre -pi et omega2
          !--nuit
          IF (omega2.LE.-omega) THEN
             zfrac2=0.0
             z2_mu =0.0
             !--jour+nuit
          ELSE
             omegadeb=-omega
             omegafin=MIN(omega,omega2)
             zfrac2=omegafin-omegadeb
             z2_mu =SIN(latr)*SIN(lat_sun) + COS(latr)*COS(lat_sun)*     &
                  (SIN(omegafin)-SIN(omegadeb))/ (omegafin-omegadeb)
             !
          ENDIF
          !-----------------------moyenne
          frac(i)=(zfrac1+zfrac2)/(omega2+deux_pi_local-omega1)
          pmu0(i)=(zfrac1*z1_mu+zfrac2*z2_mu)/MAX(zfrac1+zfrac2,1.E-10)
          !
          !---comparaison omega1 et omega2
       ENDIF
       !
    ENDDO
    !
  END SUBROUTINE zenang

  PURE SUBROUTINE mucorr(npts,pdeclin, plat, pmu, pfract,phaut,prad)

    !=======================================================================
    !
    !   Calcul of equivalent solar angle and and fraction of day whithout
    !   diurnal cycle.
    !
    !   parmeters :
    !   -----------
    !
    !      Input :
    !      -------
    !         npts             number of points
    !         pdeclin          solar declinaison
    !         plat(npts)        latitude
    !         phaut            hauteur typique de l atmosphere
    !         prad             rayon planetaire
    !
    !      Output :
    !      --------
    !         pmu(npts)          equivalent cosinus of the solar angle
    !         pfract(npts)       fractionnal day
    !
    !=======================================================================

    !-----------------------------------------------------------------------
    !
    !    0. Declarations :
    !    -----------------

    !     Arguments :
    !     -----------
    INTEGER, INTENT(IN) :: npts
    REAL, INTENT(IN)    :: phaut, prad, pdeclin, plat(npts)
    REAL, INTENT(OUT)   :: pmu(npts), pfract(npts)
    !
    !     Local variables :
    !     -----------------
    INTEGER j
    REAL z,cz,sz,tz,phi,cphi,sphi,tphi
    REAL ap,a,t,b
    REAL alph

    REAL, PARAMETER :: pi=2.*asin(1.)

    !-----------------------------------------------------------------------

    z = pdeclin
    cz = cos (z)
    sz = sin (z)

    DO j = 1, npts

       phi = plat(j)
       cphi = cos(phi)
       if (cphi.le.1.e-9) cphi=1.e-9
       sphi = sin(phi)
       tphi = sphi / cphi
       b = cphi * cz
       t = -tphi * sz / cz
       a = 1.0 - t*t
       ap = a

       IF(t.eq.0.) then
          t=0.5*pi
       ELSE
          IF (a.lt.0.) a = 0.
          t = sqrt(a) / t
          IF (t.lt.0.) then
             t = -atan (-t) + pi
          ELSE
             t = atan(t)
          ENDIF
       ENDIF

       pmu(j) = (sphi*sz*t) / pi + b*sin(t)/pi
       pfract(j) = t / pi
       IF (ap .lt.0.) then
          pmu(j) = sphi * sz
          pfract(j) = 1.0
       ENDIF
       IF (pmu(j).le.0.0) pmu(j) = 0.0
       pmu(j) = pmu(j) / pfract(j)
       IF (pmu(j).eq.0.) pfract(j) = 0.

    END DO

    !-----------------------------------------------------------------------
    !   correction de rotondite:
    !   ------------------------

    alph=phaut/prad
    DO j=1,npts
       ! !!!!!!
       pmu(j)=sqrt(1224.*pmu(j)*pmu(j)+1.)/35.
       !    $          (sqrt(alph*alph*pmu(j)*pmu(j)+2.*alph+1.)-alph*pmu(j))
    END DO

  END SUBROUTINE mucorr

END MODULE solar
