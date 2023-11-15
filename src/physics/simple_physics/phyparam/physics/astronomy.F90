MODULE astronomy

#include "use_logging.h"

  IMPLICIT NONE
  SAVE

  REAL :: aphelie, periheli, year_day, peri_day, obliquit, &
       timeperi, e_elips,p_elips
  REAL, PARAMETER :: unitastr=149.597927, & ! millions of km
       pi=2.*ASIN(1.)

CONTAINS

  SUBROUTINE solarlong(pday,psollong)
    REAL, INTENT(IN) :: pday               ! jour de l annee (le jour 0 correspondant a l equinoxe)
    REAL, INTENT(OUT) :: psollong          ! solar longitude
    LOGICAL, PARAMETER ::  lwrite=.TRUE.

    ! Local:
    ! ------
    REAL zanom,xref,zx0,zdx,zteta,zz
    INTEGER iter

    !--------------------------------------------------------
    ! calcul de l angle polaire et de la distance au soleil :
    ! -------------------------------------------------------

    !  calcul de l zanomalie moyenne

    zz=(pday-peri_day)/year_day
    zanom=2.*pi*(zz-nint(zz))
    xref=abs(zanom)

    !  resolution de l equation horaire  zx0 - e * sin (zx0) = xref
    !  methode de Newton

    zx0=xref+e_elips*sin(xref)
    DO iter=1,10
       zdx=-(zx0-e_elips*sin(zx0)-xref)/(1.-e_elips*cos(zx0))
       zx0=zx0+zdx
    END DO
    zx0=zx0+zdx
    if(zanom.lt.0.) zx0=-zx0

    ! zteta est la longitude solaire

    zteta=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))

    psollong=zteta-timeperi

    IF(psollong.LT.0.) psollong=psollong+2.*pi
    IF(psollong.GT.2.*pi) psollong=psollong-2.*pi
    !-----------------------------------------------------------------------
    !   sorties eventuelles:
    !   ---------------------

    IF (lwrite) THEN
       WRITELOG(*,*) 'day of the year  :',pday
       WRITELOG(*,*) 'solar longitude : ',psollong
       LOG_DBG('solarlong')
    ENDIF

  END SUBROUTINE solarlong

  SUBROUTINE iniorbit
    !=======================================================================
    !
    !   Auteur:
    !   -------
    !     Frederic Hourdin      22 Fevrier 1991
    !
    !   Objet:
    !   ------
    !    Initialisation du sous programme orbite qui calcule
    !    a une date donnee de l annee de duree year_day commencant
    !    a l equinoxe de printemps et dont le perihelie se situe
    !    a la date peri_day, la distance au soleil et la declinaison.
    !
    !   Interface:
    !   ----------
    !   - initialise certaines variables de ce module
    !   - Doit etre appele avant d utiliser orbite.
    !
    !   Arguments:
    !   ----------
    !
    !   Input:
    !   ------
    !   aphelie       \   aphelie et perihelie de l orbite
    !   periheli      /   en millions de kilometres.
    !
    !=======================================================================

    !-----------------------------------------------------------------------

    !   Local:
    !   ------
    REAL zxref,zanom,zz,zx0,zdx
    INTEGER iter

    !-----------------------------------------------------------------------

    WRITELOG(*,*) 'Perihelie en Mkm  ',periheli
    WRITELOG(*,*) 'Aphelise  en Mkm  ',aphelie
    WRITELOG(*,*) 'obliquite en degres  :',obliquit

    e_elips=(aphelie-periheli)/(periheli+aphelie)
    p_elips=0.5*(periheli+aphelie)*(1-e_elips*e_elips)/unitastr

    WRITELOG(*,*) 'e_elips',e_elips
    WRITELOG(*,*) 'p_elips',p_elips

    !-----------------------------------------------------------------------
    ! calcul de l angle polaire et de la distance au soleil :
    ! -------------------------------------------------------

    !  calcul de l zanomalie moyenne

    zz=(year_day-peri_day)/year_day
    zanom=2.*pi*(zz-nint(zz))
    zxref=abs(zanom)
    WRITELOG(*,*) 'zanom  ',zanom

    !  resolution de l equation horaire  zx0 - e * sin (zx0) = zxref
    !  methode de Newton

    zx0=zxref+e_elips*sin(zxref)
    DO  iter=1,100
       zdx=-(zx0-e_elips*sin(zx0)-zxref)/(1.-e_elips*cos(zx0))
       zx0=zx0+zdx
    END DO

    zx0=zx0+zdx
    if(zanom.lt.0.) zx0=-zx0
    WRITELOG(*,*) 'zx0   ',zx0

    ! zteta est la longitude solaire

    timeperi=2.*atan(sqrt((1.+e_elips)/(1.-e_elips))*tan(zx0/2.))
    WRITELOG(*,*) 'longitude solaire du perihelie timeperi = ',timeperi

    LOG_INFO('iniorbit')

  END SUBROUTINE iniorbit

  PURE SUBROUTINE orbite(pls,pdist_sol,pdecli)
    !=======================================================================
    !
    !   Objet:
    !   ------
    !
    !   Distance from sun and declimation as a function of the solar
    !   longitude Ls
    !
    !   Arguments:
    !   ----------
    !
    !   Input:
    !   ------
    !   pls          Ls
    !
    !   Output:
    !   -------
    !   pdist_sol     Distance Sun-Planet in UA
    !   pdecli        declinaison ( en radians )
    !
    !=======================================================================
    !-----------------------------------------------------------------------
    !   Declarations:
    !   -------------

    ! arguments:
    ! ----------

    REAL, INTENT(IN) :: pls
    REAL, INTENT(OUT) :: pdist_sol,pdecli

    !-----------------------------------------------------------------------

    ! Distance Sun-Planet

    pdist_sol=p_elips/(1.+e_elips*cos(pls+timeperi))

    ! Solar declination

    pdecli= asin (sin(pls)*sin(obliquit*pi/180.))

  END SUBROUTINE orbite

END MODULE astronomy
