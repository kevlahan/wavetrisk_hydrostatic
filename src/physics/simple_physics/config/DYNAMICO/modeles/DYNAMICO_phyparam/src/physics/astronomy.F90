module astronomy
  implicit none
  save
  real            :: aphelie, periheli, year_day, peri_day, obliquit, timeperi, e_elips, p_elips
  real, parameter :: UnitAstr = 149.597927 ! millions of km
  real, parameter :: pi = 2.0 * asin (1.0)
contains
  subroutine SolarLon (pDay, pSolLon)
    ! Returns solar longitude given day (from equinox)
    real, intent(in)   :: pDay     ! jour de l annee (le jour 0 correspondant a l equinoxe)
    real, intent(out)  :: pSolLon  ! solar longitude

    real    :: zAnom, xref, zx0, zdx, zTheta, zz
    integer :: iter

    ! Calcul de l angle polaire et de la distance au soleil (zAnomalie moyenne)
    zz    = (pDay - peri_day) / year_day
    zAnom = 2.0*pi * (zz - nint (zz))
    xref  = abs (zAnom)

    ! Resolution de l equation horaire  zx0 - e * sin (zx0) = xref
    !  (methode de Newton)
    zx0 = xref + e_elips * sin (xref)
    do iter = 1, 10
       zdx = - (zx0 - e_elips * sin (zx0) - xref) / (1.0 - e_elips * cos (zx0))
       zx0 = zx0 + zdx
    end do

    zx0 = zx0 + zdx
    if (zAnom < 0.0) zx0 = -zx0

    ! zTheta est la longitude solaire
    zTheta = 2.0 * atan (sqrt ((1.0 + e_elips) / (1.0 - e_elips)) * tan (zx0/2.0))

    pSolLon = zTheta - timeperi

    if (pSolLon < 0.0)     pSolLon = pSolLon + 2.0*pi
    if (pSolLon >  2.0*pi) pSolLon = pSolLon - 2.0*pi
  end subroutine SolarLon

  subroutine iniorbit
    !=======================================================================
    !
    !     Frederic Hourdin      22 fevrier 1991
    !
    !    Initialisation du sous programme orbite qui calcule
    !    a une date donnee de l annee de duree year_day commencant
    !    a l equinoxe de printemps et dont le perihelie se situe
    !    a la date peri_day, la distance au soleil et la declinaison.
    !
    !   - initialise certaines variables de ce module
    !   - doit etre appele avant d utiliser orbite.
    !
    !
    !=======================================================================
    real    :: zxref, zAnom, zz, zx0, zdx
    integer :: iter

    e_elips = (aphelie - periheli) / (periheli + aphelie)
    p_elips = 0.5 * (periheli + aphelie)*(1.0 - e_elips**2) / UnitAstr

    !-----------------------------------------------------------------------
    ! calcul de l angle polaire et de la distance au soleil :
    ! -------------------------------------------------------

    !  calcul de l zAnomalie moyenne

    zz    = (year_day - peri_day) / year_day
    zAnom = 2.0*pi * (zz - nint (zz))
    zxref = abs (zAnom)

    ! Resolution de l equation horaire  zx0 - e * sin (zx0) = zxref
    ! (methode de newton)
    zx0 = zxref + e_elips * sin (zxref)
    do  iter = 1, 100
       zdx = - (zx0 - e_elips * sin (zx0) - zxref) / (1.0 - e_elips * cos (zx0))
       zx0 = zx0 + zdx
    end do

    zx0 = zx0 + zdx
    if (zAnom < 0.0) zx0 = -zx0

    timeperi = 2.0 * atan (sqrt ((1.0 + e_elips) / (1.0 - e_elips)) * tan (zx0/2.0))
  end subroutine iniorbit

  pure subroutine Orbite (pSolLon, pDist_sol, pDecli)
    !==============================================================================
    !
    !   Distance from sun and declimation as a function of the solar longitude ls
    !
    !==============================================================================
    real, intent(in)  :: pSolLon    ! solar longitude

    real, intent(out) :: pDist_sol  ! distance between sun and planet
    real, intent(out) :: pDecli     ! solar declination angle

    pDist_sol = p_Elips / (1.0 + e_Elips * cos (pSolLon + TimePeri))

    pDecli    = asin (sin (pSolLon) * sin (obliquit * pi/180.0))
  end subroutine Orbite
end module astronomy
