MODULE astronomy
  implicit none
  save
  real            :: aphelie, periheli, year_day, peri_day, obliquit, timeperi, e_elips, p_elips
  real, parameter :: unitastr = 149.597927 ! millions of km
  real, parameter :: pi = 2.0 * asin (1.0)
contains
  subroutine SolarLong (pDay, pSolLong)
    ! Called by radiation modules
    real, intent(in)   :: pDay      ! jour de l annee (le jour 0 correspondant a l equinoxe)
    real, intent(out)  :: pSolLong  ! solar longitude

    real    :: zanom, xref, zx0, zdx, zteta, zz
    integer :: iter

    ! Calcul de l angle polaire et de la distance au soleil :

    ! Calcul de l zanomalie moyenne

    zz    = (pDay - peri_day) / year_day
    zanom = 2.0*pi * (zz - nint (zz))
    xref  = abs (zanom)

    !  Resolution de l equation horaire  zx0 - e * sin (zx0) = xref
    !  (methode de Newton)

    zx0 = xref + e_elips * sin (xref)
    do iter = 1, 10
       zdx = - (zx0 - e_elips * sin (zx0) - xref) / (1.0 - e_elips * cos (zx0))
       zx0 = zx0 + zdx
    end do
    
    zx0 = zx0 + zdx
    if (zanom < 0.0) zx0 = -zx0

    ! zteta est la longitude solaire

    zteta = 2.0 * atan (sqrt ((1.0 + e_elips) / (1.0 - e_elips)) * tan (zx0/2.0))

    pSolLong =zteta - timeperi

    if (pSolLong < 0.0)     pSolLong = pSolLong + 2.0*pi
    if (pSolLong >  2.0*pi) pSolLong = pSolLong - 2.0*pi
  end subroutine SolarLong

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
    real    :: zxref, zanom, zz, zx0, zdx
    integer :: iter

    e_elips = (aphelie - periheli) / (periheli + aphelie)
    p_elips = 0.5 * (periheli + aphelie)*(1.0 - e_elips**2) / unitastr

    !-----------------------------------------------------------------------
    ! calcul de l angle polaire et de la distance au soleil :
    ! -------------------------------------------------------

    !  calcul de l zanomalie moyenne

    zz    = (year_day - peri_day) / year_day
    zanom = 2.0*pi * (zz - nint (zz))
    zxref = abs (zanom)

    !  resolution de l equation horaire  zx0 - e * sin (zx0) = zxref
    !  methode de newton

    zx0 = zxref + e_elips * sin (zxref)
    do  iter = 1, 100
       zdx = - (zx0 - e_elips * sin (zx0) - zxref) / (1.0 - e_elips * cos (zx0))
       zx0 = zx0 + zdx
    end do

    zx0 = zx0 + zdx
    if (zanom < 0.0) zx0 = -zx0

    ! zteta est la longitude solaire

    timeperi = 2.0 * atan (sqrt ((1.0 + e_elips) / (1.0 - e_elips)) * tan (zx0/2.0))
  end subroutine iniorbit

  pure subroutine orbite (pls, pdist_sol, pdecli)
    !==============================================================================
    !
    !   Distance from sun and declimation as a function of the solar longitude ls
    !
    !==============================================================================
    real, intent(in)  :: pls
    real, intent(out) :: pdist_sol,pdecli

    ! Distance sun-planet
    pdist_sol = p_elips / (1.0 + e_elips  *cos (pls+timeperi))

    ! Solar declination
    pdecli = asin (sin (pls) * sin (obliquit * pi/180.0))
  end subroutine orbite
end module astronomy
