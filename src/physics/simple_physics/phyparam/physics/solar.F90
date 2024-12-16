module solar
  use phys_const, only : pi
  implicit none
  private
  public :: SolAng, ZenAng, MuCorr
contains
  pure subroutine SolAng (ngrid, pSiLon, pCoLon, pSiLat, pCoLat, pTim1, pTim2, pTim3, pMu0, pFract)
    ! -----------------------------------------------------------------------------------------------
    !   Computes solar angle pMu0 and day fraction of time interval pFract
    !   day fraction is set to 1 if solar angle > 0 (day), otherwise
    !   day fraction and solar angle are set to 0 (night)
    !
    !     ==== Inputs  ===
    !
    ! pSiLon(ngrid)   : sine   of longitude
    ! pCoLon(ngrid)   : cosine of longitude
    ! pSiLat(ngrid)   : sine   of latitude
    ! pCoLat(ngrid)   : cosine of latitude
    !
    ! pTim1           : sin (decli)
    ! pTim2           : cos (decli) * cos(Time)
    ! pTim3           : sin (decli) * sin(Time)
    !
    !     ==== Outputs ===
    !
    ! pMu0  (ngrid)   : cosine of solar zenith angle
    ! pFract(ngrid)   : day fraction of time interval
    !        
    !
    !      Modifications
    !     --------------
    !        Original : 1990-01-14 (F Hourdin)
    !                   1992-02-14 calculations done over entire grid       (J Polcher)
    !                   2024-12-09 simplification and modernization of code (N Kevlahan)
    ! -----------------------------------------------------------------------------------------------

    integer,                intent(in)  :: ngrid
    real,                   intent(in)  :: pTim1, pTim2, pTim3
    real, dimension(ngrid), intent(in)  :: pSiLon, pCoLon, pSiLat, pCoLat
    real, dimension(ngrid), intent(out) :: pMu0, pFract

    integer :: ig
    real    :: zTim1, zTim2, zTim3
  
    ! Computation of solar angle and day fraction
    do ig = 1, ngrid
       zTim1 = pSiLat(ig) * pTim1
       zTim2 = pCoLat(ig) * pTim2
       zTim3 = pCoLat(ig) * pTim3
       
       pMu0(ig) = zTim1 + zTim2 * pCoLon(ig) + zTim3 * pSiLon(ig)  ! cosine of solar zenith angle

       if (pMu0(ig) >  0.0) then ! day
          pFract(ig) = 1.0
       else                      ! night
          pMu0(ig)   = 0.0 
          pFract(ig) = 0.0
       end if
    enddo
  end subroutine SolAng

  subroutine ZenAng (klon, longi, gmTime, pdtrad, lat, long, pMu0, frac)
    use astronomy
    !=============================================================
    ! auteur : O. Boucher (LMD/CNRS)
    !          d apres les routines zenith et angle de z.x. li
    ! objet  : calculer les valeurs moyennes du cos de l angle zenithal
    !          et l ensoleillement moyen entre gmTime1 et gmTime2
    !          connaissant la declinaison, la latitude et la longitude.
    ! rque   : different de la routine angle en ce sens que zenang
    !          fournit des moyennes de pMu0 et non des valeurs
    !          instantanees, du coup frac prend toutes les valeurs
    !          entre 0 et 1.
    ! date   : premiere version le 13 decembre  1994
    !          revu pour gcm    le 30 septembre 1996
    !===============================================================
    ! longi----input : la longitude vraie de la terre dans son plan
    !                  solaire a partir de l equinoxe de printemps (degre)
    ! gmTime---input : temps universel en fraction de jour
    ! pdtrad---input : pas de temps du rayonnement (secondes)
    ! lat------input : latitude en degres
    ! long-----input : longitude en degres
    ! pMu0-----output: angle zenithal moyen entre gmTime et gmTime+pdtrad
    ! frac-----output: ensoleillement moyen entre gmTime et gmTime+pdtrad
    !================================================================
    ! omega1, omega2 : temps 1 et 2 exprime en radian avec 0 a midi.
    ! omega : heure en radian du coucher de soleil
    ! -omega est donc l heure en radian de lever du soleil
    
    integer               :: i, klon
    real                  :: gmTime, gmTime1, gmTime2, incl, longi, omega, omega1, omega2, omegadeb, omegafin, pdTrad
    real                  :: zfrac1, zfrac2, z1_mu, z2_mu
    
    real                  :: lat_sun     ! declinaison en radian
    real                  :: latr        ! latitude du pt de grille en radian
    real                  :: lon_sun     ! longitude solaire en radian
    
    real, dimension(klon) :: lat, long, pMu0, frac

    incl = obliquit * pi / 180.0

    lon_sun = longi
    lat_sun = asin (sin (lon_sun) * sin(incl))
    !
    gmTime1 = gmTime * 86400.0
    gmTime2 = gmTime * 86400.0 + pdTrad
    !
    do i = 1, klon
       latr = lat(i)

       !--pose probleme quand lat=+/-90 degres
       !
       !      omega = -tan(latr)*tan(lat_sun)
       !      omega = acos(omega)
       !      if (latr.ge.(pi/2.+lat_sun)
       !     .    .or. latr.le.(-pi/2.+lat_sun)) then
       !         omega = 0.0       ! nuit polaire
       !      endif
       !      if (latr.ge.(pi/2.-lat_sun)
       !     .          .or. latr.le.(-pi/2.-lat_sun)) then
       !         omega = pi  ! journee polaire
       !      endif
       !
       !--remplace par cela (le cas par defaut est different)
       !
       !--nuit polaire
       omega = 0.0
       if (latr >= (pi/2.0 - lat_sun) .or. latr <= (-pi/2.0 - lat_sun)) then ! journee polaire
          omega = pi
       endif
       if (latr < (pi/2.0 + lat_sun) .and. latr > (-pi/2.0 + lat_sun) .and. &
           latr < (pi/2.0 - lat_sun) .and. latr > (-pi/2.0 - lat_sun)) then
          omega = acos (- tan (latr) * tan (lat_sun))
       endif

       omega1 = gmTime1 + long(i) * 86400.0 / 360.0
       omega1 = omega1 / 86400.0 * 2.0*pi
       omega1 = mod (omega1 + 2.0*pi, 2.0*pi)
       
       omega1 = omega1 - pi

       omega2 = gmTime2 + long(i) * 86400.0 / 360.0
       omega2 = omega2 / 86400.0 * 2.0*pi
       omega2 = mod (omega2 + 2.0*pi, 2.0*pi)
       
       omega2 = omega2 - pi

       ! On est dans la meme journee locale
       if (omega1 <= omega2) then
          if (omega2 <= -omega .or. omega1 >= omega .or. omega < 1e-5) then !--nuit
             frac(i) = 0.0
             pMu0(i) = 0.0
          else !--jour + nuit/jour
             omegadeb = max (-omega, omega1)
             omegafin = min ( omega, omega2)
             
             frac(i) = (omegafin - omegadeb) / (omega2 - omega1)
             pMu0(i) = sin (latr) * sin (lat_sun) + cos (latr) * cos (lat_sun) * &
                  (sin (omegafin) - sin (omegadeb)) / (omegafin - omegadeb)
          endif
       else  !---omega1 gt omega2 -- a cheval sur deux journees
          if (omega1 >= omega) then !--nuit
             zfrac1 = 0.0
             z1_mu  = 0.0
          else !--jour+nuit
             omegadeb = max (-omega, omega1)
             omegafin = omega
             zfrac1 = omegafin - omegadeb
             z1_mu  = sin (latr) * sin(lat_sun) + cos (latr) * cos (lat_sun) *     &
                  (sin(omegafin) - sin (omegadeb)) / (omegafin - omegadeb)
          endif

          !---------------------entre -pi et omega2
          if (omega2 <= -omega) then  !--nuit
             zfrac2 = 0.0
             z2_mu  = 0.0
          else  !--jour+nuit
             omegadeb = -omega
             omegafin = min (omega, omega2)
             zfrac2 = omegafin - omegadeb
             z2_mu = sin (latr) * sin (lat_sun) + cos (latr) * cos (lat_sun) *     &
                  (sin (omegafin) - sin (omegadeb)) / (omegafin - omegadeb)
          endif

          !-----------------------moyenne
          frac(i) = (zfrac1 + zfrac2) / (omega2 + 2.0*pi - omega1)
          pMu0(i) = (zfrac1 * z1_mu + zfrac2 * z2_mu)   /max (zfrac1 + zfrac2, 1e-10)
       endif
    enddo
  end subroutine ZenAng

  pure subroutine MuCorr (ngrid, pDeclin, pLat, pMu, pFract, pHaut, pRad)
    !   Equivalent solar angle and fraction without diurnal cycle
    !
    !      Input :
    !      -------
    !         ngrid             number of points
    !         pDeclin          solar declination
    !         pLat(ngrid)       latitude
    !         pHaut            hauteur typique de l atmosphere
    !         pRad             rayon planetaire
    !
    !      Output :
    !      --------
    !         pMu   (ngrid)     equivalent cosine of solar angle
    !         pFract(ngrid)     fractional day
    
    integer,               intent(in) :: ngrid
    real,                  intent(in) :: pHaut, pRad, pDeclin
    real, dimension(ngrid), intent(in) :: pLat
    
    real, dimension(ngrid), intent(out):: pMu, pFract

    integer :: ig
    real    :: alph, ap, a, t, b, z, cz, sZ, tZ, Phi, cPhi, sPhi, tPhi

    z  = pdeclin
    cz = cos (z)
    sZ = sin (z)

    do ig = 1, ngrid
       Phi = plat(ig)
       cPhi = cos (Phi)
       if (cPhi <= 1e-9) cPhi = 1e-9
       
       sPhi = sin (Phi)
       tPhi = sPhi / cPhi
       
       b  =   cPhi * cz
       t  = - tPhi * sZ / cz
       a  = 1.0 - t**2
       ap = a

       if (t == 0.0) then
          t = 0.5 * pi
       else
          if (a < 0.0) a = 0.0
          t = sqrt (a) / t
          if (t < 0.0) then
             t = - atan (-t) + pi
          else
             t = atan(t)
          end if
       end if

       pMu(ig)    = (sPhi * sZ * t) / pi + b * sin (t) / pi
       pFract(ig) = t / pi
       
       if (ap < 0.0) then
          pMu(ig)    = sPhi * sZ
          pfract(ig) = 1.0
       end if
       
       if (pMu(ig) <= 0.0) pMu(ig) = 0.0

       pMu(ig) = pMu(ig) / pFract(ig)
       if (pMu(ig) == 0.0) pFract(ig) = 0.0
       
       ! Correction de rotondite
       pMu(ig) = sqrt (1224.0 * pMu(ig)**2 + 1.0) / 35.0  
    end do
  end subroutine MuCorr
end module solar
