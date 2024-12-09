module solar
  use phys_const, only : pi
  implicit none
  private
  public :: SolAng, ZenAng, MuCorr
contains
  pure subroutine SolAng (kgrid, pSiLon, pCoLon, pSiLat, pCoLat, pTim1, pTim2, pTim3, pMu0, pFract)
    ! -----------------------------------------------------------------------
    !          Calculates the solar angle for all the points of the grid (Frederic Hourdin)
    !
    !     ==== inputs  ===
    !
    ! psilon(kgrid)   : sine   of the longitude
    ! pcolon(kgrid)   : cosine of the longitude
    ! psilat(kgrid)   : sine   of the latitude
    ! pcolat(kgrid)   : cosine of the latitude
    !
    ! pTim1           : sin (decli)
    ! pTim2           : cos (decli) * cos(Time)
    ! pTim3           : sin (decli) * sin(Time)
    !
    !     ==== outputs ===
    !
    ! pMu0  (kgrid)   : solar angle
    ! pFract(kgrid)   : day fraction of the time interval
    !        
    !
    !     modifications.
    !     --------------
    !        original :90-01-14
    !                  92-02-14 calculations done the entier grid (j.polcher)
    ! -----------------------------------------------------------------------

    integer,                intent(in)  :: kgrid
    real,                   intent(in)  :: pTim1, pTim2, pTim3
    real, dimension(kgrid), intent(in)  :: pSiLon, pCoLon, pSiLat, pCoLat
    real, dimension(kgrid), intent(out) :: pMu0, pFract

    integer :: jl
    real    :: Mu0, zTim1, zTim2, zTim3
  
    do jl = 1, kgrid
       pmu0(jl)   = 0.0
       pfract(jl) = 0.0
    enddo

    ! Computation of solar angle Mu0
    do jl = 1, kgrid
       zTim1 = pSiLat(jl) * pTim1
       zTim2 = pCoLat(jl) * pTim2
       zTim3 = pCoLat(jl) * pTim3
       
       Mu0 = zTim1 + zTim2 * pCoLon(jl) + zTim3 * pSiLon(jl)

       if (Mu0 > 0.0) then ! day
          pFract(jl) = 1.0
       else                ! night
          pMu0(jl)   = 0.0
          pFract(jl) = 0.0
       endif
    enddo
  end subroutine SolAng

  subroutine ZenAng (klon, longi, gmTime, pdtrad, lat, long, pMu0, frac)
    use astronomy
    !=============================================================
    ! auteur : O. Boucher (lmd/cnrs)
    !          d apres les routines zenith et angle de z.x. li
    ! objet  : calculer les valeurs moyennes du cos de l angle zenithal
    !          et l ensoleillement moyen entre gmTime1 et gmTime2
    !          connaissant la declinaison, la latitude et la longitude.
    ! rque   : different de la routine angle en ce sens que zenang
    !          fournit des moyennes de pMu0 et non des valeurs
    !          instantanees, du coup frac prend toutes les valeurs
    !          entre 0 et 1.
    ! date   : premiere version le 13 decembre 1994
    !          revu pour  gcm  le 30 septembre 1996
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
    gmTime1 =gmTime * 86400.0
    gmTime2 =gmTime * 86400.0 + pdTrad
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
  end subroutine zenang

  pure subroutine MuCorr (npts, pDeclin, pLat, pMu, pFract, pHaut, pRad)
    !   Equivalent solar angle and fraction without diurnal cycle
    !
    !      Input :
    !      -------
    !         npts             number of points
    !         pDeclin          solar declination
    !         pLat(npts)       latitude
    !         pHaut            hauteur typique de l atmosphere
    !         pRad             rayon planetaire
    !
    !      Output :
    !      --------
    !         pMu   (npts)     equivalent cosine of solar angle
    !         pFract(npts)     fractional day
    
    integer,               intent(in) :: npts
    real,                  intent(in) :: pHaut, pRad, pDeclin
    real, dimension(npts), intent(in) :: pLat
    
    real, dimension(npts), intent(out):: pMu, pFract

    integer :: j
    real    :: alph, ap, a, t, b, z, cz, sZ, tZ, Phi, cPhi, sPhi, tPhi

    z  = pdeclin
    cz = cos (z)
    sZ = sin (z)

    do j = 1, npts
       Phi = plat(j)
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

       pMu(j)    = (sPhi * sZ * t) / pi + b * sin (t) / pi
       pFract(j) = t / pi
       
       if (ap < 0.0) then
          pMu(j)    = sPhi * sZ
          pfract(j) = 1.0
       end if
       
       if (pMu(j) <= 0.0) pMu(j) = 0.0

       pMu(j) = pMu(j) / pFract(j)
       if (pMu(j) == 0.0) pFract(j) = 0.0
       
       ! Correction de rotondite
       pMu(j) = sqrt (1224.0 * pMu(j)**2 + 1.0) / 35.0  
       !alph = pHaut / pRad
       !pMu(j) = sqrt ((alph*pMu(j))**2 + 2.0 * alph + 1.0) - alph * pMu(j)
    end do
  end subroutine MuCorr
end module solar
