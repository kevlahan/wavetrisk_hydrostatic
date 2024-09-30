module radiative_mod
#include "use_logging.h"
  use comgeomfi
  use callkeys, only : diurnal, lverbose
  implicit none
  private
  save
  real, parameter :: pi = 2.0 * asin(1.0), solarcst = 1370.0, stephan = 5.67e-08, height_scale = 10000.0, ps_rad = 1.e5
  public          :: radiative_tendencies
contains
  subroutine radiative_tendencies (lwrite, ngrid, igout, nlayer, gmTime, zdTime, zDay, pPlev, pPlay, pT, pdT, FluxRad)
    USE planet
    use phys_const,     only : planet_rad, unjours
    use astronomy,      only : orbite, solarlong
    USE solar,          only : solang, zenang, mucorr
    use soil_mod,       only : albedo, emissiv, tsurf
    USE radiative_sw,   only : sw
    USE radiative_lw,   only : lw
    use writefield_mod, only : writefield

    ! Input variables
    integer,                         intent(in)    :: ngrid, igout, nlayer
    real,                            intent(in)    :: gmTime            ! fraction of the day
    real,                            intent(in)    :: zdTime            ! time step [s]
    real,                            intent(in)    :: zDay              ! elapsed days (and fraction thereof)
    real, dimension(ngrid,nlayer),   intent(in)    :: pPlay
    real, dimension(ngrid,nlayer+1), intent(in)    :: pPlev
    real, dimension(ngrid,nlayer),   intent(in)    :: pT                ! temperature [K]
    logical,                         intent(in)    :: lwrite            ! true for output

    real, dimension(ngrid,nlayer),   intent(inout) :: pdT               ! tendency of temperature [K/s]
    real, dimension(ngrid),          intent(out)   :: fluxrad           ! net surface flux

    ! Local variables
    integer                       :: ig, l
    real                          :: zls, zInsol, ztim1,ztim2,ztim3, dist_sol, declin
    real, dimension(ngrid)        :: fract     ! day fraction
    real, dimension(ngrid)        :: zFluxSW   ! short-wave flux at surface
    real, dimension(ngrid)        :: zFluxlw   ! short-wave flux at surface
    real, dimension(ngrid)        :: Mu0       ! cosine of zenithal angle
    real, dimension(ngrid)        :: zPlanck
    real, dimension(ngrid,nlayer) :: zdTsw     ! short-wave temperature tendency
    real, dimension(ngrid,nlayer) :: zdTlw     ! long-wave temperature tendency

    ! ----------------------------------------------------
    !   2.1 Insolation
    ! ----------------------------------------------------
    call solarlong (zday, zls)
    call orbite (zls, dist_sol, declin)

    if (diurnal) THEN
       ztim1 =   SIN (declin)
       ztim2 =   COS (declin) * COS (2.0 * pi * (zday - 0.5))
       ztim3 = - COS (declin) * SIN (2.0 * pi * (zday - 0.5))

       call solang (ngrid, sinlon, coslon, sinlat, coslat,  ztim1, ztim2, ztim3, mu0, fract)
       if (lverbose) THEN
          WRITELOG (*,*) 'day, declin, sinlon,coslon, sinlat, coslat'
          WRITELOG (*,*) zday, declin, sinlon(igout), coslon(igout), sinlat(igout), coslat(igout)
          LOG_DBG ('radiative_tendencies')
       end if
    else
       WRITELOG (*,*) 'declin,ngrid,planet_rad', declin, ngrid, planet_rad
       LOG_DBG ('radiative_tendencies')
       call mucorr (ngrid, declin, lati, mu0,fract, height_scale, planet_rad)
    end if

    zinsol = solarcst / dist_sol**2

    ! ----------------------------------------------------
    !    2.2 Radiative tendencies and fluxes:
    ! ----------------------------------------------------
    call sw (ngrid, nlayer, diurnal, CoefVis, Albedo, pPlev, pS_Rad, mu0, fract, zInsol, zFluxSW, zdTsw, lverbose, lwrite)
    call lw (ngrid, nlayer, coefir, Emissiv, pPlev, pS_rad, Tsurf, pT, zFluxlw, zdTlw, lverbose, lwrite)


    ! ----------------------------------------------------
    !    2.4 Surface fluxes
    ! ----------------------------------------------------
    Fluxrad = Emissiv * zfluxlw + zFluxsw * (1.0 - Albedo)
    zPlanck = Emissiv * Stephan * Tsurf**4
    FluxRad = Fluxrad - zPlanck


    ! ----------------------------------------------------
    !    2.5 Temperature tendencies
    ! ----------------------------------------------------
    pdT = pdT + zdTsw + zdTlw


    if (lverbose) THEN
       WRITELOG (*,*) 'Diagnostics for radiation'
       WRITELOG (*,*) 'albedo, emissiv, mu0,fract,Frad,Planck'
       WRITELOG (*,*) Albedo(igout), Emissiv(igout), mu0(igout), fract(igout), FluxRad(igout), zPlanck(igout)
       WRITELOG (*,*) 'Tlay Play Plev dT/dt SW dT/dt lw (K/day)'
       WRITELOG (*,*) 'unjours',unjours
       DO l = 1,nlayer
          WRITELOG (*,'(3f15.5,2e15.2)') pT(igout,l), pPlay(igout,l), pPlev(igout,l), zdTsw(igout,l), zdTlw(igout,l)
       endDO
       LOG_DBG ('radiative_tendencies')
    end if

    if (lwrite) THEN
       call writefield ('mu0',   'Cosine zenithal angle', '',     mu0)
       call writefield ('swsurf','SW surf',               'W/m2', zFluxsw)
       call writefield ('lwsurf','LW surf',               'W/m2', zFluxlw)
       call writefield ('dtsw',  'dtsw',                  ' ',    zdTsw)
       call writefield ('dtlw',  'dtlw',                  ' ',    zdTsw)
    end if
  end subroutine radiative_tendencies
end MODULE radiative_mod
