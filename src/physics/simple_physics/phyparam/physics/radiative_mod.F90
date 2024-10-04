module radiative_mod
  use comgeomfi
  use callkeys, only : diurnal
  implicit none
  private
  save
  real, parameter :: pi = 2.0 * asin (1.0), Solarcst = 1370.0, Stephan = 5.67e-08, height_scale = 10000.0, Ps_rad = 1e5
  public          :: radiative_tendencies
contains
  subroutine radiative_tendencies (ngrid, nlayer, gmTime, zdTime, zDay, pPlev, pPlay, pT, FluxRad)
    USE planet
    use phys_const,     only : planet_rad, unjours
    use astronomy,      only : orbite, solarlong
    USE solar,          only : solang, zenang, mucorr
    use soil_mod,       only : albedo, emissiv, tsurf
    USE radiative_sw,   only : sw
    USE radiative_lw,   only : lw
    use writefield_mod, only : writefield

    ! Input variables
    integer,                         intent(in)    :: ngrid   ! number of columns
    integer,                         intent(in)    :: nlayer  ! number of layers
    real,                            intent(in)    :: gmTime  ! fraction of the day
    real,                            intent(in)    :: zdTime  ! time step [s]
    real,                            intent(in)    :: zDay    ! elapsed days (and fraction thereof)
    real, dimension(ngrid,nlayer),   intent(in)    :: pPlay
    real, dimension(ngrid,nlayer+1), intent(in)    :: pPlev

    ! Output
    real, dimension(ngrid,nlayer),   intent(inout) :: pT       ! temperature [K] (advanced from t -> t+dt)
    real, dimension(ngrid),          intent(out)   :: FluxRad  ! net surface flux

    ! Local variables
    integer                       :: ig, l
    real                          :: zls, zInsol, ztim1, ztim2, ztim3, dist_sol, declin
    real, dimension(ngrid)        :: fract     ! day fraction
    real, dimension(ngrid)        :: zFluxSW   ! short-wave flux at surface
    real, dimension(ngrid)        :: zFluxlw   ! short-wave flux at surface
    real, dimension(ngrid)        :: Mu0       ! cosine of zenithal angle
    real, dimension(ngrid)        :: zPlanck
    real, dimension(ngrid,nlayer) :: zdTsw     ! short-wave temperature tendency
    real, dimension(ngrid,nlayer) :: zdTlw     ! long-wave temperature tendency

    !  Insolation
    call solarlong (zday, zls)
    call orbite (zls, dist_sol, declin)

    if (diurnal) then
       ztim1 =   sin (declin)
       ztim2 =   cos (declin) * cos (2.0 * pi * (zday - 0.5))
       ztim3 = - cos (declin) * sin (2.0 * pi * (zday - 0.5))

       call solang (ngrid, sinlon, coslon, sinlat, coslat,  ztim1, ztim2, ztim3, mu0, fract)
    else
       call mucorr (ngrid, declin, lati, mu0, fract, height_scale, planet_rad)
    end if

    zInsol = solarcst / dist_sol**2

    ! Radiative tendencies and fluxes:
    call sw (ngrid, nlayer, diurnal, CoefVis, Albedo, pPlev, Ps_rad, mu0, fract, zInsol, zFluxSW, zdTsw)
    call lw (ngrid, nlayer, coefir, Emissiv, pPlev, Ps_rad, Tsurf, pT, zFluxlw, zdTlw)

    ! 2.4 Surface fluxes
    FluxRad = Emissiv * zFluxlw + zFluxsw * (1.0 - Albedo)
    zPlanck = Emissiv * Stephan * Tsurf**4
    FluxRad = Fluxrad - zPlanck

    ! Temperature at t+dt
    pT = pT + (zdTsw + zdTlw) * zdTime
  end subroutine radiative_tendencies
end MODULE radiative_mod
