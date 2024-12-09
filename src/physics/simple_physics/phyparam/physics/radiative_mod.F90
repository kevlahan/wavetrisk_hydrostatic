module radiative_mod
  use comgeomfi
  use callkeys,   only : diurnal
  use phys_const, only : pi, Stefan
  implicit none
  private
  save
  real, parameter :: SolarCst     = 1.37e3
  real, parameter :: height_scale = 1.00e4
  real, parameter :: Ps_rad       = 1.00e5
  public :: radiative_tendencies
contains
  subroutine radiative_tendencies (ngrid, nlayer, gmTime, pTimestep, zDay, pPint, pT, FluxRad)
    USE planet
    use phys_const,     only : planet_rad
    use astronomy,      only : Orbite, SolarLong
    USE solar,          only : SolAng, ZenAng, MuCorr
    use soil_mod,       only : Albedo, Emissiv, Tsurf
    USE radiative_sw,   only : sw
    USE radiative_lw,   only : lw

    ! Input variables
    integer,                         intent(in)    :: ngrid     ! number of columns
    integer,                         intent(in)    :: nlayer    ! number of layers
    real,                            intent(in)    :: gmTime    ! fraction of the day
    real,                            intent(in)    :: pTimestep ! time step [s]
    real,                            intent(in)    :: zDay      ! elapsed days (and fraction thereof)
    real, dimension(ngrid,nlayer+1), intent(in)    :: pPint     ! pressure at interfaces

    ! Output
    real, dimension(ngrid,nlayer),   intent(inout) :: pT        ! temperature [K] (advanced from t -> t+dt)
    real, dimension(ngrid),          intent(out)   :: FluxRad   ! net surface flux

    ! Local variables
    integer                       :: ig, l
    real                          :: zls, zInsol, zTim1, zTim2, zTim3, dist_sol, declin
    real, dimension(ngrid)        :: Fract                      ! day fraction
    real, dimension(ngrid)        :: zFluxSW                    ! short-wave flux at surface
    real, dimension(ngrid)        :: zFluxlw                    ! short-wave flux at surface
    real, dimension(ngrid)        :: Mu0                        ! cosine of zenithal angle
    real, dimension(ngrid,nlayer) :: zdTsw                      ! short-wave temperature tendency
    real, dimension(ngrid,nlayer) :: zdTlw                      ! long-wave temperature tendency

    !  Insolation
    call SolarLong (zDay, zls)
    call Orbite (zls, dist_sol, declin)

    if (diurnal) then
       zTim1 =   sin (declin)
       zTim2 =   cos (declin) * cos (2.0 * pi * (zDay - 0.5))
       zTim3 = - cos (declin) * sin (2.0 * pi * (zDay - 0.5))

       call SolAng (ngrid, SinLon, CosLon, SinLat, CosLat,  zTim1, zTim2, zTim3, Mu0, Fract)
    else
       call MuCorr (ngrid, declin, lati, Mu0, Fract, height_scale, planet_rad)
    end if

    zInsol = SolarCst / dist_sol**2

    ! Radiative tendencies and fluxes:
    call sw (ngrid, nlayer, diurnal, CoefVis, Albedo, pPint, Ps_rad, Mu0, Fract, zInsol, zFluxSW, zdTsw)
    call lw (ngrid, nlayer, coefir, Emissiv, pPint, Ps_rad, Tsurf, pT, zFluxlw, zdTlw)

    ! Surface fluxes
    FluxRad = Emissiv * (zFluxlw + (1.0 - Albedo) * zFluxsw - Stefan * Tsurf**4)

    ! Temperature at t+dt
    pT = pT + pTimestep * (zdTsw + zdTlw)
  end subroutine radiative_tendencies
end MODULE radiative_mod
