module phyparam_mod
#include "use_logging.h"
  use callkeys
  use comgeomfi
  implicit none
  private
  save
  integer         :: icount
  real            :: zday_last
  logical         :: firstcall_alloc=.true.

  real, parameter :: capcal_nosoil = 1e5
  real, parameter :: ref_temp      = 285.0
  real, parameter :: Tsoil_init    = 300.0
  public :: phyparam, alloc, precompute, zday_last, icount
contains
  subroutine phyparam (              &
       ngrid, nlayer,                &
       firstcall, lastcall,          &
       rjourvrai, gmtime, ptimestep, &
       pplev, pplay, pPhi,           &
       pu, pv, pt,                   &
       pdU, pdV, pdt, pdPsrf )       &
       bind(c, name='phyparam_phyparam')
    use,          intrinsic :: iso_c_binding
    use phys_const,     only : g, Rcp, r, unjours
    use soil_mod,       only : soil_forward, soil_backward
    use soil_mod,       only : z0, inertie, Emissiv, Albedo  ! precomputed
    use soil_mod,       only : Tsurf, Tsoil                  ! state variables
    use turbulence,     only : Vdif
    use convection,     only : convadj
    USE radiative_mod,  only : radiative_tendencies
    use writefield_mod, only : writefield
    use accelerator,    only : nvtxstartrange, nvtxendrange
    use profiling,      only : profile_enter, profile_exit, id_phyparam
    !=======================================================================
    !
    !   Physical parametrisations of LMDZ  20 parameter climate model
    !   for planetary atmospheres.
    !
    !   Includes following sub-models:
    !
    !   1. Radiative transfer (long and shortwave) for CO2 and dust.
    !   3. Convective adjustment.
    !   4. Vertical heat diffusion in soil column.
    !   2. Vertical turbulent mixing
    !
    !   Author:  Frederic Hourdin  1993-10-15
    !   Revised: Nicholas Kevlahan 2024-09-27
    !
    !=======================================================================

    ! Input variables
    integer, intent(in), value ::  & !
         ngrid,                    & ! size of the horizontal grid.
         nlayer                      ! number of vertical layers.
    logical(kind=c_bool), intent(in), value  :: & !
         firstcall,                & ! true at the first call
         lastcall                    ! true at the last call
    real, intent(in), value     :: & !
         rjourvrai,                & ! number of days counted from the north. spring equinox
         gmtime,                   & ! fraction of the day (ranges from 0 to 1)
         ptimestep                   ! timestep [s]
    real, intent(in) ::            & !
         pPlev(ngrid,nlayer+1),    & ! pressure at interfaces between layers [Pa]
         pPlay(ngrid,nlayer),      & ! pressure at the middle of the layers [Pa]
         pPhi(ngrid,nlayer),       & ! geopotential at the middle of the layers [m^2/s^2]
         pU(ngrid,nlayer),         & ! zonal velocity [m/s]
         pV(ngrid,nlayer),         & ! meridional velocity [m/s]
         pT(ngrid,nlayer)            ! temperature [K]

    ! Output variables
    real, intent(out)   ::         & !
         pdU(ngrid,nlayer),        & ! tendency of zonal velocity[m/s^2]
         pdV(ngrid,nlayer),        & ! tendency of meridional velocity [m/s^2]
         pdT(ngrid,nlayer),        & ! tendency of temperature [K/s]
         pdPsrf(ngrid)               ! tendency of surface pressure [Pa/s]

    ! Local variables
    real ::                        & !
         zH(ngrid,nlayer),         & ! potential temperature
         zpopsk(ngrid,nlayer),     & ! Exner function
         zZlev(ngrid,nlayer+1),    & ! height of layer interfaces
         zZlay(ngrid,nlayer),      & ! height of layer centres

         zc(ngrid,nsoilmx-1),      & ! Lu coefficients for soil implicit solve
         zd(ngrid,nsoilmx-1),      &

         FluxRad(ngrid),           & ! radiative flux at surface
         FluxGrd(ngrid),           & ! heat flux from deep soil
         capcal(ngrid),            & ! effective heat capacity of soil

         zdUfr(ngrid,nlayer),      & ! partial tendencies for zonal velocity,
         zdVfr(ngrid,nlayer),      & ! partial tendencies for meridional velocity,
         zdHfr(ngrid,nlayer),      & ! partial tendencies for potential temperature,

         zdTsrfr(ngrid),           & ! surface temperature
         zdTsrf(ngrid),            & ! total tendency of surface temperature

         zflubid(ngrid),           & ! radiative + deep soil fluxes
         zPmer(ngrid)                ! sea-level pressure

    integer                       :: j, l, ig, igout
    real                          :: zday, zdtime, z1, z2
    real, dimension(ngrid,nlayer) :: zdum1, zdum2, zdum3
    logical                       :: lwrite

    call nvtxstartrange ("physics")

    if (ngrid /= ngridmax) then
       print*,'STOP in inifis'
       print*,'Probleme de dimensions :'
       print*,'ngrid     = ', ngrid
       print*,'ngridmax  = ', ngridmax
       stop
    end if

    igout = ngrid/2 + 1
    zday  = rjourvrai + gmtime


    !-----------------------------------------------------------------------
    !    0. Allocate and initialize at first call
    !-----------------------------------------------------------------------
    if (firstcall) then
       zday_last = zday - ptimestep / unjours

       ! Allocate only if very first call of entire run, especially in single column cases
       if (firstcall_alloc) call alloc (ngrid, nlayer)
       call precompute
       call coldstart (ngrid)
       firstcall_alloc = .false. ! turn flag false, such that allocation occurs once
    end if

    call profile_enter (id_phyparam)


    !-----------------------------------------------------------------------
    !    1. Initialisations
    !-----------------------------------------------------------------------
    icount = icount + 1

    lwrite = .false.

    pdV     = 0.0
    pdU     = 0.0
    pdT     = 0.0
    pdPsrf  = 0.0
    zflubid = 0.0
    zdTsrf  = 0.0

    ! Calcul du geopotentiel aux niveaux intercouches ponderation des altitudes au niveau des couches en dP/P
    zZlay = pPhi / g
    zZlev(:,1) = 0.0

    do l = 2, nlayer
       do ig = 1, ngrid
          z1 = (pPlay(ig,l-1) + pPlev(ig,l)) / (pPlay(ig,l-1) - pPlev(ig,l))
          z2 = (pPlev(ig,l)   + pPlay(ig,l)) / (pPlev(ig,l)   - pPlay(ig,l))

          zZlev(ig,l) = (z1 * zZlay(ig,l-1) + z2 * zZlay(ig,l)) / (z1 + z2)
       end do
    end do

    ! Transformation de la temperature en temperature potentielle
    zpopsk = (pPlay / pPlev(:,1:nlayer))**Rcp ! surface pressure is used as reference pressure
    zH     = pT / zpopsk


    !-----------------------------------------------------------------------------------------------
    !  2. Vertical diffusion of heat in soil column
    !
    !  First half split-step of implicit time integration forward sweep from
    !  deep ground to surface.
    !
    !  Returns Lu coefficients zc, zd and CapCal, FluxGrd
    !
    !-----------------------------------------------------------------------------------------------
    if (callsoil) then
       call soil_forward (ngrid, nsoilmx, pTimestep, inertie, Tsurf, Tsoil, zc, zd, CapCal, FluxGrd)
    else
       CapCal  = CapCal_nosoil
       FluxGrd = 0.0
    end if

    if (lverbose) then
       WRITELOG (*,*) 'Surface heat capacity, conduction flux, Ts'
       WRITELOG (*,*) CapCal(igout), FluxGrd(igout), Tsurf(igout)
       LOG_DBG ('phyparam')
    end if


    !-----------------------------------------------------------------------------------------------
    !    2. Radiative tendencies
    !-----------------------------------------------------------------------------------------------
    if (callrad) &
         call radiative_tendencies (lwrite, ngrid, igout, nlayer, gmtime, pTimestep * float(iradia), zDay, pPlev, pPlay, pT,  &
         pdT, FluxRad)


    !-----------------------------------------------------------------------------------------------
    !    3. Vertical diffusion (turbulent mixing)
    !
    ! Second-order Strang splitting consisting of:
    !
    ! 1. Explicit Forward Euler  split-step for physics tendencies computed above (soil and radiative)
    !
    ! 2. Implicit Backwards Euler split-step for vertical turbulent diffusion:
    ! kz is computed and then vertical diffusion is integrated in time implicitly
    ! using a linear relationship between surface heat flux and air temperature
    ! in lowest level (Robin-type BC).
    !
    !-----------------------------------------------------------------------------------------------
    if (calldifv) then
       zflubid = FluxRad + FluxGrd

       zdum1 = 0.0
       zdum2 = 0.0
       zdum3 = pdT / zpopsk

       call Vdif (                        &
            ngrid, nlayer, zday,          &
            ptimestep, capcal, z0,        &
            pPlay, pPlev, zZlay, zZlev,   &
            pU, pV, zH, Tsurf, Emissiv,   &
            zdum1, zdum2, zdum3, zflubid, &
            zdUfr, zdVfr, zdHfr, zdTsrfr, &
            lverbose)

       pdV = pdV + zdVfr                              ! zonal velocity
       pdU = pdU + zdUfr                              ! meridional velocity

       pdT    = pdT    + zdHfr * zpopsk               ! temperature
       zdTsrf = zdTsrf + zdTsrfr                      ! surface temperature
    else ! no vertical diffusion
       zdTsrf = zdTsrf + (FluxRad + FluxGrd) / CapCal ! surface temperature
    end if

    !-------------------------------------------------------------
    !   Soil temperatures : 2nd half of implicit time integration
    !   using updated Tsurf as input
    !-------------------------------------------------------------
    Tsurf = Tsurf + ptimestep * zdTsrf

    if (callsoil) then
       call soil_backward (ngrid, nsoilmx, zc, zd, Tsurf, Tsoil)
       if (lverbose) then
          WRITELOG (*,*) 'surface ts, dts, dt'
          WRITELOG (*,*) Tsurf(igout), zdTsrf(igout), ptimestep
          LOG_DBG ('phyparam')
       end if
    end if


    !-----------------------------------------------------------------------
    !   4. Dry convective adjustment
    !-----------------------------------------------------------------------
    if (calladj) then
       zdum1 = pdT / zpopsk

       zdUfr = 0.0
       zdVfr = 0.0
       zdHfr = 0.0

       call convadj (                 &
            ngrid, nlayer, ptimestep, &
            pPlay, pPlev, zpopsk,     &
            pU, pV, zH,               &
            pdU, pdV, zdum1,          &
            zdUfr, zdVfr, zdHfr)

       pdU = pdU + zdUfr
       pdV = pdV + zdVfr
       pdT = pdT + zdHfr * zpopsk
    end if


    !-----------------------------------------------------------------------
    !   Sorties
    !-----------------------------------------------------------------------
    WRITELOG (*,*) 'zday, zday_last ', zday, zday_last, icount
    LOG_DBG ('phyparam')

    if (lwrite) then
       zPmer = pPlev(:,1) * exp (pPhi(:,1) / (r * ref_temp))

       call writefield ('U',     'vent zonal moy',        'm/s',   pU)
       call writefield ('V',     'vent meridien moy',     'm/s',   pV)
       call writefield ('temp',  'temperature',           'K',     pT)
       call writefield ('theta', 'potential temperature', 'K',     zH)
       call writefield ('geop',  'geopotential',          'm2/s2', pPhi)
       call writefield ('plev',  'plev',                  'Pa' ,   pPlev(:,1:nlayer))

       call writefield ('dU', 'dU',' ', pdU)
       call writefield ('dV', 'dV',' ', pdV)
       call writefield ('dT', 'dT',' ', pdT)

       call writefield ('Ts',     'Surface temperature', 'K',  Tsurf)
       call writefield ('coslon', 'coslon',              ' ',  coslon)
       call writefield ('sinlon', 'sinlon',              ' ',  sinlon)
       call writefield ('coslat', 'coslat',              ' ',  coslat)
       call writefield ('sinlat', 'sinlat',              ' ',  sinlat)
       call writefield ('alb',    'Albedo',              ' ',  Albedo)
       call writefield ('Ps',     'Surface pressure',    'Pa', pPlev(:,1))
       call writefield ('slp',    'Sea level pressure',  'Pa', zPmer)
    end if

    call profile_exit (id_phyparam)
    call nvtxendrange
  end subroutine phyparam

  subroutine alloc (ngrid, nlayer) bind(c, name='phyparam_alloc')
    use astronomy, only : iniorbit
    use soil_mod,  only : Tsurf, Tsoil, z0, inertie, rnatur, Albedo, Emissiv
    integer, intent(in), value :: ngrid, nlayer

    ! Allocate precomputed arrays
    allocate (rnatur(ngrid), Albedo(ngrid), Emissiv(ngrid))
    allocate (z0(ngrid),inertie(ngrid))

    ! Allocate arrays for internal state
    allocate (Tsurf(ngrid))
    allocate (Tsoil(ngrid, nsoilmx))

    call iniorbit
  end subroutine alloc

  subroutine precompute() bind(c, name='phyparam_precompute')
    ! Precompute time-independent arrays
    use soil_mod, only: rnatur, inertie, z0, Emissiv, Albedo,  i_mer, i_ter, Cd_mer, Cd_ter,  alb_mer, alb_ter, emi_mer, emi_ter

    rnatur(:)  = 1.0
    inertie(:) = (1.0 - rnatur) * i_mer   + rnatur * i_ter
    z0(:)      = (1.0 - rnatur) * Cd_mer  + rnatur * Cd_ter
    Emissiv(:) = (1.0 - rnatur) * emi_mer + rnatur * emi_ter
    Albedo(:)  = (1.0 - rnatur) * alb_mer + rnatur * alb_ter
  end subroutine precompute

  subroutine coldstart(ngrid) bind(c, name='phyparam_coldstart')
    ! Create internal state to start a run without a restart file
    use soil_mod, only : Tsurf, Tsoil
    integer, intent(in), value :: ngrid

    Tsurf  = Tsoil_init
    Tsoil  = Tsoil_init
    icount = 0
    if (.not. callsoil .and. firstcall_alloc .and. lverbose) then ! write only on first call (lverbose set by interface)
       WRITELOG (*,*) 'WARNING!! thermal conduction in the soil turned off'
       LOG_WARN ('coldstart')
    end if
  end subroutine coldstart
end module phyparam_mod
