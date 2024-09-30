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
  subroutine phyparam (ngrid, nlayer, firstcall, lastcall, rjourvrai, gmtime, pTimestep, &
       pPlev, pPlay, pPhi, pU, pV, pT, &
       pdU, pdV, pdt, pdPsrf) bind (c, name='phyparam_phyparam')

    use,          intrinsic :: iso_c_binding
    use phys_const,     only : g, Rcp, r, unjours
    use soil_mod,       only : soil_forward, soil_backward
    use soil_mod,       only : z0, pThermal_inertia, Emissiv, Albedo  ! precomputed
    use soil_mod,       only : Tsurf, Tsoil                           ! state variables
    use turbulence,     only : Vdif
    use convection,     only : ConvAdj
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
         zExner(ngrid,nlayer),     & ! Exner function

         zZlev(ngrid,nlayer+1),    & ! height of layer interfaces
         zZlay(ngrid,nlayer),      & ! height of layer centres

         FluxRad(ngrid),           & ! radiative flux at surface
         FluxGrd(ngrid),           & ! heat flux from deep soil
         capcal(ngrid),            & ! effective heat capacity of soil

         zdUfr(ngrid,nlayer),      & ! partial tendencies for zonal velocity,
         zdVfr(ngrid,nlayer),      & ! partial tendencies for meridional velocity,
         zdHfr(ngrid,nlayer),      & ! partial tendencies for potential temperature,

         zdTsrf(ngrid),            & ! total tendency of surface temperature
         zdTsrfr(ngrid),           & ! surface temperature

         zFluxid(ngrid),           & ! radiative + deep soil fluxes
         zPmer(ngrid)                ! sea-level pressure


    real, dimension(:,:), allocatable ::   zc, zd ! Lu coefficients for soil implicit solve

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
    zFluxid = 0.0
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
    zExner = (pPlay / pPlev(:,1:nlayer))**Rcp ! surface pressure is used as reference pressure
    zH     = pT / zExner


    !-----------------------------------------------------------------------------------------------
    !  2. Soil temperatures
    !
    !  First split-step of implicit time integration forward sweep from deep ground to surface.
    !
    !  Returns Lu coefficients zc, zd and CapCal, FluxGrd
    !
    !-----------------------------------------------------------------------------------------------
    if (callsoil) then
       allocate (zc(ngrid,nsoilmx-1), zd(ngrid,nsoilmx-1))

       call soil_forward (ngrid, nsoilmx, pTimestep, pThermal_inertia, Tsurf, Tsoil, zc, zd, CapCal, FluxGrd)
    else ! no diffusion of heat in soil column
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
    !    3. Vertical turbulent diffusion
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
       zFluxid = FluxRad + FluxGrd

       zdum1 = 0.0
       zdum2 = 0.0
       zdum3 = pdT / zExner

       call Vdif (ngrid, nlayer, zday, pTimestep, CapCal, z0, pPlay, pPlev, zZlay, zZlev, pU, pV, zH, Tsurf, Emissiv, &
            zdum1, zdum2, zdum3, zFluxid, &
            zdUfr, zdVfr, zdHfr, zdTsrfr, &
            lverbose)

       ! Update pseudo-tendencies
       pdV    = pdV    + zdVfr           ! zonal velocity
       pdU    = pdU    + zdUfr           ! meridional velocity
       pdT    = pdT    + zdHfr * zExner  ! temperature from potential temperature tendency
       zdTsrf = zdTsrf + zdTsrfr         ! surface temperature
    else ! no vertical turbulent diffusion
       zdTsrf = zdTsrf + (FluxRad + FluxGrd) / CapCal
    end if

    ! Surface temperature adjusted by vertical turbulent diffusion
    Tsurf = Tsurf + ptimestep * zdTsrf


    !-----------------------------------------------------------------------------------
    !   4. Soil temperatures
    !
    !   Second split step of implicit time integration using updated Tsurf as input
    !
    !-----------------------------------------------------------------------------------
    if (callsoil) then
       call soil_backward (ngrid, nsoilmx, zc, zd, Tsurf, Tsoil)

       if (lverbose) then
          WRITELOG (*,*) 'Surface Ts, dTs, dt'
          WRITELOG (*,*) Tsurf(igout), zdTsrf(igout), pTimestep
          LOG_DBG ('phyparam')
       end if
    end if


    !-----------------------------------------------------------------------
    !   4. Dry convective adjustment
    !-----------------------------------------------------------------------
    if (calladj) then
       zdum1 = pdT / zExner

       call ConvAdj (ngrid, nlayer, ptimestep, pPlay, pPlev, zExner, pU, pV, zH, pdU, pdV, zdum1,  zdUfr, zdVfr, zdHfr)

       ! Update pseudo-tendencies
       pdU = pdU + zdUfr
       pdV = pdV + zdVfr
       pdT = pdT + zdHfr * zExner
    end if


    !-----------------------------------------------------------------------
    !   Data logging
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

  subroutine alloc (ngrid, nlayer) bind (c, name='phyparam_alloc')
    use astronomy, only : iniorbit
    use soil_mod,  only : Tsurf, Tsoil, z0, pThermal_inertia, Rnatur, Albedo, Emissiv
    integer, intent(in), value :: ngrid, nlayer

    ! Allocate precomputed arrays
    allocate (Rnatur(ngrid), Albedo(ngrid), Emissiv(ngrid))
    allocate (z0(ngrid), pThermal_inertia(ngrid))

    ! Allocate arrays for internal state
    allocate (Tsurf(ngrid))
    if (callsoil) allocate (Tsoil(ngrid,nsoilmx))

    call iniorbit
  end subroutine alloc

  subroutine precompute () bind (c, name='phyparam_precompute')
    ! Precompute time-independent arrays
    use soil_mod, only: Rnatur, pThermal_inertia, z0, Emissiv, Albedo,  I_mer, I_ter, Cd_mer, Cd_ter,  Alb_mer, Alb_ter, Emi_mer, Emi_ter

    Rnatur           = 1.0
    pThermal_inertia = (1.0 - Rnatur) * I_mer   + Rnatur * i_ter
    z0               = (1.0 - Rnatur) * Cd_mer  + Rnatur * Cd_ter
    Emissiv          = (1.0 - Rnatur) * Emi_mer + Rnatur * Emi_ter
    Albedo           = (1.0 - Rnatur) * Alb_mer + Rnatur * Alb_ter
  end subroutine precompute

  subroutine coldstart (ngrid) bind (c, name='phyparam_coldstart')
    ! Create internal state to start a run without a restart file
    use soil_mod, only : Tsurf, Tsoil
    integer, intent(in), value :: ngrid

    Tsurf = Tsoil_init
    if (callsoil) Tsoil  = Tsoil_init
    icount = 0
    if (.not. callsoil .and. firstcall_alloc .and. lverbose) then ! write only on first call (lverbose set by interface)
       WRITELOG (*,*) 'WARNING!! thermal conduction in the soil turned off'
       LOG_WARN ('coldstart')
    end if
  end subroutine coldstart
end module phyparam_mod
