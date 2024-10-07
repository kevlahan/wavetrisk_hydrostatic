module phyparam_mod
  use callkeys
  use comgeomfi
  implicit none
  private
  save
  integer         :: icount
  real            :: zday_last
  logical         :: firstcall_alloc = .true.

  real, parameter :: CapCal_nosoil = 1e5
  real, parameter :: ref_temp      = 285.0
  real, parameter :: Tsoil_init    = 300.0

  public :: phyparam, alloc, precompute, zDay_last, icount
contains
  subroutine phyparam (ngrid, nlayer, mask, firstcall, lastcall, rJourvrai, gmTime, pTimestep,  pPlev, pPlay, &
       pPhi, pPhi_surf, pU, pV, pT) bind (c, name='phyparam_phyparam')

    use,          intrinsic :: iso_c_binding
    use phys_const,     only : g, Rcp, r, unjours
    use soil_mod,       only : soil_forward, soil_backward
    use soil_mod,       only : Z0, pThermal_inertia, Emissiv, Albedo  ! precomputed
    use soil_mod,       only : Tsurf, Tsoil                           ! surface temperature, soil column temperatures
    use turbulence,     only : Vdif, p_0
    use convection,     only : ConvAdj
    USE radiative_mod,  only : radiative_tendencies
    use writefield_mod, only : writefield
    use accelerator,    only : nvtxstartrange, nvtxendrange
    use profiling,      only : profile_enter, profile_exit, id_phyparam
    !==============================================================================================
    !
    !   Physical parametrisations of LMDZ  20 parameter climate model for planetary atmospheres
    !
    !   Includes following sub-models:
    !
    !   1. radiative transfer (long and shortwave) for CO2 and dust
    !   3. convective adjustment
    !   4. vertical heat diffusion in soil column
    !   2. vertical turbulent mixing
    !
    !   Author:  Frederic Hourdin  1993-10-15
    !   Revised: Nicholas Kevlahan 2024-10-04
    !
    !==============================================================================================

    ! Input variables
    integer, value,                  intent(in)    :: ngrid     ! number of columns
    integer, value,                  intent(in)    :: nlayer    ! number of vertical layers
    integer, value,                  intent(in)    :: mask      ! adaptive grid mask value
    real,    value,                  intent(in)    :: rJourvrai ! number of days, counted from northern spring equinox
    real,    value,                  intent(in)    :: gmTime    ! fraction of  day (ranges from 0 to 1)
    real,    value,                  intent(in)    :: pTimestep ! timestep [s]
    real, dimension(ngrid,nlayer),   intent(in)    :: pPlay     ! pressure at layer centres    [Pa]
    real, dimension(ngrid,nlayer+1), intent(in)    :: pPlev     ! pressure at layer interfaces [Pa]
    real, dimension(ngrid,nlayer),   intent(in)    :: pPhi      ! geopotential atlayer centres [m^2/s^2]
    real, dimension(ngrid),          intent(in)    :: pPhi_surf ! surface geopotential         [m^2/s^2]
    logical(kind=c_bool), value,     intent(in)    :: firstcall ! true at the first call
    logical(kind=c_bool), value,     intent(in)    :: lastcall  ! true at the last call

    ! Ouput: velocities and temperature from implicit Euler step: t -> t+dt
    real, dimension(ngrid,nlayer),   intent(inout) :: pU        ! zonal velocity      [m/s]
    real, dimension(ngrid,nlayer),   intent(inout) :: pV        ! meridional velocity [m/s]
    real, dimension(ngrid,nlayer),   intent(inout) :: pT        ! temperature         [K]

    ! Local variables
    real                                           :: zDay       ! day
    real, dimension(ngrid)                         :: zdTsrf     ! total tendency of surface temperature     [K/s]
    real, dimension(ngrid)                         :: zdTsrfr    ! intermediate surface temperature tendency [K/s]
    real, dimension(ngrid)                         :: zFluxId    ! surface flux
    real, dimension(ngrid)                         :: FluxRad    ! radiative flux at surface
    real, dimension(ngrid)                         :: FluxGrd    ! heat flux from deep soil
    real, dimension(ngrid)                         :: CapCal     ! effective heat capacity of soil
    real, dimension(ngrid)                         :: zPmer      ! sea-level pressure
    real, dimension(ngrid,nlayer)                  :: pH         ! potential temperature [K]
    real, dimension(ngrid,nlayer)                  :: zExner     ! Exner function
    real, dimension(ngrid,nlayer)                  :: zZlay      ! height at layer centres    [m]
    real, dimension(ngrid,2:nlayer)                :: z1, z2
    real, dimension(ngrid,nlayer+1)                :: zZlev      ! height at layer interfaces [m]
    real, dimension(:,:), allocatable              :: zc, zd     ! LU coefficients for soil implicit solve

    call nvtxstartrange ("physics")

    if (ngrid /= ngridmax) then
       print*,'STOP in inifis'
       print*,'Probleme de dimensions :'
       print*,'ngrid     = ', ngrid
       print*,'ngridmax  = ', ngridmax
       stop
    end if

    zDay = rJourvrai + gmTime

    ! Allocate and initialize at first call
    if (firstcall) then
       zDay_last = zDay - pTimestep / unjours

       ! Allocate only if very first call of entire run, especially in single column cases
       if (firstcall_alloc) call alloc (ngrid, nlayer)
       call precompute
       call coldstart (ngrid)
       firstcall_alloc = .false. ! turn flag false, such that allocation occurs once
    end if

    call profile_enter (id_phyparam)

    ! Initialisations
    icount  = icount + 1
    zFluxId = 0.0

    ! Height at middle of layers
    zZlay = pPhi / g

    ! Height at layer interfaces
    zZlev(:,1) = pPhi_surf / g ! surface height

    z1 = (pPlay(:,1:nlayer-1) + pPlev(:,2:nlayer)) / (pPlay(:,1:nlayer-1) - pPlev(:,2:nlayer))
    z2 = (pPlev(:,2:nlayer)   + pPlay(:,2:nlayer)) / (pPlev(:,2:nlayer)   - pPlay(:,2:nlayer))

    zZlev(:,2:nlayer) = (z1 * zZlay(:,1:nlayer-1) + z2 * zZlay(:,2:nlayer)) / (z1 + z2)
    
    ! Exner function for converting between temperature and potential temperature
    zExner = (pPlay/p_0)**Rcp

    !  Soil temperatures
    !  First split step of implicit time integration forward sweep from deep ground to surface
    !  (returns LU coefficients zc, zd and CapCal, FluxGrd)
    if (callsoil) then
       allocate (zc(ngrid,nsoilmx-1), zd(ngrid,nsoilmx-1))

       call soil_forward (ngrid, nsoilmx, pTimestep, pThermal_inertia, Tsurf, Tsoil, zc, zd, CapCal, FluxGrd)
    else ! no diffusion of heat in soil column
       CapCal  = CapCal_nosoil
       FluxGrd = 0.0
    end if

    !  Radiative transfer split step
    if (callrad) call radiative_tendencies (ngrid, nlayer, gmTime, pTimestep, zDay, pPlev, pPlay, pT, FluxRad)

    ! Find potential temperature from temperature
    pH = pT / zExner

    ! Vertical turbulent diffusion split step
    !
    ! Second-order Strang splitting consisting of:
    !
    ! 1. Explicit Forward Euler  split-step for physics tendencies computed above (soil and radiative)
    ! 2. Implicit Backwards Euler split-step for vertical turbulent diffusion.
    !
    ! Turbulent diffusivity Kz is computed, and then vertical diffusion is integrated in time implicitly using
    ! a linear relationship between surface heat flux and air temperature in lowest level (Robin-type BC).
    if (calldifv) then
       zFluxId = FluxRad + FluxGrd

       call Vdif (ngrid, nlayer, mask, zDay, pTimestep, CapCal, Emissiv, zFluxId, Z0, pPlay, pPlev, zZlay, zZlev, &
            pU, pV, pH, Tsurf)
    else ! no vertical turbulent diffusion
       Tsurf = Tsurf + pTimestep * (FluxRad + FluxGrd) / CapCal
    end if

    ! Soil column temperatures at t+dt
    ! (second split step of implicit time integration using updated Tsurf at d+dt as input)
    if (callsoil) call soil_backward (ngrid, nsoilmx, zc, zd, Tsurf, Tsoil)

    ! Dry convective adjustment
    ! (final values of velocities and potential temperature at t+dt)
    if (calladj) call ConvAdj (ngrid, nlayer, pTimestep, pPlay, pPlev, zExner, pU, pV, pH)

    ! Temperature at t+dt for output
    pT = pH * zExner

    call profile_exit (id_phyparam)
    call nvtxendrange
  end subroutine phyparam

  subroutine alloc (ngrid, nlayer) bind (c, name='phyparam_alloc')
    use astronomy, only : iniorbit
    use soil_mod,  only : Tsurf, Tsoil, Z0, pThermal_inertia, land, Albedo, Emissiv
    integer, intent(in), value :: ngrid, nlayer

    ! Allocate precomputed arrays
    allocate (land(ngrid), Albedo(ngrid), Emissiv(ngrid))
    allocate (Z0(ngrid), pThermal_inertia(ngrid))

    ! Allocate arrays for internal state
    allocate (Tsurf(ngrid))
    if (callsoil) allocate (Tsoil(ngrid,nsoilmx))

    call iniorbit
  end subroutine alloc

  subroutine precompute () bind (c, name='phyparam_precompute')
    ! Precompute time-independent arrays
    use soil_mod, only: land, pThermal_inertia, Z0, Emissiv, Albedo,  I_mer, I_ter, Cd_mer, Cd_ter,  Alb_mer, Alb_ter, Emi_mer, Emi_ter

    land             = 1.0                                      ! proportion of land compared to sea, 0 <= land <= 1

    pThermal_inertia = (1.0 - land) * I_mer   + land * I_ter    ! thermal inertia
    Z0               = (1.0 - land) * Cd_mer  + land * Cd_ter   ! roughness length
    Emissiv          = (1.0 - land) * Emi_mer + land * Emi_ter  ! emissivity
    Albedo           = (1.0 - land) * Alb_mer + land * Alb_ter  ! zlbedo
  end subroutine precompute

  subroutine coldstart (ngrid) bind (c, name='phyparam_coldstart')
    ! Create internal state to start a run without a restart file
    use soil_mod, only : Tsurf, Tsoil
    integer, intent(in), value :: ngrid

    Tsurf = Tsoil_init
    if (callsoil) Tsoil  = Tsoil_init
    icount = 0
  end subroutine coldstart
end module phyparam_mod
