!-------------------------------------------------------------------
!
!  Author: Gabrielle Ching-Johnson
!  Revised: 2023-05-18
!  Revised: Nicholas Kevlahan 2024-10-01
!
!     module containing functions to be utilized when the dynamics
!     is sending single columns to the physics.
!
!-------------------------------------------------------------------
module single_column_mod
  use,         intrinsic :: iso_c_binding
  use comgeomfi,    only :  nsoilmx, ngridmax, long, lati, sinlon, coslon, sinlat, coslat
  use callkeys,     only :  soil_flag => callsoil
  use phyparam_mod, only :  phyparam
  use soil_mod,     only :  Tsurf, Tsoil ! surface temperature and soil layer temperatures updated in simple physics
  implicit none
  private
  integer :: extra_levels = 1   ! default value (surface layer only)
  public :: initialize_extra_levels, get_extra_levels, physics_call_single_col, change_latitude_longitude
contains
  subroutine physics_call_single_col (ngrid, nlayer, mask, firstcall, lastcall, rJourvrai, &
       gmTime, pTimestep, pPlev, pPlay, pPhi, pPhi_surf, pU, pV, pTheta, Tsurf_soil)
    !----------------------------------------------------------------
    !
    !   WrapPer routine dynamics will use to call the physics,
    !   -----------  for a single column. it updates the surface and soil temps
    !   -----------  for the column before the call and send back the newly
    !   -----------  update temperatures.
    !
    !   extra notes: at the 1st call for each column, the dynamics wont know
    !                what the soil and surface is set to, can send random nums
    !                as phyparam() does cal coldstart() call which sets surface temperature Tsurf
    !                and soil column temperatures Tsoil to 300 K.
    !
    !
    !   author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    ! Input:
    integer, value,                      intent(in)    :: ngrid      ! number of columns
    integer, value,                      intent(in)    :: nlayer     ! number of vertical layers in atmosphere
    integer, value,                      intent(in)    :: mask       ! adaptive grid mask value
    real,    value,                      intent(in)    :: rJourvrai  ! number of days counted from the north. spring equinox
    real,    value,                      intent(in)    :: gmTime     ! fraction of the day (ranges from 0 to 1)
    real,    value,                      intent(in)    :: pTimestep  ! timestep [s]
    real, dimension(ngrid,nlayer+1),     intent(in)    :: pPlev      ! pressure at interfaces [Pa]
    real, dimension(ngrid,nlayer),       intent(in)    :: pPlay      ! pressure at layers     [Pa]
    real, dimension(ngrid,nlayer),       intent(in)    :: pPhi       ! geopotential at layers [m^2/s^2]
    real, dimension(ngrid),              intent(in)    :: pPhi_surf  ! surface geopotential   [m^2/s^2]
    logical(kind=c_bool), value,         intent(in)    :: firstcall  ! true at first call
    logical(kind=c_bool), value,         intent(in)    :: lastcall   ! true at last call

    real, dimension(ngrid,nlayer),       intent(inout) :: pU         ! zonal velocity         [m/s]
    real, dimension(ngrid,nlayer),       intent(inout) :: pV         ! meridional velocity    [m/s]
    real, dimension(ngrid,nlayer),       intent(inout) :: pTheta     ! potential temperature  [K]

    real, dimension(ngrid,extra_levels), intent(inout) :: Tsurf_soil ! temperature for surface and soil layers [K]

    ! Set current physics soil layers temperature from dynamics (shared with soil_mod)
    if (.not. firstcall) then
       Tsurf = Tsurf_soil(ngrid,1) ! surface temperature
       if (soil_flag) Tsoil = Tsurf_soil(:,2:extra_levels)
    end if

    ! Call simple physics for this column
    call phyparam (ngrid, nlayer, mask, firstcall, lastcall, rJourvrai, gmTime, pTimestep, &
         pPlev, pPlay, pPhi, pPhi_surf, pU, pV, pTheta)

    ! Update with physics surface temperature and soil column termperatures (from soil_mod)
    Tsurf_soil(:,1) = Tsurf                             ! surface temperature
    if (soil_flag) Tsurf_soil(:,2:extra_levels) = Tsoil ! soil column temperatures
  end subroutine physics_call_single_col

  subroutine initialize_extra_levels (levels)
    !----------------------------------------------------------------
    !
    !   Sets  value of private variable extra_levels
    !
    !   levels = nsoilmx + 1 if soil model is turned on
    !   levels = 1           if soil model is not used
    !
    !   author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------
    integer :: levels

    ! Check to ensure the physics was set for single column calls
    if (ngridmax /= 1) then
       print*, 'STOP initialization for single column call in physics: max columns in physics and ngrid should both be set to 1'
       stop
    end if

    if (soil_flag) then ! surface layer + soil column layers
       extra_levels = nsoilmx + 1
    else                ! surface layer only
       extra_levels = 1
    end if

    ! Check to ensure the number of levels set in dynamics is same as physics
    if (levels /= extra_levels) then
       print*, 'STOP!! the number of extra levels set in physics is not same in dynamics'
       print*, 'the levels set in the physics is', extra_levels
       print*, 'the levels set in dynamics is', levels
       stop
    end if

  end subroutine initialize_extra_levels

  integer function get_extra_levels ()
    !----------------------------------------------------------------
    !
    !  Retreive value of extra levels.
    !
    !   author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    get_extra_levels = extra_levels
  end function get_extra_levels

  subroutine change_latitude_longitude (latitude, longitude)
    !----------------------------------------------------------------
    !
    !   Changes  value of latitude and longitude in the case of single
    !   -----------   column for max number of grid points (ngridmax)
    !
    !   author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------
    real, intent(in) :: latitude, longitude ! in radians

    ! Need to ensure that only a single column, either a check for a flag, or size of ngridmax
    long(1) = longitude
    lati(1) = latitude

    sinlat(1) = sin (lati(1))
    coslat(1) = cos (lati(1))
    sinlon(1) = sin (long(1))
    coslon(1) = cos (long(1))
  end subroutine change_latitude_longitude
end module single_column_mod
