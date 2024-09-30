!-------------------------------------------------------------------
!
!  file: single_column_mod.f90
!  author: gabrielle ching-johnson
!  date revised: may 18th 2023
!
!  description:
!     module containing functions to be utilized when the dynamics
!     will  be sending single columns to the physics.
!
!-------------------------------------------------------------------

module single_column_mod
  use,         intrinsic :: iso_c_binding
  use comgeomfi,    only :  nsoilmx, ngridmax, long, lati, sinlon, coslon, sinlat, coslat
  use callkeys,     only :  soil_flag => callsoil
  use phyparam_mod, only :  phyparam
  use soil_mod,     only :  tsurf, tsoil
  implicit none
  private
  integer :: extra_levels = 1   ! soil levels + surface level
  public :: initialize_extra_levels, get_extra_levels, physics_call_single_col, change_latitude_longitude
contains

  subroutine initialize_extra_levels (levels)
    !----------------------------------------------------------------
    !
    !   Sets  value of private variable extra_levels
    !
    !   levels = nsoilmx + 1 if soil model is turned on
    !   levels = 1           if soil model is not used
    !
    !   author: gabrielle ching-johnson
    !
    !----------------------------------------------------------------
    integer :: levels

    ! Check to ensure the physics was set for single column calls
    if (ngridmax /= 1) then
       print*, 'STOP initialization for single column call in physics: max columns in physics and ngrid should both be set to 1'
       stop
    end if

    if (soil_flag) then
       extra_levels = nsoilmx + 1
    else
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

  subroutine physics_call_single_col (ngrid, nlayer, firstcall, lastcall, rjourvrai, &
       gmtime, ptimestep, pplev, pplay, pphi, pu, pv, pt, extra_temp, &
       pdu, pdv, pdt, pdpsrf)
    !----------------------------------------------------------------
    !
    !   Wrapper routine dynamics will use to call the physics,
    !   -----------  for a single column. it updates the surface and soil temps
    !   -----------  for the column before the call and send back the newly
    !   -----------  update temperatures.
    !
    !   extra notes: at the 1st call for each column, the dynamics wont know
    !                what the soil and surface is set to, can send random nums
    !                as phyparam() does cal coldstart() call which sets the tsurf
    !                and tsoil set to 300k.
    !
    !
    !   author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    ! Input:
    integer,                         intent(in), value ::  &
         ngrid,         & ! number of columns
         nlayer           ! number of vertical layers
    real,                    intent(in), value :: &
         rjourvrai,     & ! number of days counted from the north. spring equinox
         gmtime,        & ! fraction of the day (ranges from 0 to 1)
         ptimestep        ! timestep [s]

    real, dimension(ngrid,nlayer+1), intent(in) :: pplev ! pressure at  layer interfaces [Pa]

    real, dimension(ngrid,nlayer),   intent(in) :: &
         pplay,         & ! pressure at the middle of the layers (pa)
         pphi,          & ! geopotential at the middle of the layers (m2s-2)
         pu,            & ! u component of the wind (ms-1)
         pv,            & ! v component of the wind (ms-1)
         pt               ! temperature (k)

    real, dimension(ngrid,extra_levels), intent(inout) :: extra_temp  ! surface and soil temp (k) of the column

    ! Output : physical tendencies
    real, dimension(ngrid,nlayer), intent(out) ::   &
         pdu,           &                                 ! tendency of zonal velocity      [m/s^2]
         pdv,           &                                 ! tendency of meridional velocity [m/s^2]
         pdt                                              ! tendency of temperature         [K/s]
    real, dimension(ngrid),        intent(out) ::  pdpsrf ! tendency of surface pressure    [Pa/s]

    logical(kind=c_bool),   intent(in), value :: &
         firstcall,                   & ! true at first call
         lastcall                       ! true at last call

    ! Check if first call or not, if not can't update tsoil/tsurf as they haven't be allocated: will be initialized in phyparm
    if (.not. firstcall) then
       if (soil_flag) tsoil = extra_temp(:, 2:) ! update physics soil temp if turned on
       tsurf = extra_temp(ngrid,1)              ! update physics surface temp
    end if

    ! call phyparam physics
    call phyparam (ngrid, nlayer, firstcall, lastcall, rjourvrai, gmtime, ptimestep, pplev, pplay, pphi, pu, pv, pt, &
         pdu, pdv, pdt, pdpsrf)

    ! update output extra_temp with physics soil and surface temp
    extra_temp(ngrid,1) = tsurf(ngrid)

    if (soil_flag) extra_temp(:,2:) = tsoil
  end subroutine physics_call_single_col
end module single_column_mod
