!-------------------------------------------------------------------
!
!  File: single_column_mod.F90
!  Author: Gabrielle Ching-Johnson
!  Date Revised: May 18th 2023
!
!  Description:
!     Module containing functions to be utilized when the dynamics
!     will  be sending single columns to the physics.
!
!-------------------------------------------------------------------

MODULE single_column_mod

  USE, INTRINSIC :: ISO_C_BINDING
  USE comgeomfi, ONLY : nsoilmx, ngridmax, long, lati, sinlon, coslon, sinlat, coslat
  USE callkeys, ONLY : soil_flag => callsoil
  USE phyparam_mod, ONLY : phyparam
  USE soil_mod, ONLY : tsurf, tsoil

  IMPLICIT NONE
  PRIVATE

  INTEGER :: extra_levels = 1   ! states how many soil levels + surface level
  PUBLIC :: initialize_extra_levels, get_extra_levels, physics_call_single_col, change_latitude_longitude

CONTAINS

  SUBROUTINE initialize_extra_levels(levels)
    !----------------------------------------------------------------
    !
    !   Description: Sets the value of private variable extra_levels
    !   ----------   to number soil levels (nsoilmx) + 1 if soil model
    !   ----------   is turned on, otherwise to 1.
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------
    INTEGER :: levels

    ! Check to ensure the physics was set for single column calls
    IF (ngridmax .NE. 1) THEN
       PRINT*, 'STOP Initialization for single column call in physics'
       PRINT*, 'Check for number of max columns in physics to be set to 1 FALSE'
       PRINT*, 'Check previous physics initialization call to inicomgeofi that ngrid is 1'
       STOP
    END IF

    IF (soil_flag) THEN
       extra_levels = nsoilmx + 1
    ELSE
       extra_levels = 1
    END IF

    ! Check to ensure the number of levels set in dynamics is same as physics
    IF (levels .NE. extra_levels) THEN
       PRINT*, 'STOP!! The number of extra levels set in physics is not same in dynamics'
       PRINT*, 'The levels set in the physics is', extra_levels
       PRINT*, 'The levels set in dynamics is', levels
       STOP
    END IF

  END SUBROUTINE initialize_extra_levels


  INTEGER FUNCTION get_extra_levels()
    !----------------------------------------------------------------
    !
    !   Description: Retreive value of extra levels.
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    get_extra_levels = extra_levels

  END FUNCTION get_extra_levels

  SUBROUTINE change_latitude_longitude(latitude, longitude)
    !----------------------------------------------------------------
    !
    !   Description: Changes the value of latitude and longitude in the case of single
    !   -----------   column for max number of grid points (ngridmax)
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    REAL, INTENT(IN) :: latitude, longitude ! in radians

    ! Need to ensure that only a single column, either a check for a flag, or size of ngridmax
    long(1) = longitude
    lati(1) = latitude
    sinlat(1) = sin(lati(1))
    coslat(1) = cos(lati(1))
    sinlon(1) = sin(long(1))
    coslon(1) = cos(long(1))
  END SUBROUTINE change_latitude_longitude

  SUBROUTINE physics_call_single_col(ngrid, nlayer, firstcall, lastcall, rjourvrai, &
       gmtime, ptimestep, pplev, pplay, pphi, pu, pv, pt, extra_temp, &
       pdu, pdv, pdt, pdpsrf)
    !----------------------------------------------------------------
    !
    !   Description: Wrapper routine dynamics will use to call the physics,
    !   -----------  for a single column. It updates the surface and soil temps
    !   -----------  for the column before the call and send back the newly
    !   -----------  update temperatures.
    !
    !   Extra Notes: At the 1st call for each column, the dynamics wont know
    !                what the soil and surface is set to, can send random nums
    !                as phyparam() does cal coldstart() call which sets the tsurf
    !                and tsoil set to 300K.
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------

    INTEGER, INTENT(IN), VALUE :: &
         ngrid,                       & ! Size of the horizontal grid.
         nlayer                         ! Number of vertical layers.
    LOGICAL(KIND=C_BOOL), INTENT(IN), VALUE  :: &
         firstcall,                   & ! True at the first call
         lastcall                       ! True at the last call
    REAL, INTENT(IN), VALUE   ::      &
         rjourvrai,                   & ! Number of days counted from the North. Spring equinox
         gmtime,                      & ! fraction of the day (ranges from 0 to 1)
         ptimestep                      ! timestep (s)
    REAL, INTENT(IN) :: &
         pplev(ngrid,nlayer+1),       & ! Pressure at interfaces between layers (pa)
         pplay(ngrid,nlayer),         & ! Pressure at the middle of the layers (Pa)
         pphi(ngrid,nlayer),          & ! Geopotential at the middle of the layers (m2s-2)
         pu(ngrid,nlayer),            & ! u component of the wind (ms-1)
         pv(ngrid,nlayer),            & ! v component of the wind (ms-1)
         pt(ngrid,nlayer)               ! Temperature (K)
    REAL, INTENT(INOUT) :: &
         extra_temp(ngrid,extra_levels) ! Surface and soil temp (K) of the column
    REAL, INTENT(OUT)   ::            & ! output : physical tendencies
         pdu(ngrid,nlayer),           & ! tendency on zonal wind (m.s-2)
         pdv(ngrid,nlayer),           & ! tendency on meridional wind (m.s-2)
         pdt(ngrid,nlayer),           & ! tendency on temperature (K/s)
         pdpsrf(ngrid)                  ! tendency on surface pressure (Pa/s)

    ! Check if first call or not, if not can't update tsoil/tsurf as they haven't be allocated
    ! They will be initialized in phyparm!
    IF (.not.firstcall) THEN
       ! Update physics soil temp if turned on
       IF (soil_flag) tsoil(:,:) = extra_temp(:, 2:)

       ! Update physics surface temp
       tsurf(:) = extra_temp(ngrid,1)
    END IF

    ! Call phyparam physics
    call phyparam (ngrid, nlayer, firstcall, lastcall, rjourvrai, gmtime, ptimestep, pplev, pplay, &
         pphi, pu, pv, pt, pdu, pdv, pdt, pdpsrf)

    ! Update output extra_temp with physics soil and surface temp
    extra_temp(ngrid,1) = tsurf(ngrid)
    IF (soil_flag) extra_temp(:,2:) = tsoil(:,:)
  END SUBROUTINE physics_call_single_col

END MODULE single_column_mod
