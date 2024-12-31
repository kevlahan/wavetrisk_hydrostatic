!===============================================================================
! CVS $Id: calendar_mod.F90,v 1.1 2000/01/03 22:56:15 kauff Exp $
! CVS $Source: /fs/cgd/csm/models/CVS.REPOS/cpl/cpl5/calendar_mod.F90,v $
! CVS $Name: ccsm2_0_beta47 $
!===============================================================================

MODULE calendar_mod

   !============================================================================
   ! Purpose:
   !   these calendar routines do conversions between...
   !   o the integer number of elapsed days 
   !   o the integer triple (year,month,day)
   !   o the integer coded calendar date (yyyymmdd)
   !
   ! Assumptions:
   !   o there is a year 0
   !   o all years have 365 days (no leap years)
   !   o elapsed days = 0 <=> start of  1 Jan, year 0
   !
   !============================================================================

   integer,private :: dsm(12)   ! elapsed Days on Start of Month
   integer,private :: dpm(12)   ! Days Per Month
   data     dsm  / 0,31,59, 90,120,151, 181,212,243, 273,304,334/
   data     dpm  /31,28,31, 30, 31, 30,  31, 31, 30,  31, 30, 31/
   save


CONTAINS

!===============================================================================

SUBROUTINE eday2date(eday,date)

   implicit none

   integer :: eday,date

   integer :: k,year,month,day

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   compute the calendar date: year/month/day
   ! INPUT:
   !   an integer :: number of elapsed days
   ! OUTPUT:
   !   coded (yyyymmdd) calendar date
   ! NOTE:
   !   this calendar has a year zero (but no day or month zero)
   !----------------------------------------------------------------------------

   year = eday/365       ! calandar year (note: Fortran truncation)
   day  = mod(eday,365)  ! elapsed days within current year
   DO k=1,12
     IF (day .ge. dsm(k)) month=k   ! calendar month
   END DO
   day = day-dsm(month) + 1         ! calendar day
  
   date = year*10000 + month*100 + day  ! coded calendar date

END SUBROUTINE eday2date

!===============================================================================

SUBROUTINE eday2ymd (eday,year,month,day)

   implicit none

   integer :: eday,year,month,day

   integer :: k

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   compute the calendar date: year/month/day
   ! INPUT:
   !   an integer :: number of elapsed days
   ! OUTPUT:
   !   uncoded calendar date, integer :: year, month, & day
   ! NOTE:
   !   this calendar has a year zero (but no day or month zero)
   !----------------------------------------------------------------------------

   year = eday/365       ! calandar year (note: Fortran truncation)
   day  = mod(eday,365)  ! elapsed days within current year
   DO k=1,12
     IF (day .ge. dsm(k)) month=k   ! calendar month
   END DO
   day = day-dsm(month) + 1         ! calendar day

END SUBROUTINE eday2ymd 

!===============================================================================

SUBROUTINE date2ymd (date,year,month,day)

   implicit none

   integer :: date,year,month,day

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   decode the calendar date
   ! INPUT:
   !   calendar date in (integer) yyyymmdd format
   ! OUTPUT:
   !   calendar year,month,day
   !----------------------------------------------------------------------------

   if (.not. valid_date(date)) stop

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

END SUBROUTINE date2ymd 

!===============================================================================

SUBROUTINE date2eday(date,eday)

   implicit none

   integer :: date,eday

   integer :: year,month,day

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   derive elapsed days from the calendar date
   ! INPUT:
   !   calendar date in (integer) yyyymmdd format
   ! OUTPUT:
   !   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
   !----------------------------------------------------------------------------

   if (.not. valid_date(date)) stop 

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

   eday = year*365 + dsm(month) + (day-1)

END SUBROUTINE date2eday

!===============================================================================

SUBROUTINE  ymd2date(year,month,day,date)

   implicit none

   integer :: year,month,day,date

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   encode the calendar date
   ! INPUT:
   !   year, month, & date
   ! OUTPUT:
   !   coded (yyyymmdd) calendar date
   ! NOTE:
   !   this calendar has a year zero (but no day or month zero)
   !----------------------------------------------------------------------------

   if (.not. valid_ymd(year,month,day)) stop 

   date = year*10000 + month*100 + day  ! coded calendar date

END SUBROUTINE  ymd2date

!===============================================================================

SUBROUTINE  ymd2eday(year,month,day,eday)

   implicit none

   integer :: year,month,day,eday

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   derive elapsed days from the calendar date
   ! INPUT:
   !   calendara year, month, & date
   ! OUTPUT:
   !   elapsed days since yy-mm-dd = 00-01-01, with 0 elapsed seconds
   !----------------------------------------------------------------------------

   if (.not. valid_ymd(year,month,day)) stop 

   eday = year*365 + dsm(month) + (day-1)

END SUBROUTINE  ymd2eday

!===============================================================================

FUNCTION valid_date(date) 

   implicit none

   logical :: valid_date
   integer :: date

   integer :: year,month,day

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   determine if a valid year, month & day can be decoded 
   !   from the coded calendar date
   ! INPUT:
   !   calendar date in (integer) yyyymmdd format
   ! RETURNS:
   !   true or false
   !----------------------------------------------------------------------------

   year =int(     date       /10000)
   month=int( mod(date,10000)/  100)
   day  =     mod(date,  100) 

   valid_date = .true.
   if (year  .lt.0) valid_date = .false.
   if (month.lt. 1) valid_date = .false.
   if (month.gt.12) valid_date = .false.
   if (day  .lt. 1) valid_date = .false.
   if (valid_date ) then
     if (day .gt. dpm(month)) valid_date = .false.
   endif

END FUNCTION valid_date

!===============================================================================

FUNCTION valid_ymd(year,month,day)

   implicit none

   integer :: year,month,day
   logical :: valid_ymd

   !----------------------------------------------------------------------------
   ! PURPOSE:
   !   determine if a given year, month & day constitute a valid date
   ! INPUT:
   !   calendar year, month, and day
   ! RETURNS:
   !   true or false
   !----------------------------------------------------------------------------

   valid_ymd = .true.
   if (year  .lt.0) valid_ymd = .false.
   if (month.lt. 1) valid_ymd = .false.
   if (month.gt.12) valid_ymd = .false.
   if (day  .lt. 1) valid_ymd = .false.
   if (valid_ymd) then
      if (day .gt. dpm(month)) valid_ymd = .false.
   endif

END FUNCTION valid_ymd

!===============================================================================

END MODULE calendar_mod
