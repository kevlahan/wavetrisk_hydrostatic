module plugins_mod
  implicit none
contains
  SUBROUTINE mywritefield2(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:,:)
  END SUBROUTINE mywritefield2

  SUBROUTINE mywritefield1(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:)
  end subroutine mywritefield1
end module plugins_mod

program simple_physics
  implicit none
  integer, parameter :: nday = 30 ! n-day simulation
  integer, parameter :: step_per_day = 24 ! 1h time step

  integer, parameter :: ngrid=1000, llm=30, mask = 8

  real, parameter :: &
       unjours = 86400.0,   & ! solar day in seconds
       radius  = 6.4e6,     & ! planetary radius
       g       = 9.8,       & ! gravity
       cpp     = 1004.0,    & ! cp
       kappa   = 287/Cpp      ! kappa =2 96.945007/cpp ! r/cp
  real, parameter :: timestep=unjours/step_per_day ! physics time step (s)

  real :: psurf=1e5, ptop=1e-2, temp=250. ! initial values of surface pressure, temperature

  real :: &
       lon(ngrid), lat(ngrid), & ! in radian
       pPint(ngrid,llm+1), &                           ! pressure at interfaces
       pPlay(ngrid,llm),   &                           ! pressure at full levels
       pPhi(ngrid,llm),    &                           ! geopotential at full levels
       pPhi_surf(ngrid),   &                           ! surface geopotential
       pdU2dZ(ngrid,llm+1), pUmag(ngrid,llm), &
       pTheta(ngrid, llm), &                           ! potential temperature
       pU(ngrid, llm), pV(ngrid, llm), pW(ngrid, llm)  ! velocities


  CALL init_latlon(lat, lon)
  print *, 'LAT = ' , minval(lat), maxval(lat)
  print *, 'LON = ' , minval(lon), maxval(lon)

  CALL init_physics

  pU = 0.0
  pV = 0.0
  pW = 0.0

  CALL init_vitesse (kappa * cpp, Temp, pSurf, pTop, pPint, pPlay, pPhi, pPhi_surf, pTheta, pU, pV, pW)
  CALL timeloop
  CALL output_function
contains
!!!!!!!!!!!!!!! Added by Gabrielle Ching-Johnson !!!!!!!!!!!!!!!!
  SUBROUTINE output_each_step (day, hour, tot_hour)
    !-----------------------------------------------------------------------
    !
    !   Description: Outputs mean values at desired step when called.
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------
    integer, INTENT(IN) :: day, hour, tot_hour
    integer :: i, ii
    character(40) :: format
    character(4) :: s_time

    format = "(I5,A1,I2,A1,I6,A1,1P 79E16.8)"

    write (s_time, '(i4.4)') day
    open(unit=4, file = '../outputs/output_means_cells.'//s_time)
    write(4,*) day, hour, tot_hour
    write(4,*) 'i, mean_pplay, mean_pphi, mean_pTheta, mean_pu, mean_pv'

    do ii = 1, llm
       write(4,*) ii, ', ' , sum(pplay(:,ii)) /(max(1,size(pplay(:,ii)))) , ', ' , &
            sum(pphi(:,ii)) /(max(1,size(pphi(:,ii)))),   ', ' , &
            sum(pTheta(:,ii)) /(max(1,size(pTheta(:,ii)))), ', ' , &
            sum(pu(:,ii)) /(max(1,size(pu(:,ii)))), ', ' , &
            sum(pv(:,ii))  /(max(1,size(pv(:,ii))))
    enddo

    close(4)

  END SUBROUTINE output_each_step

  SUBROUTINE output_lat_long
    !-----------------------------------------------------------------------
    !
    !   Description: Outputs all latitude and longitude coordinates to
    !                 file called output_lat_long.
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------
    integer :: col
    character(40) :: format

    format = "(F15.6,A2,F15.6)"

    open(unit=1, file = "../outputs/output_lat_long")
    write(1,*) 'Col, Lat,   Long'

    do col= 1, ngrid
       write(1,*) lat(col), ', ', lon(col)
    enddo

    close(1)
  END SUBROUTINE output_lat_long


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  SUBROUTINE output_function
    !-----------------------------------------------------------------------
    !
    !   Description: Outputs mean values. Called at end of simulation.
    !
    !-----------------------------------------------------------------------
    integer :: ii, i
    character(60) :: format

    format = "(I4,A2,1PE15.8,A2,1PE15.8,A2,1PE15.8,A2,1PE15.8)"

    open(unit=1, file = "../outputs/output_pu")
    open(unit=2, file = "../outputs/output_pv")
    open(unit=3, file = "../outputs/output_pTheta")
    open(unit=4, file = "../outputs/output_means_cells")

    write(1,*) 'i,        pu(1;i),        min_pu,        max_pu,       mean_pu'
    write(2,*) 'i,      pv(1;i),      min_pv,      max_pv,       mean_p'
    write(3,*) 'i,      pTheta(1;i),      min_pTheta,      max_pTheta,       mean_pTheta'
    write(4,*) 'i, mean_pplay, mean_pphi, mean_pTheta, mean_pu, mean_pv'

    do ii = 1, llm
       write(unit=1,fmt=format) ii, ', ' , pu(1,ii) , ', ' , minval(pu(:,ii)),   ', ' , &
            maxval(pu(:,ii)), ', ' , sum(pu(:,ii)) /(max(1,size(pu(:,ii))))
       write(unit=2,fmt=format) ii, ', ' , pv(1,ii) , ', ' , minval(pv(:,ii)),   ', ' , &
            maxval(pv(:,ii)), ', ' , sum(pv(:,ii)) /(max(1,size(pv(:,ii))))
       write(unit=3,fmt=format) ii, ', ' , pTheta(1,ii) , ', ' , minval(pTheta(:,ii)),   ', ' , &
            maxval(pTheta(:,ii)), ', ' , sum(pTheta(:,ii))/(max(1,size(pTheta(:,ii))))
       write(4,*) ii, ', ' , sum(pplay(:,ii)) /(max(1,size(pplay(:,ii)))) , ', ' , &
            sum(pphi(:,ii)) /(max(1,size(pphi(:,ii)))),   ', ' , &
            sum(pTheta(:,ii)) /(max(1,size(pTheta(:,ii)))), ', ' , &
            sum(pu(:,ii)) /(max(1,size(pu(:,ii)))), ', ' , &
            sum(pv(:,ii))  /(max(1,size(pv(:,ii))))
    enddo

    close(1)
    close(2)
    close(3)
    close(4)
  END SUBROUTINE output_function


  SUBROUTINE init_vitesse (R, Temp, pSurf, pTop, pPint, pPlay, pPhi, pPhi_surf, pTheta, pU, pV, pW)
    REAL, INTENT(IN)  :: R, Temp, psurf, ptop   ! perfect gas constant, temperature, surface and top pressure
    REAL, INTENT(OUT) ::     &
         pPint(ngrid,llm+1), &    ! pressure at interfaces
         pPlay(ngrid,llm),  &    ! pressure at full levels
         pPhi(ngrid,llm),   &    ! geopotential phi=gz at full levels
         pPhi_surf(ngrid),   &    ! surface geopotential
         pTheta(ngrid, llm), &    ! potential temperature at full levels
         pU(ngrid, llm),  pV(ngrid, llm), pW(ngrid, llm) ! velocities
    REAL :: sigma
    INTEGER :: l
    REAL, PARAMETER :: d=10000 ! thickness of the eckman layer (in meter)
    REAL, PARAMETER :: U0=30 ! geostrophic wind speed (unit : m/s)

    pPhi = 0.0
    DO l=1,llm
       pU(:,l) = 0.0
       pV(:,l) = 0.0
       pW(:,l) = 0.0
    END DO
  END SUBROUTINE init_vitesse

  SUBROUTINE init_latlon (lat, lon)
    REAL, INTENT(OUT)  :: lat(ngrid), lon(ngrid)
    integer :: ii
    REAL  :: rand
    REAL, PARAMETER :: PI = 2.*asin(1.)
    INTEGER :: seed_size, seed(100)
    ! initialize seed to obtain a reproducible sequence of random numbers
    CALL RANDOM_SEED( SIZE = seed_size)
    PRINT *, 'Size of random seed array : ', seed_size
    DO ii=1, seed_size
       seed(ii) = ii
    END DO
    CALL RANDOM_SEED( PUT = seed )

    do ii=1, ngrid
       call random_number(rand)
       lat(ii) = asin(-1 + 2*rand)

       call random_number(rand)
       lon(ii) = 2*PI*rand

    enddo
!!!!!!!!!!!!!!!!!!!!! Added by Gabrielle !!!!!!!!!!!!!!!!!!!!!!!!!!!
    !CALL output_lat_long
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  END SUBROUTINE init_latlon

  SUBROUTINE readin_latlong(lat, lon)
    !-----------------------------------------------------------------------
    !
    !   Description: Read in desired latitude and longitude coordinates
    !                 from file coordinate_wavetrisk. (Can change file name)
    !
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !----------------------------------------------------------------------
    REAL, INTENT(OUT)  :: lat(ngrid), lon(ngrid)
    integer :: i

    !open file
    open(unit=5, file='coordinate_wavetrisk')

    !read from file
    do i = 1, ngrid
       read(5,*) lat(i), lon(i)
    enddo
    close(5)
  END SUBROUTINE readin_latlong

  SUBROUTINE init_physics
    USE comgeomfi
    USE iniphyparam_mod
    USE writefield_mod, ONLY : writefield1_plugin, writefield2_plugin
    USE plugins_mod
    INTEGER :: numdev

    CALL init_comgeomfi(ngrid, llm, lon, lat)
    CALL iniphyparam(timestep, unjours, radius, g, cpp*kappa, cpp)
    writefield1_plugin => mywritefield1
    writefield2_plugin => mywritefield2
  END SUBROUTINE init_physics

  subroutine timeloop
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_BOOL
    USE phyparam_mod, ONLY : phyparam
    USE callkeys,     ONLY : lverbose
    real, dimension(ngrid)     :: pPhi_surf! surface geopotential

    integer                    :: istep, iday
    integer, parameter         :: desired_write_day = 7
    real                       :: ptimestep, day_fraction, nth_day
    real(8)                    :: total_time, day
    logical(kind=c_bool)       :: firstcall, lastcall

    firstcall = .true.
    lastcall  = .false.
    ptimestep = unJours / step_per_day

    total_time = 0
    DO iday = 0, nday-1
       PRINT *, 'Potential temperature at first level', pTheta(1,1), iday, nday-1
       DO istep = 0, step_per_day-1
          call phyparam (ngrid, llm, mask, firstcall,lastcall, 1.0*iday, timestep * istep/unjours, timestep, &
               pPint, pPlay, pPhi, pPhi_surf, pUmag, &
               pU, pV, pW, pTheta)

          total_time = total_time + timestep
          firstcall = .false.
          lverbose  = .false.
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!! Added by Gabrielle Ching-Johnoson !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       ! Write at the end of every desired_write_day day mean values
       if ((mod(iday+1,desired_write_day) == 0 .or. iday == 0)) then
          call output_each_step(iday+1, (istep)*(24/step_per_day), (24*(iday+1)*3600))
       end if
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    end do
  end subroutine timeloop
end program simple_physics
