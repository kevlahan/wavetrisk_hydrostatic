MODULE plugins_mod
  IMPLICIT NONE

CONTAINS

  SUBROUTINE mywritefield2(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:,:)
  END SUBROUTINE mywritefield2

  SUBROUTINE mywritefield1(name, longname, unit, var)
    CHARACTER(*), INTENT(IN) :: name, longname, unit
    REAL, INTENT(IN)         :: var(:)
  END SUBROUTINE mywritefield1

END MODULE plugins_mod


PROGRAM simple_physics
  IMPLICIT NONE

  INTEGER, PARAMETER :: nday = 30 ! n-day simulation
  INTEGER, PARAMETER :: step_per_day = 24 ! 1h time step

  INTEGER, PARAMETER :: ngrid=1000, llm=30

  REAL, PARAMETER :: unjours=86400.,    & ! solar day in seconds
       &             radius=6.4e6,      & ! planetary radius
       &             g=9.8,             & ! gravity
       &             cpp=1004.,         & ! Cp
       &             kappa=287/cpp        !kappa=296.945007/cpp ! R/Cp
  REAL, PARAMETER :: timestep=unjours/step_per_day ! physics time step (s)

  REAL :: psurf=1e5, ptop=1e-2, Temp=250. ! initial values of surface pressure, temperature

  REAL :: lon(ngrid), lat(ngrid), & ! in radian
       &  pplev(ngrid,llm+1), &  ! pressure at interfaces
       &  pplay(ngrid, llm),  &  ! pressure at full levels
       &  pphi(ngrid, llm),   &  ! geopotential at full levels
       &  pt(ngrid, llm), &      ! temperature at full levels
       &  pu(ngrid, llm), &      ! zonal wind at full levels
       &  pv(ngrid, llm)         ! meridional wind at full levels

  !$acc data create(pplev, pplay, pphi, pt, pu, pv)

  CALL init_latlon(lat, lon)
  !CALL readin_latlong(lat, lon) ! If want to read latitude and longitude coords from a file
  print *, 'LAT = ' , minval(lat), maxval(lat)
  print *, 'LON = ' , minval(lon), maxval(lon)

  CALL init_physics

  CALL isothermal(kappa*cpp, temp, psurf, ptop, pplev, pplay, pphi, pt)

  pu(:,:)=0.
  pv(:,:)=0.

  CALL init_vitesse(kappa*cpp, temp, psurf, ptop, pplev, pplay, pphi, pt, pu, pv)
  CALL geopot(kappa*cpp, pplev, pplay, pt, pphi)  !! Added by Gabrielle to mimick wavetrisk
  CALL timeloop
  CALL output_function

  !$acc end data

CONTAINS

!!!!!!!!!!!!!!! Added by Gabrielle Ching-Johnson !!!!!!!!!!!!!!!!
  SUBROUTINE output_each_step(day, hour, tot_hour)
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
    write(4,*) 'i, mean_pplay, mean_pphi, mean_pt, mean_pu, mean_pv'

    do ii = 1, llm
       write(4,*) ii, ', ' , sum(pplay(:,ii)) /(max(1,size(pplay(:,ii)))) , ', ' , &
            sum(pphi(:,ii)) /(max(1,size(pphi(:,ii)))),   ', ' , &
            sum(pt(:,ii)) /(max(1,size(pt(:,ii)))), ', ' , &
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
    open(unit=3, file = "../outputs/output_pt")
    open(unit=4, file = "../outputs/output_means_cells")

    write(1,*) 'i,        pu(1;i),        min_pu,        max_pu,       mean_pu'
    write(2,*) 'i,      pv(1;i),      min_pv,      max_pv,       mean_p'
    write(3,*) 'i,      pt(1;i),      min_pt,      max_pt,       mean_pt'
    write(4,*) 'i, mean_pplay, mean_pphi, mean_pt, mean_pu, mean_pv'

    !$acc update host(pu, pv, pt)

    do ii = 1, llm
       write(unit=1,fmt=format) ii, ', ' , pu(1,ii) , ', ' , minval(pu(:,ii)),   ', ' , &
            maxval(pu(:,ii)), ', ' , sum(pu(:,ii)) /(max(1,size(pu(:,ii))))
       write(unit=2,fmt=format) ii, ', ' , pv(1,ii) , ', ' , minval(pv(:,ii)),   ', ' , &
            maxval(pv(:,ii)), ', ' , sum(pv(:,ii)) /(max(1,size(pv(:,ii))))
       write(unit=3,fmt=format) ii, ', ' , pt(1,ii) , ', ' , minval(pt(:,ii)),   ', ' , &
            maxval(pt(:,ii)), ', ' , sum(pt(:,ii))/(max(1,size(pt(:,ii))))
       write(4,*) ii, ', ' , sum(pplay(:,ii)) /(max(1,size(pplay(:,ii)))) , ', ' , &
            sum(pphi(:,ii)) /(max(1,size(pphi(:,ii)))),   ', ' , &
            sum(pt(:,ii)) /(max(1,size(pt(:,ii)))), ', ' , &
            sum(pu(:,ii)) /(max(1,size(pu(:,ii)))), ', ' , &
            sum(pv(:,ii))  /(max(1,size(pv(:,ii))))
    enddo

    close(1)
    close(2)
    close(3)
    close(4)
  END SUBROUTINE output_function


  SUBROUTINE init_vitesse(R, Temp, psurf, ptop, pplev, pplay, pphi, pt, pu, pv)
    REAL, INTENT(IN)  :: R, Temp, psurf, ptop   ! perfect gas constant, temperature, surface and top pressure
    REAL, INTENT(OUT) :: pplev(ngrid,llm+1), &  ! pressure at interfaces
         &               pplay(ngrid, llm), &   ! pressure at full levels
         &               pphi(ngrid, llm), &    ! geopotential phi=gz at full levels
         &               pt(ngrid, llm), &      ! temperature at full levels
         &               pu(ngrid, llm), &      ! zonal wind speed at full levels
         &               pv(ngrid, llm)         ! meridional wind speed at full levels
    REAL :: sigma
    INTEGER :: l
    REAL, PARAMETER :: d=10000 ! thickness of the eckman layer (in meter)
    REAL, PARAMETER :: U0=30 ! geostrophic wind speed (unit : m/s)

    !$acc kernels default(present)
    DO l=1,llm
       pu(:, l) = U0 * (1 - (exp(-pphi(:,l)/ (d*g)) * cos(pphi(:,l)/ (d*g))))
       pv(:, l) = U0 * (1 - (exp(-pphi(:,l)/ (d*g)) * sin(pphi(:,l)/ (d*g))))
    END DO

!!!! Added smoothing to vels to mimick wavetrisk, by Gabrielle !!!!
    DO l = 1,ngrid
       pu(l,:) = cos(lat(l))*pu(l,:)
       pv(l,:) = cos(lat(l))*pv(l,:)
    END DO
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$acc end kernels

  END SUBROUTINE init_vitesse

  SUBROUTINE init_latlon(lat, lon)
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
#ifdef _OPENACC
    USE openacc
    INTEGER :: numdev

    !$acc init

    ! get the number of NVIDIA devices on this node
    numdev = acc_get_num_devices(ACC_DEVICE_NVIDIA)
    PRINT *, 'init_physics : number of available GPU devices', numdev
    IF (numdev < 1) STOP "Error: there are no devices available on this host. ABORTING"
    CALL acc_init(ACC_DEVICE_NVIDIA)
    numdev = acc_get_device_num(ACC_DEVICE_NVIDIA)
    PRINT *, 'init_physics : GPU device id = ', numdev
#endif

    CALL init_comgeomfi(ngrid, llm, lon, lat)
    CALL iniphyparam(timestep, unjours, radius, g, cpp*kappa, cpp)
    writefield1_plugin => mywritefield1
    writefield2_plugin => mywritefield2
  END SUBROUTINE init_physics

  PURE SUBROUTINE isothermal(R, Temp, psurf, ptop, pplev, pplay, pphi, pt)
    ! isothermal profile with sigma-coordinate
    REAL, INTENT(IN) ::  R, Temp, psurf, ptop   ! perfect gas constant, temperature, surface and top pressure
    REAL, INTENT(OUT) :: pplev(ngrid,llm+1), &  ! pressure at interfaces
         &               pplay(ngrid, llm), &   ! pressure at full levels
         &               pphi(ngrid, llm), &    ! geopotential phi=gz at full levels
         &               pt(ngrid, llm)         ! temperature at full levels
    REAL :: sigma
    INTEGER :: l
    ! p = ps * exp(-gz/RT)
    ! => phi = gz = -RT * log(p/ps)
    ! sigma-coordinate : p = ptop + sigma*(ps-ptop), sigma = 1-(l/llm)
    !$acc kernels default(present)
    !$acc loop private(sigma)
    DO l=1,llm
       sigma = 1. - (l-0.5) / REAL(llm)
       pplay(:, l) = ptop + sigma*(psurf-ptop)
       pphi(:, l) = -R*temp*LOG(pplay(:,l)/psurf)
       pt(:, l)    = temp
    END DO
    !$acc loop private(sigma)
    DO l=1,llm+1
       sigma  = 1.-REAL(l-1)/REAL(llm)
       pplev(:, l) = ptop + sigma*(psurf-ptop)
    END DO
    !$acc end kernels
  END SUBROUTINE isothermal

  PURE SUBROUTINE geopot(R, pplev, pplay, pt, pphi)
    REAL, INTENT(IN) ::  R                      ! perfect gas constant
    REAL, INTENT(IN) ::  pplev(ngrid,llm+1), &  ! pressure at interfaces
         &               pplay(ngrid, llm), &   ! pressure at full levels
         &               pt(ngrid, llm)         ! temperature at full levels
    REAL, INTENT(OUT) :: pphi(ngrid, llm)       ! geopotential phi=gz at full levels
    REAL :: rho, gdz, phi(ngrid)
    INTEGER :: l, ig
    !$acc kernels create(phi) default(present)
    phi(:)=0.
    !$acc loop private(rho, gdz)
    DO l=1,llm
       DO ig=1, ngrid
          ! rho = p/RT
          rho = pplay(ig,l)/(R*pt(ig,l))
          gdz = pplev(ig,l)-pplev(ig,l+1)
          gdz  = (pplev(ig,l)-pplev(ig,l+1))/rho ! layer thickness * gravity
          pphi(ig,l) = phi(ig)+.5*gdz            ! geopotential at layer midpoint
          phi(ig) = phi(ig)+gdz                  ! geopotential at upper interface of layer
       END DO
    END DO
    !$acc end kernels
  END SUBROUTINE geopot

  SUBROUTINE timeloop
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_BOOL
    USE phyparam_mod, ONLY : phyparam
    USE callkeys,     ONLY : lverbose
    REAL :: &        ! arrays for tendencies of
         &  pdt(ngrid, llm), &  ! temperature at full levels
         &  pdu(ngrid, llm), &  ! zonal wind at full levels
         &  pdv(ngrid, llm), &  ! meridional wind at full levels
         &  pdpsrf(ngrid)       ! surface pressure

    LOGICAL(KIND=C_BOOL) :: firstcall, lastcall
    REAL :: ptimestep, day_fraction, nth_day
    INTEGER :: istep, iday
    INTEGER, PARAMETER :: desired_write_day=7
    real(8) :: total_time, day

    !day = 86400.0


    !$acc data create(pdt, pdu, pdv, pdpsrf)

    firstcall=.TRUE.
    lastcall=.FALSE.
    ptimestep = unjours/step_per_day

    total_time = 0
    DO iday = 0, nday-1
       PRINT *, 'Temperature at first level', pt(1,1), iday, nday-1
       DO istep = 0, step_per_day-1

!!! Physics call with fraction of the day calculated the same way as in wavetrisk
          !  day_fraction = real(total_time/unjours)
          !  nth_day = floor(day_fraction)
          !  day_fraction = day_fraction - nth_day
          !  CALL phyparam(ngrid,llm,                        &
          !        firstcall,lastcall,                        &
          !        nth_day, day_fraction, timestep, &
          !        pplev,pplay,pphi,                          &
          !        pu,pv,pt,                                  &
          !        pdu,pdv,pdt,pdpsrf)

          CALL phyparam(ngrid,llm,                        &
               firstcall,lastcall,                        &
               1.*iday, timestep*istep/unjours, timestep, &
               pplev,pplay,pphi,                          &
               pu,pv,pt,                                  &
               pdu,pdv,pdt,pdpsrf)

          total_time = total_time + timestep
          firstcall=.FALSE.
          lverbose=.FALSE.

          !$acc kernels default(present)
          pt = pt + pdt*timestep
          pu = pu + pdu*timestep
          pv = pv + pdv*timestep

          !$acc end kernels
          CALL geopot(kappa*cpp, pplev, pplay, pt, pphi)
       END DO
!!!!!!!!!!!!!!!!!!!!!!!!! Added by Gabrielle Ching-Johnoson !!!!!!!!!!!!!!!!!!!!!!!!!!!!
       !Write at the end of every desired_write_day day mean values
       IF ((mod(iday+1,desired_write_day) == 0 .or. iday == 0)) THEN
          CALL output_each_step(iday+1, (istep)*(24/step_per_day), (24*(iday+1)*3600))
       END IF
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    END DO

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !$acc end data

    ! INTEGER, INTENT(IN), VALUE :: &
    !      ngrid,                 & ! Size of the horizontal grid.
    !      nlayer                   ! Number of vertical layers.
    ! LOGICAL(KIND=C_BOOL), INTENT(IN), VALUE  :: &
    !      firstcall,             & ! True at the first call
    !      lastcall                 ! True at the last call
    ! REAL, INTENT(IN), VALUE     ::      &
    !      rjourvrai,             & ! Number of days counted from the North. Spring equinox
    !      gmtime,                & ! fraction of the day (ranges from 0 to 1)
    !      ptimestep                ! timestep (s)
    ! REAL, INTENT(IN) :: &
    !      pplev(ngrid,nlayer+1), & ! Pressure at interfaces between layers (pa)
    !      pplay(ngrid,nlayer),   & ! Pressure at the middle of the layers (Pa)
    !      pphi(ngrid,nlayer),    & ! Geopotential at the middle of the layers (m2s-2)
    !      pu(ngrid,nlayer),      & ! u component of the wind (ms-1)
    !      pv(ngrid,nlayer),      & ! v component of the wind (ms-1)
    !      pt(ngrid,nlayer)         ! Temperature (K)
    ! REAL, INTENT(OUT)   ::      & ! output : physical tendencies
    !      pdu(ngrid,nlayer),     &
    !      pdv(ngrid,nlayer),     &
    !      pdt(ngrid,nlayer),     &
    !      pdpsrf(ngrid)

    ! NB : pdpsrf = 0.

  END SUBROUTINE timeloop

END PROGRAM simple_physics
