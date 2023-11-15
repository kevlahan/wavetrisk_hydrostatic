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

  INTEGER, PARAMETER :: nday = 30 ! 30-day simulation
  INTEGER, PARAMETER :: step_per_day = 24 ! 1h time step

  INTEGER, PARAMETER :: ngrid=100, llm=39

  REAL, PARAMETER :: unjours=86400.,    & ! solar day in seconds
       &             radius=6.4e6,      & ! planetary radius
       &             g=9.8,             & ! gravity
       &             cpp=1004.,         & ! Cp
       &             kappa=296.945007/cpp ! R/Cp
  REAL, PARAMETER :: timestep=unjours/step_per_day ! physics time step (s)

  REAL :: psurf=1e5, ptop=1e4, Temp=250. ! initial values of surface pressure, temperature

  REAL :: lon(ngrid), lat(ngrid), & ! in radian
       &  pplev(ngrid,llm+1), &  ! pressure at interfaces
       &  pplay(ngrid, llm),  &  ! pressure at full levels
       &  pphi(ngrid, llm),   &  ! geopotential at full levels
       &  pt(ngrid, llm), &      ! temperature at full levels
       &  pu(ngrid, llm), &      ! zonal wind at full levels
       &  pv(ngrid, llm)         ! meridional wind at full levels

  !$acc data create(pplev, pplay, pphi, pt, pu, pv)

  CALL init_latlon(lat, lon)
  print *, 'LAT = ' , minval(lat), maxval(lat)
  print *, 'LON = ' , minval(lon), maxval(lon)
  print *, 'LAT = ' , lat
  print *, 'LON = ' , lon

  CALL init_physics

  CALL isothermal(kappa*cpp, temp, psurf, ptop, pplev, pplay, pphi, pt)

  pu(:,:)=0.
  pv(:,:)=0.

  CALL init_vitesse(kappa*cpp, temp, psurf, ptop, pplev, pplay, pphi, pt, pu, pv)

  CALL timeloop
  CALL output_function

  !$acc end data

CONTAINS

  SUBROUTINE output_function
    integer :: ii
    character(60) :: format

    format = "(I2,A2,1PE15.8,A2,1PE15.8,A2,1PE15.8,A2,1PE15.8)"

    open(unit=1, file = "./outputs/output_pu")
    open(unit=2, file = "./outputs/output_pv")
    open(unit=3, file = "./outputs/output_pt")
    open(unit=4, file = "./outputs/output_means_cells")

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
       pu(:, l) = U0 * (1 - exp(-pphi(:,l)/ (d*g))) * cos(pphi(:,l)/ (d*g))
       pv(:, l) = U0 * (1 - exp(-pphi(:,l)/ (d*g))) * sin(pphi(:,l)/ (d*g))
    END DO
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
       lat(ii) = (rand-0.5) * PI

       call random_number(rand)
       lon(ii) = (2*rand -1) * PI

    enddo

  END SUBROUTINE init_latlon

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
          rho = R*pplay(ig,l)/pt(ig,l)
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
    REAL :: ptimestep
    INTEGER :: istep, iday

    !$acc data create(pdt, pdu, pdv, pdpsrf)

    firstcall=.TRUE.
    lastcall=.FALSE.
    ptimestep = unjours/step_per_day

    DO iday = 0, nday-1
       PRINT *, 'Temperature at first level', pt(1,1), iday, nday-1
       DO istep = 0, step_per_day-1

          CALL phyparam(ngrid,llm,                        &
               firstcall,lastcall,                        &
               1.*iday, timestep*istep/unjours, timestep, &
               pplev,pplay,pphi,                          &
               pu,pv,pt,                                  &
               pdu,pdv,pdt,pdpsrf)
          firstcall=.FALSE.
          lverbose=.FALSE.

          !$acc kernels default(present)
          pt = pt + pdt*timestep
          pu = pu + pdu*timestep
          pv = pv + pdv*timestep
          !$acc end kernels
          CALL geopot(kappa*cpp, pplev, pplay, pt, pphi)
       END DO
    END DO

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
