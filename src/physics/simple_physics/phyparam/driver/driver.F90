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

  real, parameter :: unjours=86400.,    & ! solar day in seconds
       &             radius=6.4e6,      & ! planetary radius
       &             g=9.8,             & ! gravity
       &             cpp=1004.,         & ! cp
       &             kappa=287/cpp        !kappa=296.945007/cpp ! r/cp
  real, parameter :: timestep=unjours/step_per_day ! physics time step (s)

  real :: psurf=1e5, ptop=1e-2, temp=250. ! initial values of surface pressure, temperature

  real :: lon(ngrid), lat(ngrid), & ! in radian
       &  pPlev(ngrid,llm+1), &  ! pressure at interfaces
       &  pPlay(ngrid, llm),  &  ! pressure at full levels
       &  pPhi(ngrid, llm),   &  ! geopotential at full levels
       &  pPhi_surf(ngrid),    &  ! surface geopotential
       &  pT(ngrid, llm), &      ! temperature at full levels
       &  pU(ngrid, llm), &      ! zonal wind at full levels
       &  pV(ngrid, llm)         ! meridional wind at full levels


  CALL init_latlon(lat, lon)
  print *, 'LAT = ' , minval(lat), maxval(lat)
  print *, 'LON = ' , minval(lon), maxval(lon)

  CALL init_physics

  CALL isothermal (kappa*cpp, temp, psurf, ptop, pplev, pplay, pphi, pt)

  pU = 0.0
  pV = 0.0

  CALL init_vitesse (kappa * cpp, Temp, pSurf, pTop, pPlev, pplay, pPhi, pPhi_surf, pT, pU, pV)
  CALL geopot (kappa * cpp, pPlev, pPlay, pT, pPhi, pPhi_surf)  !! Added by Gabrielle to mimick wavetrisk
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


  SUBROUTINE init_vitesse (R, Temp, pSurf, pTop, pPlev, pPlay, pPhi, pPhi_surf, pT, pU, pV)
    REAL, INTENT(IN)  :: R, Temp, psurf, ptop   ! perfect gas constant, temperature, surface and top pressure
    REAL, INTENT(OUT) :: pplev(ngrid,llm+1), &  ! pressure at interfaces
         &               pplay(ngrid, llm), &   ! pressure at full levels
         &               pphi(ngrid, llm), &    ! geopotential phi=gz at full levels
         &               pPhi_surf(ngrid), &    ! surface geopotential
         &               pt(ngrid, llm), &      ! temperature at full levels
         &               pu(ngrid, llm), &      ! zonal wind speed at full levels
         &               pv(ngrid, llm)         ! meridional wind speed at full levels
    REAL :: sigma
    INTEGER :: l
    REAL, PARAMETER :: d=10000 ! thickness of the eckman layer (in meter)
    REAL, PARAMETER :: U0=30 ! geostrophic wind speed (unit : m/s)

    DO l=1,llm
       pU(:, l) = U0 * (1 - (exp(-pPhi(:,l)/ (d*g)) * cos(pphi(:,l)/ (d*g))))
       pV(:, l) = U0 * (1 - (exp(-pPhi(:,l)/ (d*g)) * sin(pphi(:,l)/ (d*g))))
    END DO

!!!! Added smoothing to vels to mimick wavetrisk, by Gabrielle !!!!
    DO l = 1,ngrid
       pu(l,:) = cos(lat(l)) * pU(l,:)
       pv(l,:) = cos(lat(l)) * pV(l,:)
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

  PURE SUBROUTINE geopot (R, pplev, pPlay, pt, pPhi, pPhi_surf)
    REAL, INTENT(IN) ::  R                      ! perfect gas constant
    REAL, INTENT(IN) ::  pplev(ngrid,llm+1), &  ! pressure at interfaces
         &               pplay(ngrid, llm), &   ! pressure at full levels
         &               pt(ngrid, llm)         ! temperature at full levels
    REAL, INTENT(OUT) :: pPhi(ngrid, llm)       ! geopotential phi=gz at full levels
    REAL, INTENT(OUT) :: pPhi_surf(ngrid)       ! geopotential phi=gz at full levels
    REAL :: rho, gdz, phi(ngrid)
    INTEGER :: l, ig

    pPhi_surf = 0d0 ! assume zero surface geopotential
    phi(:) = pPhi_surf

    do l=1,llm
       do ig=1, ngrid
          ! rho = p/rt
          rho = pplay(ig,l)/(r*pt(ig,l))
          gdz = pplev(ig,l)-pplev(ig,l+1)
          gdz  = (pplev(ig,l)-pplev(ig,l+1))/rho ! layer thickness * gravity
          pphi(ig,l) = phi(ig)+.5*gdz            ! geopotential at layer midpoint
          phi(ig) = phi(ig)+gdz                  ! geopotential at upper interface of layer
       end do
    end do
  end subroutine geopot

  subroutine timeloop
    USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_BOOL
    USE phyparam_mod, ONLY : phyparam
    USE callkeys,     ONLY : lverbose
    real, dimension(ngrid,llm) :: pdU    ! zonal velocity tendency
    real, dimension(ngrid,llm) :: pdV    ! meridional velocity tendency
    real, dimension(ngrid,llm) :: pdT    ! temperature tendency
    real, dimension(ngrid)     :: pdPsrf ! surface pressure tendency
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
       PRINT *, 'Temperature at first level', pT(1,1), iday, nday-1
       DO istep = 0, step_per_day-1
          call phyparam (ngrid, llm, mask, firstcall,lastcall, 1.0*iday, timestep * istep/unjours, timestep, &
               pPlev, pPlay, pPhi, pPhi_surf, &
               pU, pV, pT)

          total_time = total_time + timestep
          firstcall = .false.
          lverbose  = .false.

          pT = pT + pdT * Timestep
          pU = pU + pdU * Timestep
          pV = pV + pdV * Timestep

          CALL geopot (kappa*cpp, pplev, pplay, pt, pPhi, pPhi_surf)
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
