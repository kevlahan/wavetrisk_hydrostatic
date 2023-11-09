MODULE phyparam_mod
#include "use_logging.h"
  USE callkeys
  USE comgeomfi
  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL, PARAMETER :: ref_temp=285., capcal_nosoil=1e5, tsoil_init=300.

  INTEGER :: icount
  REAL    :: zday_last
  LOGICAL :: firstcall_alloc=.true.
  !$OMP THREADPRIVATE( icount,zday_last)

  PUBLIC :: phyparam, alloc, precompute, zday_last, icount

CONTAINS

  SUBROUTINE phyparam(ngrid,nlayer,             &
       &            firstcall,lastcall,         &
       &            rjourvrai,gmtime,ptimestep, &
       &            pplev,pplay,pphi,           &
       &            pu,pv,pt,                   &
       &            pdu,pdv,pdt,pdpsrf) BIND(C, name='phyparam_phyparam')
    USE, INTRINSIC :: ISO_C_BINDING
    USE phys_const, ONLY : g, rcp, r, unjours
    USE soil_mod,   ONLY : soil_forward, soil_backward
    USE soil_mod,   ONLY : z0, inertie, emissiv, albedo  ! precomputed
    USE soil_mod,   ONLY : tsurf, tsoil ! state variables
    USE turbulence, ONLY : vdif
    USE convection, ONLY : convadj
    USE radiative_mod,  ONLY : radiative_tendencies
    USE writefield_mod, ONLY : writefield
    USE accelerator,    ONLY : nvtxStartRange, nvtxEndRange
    USE profiling,      ONLY : profile_enter, profile_exit, id_phyparam
    !=======================================================================
    !   Top routine of the physical parametrisations of the LMD
    !   20 parameters GCM for planetary atmospheres.
    !   It includes:
    !   radiative transfer (long and shortwave) for CO2 and dust.
    !   vertical turbulent mixing
    !   convective adjsutment
    !   heat diffusion in the soil
    !
    !   author: Frederic Hourdin 15 / 10 /93
    !=======================================================================

    INTEGER, INTENT(IN), VALUE :: &
         ngrid,                 & ! Size of the horizontal grid.
         nlayer                   ! Number of vertical layers.
    LOGICAL(KIND=C_BOOL), INTENT(IN), VALUE  :: &
         firstcall,             & ! True at the first call
         lastcall                 ! True at the last call
    REAL, INTENT(IN), VALUE     ::      &
         rjourvrai,             & ! Number of days counted from the North. Spring equinox
         gmtime,                & ! fraction of the day (ranges from 0 to 1)
         ptimestep                ! timestep (s)
    REAL, INTENT(IN) :: &
         pplev(ngrid,nlayer+1), & ! Pressure at interfaces between layers (pa)
         pplay(ngrid,nlayer),   & ! Pressure at the middle of the layers (Pa)
         pphi(ngrid,nlayer),    & ! Geopotential at the middle of the layers (m2s-2)
         pu(ngrid,nlayer),      & ! u component of the wind (ms-1)
         pv(ngrid,nlayer),      & ! v component of the wind (ms-1)
         pt(ngrid,nlayer)         ! Temperature (K)
    REAL, INTENT(OUT)   ::      & ! output : physical tendencies
         pdu(ngrid,nlayer),     & ! tendency on zonal wind (m.s-2)
         pdv(ngrid,nlayer),     & ! tendency on meridional wind (m.s-2)
         pdt(ngrid,nlayer),     & ! tendency on temperature (K/s)
         pdpsrf(ngrid)            ! tendency on surface pressure (Pa/s)

    !    Local variables :
    REAL :: zh(ngrid,nlayer),      & ! potential temperature
         &  zpopsk(ngrid,nlayer),  & ! Exner function
         &  zzlev(ngrid,nlayer+1), & ! altitude of interfaces
         &  zzlay(ngrid,nlayer),   & ! altitude of full levels
         &  fluxrad(ngrid),        & ! radiative flux at surface
         &  zc(ngrid, nsoilmx-1),    & ! LU coefficients for soil implicit solve
         &  zd(ngrid, nsoilmx-1),    &
         &  fluxgrd(ngrid),        & ! heat flux from deep soil
         &  capcal(ngrid),         & ! effective heat capacity of soil
         &  zdufr(ngrid,nlayer),   & ! partial tendencies for zonal velocity,
         &  zdvfr(ngrid,nlayer),   & !   meridional velocity,
         &  zdhfr(ngrid,nlayer),   & !   potential temperature,
         &  zdtsrfr(ngrid),        & !   surface temperature
         &  zdtsrf(ngrid),         & ! total tendency of surface temperature
         &  zflubid(ngrid),        & ! radiative + deep soil fluxes
         &  zpmer(ngrid)             ! sea-level pressure
    REAL zdum1(ngrid,nlayer)
    REAL zdum2(ngrid,nlayer)
    REAL zdum3(ngrid,nlayer)

    INTEGER :: j,l,ig,igout
    LOGICAL :: lwrite
    REAL    :: zday, zdtime, z1, z2

    call nvtxStartRange("Physics")

   ! WRITELOG(*,*) 'latitude0', ngrid, lati(1:2), lati(ngrid-1:ngrid)
   ! WRITELOG(*,*) 'nlayer',nlayer
   ! LOG_DBG('phyparam')

    IF (ngrid.NE.ngridmax) THEN
       PRINT*,'STOP in inifis'
       PRINT*,'Probleme de dimensions :'
       PRINT*,'ngrid     = ',ngrid
       PRINT*,'ngridmax   = ',ngridmax
       STOP
    ENDIF

    igout=ngrid/2+1
    zday=rjourvrai+gmtime

    !-----------------------------------------------------------------------
    !    0. Allocate and initialize at first call
    !    --------------------
        
    IF(firstcall) THEN
       !         zday_last=rjourvrai
       zday_last=zday-ptimestep/unjours
       ! Alloc only if very first call of entire sim, especially in single col cases
       IF (firstcall_alloc) THEN
          CALL alloc(ngrid, nlayer)
       END IF
       CALL precompute
       CALL coldstart(ngrid)
       firstcall_alloc = .FALSE. !turn flag false, such tht allocation occurs once
    ENDIF

    !$acc data copyin(pplev(:,:),pplay(:,:),pphi(:,:),pu(:,:),pv(:,:),pt(:,:)) &
    !$acc copyout(pdu(:,:),pdv(:,:),pdt(:,:),pdpsrf(:)) &
    !$acc present(z0, inertie, emissiv, albedo, tsurf, tsoil) &
    !$acc create (zh,zpopsk,zzlev,zzlay,fluxrad,zc,zd,fluxgrd, &
    !$acc              capcal,zdufr,zdvfr,zdhfr,zdtsrfr,zdtsrf,zflubid,zpmer, &
    !$acc              zdum1,zdum2,zdum3)

    CALL profile_enter(id_phyparam)

    !-----------------------------------------------------------------------
    !    1. Initialisations :
    !    --------------------

    icount=icount+1

    IF(abs(zday-zday_last-period_sort)<=ptimestep/unjours/10.) THEN
       WRITELOG(*,*) 'zday, zday_last SORTIE ', zday, zday_last
       LOG_INFO('phyparam')
       zday_last=zday
       lwrite = .TRUE.
    ELSE
       lwrite = .FALSE.
    END IF

    !$acc kernels default(none)
    pdv(:,:)  = 0.
    pdu(:,:)  = 0.
    pdt(:,:)  = 0.
    pdpsrf(:) = 0.
    zflubid(:)= 0.
    zdtsrf(:) = 0.

    !-----------------------------------------------------------------------
    !   calcul du geopotentiel aux niveaux intercouches
    !   ponderation des altitudes au niveau des couches en dp/p

    DO l=1,nlayer
       DO ig=1,ngrid
          zzlay(ig,l)=pphi(ig,l)/g
       ENDDO
    ENDDO
    DO ig=1,ngrid
       zzlev(ig,1)=0.
    ENDDO
    DO l=2,nlayer
       DO ig=1,ngrid
          z1=(pplay(ig,l-1)+pplev(ig,l))/(pplay(ig,l-1)-pplev(ig,l))
          z2=(pplev(ig,l)+pplay(ig,l))/(pplev(ig,l)-pplay(ig,l))
          zzlev(ig,l)=(z1*zzlay(ig,l-1)+z2*zzlay(ig,l))/(z1+z2)
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   Transformation de la temperature en temperature potentielle
    DO l=1,nlayer
       DO ig=1,ngrid
          zpopsk(ig,l)=(pplay(ig,l)/pplev(ig,1))**rcp ! surface pressure is used as reference pressure
          zh(ig,l)=pt(ig,l)/zpopsk(ig,l)
       ENDDO
    ENDDO
    !$acc end kernels

    !-------------------------------------------------------------
    !  soil temperatures : 1st half of implicit time integration
    !  forward sweep from deep ground to surface
    !  yields LU coefficients zc,zd and capcal, fluxgrd
    !   ----------------------------------------------------------
    IF (callsoil) THEN

       CALL soil_forward(ngrid,nsoilmx, ptimestep, inertie, tsurf, tsoil, &
            &            zc, zd, capcal, fluxgrd)

       !       CALL soil_new(ngrid,nsoilmx,ptimestep,inertie, &
       !            tsurf, tsoil, capcal,fluxgrd)
       !       CALL soil(ngrid,nsoilmx,.false.,inertie, &
       !            &          ptimestep,tsurf,tsoil,capcal,fluxgrd)
    ELSE
       !$acc kernels default(none)
       capcal(:)  = capcal_nosoil
       fluxgrd(:) = 0.
       !$acc end kernels
    END IF

    IF(lverbose) THEN
       !$acc update host (capcal,fluxgrd,tsurf)
       WRITELOG(*,*) 'Surface Heat capacity, conduction Flux, Ts'
       WRITELOG(*,*) capcal(igout), fluxgrd(igout), tsurf(igout)
       LOG_DBG('phyparam')
    END IF

    !-----------------------------------------------------------------------
    !    2. Compute radiative tendencies :
    !    ---------------------------------------

    IF(callrad) CALL radiative_tendencies(lwrite, ngrid, igout, nlayer, &
         gmtime, ptimestep*float(iradia), zday, pplev, pplay, pt, &
         pdt, fluxrad)

    !-----------------------------------------------------------------------
    !    3. Vertical diffusion (turbulent mixing):
    ! Kz is computed then vertical diffusion is integrated in time implicitly
    ! using a linear relationship between surface heat flux and air temperature
    ! in lowest level (Robin-type BC)
    !    -------------------------------------------------------------------
    !
    IF(calldifv) THEN

       !$acc kernels default(none)
       DO ig=1,ngrid
          zflubid(ig)=fluxrad(ig)+fluxgrd(ig)
       ENDDO

       zdum1(:,:)=0.
       zdum2(:,:)=0.

       do l=1,nlayer
          do ig=1,ngrid
             zdum3(ig,l)=pdt(ig,l)/zpopsk(ig,l)
          enddo
       enddo

       !$acc end kernels
       CALL vdif(ngrid,nlayer,zday,        &
            &        ptimestep,capcal,z0,       &
            &        pplay,pplev,zzlay,zzlev,   &
            &        pu,pv,zh,tsurf,emissiv,    &
            &        zdum1,zdum2,zdum3,zflubid, &
            &        zdufr,zdvfr,zdhfr,zdtsrfr, &
            &        lverbose)

       !$acc kernels default(none)
       DO l=1,nlayer
          DO ig=1,ngrid
             pdv(ig,l)=pdv(ig,l)+zdvfr(ig,l)
             pdu(ig,l)=pdu(ig,l)+zdufr(ig,l)
             pdt(ig,l)=pdt(ig,l)+zdhfr(ig,l)*zpopsk(ig,l)
          ENDDO
       ENDDO

       DO ig=1,ngrid
          zdtsrf(ig)=zdtsrf(ig)+zdtsrfr(ig)
       ENDDO
       !$acc end kernels
    ELSE
       !$acc kernels default(none)
       DO ig=1,ngrid
          zdtsrf(ig)=zdtsrf(ig)+ &
               &      (fluxrad(ig)+fluxgrd(ig))/capcal(ig)
       ENDDO
       !$acc end kernels
    ENDIF

    !-------------------------------------------------------------
    !   soil temperatures : 2nd half of implicit time integration
    !   using updated tsurf as input
    !   ----------------------------------------------------------

    !$acc kernels default(none)
    DO ig=1,ngrid
       tsurf(ig)=tsurf(ig)+ptimestep*zdtsrf(ig)
    ENDDO
    !$acc end kernels

    WRITE(55,'(2e15.5)') zday,tsurf(ngrid/2+1)

    IF (callsoil) THEN
       CALL soil_backward(ngrid,nsoilmx, zc,zd, tsurf,tsoil)
       IF(lverbose) THEN
          !$acc update host (tsurf,zdtsrf)
          WRITELOG(*,*) 'Surface Ts, dTs, dt'
          WRITELOG(*,*) tsurf(igout), zdtsrf(igout), ptimestep
          LOG_DBG('phyparam')
       ENDIF
    END IF

    !
    !-----------------------------------------------------------------------
    !   4. Dry convective adjustment:
    !   -----------------------------

    IF(calladj) THEN

       !$acc kernels default(none)
       DO l=1,nlayer
          DO ig=1,ngrid
             zdum1(ig,l)=pdt(ig,l)/zpopsk(ig,l)
          ENDDO
       ENDDO

       zdufr(:,:)=0.
       zdvfr(:,:)=0.
       zdhfr(:,:)=0.
       !$acc end kernels

       CALL convadj(ngrid,nlayer,ptimestep, &
            &                pplay,pplev,zpopsk, &
            &                pu,pv,zh,           &
            &                pdu,pdv,zdum1,      &
            &                zdufr,zdvfr,zdhfr)

       !$acc kernels default(none)
       DO l=1,nlayer
          DO ig=1,ngrid
             pdu(ig,l)=pdu(ig,l)+zdufr(ig,l)
             pdv(ig,l)=pdv(ig,l)+zdvfr(ig,l)
             pdt(ig,l)=pdt(ig,l)+zdhfr(ig,l)*zpopsk(ig,l)
          ENDDO
       ENDDO
       !$acc end kernels

    ENDIF

    !-----------------------------------------------------------------------
    !   sorties:
    !   --------

    WRITELOG(*,*) 'zday, zday_last ',zday,zday_last,icount
    LOG_DBG('phyparam')

    IF(lwrite) THEN

       do ig=1,ngridmax
          zpmer(ig)=pplev(ig,1)*exp(pphi(ig,1)/(r*ref_temp))
       enddo
       !$acc update host (zh)
       call writefield('u','Vent zonal moy','m/s',pu)
       call writefield('v','Vent meridien moy','m/s',pv)
       call writefield('temp','Temperature','K',pt)
       call writefield('theta','Potential temperature','K',zh)
       call writefield('geop','Geopotential','m2/s2',pphi)
       call writefield('plev','plev','Pa',pplev(:,1:nlayer))

       !$acc update host (pdu,pdv,pdt)
       call writefield('du','du',' ',pdu)
       call writefield('dv','dv',' ',pdv)
       call writefield('dt','dt',' ',pdt)

       !$acc update host (tsurf,coslon,sinlon,coslat,sinlat,albedo)
       call writefield('ts','Surface temper','K',tsurf)
       call writefield('coslon','coslon',' ',coslon)
       call writefield('sinlon','sinlon',' ',sinlon)
       call writefield('coslat','coslat',' ',coslat)
       call writefield('sinlat','sinlat',' ',sinlat)
       call writefield('alb','alb',' ',albedo)
       call writefield('ps','Surface pressure','Pa',pplev(:,1))
       call writefield('slp','Sea level pressure','Pa',zpmer)
    END IF

    CALL profile_exit(id_phyparam)

    !$acc end data

    call nvtxEndRange

  END SUBROUTINE phyparam

  SUBROUTINE alloc(ngrid, nlayer) BIND(C, name='phyparam_alloc')
    !$cython header void phyparam_alloc(int, int);
    !$cython wrapper def alloc(ngrid, nlayer) : phy.phyparam_alloc(ngrid, nlayer)
    USE astronomy, ONLY : iniorbit
    USE soil_mod, ONLY: tsurf, tsoil, z0, inertie, rnatur, albedo, emissiv
    INTEGER, INTENT(IN), VALUE :: ngrid, nlayer
    ! allocate precomputed arrays
    ALLOCATE(rnatur(ngrid), albedo(ngrid), emissiv(ngrid))

    ALLOCATE(z0(ngrid),inertie(ngrid))

    ! allocate arrays for internal state
    ALLOCATE(tsurf(ngrid))
    ALLOCATE(tsoil(ngrid,nsoilmx))

    CALL iniorbit
  END SUBROUTINE alloc

  SUBROUTINE precompute() BIND(C, name='phyparam_precompute')
    !$cython header void phyparam_precompute();
    !$cython wrapper def precompute() : phy.phyparam_precompute()
    ! precompute time-independent arrays
    USE soil_mod, ONLY: rnatur, inertie, z0, emissiv, albedo, &
         I_mer,I_ter,Cd_mer,Cd_ter, &
         alb_mer,alb_ter,emi_mer,emi_ter
    !$acc kernels default(none)
    rnatur(:)  = 1.
    inertie(:) = (1.-rnatur(:))*I_mer+rnatur(:)*I_ter
    z0(:)      = (1.-rnatur(:))*Cd_mer+rnatur(:)*Cd_ter
    emissiv(:) = (1.-rnatur(:))*emi_mer+rnatur(:)*emi_ter
    albedo(:)  = (1.-rnatur(:))*alb_mer+rnatur(:)*alb_ter
    !$acc end kernels
  END SUBROUTINE precompute

  SUBROUTINE coldstart(ngrid) BIND(C, name='phyparam_coldstart')
    !$cython header void phyparam_coldstart(int);
    !$cython wrapper def coldstart (ngrid): phy.phyparam_coldstart(ngrid)
    ! create internal state to start a run without a restart file
    USE soil_mod, ONLY : tsurf, tsoil
    INTEGER, INTENT(IN), VALUE :: ngrid
    !$acc kernels default(none)
    tsurf(:)   = tsoil_init
    tsoil(:,:) = tsoil_init
    !$acc end kernels
    icount=0
    IF(.NOT. callsoil .AND. firstcall_alloc .AND. lverbose) THEN !write only on first call (lverbose set by interface)
       WRITELOG(*,*) 'WARNING!!! Thermal conduction in the soil turned off'
       LOG_WARN('coldstart')
    ENDIF
  END SUBROUTINE coldstart

END MODULE phyparam_mod
