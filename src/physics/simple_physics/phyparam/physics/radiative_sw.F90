MODULE radiative_sw

#include "use_logging.h"

  IMPLICIT NONE
  SAVE

  PRIVATE

  PUBLIC :: sw

CONTAINS

  PURE SUBROUTINE monGATHER(ngrid, n,index, a,b)
    INTEGER, INTENT(IN) :: ngrid, n, index(n)
    REAL, INTENT(IN)    :: a(ngrid)
    REAL, INTENT(OUT)   :: b(n)
    INTEGER :: i
    IF(n<ngrid) THEN
       DO i=1,n
          b(i)=a(index(i))
       END DO
    ELSE
       b(:)=a(:)
    END IF
  END SUBROUTINE monGATHER

  PURE subroutine monscatter(ngrid, n,index, b,a)
    INTEGER, INTENT(IN) :: ngrid, n,index(n)
    REAL, INTENT(IN)    :: b(n)
    REAL, INTENT(OUT)   :: a(ngrid)
    INTEGER :: i
    IF(n<ngrid) THEN
       a(:)=0.
       DO i=1,n
          a(index(i))=b(i)
       END DO
    ELSE
       a(:)=b(:)
    END IF
  end subroutine monscatter

  SUBROUTINE sw(ngrid,nlayer,ldiurn, coefvis,albedo, &
       &        plevel,ps_rad,pmu,pfract,psolarf0, &
       &        fsrfvis,dtsw, lverbose, lwrite)
    USE phys_const, ONLY : cpp, g
    USE writefield_mod, ONLY : writefield

    !=======================================================================
    !
    !   Rayonnement solaire en atmosphere non diffusante avec un
    !   coefficient d absorption gris.
    !
    !=======================================================================

    INTEGER, INTENT(IN) :: ngrid, nlayer
    LOGICAL, INTENT(IN) :: ldiurn, lverbose, lwrite
    REAL, INTENT(IN)    :: &
         psolarf0,            & ! solar constant
         ps_rad, coefvis,     & ! coefvis = attenuation at p=ps_rad
         albedo(ngrid),       & ! albedo
         pmu(ngrid),          & ! cosine zenithal angle
         pfract(ngrid),       & ! day fraction
         plevel(ngrid,nlayer+1) ! pressure at interfaces
    REAL, INTENT(OUT) :: &
         fsrfvis(ngrid),    & ! net surface flux
         dtsw(ngrid,nlayer)   ! temperature tendency

    REAL :: buf1(ngrid), buf2(ngrid, nlayer+1) ! buffers for output
    ! fluxes are non-zero only on those points where the sun shines (mu0>0)
    ! We compute only on those ncount points and gather them to vectorize
    INTEGER :: ncount, index(ngrid)
    ! In the work arrays below, ngrid should be ncount but ncount is not known yet
    REAL :: zalb(ngrid),                & ! albedo
         &  zmu(ngrid),                 & ! cosine zenithal angle
         &  zfract(ngrid),              & ! day fraction
         &  flux_in(ngrid),             & ! incoming solar flux
         &  flux_down(ngrid, nlayer+1), & ! downward flux
         &  flux_up(ngrid, nlayer+1),   & ! upward flux
         &  zplev(ngrid,nlayer+1),      & ! pressure at interfaces
         &  zflux(ngrid),               & ! net surface flux
         &  zdtsw(ngrid,nlayer),        & ! temperature tendency
         &  zu(ngrid,nlayer+1)
    INTEGER :: ig,l,igout
    REAL :: tau0

    !$acc data copyin(albedo, pmu, pfract, plevel)      &
    !$acc      copyout(fsrfvis, dtsw)                   &
    !$acc      create(flux_in, flux_down, flux_up, zu)

#ifdef _OPENACC

    ! With OpenACC we skip the gather/scatter phase
    ! therefore we read/write directly from/to arrays
    !      zfract, zmu, zalb, dtsw
    ! instead of buffers
    !      zfract, zmu, zalb, zdtsw
#define NCOUNT ngrid

#define ZFRACT pfract
#define ZMU    pmu
#define ZALB   albedo
#define ZPLEV  plevel

#define ZDTSW  dtsw
#define ZFLUX  fsrfvis

#else
    ! on the CPU we compute only on sunny points (pfract(ig)>0)
    ! we list their indices and pack the input arrays into buffers (mongather)
    ! on exit we unpack the output flux and temperature tendencies from buffers to full-size arrays (monscatter)

    IF (ldiurn) THEN
       ncount=0
       DO ig=1,ngrid
          index(ig)=0
       ENDDO
       DO ig=1,ngrid
          IF(pfract(ig).GT.1.e-6) THEN
             ncount=ncount+1
             index(ncount)=ig
          ENDIF
       ENDDO
    ELSE
       ncount=ngrid
    ENDIF

    CALL monGATHER(ngrid,ncount,index, pfract,zfract)
    CALL monGATHER(ngrid,ncount,index, pmu,   zmu)
    CALL monGATHER(ngrid,ncount,index, albedo,zalb)
    DO l=1,nlayer+1
       CALL monGATHER(ngrid,ncount,index, plevel(:,l),zplev(:,l))
    ENDDO

#endif

    !$acc kernels default(none)

    flux_in(:)=0.
    flux_down(:,:)=0.
    flux_up(:,:)=0.
    ZDTSW(:,:)=0.
    zu(:,:)=0.


    !-----------------------------------------------------------------------
    !   calcul des profondeurs optiques integres depuis p=0:
    !----------------------------------------------------

    !-----------------------------------------------------------------------
    !   calcul de la transmission depuis le sommet de l atmosphere:
    !   -----------------------------------------------------------

    DO ig=1,NCOUNT
       flux_in(ig) = psolarf0*ZFRACT(ig)*ZMU(ig)
    ENDDO

    ! calcul de la partie homogene de l opacite
    tau0=-.5*log(coefvis)/ps_rad

    DO l=1,nlayer+1
       DO ig=1,NCOUNT
          zu(ig,l)=tau0*ZPLEV(ig,l)
          !          flux_down(ig,l) = flux_in(ig)*exp(-zu(ig,l)/ZMU(ig))
          flux_down(ig,l) = flux_in(ig)*exp(-zu(ig,l)/ABS(ZMU(ig)+1e-13))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   4. calcul du flux solaire arrivant sur le sol:
    !   ----------------------------------------------

    DO ig=1,NCOUNT
       ZFLUX(ig)     = (1.-ZALB(ig))*flux_down(ig,1) ! absorbed (net)
       flux_up(ig,1) =      ZALB(ig)*flux_down(ig,1) ! reflected (up)
    ENDDO

    !-----------------------------------------------------------------------
    !   5.calcul des transmissions depuis le sol, cas diffus:
    !   ------------------------------------------------------

    DO l=2,nlayer+1
       DO ig=1,NCOUNT
          flux_up(ig,l)=flux_up(ig,1)*exp(-(zu(ig,1)-zu(ig,l))*1.66)
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   3. taux de chauffage, ray. solaire direct:
    !   ------------------------------------------

    DO l=1,nlayer
       DO ig=1,NCOUNT
          ! m.cp.dT = dflux/dz
          ! m = -(dp/dz)/g
          ZDTSW(ig,l)=(g/cpp) &
               &     * (flux_down(ig,l+1)-flux_down(ig,l)) &
               &     / (ZPLEV(ig,l)-ZPLEV(ig,l+1))
       ENDDO
    ENDDO

    !-----------------------------------------------------------------------
    !   6.ajout a l echauffement de la contribution du ray. sol. reflechit:
    !   -------------------------------------------------------------------

    DO l=1,nlayer
       DO ig=1,NCOUNT
          ZDTSW(ig,l) = ZDTSW(ig,l) +                           &
               &        (g/cpp)*(flux_up(ig,l)-flux_up(ig,l+1)) &
               &        /(ZPLEV(ig,l)-ZPLEV(ig,l+1))
       ENDDO
    ENDDO

    !$acc end kernels

    !------------------------
    !   Diagnostics
    !------------------------

    IF (lverbose) THEN
       igout=NCOUNT/2+1

       !$acc update host(ZFRACT, ZMU, ZALB, ZFLUX, ZDTSW, flux_in, flux_down, flux_up)
       WRITELOG(*,*) 'Diagnostique des transmission dans le spectre solaire'
       WRITELOG(*,*) 'zfract, zmu, zalb, flux_in'
       WRITELOG(*,*) ZFRACT(igout), ZMU(igout), ZALB(igout), flux_in(igout)
       IF (flux_in(igout)>0.) THEN
          WRITELOG(*,*) 'Pression, quantite d abs, transmission'
          DO l=1,nlayer+1
             WRITELOG(*,*) zplev(igout,l),zu(igout,l),flux_down(igout,l)/flux_in(igout)
          ENDDO
       END IF
       LOG_INFO('rad_sw')

       WRITELOG(*,*) 'Diagnostique des taux de chauffage solaires:'
       WRITELOG(*,*) ' 2 flux solaire net incident sur le sol'
       WRITELOG(*,*) zflux(igout)
       LOG_INFO('rad_sw')

       IF (flux_up(igout,1)>0.) THEN
          WRITELOG(*,*) 'Diagnostique des taux de chauffage solaires'
          WRITELOG(*,*) ' 3 transmission avec les sol'
          WRITELOG(*,*) 'niveau     transmission'
          DO l=1,nlayer+1
             WRITELOG(*,*) l, flux_up(igout,l)/flux_up(igout,1)
          ENDDO
          LOG_INFO('rad_sw')
       END IF

       WRITELOG(*,*) 'Diagnostique des taux de chauffage solaires:'
       WRITELOG(*,*) ' 3 taux de chauffage total'
       DO l=1,nlayer
          WRITELOG(*,*) ZDTSW(igout,l)
       ENDDO
       LOG_INFO('rad_sw')
    ENDIF

    !----------------------------
    !   Outputs
    !----------------------------

#ifdef _OPENACC
    ! with OpenACC we work directly in fsrfvis and dtsw

    IF(lwrite) THEN
       !$acc update host(flux_in, flux_down,flux_up)
       CALL writefield('swtop','SW down TOA','W/m2',flux_in)
       CALL writefield('swflux_down','Downward SW flux','W/m2',flux_down)
       CALL writefield('swflux_up','Upward SW flux','W/m2',flux_up)
    END IF

#else
    IF(lwrite) THEN
       CALL monscatter(ngrid,ncount,index, flux_in,buf1)
       CALL writefield('swtop','SW down TOA','W/m2',buf1)

       DO l=1,nlayer+1
          CALL monscatter(ngrid,ncount,index, flux_down(:,l),buf2(:,l))
       ENDDO
       CALL writefield('swflux_down','Downward SW flux','W/m2',buf2)

       DO l=1,nlayer+1
          CALL monscatter(ngrid,ncount,index, flux_up(:,l),buf2(:,l))
       ENDDO
       CALL writefield('swflux_up','Upward SW flux','W/m2',buf2)
    END IF

    CALL monscatter(ngrid,ncount,index, zflux,fsrfvis)
    DO l=1,nlayer
       CALL monscatter(ngrid,ncount,index, zdtsw(:,l),dtsw(:,l))
    ENDDO
#endif

    !$acc end data

  END SUBROUTINE sw

END MODULE radiative_sw
