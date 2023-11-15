MODULE radiative_lw

#include "use_logging.h"

  IMPLICIT NONE
  SAVE

  PRIVATE

  PUBLIC :: lw

  LOGICAL, PARAMETER :: lstrong=.TRUE.
  REAL,    PARAMETER :: stephan=5.67e-08

  ! tiling
  INTEGER, PARAMETER :: tile = 128
#define TILE_OUTER(ig_start,ig_end, ngrid,tile) ig_start=1,(ngrid),(tile) ; ig_end=min(ig_start+(tile)-1, (ngrid))
#define TILE_INNER(ig,ig_loc, ig_start,ig_end) ig=ig_start,ig_end ; ig_loc=ig-ig_start

CONTAINS

  SUBROUTINE lw(ngrid,nlayer,coefir,emissiv, &
       pp,ps_rad,ptsurf,pt,              &
       pfluxir,pdtlw,                    &
       lverbose, lwrite)
    USE phys_const, ONLY : cpp, g
    USE writefield_mod, ONLY : writefield
    !=======================================================================
    !
    !   calcul de l evolution de la temperature sous l effet du rayonnement
    !   infra-rouge.
    !   Pour simplifier, les transmissions sont precalculees et ne
    !   dependent que de l altitude.
    !
    !   arguments:
    !   ----------
    !
    !   entree:
    !   -------
    !      ngrid             nombres de points de la grille horizontale
    !      nlayer              nombre de couches
    !      ptsurf(ngrid)     temperature de la surface
    !      pt(ngrid,nlayer)    temperature des couches
    !      pp(ngrid,nlayer+1)  pression entre les couches
    !      lwrite            variable logique pour sorties
    !
    !   sortie:
    !   -------
    !      pdtlw(ngrid,nlayer) taux de refroidissement
    !      pfluxir(ngrid)    flux infrarouge sur le sol
    !
    !=======================================================================

    INTEGER, INTENT(IN)  :: ngrid,nlayer
    REAL,    INTENT(IN)  :: coefir, ps_rad, &
         &                  emissiv(ngrid), ptsurf(ngrid), &
         &                  pt(ngrid,nlayer), pp(ngrid,nlayer+1)
    REAL,    INTENT(OUT) :: pdtlw(ngrid,nlayer), pfluxir(ngrid)
    LOGICAL, INTENT(IN)  :: lwrite, lverbose

    CHARACTER(6), PARAMETER :: tag='rad/lw' ! tag used in LOG_XXX

    ! local arrays => !$acc data create
    REAL :: zplanck(ngrid,nlayer+1), &
         &  zfluxup(ngrid,nlayer+1), &
         &  zfluxdn(ngrid,nlayer+1), &
         &  zflux(ngrid,nlayer+1),   &
         &  zup(ngrid,nlayer+1)

    ! tiles storing intermediate values in deep loop nests => !$acc private
    REAL :: dup(0:tile-1), lwtr1(0:tile-1), lwtr2(0:tile-1), flx(0:tile-1)

    ! scalars
    INTEGER :: nlevel, ilev, il
    INTEGER :: ig, ig_start, ig_end, ig_loc
    REAL :: zcoef

    !-----------------------------------------------------------------------
    !   initialisations:
    !   ----------------

    nlevel=nlayer+1

    !$acc data copyin(emissiv, ptsurf, pt, pp) &
    !$acc      create(zplanck, zfluxup, zfluxdn, zflux, zup) &
    !$acc      copyout(pdtlw, pfluxir)

    !-----------------------------------------------------------------------
    !   2. calcul des quantites d absorbants:
    !   -------------------------------------

    IF(lstrong) THEN  !   absorption forte

       !$acc parallel loop collapse(2) default(none) async
       DO ilev=1,nlevel
          DO ig=1,ngrid
             zup(ig,ilev)=pp(ig,ilev)*pp(ig,ilev)/(2.*g)
          ENDDO
       ENDDO
       !$acc end parallel loop

       IF(lverbose) THEN
          !$acc wait
          !$acc update self(zup)
          DO ilev=1,nlayer
             WRITELOG(*,*) ' up(',ilev,')  =  ',zup(ngrid/2+1,ilev)
          ENDDO
          LOG_DBG(tag)
       ENDIF

       zcoef=-log(coefir)/sqrt(ps_rad*ps_rad/(2.*g))

    ELSE  !   absorption faible

       !$acc parallel loop collapse(2) async
       DO ilev=1,nlevel
          DO ig=1,ngrid
             zup(ig,ilev)=pp(ig,ilev)
          ENDDO
       ENDDO
       !$acc end parallel loop

       zcoef=-log(coefir)/ps_rad
    ENDIF


    !-----------------------------------------------------------------------
    !   2. calcul de la fonction de corps noir:
    !   ---------------------------------------

    !$acc parallel loop collapse(2) async

    DO ilev=1,nlayer
       DO ig=1,ngrid
          zplanck(ig,ilev)=pt(ig,ilev)*pt(ig,ilev)
          zplanck(ig,ilev)=stephan*zplanck(ig,ilev)*zplanck(ig,ilev)
       ENDDO
    ENDDO
    !$acc end parallel loop


    !-----------------------------------------------------------------------
    !   4. flux descendants:
    !   --------------------

    ! Downwards flux at interface ilev is a sum of contributions from layers above it
    ! each contribution depends on the layer itself (zplanck) and its lower (lwtr1) and upper (lwtr2) interfaces

    !$acc parallel default(none) vector_length(tile) async
    !$acc loop gang collapse(2) private(dup, lwtr1, lwtr2, flx)
    DO ilev=1,nlayer
       DO TILE_OUTER(ig_start,ig_end, ngrid,tile)
          !DIR$ SIMD
          !$acc loop vector
          DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
             flx(ig_loc)=0.
             dup(ig_loc)   = zup(ig,ilev)-zup(ig,nlayer)
             lwtr1(ig_loc) = exp(-zcoef*sqrt(dup(ig_loc)))
          END DO

          !$acc loop seq
          DO il=nlayer-1,ilev,-1
             !DIR$ SIMD
             !$acc loop vector
             DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
                lwtr2(ig_loc) = lwtr1(ig_loc)
                dup(ig_loc)   = zup(ig,ilev)-zup(ig,il)
                lwtr1(ig_loc) = exp(-zcoef*sqrt(dup(ig_loc)))
                flx(ig_loc)   = flx(ig_loc) + zplanck(ig,il)*(lwtr1(ig_loc)-lwtr2(ig_loc))
             ENDDO
          ENDDO

          !DIR$ SIMD
          !$acc loop vector
          DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
             zfluxdn(ig,ilev) = flx(ig_loc)
          END DO
       ENDDO
    ENDDO
    !$acc end parallel

    !$acc parallel loop async
    DO ig=1,ngrid
       zfluxdn(ig,nlevel)=0.
       pfluxir(ig)=emissiv(ig)*zfluxdn(ig,1)
    ENDDO
    !$acc end parallel loop


    !-----------------------------------------------------------------------
    !   3. flux montants:
    !   ------------------

    ! Upwards lw flux at the surface (ilev=1)
    !$acc parallel loop async
    DO ig=1,ngrid
       zfluxup(ig,1)=ptsurf(ig)*ptsurf(ig)
       zfluxup(ig,1)=emissiv(ig)*stephan*zfluxup(ig,1)*zfluxup(ig,1) &
            +(1.-emissiv(ig))*zfluxdn(ig,1)
    ENDDO
    !$acc end parallel loop

    ! Upwards flux at interface ilev>1 equals the surface flux times an absorption coefficient
    ! plus a sum of contributions from layers below it (il<ilev).
    ! Each contribution depends on the layer itself (zplanck) and its lower (lwtr2) and upper (lwtr1) interfaces

    !$acc parallel default(none) vector_length(tile) async
    !$acc loop gang collapse(2) private(dup, lwtr1, lwtr2, flx)
    DO ilev=2,nlevel
       DO TILE_OUTER(ig_start,ig_end, ngrid,tile)
          !DIR$ SIMD
          !$acc loop vector
          DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
             dup(ig_loc)   = zup(ig,1)-zup(ig,ilev)
             lwtr1(ig_loc) = exp(-zcoef*sqrt(dup(ig_loc)))
             flx(ig_loc)   = zfluxup(ig,1)*lwtr1(ig_loc)
          END DO
          !$acc loop seq
          DO il=1,ilev-1
             !DIR$ SIMD
             !$acc loop vector
             DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
                lwtr2(ig_loc) = lwtr1(ig_loc)
                dup(ig_loc)   = zup(ig,il+1)-zup(ig,ilev)
                lwtr1(ig_loc) = exp(-zcoef*sqrt(dup(ig_loc)))
                flx(ig_loc)   = flx(ig_loc) + zplanck(ig,il)*(lwtr1(ig_loc)-lwtr2(ig_loc))
             ENDDO
          ENDDO
          !DIR$ SIMD
          !$acc loop vector
          DO TILE_INNER(ig,ig_loc, ig_start,ig_end)
             zfluxup(ig,ilev) = flx(ig_loc)
          END DO
       ENDDO
    ENDDO
    !$acc end parallel

    !-----------------------------------------------------------------------
    !   5. calcul des flux nets:
    !   ------------------------

    !$acc parallel loop collapse(2) async
    DO ilev=1,nlevel
       DO ig=1,ngrid
          zflux(ig,ilev)=zfluxup(ig,ilev)-zfluxdn(ig,ilev)
       ENDDO
    ENDDO
    !$acc end parallel loop

    !-----------------------------------------------------------------------
    !   6. Calcul des taux de refroidissement:
    !   --------------------------------------

    !$acc parallel loop collapse(2) async
    DO ilev=1,nlayer
       DO ig=1,ngrid
          pdtlw(ig,ilev)=(zflux(ig,ilev+1)-zflux(ig,ilev))* &
               g/(cpp*(pp(ig,ilev+1)-pp(ig,ilev)))
       ENDDO
    ENDDO
    !$acc end parallel loop

    !-----------------------------------------------------------------------
    !   10. sorties eventuelles:
    !   ------------------------

    !$acc wait

    IF (lverbose) THEN
       !$acc update self(ptsurf, pt, pdtlw, zfluxup, zfluxdn)
       WRITELOG(*,*) 'Diagnostique rayonnement thermique'
       WRITELOG(*,*) 'temperature     ', &
            'flux montant    flux desc.     taux de refroid.'
       ig=ngrid/2+1
       WRITELOG(6,'(4e18.4)') ptsurf(ig)
       DO ilev=1,nlayer
          WRITELOG(6,'(i4,4e18.4)') ilev,pt(ig,ilev), &
               zfluxup(ig,ilev),zfluxdn(ig,ilev),pdtlw(ig,ilev)
       ENDDO
       WRITELOG(6,'(4e18.4)') zfluxup(ig,nlevel),zfluxdn(ig,nlevel)
       LOG_DBG(tag)
    ENDIF

    IF(lwrite) THEN
       !$acc update self(zfluxup, zfluxdn)
       CALL writefield('lwflux_up', 'Upward LW flux', 'W/m2', zfluxup)
       CALL writefield('lwflux_down', 'Downward LW flux', 'W/m2', zfluxdn)
    END IF
    !-----------------------------------------------------------------------

    !$acc end data

  END SUBROUTINE lw

  PURE SUBROUTINE lwtr(ngrid,coef,lstrong,dup,transm)
    INTEGER, INTENT(IN) :: ngrid
    REAL,    INTENT(IN) :: coef
    LOGICAL, INTENT(IN) :: lstrong
    REAL,    INTENT(IN) :: dup(ngrid)
    REAL,    INTENT(OUT) :: transm(ngrid)
    INTEGER ig
    IF(lstrong) THEN
       DO ig=1,ngrid
          transm(ig)=exp(-coef*sqrt(dup(ig)))
       ENDDO
    ELSE
       DO ig=1,ngrid
          transm(ig)=exp(-coef*dup(ig))
       ENDDO
    ENDIF

  END SUBROUTINE lwtr

END MODULE radiative_lw
