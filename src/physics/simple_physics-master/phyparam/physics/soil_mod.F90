MODULE soil_mod

#include "use_logging.h"

  IMPLICIT NONE
  PRIVATE
  SAVE

  REAL, PARAMETER :: min_period=20000., dalph_soil=2.,  &
       pi=2.*ASIN(1.), fz1=SQRT(min_period/pi)

  ! common variables
  REAL, PUBLIC ::  I_mer,I_ter,Cd_mer,Cd_ter, &
       &           alb_mer,alb_ter,emi_mer,emi_ter

  ! precomputed variables
  REAL :: lambda
  REAL, ALLOCATABLE :: dz1(:),dz2(:)
  !$OMP THREADPRIVATE(dz1,dz2)
  !$acc declare create(dz1,dz2)

  REAL, ALLOCATABLE :: rnatur(:), albedo(:),emissiv(:), z0(:), inertie(:)
  !$OMP THREADPRIVATE( rnatur, albedo, emissiv, z0, inertie)
  !$acc declare create(rnatur, albedo, emissiv, z0, inertie)

  ! internal state, written to / read from disk at checkpoint / restart
  REAL, ALLOCATABLE :: tsurf(:), tsoil(:,:)
  !$OMP THREADPRIVATE(tsurf, tsoil)
  !$acc declare create(tsurf,tsoil)

  PUBLIC :: init_soil, soil_forward, soil_backward, &
       rnatur, albedo, emissiv, z0, inertie, &
       tsurf, tsoil

CONTAINS

  FUNCTION fz(rk) RESULT(val)
    !$acc routine seq
    REAL :: val, rk
    val = fz1*(dalph_soil**rk-1.)/(dalph_soil-1.)
  END FUNCTION fz

  SUBROUTINE init_soil(nsoil)
    INTEGER, INTENT(IN) :: nsoil
    REAL :: rk, rk1, rk2
    INTEGER :: jk

    !-----------------------------------------------------------------------
    !   ground levels
    !   grnd=z/l where l is the skin depth of the diurnal cycle:
    !   --------------------------------------------------------

    WRITELOG(*,*) 'nsoil,firstcall=',nsoil, .TRUE.

    ALLOCATE(dz1(nsoil-1),dz2(nsoil))

    !$acc kernels default(none) present(dz1, dz2)
    !$acc loop private(rk1, rk2)
    DO jk=1,nsoil
       rk1=jk
       rk2=jk-1
       dz2(jk)=fz(rk1)-fz(rk2)  !numerator of c_k+0.5 A.11
    ENDDO
    !$acc loop private(rk1, rk2)
    DO jk=1,nsoil-1
       rk1=jk+.5
       rk2=jk-.5
       dz1(jk)=1./(fz(rk1)-fz(rk2))  ! d_k A.12 in Frederics these
    ENDDO
    !$acc end kernels

    !$acc update host(dz1, dz2)
    lambda=fz(.5)*dz1(1)  ! mu (A.28) in Frederics these

    WRITELOG(*,*) 'full layers, intermediate layers (secoonds)'
    DO jk=1,nsoil
       rk=jk
       rk1=jk+.5
       rk2=jk-.5
       WRITELOG(*,*) fz(rk1)*fz(rk2)*pi,        &
            &        fz(rk)*fz(rk)*pi
    ENDDO

    LOG_INFO('init_soil')

  END SUBROUTINE init_soil

  !=======================================================================
  !
  !   Auteur:  Frederic Hourdin     30/01/92
  !   -------
  !
  !   objet:  computation of : the soil temperature evolution
  !   ------                   the surfacic heat capacity "Capcal"
  !                            the surface conduction flux pcapcal
  !
  !
  !   Method: implicit time integration
  !   -------
  !   Consecutive ground temperatures are related by:
  !           T(k+1) = C(k) + D(k)*T(k)  (1)
  !   the coefficients C and D are computed at the t-dt time-step.
  !   Routine structure:
  !   1)new temperatures are computed  using (1)
  !   2)C and D coefficients are computed from the new temperature
  !     profile for the t+dt time-step
  !   3)the coefficients A and B are computed where the diffusive
  !     fluxes at the t+dt time-step is given by
  !            Fdiff = A + B Ts(t+dt)
  !     or     Fdiff = F0 + Capcal (Ts(t+dt)-Ts(t))/dt
  !            with F0 = A + B (Ts(t))
  !                 Capcal = B*dt
  !

  PURE SUBROUTINE soil_backward(ngrid,nsoil, zc,zd, ptsrf,ptsoil)
    INTEGER, INTENT(IN) :: ngrid, nsoil         ! number of columns, of soil layers
    REAL, INTENT(IN)    :: zc(ngrid, nsoil-1), zd(ngrid, nsoil-1) ! LU factorization
    REAL, INTENT(IN)    :: ptsrf(ngrid)         ! new surface temperature
    REAL, INTENT(INOUT) :: ptsoil(ngrid,nsoil)  ! soil temperature
    INTEGER :: ig, jk

    !-----------------------------------------------------------------------
    !  Computation of the soil temperatures using the Cgrd and Dgrd
    !  coefficient computed during the forward sweep
    !  -----------------------------------------------

    !  surface temperature => temperature in first soil layer

    !$acc data copyin(zc(:,:),zd(:,:),ptsrf(:))        &
    !$acc &    copyout(ptsoil(:,:))                    &
    !$acc &

    !$acc kernels default(none)
    DO ig=1,ngrid
       ptsoil(ig,1)=(lambda*zc(ig,1)+ptsrf(ig))/                   &
             &      (lambda*(1.-zd(ig,1))+1.)   ! A.27 re-arragend to solve for T_0.5 in Frederics these
    ENDDO

    !   other temperatures
    DO jk=1,nsoil-1
       DO ig=1,ngrid
          ptsoil(ig,jk+1)=zc(ig,jk)+zd(ig,jk)*ptsoil(ig,jk) ! A.15 in Frederics these
       ENDDO
    ENDDO
    !$acc end kernels
    !$acc end data

  END SUBROUTINE Soil_backward

  PURE SUBROUTINE soil_forward(ngrid, nsoil, ptimestep, ptherm_i, ptsrf, ptsoil, &
       &                       zc, zd, pcapcal, pfluxgrd)

    INTEGER, INTENT(IN) :: ngrid, nsoil         ! number of columns, of soil layers
    REAL, INTENT(IN)    :: ptimestep,         & ! time step
         &                 ptherm_i(ngrid),   & ! thermal inertia ??
         &                 ptsrf(ngrid),      & ! surface temperature before heat conduction
         &                 ptsoil(ngrid, nsoil) ! soil temperature before heat conduction
    REAL, INTENT(OUT)   :: zc(ngrid,nsoil-1),   &
         &                 zd(ngrid, nsoil-1),  & ! LU factorization for backward sweep
         &                 pcapcal(ngrid),    & ! effective calorific capacity
         &                 pfluxgrd(ngrid)      ! conductive heat flux at the ground
    REAL :: z1, zdz2(nsoil)
    INTEGER :: jk, ig
   
    !-----------------------------------------------------------------------
    !   Computation of the Cgrd and Dgrd coefficients the backward sweep :
    !   ---------------------------------------------------------------

    !$acc data copyin(dz1, dz2, ptsrf(:),ptherm_i(:),ptsoil(:,:))     &
    !$acc &    copyout(zc(:,1:nsoil-1), zd(:,1:nsoil-1), pcapcal(:), pfluxgrd(:))               &
    !$acc &    create(zdz2(:))

    !$acc kernels default(none)

    DO jk=1,nsoil
       zdz2(jk)=dz2(jk)/ptimestep   ! c_k+0.5 A.11 in Frederics these
    ENDDO

    !$acc loop private(z1)
    DO ig=1,ngrid
       z1=zdz2(nsoil)+dz1(nsoil-1)
       zc(ig,nsoil-1)=zdz2(nsoil)*ptsoil(ig,nsoil)/z1 ! B__N-1 (A.17)
       zd(ig,nsoil-1)=dz1(nsoil-1)/z1                 ! a_N-1 (A.16) in Frederics these
    ENDDO

    !$acc loop private(z1)
    DO jk=nsoil-1,2,-1
       DO ig=1,ngrid
          z1=1./(zdz2(jk)+dz1(jk-1)+dz1(jk)*(1.-zd(ig,jk)))
          zc(ig,jk-1)=                                                &
          &      (ptsoil(ig,jk)*zdz2(jk)+dz1(jk)*zc(ig,jk))*z1  ! all the B_k ()
          zd(ig,jk-1)=dz1(jk-1)*z1 ! all the a_ks
       ENDDO
    ENDDO
    !-----------------------------------------------------------------------
    !   computation of the surface diffusive flux from ground and
    !   calorific capacity of the ground:
    !   ---------------------------------

    !$acc loop private(z1)
    DO ig=1,ngrid
       pfluxgrd(ig)=ptherm_i(ig)*dz1(1)*                              &
            &   (zc(ig,1)+(zd(ig,1)-1.)*ptsoil(ig,1))               ! F* A.25 in Frederics these
       z1=lambda*(1.-zd(ig,1))+1.
       pcapcal(ig)=ptherm_i(ig)*                                      &
           &   ptimestep*(zdz2(1)+(1.-zd(ig,1))*dz1(1))/z1         ! C_s A.30
       pfluxgrd(ig)=pfluxgrd(ig)                                      &
           &   +pcapcal(ig)*(ptsoil(ig,1)*z1-lambda*zc(ig,1)-ptsrf(ig))   &
           &   /ptimestep                                           ! F_s A.31
    ENDDO
    !$acc end kernels
    !$acc end data

  END SUBROUTINE soil_forward

END MODULE soil_mod
