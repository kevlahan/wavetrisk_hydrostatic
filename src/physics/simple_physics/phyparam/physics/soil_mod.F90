module soil_mod
  ! Vertical diffusion of heat in soil column
  !
  ! References refer to Frederic Hourdin's thesis
#include "use_logging.h"
  implicit none
  private
  save
  real, parameter :: min_period = 20000.0
  real, parameter :: dalph_soil = 2.0
  real, parameter :: pi         = 2.0 * asin (1.0)
  real, parameter :: fz1        = sqrt (min_period/pi)

  ! Common variables
  real, public ::  i_mer, i_ter, cd_mer, cd_ter, alb_mer, alb_ter, emi_mer, emi_ter

  ! Precomputed variables
  real                              :: lambda
  real, dimension(:),  allocatable  :: dz1, dz2
  real, dimension(:),  allocatable  :: rnatur, Albedo, emissiv, z0, inertie

  ! Internal state, written to / read from disk at checkpoint / restart
  real, dimension(:),   allocatable :: Tsurf
  real, dimension(:,:), allocatable :: tsoil

  public :: init_soil, soil_forward, soil_backward, rnatur, Albedo, Emissiv, z0, inertie, Tsurf, Tsoil
contains
  
  ! -----------------------------------------------------------------------------------------------------
  !
  !   Auteur:  Frederic Hourdin     30/01/92
  !   -------
  !
  !   Objet:  computation of : the soil temperature evolution
  !   ------                   the surfacic heat capacity "capcal"
  !                            the surface conduction flux pcapcal
  !
  !
  !   Method: implicit time integration
  !   -------
  !   consecutive ground temperatures are related by:
  !
  !           T(k+1)  =  c(k) + d(k) * T(k)  (1)
  !
  !   the coefficients c and d are computed at the t-dt time-step.
  !   routine structure:
  !
  !   1) New temperatures are computed  using (1)
  !
  !   2) C and d coefficients are computed from the new temperature
  !     profile for the t+dt time-step
  !
  !   3) The coefficients a and b are computed where the diffusive
  !     fluxes at the t+dt time-step is given by
  !
  !            fdiff  =  a + b ts(t+dt)
  !     or     fdiff  =  f0 + capcal (ts(t+dt)-ts(t))/dt
  !
  !          with f0  =  a + b (ts(t))
  !           CapCal  =  b * dt
  !
  ! -----------------------------------------------------------------------------------------------------
   pure subroutine soil_forward (ngrid, nsoil, ptimestep, ptherm_i, pTsrf, pTsoil, zc, zd, pcapcal, pfluxgrd)
    integer,                        intent(in)  :: ngrid     ! number of columns, of soil layers
    integer,                        intent(in)  :: nsoil     ! number of columns, of soil layers
    real,                           intent(in)  :: ptimestep ! time step
    real, dimension(ngrid),         intent(in)  :: ptherm_i  ! thermal inertia ??
    real, dimension(ngrid),         intent(in)  :: pTsrf     ! surface temperature before heat conduction
    real, dimension(ngrid,nsoil),   intent(in)  :: pTsoil    ! soil temperature before heat conduction

    real, dimension(ngrid,nsoil-1), intent(out) :: zc, zd    ! Lu factorization for backward sweep
    real, dimension(ngrid),         intent(out) :: pCapCal   ! effective calorific capacity
    real, dimension(ngrid),         intent(out) :: pFluxGrd  ! conductive heat flux at the ground
   
    integer                :: ig, jk
    real                   :: z1
    real, dimension(nsoil) :: zdz2

    ! Computation of the cGrd and dGrd coefficients the backward sweep :
    zdz2 = dz2 / ptimestep   ! c_k + 0.5 A.11 

    z1 = zdz2(nsoil) + dz1(nsoil-1)
    
    zc(:,nsoil-1) = zdz2(nsoil) * pTsoil(:,nsoil) / z1 ! b_n - 1 (A.17)
    zd(:,nsoil-1) = dz1(nsoil-1)                  / z1 ! a_n - 1 (A.16) 

    do jk = nsoil-1, 2, -1
       do ig = 1, ngrid
          z1 = 1.0 / (zdz2(jk) + dz1(jk-1) + dz1(jk) * (1.0 - zd(ig,jk)))
          
          zc(ig,jk-1) = z1 * (pTsoil(ig,jk) * zdz2(jk) + dz1(jk) * zc(ig,jk)) ! b_k
          zd(ig,jk-1) = z1 * dz1(jk-1)                                        ! a_k
       end do
    end do
    
    ! Surface diffusive flux and calorific capacity of ground
    do ig = 1, ngrid
       pFluxGrd(ig) = pTherm_i(ig) * dz1(1) * (zc(ig,1) + (zd(ig,1)-1.) * pTsoil(ig,1))                             ! f *  A.25
       
       z1 = lambda * (1.0 - zd(ig,1)) + 1.0
       
       pCapCal(ig)  = pTherm_i(ig) * ptimestep * (zdz2(1) + (1.0 - zd(ig,1)) * dz1(1)) / z1                         ! c_s A.30
       pFluxGrd(ig) = pFluxGrd(ig) + pCapCal(ig) * (pTsoil(ig,1) * z1 - lambda * zc(ig,1) - pTsrf(ig)) / ptimestep  ! f_s A.31
    end do
  end subroutine soil_forward

  pure subroutine soil_backward (ngrid, nsoil, zc, zd, pTsrf, pTsoil)
    ! Soil temperatures using cgrd and dgrd  coefficient computed during the forward sweep
    integer,                         intent(in)    :: ngrid   ! number of columns
    integer,                         intent(in)    :: nsoil   ! number of soil layers
    real, dimension(ngrid, nsoil-1), intent(in)    :: zc, zd  ! Lu factorization
    real, dimension(ngrid),          intent(in)    :: pTsrf   ! new surface temperature
    real, dimension(ngrid,nsoil),    intent(inout) :: pTsoil  ! soil temperature

    integer :: ig, jk

    pTsoil(:,1) = (lambda * zc(:,1) + pTsrf) / (lambda * (1.0 - zd(:,1)) + 1.0) ! A.27 re-arragend to solve for t_0.5
    pTsoil(:,2:nsoil) = zc(:,1:nsoil-1) + zd(:,1:nsoil-1) * pTsoil(:,1:nsoil-1) ! A.15 
  end subroutine soil_backward

  subroutine init_soil (nsoil)
    !   Ground levels
    !   grnd = z/l where l is the skin depth of the diurnal cycle:
    integer, intent(in) :: nsoil
    
    real    :: rk, rk1, rk2
    integer :: jk

    WRITELOG (*,*) 'nsoil, firstcall  =  ', nsoil, .true.

    allocate (dz1(nsoil-1), dz2(nsoil))

    do jk = 1, nsoil
       rk1 = jk
       rk2 = jk-1
       dz2(jk) = fz(rk1) - fz(rk2)          ! numerator of c_k + 0.5 A.11
    end do
    do jk = 1, nsoil-1
       rk1 = jk + 0.5
       rk2 = jk - 0.5
       dz1(jk) = 1.0 / (fz(rk1) - fz(rk2))  ! d_k A.12 
    end do
    lambda = fz(0.5) * dz1(1)               ! mu A.28

    WRITELOG (*,*) 'Full layers, intermediate layers (seconds)'
    
    do jk = 1, nsoil
       rk  = jk
       rk1 = jk + 0.5
       rk2 = jk - 0.5
       WRITELOG (*,*) fz(rk1) * fz(rk2) * pi, fz(rk) * fz(rk) * pi
    end do
    LOG_INFO ('init_soil')
  end subroutine init_soil

  function fz (rk) result (val)
    real :: val, rk

    val = fz1 * (dalph_soil**rk - 1.0) / (dalph_soil - 1.0)
  end function fz
end module soil_mod
