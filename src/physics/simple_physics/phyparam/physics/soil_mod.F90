module soil_mod
  ! Vertical diffusion of heat in soil column
  !
  ! References refer to Frederic Hourdin thesis
  !
  !   Auteur:  Frederic Hourdin  1992-01-30
  !   -------
  !
  !   Computation of : 1. soil temperature evolution
  !                    2. surface heat capacity CapCal
  !                    3. surface conduction flux pCapCal
  !
  !
  !   Method: implicit time integration
  !   -------
  !
  !   Consecutive ground temperatures are related by:
  !
  !           T(k+1)  =  c(k) + d(k) * T(k)  (1)
  !
  !   the coefficients c and d are computed at the t-dt time-step.
  !   routine structure:
  !
  !   1. New temperatures are computed  using (1)
  !
  !   2. C and d coefficients are computed from the new temperature
  !     profile for the t+dt time-step
  !
  !   3. Coefficients a and b are computed where the diffusive
  !     fluxes at the t+dt time-step is given by
  !
  !            Fdiff  =  a + b Ts(t+dt)
  !     or     Fdiff  =  f0 + CapCal (Ts(t+dt) - Ts(t)) / dt
  !
  !          with f0  =  a + b (Ts(t))
  !           CapCal  =  b * dt
  !
  ! -----------------------------------------------------------------------------------------------------
  implicit none
  private
  save
  real, parameter :: min_period = 20000.0
  real, parameter :: dAlph_soil = 2.0
  real, parameter :: Pi         = 2.0 * asin (1.0)
  real, parameter :: fz1        = sqrt (min_period / Pi)

  ! Common variables
  real, public ::  I_mer, I_ter, Cd_mer, Cd_ter, Alb_mer, Alb_ter, Emi_mer, Emi_ter

  ! Precomputed variables
  real                              :: lambda
  real, dimension(:),  allocatable  :: dz1, dz2
  real, dimension(:),  allocatable  :: land, Albedo, Emissiv, Z0, pThermal_inertia

  ! Internal state, written to / read from disk at checkpoint / restart
  real, dimension(:),   allocatable :: Tsurf
  real, dimension(:,:), allocatable :: tsoil

  public :: init_soil, soil_forward, soil_backward, land, Albedo, Emissiv, Z0, pThermal_inertia, Tsurf, Tsoil
contains
  pure subroutine soil_forward (ngrid, nsoil, pTimestep, pThermal_inertia, pTsrf, pTsoil, zc, zd, pCapCal, pFluxGrd)
    ! Input
    integer,                        intent(in)  :: ngrid             ! number of columns
    integer,                        intent(in)  :: nsoil             ! number of soil layers
    real,                           intent(in)  :: pTimestep         ! time step
    real, dimension(ngrid),         intent(in)  :: pThermal_Inertia  ! thermal inertia
    real, dimension(ngrid),         intent(in)  :: pTsrf             ! surface temperature before heat conduction
    real, dimension(ngrid,nsoil),   intent(in)  :: pTsoil            ! soil temperature before heat conduction

    real, dimension(ngrid,nsoil-1), intent(out) :: zc, zd            ! LU factorization for backward sweep

    ! Output
    real, dimension(ngrid),         intent(out) :: pCapCal           ! effective calorific capacity
    real, dimension(ngrid),         intent(out) :: pFluxGrd          ! conductive heat flux at the ground

    integer                :: ig, jk
    real, dimension(nsoil) :: z1, zdz2

    ! Computation of the cgrd and dgrd coefficients the backward sweep :
    zdz2 = dz2 / ptimestep                             ! c_k + 0.5 (A.11)

    z1 = zdz2(nsoil) + dz1(nsoil-1)

    zc(:,nsoil-1) = zdz2(nsoil) * pTsoil(:,nsoil) / z1 ! b_n - 1 (A.17)
    zd(:,nsoil-1) = dz1(nsoil-1)                  / z1 ! a_n - 1 (A.16)
    do jk = nsoil-1, 2, -1
       z1 = 1.0 / (zdz2(jk) + dz1(jk-1) + dz1(jk) * (1.0 - zd(:,jk)))

       zc(:,jk-1) = z1 * (pTsoil(:,jk) * zdz2(jk) + dz1(jk) * zc(:,jk)) ! b_k
       zd(:,jk-1) = z1 * dz1(jk-1)                                      ! a_k
    end do

    ! Surface diffusive flux and calorific capacity of ground
    z1 = 1.0 + lambda * (1.0 - zd(:,1))

    pCapCal  = pThermal_inertia * pTimestep * (zdz2(1) + (1.0 - zd(:,1)) * dz1(1)) / z1        ! c_s (A.30)

    pFluxGrd = pThermal_inertia * dz1(1) * (zc(:,1) + (zd(:,1) - 1.0) * pTsoil(:,1))           ! f * (A.25)
    pFluxGrd = pFluxGrd + pCapCal * (pTsoil(:,1) * z1 - lambda * zc(:,1) - pTsrf) / pTimestep  ! f_s (A.31)
  end subroutine soil_forward

  pure subroutine soil_backward (ngrid, nsoil, zc, zd, pTsrf, pTsoil)
    ! Soil column temperatures using cgrd and dgrd  coefficient computed during the forward sweep

    ! Input
    integer,                         intent(in)    :: ngrid   ! number of columns
    integer,                         intent(in)    :: nsoil   ! number of soil layers
    real, dimension(ngrid, nsoil-1), intent(in)    :: zc, zd  ! Lu factorization
    real, dimension(ngrid),          intent(in)    :: pTsrf   ! new surface temperature

    ! Output
    real, dimension(ngrid,nsoil),    intent(inout) :: pTsoil  ! soil temperature updated from t -> t+dt

    integer :: ig, jk

    ! Soil column temperatures
    pTsoil(:,1) = (lambda * zc(:,1) + pTsrf) / (lambda * (1.0 - zd(:,1)) + 1.0) ! (A.27) re-arragend to solve for t_0.5
    do jk = 2, nsoil
       pTsoil(:,jk) = zc(:,jk-1) + zd(:,jk-1) * pTsoil(:,jk-1)                  ! (A.15)
    end do
  end subroutine soil_backward

  subroutine init_soil (nsoil)
    !   Ground levels
    !   grnd = z/l where l is the skin depth of the diurnal cycle:
    integer, intent(in) :: nsoil

    real    :: rk, rk1, rk2
    integer :: jk

    allocate (dz1(1:nsoil-1), dz2(1:nsoil))

    do jk = 1, nsoil-1
       rk1 = float (jk) + 0.5
       rk2 = float (jk) - 0.5
       dz1(jk) = 1.0 / (fz (rk1) - fz (rk2))  ! d_k (A.12)
    end do
    lambda = fz (0.5) * dz1(1)                ! mu (A.28)

    do jk = 1, nsoil
       rk1 = float (jk)
       rk2 = float (jk) - 1.0
       dz2(jk) = fz (rk1) - fz (rk2)          ! numerator of c_k + 0.5 (A.11)
    end do
  end subroutine init_soil

  real(8) function fz (rk)
    real :: rk

    fz = fz1 * (dAlph_soil**rk - 1.0) / (dAlph_soil - 1.0)
  end function fz
end module soil_mod
