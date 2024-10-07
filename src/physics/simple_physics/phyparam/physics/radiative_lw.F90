module radiative_lw
  implicit none
  save
  private
  public :: lw
  logical, parameter :: lstrong =.true.
  real,    parameter :: Stephan = 5.67e-08
contains
  subroutine lw (ngrid, nlayer, CoefIR, Emissiv, pP, Ps_rad, pTsurf, pT, pFluxIR, pdTlw)
    use phys_const,     only : Cpp, G
    !===============================================================================================
    !  Cooling rate due to long wave (infrared) radiation
    !
    !  To simplify, the absorpotion coefficients only depend only on height and are pre-calculated.
    !
    !===============================================================================================

    ! Input
    integer,                         intent(in) :: ngrid    ! number of columns
    integer,                         intent(in) :: nlayer   ! number of layers
    real,                            intent(in) :: CoefIR
    real,                            intent(in) :: Ps_rad
    real, dimension(ngrid),          intent(in) :: Emissiv  ! Emissivity of surface
    real, dimension(ngrid),          intent(in) :: pTsurf   ! surface temperature
    real, dimension(ngrid,nlayer),   intent(in) :: pT       ! layer temperatures
    real, dimension(ngrid,nlayer+1), intent(in) :: pP       ! interface pressures

    ! Output
    real, dimension(ngrid),          intent(out) :: pFluxIR ! infrared flux through surface
    real, dimension(ngrid,nlayer),   intent(out) :: pdTlw   ! cooling rate

    integer                         :: ig, il, ilev
    real                            :: zcoef
    real, dimension(ngrid)          :: dup, lwtr1, lwtr2, Flux
    real, dimension(ngrid,nlayer+1) :: zPlanck
    real, dimension(ngrid,nlayer+1) :: zFluxup, zFluxdn, zFlux, zup

    ! Absorption coefficients
    if (lstrong) then ! strong absorption
       zup = pP**2 / (2.0 * g)
       zcoef = - log (CoefIR) / sqrt (Ps_rad**2 / (2.0 * g))
    else              ! weak absorption
       zup = pP
       zcoef = - log (CoefIR) / Ps_rad
    end if

    ! Black body function
    zPlanck(:,1:nlayer) = Stephan * pT**4

    ! Downwards fluxes
    !
    ! Downwards flux at interface ilev is a sum of contributions from layers above it
    ! each contribution depends on the layer itself (zPlanck) and its lower (lwtr1) and upPer (lwtr2) interfaces
    do ilev = 1, nlayer
       Flux  = 0.0
       dup   = zup(:,ilev) - zup(:,nlayer)
       lwtr1 = exp ( - zcoef * sqrt (dup))
       do il = nlayer-1, ilev, -1
          dup  = zup(:,ilev) - zup(:,il)
          lwtr2 = lwtr1
          lwtr1 = exp ( - zcoef * sqrt (dup))
          Flux  = Flux + zPlanck(:,il) * (lwtr1 - lwtr2)
       end do
       zFluxdn(:,ilev) = Flux
    end do
    zfluxdn(:,nlayer+1) = 0.0
    
    pFluxIr = Emissiv * zfluxdn(:,1)

    ! Upwards fluxes
    
    ! Upwards lw flux at the surface (ilev=1)
    zFluxup(:,1) = Emissiv * Stephan *  pTsurf**4 + (1.0 - Emissiv) * zFluxdn(:,1)

    ! Upwards flux at interface ilev>1 equals the surface flux times an absorption coefficient
    ! plus a sum of contributions from layers below it (il<ilev).
    ! each contribution depends on the layer itself (zPlanck) and its lower (lwtr2) and upPer (lwtr1) interfaces
    do ilev = 2, nlayer + 1
       dup   = zup(:,1) - zup(:,ilev)
       lwtr1 = exp ( - zcoef * sqrt (dup))
       Flux  = zfluxup(:,1) * lwtr1
       do il = 1, ilev - 1
          dup   = zup(:,il+1 ) - zup(:,ilev)
          lwtr2 = lwtr1
          lwtr1 = exp ( - zcoef * sqrt (dup))
          Flux  = Flux + zPlanck(:,il) * (lwtr1 - lwtr2)
       end do
       zFluxup(:,ilev) = Flux
    end do

    ! Net fluxes
    zFlux = zFluxup - zFluxdn

    ! Cooling rate
    pdTlw(:,1:nlayer) = (zFlux(:,2:nlayer+1) - zFlux(:,1:nlayer)) * G / (Cpp * (pP(:,2:nlayer+1) - pP(:,1:nlayer)))
  end subroutine lw

  pure subroutine lwtr (ngrid, coef, lstrong, dup, transm)
    ! Input
    integer,                intent(in) :: ngrid
    real,                   intent(in) :: coef
    logical,                intent(in) :: lstrong
    real, dimension(ngrid), intent(in) :: dup

    ! Output
    real, dimension(ngrid), intent(out) :: transm

    if (lstrong) then
       transm = exp( - coef* sqrt (dup))
    else
       transm = exp( - coef * dup)
    end if
  end subroutine lwtr
end module radiative_lw
