module radiative_sw
  implicit none
  save
  private
  public :: sw
contains
  subroutine sw (ngrid, nlayer, Diurn, CoefVis, Albedo, Pint, pS_rad,pMu, pFract, pSolarf0, FsrFvis, dTsw)
    !=======================================================================
    !
    !   Rayonnement solaire en atmosphere non diffusante avec un
    !   coefficient d absorption gris.
    !
    !=======================================================================
    USE phys_const,     only : Cpp, g
    USE writefield_mod, only : writefield
    integer,                         intent(in) :: ngrid    ! number of columns
    integer,                         intent(in) :: nlayer   ! number of vertical layers
    real,                            intent(in) :: pSolarf0 ! solar constant
    real,                            intent(in) :: CoefVis  ! CoefVis = attenuation at p = Ps_rad
    real,                            intent(in) :: Ps_rad   !
    real, dimension(ngrid),          intent(in) :: Albedo   ! Albedo
    real, dimension(ngrid),          intent(in) :: pMu      ! cosine of solar zenithal angle
    real, dimension(ngrid),          intent(in) :: pFract   ! day fraction
    real, dimension(ngrid,nlayer+1), intent(in) :: Pint     ! interface pressures
    logical,                         intent(in) :: Diurn    ! diurnal cycle

    ! Output
    real, dimension(ngrid),          intent(out) :: FsrfVis ! net surface flux
    real, dimension(ngrid,nlayer),   intent(out) :: dTsw    ! temperature tendency

    ! Fluxes are non-zero only on those ncount points where sun shines (Mu0 > 0)
    integer                            :: ncount
    integer, dimension(ngrid)          :: index

    integer                         :: ig, l
    real                            :: tau0
    real, dimension(ngrid)          :: zAlb      ! Albedo
    real, dimension(ngrid)          :: zMu       ! cosine zenithal angle
    real, dimension(ngrid)          :: zFract    ! day fraction
    real, dimension(ngrid)          :: Flux_in   ! incoming solar flux
    real, dimension(ngrid)          :: zFlux     ! net surface flux
    real, dimension(ngrid,nlayer)   :: zdTsw     ! temperature tendency
    real, dimension(ngrid,nlayer+1) :: Flux_down ! downward flux
    real, dimension(ngrid,nlayer+1) :: Flux_up   ! upward flux
    real, dimension(ngrid,nlayer+1) :: zPint     ! pressure at interfaces
    real, dimension(ngrid,nlayer+1) :: zU
    
    ! Count number day cells 
    if (Diurn) then
       ncount = 0
       index  = 0
       do ig = 1, ngrid
          if (pFract(ig) > 1e-6) then
             ncount = ncount + 1
             index (ncount) = ig
          end if
       end do
    else
       ncount = ngrid
    end  if

    call mongather (ngrid, ncount, index, pFract, zFract)
    call mongather (ngrid, ncount, index, pMu,    zMu)
    call mongather (ngrid, ncount, index, Albedo, zAlb)
    do l=  1, nlayer+1
       call mongather (ngrid, ncount, index, Pint(:,l), zPint(:,l))
    end do

    Flux_in   = 0.0
    Flux_down = 0.0
    Flux_up   = 0.0
    zdTsw     = 0.0
    zU        = 0.0

    ! Profondeurs optiques integres depuis p=0:

    ! Transmission depuis le sommet de l atmosphere:
    Flux_in(1:ncount) = pSolarf0 * zFract(1:ncount) * zMu(1:ncount)

    ! Partie homogene de l opacite
    tau0 = -0.5 * log (CoefVis) / Ps_rad

    zU(1:ncount,:) = tau0 * zPint(1:ncount,:)

    do l = 1, nlayer+1
       Flux_down(1:ncount,l) = Flux_in(1:ncount) * exp (- zU(1:ncount,l) / abs (zMu(1:ncount) + 1e-13))
    end do

    ! Flux solaire arrivant sur le sol:
    zFlux(1:ncount)     = (1.0 - zAlb(1:ncount)) * Flux_down(1:ncount,1) ! absorbed (net)
    Flux_up(1:ncount,1) =        zAlb(1:ncount)  * Flux_down(1:ncount,1) ! reflected (up)

    ! Transmissions depuis le sol (cas diffus)
    do l = 2, nlayer+1
       flux_up(1:ncount,l) = Flux_up(1:ncount,1) * exp (- (zU(1:ncount,1) - zU(1:ncount,l)) * 1.66)
    end do

    ! Taux de chauffage, ray. solaire direct
    ! m cp dt = dflux/dz, m = -(dp/dz)/g
    do l = 1, nlayer
       zdTsw(1:ncount,l) = (g/cpp) * (Flux_down(1:ncount,l+1) - Flux_down(1:ncount,l))  &
            / (zPint(1:ncount,l) - zPint(1:ncount,l+1))
    end do
    
    ! Ajout l echauffement de la contribution du ray. sol. reflechit:
    do l = 1, nlayer
       zdTsw(1:ncount,l) = zdTsw(1:ncount,l) &
            + (g/Cpp) * (Flux_up(1:ncount,l) - Flux_up(1:ncount,l+1)) / (zPint(1:ncount,l) - zPint(1:ncount,l+1))
    end do
    
    call monscatter (ngrid, ncount, index, zflux, FsrfVis)
    do l = 1, nlayer
       call monscatter (ngrid, ncount, index, zdTsw(:,l), dTsw(:,l))
    end do
  end subroutine sw
  
  pure subroutine mongather (ngrid, n, index, a, b) 
    ! Input:
    integer,                   intent(in)  :: ngrid, n
    integer, dimension(n),     intent(in)  :: index
    real,    dimension(ngrid), intent(in)  :: a

    ! Output:
    real,    dimension(n),     intent(out) :: b(n)

    integer :: i

    if (n < ngrid) then
       do i = 1, n
          b(i) = a(index(i))
       end do
    else
       b = a
    end if
  end subroutine mongather

  pure subroutine monscatter (ngrid, n, index, b, a) 
    ! Input
    integer,                   intent(in)  :: ngrid, n
    integer, dimension(n),     intent(in)  :: index
    real, dimension(n),        intent(in)  :: b

    ! Output
    real, dimension(ngrid),    intent(out) :: a

    integer :: i

    if (n < ngrid) then
       a = 0.0
       do i=1,n
          a(index(i)) = b(i)
       end do
    else
       a = b
    end if
  end subroutine monscatter
end module radiative_sw
