module radiative_sw
#include "use_logging.h"
  implicit none
  save
  private
  public :: sw
contains
  subroutine sw (ngrid, nlayer,ldiurn, CoefVis, Albedo, pLevel, pS_rad,pMu, pfract, pSolarf0, FsrFvis, dTsw, lverbose, lwrite)
    !=======================================================================
    !
    !   Rayonnement solaire en atmosphere non diffusante avec un
    !   coefficient d absorption gris.
    !
    !=======================================================================
    USE phys_const,     only : Cpp, g
    USE writefield_mod, only : writefield
    integer,                         intent(in) :: ngrid            ! number of columns
    integer,                         intent(in) :: nlayer           ! number of vertical layers
    real,                            intent(in) :: pSolarf0         ! solar constant
    real,                            intent(in) :: pS_rad, CoefVis  ! CoefVis = attenuation at p = pS_rad
    real, dimension(ngrid),          intent(in) :: Albedo           ! Albedo
    real, dimension(ngrid),          intent(in) :: pMu              ! cosine zenithal angle
    real, dimension(ngrid),          intent(in) :: pFract           ! day fraction
    real, dimension(ngrid,nlayer+1), intent(in) :: pLevel           ! pressure at interfaces
    logical,                         intent(in) :: ldiurn, lverbose, lwrite

    ! Output
    real, dimension(ngrid),          intent(out) :: FsrFvis         ! net surface flux
    real, dimension(ngrid,nlayer),   intent(out) :: dTsw            ! temperature tendency

    ! Fluxes are non-zero only on those points where the sun shines (Mu0 > 0)
    ! we compute only on those ncount points
    integer                            :: ncount
    integer, dimension(ngrid)          :: index
    real,    dimension(ngrid)          :: buf1
    real,    dimension(ngrid,nlayer+1) :: buf2

    integer                         :: ig,l,igout
    real                            :: tau0
    real, dimension(ngrid)          :: zAlb      ! Albedo
    real, dimension(ngrid)          :: zMu       ! cosine zenithal angle
    real, dimension(ngrid)          :: zFract    ! day fraction
    real, dimension(ngrid)          :: Flux_in   ! incoming solar flux
    real, dimension(ngrid)          :: zFlux     ! net surface flux
    real, dimension(ngrid,nlayer)   :: zdtsw     ! temperature tendency
    real, dimension(ngrid,nlayer+1) :: Flux_down ! downward flux
    real, dimension(ngrid,nlayer+1) :: Flux_up   ! upward flux
    real, dimension(ngrid,nlayer+1) :: zPlev     ! pressure at interfaces
    real, dimension(ngrid,nlayer+1) :: zU        !

    if (ldiurn) then
       ncount = 0
       index = 0
       do ig = 1, ngrid
          if (pfract(ig) > 1e-6) then
             ncount = ncount+1
             index (ncount) = ig
          end if
       end do
    else
       ncount = ngrid
    end  if

    call mongather (ngrid, ncount, index, pFract, zfract)
    call mongather (ngrid, ncount, index, pMu,    zMu)
    call mongather (ngrid, ncount, index, Albedo, zAlb)
    do l=  1, nlayer+1
       call mongather (ngrid, ncount, index, pLevel(:,l), zPlev(:,l))
    end do

    Flux_in   = 0.0
    Flux_down = 0.0
    Flux_up   = 0.0
    zdTsw     = 0.0
    zU        = 0.0

    !-----------------------------------------------------------------------
    !   calcul des profondeurs optiques integres depuis p=0:
    !----------------------------------------------------

    !-----------------------------------------------------------------------
    !   calcul de la transmission depuis le sommet de l atmosphere:
    !   -----------------------------------------------------------
    Flux_in(1:ncount) = pSolarf0 * zfract(1:ncount) * zMu(1:ncount)

    ! calcul de la partie homogene de l opacite
    tau0 = -0.5 * log (CoefVis) / pS_rad

    zU(1:ncount,:) = tau0 * zPlev(1:ncount,:)

    do l = 1, nlayer+1
       Flux_down(1:ncount,l) = Flux_in(1:ncount) * exp (- zU(1:ncount,l) / abs (zMu(1:ncount) + 1e-13))
    end do

    !-----------------------------------------------------------------------
    !   4. calcul du flux solaire arrivant sur le sol:
    !   ----------------------------------------------
    zFlux(1:ncount)     = (1.0 - zAlb(1:ncount)) * Flux_down(1:ncount,1) ! absorbed (net)
    Flux_up(1:ncount,1) =        zAlb(1:ncount)  * Flux_down(1:ncount,1) ! reflected (up)

    !-----------------------------------------------------------------------
    !   5.calcul des transmissions depuis le sol, cas diffus:
    !   ------------------------------------------------------
    do l = 2, nlayer+1
       flux_up(1:ncount,l) = Flux_up(1:ncount,1) * exp (- (zU(1:ncount,1) - zU(1:ncount,l)) * 1.66)
    end do

    !-----------------------------------------------------------------------
    !   3. taux de chauffage, ray. solaire direct:
    !   ------------------------------------------
    ! m.cp.dt = dflux/dz
    ! m = -(dp/dz)/g
    zdTsw(1:ncount,1:nlayer) = (g/cpp) * (Flux_down(1:ncount,2:nlayer+1) - Flux_down(1:ncount,1:nlayer))  / (zplev(1:ncount,1:nlayer) - zplev(1:ncount,2:nlayer+1))

    !-----------------------------------------------------------------------
    !   6. ajout a l echauffement de la contribution du ray. sol. reflechit:
    !   -------------------------------------------------------------------
    zdTsw(1:ncount,1:nlayer) = zdTsw(1:ncount,1:nlayer) + (g/Cpp) * (Flux_up(1:ncount,1:nlayer) - Flux_up(1:ncount,2:nlayer+1)) /( zPlev(1:ncount,1:nlayer) - zPlev(1:ncount,2:nlayer+1))


    !------------------------
    !   Diagnostics
    !------------------------
    if (lverbose) then
       igout = ncount/2 + 1

       WRITELOG (*,*) 'Diagnostique des transmission dans le spectre solaire'
       WRITELOG (*,*) 'zfract, zmu, zalb, flux_in'
       WRITELOG (*,*) ZFRACT(igout), ZMU(igout), ZALB(igout), flux_in(igout)

       if (flux_in(igout) > 0.0) then
          WRITELOG (*,*) 'Pression, quantite d abs, transmission'
          do l = 1, nlayer+1
             WRITELOG (*,*) zplev(igout,l),zu(igout,l),flux_down(igout,l)/flux_in(igout)
          end do
       end if
       LOG_INFO ('rad_sw')

       WRITELOG (*,*) 'Diagnostique des taux de chauffage solaires:'
       WRITELOG (*,*) ' 2 flux solaire net incident sur le sol'
       WRITELOG (*,*) zflux(igout)
       LOG_INFO ('rad_sw')

       if (flux_up(igout,1) > 0.0) then
          WRITELOG (*,*) 'Diagnostique des taux de chauffage solaires'
          WRITELOG (*,*) ' 3 transmission avec les sol'
          WRITELOG (*,*) 'niveau     transmission'
          do l=  1, nlayer+1
             WRITELOG (*,*) l, flux_up(igout,l) / flux_up(igout,1)
          end do
          LOG_INFO ('rad_sw')
       end if

       WRITELOG (*,*) 'Diagnostique des taux de chauffage solaires:'
       WRITELOG (*,*) ' 3 taux de chauffage total'
       do l = 1, nlayer
          WRITELOG(*,*) ZDTSW(igout,l)
       end do
       LOG_INFO( 'rad_sw')
    end if

    if (lwrite) then
       call monscatter (ngrid, ncount, index, flux_in, buf1)
       call writefield ('swtop','sw down toa','w/m2',buf1)

       do l = 1, nlayer+1
          call monscatter (ngrid, ncount, index, flux_down(:,l), buf2(:,l))
       end do
       call writefield ('swflux_down','downward sw flux','w/m2', buf2)

       do l = 1, nlayer+1
          call monscatter (ngrid, ncount, index, flux_up(:,l), buf2(:,l))
       end do
       call writefield ('swflux_up','upward sw flux','w/m2', buf2)
    end if

    call monscatter (ngrid, ncount, index, zflux,fsrfvis)
    do l = 1, nlayer
       call monscatter (ngrid, ncount, index, zdTsw(:,l), dTsw(:,l))
    end do
  end subroutine sw

  pure subroutine mongather (ngrid, n,index, a,b)
    integer,                   intent(in)  :: ngrid, n
    integer, dimension(n),     intent(in)  :: index
    real,    dimension(ngrid), intent(in)  :: a

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
    integer,                   intent(in)  :: ngrid, n
    integer, dimension(n),     intent(in)  :: index
    real, dimension(n),        intent(in)  :: b

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
