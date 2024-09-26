module turbulence
#include "use_logging.h"
  implicit none
  save
  private
  real, parameter :: karman = 0.4
  real            :: lmixmin = 100.0, emin_turb = 1e-8
  public          :: vdif, lmixmin, emin_turb
contains
  pure subroutine vdif_cd (ngrid, pz0, pg, pz, pu, pv, pts, ph, pcdv, pcdh)
    !=======================================================================
    !
    !   Subject: computation of the surface drag coefficient using the
    !   -------  approch developed by Louis for ECMWF.
    !
    !   Author: Frederic Hourdin  15 /10 /93
    !   -------
    !=======================================================================
    integer,                intent(in) :: ngrid ! number of atmospheric columns
    
    real,                   intent(in) :: pg ! gravity (m s -2)
    
    real, dimension(ngrid), intent(in) :: pz0   ! surface roughness length scale (m)
    real, dimension(ngrid), intent(in) :: pz    ! height of the first atmospheric layer (m)
    real, dimension(ngrid), intent(in) :: pu    ! zonal wind (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pv    ! meridional wind (m/s) in the layer
    real, dimension(ngrid), intent(in) :: pts   ! surface temperature (K)
    real, dimension(ngrid), intent(in) :: ph    ! potential temperature T*(p/ps)^kappa
    
    real, dimension(ngrid), intent(out) :: pcdv ! Cd for the wind
    real, dimension(ngrid), intent(out) :: pcdh ! Cd for potential temperature

    integer         :: ig
    real            :: zu2, z1, zri, zcd0, zz
    real, parameter :: b = 5.0, c = 5.0, d = 5.0, umin2 = 1e-12, c2b = 2.0 * b, c3bc = 3.0 * b * c, c3b = 3.0 * b
    
    !-----------------------------------------------------------------------
    !   couche de surface:
    !   ------------------

    !     do ig=1,ngrid
    !        zu2=pu(ig)*pu(ig)+pv(ig)*pv(ig)+umin2
    !        pcdv(ig)=pz0(ig)*(1.+sqrt(zu2))
    !        pcdh(ig)=pcdv(ig)
    !     enddo
    !     RETURN
!!!! WARNING, verifier la formule originale de Louis!
    do ig = 1, ngrid
       zu2 = pu(ig)**2 +  pv(ig)**2 + umin2
       zri = pg * pz(ig) * (ph(ig) - pts(ig)) / (ph(ig) * zu2)
       z1 = 1.0 + pz(ig) / pz0(ig)
       
       zcd0 = karman / log (z1)
       zcd0 = zcd0**2 * sqrt (zu2)
       if (zri < 0.0) then
          z1 = b * zri / (1.0 + c3bc * zcd0 * sqrt (- z1 * zri))
          pcdv(ig) = zcd0 * (1.0 - 2.0 * z1)
          pcdh(ig) = zcd0 * (1.0 - 3.0 * z1)
       else
          zz = sqrt (1.0 + d * zri)
          pcdv(ig) = zcd0 / (1.0 + c2b * zri / zz)
          pcdh(ig) = zcd0 / (1.0 + c3b * zri * zz)
       end if
    end do
  end subroutine vdif_cd

  pure subroutine vdif_k (ngrid, nlay, ptimestep, pg, pzlev, pzlay, pz0, pu, pv, ph, pkv, pkh)
    ! FIXME : pkh := pkv
    integer,                       intent(in)  :: ngrid, nlay
    real,                          intent(in)  :: ptimestep, pg 
    real, dimension(ngrid),        intent(in)  :: pz0
    real, dimension(ngrid,nlay),   intent(in)  :: ph, pu, pv, pzlay
    real, dimension(ngrid,nlay+1), intent(in)  :: pzlev
    real, dimension(ngrid,nlay+1), intent(out) :: pkh, pkv

    integer :: ig, il
    real    :: lmix, zdu, zdv, zri, zdvodz2, zdz, z1

    pkv(:,1)      = 0.0
    pkh(:,1)      = 0.0
    pkv(:,nlay+1) = 0.0
    pkh(:,nlay+1) = 0.0

    do il = 2 , nlay
       do ig = 1, ngrid
          z1 = pzlev(ig,il) + pz0(ig)
          lmix = karman * z1 / (1.0 + karman * z1 / lmixmin)
          !           lmix=lmixmin
          ! WARNING test lmix=lmixmin
          zdu = pu(ig,il)    -    pu(ig,il-1)
          zdv = pv(ig,il)    -    pv(ig,il-1)
          zdz = pzlay(ig,il) - pzlay(ig,il-1)
          
          zdvodz2 = (zdu**2 + zdv**2) / zdz**2
          if (zdvodz2 < 1.0e-5) then 
             pkv(ig,il) = lmix * sqrt (emin_turb)
          else
             zri = 2.0 * pg * (ph(ig,il) - ph(ig,il-1)) / (zdz * (ph(ig,il) + ph(ig,il-1)) * zdvodz2)
             pkv(ig,il)= lmix * sqrt (max (lmix * lmix*zdvodz2 * (1.0 - zri/0.4), emin_turb))
          end if
          pkh(ig,il) = pkv(ig,il)
       end do
    end do
  end subroutine vdif_k

  subroutine vdif (ngrid, nlay, ptime, ptimestep, pcapcal, pz0, pplay, pplev, pzlay, pzlev, pu, pv, ph, ptsrf, pemis, &
       pdufi, pdvfi, pdhfi, pfluxsrf, pdudif, pdvdif, pdhdif, pdtsrf, lwrite)
    use phys_const, only: g,r,rcp,cpp
    !=======================================================================
    !
    !   Vertical diffusion
    !   Implicit scheme
    !   We start by adding physical tendencies to variables x
    !   and then solve:
    !      x(t+1) =  x(t) + dt * (dx/dt)phys(t)  +  dt * (dx/dt)difv(t+1)
    !
    !=======================================================================
    integer,                       intent(in) :: ngrid,nlay
    real,                          intent(in) :: ptime,ptimestep
    real, dimension(ngrid),        intent(in) :: pcapcal, pemis, pfluxsrf, pz0, ptsrf
    real, dimension(ngrid,nlay),   intent(in) :: pplay, pzlay, ph, pu, pv, pdufi, pdvfi, pdhfi
    real, dimension(ngrid,nlay+1), intent(in) :: pplev, pzlev
    logical,                       intent(in) :: lwrite

    real, dimension(ngrid),        intent(out) :: pdtsrf
    real, dimension(ngrid,nlay),   intent(out) :: pdhdif, pdudif, pdvdif

    
    integer                       :: ilev, ig, ilay
    integer                       :: unit, ierr, it1, it2
    integer                       :: cluvdb, putdat, putvdim, setname, setvdim

    real                          :: z4st, zcst1
    real, dimension(ngrid)        :: z1, z2, zdplanck, zcdh, zcdv, ztsrf2, zu2
    real, dimension(ngrid,nlay)   :: za, zb, zb0, zc, zd, zh, zu, zv
    real, dimension(ngrid,nlay+1) :: zkv, zkh
    !
    !------------------------------------------------------------------------------------------------------
    !   Initializations:
    !------------------------------------------------------------------------------------------------------
    !   Computation of rho*dz and dt*rho/dz=dt*rho**2 g/dp with rho = p/RT  =p/ (R Theta) (p/ps)**kappa
    !------------------------------------------------------------------------------------------------------
    za(:,1:nlay) = (pplev(:,1:nlay) - pplev(:,2:nlay+1)) / g

    zcst1 = 4.0 * g * ptimestep / r**2
    do ilev = 2, nlay
       zb0(:,ilev) = pplev(:,ilev) * (pplev(:,1) / pplev(:,ilev))**rcp / (ph(:,ilev-1) + ph(:,ilev))
       zb0(:,ilev) = zcst1 * zb0(:,ilev) * zb0(:,ilev) / (pplay(:,ilev-1) - pplay(:,ilev))
    end do
    zb0(:,1) = ptimestep * pplev(:,1) / (r*ptsrf)

    if (lwrite) then
       ig = ngrid/2 + 1
       WRITELOG(*,*) 'Pression (mbar) altitude (km),u,v,theta, rho dz'
       do ilay = 1, nlay
          WRITELOG(*,*) 0.01 * pplay(ig,ilay), 0.001 * pzlay(ig,ilay), pu(ig,ilay), pv(ig,ilay), ph(ig,ilay), za(ig,ilay)
       end do
       WRITELOG(*,*) 'Pression (mbar) altitude (km),zb'
       do ilev = 1, nlay
          WRITELOG(*,*) 0.01 * pplev(ig,ilev), 0.001 * pzlev(ig,ilev), zb0(ig,ilev)
       end do
       LOG_DBG('vdif')
    endif 

    !-----------------------------------------------------------------------
    !   2. Add physical tendencies
    !-----------------------------------------------------------------------
    zu = pu + pdufi * ptimestep
    zv = pv + pdvfi * ptimestep
    zh = ph + pdhfi * ptimestep

    !-----------------------------------------------------------------------
    !   3. Compute  drag coefficients Cd
    !-----------------------------------------------------------------------
    ! Compute surface drag coefficents zcdv() and zcdh()
    CALL vdif_cd (ngrid, pz0, g, pzlay, pu, pv, ptsrf, ph, zcdv, zcdh)

    ! Compute zkv() and zkh()
    CALL vdif_k (ngrid, nlay, ptimestep, g, pzlev, pzlay, pz0, pu, pv, ph, zkv, zkh)

    zu2 = pu(:,1)**2 + pv(:,1)**2
    zcdv = zcdv * sqrt (zu2)
    zcdh = zcdh * sqrt (zu2)

    if (lwrite) then
       WRITELOG(*,*) 'Diagnostique diffusion verticale'
       WRITELOG(*,*) 'LMIXMIN',lmixmin
       
       WRITELOG(*,*) 'Coefficients Cd pour v et h'
       WRITELOG(*,*) zcdv(ngrid/2+1), zcdh(ngrid/2+1)
       
       WRITELOG(*,*) 'Coefficients K pour v et h'
       do ilev = 1, nlay
          WRITELOG(*,*) zkv(ngrid/2+1,ilev), zkh(ngrid/2+1,ilev)
       end do
       LOG_DBG ('vdif')
    endif 

    !-----------------------------------------------------------------------
    !   Vertical integration for u
    !-----------------------------------------------------------------------
    zb(:,2:nlay) = zkv(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1) = zcdv * zb0(:,1)

    z1 = 1.0/ (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zu(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zu(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    zu(:,1) = zc(:,1)
    zu(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zu(:,1:nlay-1)

    !-----------------------------------------------------------------------
    !   Vertical integration for v
    !-----------------------------------------------------------------------
    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zv(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zv(:,ilay) + zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    zv(:,1)      = zc(:,1)
    zv(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zv(:,1:nlay-1)

    !------------------------------------------------------------------------
    !   Vertical integration for h
    !------------------------------------------------------------------------
    zb(:,2:nlay) = zkh(:,2:nlay) * zb0(:,2:nlay)
    zb(:,1)      = zcdh * zb0(:,1)

    z1 = 1.0 / (za(:,nlay) + zb(:,nlay))
    zc(:,nlay) = za(:,nlay) * zh(:,nlay) * z1
    zd(:,nlay) = zb(:,nlay) * z1

    do ilay = nlay-1, 1, -1
       z1 = 1.0 / (za(:,ilay) + zb(:,ilay) + zb(:,ilay+1) * (1.0 - zd(:,ilay+1)))
       zc(:,ilay) = (za(:,ilay) * zh(:,ilay)+ zb(:,ilay+1) * zc(:,ilay+1)) * z1
       zd(:,ilay) = zb(:,ilay) * z1
    end do

    !------------------------------------------------------------------------
    !   Add Planck contribution to the implicit scheme
    !------------------------------------------------------------------------
    z4st = 4.0 * 5.67e-8 * ptimestep
    zdplanck = z4st * pemis * ptsrf**3

    !-----------------------------------------------------------------------
    !   Compute the evolution of the surface temperature
    !-----------------------------------------------------------------------
    z1 = pcapcal * ptsrf + cpp * zb(:,1) * zc(:,1)         + zdplanck * ptsrf + pfluxsrf * ptimestep
    z2 = pcapcal         + cpp * zb(:,1) * (1.0 - zd(:,1)) + zdplanck

    ztsrf2 = z1 / z2
    zh(:,1) = zc(:,1) + zd(:,1) * ztsrf2
    pdtsrf = (ztsrf2 - ptsrf) / ptimestep

    !-----------------------------------------------------------------------
    !   Final vertical integration
    !-----------------------------------------------------------------------
    zh(:,2:nlay) = zc(:,2:nlay) + zd(:,2:nlay) * zh(:,1:nlay-1)

    !-----------------------------------------------------------------------
    !   Compute the vertical diffusion tendencies
    !-----------------------------------------------------------------------
    pdudif = (zu - (pu + pdufi * ptimestep)) / ptimestep
    pdvdif = (zv - (pv + pdvfi * ptimestep)) / ptimestep
    pdhdif = (zh - (ph + pdhfi * ptimestep)) / ptimestep

    if (lwrite) then
       WRITELOG(*,*)
       WRITELOG(*,*) 'Diagnostique de la diffusion verticale'
       WRITELOG(*,*) 'h avant et apres diffusion verticale'
       WRITELOG(*,*) ptsrf(ngrid/2+1), ztsrf2(ngrid/2+1)
       do  ilev = 1, nlay
          WRITELOG(*,*) ph(ngrid/2+1,ilev), zh(ngrid/2+1,ilev)
       end do
       LOG_DBG('vdif')
    end if 
  end subroutine vdif
end module turbulence
