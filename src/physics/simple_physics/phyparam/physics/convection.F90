module convection
  implicit none
contains
  subroutine ConvAdj (ngrid, nlay, pTimestep, pPlay, pPlev, zExner, pU, pV, pTheta)
    !====================================================================================================================
    !
    !   Dry convective adjustment
    !
    !   Layer is statically unstable if potential temperatures pTheta(above) < pTheta(below).
    !
    !====================================================================================================================

    integer,                       intent(in)  :: ngrid     ! number of columns
    integer,                       intent(in)  :: nlay      ! number of vertical layers
    real,                          intent(in)  :: pTimestep ! timestep
    real, dimension(ngrid,nlay),   intent(in)  :: pPlay     ! pressure at layer centres
    real, dimension(ngrid,nlay+1), intent(in)  :: pPlev     ! pressure at interfaces
    real, dimension(ngrid,nlay),   intent(in)  :: zExner    ! Exner function at layer centres

    ! Output
    real, dimension(ngrid,nlay),   intent(inout)  :: pU     ! zonal velocity
    real, dimension(ngrid,nlay),   intent(inout)  :: pV     ! meridional velocity
    real, dimension(ngrid,nlay),   intent(inout)  :: pTheta ! potential temperature

    integer                        :: Icol, ig, jcnt, l, jj
    integer, dimension(ngrid)      :: jadrs
    real,    dimension(nlay)       :: sdSig, dsig
    real,    dimension(nlay+1)     :: Sig
    logical, dimension(ngrid)      :: vtest

    ! Label unstable layers
    vtest = .false.
    do l = 2, nlay
       do ig = 1, ngrid
          if (pTheta(ig,l) < pTheta(ig,l-1)) vtest(ig) = .true.
       end do
    end do

    ! Adjust unstable layers in each column
    do ig = 1, ngrid
       if (vtest(ig)) then
          Icol = ig
          call sigma_levels     (ngrid, nlay, Icol, pPlev, zExner, sig, dsig, sdsig) ! compute sigma coordinates of unstable column Icol
          call adjust_onecolumn (ngrid, nlay, Icol, sig, dsig, sdsig, pU, pV, pTheta)    ! convective adjustment of unstable column Icol
       endif
    end do
  end subroutine ConvAdj

  pure subroutine sigma_levels (ngrid, nlay, Icol, pPlev, zExner, Sig, dSig, sdSig)
    ! Compute sigma coordinates for column Icol

    ! Input
    integer,                       intent(in) :: ngrid, nlay, Icol
    real, dimension(ngrid,nlay),   intent(in) :: zExner
    real, dimension(ngrid,nlay+1), intent(in) :: pPlev

    ! Output
    real, dimension(nlay),        intent(out) :: dSig, sdSig
    real, dimension(nlay+1),      intent(out) :: Sig

    integer :: l

    Sig = pPlev(Icol,:) / pPlev(Icol,1)

    dSig(1:nlay) = Sig(1:nlay) - Sig(2:nlay+1)

    sdSig = zExner(Icol,:) * dSig
  end subroutine sigma_levels

  subroutine adjust_onecolumn (ngrid, nlay, Icol, Sig, dSig, sdSig, zU, zV, zH)
    ! Input
    integer,                     intent(in)    :: ngrid, nlay, Icol
    real, dimension(nlay),       intent(in)    :: dSig, sdSig
    real, dimension(nlay+1),     intent(in)    :: Sig

    ! Output
    real, dimension(ngrid,nlay), intent(inout) :: zU, zV, zH  ! adjusted velocities and potential temperature

    real    :: zHm, zSm, zUm, zVm, zAlpha
    integer :: l, l1, l2
    logical :: extend

    l2 = 1
    do while (l2 < nlay)
       l2 = l2 + 1
       ! Starting from l1 = l2, we extend l1..l2 downwards and upwards
       ! until zHm = potential temperature mass-weighted over l1..l2 is higher than potential temperature at level l1 - 1
       ! and lower than pot. temp. at level l2+1.
       !
       ! We incrementally compute zSm = mass of layers l1 .. l2 and zHm
       !
       zSm = sdSig(l2)
       zHm = zH(Icol,l2)
       l1 = l2
       extend = .true.
       do while (extend)
          ! Extend downwards if l1>1 and zHm is lower than layer l1-1
          extend = .false.
          if (l1 > 1) then
             if (zHm < zH(Icol,l1-1)) extend = .true.
          end if

          if (extend) then
             l1 = l1 - 1
             zSm = zSm + sdsig(l1)
             zHm = zHm + sdsig(l1) * (zH(Icol,l1) - zHm) / zSm
             cycle
          end if

          ! Extend upwards if l2 < nlay and zHm is higher than layer l2+1
          extend = .false.
          if (l2 < nlay) then
             if (zH(Icol,l2+1) < zHm) extend = .true.
          end if

          if (extend) then
             l2 = l2 + 1
             zSm = zSm + sdSig(l2)
             zHm = zHm + sdSig(l2) * (zH(Icol,l2) - zHm) / zSm
          end if
       end do

       ! At this point the profile l1-1 (l1 ...l2) l2+1 is stable

       if (l1 < l2) then
          ! Perform the mixing of layers l1...l2 :
          ! * potential temperature is fully mixed, ie replaced by mass-weighted average
          ! * momentum u,v is mixed with weight zAlpha, i.e. u: = zAlpha * mean(u) + (1 - zAlpha) * u
          !   where zAlpha = sum( abs(theta - mean(theta))*  dsig) / mean(theta) * sum(dsig)
          !   for large deviations of theta from its mean, this formula could in principe yield zAlpha>1.
          zAlpha = 0.0
          zUm    = 0.0
          zVm    = 0.0
          do l = l1, l2
             zAlpha = zAlpha + abs (zH(Icol,l) - zHm) * dSig(l)
             zH(Icol,l) = zHm

             ! Must use compute zUm, zVm from zU, zV to conserve momentum (see below)
             zUm = zUm + dSig(l) * zU(Icol,l)
             zVm = zVm + dSig(l) * zV(Icol,l)
          end do
          zAlpha = zAlpha / (zHm * (sig(l1)  - sig(l2+1)))
          zUm    = zUm    / (sig(l1) - sig(l2+1))
          zVm    = zVm    / (sig(l1) - sig(l2+1))
          if (zAlpha < 1e-5) zAlpha = 1e-5

          ! Since zUm, zVm are mass-weighted averages of zU, zV => mixing preserves total momentum
          zU(Icol,l1:l2) = zU(Icol,l1:l2) + zAlpha * (zUm - zU(Icol,l1:l2))
          zV(Icol,l1:l2) = zV(Icol,l1:l2) + zAlpha * (zVm - zV(Icol,l1:l2))
       end if
    end do
  end subroutine adjust_onecolumn
end module convection
