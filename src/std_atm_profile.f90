module std_atm_profile_mod

  ! -------------------------------------------------------------------------------
  !
  ! The barometric formula for U.S. Standard Atmosphere is valid up to 86 km.
  ! see https://en.wikipedia.org/wiki/Barometric_formula.
  !
  ! N.B.  The extension above 86 km is using data from Hanli.  It is not complete
  !       since the hardcoded parameter (c1) needs adjustment above 86 km.
  !
  ! -------------------------------------------------------------------------------
  use kind_mod
  implicit none
  private
  save

  public :: &
       std_atm_pres,   & ! compute pressure at a vector of given heights
       std_surf_pres,  & ! compute pressure at a single height
       std_atm_height, & ! compute height given pressure
       std_atm_temp      ! compute temperature given height

  ! Parameters for barometric formula for U.S. Standard Atmosphere.

  integer, parameter  :: nreg = 15  ! number of regions

  real(dp), parameter :: hb(nreg) = & ! height at bottom of layer (m)
       (/ 0.0_dp, 1.1e4_dp, 2.0e4_dp, 3.2e4_dp, 4.7e4_dp, 5.1e4_dp, 7.1e4_dp, 8.6e4_dp, &
       9.1e4_dp, 1.1e5_dp, 1.2e5_dp, 1.5e5_dp, 2e5_dp, 3e5_dp, 7e5_dp /)

  real(dp), parameter :: pb(nreg) = & ! standard pressure (Pa)
       (/ 101325.0_dp, 22632.1_dp, 5474.89_dp, 868.02_dp, 110.91_dp, 66.94_dp, 3.96_dp, 3.7e-1_dp,  &
       1.5e-1_dp, 7.1e-3_dp, 2.5e-3_dp, 4.5e-4_dp, 8.47e-5_dp, 8.77e-6_dp, 3.19e-8_dp /)

  real(dp), parameter :: tb(nreg) = & ! standard temperature (K)
       (/288.15_dp, 216.65_dp, 216.65_dp, 228.65_dp, 270.65_dp, 270.65_dp, 214.65_dp, 186.87_dp,  &
       186.87_dp, 240._dp, 360._dp, 634.39_dp, 854.56_dp, 976.01_dp, 1000.0_dp/)

  real(dp), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
       (/-0.0065_dp, 0.0_dp, 0.001_dp, 0.0028_dp, 0.0_dp, -0.0028_dp, -0.001852_dp, 0.0_dp,       &
       2.796e-3_dp, 0.012_dp, 9.15e-3_dp, 4.4e-3_dp, 1.21e-3_dp, 6e-5_dp, 0.0_dp/)

  real(dp), parameter :: rg = 8.3144598_dp ! universal gas constant (J/mol/K)
  real(dp), parameter :: g0 = 9.80665_dp   ! gravitational acceleration (m/s^2)
  real(dp), parameter :: mw = 0.0289644_dp ! molar mass of dry air (kg/mol)
  real(dp), parameter :: c1 = g0*mw/rg

  !=========================================================================================
CONTAINS
  !=========================================================================================

  subroutine std_atm_pres(height, pstd)

    ! arguments
    real(dp), intent(in)  :: height(:) ! height above sea level in meters
    real(dp), intent(out) :: pstd(:)   ! std pressure in Pa

    integer :: i, ii, k, nlev
    character(len=*), parameter :: routine = 'std_atm_pres'
    !----------------------------------------------------------------------------

    nlev = size(height)
    do k = 1, nlev
       if (height(k) < 0.0_dp) then
          ! Extrapolate below mean sea level using troposphere lapse rate.
          ii = 1
       else
          ! find region containing height
          find_region: do i = nreg, 1, -1
             if (height(k) >= hb(i)) then
                ii = i
                exit find_region
             end if
          end do find_region
       end if

       if (lb(ii) /= 0.0_dp) then
          pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
       else
          pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
       end if

    end do

  end subroutine std_atm_pres

  subroutine std_surf_pres (z_s, pstd)
    ! Finds surface pressure (i.e. pressure at a single height z_s)

    ! arguments
    real(dp), intent(in)  :: z_s  ! height of surface above sea level in meters
    real(dp), intent(out) :: pstd ! std pressure in Pa

    integer                     :: i, ii, k, nlev
    character(len=*), parameter :: routine = 'std_surf_pres'
    !----------------------------------------------------------------------------

    if (z_s < 0.0_dp) then ! extrapolate below mean sea level using troposphere lapse rate.
       ii = 1
    else
       ! find region containing z_s
       find_region: do i = nreg, 1, -1
          if (z_s >= hb(i)) then
             ii = i
             exit find_region
          end if
       end do find_region
    end if

    if (lb(ii) /= 0.0_dp) then
       pstd = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(z_s - hb(ii)) ) )**(c1/lb(ii))
    else
       pstd = pb(ii) * exp( -c1 * (z_s - hb(ii))/tb(ii) )
    end if

  end subroutine std_surf_pres

  !=========================================================================================

  subroutine std_atm_height(pstd, height)

    ! arguments
    real(dp), intent(in)   :: pstd(:)   ! std pressure in Pa
    real(dp), intent(out)  :: height(:) ! height above sea level in meters

    integer :: i, ii, k, nlev
    logical :: found_region
    character(len=*), parameter :: routine = 'std_atm_height'
    !----------------------------------------------------------------------------

    nlev = size(height)
    do k = 1, nlev

       if (pstd(k) <= pb(nreg)) then
          ii = nreg
       else if (pstd(k) > pb(1)) then
          ii = 1
       else
          ! find region containing pressure
          find_region: do i = 2, nreg
             if (pstd(k) > pb(i)) then
                ii = i - 1
                exit find_region
             end if
          end do find_region
       end if

       if (lb(ii) /= 0.0_dp) then
          height(k) = hb(ii) + (tb(ii)/lb(ii)) * ( (pb(ii)/pstd(k))**(lb(ii)/c1) - 1.0_dp)
       else
          height(k) = hb(ii) + (tb(ii)/c1)*log(pb(ii)/pstd(k))
       end if
    end do

  end subroutine std_atm_height

  !=========================================================================================

  subroutine std_atm_temp(height, temp)

    ! arguments
    real(dp), intent(in)  :: height(:) ! std pressure in Pa
    real(dp), intent(out) :: temp(:)   ! temperature

    ! local vars
    integer :: i, ii, k, nlev
    character(len=*), parameter :: routine = 'std_atm_temp'
    !----------------------------------------------------------------------------

    nlev = size(height)
    do k = 1, nlev
       if (height(k) < 0.0_dp) then
          ii = 1
       else
          ! find region containing height
          find_region: do i = nreg, 1, -1
             if (height(k) >= hb(i)) then
                ii = i
                exit find_region
             end if
          end do find_region
       end if

       if (lb(ii) /= 0.0_dp) then
          temp(k) = tb(ii) + lb(ii)*(height(k) - hb(ii))
       else
          temp(k) = tb(ii)
       end if

    end do

  end subroutine std_atm_temp

end module std_atm_profile_mod
