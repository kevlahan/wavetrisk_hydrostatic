module std_atm_profile_mod

!-------------------------------------------------------------------------------
!
! The barometric formula for U.S. Standard Atmosphere is valid up to 86 km.
! see https://en.wikipedia.org/wiki/Barometric_formula.
!
! N.B.  The extension above 86 km is using data from Hanli.  It is not complete
!       since the hardcoded parameter (c1) needs adjustment above 86 km.
!
!-------------------------------------------------------------------------------
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

real(8), parameter :: hb(nreg) = & ! height at bottom of layer (m)
     (/0d0, 1.1d4, 2.0d4, 3.2d4, 4.7d4, 5.1d4, 7.1d4, 8.6d4, &
     9.1d4, 1.1d5, 1.2d5, 1.5d5, 2d5, 3d5, 7d5/)

real(8), parameter :: pb(nreg) = & ! standard pressure (Pa)
     (/101325.d0, 22632.1d0, 5474.89d0, 868.02d0, 110.91d0, 66.94d0, 3.96d0, 3.7d-1,  &
     1.5d-1, 7.1d-3, 2.5d-3, 4.5d-4, 8.47d-5, 8.77d-6, 3.19d-8/)

real(8), parameter :: tb(nreg) = & ! standard temperature (K)
     (/288.15d0, 216.65d0, 216.65d0, 228.65d0, 270.65d0, 270.65d0, 214.65d0, 186.87d0,  &
     186.87d0, 240.d0, 360.d0, 634.39d0, 854.56d0, 976.01d0, 1d3/)

real(8), parameter :: lb(nreg) = & ! temperature lapse rate (K/m)
     (/-0.0065d0, 0.0d0, 0.001d0, 0.0028d0, 0.0d0, -0.0028d0, -0.001852d0, 0.0d0,       &
     2.796d-3, 0.012d0, 9.15d-3, 4.4d-3, 1.21d-3, 6d-5, 0.0d0/)

real(8), parameter :: rg = 8.3144598d0 ! universal gas constant (J/mol/K)
real(8), parameter :: g0 = 9.80665d0   ! gravitational acceleration (m/s^2)
real(8), parameter :: mw = 0.0289644d0 ! molar mass of dry air (kg/mol)
real(8), parameter :: c1 = g0*mw/rg
  
!=========================================================================================
CONTAINS
!=========================================================================================

subroutine std_atm_pres(height, pstd)
    
   ! arguments
   real(8), intent(in)  :: height(:) ! height above sea level in meters
   real(8), intent(out) :: pstd(:)   ! std pressure in Pa
    
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_pres'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0.0d0) then
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
      
      if (lb(ii) /= 0.d0) then
         pstd(k) = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(height(k) - hb(ii)) ) )**(c1/lb(ii))
      else
         pstd(k) = pb(ii) * exp( -c1*(height(k) - hb(ii))/tb(ii) )
      end if
      
   end do

 end subroutine std_atm_pres

 subroutine std_surf_pres (z_s, pstd)
   ! Finds surface pressure (i.e. pressure at a single height z_s)
   
   ! arguments
   real(8), intent(in)  :: z_s  ! height of surface above sea level in meters
   real(8), intent(out) :: pstd ! std pressure in Pa
    
   integer                     :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_surf_pres'
   !----------------------------------------------------------------------------

   if (z_s < 0d0) then ! extrapolate below mean sea level using troposphere lapse rate.
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

   if (lb(ii) /= 0d0) then
      pstd = pb(ii) * ( tb(ii) / (tb(ii) + lb(ii)*(z_s - hb(ii)) ) )**(c1/lb(ii))
   else
      pstd = pb(ii) * exp( -c1 * (z_s - hb(ii))/tb(ii) )
   end if

end subroutine std_surf_pres

!=========================================================================================

subroutine std_atm_height(pstd, height)
    
   ! arguments
   real(8), intent(in)   :: pstd(:)   ! std pressure in Pa
   real(8), intent(out)  :: height(:) ! height above sea level in meters
    
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

      if (lb(ii) /= 0.d0) then
         height(k) = hb(ii) + (tb(ii)/lb(ii)) * ( (pb(ii)/pstd(k))**(lb(ii)/c1) - 1d0)
      else
         height(k) = hb(ii) + (tb(ii)/c1)*log(pb(ii)/pstd(k))
      end if
   end do

end subroutine std_atm_height

!=========================================================================================

subroutine std_atm_temp(height, temp)
    
   ! arguments
   real(8), intent(in)   :: height(:) ! std pressure in Pa
   real(8), intent(out)  :: temp(:)   ! temperature
    
   ! local vars
   integer :: i, ii, k, nlev
   character(len=*), parameter :: routine = 'std_atm_temp'
   !----------------------------------------------------------------------------
    
   nlev = size(height)
   do k = 1, nlev
      if (height(k) < 0d0) then
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

      if (lb(ii) /= 0d0) then
         temp(k) = tb(ii) + lb(ii)*(height(k) - hb(ii))
      else
         temp(k) = tb(ii)
      end if
      
   end do

end subroutine std_atm_temp

end module std_atm_profile_mod
