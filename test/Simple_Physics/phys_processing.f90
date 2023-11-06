!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: phys_processing.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Module contianing all subroutines to save the mean_values of the prognostic variables
!               during a write to a file and routines to save grid coordinates to a file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module phys_processing_mod
   ! Use cases
   use test_case_mod
   use physics_call_mod
   implicit none

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Post Processing Calculations - Mean Values !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine mean_values (iwrt)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Saves averages over the sphere (assumes non-adaptive grid) of key
      !                 variables to a file.
      !
      !   Files Saved to: "run_id".6."iwrt"
      !   Note: iwrt is only a 4 character variable, if it represents the number of days
      !          then need to update this if simulation time is gerater than 4 chars
      !          (ie greater than 27 years). (iwrt is called in Simple_Physics.f90 and
      !           represents number of days)
      !
      !   Variables Saved:
      !        - Temperature
      !        - Zonal and Meridional Velocity
      !        - Zonal and Meridional Kinetic Energy (included density)
      !        - Low (0-23.5 deg), Mid (23.5-66.5) and High (66.5-90) Latitude Temperatures
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      use io_mod
      implicit none
      integer :: iwrt

      integer                       :: k
      real(8)                       :: area, high_lat_area, mid_lat_area, low_lat_area, low_lat_upper_lim, high_lat_lower_lim
      real(8), dimension(1:zlevels) :: dz, T_avg, Total_Tavg, Pressure_avg, Geopot_avg, Zonal_vel_avg, Meridional_avg, &
         Zonal_KE_avg, Merid_KE_avg, Low_Lat_avg, Mid_lat_avg, High_lat_avg
      real(8), dimension(1:zlevels) ::  z
      real(8), dimension(2) :: mid_lat_range
      real(8), dimension(0:Nsoil) :: tsurf_soil_avg
      character(4)                  :: s_time

      ! Set the ranges
      low_lat_upper_lim = (23.5*MATH_PI/180) ! Upper limit of low latitude range
      mid_lat_range = (/(23.5*MATH_PI/180), (66.5*MATH_PI/180)/) ! Mid latitude range
      high_lat_lower_lim = (66.5*MATH_PI/180) ! Lower limit of high latitude range

      !Set latatide averages and area to zero
      Low_Lat_avg = 0
      Mid_lat_avg = 0
      High_lat_avg = 0
      high_lat_area = 0
      mid_lat_area = 0
      low_lat_area = 0
      tsurf_soil_avg = 0

      ! Calculate surface pressure, as pressure is calculated in temp_fun to get temperature
      call cal_surf_press_phys(sol(1:N_VARIABLE,1:zmax))

      ! Calucate the sums of the area over sphere and over low, mid and high regions (high/mid/low_lat_area updated in area)
      area = integrate_hex (area_fun, 1, .true.)
      high_lat_area = sum_real(high_lat_area) !Get area summed over all ranks
      mid_lat_area = sum_real(mid_lat_area)   !Get area summed over all ranks
      low_lat_area = sum_real(low_lat_area)   !Get area summed over all ranks

      !Calculate Sums of temp, geopotential, and velocities, Kinetic Energies, Temp of Latitude zones
      do k = 1, zlevels
         T_avg(k)  = integrate_hex (temp_fun, k, .true.) ! This will take care if mpi is used
         Pressure_avg(k) = integrate_hex(pressure_fun, k, .true.)
         Geopot_avg(k) = integrate_hex (geopot_fun, k, .true.)
         Zonal_vel_avg(k) = integrate_hex(zonal_fun, k, .true.)
         Meridional_avg(k) = integrate_hex(merid_fun, k, .true.)
         Zonal_KE_avg(k) = integrate_hex(zonal_KE, k, .true.)
         Merid_KE_avg(k) = integrate_hex(merid_KE, k, .true.)
         ! Temps of zones on each rank summed in temp_fun call for T_avg
         High_lat_avg(k) = sum_real(High_lat_avg(k)) ! Get temp of high lat zones from all ranks
         Mid_lat_avg(k) = sum_real(Mid_lat_avg(k)) ! Get temp of mid lat zones from all ranks
         Low_Lat_avg(k) = sum_real(Low_Lat_avg(k)) ! Get temp of low lat zones from all ranks
      end do
      ! Calculate sums of the surface temp and soil temps
      do k = 0,zmin,-1
         tsurf_soil_avg(abs(k)) = integrate_hex(surf_soil_temp_fun,k,.true.)
      end do

      if (rank == 0) then
         ! Get average
         Total_Tavg  = T_avg  / area
         Pressure_avg = Pressure_avg / area
         Geopot_avg = Geopot_avg / area
         Zonal_vel_avg = Zonal_vel_avg / area
         Meridional_avg = Meridional_avg / area
         Zonal_KE_avg = 0.5 * Zonal_KE_avg / area
         Merid_KE_avg = 0.5 * Merid_KE_avg / area
         Low_Lat_avg = Low_Lat_avg / low_lat_area
         Mid_lat_avg = Mid_lat_avg / mid_lat_area
         High_lat_avg = High_lat_avg / high_lat_area
         tsurf_soil_avg = tsurf_soil_avg / area

         ! Write values to terminal and file
         write (s_time, '(i4.4)') iwrt
         open (unit=20, file=trim(run_id)//'.6.'//s_time, form="FORMATTED", action='WRITE', status='REPLACE')

         write (6,'(a, f4.1, a)') "Temperature profile at time ", time/HOUR, " h"
         write (6,'(a)') 'Level    z_k     pres_k     geopot_k     T(z_k)       u(z_k)       v(z_k)        1/2u(z_k)^2 '// &
            '  1/2v(z_k)^2  Low_lat_T(z_k) Mid_lat_T(z_k) High_lat_T(z_k)'

         write(20,*) time ! write the time
         do k = 1, zlevels
            z(k) = log(Pressure_avg(k)/p_0)*(-R_d*Total_Tavg(k)/grav_accel)
            write (6,'(i3, 3x, f9.2, 1x, f9.2, 1x, f9.2, 9(1x, es13.4))') &
               k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), Zonal_vel_avg(k), Meridional_avg(k),&
               Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k), High_lat_avg(k)
            ! Write unformatated to file
            write (20,*) k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
               Zonal_vel_avg(k), Meridional_avg(k), Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k),&
               High_lat_avg(k)
            ! Write formatted to file
            ! write (20,'(i3, 1x, 11(es13.6,1x))') k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
            !    Zonal_vel_avg(k), Meridional_avg(k), Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k),&
            !    High_lat_avg(k)
         end do
         ! Write the surface temperature and soil temp (if soil included)
         do k = 0,zmin,-1
            write(20,*) k, tsurf_soil_avg(abs(k))
         end do

         close (20)
      end if
   contains
      real(8) function temp_fun (dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the temperature at center of zlev layer of an element
         !                 in a domain and the zonal latitude temperatures.
         !                 Requires integration of the pressure.
         !                 Saves the density at the element node in dom%ke%elts for
         !                 the Kintic Energies.
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i
         real(8) :: full_mass, full_theta, potential_temp, temperature, lat, long

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1
         call cart2sph(dom%node%elts(id_i), long, lat)
         lat = abs(lat)

         ! Calculate the pressure
         call cal_press_geopot_layer(dom, i, j, zlev, offs, dims)

         full_mass  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)
         full_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)
         potential_temp = full_theta/full_mass

         ! Calculate the temperature at the center of zlayer
         temperature = potential_temp * ((dom%press%elts(id_i)/dom%surf_press%elts(id_i))**kappa)

         ! Gather the zonal latitiude temps
         if (lat .ge. mid_lat_range(1) .and. lat .le. mid_lat_range(2)) then
            Mid_lat_avg(zlev) = Mid_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)
         end if
         if (lat .ge. high_lat_lower_lim) High_lat_avg(zlev) = High_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)
         if (lat .le. low_lat_upper_lim) Low_Lat_avg(zlev) = Low_lat_avg(zlev) + (temperature/dom%areas%elts(id_i)%hex_inv)

         ! Save density in Ke dom arg for KE
         dom%ke%elts(id_i) = dom%press%elts(id_i)/(temperature*R_d)

         !Return temp
         temp_fun = temperature

      end function temp_fun

      real(8) function zonal_fun(dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the velocity of an element at the center of a zlev layer
         !                 in a domain, but returns zonal velocity.
         !                 The meridional velocity is saved in the dom%v_merid%elts for when
         !                    meridional integration is called.
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1

         velo => sol(S_VELO,zlev)%data(d)%elts
         velo1 => grid(d)%u_zonal%elts
         velo2 => grid(d)%v_merid%elts
         call interp_UVW_latlon(dom, i, j, zlev, offs, dims)
         zonal_fun = velo1(id_i)
         nullify (velo, velo1, velo2)
      end function zonal_fun

      real(8) function zonal_KE(dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the Zonal Kinetic energy of an element at the center of
         !                 a zlev layer in a domain
         !
         !   Assumptions:
         !        - velocity was already calculated and stored in dom%u_zonal, from zonal vel call
         !        - density was already calculated and stored in dom%ke, from the temp call
         !
         !-----------------------------------------------------------------------------------

         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: d, id_i

         d = dom%id + 1
         id_i = idx (i, j, offs, dims) + 1

         zonal_KE = dom%ke%elts(id_i)*((dom%u_zonal%elts(id_i))**2)

      end function zonal_KE

      real(8) function merid_fun(dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the meridional velocity of an element at the center of
         !                 a zlev layer in a domain
         !
         !   Assumptions:
         !        - velocity was already calculated and stored in dom%v_merid, from zonal velocity call
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         merid_fun = dom%v_merid%elts(id_i)

      end function merid_fun

      real(8) function merid_KE(dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the meridional kinetic energy of an element at the center of
         !                 a zlev layer in a domain.
         !
         !   Assumptions:
         !        - velocity was already calculated and stored in dom%v_merid, from zonal velocity call
         !        - density was already calculated and stored in dom%ke, from the temp call
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         merid_KE = dom%ke%elts(id_i)*((dom%v_merid%elts(id_i))**2)

      end function merid_KE

      real(8) function pressure_fun (dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the pressure of an element at the center of
         !                 a zlev layer in a domain.
         !
         !   Assumptions:
         !        - pressure at the zlevel interfaces are already calculated and stored in dom%press
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         pressure_fun = dom%press%elts(id_i)

      end function pressure_fun

      real(8) function geopot_fun (dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Calculates the geopotential of an element at the center of
         !                 a zlev layer in a domain.
         !
         !   Assumptions:
         !        - geopot at the zlevel interfaces above and below already calculated
         !              and stored in dom%geopot and dom%geopot_lower.
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i

         id_i = idx (i, j, offs, dims) + 1

         geopot_fun = 0.5*(dom%geopot_lower%elts(id_i)+dom%geopot%elts(id_i))

      end function geopot_fun

      real(8) function area_fun (dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Defines mass for total mass integration and for zonal latitude regions.
         !
         !-----------------------------------------------------------------------------------

         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i
         real(8) :: lat, long

         id_i = idx (i, j, offs, dims) + 1
         call cart2sph(dom%node%elts(id_i), long, lat)
         lat = abs(lat)

         if (lat .ge. mid_lat_range(1) .and. lat .le. mid_lat_range(2)) then
            mid_lat_area = mid_lat_area + (1/dom%areas%elts(id_i)%hex_inv)
         end if
         if (lat .ge. high_lat_lower_lim) high_lat_area = high_lat_area + (1/dom%areas%elts(id_i)%hex_inv)
         if (lat .le. low_lat_upper_lim) low_lat_area = low_lat_area + (1/dom%areas%elts(id_i)%hex_inv)

         area_fun = 1.0_8

      end function area_fun

      real(8) function surf_soil_temp_fun(dom, i, j, zlev, offs, dims)
         !-----------------------------------------------------------------------------------
         !
         !   Description: Retrieve the temperature of the surface and/or soil layer of an element.
         !
         !   Notes: zleve will be 0 or a negative number indicating soil layer or surface.
         !
         !-----------------------------------------------------------------------------------
         implicit none
         type(Domain)                   :: dom
         integer                        :: i, j, zlev
         integer, dimension(N_BDRY+1)   :: offs
         integer, dimension(2,N_BDRY+1) :: dims

         integer :: id_i,d

         id_i = idx (i, j, offs, dims) + 1
         d = dom%id + 1

         surf_soil_temp_fun = sol(S_TEMP,zlev)%data(d)%elts(id_i)

      end function surf_soil_temp_fun

   end subroutine mean_values

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Save grid coords !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   subroutine get_coordinates()
      !-----------------------------------------------------------------------------------
      !
      !   Description: Retrieve all latitude and longitude coordinates and save in
      !                 a file named: run_id_coordinates
      !
      !-----------------------------------------------------------------------------------
      !Arguments
      integer   :: d, p
      real(8)  :: lat, long
      type(Coord) :: x_i

      ! Open a file to write the latitude and longitude
      open(90, file=trim(run_id)//'_coordinates', form="FORMATTED", action='WRITE', status='REPLACE')

      !Retrieve longitude and latitude elements of each domain and combine into single array
      do d = 1, size(grid)
         do p = 3, grid(d)%patch%length
            call apply_onescale_to_patch(get_lat_long, grid(d), p-1, z_null, 0, 0)
         end do
      end do

      close(90)
   end subroutine

   subroutine get_lat_long(dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Retrieve latitude and longitude of an element
      !
      !-----------------------------------------------------------------------------------
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims
      integer :: id_i
      real(8) :: lat, long

      id_i = idx (i, j, offs, dims) + 1

      call cart2sph(dom%node%elts(id_i), long, lat)

      write(90, '(2(es13.6,1x))') lat, long

   end subroutine get_lat_long


end module phys_processing_mod
