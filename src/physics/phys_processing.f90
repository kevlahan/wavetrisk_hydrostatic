!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: phys_processing.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Module contianing all subroutines to save the mean_values of the prognostic variables
!               during a write to a file and routines to save grid coordinates to a file.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module phys_processing_mod
  ! Use cases
  use ops_mod
  use physics_call_mod
  implicit none
contains
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
    integer                       :: iwrt
    integer                       :: k
    
    real(8)                       :: area, high_lat_area, mid_lat_area, low_lat_area, low_lat_upper_lim, high_lat_lower_lim
    real(8), dimension(1:zlevels) :: z, dz, T_avg, Total_Tavg, Pressure_avg, Geopot_avg
    real(8), dimension(1:zlevels) :: Zonal_vel_avg, Meridional_avg, Zonal_KE_avg, Merid_KE_avg
    real(8), dimension(1:zlevels) :: Low_Lat_avg, Mid_lat_avg, High_lat_avg
    real(8), dimension(2)         :: mid_lat_range
    real(8), dimension(0:Nsoil)   :: tsurf_soil_avg
    
    character(4)                  :: s_time

    ! Set the ranges
    low_lat_upper_lim  = 23.5d0 * MATH_PI/180d0                       ! upper limit of low latitude range
    high_lat_lower_lim = 66.5d0 * MATH_PI/180d0                       ! lower limit of high latitude range
    mid_lat_range      = (/ low_lat_upper_lim , high_lat_lower_lim /) ! mid latitude range

    ! Set latitude averages and area to zero
    high_lat_area  = 0
    mid_lat_area   = 0
    low_lat_area   = 0

    Low_Lat_avg    = 0
    Mid_lat_avg    = 0
    High_lat_avg   = 0
    tsurf_soil_avg = 0

    ! Calculate surface pressure, as pressure is calculated in temp_fun to get temperature
    call cal_surf_press (sol)

!!! change to integrate_tri for adaptive grid !!!
    
    ! Calculate sums of the area over sphere and over low, mid and high regions (high/mid/low_lat_area updated in area)
    area = integrate_hex (area_fun, 1, level_start)
    
    high_lat_area = sum_real (high_lat_area)   
    mid_lat_area  = sum_real (mid_lat_area)   
    low_lat_area  = sum_real (low_lat_area)   

    ! Sums of temp, geopotential, and velocities, Kinetic Energies, Temp of Latitude zones
    do k = 1, zlevels
       T_avg(k)          = integrate_hex (temp_fun, k, level_start) ! this will take care if mpi is used

       Pressure_avg(k)   = integrate_hex (pressure_fun, k, level_start)
       Geopot_avg(k)     = integrate_hex (geopot_fun, k, level_start)
       
       Zonal_vel_avg(k)  = integrate_hex (zonal_fun, k, level_start)
       Meridional_avg(k) = integrate_hex (merid_fun, k, level_start)
       
       Zonal_KE_avg(k)   = integrate_hex (zonal_KE, k, level_start)
       Merid_KE_avg(k)   = integrate_hex (merid_KE, k, level_start)
       
       ! Temps of zones on each rank summed in temp_fun call for T_avg
       High_lat_avg(k) = sum_real (High_lat_avg(k)) 
       Mid_lat_avg(k)  = sum_real (Mid_lat_avg(k)) 
       Low_Lat_avg(k)  = sum_real (Low_Lat_avg(k)) 
    end do
    
    ! Calculate sums of the surface temp and soil temps
    do k = 0, zmin, -1
       tsurf_soil_avg(abs(k)) = integrate_hex (surf_soil_temp_fun, k, level_start)
    end do

    if (rank == 0) then
       ! Averages
       Total_Tavg     = T_avg        / area
       Pressure_avg   = Pressure_avg / area
       Geopot_avg     = Geopot_avg   / area
       
       Zonal_vel_avg  = Zonal_vel_avg / area
       Meridional_avg = Meridional_avg / area
       
       Zonal_KE_avg   = 0.5d0 * Zonal_KE_avg / area
       Merid_KE_avg   = 0.5d0 * Merid_KE_avg / area
       
       Low_Lat_avg    = Low_Lat_avg  / low_lat_area
       Mid_lat_avg    = Mid_lat_avg  / mid_lat_area
       High_lat_avg   = High_lat_avg / high_lat_area
       
       tsurf_soil_avg = tsurf_soil_avg / area

       ! Write values to terminal and file
       write (s_time, '(i4.4)') iwrt
       open (unit=20, file=trim(run_id)//'.6.'//s_time, form="FORMATTED", action='WRITE', status='REPLACE')

       write (6,'(a, f4.1, a)') "Temperature profile at time ", time/HOUR, " h"
       write (6,'(a)') 'Level    z_k     pres_k     geopot_k     T(z_k)       u(z_k)       v(z_k)        1/2u(z_k)^2 '// &
            '  1/2v(z_k)^2  Low_lat_T(z_k) Mid_lat_T(z_k) High_lat_T(z_k)'

       write(20,*) time ! write the time
       do k = 1, zlevels
          z(k) = log (Pressure_avg(k)/p_0)*(-R_d*Total_Tavg(k)/grav_accel)
          
          write (6,'(i3, 3x, f9.2, 1x, f9.2, 1x, f9.2, 9(1x, es13.4))') &
               k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), Zonal_vel_avg(k), Meridional_avg(k),&
               Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k), High_lat_avg(k)
          
          ! Write unformatted to file
          write (20,*) k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
               Zonal_vel_avg(k), Meridional_avg(k), Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k),&
               High_lat_avg(k)
          
          ! Write formatted to file
          ! write (20,'(i3, 1x, 11(es13.6,1x))') k, z(k), Pressure_avg(k), Geopot_avg(k), Total_Tavg(k), &
          !    Zonal_vel_avg(k), Meridional_avg(k), Zonal_KE_avg(k), Merid_KE_avg(k), Low_Lat_avg(k), Mid_lat_avg(k),&
          !    High_lat_avg(k)
       end do
       
       ! Write  surface temperature and soil temp (if soil included)
       do k = 0, zmin,-1
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
      real(8) :: rho_dz, rho_dz_theta, potential_temp, temperature, lat, lon

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1
      
      call cart2sph (dom%node%elts(id_i), lon, lat)
      lat = abs (lat)

      mass      =>      sol(S_MASS,k)%data(d)%elts
      temp      =>      sol(S_TEMP,k)%data(d)%elts
      mean_m    => sol_mean(S_MASS,k)%data(d)%elts
      mean_t    => sol_mean(S_TEMP,k)%data(d)%elts
      exner     =>       exner_fun(k)%data(d)%elts

      ! Pressure at layer centers and interfaces of the column and geopotential at next interface 
      call integrate_pressure_up (dom, i, j, zlev, offs, dims)

      rho_dz       = mean_m(id_i) + mass(id_i)
      rho_dz_theta = mean_t(id_i) + temp(id_i)
      potential_temp = rho_dz_theta / rho_dz

      ! Temperature at the center of zlayer
      temperature = theta2temp (potential_temp, dom%press%elts(id_i))

      ! Gather the zonal latitude temps
      if (lat >= mid_lat_range(1) .and. lat <= mid_lat_range(2)) then
         Mid_lat_avg(zlev) = Mid_lat_avg(zlev) + temperature / dom%areas%elts(id_i)%hex_inv
      end if
      
      if (lat >= high_lat_lower_lim) High_lat_avg(zlev) = High_lat_avg(zlev) + temperature / dom%areas%elts(id_i)%hex_inv
      if (lat <= low_lat_upper_lim ) Low_Lat_avg(zlev)  = Low_lat_avg(zlev)  + temperature / dom%areas%elts(id_i)%hex_inv

      ! Save density in Ke dom arg for KE
      dom%ke%elts(id_i) = dom%press%elts(id_i) / (temperature * R_d)

      ! Return temperature
      temp_fun = temperature

      nullify (mass, temp, mean_m, mean_t, exner)
    end function temp_fun

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

      geopot_fun = interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i))
    end function geopot_fun

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

    real(8) function zonal_fun (dom, i, j, zlev, offs, dims)
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

      velo  => sol(S_VELO,zlev)%data(d)%elts
      velo1 => grid(d)%u_zonal%elts
      velo2 => grid(d)%v_merid%elts
      
      call interp_UVW_latlon (dom, i, j, zlev, offs, dims)
      
      zonal_fun = velo1(id_i)
      
      nullify (velo, velo1, velo2)
    end function zonal_fun

    real(8) function zonal_KE (dom, i, j, zlev, offs, dims)
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

      zonal_KE = dom%ke%elts(id_i) * dom%u_zonal%elts(id_i)**2
    end function zonal_KE

    real(8) function merid_fun (dom, i, j, zlev, offs, dims)
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

    real(8) function merid_KE (dom, i, j, zlev, offs, dims)
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

      merid_KE = dom%ke%elts(id_i) * dom%v_merid%elts(id_i)**2
    end function merid_KE

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
      real(8) :: lat, lon

      id_i = idx (i, j, offs, dims) + 1
      
      call cart2sph (dom%node%elts(id_i), lon, lat)
      
      lat = abs (lat)

      if (lat >= mid_lat_range(1) .and. lat <= mid_lat_range(2)) then
         mid_lat_area = mid_lat_area + 1d0 / dom%areas%elts(id_i)%hex_inv
      end if
      
      if (lat >= high_lat_lower_lim) high_lat_area = high_lat_area + 1d0 / dom%areas%elts(id_i)%hex_inv
      if (lat <= low_lat_upper_lim ) low_lat_area  = low_lat_area +  1d0 / dom%areas%elts(id_i)%hex_inv

      area_fun = 1d0
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

      integer :: d, id_i

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1

      surf_soil_temp_fun = sol(S_TEMP,zlev)%data(d)%elts(id_i)
    end function surf_soil_temp_fun
  end subroutine mean_values

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Save grid coords !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine get_coordinates ()
    !-----------------------------------------------------------------------------------
    !
    !   Description: Retrieve all latitude and longitude coordinates and save in
    !                 a file named: run_id_coordinates
    !
    !-----------------------------------------------------------------------------------
    ! Open a file to write the latitude and longitude
    open (90, file=trim(run_id)//'_coordinates', form="FORMATTED", action='WRITE', status='REPLACE')

    ! Retrieve longitude and latitude elements of each domain and combine into single array
    call apply_bdry (get_lat_lon, z_null, 0, 0)

    close (90)
  end subroutine get_coordinates

  subroutine get_lat_lon (dom, i, j, zlev, offs, dims)
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
    real(8) :: lat, lon

    id_i = idx (i, j, offs, dims) + 1

    call cart2sph (dom%node%elts(id_i), lon, lat)

    write (90, '(2(es13.6,1x))') lat, lon
  end subroutine get_lat_lon
end module phys_processing_mod
