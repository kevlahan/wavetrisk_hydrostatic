!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_simple.f90
! Author: Gabrielle Ching-Johnson, Nicholas Kevlahan
!
! Date Revised: Sept 11 2024
!
! Contains all subroutines needed to compute trend using simple physics package
!
! Assumes  trend will be used in an Euler step only (uses sol, not input q)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physics_simple_mod
  use utils_mod
  use ops_mod
  use init_physics_mod
  implicit none
  type(Float_Field), dimension(:), allocatable, target :: dzonal, dmerid
contains
  subroutine trend_physics_simple (q, dq)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Trend used to call the physics for each column on each domain. Also
    !                updates velocity edge tendencies once all elements on a domain have
    !                been set in the placeholder trends.
    !
    !   Notes: This routine is input arguement to the physics step in the main program
    !           Simple_Physics.f90, (Euler subroutine).
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: q 
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: dq
    
    integer :: d, k

    call update_bdry (sol, NONE)
    
    call cal_surf_press (sol(:,1:zlevels))
    
    ! Initialize zonal and meridional velocity trends to zero
    dzonal = dq(S_MASS,1:zlevels); call zero_float (dzonal)
    dmerid = dzonal

    call apply_bdry (physics_call, z_null, 0, 1)
    physics_firstcall_flag = .false. ! update flag to false, once 1st call for all columns finished

    ! Assign velocity tendencies at edges
    do k = 1, zlevels
       do d = 1, size(grid)
          grid(d)%u_zonal%elts = dzonal(k)%data(d)%elts
          grid(d)%v_merid%elts = dmerid(k)%data(d)%elts
       end do
       call apply_bdry (save_velocity_tendencies, k, 0, 0)
    end do
    trend%bdry_uptodate = .false.

    sol(S_TEMP,zmin:0)%bdry_uptodate = .false.
    call update_bdry (sol(S_TEMP,zmin:0), NONE)
    
    q(S_TEMP,zmin:0) = sol(S_TEMP,zmin:0)
    
    dq = trend; dq%bdry_uptodate = .false.
  end subroutine trend_physics_simple

  subroutine physics_call (dom, i, j, zlev, offs, dims)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Subroutine used to update calculate the trend for a single element/column.
    !                Converts the dynamics progrnostic vars to needed physics structure, calls
    !                the physics package to get tendencies for the column and sets temp
    !                tendencies back to hybrid structure. With the physics package, the surface
    !                temp and soil temp (if on) will be saved in the pot temp hybrid structure,
    !                for when load balancing occurs.
    !
    !   Expectation: Expected that this subroutine is called/used by apply_one_scale_to_patch
    !                (or similar routines) routine of Wavetrisk to calculate trend for the column.
    !
    !   Called by: trend_physics subroutine
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    USE callkeys,          only : lverbose ! print physics model parameters
    USE single_column_mod, only : change_latitude_longitude, physics_call_single_col
    implicit none
    type(Domain)                     :: dom                  
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims

    integer :: d, id, id_i

    real(8), dimension(1:zlevels) :: &
         phys_temp,&             ! temperature
         phys_zonal,&            ! zonal velocity
         phys_meridional, &      ! meridional velocity
         phys_geopot,&           ! geopotential
         phys_pres_lay, &        ! layer pressures
         phys_dtemp,&            ! temperature tendencies
         phys_dzonal,&           ! zonal velocity tendencies
         phys_dmeridional, &     ! meridonal velcity tendencies
         phys_dsurfpres          ! surface pressure tendency

    real(8), dimension(1:Nsoil+1) :: surf_soil_temp !  surface temp and soil temperatures
    real(8), dimension(0:zlevels) :: phys_pres_lev  !  interfaces pressure

    real(8) :: latitude, longitude   ! coordinates of the 
    real(8) :: nth_day, day_fraction ! day in simulation, fraction of the day

    logical(kind=C_BOOL) :: lastcall_flag=.false.

    d = dom%id +1

    id = idx (i, j, offs, dims)
    id_i = id + 1

    day_fraction = (time - dt) / DAY
    nth_day      = floor (day_fraction)
    day_fraction = day_fraction - nth_day

    ! Save prognostic variables in physics data structure
    call pack_physics_vars

    ! Get latitude and longitude of the column
    call cart2sph (dom%node%elts(id_i), longitude, latitude)

    ! Update physics latitude and longitude
    call change_latitude_longitude (latitude, longitude)

    ! Call physics for current column
    phys_dtemp       = 0d0
    phys_dzonal      = 0d0
    phys_dmeridional = 0d0
    phys_dsurfpres   = 0d0
    
    call physics_call_single_col (1, zlevels, &
         physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, dt, &
         phys_pres_lev, phys_pres_lay, phys_geopot, phys_zonal,  phys_meridional, phys_temp, surf_soil_temp, &
         phys_dzonal, phys_dmeridional, phys_dtemp, phys_dsurfpres)

    lverbose = .false.
    
    ! Convert data received from physics back to hybrid structure
    call save_tendencies
  contains
    subroutine pack_physics_vars
      ! Gathers all prognostic variables for all levels of the column into physics 2D data structure from the dynamics hybrid data structure
      integer :: k
      real(8) :: rho_dz, rho_dz_theta, theta

      phys_pres_lev(0) = dom%surf_press%elts(id_i)
      
      do k = 1, zlevels
         mass   =>      sol(S_MASS,k)%data(d)%elts
         temp   =>      sol(S_TEMP,k)%data(d)%elts
         mean_m => sol_mean(S_MASS,k)%data(d)%elts
         mean_t => sol_mean(S_TEMP,k)%data(d)%elts
         exner  =>       exner_fun(k)%data(d)%elts

         ! Get pressure at layer centers and interfaces of the column and geopotential at next interface (set in dom%geopot)
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)

         phys_pres_lev(k) = dom%press_lower%elts(id_i)
         phys_pres_lay(k) = dom%press%elts(id_i)
         phys_geopot(k)   = interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i))

         rho_dz       = mass(id_i) + mean_m(id_i)
         rho_dz_theta = temp(id_i) + mean_t(id_i)

         theta = rho_dz_theta / rho_dz

         phys_temp(k) = theta2temp (theta, phys_pres_lay(k))

         ! Convert the edge velocities to zonal and meridional velocities
         velo  => sol(S_VELO,k)%data(d)%elts
         velo1 =>       grid(d)%u_zonal%elts
         velo2 =>       grid(d)%v_merid%elts
         
         call interp_UVW_latlon (dom, i, j, k, offs, dims)
         phys_zonal(k)      = velo1(id_i)
         phys_meridional(k) = velo2(id_i)

         nullify (mass, temp, mean_m, mean_t, exner, velo, velo1, velo2)
      end do

      ! Retrieve column surface and soil temperatures from hybrid structure
      do k = zmin, 0
         surf_soil_temp(abs(k)+1) = sol(S_TEMP,k)%data(d)%elts(id_i)
      end do
    end subroutine pack_physics_vars

    subroutine save_tendencies
      ! Saves tendencies calculated by the physics to dynamics hybrid data structure.
      integer :: k
      real(8) :: dtheta, rho_dz   

      do k = 1, zlevels
         ! Save zonal and meridional velocities at nodes (interpolated to edges once entire domain finished)
         dzonal(k)%data(d)%elts(id_i) = phys_dzonal(k)
         dmerid(k)%data(d)%elts(id_i) = phys_dmeridional(k)
         
         rho_dz = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
         
         dtheta = temp2theta (phys_dtemp(k), phys_pres_lay(k))

         trend(S_MASS,k)%data(d)%elts(id_i) = 0d0             ! no physics mass tendency due to hydrostatic approximation
         trend(S_TEMP,k)%data(d)%elts(id_i) = dtheta * rho_dz
      end do

      ! Save column surface soil temp in temperature hybrid structure
      do k = zmin, 0
         sol(S_TEMP,k)%data(d)%elts(id_i) = surf_soil_temp(abs(k)+1)
      end do
    end subroutine save_tendencies
  end subroutine physics_call

  subroutine save_velocity_tendencies (dom, i, j, zlev, offs, dims)
    !-----------------------------------------------------------------------------------
    !
    !   Description: Subroutine used to update the trend hybrid structure for a single element/column
    !                with the converted edge tendencies from the zonal and meridional tendencies.
    !
    !   Expectation: Expected that this subroutine is called/used by apply_one_scale_to_patch
    !                (or similar routines) routine of Wavetrisk to update the hybrid structure.
    !
    !   Assumptions: dom%u_zonal and dom%v_merid have been populated at all gird points for
    !                domain dom.
    !
    !   Author: Gabrielle Ching-Johnson
    !
    !-----------------------------------------------------------------------------------
    type(Domain)                         :: dom 
    integer                              :: i, j, zlev
    integer, dimension(N_BDRY+1)         :: offs
    integer, dimension(2,N_BDRY+1)       :: dims

    integer                    :: d, id
    real(8), dimension(1:EDGE) :: uvw

    d  = dom%id +1  
    id = idx (i, j, offs, dims)

    call interp_latlon_UVW (dom, i, j, zlev, offs, dims, uvw)

    trend(S_VELO,zlev)%data(d)%elts(id_edge(id)) = uvw
  end subroutine save_velocity_tendencies
end module physics_simple_mod
