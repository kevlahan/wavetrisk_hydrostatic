!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_call_module.f90
! Author: Gabrielle Ching-Johnson
! Date Revised: Nov 5 2023
! Description: Module contianing all subroutines for a simple physics call of a single time step.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module physics_call_mod
   ! Use cases
   use init_physics_mod
   implicit none

   ! Variables
   type(Float_Field), dimension(:), pointer :: dzonal, dmerid ! pointers for the tendencies

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!! Routines Needed to Call the Physics !!!!!!!!!!!!!!!!!!!!
   subroutine trend_physics(q, dq)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Trend used to call the physics for each column on each domain. Also
      !                updates velocity edge tendencies once all elements on a domain have
      !                been set in the placeholder trends.
      !
      !   Notes: This routine is input arguement to the physics step in the main program
      !           Simple_Physics.f90, (euler subroutine).
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------

      implicit none
      type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: q !solutions
      type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: dq !tendencies
      integer :: p, d, dom_length, k
      type(Float_Field), dimension(1:zlevels), target :: trend_zonal, trend_merid

      ! Initialize zonal and meridional velocity placeholder trends
      trend_zonal = dq(S_MASS,1:zlevels)
      trend_merid = dq(S_MASS,1:zlevels)
      do k = 1, zlevels
         call zero_float_field (trend_zonal(k), S_MASS)
         call zero_float_field (trend_merid(k), S_MASS)
      end do
      dzonal => trend_zonal
      dmerid => trend_merid

      call update_array_bdry (sol, NONE, 27)

      !Get Surface pressure of all domains on rank
      call cal_surf_press_phys(q)

      ! Loop through each domain on rank to call physics
      do d = 1,size(grid)
         ! Loop through each patch on a domain
         do p = 3, grid(d)%patch%length
            call apply_onescale_to_patch(physics_call, grid(d), p-1, z_null, 0, 1)
         end do
      end do
      ! Update flags to false, once 1st call for all columns finished
      physics_firstcall_flag = .false.

      ! Update the edge tendencies of a domain (Can only occur once all velocities on a domain set)
      do k=1,zlevels
         do d = 1,size(grid)
            grid(d)%u_zonal%elts = dzonal(k)%data(d)%elts
            grid(d)%v_merid%elts = dmerid(k)%data(d)%elts
            !once all columns on domain has been updated, go patch by patch to change tendencies and update them
            do p = 3, grid(d)%patch%length
               call apply_onescale_to_patch(update_velocity_tendencies, grid(d), p-1, k, 0, 0)
            end do
         end do
      end do

      !nullify pointers
      nullify (dzonal, dmerid)

      ! Boundary update ! setting only the zleveld 1:zlevs float fields boundary to be updated
      dq%bdry_uptodate = .false.
   end subroutine trend_physics

   subroutine update_velocity_tendencies(dom, i, j, zlev, offs, dims)
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

      ! Arguments
      type(Domain)                         :: dom ! domain the column is on
      integer                              :: i, j, zlev
      integer, dimension(N_BDRY+1)         :: offs
      integer, dimension(2,N_BDRY+1)       :: dims
      ! Local variables
      integer :: id, id_i
      real(8), dimension (1:EDGE)     :: uvw

      ! Get id of element
      id = idx(i, j, offs, dims)
      id_i = id + 1

      !Convert zonal and meridional tendencies to edge tendencies (saved in uvw)
      call interp_latlon_UVW(dom, i, j, zlev, offs, dims, uvw)

      !Update the trend
      trend(S_VELO, zlev)%data(dom%id+1)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = uvw

   end subroutine update_velocity_tendencies

   subroutine physics_call(dom, i, j, zlev, offs, dims)
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

      ! Physics use cases
      USE callkeys,     ONLY : lverbose
      USE single_column_mod, ONLY : change_latitude_longitude, physics_call_single_col
      implicit none

      ! Arguments
      type(Domain)                     :: dom                  ! domain the column is on
      integer                          :: i, j, zlev
      integer, dimension(N_BDRY+1)     :: offs
      integer, dimension(2,N_BDRY+1)   :: dims

      ! Local variables for physics
      integer :: d,&   ! domain id
         id, &  ! id of the column in the domain -1
         id_i   ! id of column
      real(8), dimension(1:zlevels) :: phys_temp,&             ! Column temperature
         phys_zonal,&            ! Column zonal velocity
         phys_meridional, &      ! Column meridional velocity
         phys_geopot,&           ! Column geopotential
         phys_pres_lay, &        ! Column layer pressures
         phys_dtemp,&            ! Column temperature tendencies
         phys_dzonal,&           ! Column zonal velolcity tendencies
         phys_dmeridional, &     ! Column meridonal velcity tendencies
         phys_dsurfpres,&        ! Column surface pressure tendency
         phys_dtheta             ! Column potential temperature tendency
      real(8), dimension(1:Nsoil+1) :: surf_soil_temp          ! Column surface temp and soil temperatures
      real(8), dimension (0:zlevels) :: phys_pres_lev          ! Column interfaces pressure
      real(8) :: latitude, longitude                           ! Coordinates of the column
      real(8) :: nth_day, day_fraction                         ! Day in simulation, fraction of the day
      LOGICAL(KIND=C_BOOL) :: lastcall_flag=.false.

      ! Get domain id
      d = dom%id +1
      ! Get the id of the column in the domain
      id = idx (i, j, offs, dims)
      id_i = id + 1

      !Get the number of days and fraction of day
      day_fraction = (time-dt)/DAY
      nth_day = FLOOR(day_fraction)
      day_fraction = day_fraction-nth_day

      ! Set Surface pressure to what was initially calculated
      phys_pres_lev(0) = dom%press_lower%elts(id_i)

      ! Retrieve prognositc variables and save in physics data structure
      call retrieve_prog_vars

      ! Intialize tendencies to 0
      phys_dtemp = 0
      phys_dzonal = 0
      phys_dmeridional = 0
      phys_dsurfpres = 0

      ! Get latitude and longitude of the column
      call cart2sph(dom%node%elts(id_i), longitude, latitude)

      ! Update physics latitude and longitude
      call change_latitude_longitude(latitude, longitude)

      ! Call physics for current column
      call physics_call_single_col(1, zlevels, physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, dt, &
         phys_pres_lev, phys_pres_lay, phys_geopot, phys_zonal,  phys_meridional, phys_temp, surf_soil_temp, &
         phys_dzonal, phys_dmeridional, phys_dtemp, phys_dsurfpres)

      lverbose = .false.
      ! Convert data received from physics back to hybrid structure
      call save_tendencies

   contains
      subroutine retrieve_prog_vars
         ! Description: Gathers all prognostic variables for all levels of the column into physics 2D data structure
         !              from the dynamics hybrid data structure
         !Arguments
         integer :: k
         real(8) :: theta, full_mass_theta, full_mass

         do k= 1,zlevels
            ! Get pressure at layer centers and interfaces of the column and geopotential at next interface (set in dom%geopot)
            call cal_press_geopot_layer (dom, i, j, k, offs, dims)
            phys_pres_lev(k) = dom%press_lower%elts(id_i)
            phys_pres_lay(k) = dom%press%elts(id_i)

            ! Get geopotential at layer center
            phys_geopot(k) = 0.5*(dom%geopot_lower%elts(id_i)+dom%geopot%elts(id_i))

            ! Convert the dynamics temp (pot temp ie theta) to temperature  (for entire column)
            ! use T = theta(p/p0)^kappa
            full_mass = sol(S_MASS, k)%data(d)%elts(id_i) + sol_mean(S_MASS, k)%data(d)%elts(id_i)
            full_mass_theta = sol(S_TEMP, k)%data(d)%elts(id_i) + sol_mean(S_TEMP, k)%data(d)%elts(id_i)
            theta = full_mass_theta/full_mass
            phys_temp(k) = (((phys_pres_lay(k)/dom%surf_press%elts(id_i))**kappa)*theta)

            ! Convert the edge velocities to zonal and meridional velocities
            velo => sol(S_VELO,k)%data(d)%elts
            velo1 => grid(d)%u_zonal%elts
            velo2 => grid(d)%v_merid%elts
            call interp_UVW_latlon(dom, i, j, k, offs, dims)
            phys_zonal(k) = velo1(id_i)
            phys_meridional(k) = velo2(id_i)
            nullify (velo, velo1, velo2)
         end do

         ! Retrieve column surface and soil temperatures from hybrid structure
         do k = 0,zmin,-1
            surf_soil_temp(abs(k) + 1) = sol(S_TEMP, k)%data(d)%elts(id_i)
         end do
      end subroutine retrieve_prog_vars

      subroutine save_tendencies
         ! Description: Saves tendencies calculated by the physics to dynamics hybrid data structure.
         ! Arguments
         integer :: k
         real(8) :: full_mass

         do k= 1, zlevels
            ! updated mass trend - none due to hydrostatic balance and get full_mass
            trend(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
            full_mass = sol(S_MASS, k)%data(d)%elts(id_i) + sol_mean(S_MASS, k)%data(d)%elts(id_i)

            !Convert temperature tendency
            phys_dtheta(k) = phys_dtemp(k)*(dom%surf_press%elts(id_i)/phys_pres_lay(k))**kappa
            trend(S_TEMP,k)%data(d)%elts(id_i) = phys_dtheta(k)*full_mass

            !Save Zonal and Meridional Velocities in place holder, will be converted to edges once entire domain finished
            dzonal(k)%data(d)%elts(id_i) = phys_dzonal(k)
            dmerid(k)%data(d)%elts(id_i) = phys_dmeridional(k)
         end do

         ! Save column surface soil temp in temperature hybrid structure
         do k = 0,zmin,-1
            sol(S_TEMP,k)%data(d)%elts(id_i) = surf_soil_temp(abs(k) + 1)
         end do
      end subroutine save_tendencies
   end subroutine physics_call

   subroutine cal_surf_press_phys (q)
      implicit none
      !-----------------------------------------------------------------------------------
      !
      !   Description: Compute surface pressure of all domains and save in domain press_lower
      !                 element for upward integration. Set geopotential to surface
      !                 geopotential for upward integration.
      !
      !   Author: Gabrielle Ching-Johnson
      !
      !-----------------------------------------------------------------------------------
      type(Float_Field), dimension(1:N_VARIABLE,1:zmax), target :: q

      integer :: d, k, mass_type, p

      ! Set surface geopotential
      call apply (set_surf_geopot_phys, z_null)

      ! Calculate surface pressure
      do d = 1, size(grid)
         grid(d)%surf_press%elts = 0.0_8
         ! Get total mass of the column
         do k = 1, zlevels
            mass   => q(S_MASS,k)%data(d)%elts
            temp   => q(S_TEMP,k)%data(d)%elts
            mean_m => sol_mean(S_MASS,k)%data(d)%elts
            mean_t => sol_mean(S_TEMP,k)%data(d)%elts
            do p = 3, grid(d)%patch%length
               call apply_onescale_to_patch (column_mass_phys, grid(d), p-1, k, 0, 1)
            end do
            nullify (mass, mean_m, mean_t, temp)
         end do

         ! using hydrostatic approx get surface pressure (P-bottom - p_top = -row*g*dz) (recall mass = reference density*dz)
         grid(d)%surf_press%elts = grav_accel*grid(d)%surf_press%elts + p_top

         grid(d)%press_lower%elts = grid(d)%surf_press%elts
      end do
   end subroutine cal_surf_press_phys

   subroutine column_mass_phys (dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Sum up mass (ie add mass to previous sum) and save in surface
      !                 pressure array
      !
      !-----------------------------------------------------------------------------------
      implicit none
      type (Domain)                  :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i
      ! Get id of column
      id_i = idx (i, j, offs, dims) + 1

      dom%surf_press%elts(id_i) = dom%surf_press%elts(id_i) + mass(id_i)
   end subroutine column_mass_phys

   subroutine set_surf_geopot_phys (dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Set initial geopotential to surface geopotential
      !
      !-----------------------------------------------------------------------------------
      implicit none
      type (Domain)                  :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id

      d = dom%id + 1
      id = idx (i, j, offs, dims) + 1

      dom%geopot%elts(id) = surf_geopot (d, id)
   end subroutine set_surf_geopot_phys

   subroutine cal_press_geopot_layer (dom, i, j, zlev, offs, dims)
      !-----------------------------------------------------------------------------------
      !
      !   Description: Integrate pressure and geopotential up from surface to top layer.
      !                 Each time routine is called, integrate up one vertical layer.
      !
      !   Outputs/Saves:
      !        dom%press%elts   -> pressure at center of the layer
      !        dom%press_lower  -> pressure at the top interface to be used at next call
      !        dom%geopot_lower -> geopotential at bottom of interface
      !        dom%geopot       -> geopotential at top interface of layer
      !        * geopotentials will be used to calculate at center of layer
      !            in retrieve_prog_vars subroutine of physics_call
      !
      !-----------------------------------------------------------------------------------
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i, d
      real(8) :: full_mass, pressure_upper, pressure_lower, exner

      id_i = idx (i, j, offs, dims) + 1
      d = dom%id + 1

      ! Calculate mass
      full_mass = sol(S_MASS,zlev)%data(d)%elts(id_i) + sol_mean(S_MASS,zlev)%data(d)%elts(id_i)

      ! Get pressure at next zlevel interface
      pressure_lower = dom%press_lower%elts(id_i)
      pressure_upper = pressure_lower - grav_accel * full_mass
      ! Get pressure at the center of current zlevel
      dom%press%elts(id_i) = 0.5 * (pressure_lower + pressure_upper)

      ! Get geopotential at next zlevel interface
      dom%geopot_lower%elts(id_i) = dom%geopot%elts(id_i)
      exner = c_p * (dom%press%elts(id_i)/p_0)**kappa
      dom%geopot%elts(id_i) = dom%geopot_lower%elts(id_i) + &
         grav_accel*kappa*sol(S_TEMP,zlev)%data(d)%elts(id_i)*exner/dom%press%elts(id_i)

      ! Update pressure lower to next interface for next integration
      dom%press_lower%elts(id_i) = pressure_upper
   end subroutine cal_press_geopot_layer

end module physics_call_mod
