!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_simple.f90
! Author: Gabrielle Ching-Johnson, Nicholas Kevlahan
!
! Date Revised: Sept 11 2024
!
! Contains all subroutines needed to compute Backwards Euler step using simple physics package
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physics_simple_mod
  use utils_mod
  use adapt_mod
  use ops_mod
  use init_physics_mod
  implicit none
  type(Float_Field), dimension(:), allocatable, target :: U, V
contains
  subroutine physics_simple_step
    ! Uses simple physics modules to take a Backwards Euler step for physics using time step dt set by dynamics
    implicit none
    integer :: d, k

    call update_bdry (sol, NONE)
    
    call cal_surf_press (sol(:,1:zlevels))
    
    ! Initialize zonal and meridional velocity to zero
    U = sol(S_MASS,1:zlevels); call zero_float (U); V = U

    call apply_bdry (physics_call, z_null, 0, 1)
    physics_firstcall_flag = .false. ! update flag to false, once 1st call for all columns finished

    U%bdry_uptodate = .false.
    V%bdry_uptodate = .false.
    call update_bdry (U, NONE)
    call update_bdry (V, NONE)

    ! Assign velocity tendencies at edges
    do k = 1, zlevels
       do d = 1, size(grid)
          grid(d)%u_zonal%elts = U(k)%data(d)%elts
          grid(d)%v_merid%elts = V(k)%data(d)%elts
       end do
       call apply_bdry (save_velocity, k, 0, 0)
    end do
    sol%bdry_uptodate = .false.

    call WT_after_step (sol, wav_coeff, level_start-1)
  end subroutine physics_simple_step

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
    type(Domain)                   :: dom                  
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i, k, mask

    real(8), dimension(1:zlevels) :: phys_U       ! zonal velocity
    real(8), dimension(1:zlevels) :: phys_V       ! meridional velocity
    real(8), dimension(1:zlevels) :: phys_T       ! temperature
    real(8), dimension(1:zlevels) :: phys_Phi     ! geopotential
    real(8), dimension(1:1)       :: phys_Phisurf ! surface geopotential
    
    real(8), dimension(1:zlevels) :: phys_Play    ! pressure at layer centres
    real(8), dimension(0:zlevels) :: phys_Plev    ! pressure at layer interfaces

    real(8), dimension(1:Nsoil+1) :: Tsoil        !  surface temp and soil temperatures


    real(8) :: rho_dz
    real(8) :: latitude, longitude                ! coordinates of the column
    real(8) :: nth_day, day_fraction              ! day in simulation, fraction of the day

    logical(kind=C_BOOL) :: lastcall_flag = .false.

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

    ! Surface geopotential for this column
    phys_Phisurf = grav_accel * topography%data(d)%elts(id_i)

    mask = grid(d)%mask_n%elts(id_i)
    
    call physics_call_single_col (1, zlevels, mask, &
         physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, dt, &
         phys_Plev, phys_Play, phys_Phi, phys_Phisurf, &
         phys_U, phys_V, phys_T, Tsoil) ! updated

    lverbose = .false.

    do k = 1, zlevels
       ! Save zonal and meridional velocities at nodes (interpolated to edges once entire domain finished)
       U(k)%data(d)%elts(id_i) = phys_U(k)
       V(k)%data(d)%elts(id_i) = phys_V(k)

       ! Assign new potential temperature to WAVETRISK data structure
       rho_dz = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)

       sol(S_TEMP,k)%data(d)%elts(id_i) = rho_dz * temp2theta (phys_T(k), phys_Play(k))
    end do

    ! Assign soil column data to WAVETRISK data structure
    do k = zmin, 0
       sol(S_TEMP,k)%data(d)%elts(id_i) = Tsoil(abs(k)+1)
    end do
  contains
    subroutine pack_physics_vars
      ! Gathers all prognostic variables for all levels of the column into physics 2D data structure from the dynamics hybrid data structure
      integer :: k
      real(8) :: rho_dz, rho_dz_theta, theta

      phys_Plev(0) = dom%surf_press%elts(id_i)
      
      do k = 1, zlevels
         mass   =>      sol(S_MASS,k)%data(d)%elts
         temp   =>      sol(S_TEMP,k)%data(d)%elts
         mean_m => sol_mean(S_MASS,k)%data(d)%elts
         mean_t => sol_mean(S_TEMP,k)%data(d)%elts
         
         exner  =>       exner_fun(k)%data(d)%elts

         ! Get pressure at layer centers and interfaces of the column and geopotential at next interface (set in dom%geopot)
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)

         phys_Plev(k) = dom%press_lower%elts(id_i)
         phys_Play(k) = dom%press%elts(id_i)
         phys_Phi(k)  = interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i))

         rho_dz       = mass(id_i) + mean_m(id_i)
         rho_dz_theta = temp(id_i) + mean_t(id_i)

         theta = rho_dz_theta / rho_dz

         ! Temperature
         phys_T(k) = theta2temp (theta, phys_Play(k))

         ! Convert the edge velocities to zonal and meridional velocities
         velo  => sol(S_VELO,k)%data(d)%elts
         velo1 =>       grid(d)%u_zonal%elts
         velo2 =>       grid(d)%v_merid%elts
         
         call interp_UVW_latlon (dom, i, j, k, offs, dims)
         phys_U(k) = velo1(id_i)
         phys_V(k) = velo2(id_i)

         nullify (mass, temp, mean_m, mean_t, exner, velo, velo1, velo2)
      end do

      ! Retrieve surface temperature and soil column temperatures from dynamics
      ! (Tsoil(1) is surface temperature)
      do k = zmin, 0
         Tsoil(abs(k)+1) = sol(S_TEMP,k)%data(d)%elts(id_i)
      end do
    end subroutine pack_physics_vars
  end subroutine physics_call

  subroutine save_velocity (dom, i, j, zlev, offs, dims)
    ! Interpolates new zonal and meridional velocities to UVW and assigns to WAVETRISK data structure
    type(Domain)                         :: dom 
    integer                              :: i, j, zlev
    integer, dimension(N_BDRY+1)         :: offs
    integer, dimension(2,N_BDRY+1)       :: dims

    integer                    :: d, id
    real(8), dimension(1:EDGE) :: uvw

    d  = dom%id + 1  
    id = idx (i, j, offs, dims)

    call interp_latlon_UVW (dom, i, j, zlev, offs, dims, uvw)
    
    sol(S_VELO,zlev)%data(d)%elts(id_edge(id)) = uvw
  end subroutine save_velocity
end module physics_simple_mod
