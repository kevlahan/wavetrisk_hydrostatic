!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! File Name: physics_simple.f90
! Author: Gabrielle Ching-Johnson, Nicholas Kevlahan
!
! Date Revised: 2024-10-10
!
! Contains all subroutines needed to compute a Backwards Euler step using simple physics package
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module physics_simple_mod
  use utils_mod
  use adapt_mod
  use ops_mod
  use init_physics_mod
  implicit none
contains
  subroutine physics_simple_step
    ! Uses simple physics modules to take a Backwards Euler step for physics using time step dt set by dynamics
    implicit none
    integer :: d, k
    
    call update_bdry (sol, NONE)

    call cal_surf_press (sol(1:N_VARIABLE,1:zlevels))
    
    call apply_bdry (physics_call, z_null, 0, 1)

    physics_firstcall_flag = .false. ! update flag to false, once 1st call for all columns finished

    sol%bdry_uptodate = .false.
  end subroutine physics_simple_step

  subroutine physics_call (dom, i, j, zlev, offs, dims)
    !-----------------------------------------------------------------------------------
    !
    !   Backwards Euler physics step on a single element/column
    !
    !   Author:     Gabrielle Ching-Johnson
    !   Revised by: Nicholas Kevlahan 2024-10
    !
    !-----------------------------------------------------------------------------------
    use callkeys,           only : lverbose ! print physics model parameters
    use io_mod,             only : kinetic_energy
    use single_column_mod,  only : change_latitude_longitude, physics_call_single_col 
    implicit none
    type(Domain)                     :: dom                  
    integer                          :: i, j, zlev
    integer, dimension(N_BDRY+1)     :: offs
    integer, dimension(2,N_BDRY+1)   :: dims

    integer :: d, id, id_i, k, mask

    ! Variables needed for simple physics
    real(8), dimension(1:zlevels) :: phys_Play               ! pressure at layer centres
    real(8), dimension(0:zlevels) :: phys_Pint               ! pressure at layer interfaces
    real(8), dimension(1:1)       :: phys_Phisurf            ! surface geopotential

    real(8), dimension(1:zlevels) :: phys_Theta              ! potential temperature
    real(8), dimension(1:zlevels) :: phys_Phi                ! geopotential
    real(8), dimension(1:zlevels) :: phys_U, phys_V, phys_W  ! velocities at edges
    real(8), dimension(1:zlevels) :: phys_Umag               ! speed at nodes
    real(8), dimension(1:Nsoil+1) :: Tsoil                   ! surface temp and soil temperatures


    real(8), dimension(1:zlevels) :: rho_dz
    real(8)                       :: latitude, longitude     ! coordinates of the column
    real(8)                       :: nth_day, day_fraction   ! day in simulation, fraction of the day
    logical(kind=C_BOOL)          :: lastcall_flag = .false.

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
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

    ! Mask for adaptive grid
    mask = grid(d)%mask_n%elts(id_i)

    ! Backwards Euler step on current column
    call physics_call_single_col (1, zlevels, mask, &
         physics_firstcall_flag, lastcall_flag, nth_day, day_fraction, dt, &
         phys_Play, phys_Pint, phys_Phi, phys_Phisurf, phys_Umag, &
         phys_U, phys_V, phys_W, phys_Theta, Tsoil) ! updated

    ! Assign solution at t+dt
    do k = 1, zlevels
       sol(S_VELO,k)%data(d)%elts(id_edge(id)) = (/ phys_U(k), phys_V(k), phys_W(k) /)
       sol(S_TEMP,k)%data(d)%elts(id_i)        = rho_dz(k) * phys_Theta(k) - sol_mean(S_TEMP,k)%data(d)%elts(id_i)
    end do

    ! Assign soil column solution at t+dt to WAVETRISK data structure
    do k = zmin, 0
       sol(S_TEMP,k)%data(d)%elts(id_i) = Tsoil(abs(k)+1)
    end do
  contains
    subroutine pack_physics_vars
      ! Gathers all required variables from wavetrisk data structure for all layers of column into physics data structure 
      integer :: k
      real(8) :: rho_dz_theta

      phys_Pint(0)  = dom%surf_press%elts(id_i)

      do k = 1, zlevels
         mass   =>      sol(S_MASS,k)%data(d)%elts
         temp   =>      sol(S_TEMP,k)%data(d)%elts
         mean_m => sol_mean(S_MASS,k)%data(d)%elts
         mean_t => sol_mean(S_TEMP,k)%data(d)%elts
         exner  =>       exner_fun(k)%data(d)%elts

         rho_dz(k)    = mass(id_i) + mean_m(id_i)
         rho_dz_theta = temp(id_i) + mean_t(id_i)

         ! Pressure at layers and interface and geopotential at next interface (set in dom%geopot)
         call integrate_pressure_up (dom, i, j, zlev, offs, dims)

         ! Input variables for simple physics module
         phys_Pint(k)  = dom%press_lower%elts(id_i)
         phys_Play(k)  = dom%press%elts(id_i)
         phys_Phi(k)   = interp (dom%geopot_lower%elts(id_i), dom%geopot%elts(id_i))
         phys_Umag(k)  = sqrt (2.0 * kinetic_energy (dom, i, j, k, offs, dims))
         
         phys_U(k)     = sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1)
         phys_V(k)     = sol(S_VELO,k)%data(d)%elts(EDGE*id+DG+1)
         phys_W(k)     = sol(S_VELO,k)%data(d)%elts(EDGE*id+UP+1)
         phys_Theta(k) = rho_dz_theta / rho_dz(k)

         nullify (mass, temp, mean_m, mean_t, exner)
      end do

      ! Retrieve surface temperature and soil column temperatures from dynamics
      ! (Tsoil(1) is surface temperature)
      do k = zmin, 0
         Tsoil(abs(k)+1) = sol(S_TEMP,k)%data(d)%elts(id_i)
      end do
    end subroutine pack_physics_vars
  end subroutine physics_call
end module physics_simple_mod
