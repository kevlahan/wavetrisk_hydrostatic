module physics_Held_Suarez_mod
  ! Module providing trend routine for Held and Suarez (1994) physics.
  ! Bulletin of the American Meteorological Society 75 (10), 1825-1830
  use ops_mod
  use sso_mod
  implicit none

  ! Held-Suarez model parameters
  real(8) :: T_0            = 300d0      * KELVIN              ! reference temperature
  real(8) :: T_mean         = 315d0      * KELVIN              ! mean temperature
  real(8) :: T_tropo        = 200d0      * KELVIN              ! tropopause temperature
  real(8) :: u_0            = 70d0       * METRE/SECOND        ! maximum velocity of zonal wind
  real(8) :: k_a            = 1d0/40d0   / DAY                 ! cooling at free surface of atmosphere
  real(8) :: k_f            = 1d0        / DAY                 ! Rayleigh friction
  real(8) :: k_s            = 1d0/4d0    / DAY                 ! cooling at surface
  real(8) :: delta_T        = 65d0       * KELVIN/METRE        ! meridional temperature gradient
  real(8) :: delta_theta    = 10d0       * KELVIN/METRE        ! vertical temperature gradient
  real(8) :: sigma_b        = 0.7d0                            ! normalized tropopause pressure height
  real(8) :: gamma_T        = 5d-3       * KELVIN/METRE        ! temperature lapse rate
  real(8) :: delta_T2       = 4.8d5      * KELVIN              ! empirical temperature difference
  real(8) :: sigma_0        = 0.252d0                          ! value of sigma at reference level (level of the jet)
  real(8) :: sigma_t        = 0.2d0                            ! value of sigma at the tropopauses
contains
  subroutine trend_physics_held_suarez (q, dq)
    ! Trend for Held-Suarez physics
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target :: q, dq 

    integer :: d, k, p

    call update_bdry (q, NONE)

    ! Assign shared variable
    call zero_float (dq)
    
    ! Current surface pressure
    call cal_surf_press (q)

    do d = 1, size(grid)
       do k = 1, zlevels
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          
          mass   =>        q(S_MASS,k)%data(d)%elts
          temp   =>        q(S_TEMP,k)%data(d)%elts
          velo   =>        q(S_VELO,k)%data(d)%elts
          
          dmass  =>       dq(S_MASS,k)%data(d)%elts
          dtemp  =>       dq(S_TEMP,k)%data(d)%elts
          dvelo  =>       dq(S_VELO,k)%data(d)%elts
          
          exner  =>       exner_fun(k)%data(d)%elts
          
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (integrate_pressure_up, grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_scalars,         grid(d), p-1, k, 0, 1)
             call apply_onescale_to_patch (trend_velo,            grid(d), p-1, k, 0, 0)
          end do
          nullify (dmass, dtemp, exner, mass, mean_m, mean_t, temp, velo)

          ! Add SSO drag
          if (sso) then
             do p = 3, grid(d)%patch%length
                call apply_onescale_to_patch (trend_velo_sso, grid(d), p-1, z_null, 0, 0)
             end do
          end if
          nullify (dvelo)
       end do
    end do

    ! Set output
    dq%bdry_uptodate = .false.
  end subroutine trend_physics_held_suarez

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Trend for physics step (relaxation to equilibrium temperature)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer :: id
    real(8) :: k_T, lat, lon, sigma, theta_equil

    id = idx (i, j, offs, dims) + 1

    call cart2sph (dom%node%elts(id), lon, lat)

    call cal_theta_eq (dom%press%elts(id), dom%surf_press%elts(id), lat, theta_equil, k_T)

    dmass(id) = 0d0
    dtemp(id) = - k_T * temp(id)
  end subroutine trend_scalars

  subroutine trend_velo (dom, i, j, zlev, offs, dims)
    ! Velocity trend for physics step (Rayleigh friction)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                    :: id, id_i
    integer, dimension(1:EDGE) :: id_e
    real(8)                    :: k_v, sigma, sigma_c

    real(8), dimension(1:EDGE) :: drag

    id   = idx (i, j, offs, dims)
    id_i = id + 1
    id_e = id_edge (id)

    sigma = (dom%press%elts(id_i) - p_top) / (dom%surf_press%elts(id_i) - p_top)
    sigma_c = 1d0 - sigma_b
    
    k_v = k_f * max (0d0, (sigma - sigma_b) / sigma_c)
    
    dvelo(id_e) = - k_v * velo(id_e)
  end subroutine trend_velo

  subroutine trend_velo_sso (dom, i, j, zlev, offs, dims)
    ! Include SSO drag velocity trend 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                              :: d, id, k
    integer, dimension(1:EDGE)           :: id_e
    real(8), dimension(1:zlevels,1:EDGE) :: drag

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_e = id_edge (id)

    drag = sso_drag (dom, i, j, z_null, offs, dims)
    do k = 1, zlevels
      dvelo(id_e) = dvelo(id_e) + drag(k,:)
    end do
  end subroutine trend_velo_sso

  subroutine cal_theta_eq (p, p_s, lat, theta_equil, k_T)
    ! Returns equilibrium potential temperature theta_equil and Newton cooling constant k_T
    use domain_mod
    implicit none
    real(8) :: p, p_s, lat, theta_equil, k_T

    real(8) :: cs2, sigma, sigma_c, theta_force, theta_tropo

    cs2 = cos (lat)**2

    sigma = (p - p_top) / (p_s - p_top)
    sigma_c = 1d0 - sigma_b

    k_T = k_a + (k_s - k_a) * max (0d0, (sigma - sigma_b) / sigma_c) * cs2**2

    theta_tropo = T_tropo * (p / p_0)**(-kappa)  ! potential temperature at tropopause

    theta_force = T_mean - delta_T * (1d0 - cs2) - delta_theta * cs2 * log (p / p_0)

    theta_equil = max (theta_tropo, theta_force) ! equilibrium temperature
  end subroutine cal_theta_eq
end module physics_Held_Suarez_mod
