module physics_Held_Suarez_mod
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! 
  !    Trend using Held and Suarez physics model
  !
  !    Bulletin of the American Meteorological Society 75 (10), 1825-1830
  !
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  use ops_mod
  use sso_mod
  implicit none

  ! Held-Suarez model parameters
  real(dp) :: T_mean         = 315.0_dp       * KELVIN              ! mean temperature
  real(dp) :: T_tropo        = 200.0_dp       * KELVIN              ! tropopause temperature
  real(dp) :: k_a            = 1.0_dp/40.0_dp / DAY                 ! cooling at free surface of atmosphere
  real(dp) :: k_f            = 1.0_dp         / DAY                 ! Rayleigh friction
  real(dp) :: k_s            = 0.25_dp        / DAY                 ! cooling at surface
  real(dp) :: delta_T        = 60.0_dp        * KELVIN/METRE        ! meridional temperature gradient
  real(dp) :: delta_theta    = 10.0_dp        * KELVIN/METRE        ! vertical temperature gradient
  real(dp) :: sigma_b        = 0.7_dp                               ! normalized tropopause pressure height
  real(dp) :: gamma_T        = 5e-3_dp        * KELVIN/METRE        ! temperature lapse rate
  real(dp) :: delta_T2       = 4.8e5_dp       * KELVIN              ! empirical temperature difference
  real(dp) :: sigma_0        = 0.252_dp                             ! value of sigma at reference level (level of the jet)
  real(dp) :: sigma_t        = 0.2_dp                               ! value of sigma at the tropopauses
contains
  subroutine trend_physics_Held_Suarez (q, dq)
    ! Trend for Held-Suarez physics
    implicit none
    type(Float_Field), dimension(1:N_VARIABLE,zmin:zmax), target  :: q  ! includes soil layers
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels), target  :: dq

    integer :: d, k, p

    call update_bdry (q, NONE, 941)

    ! Assign shared variable
    call zero_float (dq)
    
    ! Current surface pressure
    call cal_surf_press (q(1:N_VARIABLE,1:zlevels))

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
    dq%bdry_uptodate = .false.
  end subroutine trend_physics_Held_Suarez

  subroutine trend_scalars (dom, i, j, zlev, offs, dims)
    ! Trend for physics step (relaxation to equilibrium temperature)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    
    integer  :: id
    real(dp) :: k_T, lat, lon, rho_dz, rho_dz_dtheta, sigma, theta_equil

    id = idx (i, j, offs, dims) + 1

    rho_dz        = mean_m(id) + mass(id)
    rho_dz_dtheta = mean_t(id) + temp(id)

    call cart2sph (dom%node%elts(id), lon, lat)
    call cal_theta_eq (dom%press%elts(id), dom%surf_press%elts(id), lat, theta_equil, k_T)

    dmass(id) = 0.0_dp
    dtemp(id) = - k_T * (rho_dz_dtheta - rho_dz * theta_equil)
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
    real(dp)                   :: k_v, sigma, sigma_c

    id   = idx (i, j, offs, dims)
    id_i = id + 1
    id_e = id_edge (id)

    sigma = (dom%press%elts(id_i) - p_top) / (dom%surf_press%elts(id_i) - p_top)
    sigma_c = 1.0_dp - sigma_b
    
    k_v = k_f * max (0.0_dp, (sigma - sigma_b) / sigma_c)
    
    dvelo(id_e) = - k_v * velo(id_e)
  end subroutine trend_velo

  subroutine trend_velo_sso (dom, i, j, zlev, offs, dims)
    ! Include SSO drag velocity trend 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                               :: d, id, k
    integer,  dimension(1:EDGE)           :: id_e
    real(dp), dimension(1:zlevels,1:EDGE) :: drag

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
    real(dp) :: p, p_s, lat, theta_equil, k_T

    real(dp) :: cs2, sn2, sigma, sigma_c, theta_force, theta_tropo

    cs2 = cos (lat)**2
    sn2 = sin (lat)**2

    sigma = (p - p_top) / (p_s - p_top)
    sigma_c = 1.0_dp - sigma_b

    k_T = k_a + (k_s - k_a) * max (0.0_dp, (sigma - sigma_b) / sigma_c) * cs2**2

    theta_tropo = T_tropo * (p / p_0)**(-kappa)  ! potential temperature at tropopause

    theta_force = T_mean - delta_T * sn2 - delta_theta * cs2 * log (p / p_0)

    theta_equil = max (theta_tropo, theta_force) ! equilibrium temperature
  end subroutine cal_theta_eq
end module physics_Held_Suarez_mod
