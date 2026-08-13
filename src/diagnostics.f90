module diagnostics_mod

  ! Diagnostic routines and functions for physical quantities

  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, EDGE, MATH_PI, N_BDRY, N_VARIABLE, RT, DG, UP, S_MASS, S_TEMP, S_VELO, TRIAG, LORT, UPLT, &
       IJMINUS, IJPLUS, IMINUSJPLUS, IPLUSJMINUS, &
       c_p, compressible, grav_accel, H_rho, kappa, mode_split, omega, p_0, pressure_save, R_d, radius, ref_density, zlevels, &
       eps, p_top, z_null

  use domain_ops_mod, only : apply, apply_onescale_to_patch
  use init_mod,       only : surf_geopot
  use integrate_mod,  only : integrate_hex
  use patch_mod,      only : PATCH_SIZE
  use utils_mod,      only : dz_l, interp, porous_density, pressure_i, z_i

  use domain_mod,     only : Domain, Float_Field, dscalar, exner, exner_fun, grid, h_flux, ke, mass, temp, mean_m, mean_t, scalar, &
       sol, sol_mean, topography, trend, velo, velo1, velo2, vort, idx

  implicit none

  private
  public :: eta_e, free_surface
  public :: N_e, N_i
  public :: rho_dz_i
  public :: theta2temp, temp2theta
  public :: buoyancy, cal_buoyancy, density_i, density_e, potential_density, potential_buoyancy, theta_i
  public :: f_coriolis_edge, f_coriolis_node
  public :: cal_density, cal_div, cal_geopot, cal_surf_press, cal_temp, cal_vort, post_vort, integrate_pressure_up
  public :: curlv_e, div, gradi_e
  public :: kinetic_energy, pot_energy, pot_enstrophy, topo, total_ke, u_mag, umag
  public :: layer1_ke, layer2_ke, one_layer_ke, barotropic_ke, baroclinic_velocity, barotropic_velocity


contains


  function barotropic_ke (dom, i, j, zlev, offs, dims) result(val)
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer                     :: d, id, idE, idNE, idN, idW, idSW, idS
    integer,  dimension(1:EDGE) :: id_edge, id_node
    real(dp)                    :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp)                    :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S
    real(dp), dimension(1:EDGE) :: u

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idE  = idx (i+1, j,   offs, dims)
       idNE = idx (i+1, j+1, offs, dims)
       idN  = idx (i,   j+1, offs, dims)
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       id_node = [idE, idNE, idN]
       id_edge = [id,  id,   id ]
       u = barotropic_velo ()
       u_prim_RT = u(1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = u(2) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = u(3) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = u(1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = u(2) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = u(3) * dom%pedlen%elts(EDGE*id+UP+1)

       id_node = [ idW, idSW, idS ]
       id_edge = id_node
       u = barotropic_velo ()
       u_prim_RT_W  = u(1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_DG_SW = u(2) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_S  = u(3) * dom%len%elts(EDGE*idS +UP+1)

       u_dual_RT_W  = u(1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = u(2) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = u(3) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = &
            (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1) + &
            sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if

  contains

    function barotropic_velo ()
      real(dp), dimension(1:EDGE) :: barotropic_velo

      integer                     :: e, id_e, k
      real(dp), dimension(1:EDGE) :: dz

      do e = 1, EDGE
         id_e = EDGE*id_edge(e) + e

         dz = 0.0_dp
         do k = 1, 2
            dz(k) = interp (&
                 sol_mean(S_MASS,k)%data(d)%elts(id+1)         + sol(S_MASS,k)%data(d)%elts(id+1), &
                 sol_mean(S_MASS,k)%data(d)%elts(id_node(e)+1) + sol(S_MASS,k)%data(d)%elts(id_node(e)+1))
         end do

         barotropic_velo(e) = (dz(1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(2)*sol(S_VELO,2)%data(d)%elts(id_e)) / sum(dz)
      end do

    end function barotropic_velo

  end function barotropic_ke


  subroutine barotropic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate barotropic velocity in two-layer model

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer(4)                     :: d, e, id, id_e, id_i, idE, idNE, idN, k
    real(dp)                       :: dz0
    real(dp), dimension (1:EDGE,2) :: dz

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d = dom%id + 1

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    do k = 1, 2
       dz0 = sol_mean(S_MASS,k)%data(d)%elts(id+1) + sol(S_MASS,k)%data(d)%elts(id_i)
       dz(RT+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idE+1)  + sol(S_MASS,k)%data(d)%elts(idE+1))
       dz(DG+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idNE+1) + sol(S_MASS,k)%data(d)%elts(idNE+1))
       dz(UP+1,k) = interp (dz0, sol_mean(S_MASS,k)%data(d)%elts(idN+1)  + sol(S_MASS,k)%data(d)%elts(idN+1))
    end do

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = (dz(e,1)*sol(S_VELO,1)%data(d)%elts(id_e) + dz(e,2)*sol(S_VELO,2)%data(d)%elts(id_e)) / (dz(e,1) + dz(e,2))
    end do
  end subroutine barotropic_velocity


  subroutine baroclinic_velocity (dom, i, j, zlev, offs, dims)
    ! Calculate baroclinic velocity in top layer

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, e, id, id_e

    id = idx (i, j, offs, dims)
    d = dom%id + 1

    do e = 1, EDGE
       id_e = EDGE*id + e
       velo(id_e) = velo2(id_e) - velo1(id_e)
    end do
  end subroutine baroclinic_velocity


  function buoyancy (dom, i, j, zlev, offs, dims, q) result(val)
    ! at nodes
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(dp)                                  :: val

    integer  :: d, id_i
    real(dp) :: rho_dz, rho_dz_theta

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    rho_dz  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + q(S_MASS,zlev)%data(d)%elts(id_i)
    rho_dz_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + q(S_TEMP,zlev)%data(d)%elts(id_i)

    val = rho_dz_theta / rho_dz
  end function buoyancy


  subroutine cal_buoyancy (dom, i, j, zlev, offs, dims)
    ! at nodes

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id_i
    real(dp) :: rho_dz, rho_dz_theta

    id_i = idx (i, j, offs, dims) + 1

    rho_dz  = mean_m(id_i) + mass(id_i)
    rho_dz_theta = mean_t(id_i) + temp(id_i)

    scalar(id_i) = rho_dz_theta / rho_dz
  end subroutine cal_buoyancy



  subroutine cal_density (dom, i, j, zlev, offs, dims)
    ! Compute density
    ! *** compressible case requires pressure ***

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1) 

    integer  :: d, id
    real(dp) :: rho_dz_theta, exner

    id = idx (i, j, offs, dims) + 1

    if (compressible) then ! rho = P / (kappa theta pi)
       d = dom%id + 1
       rho_dz_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id) + sol(S_TEMP,zlev)%data(d)%elts(id)
       exner = c_p * (dom%press%elts(id)/p_0)**kappa

       scalar(id) = dom%press%elts(id) / (kappa * rho_dz_theta * exner) 
    else ! gravitational density (Boussinesq approximation)
       scalar(id) = ref_density * (1.0_dp - (mean_t(id) + temp(id)) / (mean_m(id) + mass(id)))
    end if
  end subroutine cal_density

  subroutine cal_div (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: id

    id = idx (i, j, offs, dims) + 1

    dscalar(id) = div (h_flux, dom, i, j, offs, dims)
  end subroutine cal_div

  subroutine cal_geopot (dom, i, j, zlev, offs, dims)
    ! Compute geopotential in compressible case
    ! Assumes that temperature has already been calculated and stored in exner_fun

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1) 

    integer :: id, d, k
    real(dp) :: rho_dz, pressure_lower, pressure_upper

    d = dom%id + 1
    id = idx(i, j, offs, dims) + 1

    ! Integrate geopotential upwards from surface
    rho_dz = sol_mean(S_MASS,1)%data(d)%elts(id) + sol(S_MASS,1)%data(d)%elts(id)

    pressure_lower = dom%surf_press%elts(id)
    pressure_upper = pressure_lower - grav_accel * rho_dz

    dom%geopot_lower%elts(id) = surf_geopot (d, id) / grav_accel

    k = 1
    do while (pressure_upper > pressure_save(1))
       dom%geopot_lower%elts(id) = dom%geopot_lower%elts(id) + &
            R_d/grav_accel * exner_fun(k)%data(d)%elts(id) * (log(pressure_lower)-log(pressure_upper))

       k = k+1
       rho_dz = sol_mean(S_MASS,k+1)%data(d)%elts(id) + sol(S_MASS,k+1)%data(d)%elts(id)

       pressure_lower = pressure_upper
       pressure_upper = pressure_lower - grav_accel * rho_dz
    end do

    ! Add additional contribution up to pressure level pressure_save(1)
    dom%geopot_lower%elts(id) = dom%geopot_lower%elts(id) &
         + R_d/grav_accel * exner_fun(k)%data(d)%elts(id) * (log(pressure_lower) - log(pressure_save(1)))
  end subroutine cal_geopot

  
  subroutine cal_surf_press (q)
    ! Compute surface pressure and save in press_lower for upward integration
    ! Set geopotential to surface geopotential for upward integration

    implicit none

    type(Float_Field), intent(inout), target :: q(N_VARIABLE,1:zlevels)

    integer :: d, k, p

    call apply (set_surf_geopot, z_null)

    do d = 1, size(grid)
       grid(d)%surf_press%elts = 0.0_dp
       do k = 1, zlevels
          mass   =>        q(S_MASS,k)%data(d)%elts
          temp   =>        q(S_TEMP,k)%data(d)%elts
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          do p = 3, grid(d)%patch%length
             call apply_onescale_to_patch (column_mass, grid(d), p-1, k, 0, 1)
          end do
          nullify (mass, mean_m, mean_t, temp)
       end do
       grid(d)%surf_press%elts  = grav_accel * grid(d)%surf_press%elts + p_top

       grid(d)%press_lower%elts = grid(d)%surf_press%elts
    end do
  end subroutine cal_surf_press


  subroutine column_mass (dom, i, j, zlev, offs, dims)
    ! Sum up mass

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id
    real(dp) :: rho_dz, rho_dz_theta

    id = idx (i, j, offs, dims) + 1

    rho_dz = mean_m(id) + mass(id)
    if (compressible) then
       dom%surf_press%elts(id) = dom%surf_press%elts(id) + rho_dz
    else
       rho_dz_theta = mean_t(id) + temp(id)
       
       dom%surf_press%elts(id) = dom%surf_press%elts(id) + (rho_dz - rho_dz_theta)
    end if
  end subroutine column_mass


  subroutine set_surf_geopot (dom, i, j, zlev, offs, dims)
    ! Set initial geopotential to surface geopotential (negative for incompressible ocean flows)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    if (compressible) then
       dom%geopot%elts(id) = surf_geopot (d, id)
    else
       dom%geopot%elts(id) = grav_accel * topography%data(d)%elts(id)
    end if
  end subroutine set_surf_geopot

   subroutine integrate_pressure_up (dom, i, j, zlev, offs, dims)
    ! Integrate pressure (compressible case)/Lagrange multiplier (incompressible case) and geopotential up from surface to top layer
    !
    ! Hydrostatic equilibrium:  dP = - g rho dz 
    ! compressible case:   rho dz = mu, rho = P / (kappa theta pi)
    ! incompressible case: rho dz = (1 - theta) mu = mu - Theta
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: d, id
    real(dp) :: dz, rho_dz, rho_dz_theta, p_upper
    
    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    rho_dz       = mean_m(id) + mass(id)
    rho_dz_theta = mean_t(id) + temp(id)
    dz           = rho_dz / porous_density (d, id, zlev)
    
    dom%geopot_lower%elts(id) = dom%geopot%elts(id)
    if (compressible) then ! compressible case:, rho = P / (kappa theta pi)
       p_upper = dom%press_lower%elts(id) - grav_accel * rho_dz

       dom%press%elts(id) = interp (dom%press_lower%elts(id), p_upper)

       exner(id) = c_p * (dom%press%elts(id)/p_0)**kappa

       dom%geopot%elts(id) = dom%geopot_lower%elts(id) + grav_accel * kappa * rho_dz_theta * exner(id) / dom%press%elts(id)
    else ! incompressible case
       p_upper = dom%press_lower%elts(id) - grav_accel * (rho_dz - rho_dz_theta)

       dom%press%elts(id) = interp (dom%press_lower%elts(id), p_upper)

       dom%geopot%elts(id) = dom%geopot_lower%elts(id) + grav_accel * dz
    end if
    dom%press_lower%elts(id) = p_upper
  end subroutine integrate_pressure_up


  subroutine cal_temp (dom, i, j, zlev, offs, dims)
    ! Compute temperature in compressible case

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1) 

    integer                        :: id, d, k
    real(dp)                       :: rho_dz, rho_dz_lower, rho_dz_theta
    real(dp), dimension(1:zlevels) :: p

    d = dom%id + 1
    id = idx(i, j, offs, dims) + 1

    ! Integrate the pressure upwards
    rho_dz = sol_mean(S_MASS,1)%data(d)%elts(id) + sol(S_MASS,1)%data(d)%elts(id)
    p(1) = dom%surf_press%elts(id) - 0.5_dp * grav_accel * rho_dz
    do k = 2, zlevels
       rho_dz_lower = rho_dz
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       p(k) = p(k-1) - grav_accel * interp (rho_dz, rho_dz_lower)
    end do

    ! Temperature at all vertical levels (saved in exner_fun)
    do k = 1, zlevels
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       rho_dz_theta = sol_mean(S_TEMP,k)%data(d)%elts(id) + sol(S_TEMP,k)%data(d)%elts(id)

       exner_fun(k)%data(d)%elts(id) = rho_dz_theta / rho_dz * (p(k)/p_0)**kappa
    end do

    ! Temperature at save levels (saved in trend)
    do k = 1, zlevels
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id) + sol(S_MASS,k)%data(d)%elts(id)
       rho_dz_theta = sol_mean(S_TEMP,k)%data(d)%elts(id) + sol(S_TEMP,k)%data(d)%elts(id)

       trend(1,k)%data(d)%elts(id) = rho_dz_theta / rho_dz * (pressure_save(k) / p_0)**kappa
    end do
  end subroutine cal_temp
  

  subroutine cal_vort (dom, i, j, zlev, offs, dims)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id, idE, idN
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_UP_E, u_prim_RT_N

    id  = idx (i,   j,   offs, dims)
    idE = idx (i+1, j,   offs, dims)
    idN = idx (i,   j+1, offs, dims)

    u_prim_RT   = velo(EDGE*id+RT +1) * dom%len%elts(EDGE*id +RT+1)
    u_prim_DG   = velo(EDGE*id+DG +1) * dom%len%elts(EDGE*id +DG+1)
    u_prim_UP   = velo(EDGE*id+UP +1) * dom%len%elts(EDGE*id +UP+1)
    u_prim_UP_E = velo(EDGE*idE+UP+1) * dom%len%elts(EDGE*idE+UP+1)
    u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)

    if (dom%triarea%elts(TRIAG*id+LORT+1) > eps (radius**2)) then
       vort(TRIAG*id+LORT+1) =   (u_prim_RT + u_prim_UP_E + u_prim_DG  ) / dom%triarea%elts(TRIAG*id+LORT+1)
    else
       vort(TRIAG*id+LORT+1) = 0.0_dp
    endif

    if (dom%triarea%elts(TRIAG*id+UPLT+1) > eps (radius**2)) then
       vort(TRIAG*id+UPLT+1) = - (u_prim_DG + u_prim_UP   + u_prim_RT_N) / dom%triarea%elts(TRIAG*id+UPLT+1)
    else
       vort(TRIAG*id+UPLT+1) = 0.0_dp
    end if
  end subroutine cal_vort
  
  
  subroutine post_vort (dom, p, c, offs, dims, zlev)
    ! Correct values for vorticity and qe at pentagon points
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: p, c, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer  :: id, idS, idW, idSW, idN
    real(dp) :: u_prim_RT, u_prim_RT_N, u_prim_RT_SW, u_prim_RT_W, u_prim_DG_SW, u_prim_UP, u_prim_UP_S, u_prim_UP_SW

    ! Parts 4, 5 of hexagon IJMINUS  (lower left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idSW+DG+1) = 0 in this case
    if (c == IJMINUS) then
       idS  = idx ( 0, -1, offs, dims)
       idSW = idx (-1, -1, offs, dims)
       idW  = idx (-1,  0, offs, dims)

       u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)
       u_prim_RT_SW = velo(EDGE*idSW+RT+1) * dom%len%elts(EDGE*idSW+RT+1) 
       u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)

       vort(TRIAG*idSW+LORT+1) = (u_prim_UP_S - u_prim_RT_W + u_prim_RT_SW) / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idSW+UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    ! Parts 5, 6 of hexagon IPLUSJMINUS (lower right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idS+UP+1) = 0 in this case
    if (c == IPLUSJMINUS) then 
       id   = idx (PATCH_SIZE,    0, offs, dims)
       idS  = idx (PATCH_SIZE,   -1, offs, dims)
       idSW = idx (PATCH_SIZE-1, -1, offs, dims)

       u_prim_RT_SW = velo(EDGE*idSW+RT+1)*dom%len%elts(EDGE*idSW+RT+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1)*dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT    = velo(EDGE*id  +RT+1)*dom%len%elts(EDGE*id  +RT+1)

       vort(TRIAG*idSW+LORT+1) = (- u_prim_RT + u_prim_RT_SW + u_prim_DG_SW) / dom%triarea%elts(TRIAG*idSW+LORT+1)
       vort(TRIAG*idS +UPLT+1) = vort(TRIAG*idSW+LORT+1)
    end if

    ! Parts 3, 4 of hexagon IMINUSJPLUS (upper left corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*idW+RT+1) = 0 in this case
    if (c == IMINUSJPLUS) then
       id   = idx ( 0, PATCH_SIZE,   offs, dims)
       idSW = idx (-1, PATCH_SIZE-1, offs, dims)
       idW  = idx (-1, PATCH_SIZE,   offs, dims)

       u_prim_UP    = velo(EDGE*id  +UP+1) * dom%len%elts(EDGE*id  +UP+1)
       u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_UP_SW = velo(EDGE*idSW+UP+1) * dom%len%elts(EDGE*idSW+UP+1)

       vort(TRIAG*idSW+UPLT+1) = (u_prim_UP - u_prim_DG_SW - u_prim_UP_SW) / dom%triarea%elts(TRIAG*idSW+UPLT+1)  
       vort(TRIAG*idW +LORT+1) = vort(TRIAG*idSW+UPLT+1)
    end if

    ! Parts 1, 2 of hexagon IJPLUS (upper right corner of lozenge) combined to form pentagon
    ! Note that pedlen(EDGE*id+DG+1) = 0 in this case
    if (c == IJPLUS) then 
       id  = idx (PATCH_SIZE, PATCH_SIZE,   offs, dims)
       idN = idx (PATCH_SIZE, PATCH_SIZE+1, offs, dims)

       u_prim_RT   = velo(EDGE*id +RT+1) * dom%len%elts(EDGE*id +RT+1)
       u_prim_RT_N = velo(EDGE*idN+RT+1) * dom%len%elts(EDGE*idN+RT+1)
       u_prim_UP   = velo(EDGE*id +UP+1) * dom%len%elts(EDGE*id +UP+1)

       vort(TRIAG*id+LORT+1) = (u_prim_RT - u_prim_RT_N - u_prim_UP) / dom%triarea%elts(TRIAG*id+LORT+1)
       vort(TRIAG*id+UPLT+1) = vort(TRIAG*id+LORT+1)
    end if
  end subroutine post_vort


  function curlv_e (curl, dom, i, j, offs, dims) result(val)
    ! Curl of vorticity given at triangle circumcentres x_v, rot(rot(u))
    ! output is at edges x_e

    implicit none

    real(dp), pointer, intent(in)    :: curl(:)
    type(Domain),      intent(inout) :: dom
    integer,           intent(in)    :: i, j
    integer,           intent(in)    :: offs(N_BDRY+1)
    integer,           intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                        :: val(3)
    
    integer :: id, idS, idW

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)

    val(RT+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*idS+UPLT+1)) / dom%pedlen%elts(EDGE*id+RT+1)
    val(DG+1) = (curl(TRIAG*id +LORT+1) - curl(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+DG+1)
    val(UP+1) = (curl(TRIAG*idW+LORT+1) - curl(TRIAG*id +UPLT+1)) / dom%pedlen%elts(EDGE*id+UP+1)
  end function curlv_e

  
  function density_i (dom, i, j, zlev, offs, dims, q) result(val)
    ! Density at nodes

    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val

    integer  :: d, id
    real(dp) :: rho_dz, rho_dz_theta, p, temp

    d = dom%id + 1
    id  = idx (i, j, offs, dims) + 1

    if (compressible) then ! rho = P / (kappa theta pi)
       rho_dz_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id) + q(S_TEMP,zlev)%data(d)%elts(id)
       rho_dz       = sol_mean(S_MASS,zlev)%data(d)%elts(id) + q(S_MASS,zlev)%data(d)%elts(id)

       p    = pressure_i (dom, i, j, zlev, offs, dims, sol)

       temp = (rho_dz_theta / rho_dz) * (p/p_0)**kappa ! temperature

       val = p / (R_d * temp) 
    else                   ! gravitational density using Boussinesq approximation
       val = ref_density * (1.0_dp - buoyancy (dom, i, j, zlev, offs, dims, q))
    end if
  end function density_i


  function density_e (dom, i, j, zlev, offs, dims, q) result(val)
    ! Density at edges
    ! *** compressible case requires pressure at zlev ***

    implicit none

    type(Domain), intent(inout)           :: dom
    integer,      intent(in)              :: i, j, zlev
    integer,      intent(in)              :: offs(N_BDRY+1) 
    integer,      intent(in)              :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val(EDGE)

    real(dp) :: rho(0:EDGE) 

    rho(0)    = density_i (dom, i,   j,   zlev, offs, dims, q)
    rho(RT+1) = density_i (dom, i+1, j,   zlev, offs, dims, q)
    rho(DG+1) = density_i (dom, i+1, j+1, zlev, offs, dims, q)
    rho(UP+1) = density_i (dom, i  , j+1, zlev, offs, dims, q)

    val= 0.5 * (rho(0) + rho(1:EDGE))
  end function density_e


  function div (hflux, dom, i, j, offs, dims) result(val)
    ! Divergence at nodes x_i given horizontal fluxes at edges x_e

    implicit none

    real(dp), pointer, intent(in)    :: hflux(:)
    type(Domain),      intent(inout) :: dom
    integer,           intent(in)    :: i, j
    integer,           intent(in)    :: offs(N_BDRY+1)
    integer,           intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                        :: val

    integer :: id, idW, idS, idSW

    id   = idx (i,   j,   offs, dims)
    idS  = idx (i,   j-1, offs, dims)
    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)

    val = (hflux(EDGE*id+RT+1)-hflux(EDGE*idW+RT+1) + hflux(EDGE*idSW+DG+1)-hflux(EDGE*id+DG+1) &
         + hflux(EDGE*id+UP+1)-hflux(EDGE*idS+UP+1)) * dom%areas%elts(id+1)%hex_inv
  end function div


  function eta_e (dom, i, j, zlev, offs, dims, q) result(val)
    ! Free surface perturbation at edges
    implicit none
    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(dp)                                  :: val(EDGE)   

    integer  :: d, id, idE, idN, idNE
    real(dp) :: eta0

    d = dom%id + 1
    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    if (mode_split) then
       val(RT+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idE+1))
       val(DG+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idNE+1))
       val(UP+1) = interp (q(S_MASS,zlevels+1)%data(d)%elts(id+1), q(S_MASS,zlevels+1)%data(d)%elts(idN+1))
    else
       eta0 = free_surface (dom, i, j, zlev, offs, dims, q)
       val(RT+1) = interp (eta0, free_surface (dom, i+1, j,   zlev, offs, dims, q))
       val(DG+1) = interp (eta0, free_surface (dom, i+1, j+1, zlev, offs, dims, q))
       val(UP+1) = interp (eta0, free_surface (dom, i,   j+1, zlev, offs, dims, q))
    end if
  end function eta_e


  function f_coriolis_edge (dom, i, j, zlev, offs, dims) result(val)
    ! Coriolis parameter at edges

    implicit none

    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i, j, zlev
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)

    real(dp)                 :: val(EDGE) 

    integer :: id

    id = idx (i, j, offs, dims)

    val = dom%midpt%elts(EDGE*id+RT+1:EDGE*id+UP+1)%z / radius * 2*omega
  end function f_coriolis_edge


  function f_coriolis_node (dom, i, j, zlev, offs, dims) result(val)
    ! Coriolis parameter at nodes

    implicit none

    type(Domain), intent(in) :: dom
    integer,      intent(in) :: i, j, zlev
    integer,      intent(in) :: offs(N_BDRY+1)
    integer,      intent(in) :: dims(2,N_BDRY+1)
    real(dp)                 :: val

    integer :: id

    id = idx (i, j, offs, dims)

    val = dom%node%elts(id+1)%z / radius * 2*omega
  end function f_coriolis_node


  function free_surface (dom, i, j, zlev, offs, dims, q) result(val)
    ! Computes free surface perturbations

    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1)
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val

    integer  :: d, id_i, k
    real(dp) :: rho_dz, total_depth

    d = dom%id + 1
    id_i  = idx (i, j, offs, dims) + 1

    total_depth = 0.0_dp
    do k = 1, zlevels
       rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id_i) + q(S_MASS,k)%data(d)%elts(id_i)
       total_depth = total_depth + rho_dz /  porous_density (d, id_i, k)
    end do

    val = total_depth + topography%data(d)%elts(id_i)
  end function free_surface


  function gradi_e (scalar, dom, i, j, offs, dims) result(val)
    ! Gradient of a scalar at nodes x_i
    ! output is at edges

    implicit none

    real(dp), pointer, intent(in) :: scalar(:)
    type(Domain), intent(inout)   :: dom
    integer,      intent(in)      :: i, j
    integer,      intent(in)      :: offs(N_BDRY+1)
    integer,      intent(in)      :: dims(2,N_BDRY+1)
    real(dp)                      :: val(3)

    integer :: id, idE, idN, idNE

    id   = idx (i,   j,   offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    val(RT+1) = (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1)
    val(DG+1) = (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1)
    val(UP+1) = (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1)
  end function gradi_e


  function kinetic_energy (dom, i, j, zlev, offs, dims) result(val)
    ! Kinetic energy u^2/2 at level zlev using approximation to TRiSK formula

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    u_prim_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_prim_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_prim_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

    u_dual_RT = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP = sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

    u_prim_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
    u_prim_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_prim_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

    u_dual_RT_W  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
    u_dual_DG_SW = sol(S_VELO,zlev)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
    u_dual_UP_S  = sol(S_VELO,zlev)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

    val = (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT + &
         u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
         * dom%areas%elts(id+1)%hex_inv/4.0_dp 
  end function kinetic_energy

  function layer1_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1))  &
            * (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4.0_dp
    else
       val = 0.0_dp
    end if
  end function layer1_ke


  function layer2_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,2)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,2)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,2)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,2)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,2)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,2)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,2)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,2)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,2)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,2)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,2)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,2)%data(d)%elts(id+1) + sol(S_MASS,2)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if
  end function layer2_ke


  function N_i (dom, i, j, zlev, offs, dims) result(val)
    ! Brunt-Vaisala frequency at nodes at layer interface above layer zlev
    ! *** need zlev /= zlevels ***

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    real(dp)                                      :: val

    integer  :: d, id
    real(dp) :: drho, dtheta, dz, rho_l, rho1, rho2, theta_l, theta1, theta2

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    dz = dz_l (dom, i, j, zlev, offs, dims, sol)

    if (compressible) then
       theta1 = theta_i (dom, i, j, zlev,   offs, dims)
       theta2 = theta_i (dom, i, j, zlev+1, offs, dims)
       theta_l = interp (theta1, theta2)
       dtheta =  theta2 - theta1

       val = sqrt ( grav_accel *  dtheta/dz / theta_l) 
    else                     ! incompressible
       rho1 = porous_density (d, id+1, zlev)
       rho2 = porous_density (d, id+1, zlev+1) 
       rho_l = interp (rho1, rho2)
       drho = rho2 - rho1

       val = sqrt ( grav_accel * drho/dz / rho_l )
    end if
  end function N_i


  function N_e (dom, i, j, zlev, offs, dims) result(val)
    ! Brunt-Vaisala frequency at edges at layer interfaces

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    real(dp)                                      :: val(EDGE)

    val = 0.5 * ( &
         N_i (dom, i,   j,   zlev, offs, dims) + [ &
         N_i (dom, i+1, j,   zlev, offs, dims),     &
         N_i (dom, i+1, j+1, zlev, offs, dims),     &
         N_i (dom, i,   j+1, zlev, offs, dims) ])
  end function N_e


  function one_layer_ke (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    id = idx (i, j, offs, dims)

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       idW  = idx (i-1, j,   offs, dims)
       idSW = idx (i-1, j-1, offs, dims)
       idS  = idx (i,   j-1, offs, dims)

       d  = dom%id + 1

       u_prim_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
       u_prim_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
       u_prim_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

       u_dual_RT = sol(S_VELO,1)%data(d)%elts(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
       u_dual_DG = sol(S_VELO,1)%data(d)%elts(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
       u_dual_UP = sol(S_VELO,1)%data(d)%elts(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

       u_prim_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
       u_prim_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
       u_prim_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

       u_dual_RT_W  = sol(S_VELO,1)%data(d)%elts(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
       u_dual_DG_SW = sol(S_VELO,1)%data(d)%elts(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
       u_dual_UP_S  = sol(S_VELO,1)%data(d)%elts(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

       val = (sol_mean(S_MASS,1)%data(d)%elts(id+1) + sol(S_MASS,1)%data(d)%elts(id+1)) * &
            (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT &
            + u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
            * dom%areas%elts(id+1)%hex_inv/4
    else
       val = 0.0_dp
    end if
  end function one_layer_ke

  function pot_energy (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: id
    real(dp) :: rho_dz

    id = idx (i, j, offs, dims)

    rho_dz = sol_mean(S_MASS,zlev)%data(dom%id+1)%elts(id+1) + sol(S_MASS,zlev)%data(dom%id+1)%elts(id+1)
    val = rho_dz**2
  end function pot_energy


  function pot_enstrophy (dom, i, j, zlev, offs, dims) result(val)
    ! Computes potential enstrophy in two layer mode split case
    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: d, id, id_i
    real(dp) :: f, h, w

    id = idx (i, j, offs, dims)
    id_i = id + 1

    d = dom%id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       ! Approximate Coriolis term
       f = dom%coriolis%elts(TRIAG*id+LORT+1)/dom%triarea%elts(TRIAG*id+LORT+1)
       ! Total vorticity
       w = dom%ke%elts(id_i)  + f
       ! Height
       if (zlev == 3) then ! barotropic
          h = (sol_mean(S_MASS,1)%data(d)%elts(id_i) + sol(S_MASS,1)%data(d)%elts(id_i) + &
               sol_mean(S_MASS,2)%data(d)%elts(id_i) + sol(S_MASS,2)%data(d)%elts(id_i)) / ref_density
       else ! single layer
          h = (sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i)) / ref_density
       end if
       val = 0.5_dp * (w / h)**2
    else
       val = 0.0_dp
    end if
  end function pot_enstrophy


  function potential_density (dom, i, j, zlev, offs, dims, q) result(val)
    ! Potential density at nodes (neglect free surface perturbation)

    implicit none

    type(Domain),      intent(in)         :: dom
    integer,           intent(in)         :: i, j, zlev
    integer,           intent(in)         :: offs(N_BDRY+1) 
    integer,           intent(in)         :: dims(2,N_BDRY+1) 
    type(Float_Field), intent(in), target :: q(:,:)
    real(dp)                              :: val

    integer  :: d, id_i
    real(dp) :: z

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims, q)

    val = density_i (dom, i, j, zlev, offs, dims, q) - ref_density * z / H_rho
  end function potential_density


  function potential_buoyancy (dom, i, j, zlev, offs, dims, q) result(val)
    ! Potential buoyancy at nodes (neglect free surface perturbation)

    implicit none

    type(Domain)                              :: dom
    integer                                   :: i, j, zlev
    integer, dimension(N_BDRY+1)              :: offs
    integer, dimension(2,N_BDRY+1)            :: dims
    type(Float_Field), dimension(:,:), target :: q
    real(dp)                                  :: val

    integer  :: d, id_i
    real(dp) :: z

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    z = z_i (dom, i, j, zlev, offs, dims, q)

    val = buoyancy (dom, i, j, zlev, offs, dims, q) + z / H_rho
  end function potential_buoyancy


  function rho_dz_i (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims)

    val = sol(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(id+1)
  end function rho_dz_i


  function theta_i (dom, i, j, zlev, offs, dims) result(val)
    ! Potential temperature at layer centre

    implicit none

    type(Domain),                   intent(inout) :: dom
    integer,                        intent(in)    :: i, j, zlev
    integer, dimension(N_BDRY+1),   intent(in)    :: offs
    integer, dimension(2,N_BDRY+1), intent(in)    :: dims
    real(dp)                                      :: val

    integer  :: d, id
    real(dp) :: rho_dz, rho_dz_theta

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    rho_dz_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id) + sol(S_TEMP,zlev)%data(d)%elts(id)
    rho_dz       = sol_mean(S_MASS,zlev)%data(d)%elts(id) + sol(S_MASS,zlev)%data(d)%elts(id) 

    val = rho_dz_theta / rho_dz
  end function theta_i

  function theta2temp (theta, p) result(val)
    ! Convert potential temperature (theta) to temperature
    ! requires pressure
    implicit none
    real(dp), intent(in) :: p, theta
    real(dp)             :: val

    val = theta * (p / p_0)**kappa
  end function theta2temp


  function temp2theta (temp, p) result(val)
    ! Convert temperature (temp) to potential temperature
    ! requires pressure
    implicit none
    real(dp), intent(in) :: p, temp
    real(dp)             :: val

    val = temp / (p / p_0)**kappa
  end function temp2theta


  function total_ke (itype) result(val)
    ! Computes total kinetic energy

    implicit none

    character(len=*), intent(in) :: itype
    real(dp)                     :: val

    integer :: k

    val = 0.0_dp
    do k = 1, zlevels
       val = val + integrate_hex (kinetic_energy, k)
    end do
  end function total_ke


  function topo (dom, i, j, zlev, offs, dims) result(val)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer :: d, id

    d = dom%id + 1
    id = idx (i, j, offs, dims) + 1

    val = topography%data(d)%elts(id)
  end function topo


  function u_mag (dom, i, j, zlev, offs, dims) result(val)
    ! Velocity magnitude using data from a single element

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: val

    integer  :: id
    real(dp) :: prim_dual(1:EDGE), u(1:EDGE)

    id = idx (i, j, offs, dims)

    prim_dual = dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+RT+1:EDGE*id+UP+1)

    u = sol(S_VELO,zlev)%data(dom%id+1)%elts(EDGE*id+RT+1:EDGE*id+UP+1)

    val = sqrt (sum (u**2 * prim_dual) * dom%areas%elts(id+1)%hex_inv)
  end function u_mag


  subroutine umag (q)
    ! Evaluate complete velocity trend by adding gradient terms to previously calculated source terms on entire grid
    !
    ! Input:  velocity field q at a single vertical layer
    ! Output: velocity magnitude stored in dom%ke

    implicit none

    type(Float_Field), target, intent(in) :: q

    integer :: k

    integer :: d, p

    do d = 1, size(grid)
       velo => q%data(d)%elts
       ke   => grid(d)%ke%elts
       do p = 3, grid(d)%patch%length
          call apply_onescale_to_patch (cal_umag, grid(d), p-1, k, 0, 1)
       end do
       nullify (ke, velo)
    end do
  end subroutine umag


  subroutine cal_umag (dom, i, j, zlev, offs, dims)
    ! Velocity magnitude: sqrt(2*ke) using approximation to TRiSK formula
    ! divide out surface area

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, idW, idSW, idS
    real(dp) :: u_prim_RT, u_prim_DG, u_prim_UP, u_prim_RT_W, u_prim_DG_SW, u_prim_UP_S
    real(dp) :: u_dual_RT, u_dual_DG, u_dual_UP, u_dual_RT_W, u_dual_DG_SW, u_dual_UP_S

    d  = dom%id + 1
    id = idx (i, j, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    u_prim_RT = velo(EDGE*id+RT+1) * dom%len%elts(EDGE*id+RT+1)
    u_prim_DG = velo(EDGE*id+DG+1) * dom%len%elts(EDGE*id+DG+1)
    u_prim_UP = velo(EDGE*id+UP+1) * dom%len%elts(EDGE*id+UP+1)

    u_dual_RT = velo(EDGE*id+RT+1) * dom%pedlen%elts(EDGE*id+RT+1)
    u_dual_DG = velo(EDGE*id+DG+1) * dom%pedlen%elts(EDGE*id+DG+1)
    u_dual_UP = velo(EDGE*id+UP+1) * dom%pedlen%elts(EDGE*id+UP+1)

    u_prim_UP_S  = velo(EDGE*idS +UP+1) * dom%len%elts(EDGE*idS +UP+1)
    u_prim_DG_SW = velo(EDGE*idSW+DG+1) * dom%len%elts(EDGE*idSW+DG+1)
    u_prim_RT_W  = velo(EDGE*idW +RT+1) * dom%len%elts(EDGE*idW +RT+1)

    u_dual_RT_W  = velo(EDGE*idW +RT+1) * dom%pedlen%elts(EDGE*idW +RT+1)
    u_dual_DG_SW = velo(EDGE*idSW+DG+1) * dom%pedlen%elts(EDGE*idSW+DG+1)         
    u_dual_UP_S  = velo(EDGE*idS +UP+1) * dom%pedlen%elts(EDGE*idS +UP+1)

    ke(id+1) = sqrt( (u_prim_UP   * u_dual_UP   + u_prim_DG    * u_dual_DG    + u_prim_RT   * u_dual_RT + &
         u_prim_UP_S * u_dual_UP_S + u_prim_DG_SW * u_dual_DG_SW + u_prim_RT_W * u_dual_RT_W) &
         * dom%areas%elts(id+1)%hex_inv/(4.0_dp*MATH_PI*radius**2)) 
  end subroutine cal_umag

  
end module diagnostics_mod
