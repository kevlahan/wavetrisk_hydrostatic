Module test_case_mod
  ! Module file for tke test case
  use shared_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use equation_of_state_mod
  implicit none

  ! Standard variables
  integer                              :: CP_EVERY, resume_init
  real(8)                              :: dt_cfl,  tau_diffusion, total_cpu_time
  real(8)                              :: dPdim, R_ddim, specvoldim, dTempdim
  real(8), target                      :: bottom_friction_case

  ! Local variables
  real(8)                              :: drho
  real(8)                              :: N_0, Q_0, r_max, r_max_loc, tau_0, u_0
  character(255)                       :: coords
  logical                              :: kato
contains
  subroutine assign_functions
    ! Assigns generic pointer functions to functions defined in test cases
    implicit none

    ! Standard functions
    apply_initial_conditions => apply_initial_conditions_case
    dump                     => dump_case
    load                     => load_case
    initialize_a_b_vert      => initialize_a_b_vert_case
    initialize_dt_viscosity  => initialize_dt_viscosity_case
    initialize_thresholds    => initialize_thresholds_case
    physics_scalar_flux      => physics_scalar_flux_case
    physics_velo_source      => physics_velo_source_case
    set_save_level           => set_save_level_case
    set_thresholds           => set_thresholds_case
    surf_geopot              => surf_geopot_case
    update                   => update_case
    z_coords                 => z_coords_case

    ! Needed for vertical diffusion
    bottom_friction  => bottom_friction_case
    bottom_buoy_flux => bottom_buoy_flux_case
    top_buoy_flux    => top_buoy_flux_case
    wind_flux        => wind_flux_case
    tau_mag          => tau_mag_case
  end subroutine assign_functions

  function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type)
    ! Additional physics for the flux term of the scalar trend
    ! In this test case we add -gradient to the flux to include a Laplacian diffusion (div grad) to the scalar trend
    !
    ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
    use domain_mod
    implicit none

    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
    type(domain)                                         :: dom
    integer                                              :: d, id, idE, idNE, idN, v, zlev
    logical, optional                                    :: type

    physics_scalar_flux_case = 0.0_8
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    ! Additional physics for the source term of the velocity trend
    ! wind stress and bottom friction are included as surface fluxes in the split eddy viscosity split step
    implicit none

    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    physics_velo_source_case = 0.0_8 
  end function physics_velo_source_case

  subroutine read_test_case_parameters
    implicit none
    integer            :: ilat, ilon, k
    integer, parameter :: fid = 500
    real(8)            :: lat, lon, press_save
    character(255)     :: filename, varname

    ! Find input parameters file name
    if (command_argument_count() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'test_case.in'
    end if
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

    open(unit=fid, file=filename, action='READ')
    read (fid,*) varname, test_case
    read (fid,*) varname, run_id
    read (fid,*) varname, kato
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, dt_init
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)

    press_save = 0.0_8
    allocate (pressure_save(1))
    pressure_save(1) = press_save
    dt_write = dt_write * HOUR
    time_end = time_end * HOUR
    resume   = resume_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none

    call set_save_level

    call cal_r_max

    if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case                      = ", trim (test_case)
       write (6,'(A,A)')      "run_id                         = ", trim (run_id)
       write (6,'(A,L1)')     "kato                           = ", kato
       write (6,'(A,L1)')     "compressible                   = ", compressible
       write (6,'(A,L1)')     "penalize                       = ", penalize
       write (6,'(A,i3)')     "min_level                      = ", min_level
       write (6,'(A,i3)')     "max_level                      = ", max_level
       write (6,'(A,i3)')     "level_fill                     = ", level_fill
       write (6,'(A,i5)')     "number of domains              = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors           = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL                   = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL                    = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels                        = ", zlevels
       write (6,'(A,L1)')     "remap                          = ", remap
       write (6,'(A,L1)')     "sigma_z                        = ", sigma_z
       write (6,'(a,a)')      "remapscalar_type               = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type                 = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap                         = ", iremap
       write (6,'(A,L1)')     "default_thresholds             = ", default_thresholds
       write (6,'(A,es10.4)') "tolerance                      = ", tol
       write (6,'(A,i1)')     "optimize_grid                  = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt                       = ", adapt_dt
       write (6,'(A,es10.4)') "dt_init                        = ", dt_init
       write (6,'(a,a)')      "timeint_type                   = ", trim (timeint_type)
       write (6,'(A,i1)')     "Laplace_order                  = ", Laplace_order_init
       write (6,'(A,i1)')     "n_diffuse                      = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write [h]                   = ", dt_write/HOUR
       write (6,'(A,i6)')     "CP_EVERY                       = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance                      = ", rebalance
       write (6,'(A,es10.4)') "time_end [h]                   = ", time_end/HOUR
       write (6,'(A,i6)')     "resume                         = ", resume_init

       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius                   [km]  = ", radius / KM
       write (6,'(A,es11.4)') "omega                 [rad/s]  = ", omega
       write (6,'(A,es10.4)') "ref density          [kg/m^3]  = ", ref_density
       write (6,'(A,es10.4)') "grav accel            [m/s^2]  = ", grav_accel

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es11.4)') "max_depth                 [m]  = ", max_depth
       write (6,'(A,es11.4)') "min_depth                 [m]  = ", min_depth
       write (6,'(A,es11.4)') "c0 wave speed           [m/s]  = ", wave_speed
       write (6,'(A,es11.4)') "max wind stress       [N/m^2]  = ", tau_0
       write (6,'(A,es11.4)') "bottom friction         [m/s]  = ", bottom_friction_case
       write (6,'(A,es11.4)') "bottom drag decay         [d]  = ", 1/bottom_friction_case / DAY
       write (6,'(A,es11.4)') "dx_max                   [km]  = ", dx_max   / KM
       write (6,'(A,es11.4)') "dx_min                   [km]  = ", dx_min   / KM
       write (6,'(a,es11.4)') "r_max                          = ", r_max
       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'

       call print_density
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)

    if (rank == 0) then
       write (6,'(a,es13.6,4(a,es8.2),a,i2,a,i12,3(a,es9.2,1x))') &
            'time [h] = ', time/HOUR, &
            ' dt [s] = ', dt, &
            '  mass threshold = ', sum (threshold(S_MASS,:))/zlevels, &
            ' temp threshold = ', sum (threshold(S_TEMP,:))/zlevels, &
            ' velo threshold = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance

       write (12,'(5(es15.9,1x),i2,1x,i12,1x,3(es15.9,1x))')  time/DAY, dt, &
            threshold(S_MASS,zlevels), threshold(S_TEMP,zlevels), threshold(S_VELO,zlevels), &
            level_end, sum (n_active), min_mass, mass_error, rel_imbalance
    end if
  end subroutine print_log

  subroutine avg_temp (iwrt)
    ! Saves  temperature averaged over the sphere
    ! (assumes non-adaptive grid)
    use io_mod
    implicit none
    integer :: iwrt

    integer                       :: k
    real(8)                       :: area, dzl, eta, Nsq, z_k, z_s
    real(8), parameter            :: lat_band = 30
    real(8), dimension(1:zlevels) :: dz, T_avg
    real(8), dimension(0:zlevels) :: Kt_avg, Kv_avg, z
    character(4)                  :: s_time

    area = integrate_hex (area_fun, 1, level_start)
    do k = 1, zlevels
       Kt_avg(k) = integrate_hex (Kt_fun, k, level_start)
       Kv_avg(k) = integrate_hex (Kv_fun, k, level_start)
    end do
    do k = 1, zlevels
       T_avg(k)  = integrate_hex (temp_fun, k, level_start)
    end do
    Kt_avg = Kt_avg / area
    Kv_avg = Kv_avg / area
    T_avg  = T_avg  / area

    if (rank == 0) then
       eta = 0.0_8
       z_s = max_depth
       if (sigma_z) then
          z = z_coords_case (eta, z_s)
       else
          z = a_vert * eta + b_vert * z_s
       end if
       dz = z(1:zlevels) - z(0:zlevels-1)

       write (s_time, '(i4.4)') iwrt
       open (unit=20, file=trim(run_id)//'.6.'//s_time, form="FORMATTED", action='WRITE', status='REPLACE')

       write (6,'(a, f4.1, a)') "Temperature profile at time ", time/HOUR, " h"
       write (6,'(a)') "Level    z_l     Kt(z_l)     Kv(z_l)     z_k    T(z_k)        rho(k)"

       z_k = interp (z(0), z(1))
       write (6,'(i3, 3x, f7.2, 1x, 2(es11.4, 1x), f7.2, 1x, es11.4, 1x, es14.7)') &
            0, z(0), Kt_avg(0), Kv_avg(0), z_k, T_avg(1), density_eos (S_ref, T_avg(1),  z_k)
       write (20,'(i3, 1x, 5(es13.6,1x))')    0, z(0), Kt_avg(0), Kv_avg(0), z_k, T_avg(1)

       do k = 1, zlevels
          z_k = interp (z(k-1), z(k))
          write (6,'(i3, 3x, f7.2, 1x, 2(es11.4, 1x), f7.2, 1x, es11.4, 1x, es14.7)') &
               k, z(k), Kt_avg(k), Kv_avg(k), z_k, T_avg(k), density_eos (S_ref, T_avg(k), z_k)
          write (20,'(i3, 1x, 5(es13.6,1x))') k, z(k), Kt_avg(k), Kv_avg(k), z_k, T_avg(k)
       end do

       write (6,'(a, es9.2)') " "
       write (6,'(a)') "Level    z_l      Nsq"
       do k = 1, zlevels-1
          drho = density_eos (S_ref, T_avg(k+1), interp (z(k), z(k+1))) - density_eos (S_ref, T_avg(k), interp (z(k-1), z(k)))
          dzl = interp (dz(k), dz(k+1))
          Nsq = - grav_accel * drho/dzl / ref_density 
          write (6, '(i3, 3x, f7.2, 1x, es11.4, 1x)') k, z(k), Nsq
       end do
       write (6,'(a, es9.2)') " "
       close (20)
    end if
  contains
    real(8) function temp_fun (dom, i, j, zlev, offs, dims)
      ! Defines mass for total mass integration
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i
      real(8) :: density, full_mass, full_theta, lat, lon, z

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1

      call cart2sph (dom%node%elts(id_i), lon, lat)

      if (abs(lat/DEG) <= lat_band) then
         full_mass  = sol_mean(S_MASS,zlev)%data(d)%elts(id_i) + sol(S_MASS,zlev)%data(d)%elts(id_i) 
         full_theta = sol_mean(S_TEMP,zlev)%data(d)%elts(id_i) + sol(S_TEMP,zlev)%data(d)%elts(id_i)

         z = z_i (dom, i, j, zlev, offs, dims, sol)
         density = ref_density * (1 - full_theta/full_mass)
         temp_fun = temperature_eos (density, S_ref, z)
      else
         temp_fun = 0.0_8
      end if
    end function temp_fun

    real(8) function Kt_fun (dom, i, j, zlev, offs, dims)
      ! Defines mass for total mass integration
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i
      real(8) :: lat, lon

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1

      call cart2sph (dom%node%elts(id_i), lon, lat)

      if (abs(lat/DEG) <= lat_band) then
         Kt_fun = Kt(zlev)%data(d)%elts(id_i)
      else
         Kt_fun = 0.0_8
      end if
    end function Kt_fun

    real(8) function Kv_fun (dom, i, j, zlev, offs, dims)
      ! Defines mass for total mass integration
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: d, id_i
      real(8) :: lat, lon

      d = dom%id + 1
      id_i = idx (i, j, offs, dims) + 1

      call cart2sph (dom%node%elts(id_i), lon, lat)

      if (abs(lat/DEG) <= lat_band) then
         Kv_fun = Kv(zlev)%data(d)%elts(id_i)
      else
         Kv_fun = 0.0_8
      end if
    end function Kv_fun

    real(8) function area_fun (dom, i, j, zlev, offs, dims)
      ! Defines mass for total mass integration
      implicit none
      type(Domain)                   :: dom
      integer                        :: i, j, zlev
      integer, dimension(N_BDRY+1)   :: offs
      integer, dimension(2,N_BDRY+1) :: dims

      integer :: id_i
      real(8) :: lat, lon

      id_i = idx (i, j, offs, dims) + 1

      call cart2sph (dom%node%elts(id_i), lon, lat)

      if (abs(lat/DEG) <= lat_band) then
         area_fun = 1.0_8
      else
         area_fun = 0.0_8
      end if
    end function area_fun
  end subroutine avg_temp

  subroutine apply_initial_conditions_case
    implicit none
    integer :: d, k, l

    do l = level_end, level_start, -1
       call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       do k = 1, zmax
          call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do

       do k = 0, zlevels
          call apply_onescale (init_eddy, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do

       do k = 1, zlevels
          call apply_onescale (init_tke, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
    end do

    do l = level_end, level_start, -1
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_sol,  l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update means, bathymetry and penalization mask
    use wavelet_mod
    implicit none
    integer :: d, k, l, p

    if (istep /= 0) then
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             call apply_onescale_to_patch (set_bathymetry, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
             do k = 1, zlevels
                call apply_onescale_to_patch (set_penal, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
             end do
          end do
       end do
       call barrier

       do k = 1, zlevels
          do d = 1, size(grid)
             do p = n_patch_old(d)+1, grid(d)%patch%length
                call apply_onescale_to_patch (init_mean, grid(d), p-1, k, -BDRY_THICKNESS, BDRY_THICKNESS)
             end do
          end do
       end do
    else ! need to set values over entire grid on restart
       do l = level_start, level_end
          call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          do k = 1, zmax
             call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do

       do l = level_start, level_end
          do k = 1, zmax
             call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    end if
  end subroutine update_case
  
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Initial perturbation to mean for an entire vertical column 
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_k, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    eta = init_free_surface (dom%node%elts(id_i))
    z_s = topography%data(d)%elts(id_i)

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (d, id_i, k)
       z_k = interp (z(k-1), z(k))

       if (k == zlevels) then
          sol(S_MASS,zlevels)%data(d)%elts(id_i) = rho * eta
       else
          sol(S_MASS,k)%data(d)%elts(id_i) = 0.0_8
       end if
       sol(S_TEMP,k)%data(d)%elts(id_i)                      = rho * dz(k) * buoyancy_init (z_k)
       sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end do

    if (mode_split) then
       phi = phi_node (d, id_i, zlevels)
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i) = phi * eta ! free surface perturbation
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0.0_8
    end if
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values for an entire vertical column
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer                       :: d, id, id_i, k 
    real(8)                       :: eta, phi, rho, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    d    = dom%id+1
    id   = idx (i, j, offs, dims) 
    id_i = id + 1

    eta = 0d0
    z_s = topography%data(d)%elts(id_i)

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       rho = porous_density (d, id_i, k)

       sol_mean(S_MASS,k)%data(d)%elts(id_i)                      = rho * dz(k)
       sol_mean(S_TEMP,k)%data(d)%elts(id_i)                      = 0.0_8
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8
    end do

    if (mode_split) then
       sol_mean(S_MASS,zlevels+1)%data(d)%elts(id_i) = 0.0_8
       sol_mean(S_TEMP,zlevels+1)%data(d)%elts(id_i) = 0.0_8
    end if
  end subroutine init_mean

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid 
    ! (a_vert, b_vert not used if sigma_z = .true.)
    implicit none
    integer               :: k
    real(8)               :: z
    real(8), dimension(6) :: p

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    b_vert(0) = 1.0_8 ; b_vert(zlevels) = 0.0_8
    if (trim (coords) == "uniform") then
       do k = 2, zlevels-1
          b_vert(k) = 1.0_8 - dble(k)/dble(zlevels)
       end do
    elseif (trim(coords) == "croco") then
       p = (/ -5.7831,  18.9754, -24.6521,  16.1698, -5.7092, 0.9972 /)
       do k = 1, zlevels-1
          z = dble(k)/dble(zlevels)
          b_vert(k) = p(1)*z**5 + p(2)*z**4 + p(3)*z**3 + p(4)*z**2 + p(5)*z + p(6)
       end do
    end if
    a_vert = 1.0_8 - b_vert

    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

  subroutine init_eddy (dom, i, j, zlev, offs, dims)
    ! Initialize eddy diffusivity and eddy viscosity
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer  :: d, id

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    Kt(zlev)%data(d)%elts(id) = Kt_0
    Kv(zlev)%data(d)%elts(id) = Kv_0
  end subroutine init_eddy

  subroutine init_tke (dom, i, j, zlev, offs, dims)
    ! Initialize TKE
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer  :: d, id

    d  = dom%id+1
    id = idx (i, j, offs, dims) + 1

    tke(zlev)%data(d)%elts(id) = e_min
  end subroutine init_tke

  real(8) function surf_geopot_case (d, id)
    ! Surface geopotential: postive if greater than mean seafloor                                                                                        
    implicit none
    integer :: d, id

    surf_geopot_case = 0.0_8
  end function surf_geopot_case

  real(8) function init_free_surface (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    init_free_surface = 0.0_8
  end function init_free_surface

  real(8) function buoyancy_init (z)
    ! Initial buoyancy at depth z
    ! buoyancy = (ref_density - density)/ref_density
    implicit none
    real(8) :: z

    real(8) :: rho

    rho = density_eos (S_ref, temp_init (z), z)
    buoyancy_init = (ref_density - rho) / ref_density
  end function buoyancy_init

  real(8) function temp_init (z)
    implicit none
    real(8) :: z

    temp_init = T_ref - N_0**2 * ref_density / (a_0 * grav_accel) * abs (z)
  end function temp_init

  subroutine print_density
    implicit none
    integer                       :: k
    real(8)                       :: bv, c_k, c1, drho, dzl, eta, rho, rho_above, z_k, z_s, z_above
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    eta = 0.0_8
    z_s = max_depth

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = a_vert * eta + b_vert * z_s
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    write (6,'(a)') " Layer        z         dz         T          rho "
    do k = 1, zlevels
       z_k = interp (z(k-1), z(k))
       rho = ref_density * (1.0_8 - buoyancy_init (z_k))
       write (6, '(2x, i4, 4x, 2(es9.2, 1x), es10.4, 1x, es13.7)') k, z_k, dz(k), temp_init (z_k), rho
    end do

    write (6,'(/,a)') " Interface     z"
    do k = 0, zlevels
       write (6, '(3x, i4, 5x, es9.2)') k, z(k)
    end do

    write (6,'(/,a)') " Interface      V        c1      CFL_c1"
    c1 = 0.0_8
    do k = 1, zlevels-1
       z_above = interp (z(k),   z(k+1))
       z_k     = interp (z(k-1), z(k))
       dzl    = z_above - z_k

       rho_above = ref_density * (1.0_8 - buoyancy_init (z_above))
       rho  = ref_density * (1.0_8 - buoyancy_init (z_k))
       drho = rho_above - rho

       bv = sqrt(- grav_accel * drho/dzl/rho)
       c_k = bv * abs(max_depth) / MATH_PI
       c1 = max (c1, c_k)

       write (6, '(3x, i4, 5x,3(es9.2,1x))') k, bv, c_k, c1*dt_init/dx_min
    end do
    write (6,'(/,a,es10.4)') "Maximum internal wave speed = ", c1
    write (6,'(A)') &
         '*********************************************************************&
         ************************************************************'
  end subroutine print_density

  subroutine set_thresholds_case
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                 :: k, v
    real(8), dimension(1:N_VARIABLE,1:zmax) :: threshold_new
    character(3), parameter                 :: order = "inf"

    if (default_thresholds) then ! initialize once
       threshold = threshold_def
    else
       call cal_lnorm (order)
       threshold_new = tol * lnorm

       ! Correct very small values
       do k = 1, zmax
          if (threshold_new(S_MASS,k) < threshold_def(S_MASS,k)/10) threshold_new(S_MASS,k) = threshold_def(S_MASS,k)
          if (threshold_new(S_TEMP,k) < threshold_def(S_TEMP,k)/10) threshold_new(S_TEMP,k) = threshold_def(S_TEMP,k)
          if (threshold_new(S_VELO,k) < threshold_def(S_VELO,k)/10) threshold_new(S_VELO,k) = threshold_def(S_VELO,k)
       end do

       if (istep >= 10) then
          threshold = 0.01*threshold_new + 0.99*threshold
       else
          threshold = threshold_new
       end if
    end if
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
    ! Set default thresholds based on dimensional scalings of norms
    implicit none
    integer                       :: k
    real(8)                       :: eta, z_s
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z
    type(Coord)                   :: x_i

    eta = 0.0_8
    z_s = max_depth

    if (sigma_z) then
       z = z_coords_case (eta, z_s)
    else
       z = eta * a_vert + z_s * b_vert
    end if
    dz = z(1:zlevels) - z(0:zlevels-1)

    do k = 1, zlevels
       lnorm(S_MASS,k) = ref_density * dz(k)
       lnorm(S_TEMP,k) = drho * dz(k)
       lnorm(S_VELO,k) = Udim
    end do

    if (mode_split) lnorm(:,zlevels+1) = lnorm(:,zlevels) ! not used

    threshold_def = tol * lnorm  
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    ! Initializes viscosity, time step and penalization parameter eta
    implicit none
    real(8) :: area, C, C_b, C_divu, C_mu, C_rotu, dlat, tau_b, tau_divu, tau_mu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = 0.891 * sqrt (4/sqrt(3.0_8) * area) ! edge length of average triangle

    area = 4*MATH_PI*radius**2/(20*4**min_level)
    dx_max = sqrt (4/sqrt(3.0_8) * area)


    C = 0.0_8 ! <= 1/2 if explicit
    C_rotu = C / 4**Laplace_order_init
    C_divu = C
    C_mu   = C
    C_b    = C

    ! Diffusion time scales
    tau_mu   = dt_cfl / C_mu
    tau_b    = dt_cfl / C_b
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr(S_MASS) = dx_min**(2*Laplace_order_init) / tau_mu
       visc_sclr(S_TEMP) = dx_min**(2*Laplace_order_init) / tau_b
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,4(a,es8.2),a,/)') &
            "dx_max  = ", dx_max/KM, " dx_min  = ", dx_min/KM, " [km] dt_cfl = ", dt_cfl, " [s] tau_mu = ", tau_mu/HOUR, " [h]"
       write (6,'(4(a,es8.2),/)') "C_mu = ", C_mu,  " C_b = ", C_mu, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2),/)') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
            " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
            " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity_case

  subroutine set_bathymetry (dom, i, j, zlev, offs, dims)
    ! Set bathymetry
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    call topo_tke1d (dom, i, j, zlev, offs, dims, 'bathymetry')
  end subroutine set_bathymetry

  subroutine set_penal (dom, i, j, zlev, offs, dims)
    ! Set penalization mask
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, id_i

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1 

    if (penalize) then
       call topo_tke1d (dom, i, j, zlev, offs, dims, "penalize")
    else
       penal_node(zlev)%data(d)%elts(id_i)                      = 0.0_8
       penal_edge(zlev)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0.0_8       
    end if
  end subroutine set_penal

  subroutine topo_tke1d (dom, i, j, zlev, offs, dims, itype)
    ! Not used
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    character(*)                   :: itype

    integer :: d, id, id_i

    d = dom%id + 1
    id = idx (i, j, offs, dims)
    id_i = id + 1 

    select case (itype)
    case ("bathymetry")
       topography%data(d)%elts(id_i) = max_depth + surf_geopot_case (d, id_i) / grav_accel
    case ("penalize") ! not used
    end select
  end subroutine topo_tke1d

  subroutine wind_stress (lon, lat, tau_zonal, tau_merid)
    implicit none
    real(8) :: lat, lon, tau_zonal, tau_merid

    tau_zonal = tau_0
    tau_merid = 0d0
  end subroutine wind_stress

  real(8) function tau_mag_case (p)
    ! Magnitude of wind stress at node p
    implicit none
    type(Coord) :: p

    real(8) :: lat, lon, tau_zonal, tau_merid

    call cart2sph (p, lon, lat)

    call wind_stress (lon, lat, tau_zonal, tau_merid)

    tau_mag_case = sqrt (tau_zonal**2 + tau_merid**2)
  end function tau_mag_case

  subroutine set_save_level_case
    ! Save top layer
    implicit none
    real(8) :: save_height

    save_height = 0.5 * (b_vert(save_zlev)+b_vert(save_zlev-1)) * max_depth

    if (rank==0) write (6,'(/,A,i2,A,es11.4,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_height, " [m])"
  end subroutine set_save_level_case

  function z_coords_case (eta_surf, z_s)
    ! Hybrid sigma-z vertical coordinates to minimize inclination of layers to geopotential
    ! near the free surface over strong bathymetry gradients.
    ! Reference: similar to Shchepetkin and McWilliams (JCP vol 228, 8985-9000, 2009)
    !
    ! Sets the a_vert parameter that depends on eta_surf (but not b_vert).
    implicit none
    real(8)                       :: eta_surf, z_s ! free surface and bathymetry
    real(8), dimension(0:zlevels) :: z_coords_case

    integer                       :: k
    real(8)                       :: cff, cff1, cff2, hc, z_0
    real(8), parameter            :: theta_b = 0d0, theta_s = 7d0
    real(8), dimension(0:zlevels) :: Cs, sc

    hc = abs(min_depth)

    cff1 = 1.0_8 / sinh (theta_s)
    cff2 = 0.5d0 / tanh (0.50 * theta_s)

    sc(0) = -1.0_8
    Cs(0) = -1.0_8
    cff = 1d0 / dble(zlevels)
    do k = 1, zlevels
       sc(k) = cff * dble (k - zlevels)
       Cs(k) = (1.0_8 - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5d0)) - 0.5d0)
    end do

    z_coords_case(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1.0_8 - z_0 / z_s
       z_coords_case(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords_case

  subroutine dump_case (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump_case

  subroutine load_case (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load_case

  subroutine cal_r_max
    ! Calculates minimum relative mass and checks diffusion stability limits
    use mpi
    implicit none
    integer :: ierror, k, l

    r_max_loc = 1d-16
    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (cal_rmax_loc, l, k, 0, 0)
       end do
    end do

    call MPI_Allreduce (r_max_loc, r_max, 1, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_COMM_WORLD, ierror)
  end subroutine cal_r_max

  subroutine cal_rmax_loc (dom, i, j, zlev, offs, dims)
    ! Calculates minimum mass and diffusion stability limits
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idE, idN, idNE, idS, idSW, idW, k
    real(8) :: dz0, dz_e, r_loc

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    d    = dom%id + 1

    if (dom%mask_n%elts(id+1) >= ADJZONE) then
       dz0  = (sol(S_MASS,zlev)%data(d)%elts(id+1) + sol_mean(S_MASS,zlev)%data(d)%elts(id+1)) &
            / porous_density (d, id+1, zlev)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idE+1)) &
            / porous_density (d, idE+1, zlev)
       r_loc = abs (dz0 - dz_e) / (dz0 + dz_e)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idNE+1) + sol_mean(S_MASS,zlev)%data(d)%elts(idNE+1)) &
            / porous_density (d, idNE+1, zlev)
       r_max_loc = max (r_max_loc, r_loc)

       dz_e = (sol(S_MASS,zlev)%data(d)%elts(idN+1)  + sol_mean(S_MASS,zlev)%data(d)%elts(idN+1)) &
            / porous_density(d, idN+1, zlev)
       r_max_loc = max (r_max_loc, r_loc)
    end if
  end subroutine cal_rmax_loc

  real(8) function bottom_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Bottom boundary flux for vertical diffusion of buoyancy 
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id_i

    d = dom%id + 1
    id_i = idx (i, j, offs, dims) + 1

    bottom_buoy_flux_case =  Kt(0)%data(d)%elts(id_i) * N_0**2 / grav_accel
  end function bottom_buoy_flux_case

  real(8) function top_buoy_flux_case (dom, i, j, z_null, offs, dims)
    ! Top boundary flux for vertical diffusion of buoyancy
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, z_null
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    if (kato) then
       top_buoy_flux_case = 0.0_8
    else
       top_buoy_flux_case = Q_0 / (ref_density*c_p)  * a_0/ref_density 
    end if
  end function top_buoy_flux_case

  function wind_flux_case (dom, i, j, zlev, offs, dims)
    ! Wind stress velocity source term evaluated at edges (top boundary condition for vertical diffusion of velocity)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
    real(8), dimension(1:EDGE)     :: wind_flux_case

    integer                     :: d, id, idE, idN, idNE
    real(8)                     :: rho
    real(8), dimension(1:EDGE)  :: tau_wind

    id = idx (i, j, offs, dims)

    d = dom%id + 1
    idE  = idx (i+1, j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    tau_wind(RT+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idE+1))
    tau_wind(DG+1) = proj_vel (wind_stress, dom%node%elts(idNE+1), dom%node%elts(id+1))
    tau_wind(UP+1) = proj_vel (wind_stress, dom%node%elts(id+1),   dom%node%elts(idN+1))

    rho = porous_density (d, id+1, zlevels)

    wind_flux_case = tau_wind / rho
  end function wind_flux_case
end module test_case_mod
