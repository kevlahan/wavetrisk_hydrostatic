program flat_projection_data
  ! Post-processing of checkpoint data to calculate flat projection
  use main_mod
  use test_case_mod
  use io_mod
  use projection_mod
  implicit none
  
  integer                                :: k, l, nt, Ncumul, p, d
  integer, parameter                     :: nvar_save = 6, nvar_drake = 12, nvar_1layer = 5
  real(8)                                :: area1, area2
  real(8), dimension(:),     allocatable :: eta_lat, eta_lon
  real(8), dimension(:,:),   allocatable :: drake_ke, drake_enstrophy
  real(8), dimension(:,:,:), allocatable :: field2d_save, lat_slice, lon_slice, zonal_av, zcoord_lat, zcoord_lon, field2d_simplephys
  
  character(2)                           :: var_file
  character(8)                           :: itype
  character(130)                         :: command

  logical, parameter                     :: welford = .true. ! use Welford's one-pass algorithm or naive two-pass algorithm
  
  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters
  
  select case (test_case)
  case ("DCMIP2012c4")
     compressible   = .true.                      ! Compressible equations
     mean_split     = .false.

     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
     omega          = 7.29212d-5                  ! Earth’s angular velocity in radians per second
     p_0            = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = p_0                         ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

     u_0            = 35.0_8                      ! maximum velocity of zonal wind
     eta_0          = 0.252_8                     ! value of eta at reference level (level of the jet)
  case ("DCMIP2008c5")
     compressible   = .true.                      ! Compressible equations
     mean_split     = .false.

     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
     omega          = 7.29211d-5                  ! Earth’s angular velocity in radians per second
     p_0            = 100145.6_8                  ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = 930.0d2                     ! reference surface pressure
     R_d            = 287.04_8                    ! ideal gas constant for dry air in joules per kilogram Kelvin
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

     d2             = 1.5d6**2                    ! square of half width of Gaussian mountain profile in meters
     h_0            = 2.0d3                       ! mountain height in meters
     lon_c          = MATH_PI/2                   ! longitude location of mountain
     lat_c          = MATH_PI/6                   ! latitude location of mountain
  case ("Held_Suarez")
     compressible   = .true.                      ! Compressible equations
     mean_split     = .false.

     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.8_8                       ! gravitational acceleration in meters per second squared
     omega          = 7.292d-5                    ! Earth’s angular velocity in radians per second
     p_0            = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = p_0                         ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     gamma          = c_p/c_v                     ! heat capacity ratio
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

     u_0            = 35.0_8                      ! maximum velocity of zonal wind
     eta_0          = 0.252_8                     ! value of eta at reference level (level of the jet)
  case ("drake")
     scale          = 6                                  ! scale factor for small planet (1/6 Earth radius)
     radius         = 6371.229/scale   * KM              ! mean radius of the small planet
     grav_accel     = 9.80616          * METRE/SECOND**2 ! gravitational acceleration 
     omega          = 7.29211d-5/scale * RAD/SECOND      ! angular velocity (scaled for small planet to keep beta constant)
     p_top          = 0.0_8            * hPa             ! pressure at free surface
     ref_density    = 1028             * KG/METRE**3     ! reference density at depth (seawater)

     max_depth   = -4000 * METRE
     halocline   = -1000 * METRE                      ! location of top (less dense) layer in two layer case
     mixed_layer = -1000 * METRE                      ! location of layer forced by surface wind stress
     drho        =    -8 * KG/METRE**3                ! density perturbation at free surface (density of top layer is rho0 + drho/2)
     density_drake = (/ ref_density, ref_density + drho/2 /)    ! densities in each layer
     height      = (/ abs(max_depth - halocline), abs(halocline) /) ! depths of each layer
     npts_penal  = 4

     mode_split     = .true.                             ! split barotropic mode if true
     mean_split     = .true.
     compressible   = .false.                            ! always run with incompressible equations
     penalize       = .true.                             ! penalize land regions
  case ("seamount") 
     scale          = 41.75d0                             
     radius         = 6371.229/scale   * KM              
     grav_accel     = 9.80616          * METRE/SECOND**2 
     omega          = 7.29211d-5       * RAD/SECOND      
     p_top          = 0.0_8            * hPa             
     ref_density    = 1000             * KG/METRE**3     

     max_depth   = -5000 * METRE
     drho        =    -3 * KG/METRE**3               

     mode_split     = .true.
     mean_split     = .true.
     compressible   = .false.                            
     penalize       = .false.

     lat_c          = 43.29 * DEG                               ! latitude of seamount
     lon_c          =     0 * DEG                               ! longitude
     h0             =  4500 * METRE                             ! height of seamount
     width          =    40 * KM                                ! radius of seamount
     delta          =   500 * METRE                             ! vertical decay of density

     sigma_z        = .false.
     coords         = "chebyshev"
     stratification = "exponential"
  case ("upwelling")
     radius         = 240      * KM             
     grav_accel     = 9.80616  * METRE/SECOND**2 
     omega          = 6d-5     * RAD/SECOND      
     p_top          = 0.0_8    * hPa             
     ref_density    = 1027     * KG/METRE**3     

     max_depth   =  -150 * METRE
     min_depth   =   -25 * METRE
     Tcline      =  - 50 * METRE
     drho        =    -2.5 * KG/METRE**3

     sigma_z        = .true.                       ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     coords         = "croco"                      ! grid type for pure sigma grid ("croco" or "uniform")
     mode_split     = .true.
     mean_split     = .true.
     compressible   = .false.                            
     penalize       = .true.
     vert_diffuse   = .true.

     a_0            = 0.28 / CELSIUS
     b_0            = 0.0_8
     mu_1           = 0.0_8
     mu_2           = 0.0_8
     T_ref          = 14   * CELSIUS

     alpha          = 1d-6    ! porosity
     npts_penal     = 2.5

     width          = 80 * KM                           ! width of channel in km
     lat_width      = (width/radius)/DEG
     lat_c          = 45                                ! centre of zonal channel (in degrees)
  case ("jet")
     soufflet           = .false.                          ! set radius to exactly match Soufflet domain
     lat_c              = 30d0                            ! centre of zonal channel (in degrees)
     
     if (soufflet) then
        width           = 2000d0 * KM
        beta            = 1.6d-11 / (METRE * SECOND)      ! beta parameter
        f0              = 1d-4    / SECOND                ! Coriolis parameter
        omega           = f0 / (2d0*sin(lat_c * DEG))     ! planet rotation
        radius          = f0 / (beta * tan (lat_c * DEG)) ! planet radius to exactly match Soufflet beta plane
        L_jet              = 0.8d0 * width                   ! width of jet transition region
     else
        lat_c           = 30d0                            ! centre of zonal channel (in degrees)
        radius          = 1000d0 * KM                     ! meridional width of zonal channel
        width           = radius                          ! zonal channel width
        L_jet           = 0.4d0 * width                   ! width of jet transition region
        f0              = 1d-4  / SECOND                  ! Coriolis parameter
        omega           = f0 / (2d0*sin(lat_c*DEG))       ! planet rotation
     end if

     grav_accel         = 9.80616d0    * METRE/SECOND**2  ! gravitational acceleration 
     ref_density        = 1027.75d0    * KG/METRE**3      ! reference density at depth (maximum density)

     sigma_z            = .true.                        ! use sigma-z Schepetkin/CROCO type vertical coordinates (pure sigma grid if false)
     coords             = "croco"                       ! grid type for pure sigma grid ("croco" or "uniform")
     max_depth          = -4000d0 * METRE               ! total depth
     min_depth          = -4000d0 * METRE               ! minimum depth
     Tcline             =  -100d0 * METRE               ! position of thermocline

     lat_width          = (width/radius)/DEG              ! width of zonal channel (in degrees)
     
     mode_split     = .true.
     mean_split     = .true.
     compressible   = .false.                            
     penalize       = .false.
     vert_diffuse   = .true.

     a_0            = 0.28 / CELSIUS
     b_0            = 0.0_8
     mu_1           = 0.0_8
     mu_2           = 0.0_8
     T_ref          = 14   * CELSIUS

     alpha          = 1d-2    ! porosity
     npts_penal     = 4.5d0
   case ("Simple_Physics")
      radius         = 6400      * KM                 ! mean radius of the Earth
      grav_accel     = 9.8       * METRE/SECOND**2    ! gravitational acceleration
      p_0            = 1000      * hPa                ! reference pressure (mean surface pressure) in Pascals
      p_top          = 0.01       * Pa                 ! pressure at the top in Pascals
      c_p            = 1004.0_8  * JOULE/(KG*KELVIN)  ! specific heat at constant pressure in joules per kilogram Kelvin
      R_d            = 287 * JOULE/(KG*KELVIN)         ! ideal gas constant for dry air in joules per kilogram Kelvin
      c_v            = c_p - R_d * JOULE/(KG*KELVIN)  ! specific heat at constant volume c_v = c_p - R_d
      ref_density    = 1.204     * KM                  ! Reference density (km/m^3)

      kappa          = R_d/c_p                    ! kappa
      gamma          = c_p/c_v  
      zmin = -10
      ref_surf_press = p_0
      climatology = .true.
   case default
     if (rank == 0) write (6,'(a)') "Case not supported ... aborting"
     call abort
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  resume = mean_beg

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize (run_id)

  ! Initialize statistics
  call initialize_projection (N)
  if (trim (test_case) == "drake") then
     call initialize_stat_drake
  elseif (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling" .or. trim (test_case) == "jet") then
     call initialize_stat_vertical
  else
     if (trim(test_case) == "Simple_Physics" .and. climatology) call init_physics_climatology
     call initialize_stat
  end if
  
  Nt = 0
  do cp_idx = mean_beg, mean_end
     Nt = Nt + 1
     resume = NONE
     call restart (run_id)

     ! Set means
     do l = level_start, level_end
        do k = 1, zmax
           call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
        end do
     end do

     if  (trim (test_case) == "drake") then
        if (zlevels == 2) then
           drake_ke(Nt,1) = time
           drake_ke(Nt,2:5) = energy_drake ('adaptive')
           drake_enstrophy(Nt,1) = time
           drake_enstrophy(Nt,2:6) = pot_enstrophy_drake ('adaptive')
           if (cp_idx == cp_2d) call latlon_drake
        elseif (zlevels == 1) then
           drake_ke(Nt,1) = time
           drake_ke(Nt,2) = energy_1layer ('adaptive')
           drake_enstrophy(Nt,1) = time
           drake_enstrophy(Nt,2) = pot_enstrophy_1layer ('adaptive')
           if (cp_idx == cp_2d) call latlon_1layer
        end if
     elseif  (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling" .or. trim (test_case) == "jet") then
        call vertical_slice
        if (rank == 0) call write_slice
     elseif (trim (test_case) == "Simple_Physics") then
         if (climatology) then
            call cal_surf_press(sol(1:N_VARIABLE,1:zmax))
            ! Add each temp & KE for each checkpoint for the climatology
            !Update the boundary for the velocities
            sol%bdry_uptodate = .false.
            call update_array_bdry (sol, NONE, 26)
            do k = 1, zlevels
               do d = 1, size(grid)
                  temp   => sol(S_TEMP,k)%data(d)%elts
                  temp1  => simple_phys_temp(k)%data(d)%elts
                  mass   => sol(S_MASS,k)%data(d)%elts
                  mean_m => sol_mean(S_MASS,k)%data(d)%elts
                  mean_t => sol_mean(S_TEMP,k)%data(d)%elts
                  velo   => sol(S_VELO,k)%data(d)%elts
                  velo_2d  => simple_phys_vels(k)%data(d)%elts
                  velo1  => grid(d)%u_zonal%elts
                  velo2  => grid(d)%v_merid%elts
                  do p = 3, grid(d)%patch%length
                     call apply_onescale_to_patch (climatology_add_temp, grid(d), p-1, k, 0, 1)
                     call apply_onescale_to_patch(climatology_add_velocities, grid(d), p-1, k, 0, 0)
                     call apply_onescale_to_patch (climatology_add_KE, grid(d), p-1, k, 0, 1)
                     if (cp_idx==cp_2d) then
                         call apply_onescale_to_patch (climatology_temp_mean, grid(d), p-1, k, 0, 1)
                         call apply_onescale_to_patch (climatology_velocity_mean, grid(d), p-1, k, 0, 0)
                         call apply_onescale_to_patch (climatology_KE_mean, grid(d), p-1, k, 0, 1)
                     end if
                  end do
                  nullify(temp, temp1, mass, mean_m, mean_t, velo, velo_2d, velo1, velo2)
               end do
            end do
         end if 
         if (welford) then
            call cal_zonal_av
         else
            call cal_zonal_average
         end if
         if (cp_idx == cp_2d) call latlon
     else
        if (welford) then
           call cal_zonal_av
        else
           call cal_zonal_average
        end if
        if (cp_idx == cp_2d) call latlon
     end if
  end do

  if (trim (test_case) == "drake") then
     if (rank==0) then
        if (zlevels == 2) then
           call write_out_drake
        elseif (zlevels == 1) then
           call write_out_1layer
        end if
     end if
  elseif (.not. (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling" .or. trim (test_case) == "jet")) then
     if (.not. welford) then
        zonal_av(:,:,1)   = zonal_av(:,:,1)   / Ncumul
        zonal_av(:,:,3:5) = zonal_av(:,:,3:5) / Ncumul
        call barrier

        do cp_idx = mean_beg, mean_end
           resume = NONE
           call restart (run_id)
           call cal_variance
        end do
        call barrier
     end if

     ! Finish covariance calculations
     zonal_av(:,:,2)   = zonal_av(:,:,2)   / (Ncumul-1)
     zonal_av(:,:,6:9) = zonal_av(:,:,6:9) / (Ncumul-1)

     if (rank==0) call write_out
  end if
  call finalize
contains
  subroutine cal_zonal_av
    ! Zonal average means and covariances over all checkpoints using stable online algorithm
    ! Uses Welford's stable onlne algorithm
    use domain_mod
    implicit none
    
    integer                                       :: d, ix, j, k
    real(8), dimension (Ny(1):Ny(2))              :: Tprime, Uprime, Vprime, Tprime_new, Uprime_new, Vprime_new
    real(8), dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level_save
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 41)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (exner_fun(k), level_save, 1d0)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       call project_array_onto_plane ("u_zonal", level_save, 0d0)
       Uproj = field2d
       call project_array_onto_plane ("v_merid", level_save, 0d0)
       Vproj = field2d

       ! Update means and covariances
       do ix = Nx(1), Nx(2)
          Ncumul = (Nt-1)*(Nx(2)-Nx(1)+1) + ix-Nx(1)+1
          
          Tprime = Tproj(ix,:) - zonal_av(k,:,1)
          Uprime = Uproj(ix,:) - zonal_av(k,:,3)
          Vprime = Vproj(ix,:) - zonal_av(k,:,4)

          ! Mean values
          zonal_av(k,:,1) = zonal_av(k,:,1) + Tprime/Ncumul
          zonal_av(k,:,3) = zonal_av(k,:,3) + Uprime/Ncumul
          zonal_av(k,:,4) = zonal_av(k,:,4) + Vprime/Ncumul
          zonal_av(k,:,5) = zonal_av(k,:,5) +  0.5 * (Uprime**2 + Vprime**2)/Ncumul

          Tprime_new = Tproj(ix,:) - zonal_av(k,:,1)
          Uprime_new = Uproj(ix,:) - zonal_av(k,:,3)
          Vprime_new = Vproj(ix,:) - zonal_av(k,:,4)

          ! Temperature variance
          zonal_av(k,:,2) = zonal_av(k,:,2) + Tprime * Tprime_new

          ! Eddy momentum flux (covariance)
          zonal_av(k,:,6) = zonal_av(k,:,6) + Uprime * Vprime_new
       
          ! Zonal velocity variance
          zonal_av(k,:,7) = zonal_av(k,:,7) + Uprime * Uprime_new

          ! Meridional velocity variance
          zonal_av(k,:,8) = zonal_av(k,:,8) + Vprime * Vprime_new

          ! Eddy heat flux (covariance)
          zonal_av(k,:,9) = zonal_av(k,:,9) + Vprime * Tprime_new
       end do
    end do
  end subroutine cal_zonal_av

  subroutine cal_zonal_average
    ! Zonal average means over all checkpoints
    use domain_mod
    implicit none
    
    integer                                       :: d, ix, j, k
    real(8), dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 42)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       call project_array_onto_plane ("u_zonal", level_save, 0d0)
       Uproj = field2d
       call project_array_onto_plane ("v_merid", level_save, 0d0)
       Vproj = field2d

       ! Update means
       do ix = Nx(1), Nx(2)
          Ncumul = (Nt-1)*(Nx(2)-Nx(1)+1) + ix-Nx(1)+1

          zonal_av(k,:,1) = zonal_av(k,:,1) + Tproj(ix,:)
          zonal_av(k,:,3) = zonal_av(k,:,3) + Uproj(ix,:)
          zonal_av(k,:,4) = zonal_av(k,:,4) + Vproj(ix,:)
          zonal_av(k,:,5) = zonal_av(k,:,5) +  0.5 * (Uproj(ix,:)**2 + Vproj(ix,:)**2)
       end do
    end do
  end subroutine cal_zonal_average

  subroutine cal_variance
    ! Zonal average variances over all checkpoints 
    use domain_mod
    implicit none
    
    integer                                       :: d, ix, j, k
    real(8), dimension (Ny(1):Ny(2))              :: Tprime, Uprime, Vprime
    real(8), dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 43)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo => sol(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo)
       end do
       call project_array_onto_plane ("u_zonal", level_save, 0d0)
       Uproj = field2d
       call project_array_onto_plane ("v_merid", level_save, 0d0)
       Vproj = field2d

       ! Update covariances
       do ix = Nx(1), Nx(2)          
          Tprime = Tproj(ix,:) - zonal_av(k,:,1)
          Uprime = Uproj(ix,:) - zonal_av(k,:,3)
          Vprime = Vproj(ix,:) - zonal_av(k,:,4)

          ! Temperature variance
          zonal_av(k,:,2) = zonal_av(k,:,2) + Tprime**2

          ! Eddy momentum flux (covariance)
          zonal_av(k,:,6) = zonal_av(k,:,6) + Uprime * Vprime
       
          ! Zonal velocity variance
          zonal_av(k,:,7) = zonal_av(k,:,7) + Uprime**2

          ! Meridional velocity variance
          zonal_av(k,:,8) = zonal_av(k,:,8) + Vprime**2

          ! Eddy heat flux (covariance)
          zonal_av(k,:,9) = zonal_av(k,:,9) + Vprime * Tprime
       end do
    end do
  end subroutine cal_variance

  subroutine latlon
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels)
    use domain_mod
    implicit none
    integer :: d, j, k

    if (rank == 0) write (6,'(A,i6)') "Saving latitude-longitude projection of checkpoint file = ", cp_2d

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))

    ! Remap to pressure_save vertical levels for saving data
    sol_save = sol(:,1:save_levels)
    call apply_onescale (interp_save, level_save, z_null, -1, 2)

    ! Calculate temperature at all vertical levels (saved in exner_fun) and temperature at interpolated saved vertical levels
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 44)

    ! Latitude-longitude projections
    do k = 1, save_levels
       ! Temperature
       call project_field_onto_plane (trend(1,k), level_save, 0.0_8)
       field2d_save(:,:,1+k-1) = field2d

       ! Calculate zonal and meridional velocities and vorticity
       do d = 1, size(grid)
          velo  => sol_save(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,         grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (velo, velo1, velo2, vort)
       end do

       ! Zonal velocity
       call project_array_onto_plane ("u_zonal", level_save, 0d0)
       field2d_save(:,:,2+k-1) = field2d
       
       ! Meridional velocity
       call project_array_onto_plane ("v_merid", level_save, 0d0)
       field2d_save(:,:,3+k-1) = field2d

       ! Geopotential
       call apply_onescale (cal_geopot, level_save, z_null, 0, 1)
       call project_array_onto_plane ("geopot", level_save, 1d0)
       field2d_save(:,:,4+k-1) = field2d

       ! Vorticity
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          nullify (vort)
       end do
       call project_array_onto_plane ("press_lower", level_save, 1d0)
       field2d_save(:,:,5+k-1) = field2d

       ! Surface pressure
       call project_array_onto_plane ("surf_press", level_save, 1d0)
       field2d_save(:,:,6+k-1) = field2d

       ! Save climatology for simple Physics
       if (trim(test_case)=="Simple_Physics" .and. climatology) then
         ! update the boundarys
         simple_phys_temp%bdry_uptodate = .false.
         call update_vector_bdry (simple_phys_temp, NONE, 44)
         simple_phys_zonal%bdry_uptodate = .false.
         call update_vector_bdry (simple_phys_zonal, NONE, 44)
         simple_phys_merid%bdry_uptodate = .false.
         call update_vector_bdry (simple_phys_merid, NONE, 44)

         ! save 2D projections
         call project_field_onto_plane(simple_phys_temp(k-1), level_save, 0.0_8)
         field2d_simplephys(:,:,1+k-1) = field2d
         call project_field_onto_plane(simple_phys_zonal(k-1), level_save, 0.0_8)
         field2d_simplephys(:,:,4+k-1) = field2d
         call project_field_onto_plane(simple_phys_merid(k-1), level_save, 0.0_8)
         field2d_simplephys(:,:,5+k-1) = field2d

         simple_phys_vels%bdry_uptodate= .false.
         call update_vector_bdry(simple_phys_vels,NONE,27)
         ! Calculate zonal and meridional velocity
         do d = 1, size(grid)
            temp  => simple_phys_temp(k-1)%data(d)%elts
            velo  => simple_phys_vels(k-1)%data(d)%elts
            velo1 => grid(d)%u_zonal%elts
            velo2 => grid(d)%v_merid%elts
            do j = 1, grid(d)%lev(level_save)%length
               call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
            end do
            nullify (temp,velo, velo1, velo2)
         end do
         
         ! Zonal Vel
         call project_array_onto_plane ("u_zonal", level_save, 0d0)
         field2d_simplephys(:,:,2+k-1) = field2d
         
         ! Meridional Vel
         call project_array_onto_plane ("v_merid", level_save, 0d0)
         field2d_simplephys(:,:,3+k-1) = field2d

       end if
    end do
  end subroutine latlon

  subroutine vertical_slice
    ! Save vertical slice along given latitude and longitude for incompressible test cases
    use barotropic_2d_mod
    implicit none
    integer                         :: d, i, j, k, l, idx_lat, idx_lon
    real(8)                         :: dz, eta, z_s
    real(8), dimension(0:zlevels)   :: z
    real(8), dimension(Ny(1):Ny(2)) :: bathy_lat
    real(8), dimension(Nx(1):Nx(2)) :: bathy_lon

    ! Find indices of latitude and longitude slices
    idx_lat = minloc (abs(lat-lat_val), DIM=1) + Ny(1) - 1
    idx_lon = minloc (abs(lon-lon_val), DIM=1) + Nx(1) - 1
    
    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_save
    call fill_up_grid_and_IWT (l)
    
    call apply_onescale (set_bathymetry, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    do k = 1, zmax
       call apply_onescale (set_penal, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do

    do k = 1, zmax
       call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do

    do k = 1, zlevels
       do d = 1, size(grid)
          mean_m => sol_mean(S_MASS,k)%data(d)%elts
          mean_t => sol_mean(S_TEMP,k)%data(d)%elts
          mass   => sol(S_MASS,k)%data(d)%elts
          temp   => sol(S_TEMP,k)%data(d)%elts
          scalar => sol(S_TEMP,zlevels+1)%data(d)%elts
          velo   => sol(S_VELO,k)%data(d)%elts
          divu   => grid(d)%divu%elts
          velo1  => grid(d)%u_zonal%elts
          velo2  => grid(d)%v_merid%elts
          vort   => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_density,      grid(d), grid(d)%lev(l)%elts(j), k,       0, 1)
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,         grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (mass, mean_m, mean_t, scalar, temp, velo, divu, velo1, velo2, vort)
       end do

       ! Velocities
       call project_array_onto_plane ("u_zonal", l, 0d0)
       if (zonal) then
          lat_slice(:,k,1) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,1) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,1) = field2d(:,idx_lat)

       call project_array_onto_plane ("v_merid", l, 0d0)
       if (zonal) then
          lat_slice(:,k,2) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,2) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,2) = field2d(:,idx_lat)

       ! Temperature
       call project_field_onto_plane (sol(S_TEMP,zlevels+1), l, 0d0)
        if (zonal) then
          lat_slice(:,k,3) = sum(field2d,1) / size(field2d,1) 
       else
          lat_slice(:,k,3) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,3) = field2d(:,idx_lat)

       ! Vorticity
       call project_array_onto_plane ("press_lower", l, 1d0)
       if (zonal) then
          lat_slice(:,k,4) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,4) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,4) = field2d(:,idx_lat)
    end do

    ! Compute and project vertical velocity
    call vertical_velocity 
    do k = 1, zlevels
       call project_field_onto_plane (trend(S_TEMP,k), l, 0.0_8)
       if (zonal) then
          lat_slice(:,k,5) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,5) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,5) = field2d(:,idx_lat)
    end do
    
    ! Set free surface
    call project_field_onto_plane (sol(S_MASS,zlevels+1), l, 0.0_8)
    eta_lat = field2d(idx_lon,:)
    eta_lon = field2d(:,idx_lat)

    ! Set bathymetry
    call project_array_onto_plane ("topo", l, 0d0)
    bathy_lat = field2d(idx_lon,:)
    bathy_lon = field2d(:,idx_lat)
    
    do i = Ny(1), Ny(2)
       z_s = bathy_lat(i)
       eta = eta_lat(i)
       if (sigma_z) then
          z = z_coords (eta_lat(i), z_s)
       else
          z = a_vert * eta + b_vert * z_s
       end if
    
       do k = 1, zlevels
          zcoord_lat(i,k,1) = z(k-1)
          zcoord_lat(i,k,2) = z(k)  
       end do
    end do

    do i = Nx(1), Nx(2)
       z_s = bathy_lon(i)
       eta = eta_lon(i)
       if (sigma_z) then
          z = z_coords (eta, z_s)
       else
          z = a_vert * eta + b_vert * z_s
       end if
       
       do k = 1, zlevels
          zcoord_lon(i,k,1) = z(k-1)
          zcoord_lon(i,k,2) = z(k)  
       end do
    end do
  end subroutine vertical_slice

  function energy_drake (itype)
    ! Calculates baroclinic and bartoropic energies for two layer case
    use io_mod
    implicit none
    real(8), dimension(4) :: energy_drake
    character(*)          :: itype
    
    energy_drake(1) = integrate_hex (layer1_ke,     z_null)
    energy_drake(2) = integrate_hex (layer2_ke,     z_null)
    energy_drake(3) = integrate_hex (barotropic_ke, z_null) 
    energy_drake(4) = energy_drake(1) + energy_drake(2) - energy_drake(3) ! baroclinic energy
    
    energy_drake = energy_drake / (4*MATH_PI*radius**2*ref_density)
  end function energy_drake
  
  real(8) function energy_1layer (itype)
    ! Calculates baroclinic and bartoropic energies for two layer case
    use io_mod
    implicit none
    character(*) :: itype
    
    energy_1layer = integrate_hex (one_layer_ke, z_null)
    energy_1layer = energy_1layer / (4*MATH_PI*radius**2*ref_density)
  end function energy_1layer
  
  function pot_enstrophy_drake (itype)
    use io_mod
    implicit none
    real(8), dimension(5) :: pot_enstrophy_drake
    character(*)          :: itype
    
    integer :: d, j, k, l

    ! Potential enstrophy in each layer
    do k = 1, 2
       do l = level_start, level_end
          do d = 1, size(grid)
             velo  => sol(S_VELO,k)%data(d)%elts
             vort  => grid(d)%vort%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
             end do
             call apply_to_penta_d (post_vort, grid(d), l, z_null)
             nullify (velo, vort)
          end do

          ! Calculate vorticity at hexagon points (stored in press_lower)
          do d = 1, size(grid)
             vort => grid(d)%press_lower%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (vort)
          end do
       end do
       pot_enstrophy_drake(k) = integrate_hex (pot_enstrophy, k)
    end do

    ! Barotropic enstrophy
    do l = level_start, level_end
       do d = 1, size(grid)
          velo => sol(S_VELO,zlevels+1)%data(d)%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (barotropic_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          nullify (velo)
       end do
       do d = 1, size(grid)
          velo  => sol(S_VELO,zlevels+1)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       do d = 1, size(grid)
          velo  => sol(S_VELO,zlevels+1)%data(d)%elts
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, z_null)
          nullify (velo, vort)
       end do
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (vort)
       end do
    end do
    pot_enstrophy_drake(3) = integrate_hex (pot_enstrophy, 3)
    
    ! Baroclinic velocity and vorticity in each layer at hexagon points
    do k = 1, 2
       do l = level_start, level_end
          do d = 1, size(grid)
             velo1  => sol(S_VELO,zlevels+1)%data(d)%elts   ! barotropic velocity
             velo2  => sol(S_VELO,k)%data(d)%elts           ! velocity in current layer
             velo   => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in current layer
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (baroclinic_velocity, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
             end do
             nullify (velo, velo1, velo2)
          end do
          do d = 1, size(grid)
             velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in currrent layer
             velo1 => trend(S_VELO,1)%data(d)%elts
             velo2 => trend(S_VELO,2)%data(d)%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
             end do
             nullify (velo, velo1, velo2)
          end do
          do d = 1, size(grid)
             velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in current layer
             vort  => grid(d)%vort%elts
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 0)
             end do
             call apply_to_penta_d (post_vort, grid(d), l, z_null)
             nullify (velo, vort)
          end do
          do d = 1, size(grid)
             vort => grid(d)%press_lower%elts ! baroclinic vorticity in current layer in 
             do j = 1, grid(d)%lev(l)%length
                call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
             end do
             nullify (vort)
          end do
       end do
       pot_enstrophy_drake(k+3) = integrate_hex (pot_enstrophy, k)
    end do
  end function pot_enstrophy_drake

  real(8) function pot_enstrophy_1layer (itype)
    use io_mod
    implicit none
    character(*) :: itype

    integer :: d, j, l

    ! Potential enstrophy in each layer
    do l = level_start, level_end
       do d = 1, size(grid)
          velo  => sol(S_VELO,1)%data(d)%elts
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, z_null)
          nullify (velo, vort)
       end do

       ! Calculate vorticity at hexagon points (stored in press_lower)
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (vort)
       end do
    end do
    
    pot_enstrophy_1layer = integrate_hex (pot_enstrophy, 1)
  end function pot_enstrophy_1layer

  subroutine latlon_drake
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels)
    ! specific to 2 layer mode split drake test case
    use domain_mod
    use io_mod
    integer                              :: d, ibeg, ibeg_m, iend, iend_m, j, k, l
    real(8), dimension(:,:), allocatable :: dz

    if (rank == 0) write (6,'(A,i6)') "Saving latitude-longitude projection of checkpoint file = ", cp_2d

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_save
    call fill_up_grid_and_IWT (l)
    trend = sol

    ! Set mean on filled grid
    do k = 1, zmax
       call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do

    ! Compute barotropic velocity and vorticity at hexagon points
    
    ! Barotropic velocity
    do d = 1, size(grid)
       ibeg   = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       ibeg_m = (1+2*(POSIT(S_MASS)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
       iend   = sol(S_VELO,1)%data(d)%length
       iend_m = sol(S_MASS,1)%data(d)%length
       
       allocate (dz(ibeg_m:iend_m,1:2))
       dz(:,1) = sol_mean(S_MASS,1)%data(d)%elts(ibeg_m:iend_m) + sol(S_MASS,1)%data(d)%elts(ibeg_m:iend_m)
       dz(:,2) = sol_mean(S_MASS,2)%data(d)%elts(ibeg_m:iend_m) + sol(S_MASS,2)%data(d)%elts(ibeg_m:iend_m)

       sol(S_VELO,zlevels+1)%data(d)%elts(ibeg:iend:3) = &
            (dz(:,1) * sol(S_VELO,1)%data(d)%elts(ibeg:iend:3) + dz(:,2) * sol(S_VELO,2)%data(d)%elts(ibeg:iend:3)) &
            / sum (dz, dim=2)

       sol(S_VELO,zlevels+1)%data(d)%elts(ibeg+1:iend:3) = &
            (dz(:,1) * sol(S_VELO,1)%data(d)%elts(ibeg+1:iend:3) + dz(:,2) * sol(S_VELO,2)%data(d)%elts(ibeg+1:iend:3)) &
            / sum (dz, dim=2)

       sol(S_VELO,zlevels+1)%data(d)%elts(ibeg+2:iend:3) = &
            (dz(:,1) * sol(S_VELO,1)%data(d)%elts(ibeg+2:iend:3) + dz(:,2) * sol(S_VELO,2)%data(d)%elts(ibeg+2:iend:3))  &
            / sum (dz, dim=2)
       deallocate (dz)
    end do
    do d = 1, size(grid)
       velo  => sol(S_VELO,zlevels+1)%data(d)%elts
       velo1 => grid(d)%u_zonal%elts
       velo2 => grid(d)%v_merid%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (velo, velo1, velo2)
    end do
    
    ! Barotropic vorticity
    do d = 1, size(grid)
       velo  => sol(S_VELO,zlevels+1)%data(d)%elts
       vort  => grid(d)%vort%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       call apply_to_penta_d (post_vort, grid(d), l, z_null)
       nullify (velo, vort)
    end do
    do d = 1, size(grid)
       vort => grid(d)%press_lower%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       nullify (vort)
    end do

    ! Project barotropic velocity onto plane
    call project_array_onto_plane ("u_zonal", l, 0d0)
    field2d_save(:,:,1) = field2d
    call project_array_onto_plane ("v_merid", l, 0d0)
    field2d_save(:,:,2) = field2d
    
    ! Project barotropic vorticity onto plane
    call project_array_onto_plane ("press_lower", l, 1d0)
    field2d_save(:,:,3) = field2d

    ! Baroclinic velocity and baroclinic vorticity in each layer at hexagon points
    do k = 1, 2
       ! Baroclinic velocity 
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = sol(S_VELO,1)%data(d)%length
          trend(S_VELO,zlevels+1)%data(d)%elts(ibeg:iend) = &
               sol(S_VELO,k)%data(d)%elts(ibeg:iend) - sol(S_VELO,3)%data(d)%elts(ibeg:iend)
       end do
       ! Baroclinic velocity on hexagons
       do d = 1, size(grid)
          velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in currrent layer
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       
       ! Baroclinic vorticity 
       do d = 1, size(grid)
          velo  => trend(S_VELO,zlevels+1)%data(d)%elts ! baroclinic velocity in current layer
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), l, z_null)
          nullify (velo, vort)
       end do
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts ! baroclinic vorticity in current layer
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          nullify (vort)
       end do
       
       ! Project layer k baroclinic velocity onto plane
       call project_array_onto_plane ("u_zonal", l, 0d0)
       field2d_save(:,:,3*k+1) = field2d
       call project_array_onto_plane ("v_merid", l, 0d0)
       field2d_save(:,:,3*k+2) = field2d
       
       ! Project layer k barotropic vorticity onto plane
       call project_array_onto_plane ("press_lower", l, 1d0)
       field2d_save(:,:,3*k+3) = field2d
    end do

    ! Free surface
    call project_field_onto_plane (sol(S_MASS,zlevels+1), l, 1d0)
    field2d_save(:,:,10) = field2d

    ! Internal free surface
    call project_field_onto_plane (sol(S_MASS,1), l, 1d0)
    field2d_save(:,:,11) = field2d

    ! Penalization
    do k = 1, 2
       do d = 1, size(grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (set_penal, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
       end do
    end do
    call project_field_onto_plane (penal_node(1), l, 1d0)
    field2d_save(:,:,12) = field2d
  end subroutine latlon_drake

  subroutine latlon_1layer
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels)
    ! specific to one layer mode split drake test case
    use domain_mod
    use io_mod
    integer                              :: d, ibeg, ibeg_m, iend, iend_m, j, k, l
    real(8), dimension(:,:), allocatable :: dz

    if (rank == 0) write (6,'(A,i6)') "Saving latitude-longitude projection of checkpoint file = ", cp_2d

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_save
    call fill_up_grid_and_IWT (l)
    trend = sol

    ! Set mean on filled grid
    do k = 1, zmax
       call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
   
    ! Compute velocity at hexagon points
    do d = 1, size(grid)
       velo  => sol(S_VELO,1)%data(d)%elts
       velo1 => grid(d)%u_zonal%elts
       velo2 => grid(d)%v_merid%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (velo, velo1, velo2)
    end do
    
    ! Vorticity at hexagon points
    do d = 1, size(grid)
       velo  => sol(S_VELO,1)%data(d)%elts
       vort  => grid(d)%vort%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       call apply_to_penta_d (post_vort, grid(d), l, z_null)
       nullify (velo, vort)
    end do
    do d = 1, size(grid)
       vort => grid(d)%press_lower%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       nullify (vort)
    end do

    ! Project velocity onto plane
    call project_array_onto_plane ("u_zonal", l, 0d0)
    field2d_save(:,:,1) = field2d
    call project_array_onto_plane ("v_merid", l, 0d0)
    field2d_save(:,:,2) = field2d

    ! Project vorticity onto plane
    call project_array_onto_plane ("press_lower", l, 1d0)
    field2d_save(:,:,3) = field2d

    ! Free surface
    call project_field_onto_plane (sol(S_MASS,zlevels+1), l, 1d0)
    field2d_save(:,:,4) = field2d

    ! Penalization
    do k = 1, 2
       do d = 1, size(grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (set_penal, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
       end do
    end do
    call project_field_onto_plane (penal_node(1), l, 1d0)
    field2d_save(:,:,5) = field2d
  end subroutine latlon_1layer

  subroutine write_out
    ! Writes out results
    integer            :: i, k, v
    integer, parameter :: funit = 400

    ! 2d projections
    do v = 1, nvar_save*save_levels
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.4.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_save(:,i,v)
       end do
       close (funit)
    end do

    if (trim(test_case)=="Simple_Physics" .and. climatology) then
      do v = 1, 5*save_levels
         write (var_file, '(i2)') v+30
         open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
         do i = Ny(1), Ny(2)
            write (funit) field2d_simplephys(:,i,v)
         end do
         close (funit)
      end do
    end if

    ! Zonal average of solution over all vertical levels
    do v = 1, nvar_zonal
       write (var_file, '(i2)') v+10
       open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do k = zlevels,1,-1
          write (funit) zonal_av(k,:,v)
       end do
       close (funit)
    end do

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
    close (funit)

    ! Non-dimensional pressure based vertical coordinates p_k/p_s
    write (var_file, '(i2)') 22
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (0.5*((a_vert(k)+a_vert(k+1))/ref_surf_press + b_vert(k)+b_vert(k+1)), k = zlevels, 1, -1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.4.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (command)
  end subroutine write_out

  subroutine write_out_drake
    ! Writes out results
    implicit none
    integer            :: i, it, v
    integer, parameter :: funit = 400

    ! 2d projections
    do v = 1, 9
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.4.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_save(:,i,v)
       end do
       close (funit)
    end do

    do v = 10, nvar_drake
       write (var_file, '(i2)') v
       open (unit=funit, file=trim(run_id)//'.4.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_save(:,i,v)
       end do
       close (funit)
    end do

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
    close (funit)
    
    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.4.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (command)

    ! Write out kinetic energies
    open (unit=funit, file=trim(run_id)//'_kinetic_energy', form="FORMATTED", status="REPLACE")
    do it = 1, mean_end-mean_beg+1
       write (funit,'(5(es11.4,1x))') drake_ke(it,1)/DAY, drake_ke(it,2:5)
    end do
    close (funit)

    ! Write out potential enstrophy
    open (unit=funit, file=trim(run_id)//'_pot_enstrophy', form="FORMATTED", status="REPLACE")
    do it = 1, mean_end-mean_beg+1
       write (funit,'(6(es11.4,1x))') drake_enstrophy(it,1)/DAY, drake_enstrophy(it,2:6)
    end do
    close (funit)
  end subroutine write_out_drake

  subroutine write_out_1layer
    ! Writes out results
    implicit none
    integer            :: i, it, v
    integer, parameter :: funit = 400

    ! 2d projections
    do v = 1, nvar_1layer
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.4.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_save(:,i,v)
       end do
       close (funit)
    end do

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
    
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.4.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (command)

    ! Write out kinetic energies
    open (unit=funit, file=trim(run_id)//'_kinetic_energy', form="FORMATTED", status="REPLACE")
    do it = 1, mean_end-mean_beg+1
       write (funit,'(5(es11.4,1x))') drake_ke(it,1)/DAY, drake_ke(it,2)
    end do
    close (funit)
    deallocate (drake_ke)

    ! Write out potential enstrophy
    open (unit=funit, file=trim(run_id)//'_pot_enstrophy', form="FORMATTED", status="REPLACE")
    do it = 1, mean_end-mean_beg+1
       write (funit,'(6(es11.4,1x))') drake_enstrophy(it,1)/DAY, drake_enstrophy(it,2)
    end do
    close (funit)
  end subroutine write_out_1layer

   subroutine write_slice
    ! Writes out results
    integer            :: i, k, v
    integer, parameter :: funit = 400
    character(4)       :: s_time

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 50
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) lon
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 51
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) lat
    close (funit)

    ! coordinates
    write (var_file, '(i2)') 52
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) xcoord_lat
    close (funit)
    
    write (var_file, '(i2)') 53
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) xcoord_lon
    close (funit)

    write (var_file, '(i2)') 54
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) zcoord_lat
    close (funit)
    
    write (var_file, '(i2)') 55
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) zcoord_lon
    close (funit)

    ! data
    write (var_file, '(i2)') 56
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) lat_slice
    close (funit)

    write (var_file, '(i2)') 57
    open (unit=funit, file=trim(run_id)//'.5.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
    write (funit) lon_slice
    close (funit)

    ! Compress files
    write (s_time, '(i4.4)') cp_idx
    command = 'ls -1 '//trim(run_id)//'.5.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.5.'//s_time//'.tgz -T tmp --remove-files &'
    call system (command)
  end subroutine write_slice

  subroutine initialize_stat
    implicit none

    allocate (zonal_av(1:zlevels,Ny(1):Ny(2),nvar_zonal))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*save_levels))
    if (trim(test_case)=="Simple_Physics" .and. climatology) allocate (field2d_simplephys(Nx(1):Nx(2),Ny(1):Ny(2),5*save_levels))
    zonal_av = 0d0
  end subroutine initialize_stat

  subroutine initialize_stat_drake
    implicit none

    if (zlevels == 2) then
       allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),1:nvar_drake))
    else
       allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),1:nvar_1layer))
    end if

    if (zlevels == 2) then
       allocate (drake_ke(1:mean_end-mean_beg+1,1:5))
       allocate (drake_enstrophy(1:mean_end-mean_beg+1,1:6))
    elseif (zlevels == 1) then
       allocate (drake_ke(1:mean_end-mean_beg+1,1:2))
       allocate (drake_enstrophy(1:mean_end-mean_beg+1,1:2))
    end if
  end subroutine initialize_stat_drake

  subroutine initialize_stat_vertical
    implicit none

    allocate (lon_slice(Nx(1):Nx(2),1:zlevels,1:5))
    allocate (lat_slice(Ny(1):Ny(2),1:zlevels,1:5))
    allocate (zcoord_lat(Ny(1):Ny(2),1:zlevels,2), zcoord_lon(Nx(1):Nx(2),1:zlevels,2))
    allocate (eta_lat(Ny(1):Ny(2)), eta_lon(Nx(1):Nx(2)))
  end subroutine initialize_stat_vertical

  subroutine interp_save (dom, i, j, zlev, offs, dims)
    ! Linear interpolation to save levels
    ! Assumes variables have been remapped to original vertical grid
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: p, i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, d, e, k, kk
    real(8) :: dpressure, pressure_lower, pressure_upper, p_s

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    p_s = 0d0
    do k = 1, zlevels
       p_s = p_s + sol(S_MASS,k)%data(d)%elts(id+1)
    end do
    p_s = p_top + grav_accel*p_s

    do kk = 1, save_levels
       ! Find pressure at current levels (not interfaces)
       pressure_lower = p_s
       pressure_upper = 0.5d0 * (a_vert(1)+a_vert(2) + (b_vert(1)+b_vert(2))*p_s)
       k = 1
       do while (pressure_upper > pressure_save(kk))
          k = k+1
          pressure_lower = pressure_upper
          pressure_upper = 0.5d0 * (a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1))*p_s)
       end do
       if (k==1) return ! Skip incorrect values for pentagons
       dpressure =  (pressure_save(kk)-pressure_upper)/(pressure_lower-pressure_upper)

       ! Linear interpolation
       sol_save(S_MASS,kk)%data(d)%elts(id+1) = sol(S_MASS,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_MASS,k)%data(d)%elts(id+1) - sol(S_MASS,k+1)%data(d)%elts(id+1))
       sol_save(S_TEMP,kk)%data(d)%elts(id+1) = sol(S_TEMP,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_TEMP,k)%data(d)%elts(id+1) - sol(S_TEMP,k+1)%data(d)%elts(id+1))
       do e = 1, EDGE
          sol_save(S_VELO,kk)%data(d)%elts(EDGE*id+e) = sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e)  + &
               dpressure * (sol(S_VELO,k)%data(d)%elts(EDGE*id+e) - sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e))
       end do
       if (trim(test_case)=="Simple_Physics" .and. climatology) then
         simple_phys_temp(kk-1)%data(d)%elts(id+1) = simple_phys_temp(k+1)%data(d)%elts(id+1) + &
            dpressure * (simple_phys_temp(k)%data(d)%elts(id+1) - simple_phys_temp(k+1)%data(d)%elts(id+1))
         simple_phys_zonal(kk-1)%data(d)%elts(id+1) = simple_phys_zonal(k+1)%data(d)%elts(id+1) + &
            dpressure * (simple_phys_zonal(k)%data(d)%elts(id+1) - simple_phys_zonal(k+1)%data(d)%elts(id+1))
         simple_phys_merid(kk-1)%data(d)%elts(id+1) = simple_phys_merid(k+1)%data(d)%elts(id+1) + &
            dpressure * (simple_phys_merid(k)%data(d)%elts(id+1) - simple_phys_merid(k+1)%data(d)%elts(id+1))
         do e = 1, EDGE
            simple_phys_vels(kk-1)%data(d)%elts(EDGE*id+e) = simple_phys_vels(k+1)%data(d)%elts(EDGE*id+e) + &
               dpressure * (simple_phys_vels(k)%data(d)%elts(EDGE*id+e) - simple_phys_vels(k+1)%data(d)%elts(EDGE*id+e))
         end do
       end if
    end do
  end subroutine interp_save
end program
