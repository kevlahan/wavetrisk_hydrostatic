program flat_projection_data
  ! Post-processing of checkpoint data to calculate flat projection
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  
  integer                                :: k, l, nt, Ncumul
  integer, parameter                     :: nvar_save = 6, nvar_drake = 12, nvar_1layer = 5
  integer, dimension(2)                  :: Nx, Ny
  
  real(4), dimension(:,:),   allocatable :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export, area1, area2
  real(8), dimension(2)                  :: lon_lat_range
  real(8), dimension(:),     allocatable :: eta_lat, eta_lon, lat, lon
  real(8), dimension(:,:),   allocatable :: drake_ke, drake_enstrophy, xcoord_lat, xcoord_lon
  real(8), dimension(:,:,:), allocatable :: field2d_save, lat_slice, lon_slice, zonal_av, zcoord_lat, zcoord_lon
  
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

     alpha          = 1d-2    ! porosity
     npts_penal     = 2.5

     width          = 80 * KM                           ! width of channel in km
     lat_width      = (width/radius)/DEG
     lat_c          = 45                                ! centre of zonal channel (in degrees)

     coords         = "croco"
  end select

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  resume = mean_beg
  
  ! Initialize variables
  call initialize (run_id)

  ! Initialize statistics
  if (trim (test_case) == "drake") then
     call initialize_stat_drake
  elseif (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling") then
     call initialize_stat_vertical
  else
     call initialize_stat
  end if
  
  ! Calculate zonal average over all check points
  if (rank == 0) write (6,'(A)') "Calculating zonal averages over all checkpoints"

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
     elseif  (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling") then
        call vertical_slice
        if (rank == 0) call write_slice
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
  elseif (.not. (trim (test_case) == "seamount" .or. trim (test_case) == "upwelling")) then
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
    call cal_surf_press (sol)
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 41)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       call project_uzonal_onto_plane (level_save, 0.0_8)
       Uproj = field2d
       call project_vmerid_onto_plane (level_save, 0.0_8)
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
    call cal_surf_press (sol)
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 42)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo, velo1, velo2)
       end do
       call project_uzonal_onto_plane (level_save, 0.0_8)
       Uproj = field2d
       call project_vmerid_onto_plane (level_save, 0.0_8)
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
    call cal_surf_press (sol)
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE, 43)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo => sol(S_VELO,k)%data(d)%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), k, 0, 1)
          end do
          nullify (velo)
       end do
       call project_uzonal_onto_plane (level_save, 0.0_8)
       Uproj = field2d
       call project_vmerid_onto_plane (level_save, 0.0_8)
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

    call cal_surf_press (sol)

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
       call project_onto_plane (trend(1,k), level_save, 0.0_8)
       field2d_save(:,:,1+k-1) = field2d

       ! Calculate zonal and meridional velocities and vorticity
       do d = 1, size(grid)
          velo  => sol_save(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (velo, velo1, velo2, vort)
       end do

       ! Zonal velocity
       call project_uzonal_onto_plane (level_save, 0.0_8)
       field2d_save(:,:,2+k-1) = field2d
       
       ! Meridional velocity
       call project_vmerid_onto_plane (level_save, 0.0_8)
       field2d_save(:,:,3+k-1) = field2d

       ! Geopotential
       call apply_onescale (cal_geopot, level_save, z_null, 0, 1)
       call project_geopot_onto_plane (level_save, 1.0_8)
       field2d_save(:,:,4+k-1) = field2d

       ! Vorticity
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          nullify (vort)
       end do
       call project_vorticity_onto_plane (level_save, 1.0_8)
       field2d_save(:,:,5+k-1) = field2d

       ! Surface pressure
       call project_surf_press_onto_plane (level_save, 1.0_8)
       field2d_save(:,:,6+k-1) = field2d
    end do
  end subroutine latlon

  subroutine vertical_slice
    ! Save vertical slicse along given latitude and longitude for incompressible test cases
    use barotropic_2d_mod
    implicit none
    integer                         :: d, i, j, k, l, idx_lat, idx_lon
    real(8)                         :: dz, eta, z_s
    real(8), dimension(0:zlevels)   :: z
    real(8), dimension(Ny(1):Ny(2)) :: bathy_lat
    real(8), dimension(Nx(1):Nx(2)) :: bathy_lon
    logical, parameter              :: zonal = .true. ! compute zonal averages

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
             call apply_onescale_to_patch (cal_density,    grid(d), grid(d)%lev(l)%elts(j), k,      -2, 3)
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (mass, mean_m, mean_t, scalar, temp, velo, divu, velo1, velo2, vort)
       end do
       call project_uzonal_onto_plane (l, 0.0_8)
        if (zonal) then
          lat_slice(:,k,1) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,1) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,1) = field2d(:,idx_lat)
       
       call project_vmerid_onto_plane (l, 0.0_8)
        if (zonal) then
          lat_slice(:,k,2) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,2) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,2) = field2d(:,idx_lat)
       
       call project_onto_plane (sol(S_TEMP,zlevels+1), l, 0.0_8)
        if (zonal) then
          lat_slice(:,k,3) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,3) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,3) = field2d(:,idx_lat)

       call project_vorticity_onto_plane (l, 1.0_8)
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
       call project_w_onto_plane (k, l, 0.0_8)
       if (zonal) then
          lat_slice(:,k,5) = sum(field2d, 1) / size(field2d,1) 
       else
          lat_slice(:,k,5) = field2d(idx_lon,:)
       end if
       lon_slice(:,k,5) = field2d(:,idx_lat)
    end do
    
    ! Set free surface
    call project_onto_plane (sol(S_MASS,zlevels+1), l, 0.0_8)
    eta_lat = field2d(idx_lon,:)
    eta_lon = field2d(:,idx_lat)

    ! Set bathymetry
    call project_topo_onto_plane (l, 0.0_8)
    bathy_lat = field2d(idx_lon,:)
    bathy_lon = field2d(:,idx_lat)
    
    do i = Ny(1), Ny(2)
       xcoord_lat(i,1) = lat(i) - dy_export/2 / DEG
       xcoord_lat(i,2) = lat(i) + dy_export/2 / DEG
       
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
       xcoord_lon(i,1) = lon(i) - dx_export/2 / DEG
       xcoord_lon(i,2) = lon(i) + dx_export/2 / DEG
       
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
    
    if (trim(itype) == 'adaptive') then
       energy_drake(1) = integrate_adaptive (layer1_ke, z_null)
       energy_drake(2) = integrate_adaptive (layer2_ke, z_null)
       energy_drake(3) = integrate_adaptive (barotropic_ke, z_null)
    elseif (trim(itype) == 'coarse') then
       energy_drake(1) = integrate_hex (layer1_ke,     level_start, z_null)
       energy_drake(2) = integrate_hex (layer2_ke,     level_start, z_null)
       energy_drake(3) = integrate_hex (barotropic_ke, level_start, z_null) 
    end if
    energy_drake(4) = energy_drake(1) + energy_drake(2) - energy_drake(3) ! baroclinic energy
    
    energy_drake = energy_drake / (4*MATH_PI*radius**2*ref_density)
  end function energy_drake
  
  real(8) function energy_1layer (itype)
    ! Calculates baroclinic and bartoropic energies for two layer case
    use io_mod
    implicit none
    character(*) :: itype
    
    if (trim(itype) == 'adaptive') then
       energy_1layer = integrate_adaptive (one_layer_ke, z_null)
    elseif (trim(itype) == 'coarse') then
       energy_1layer = integrate_hex (one_layer_ke, level_start, z_null) 
    end if
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
       if (trim(itype) == 'adaptive') then
          pot_enstrophy_drake(k) = integrate_adaptive (pot_enstrophy, k)
       elseif (trim(itype) == 'coarse') then
          pot_enstrophy_drake(k) = integrate_hex (pot_enstrophy, level_start, k)
       end if
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
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
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
    if (trim(itype) == 'adaptive') then
       pot_enstrophy_drake(3) = integrate_adaptive (pot_enstrophy, 3)
    elseif (trim(itype) == 'coarse') then
       pot_enstrophy_drake(3) = integrate_hex (pot_enstrophy, level_start, 3)
    end if
    
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
                call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
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
       if (trim(itype) == 'adaptive') then
          pot_enstrophy_drake(k+3) = integrate_adaptive (pot_enstrophy, k)
       elseif (trim(itype) == 'coarse') then
          pot_enstrophy_drake(k+3) = integrate_hex (pot_enstrophy, level_start, k)
       end if
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
    
    if (trim(itype) == 'adaptive') then
       pot_enstrophy_1layer = integrate_adaptive (pot_enstrophy, 1)
    elseif (trim(itype) == 'coarse') then
       pot_enstrophy_1layer = integrate_hex (pot_enstrophy, level_start, 1)
    end if
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
          call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
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
    call project_uzonal_onto_plane (l, 0.0_8)
    field2d_save(:,:,1) = field2d
    call project_vmerid_onto_plane (l, 0.0_8)
    field2d_save(:,:,2) = field2d
    ! Project barotropic vorticity onto plane
    call project_vorticity_onto_plane (l, 1.0_8)
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
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
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
       call project_uzonal_onto_plane (l, 0.0_8)
       field2d_save(:,:,3*k+1) = field2d
       call project_vmerid_onto_plane (l, 0.0_8)
       field2d_save(:,:,3*k+2) = field2d
       ! Project layer k barotropic vorticity onto plane
       call project_vorticity_onto_plane (l, 1.0_8)
       field2d_save(:,:,3*k+3) = field2d
    end do

    ! Free surface
    call project_freesurface_onto_plane (l, 1.0_8)
    field2d_save(:,:,10) = field2d

    ! Internal free surface
    call project_internal_freesurface_onto_plane (l, 1.0_8)
    field2d_save(:,:,11) = field2d

    ! Penalization
    do k = 1, 2
       do d = 1, size(grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (set_penal, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
       end do
    end do
    call project_penal_onto_plane (l, 1.0_8)
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
          call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
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
    call project_uzonal_onto_plane (l, 0.0_8)
    field2d_save(:,:,1) = field2d
    call project_vmerid_onto_plane (l, 0.0_8)
    field2d_save(:,:,2) = field2d
    ! Project vorticity onto plane
    call project_vorticity_onto_plane (l, 1.0_8)
    field2d_save(:,:,3) = field2d

    ! Free surface
    call project_freesurface_onto_plane (l, 1.0_8)
    field2d_save(:,:,4) = field2d

    ! Penalization
    do k = 1, 2
       do d = 1, size(grid)
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (set_penal, grid(d), grid(d)%lev(l)%elts(j), k, 0, 1)
          end do
       end do
    end do
    call project_penal_onto_plane (l, 1.0_8)
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
    integer :: i

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (zonal_av(1:zlevels,Ny(1):Ny(2),nvar_zonal))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*save_levels))
    allocate (lat(Ny(1):Ny(2)), lon(Nx(1):Nx(2)))
    zonal_av = 0.0_8

    do i = Nx(1), Nx(2)
       lon(i) = -180+dx_export*(i-Nx(1))/MATH_PI*180
    end do

    do i = Ny(1), Ny(2)
       lat(i) = -90+dy_export*(i-Ny(1))/MATH_PI*180
    end do
  end subroutine initialize_stat

  subroutine initialize_stat_drake
    implicit none
    integer :: i

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))

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
    
    allocate (lat(Ny(1):Ny(2)), lon(Nx(1):Nx(2)))

    do i = Nx(1), Nx(2)
       lon(i) = -180+dx_export*(i-Nx(1))/MATH_PI*180
    end do

    do i = Ny(1), Ny(2)
       lat(i) = -90+dy_export*(i-Ny(1))/MATH_PI*180
    end do
  end subroutine initialize_stat_drake

  subroutine initialize_stat_vertical
    implicit none
    integer :: i

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (lon_slice(Nx(1):Nx(2),1:zlevels,1:5))
    allocate (lat_slice(Ny(1):Ny(2),1:zlevels,1:5))
    allocate (lat(Ny(1):Ny(2)), lon(Nx(1):Nx(2)))
    allocate (xcoord_lat(Ny(1):Ny(2),1:2), xcoord_lon(Nx(1):Nx(2),1:2))
    allocate (zcoord_lat(Ny(1):Ny(2),1:zlevels,2), zcoord_lon(Nx(1):Nx(2),1:zlevels,2))
    allocate (eta_lat(Ny(1):Ny(2)), eta_lon(Nx(1):Nx(2)))

    do i = Nx(1), Nx(2)
       lon(i) = -180+dx_export*(i-Nx(1))/MATH_PI*180
    end do

    do i = Ny(1), Ny(2)
       lat(i) = -90+dy_export*(i-Ny(1))/MATH_PI*180
    end do
  end subroutine initialize_stat_vertical

  subroutine project_onto_plane (field, l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val
    Type(Float_field)     :: field

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = field%data(d)%elts(id+1)
                valN  = field%data(d)%elts(idN+1)
                valE  = field%data(d)%elts(idE+1)
                valNE = field%data(d)%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_onto_plane

  subroutine project_geopot_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%geopot_lower%elts(id+1)
                valN  = grid(d)%geopot_lower%elts(idN+1)
                valE  = grid(d)%geopot_lower%elts(idE+1)
                valNE = grid(d)%geopot_lower%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_geopot_onto_plane

  subroutine project_surf_press_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%surf_press%elts(id+1)/100
                valN  = grid(d)%surf_press%elts(idN+1)/100
                valE  = grid(d)%surf_press%elts(idE+1)/100
                valNE = grid(d)%surf_press%elts(idNE+1)/100

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_surf_press_onto_plane

  subroutine project_vorticity_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%press_lower%elts(id+1)
                valN  = grid(d)%press_lower%elts(idN+1)
                valE  = grid(d)%press_lower%elts(idE+1)
                valNE = grid(d)%press_lower%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_vorticity_onto_plane

  subroutine project_baroclinic_vorticity_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%geopot_lower%elts(id+1)
                valN  = grid(d)%geopot_lower%elts(idN+1)
                valE  = grid(d)%geopot_lower%elts(idE+1)
                valNE = grid(d)%geopot_lower%elts(idNE+1)

                if (abs(cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs(cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_baroclinic_vorticity_onto_plane

  subroutine project_w_onto_plane (k, l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: k, l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = trend(S_TEMP,k)%data(d)%elts(id+1)
                valN  = trend(S_TEMP,k)%data(d)%elts(idN+1)
                valE  = trend(S_TEMP,k)%data(d)%elts(idE+1)
                valNE = trend(S_TEMP,k)%data(d)%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_w_onto_plane

  subroutine project_uzonal_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%u_zonal%elts(id+1)
                valN  = grid(d)%u_zonal%elts(idN+1)
                valE  = grid(d)%u_zonal%elts(idE+1)
                valNE = grid(d)%u_zonal%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_uzonal_onto_plane

  subroutine project_baroclinic_uzonal_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = trend(S_VELO,1)%data(d)%elts(id+1)
                valN  = trend(S_VELO,1)%data(d)%elts(idN+1)
                valE  = trend(S_VELO,1)%data(d)%elts(idE+1)
                valNE = trend(S_VELO,1)%data(d)%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_baroclinic_uzonal_onto_plane

   subroutine project_baroclinic_vmerid_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = trend(S_VELO,2)%data(d)%elts(id+1)
                valN  = trend(S_VELO,2)%data(d)%elts(idN+1)
                valE  = trend(S_VELO,2)%data(d)%elts(idE+1)
                valNE = trend(S_VELO,2)%data(d)%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt (1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_baroclinic_vmerid_onto_plane

  subroutine project_vmerid_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%v_merid%elts(id+1)
                valN  = grid(d)%v_merid%elts(idN+1)
                valE  = grid(d)%v_merid%elts(idE+1)
                valNE = grid(d)%v_merid%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_vmerid_onto_plane

  subroutine project_topo_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = grid(d)%topo%elts(id+1)
                valN  = grid(d)%topo%elts(idN+1)
                valE  = grid(d)%topo%elts(idE+1)
                valNE = grid(d)%topo%elts(idNE+1)

                if (abs (cN(2) - MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_topo_onto_plane

   subroutine project_freesurface_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use ops_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = sol(S_MASS,zlevels+1)%data(d)%elts(id+1)   / phi_node (d, id+1,   zlevels)
                valN  = sol(S_MASS,zlevels+1)%data(d)%elts(idN+1)  / phi_node (d, idN+1,  zlevels)
                valE  = sol(S_MASS,zlevels+1)%data(d)%elts(idE+1)  / phi_node (d, idE+1,  zlevels)
                valNE = sol(S_MASS,zlevels+1)%data(d)%elts(idNE+1) / phi_node (d, idNE+1, zlevels)

                if (abs (cN(2) - MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_freesurface_onto_plane

  subroutine project_internal_freesurface_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use ops_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = sol(S_MASS,1)%data(d)%elts(id+1)   / (ref_density * phi_node (d, id+1,   1))
                valN  = sol(S_MASS,1)%data(d)%elts(idN+1)  / (ref_density * phi_node (d, idN+1,  1))
                valE  = sol(S_MASS,1)%data(d)%elts(idE+1)  / (ref_density * phi_node (d, idE+1,  1))
                valNE = sol(S_MASS,1)%data(d)%elts(idNE+1) / (ref_density * phi_node (d, idNE+1, 1))

                if (abs (cN(2) - MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_internal_freesurface_onto_plane

  subroutine project_penal_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use ops_mod
    use comm_mpi_mod
    integer               :: l, itype
    real(8)               :: default_val

    integer                        :: d, i, j, jj, p, c, p_par, l_cur
    integer                        :: id, idN, idE, idNE
    real(8)                        :: val, valN, valE, valNE
    real(8), dimension(2)          :: cC, cN, cE, cNE
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    field2d = default_val
    do d = 1, size(grid)
       do jj = 1, grid(d)%lev(l)%length
          call get_offs_Domain (grid(d), grid(d)%lev(l)%elts(jj), offs, dims)
          do j = 0, PATCH_SIZE-1
             do i = 0, PATCH_SIZE-1
                id   = idx(i,   j,   offs, dims)
                idN  = idx(i,   j+1, offs, dims)
                idE  = idx(i+1, j,   offs, dims)
                idNE = idx(i+1, j+1, offs, dims)

                call cart2sph2 (grid(d)%node%elts(id+1),   cC)
                call cart2sph2 (grid(d)%node%elts(idN+1),  cN)
                call cart2sph2 (grid(d)%node%elts(idE+1),  cE)
                call cart2sph2 (grid(d)%node%elts(idNE+1), cNE)

                val   = penal_node(1)%data(d)%elts(id+1)   
                valN  = penal_node(1)%data(d)%elts(idN+1)  
                valE  = penal_node(1)%data(d)%elts(idE+1)  
                valNE = penal_node(1)%data(d)%elts(idNE+1) 

                if (abs (cN(2) - MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cNE, (/cNE(1), cN(2)/), cC, (/valNE, valN, val/))
                   call interp_tri_to_2d_and_fix_bdry ((/cNE(1), cN(2)/), (/cC(1), cN(2)/), cC, (/valN, valN, val/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cNE, cN, cC, (/valNE, valN, val/))
                end if
                if (abs (cE(2) + MATH_PI/2) < sqrt(1d-15)) then
                   call interp_tri_to_2d_and_fix_bdry (cC, (/cC(1), cE(2)/), cNE, (/val, valE, valNE/))
                   call interp_tri_to_2d_and_fix_bdry ((/cC(1), cE(2)/), (/cNE(1), cE(2)/), cNE, (/valE, valE, valNE/))
                else
                   call interp_tri_to_2d_and_fix_bdry (cC, cE, cNE, (/val, valE, valNE/))
                end if
             end do
          end do
       end do
    end do
    ! Synchronize array over all processors
    sync_val = default_val
    call sync_array (field2d(Nx(1),Ny(1)), size(field2d))
  end subroutine project_penal_onto_plane

  subroutine interp_tri_to_2d (a, b, c, val)
    real(8), dimension(2) :: a, b, c
    real(8), dimension(3) :: val

    integer               :: id_x, id_y
    real(8)               :: ival, minx, maxx, miny, maxy
    real(8), dimension(2) :: ll
    real(8), dimension(3) :: bac
    logical               :: inside

    minx = min (min (a(1), b(1)), c(1))
    maxx = max (max (a(1), b(1)), c(1))
    miny = min (min (a(2), b(2)), c(2))
    maxy = max (max (a(2), b(2)), c(2))
    if (maxx-minx > MATH_PI/2) then
       write (0,'(A,i4,A)') 'ERROR (rank = ', rank, '): io-333 "export"'
       return
    end if

    do id_x = floor (kx_export*minx), ceiling (kx_export*maxx)
       if (id_x < lbound (field2d,1) .or. id_x > ubound (field2d,1)) cycle
       do id_y = floor (ky_export*miny), ceiling (ky_export*maxy)
          if (id_y < lbound (field2d,2) .or. id_y > ubound (field2d,2)) cycle
          ll = (/dx_export*id_x, dy_export*id_y/)
          call interp_tria (ll, a, b, c, val, ival, inside)
          if (inside) field2d(id_x,id_y) = ival
       end do
    end do
  end subroutine interp_tri_to_2d

  subroutine interp_tri_to_2d_and_fix_bdry (a0, b0, c0, val)
    implicit none
    real(8), dimension(2) :: a0, b0, c0
    real(8), dimension(3) :: val

    integer               :: i
    integer, dimension(3) :: fixed
    real(8), dimension(2) :: a, b, c

    a = a0
    b = b0
    c = c0
    call fix_boundary (a(1), b(1), c(1), fixed(1))
    call fix_boundary (b(1), c(1), a(1), fixed(2))
    call fix_boundary (c(1), a(1), b(1), fixed(3))
    call interp_tri_to_2d (a, b, c, val)

    if (sum(abs(fixed)) > 1) write (0,'(A)') 'ALARM'

    if (sum(fixed) /= 0) then
       a(1) = a(1) - sum(fixed) * 2*MATH_PI
       b(1) = b(1) - sum(fixed) * 2*MATH_PI
       c(1) = c(1) - sum(fixed) * 2*MATH_PI
       call interp_tri_to_2d (a, b, c, val)
    end if
  end subroutine interp_tri_to_2d_and_fix_bdry

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

    p_s = 0.0_8
    do k = 1, zlevels
       p_s = p_s + sol(S_MASS,k)%data(d)%elts(id+1)
    end do
    p_s = p_top + grav_accel*p_s

    do kk = 1, save_levels
       ! Find pressure at current levels (not interfaces)
       pressure_lower = p_s
       pressure_upper = 0.5*(a_vert(1)+a_vert(2) + (b_vert(1)+b_vert(2))*p_s)
       k = 1
       do while (pressure_upper > pressure_save(kk))
          k = k+1
          pressure_lower = pressure_upper
          pressure_upper = 0.5*(a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1))*p_s)
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
    end do
  end subroutine interp_save
end program

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Physics routines for this test case (including diffusion)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function physics_scalar_flux (dom, id, idE, idNE, idN, type)
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP,1:EDGE) :: physics_scalar_flux
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical, optional                        :: type

  physics_scalar_flux(S_MASS,:) = 0.0_8
end function physics_scalar_flux

function physics_scalar_source (dom, i, j, zlev, offs, dims)
  ! Additional physics for the source term of the scalar trend
  ! Newton cooling to equilibrium potential temperature theta_equil
  use domain_mod
  implicit none

  real(8), dimension(S_MASS:S_TEMP) :: physics_scalar_source
  type(Domain)                      :: dom
  integer                           :: i, j, zlev
  integer, dimension(N_BDRY+1)      :: offs
  integer, dimension(2,N_BDRY+1)    :: dims

  physics_scalar_source(S_MASS) = 0.0_8
  physics_scalar_source(S_TEMP) = 0.0_8
end function physics_scalar_source

function physics_velo_source (dom, i, j, zlev, offs, dims)
  use domain_mod
  implicit none

  real(8), dimension(1:EDGE)     :: physics_velo_source
  type(Domain)                   :: dom
  integer                        :: i, j, zlev
  integer, dimension(N_BDRY+1)   :: offs
  integer, dimension(2,N_BDRY+1) :: dims

  physics_velo_source = 0.0_8
end function physics_velo_source




