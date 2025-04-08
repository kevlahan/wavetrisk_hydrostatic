program save_vtk_data
  ! Reads data from checkpoints and saves in vtk files
  ! (allows different post-processing data analysis from checkpoint data)
  use main_mod
  use test_case_mod
  use io_mod
  use projection_mod
  use std_atm_profile_mod
  use io_vtk_mod
  implicit none
  
  integer                                :: i, k, l, nt, Ncumul, p, d, nvar_total
  integer, parameter                     :: funit = 400
  integer, parameter                     :: nvar_save = 7, nvar_drake = 12, nvar_1layer = 5
  real(8)                                :: area1, area2
  real(8), dimension(:),     allocatable :: eta_lat, eta_lon
  real(8), dimension(:,:),   allocatable :: drake_ke, drake_enstrophy
  real(8), dimension(:,:,:), allocatable :: field2d_av, field2d_save, field2d_incr, field2d_simplephys
  real(8), dimension(:,:,:), allocatable :: lat_slice, lon_slice, zonal_av, zcoord_lat, zcoord_lon
  
  character(2)                           :: var_file
  character(8)                           :: itype
  character(9999)                        :: bash_cmd, command, topo_filename

  logical, parameter                     :: welford = .true. ! use Welford's one-pass algorithm or naive two-pass algorithm
  
  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters
  
  select case (test_case)
  case ("DCMIP2012c4")
     compressible            = .true.                
     split_mean_perturbation = .false.
  case ("DCMIP2008c5")
     compressible            = .true.                
     split_mean_perturbation = .false.
  case ("climate") ! compile with PHYSICS=true
     compressible            = .true.                    
     split_mean_perturbation = .true.           
     Nsoil                   = 10
     physics_model           = .true.
     physics_type            = "Simple"
  case ("drake")
     compressible            = .false.
     mode_split              = .true.          
     split_mean_perturbation = .true.
  case ("seamount") 
     compressible            = .false.                            
     split_mean_perturbation = .true.
     mode_split              = .true.
     sigma_z                 = .false.
     coords                  = "chebyshev"
     stratification          = "exponential"
  case ("upwelling")
     compressible            = .false.
     split_mean_perturbation = .true.
     sigma_z                 = .true.                       
     coords                  = "croco"                      
     mode_split              = .true.
     vert_diffuse            = .true.
  case ("jet")
     split_mean_perturbation = .true.
     compressible            = .false.                            
     mode_split              = .true.
     penalize                = .false.
     vert_diffuse            = .true.
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

   Nt = 0
  do cp_idx = mean_beg, mean_end
     Nt = Nt + 1
     resume = NONE
     call restart (run_id)

     ! Set means
     do l = level_start, level_end
        do k = 1, zlevels
           call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
        end do
     end do

     ! Save data
     call write_and_export (cp_idx) 
  end do
  
  call finalize
contains
  subroutine cal_zonal_av
    ! Zonal average means and covariances over all checkpoints using stable online algorithm
    ! Uses Welford's stable onlne algorithm
    use domain_mod
    implicit none
    
    integer                                       :: d, ix, j, k
    real(8), dimension (Ny(1):Ny(2))              :: Tprime, Uprime, Vprime, Tprime_new, Uprime_new, Vprime_new, uKEprime, &
         vKEprime, rho_prime
    real(8), dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj, Dproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level_save
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_bdry (exner_fun, NONE, 41)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (exner_fun(k), level_save, 1d0)
       Tproj = field2d

       ! Simple phys Density for KEs
       if (trim (test_case) == "Simple_Physics") then
         call project_field_onto_plane(penal_node(k), level_save, 1.0_8)
         Dproj = field2d
       end if

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

          if (trim (test_case) == "Simple_Physics") then
            uKEprime = (0.5 * Dproj(ix,:)*Uproj(ix,:)**2) - zonal_av(k,:,10)
            vKEprime = (0.5 * Dproj(ix,:)*Vproj(ix,:)**2) - zonal_av(k,:,11)
            rho_prime = Dproj(ix,:) - zonal_av(k,:,12)
            zonal_av(k,:,10) = zonal_av(k,:,10) + uKEprime/Ncumul
            zonal_av(k,:,11) = zonal_av(k,:,11) + vKEprime/Ncumul
            zonal_av(k,:,12) = zonal_av(k,:,12) + rho_prime/Ncumul
          end if

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
    real(8), dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj, Dproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_bdry (exner_fun, NONE, 42)

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

       if (trim (test_case) == "Simple_Physics") then
         call project_field_onto_plane(penal_node(k), level_save, 1.0_8)
         Dproj = field2d
       end if

       ! Update means
       do ix = Nx(1), Nx(2)
          Ncumul = (Nt-1)*(Nx(2)-Nx(1)+1) + ix-Nx(1)+1

          zonal_av(k,:,1) = zonal_av(k,:,1) + Tproj(ix,:)
          zonal_av(k,:,3) = zonal_av(k,:,3) + Uproj(ix,:)
          zonal_av(k,:,4) = zonal_av(k,:,4) + Vproj(ix,:)
          zonal_av(k,:,5) = zonal_av(k,:,5) +  0.5 * (Uproj(ix,:)**2 + Vproj(ix,:)**2)
          if (trim (test_case) == "Simple_Physics") then
            zonal_av(k,:,10) = zonal_av(k,:,10) + 0.5 * Dproj(ix,:) * (Uproj(ix,:)**2)
            zonal_av(k,:,11) = zonal_av(k,:,11) + 0.5 * Dproj(ix,:) * (Vproj(ix,:)**2)
            zonal_av(k,:,12) = zonal_av(k,:,12) + Dproj(ix,:)
          end if
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
    call update_bdry (exner_fun, NONE, 43)

    ! Zonal averages
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (exner_fun(k), level_save, 1.0_8)
       Tproj = field2d

       ! Zonal and meridional velocities
       do d = 1, size(grid)
          velo => sol(S_VELO,k)%data(d)%elts
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

  subroutine latlon (field)
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels)
    use domain_mod
    implicit none
    real(8), dimension (Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*zlevels) :: field
    
    integer :: d, j, k

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    ! Compute vertical velocities at all vertical layers
    call omega_velocity 

    call cal_surf_press (sol(1:N_VARIABLE,1:zmax))

    ! Remap to pressure_save vertical levels for saving data
    call apply_onescale (interp_save, level_save, z_null, -1, 2)

    ! Calculate temperature at all vertical levels (saved in exner_fun) and temperature at interpolated saved vertical levels
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_bdry (exner_fun, NONE, 44)

    ! Latitude-longitude projections
    do k = 1, zlevels
       ! Temperature
       call project_field_onto_plane (trend(1,k), level_save, 0.0_8)
       field(:,:,1+k-1) = field2d

       ! Calculate zonal and meridional velocities and vorticity
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
          velo1 => grid(d)%u_zonal%elts
          velo2 => grid(d)%v_merid%elts
          vort  => grid(d)%vort%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_UVW_latlon, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,          grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (velo, velo1, velo2, vort)
       end do

       ! Zonal velocity
       call project_array_onto_plane ("u_zonal", level_save, 0d0)
       field(:,:,2+k-1) = field2d

       ! Meridional velocity
       call project_array_onto_plane ("v_merid", level_save, 0d0)
       field(:,:,3+k-1) = field2d

       ! Geopotential
       call apply_onescale (cal_geopot, level_save, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call project_array_onto_plane ("geopot", level_save, 1d0)
       field(:,:,4+k-1) = field2d

       ! Vorticity
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          nullify (vort)
       end do
       call project_array_onto_plane ("press_lower", level_save, 1d0)
       field(:,:,5+k-1) = field2d

       ! Surface pressure
       call project_array_onto_plane ("surf_press", level_save, 1d0)
       field(:,:,6+k-1) = field2d

       ! Project vertical velocity
       call project_field_onto_plane (trend(S_TEMP,k), level_save, 0d0)
       field(:,:,7+k-1) = field2d
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
    integer            :: i, info, k, v
    integer, parameter :: funit = 400
    character(9999)    :: s_time

    ! 2d projections
    do v = 1, nvar_save*zlevels
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.4.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_save(:,i,v)
       end do
       close (funit)
    end do


    !Save KE zonal averages for simple Physics (50 and 51)
    if (trim(test_case)=="climate")then
      do v = 10,12
         write (var_file, '(i2)') v+40
         open (unit=funit, file=trim(run_id)//'.4.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
         do k = zlevels,1,-1
            write (funit) zonal_av(k,:,v)
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
    write (funit) (0.5d0*((a_vert(k)+a_vert(k+1))/ref_surf_press + b_vert(k)+b_vert(k+1)), k = zlevels, 1, -1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.4.?? > tmp'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.4.?? > tmp'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
    command = 'gtar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (trim(command))
  end subroutine write_out

  subroutine write_out_av
    ! Writes out time averaged 2d projection results
    integer            :: i, info, k, v
    integer, parameter :: funit = 400

    ! 2d projections
    do v = 1, nvar_save * zlevels
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.6.0'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (funit) field2d_av(:,i,v)
       end do
       close (funit)
    end do

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.6.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.6.'//var_file, access="STREAM", form="UNFORMATTED", status="REPLACE") 
    write (funit) (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.6.?? > tmp'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.6.?? > tmp'
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
    command = 'gtar czf '//trim(run_id)//'.6.tgz -T tmp --remove-files &'
    call system (trim(command))
  end subroutine write_out_av

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
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
    command = 'gtar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (trim(command))

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
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
    command = 'gtar czf '//trim(run_id)//'.4.tgz -T tmp --remove-files &'
    call system (trim(command))

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
    write (bash_cmd,'(a,a,a)') 'bash -c "', trim (command), '"'
    call system (trim(bash_cmd))
    command = 'gtar czf '//trim(run_id)//'.5.'//trim(s_time)//'.tgz -T tmp --remove-files &'
    call system (trim(command))
  end subroutine write_slice

  subroutine initialize_stat
    implicit none

    allocate (zonal_av(1:zlevels,Ny(1):Ny(2),nvar_total))
    
    allocate (field2d_av  (Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*zlevels))
    allocate (field2d_incr(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*zlevels))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*zlevels))
    
    field2d_av = 0d0; field2d_incr = 0d0; field2d_save = 0d0; zonal_av = 0d0
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

    do kk = 1, zlevels
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
       sol(S_MASS,kk)%data(d)%elts(id+1) = sol(S_MASS,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_MASS,k)%data(d)%elts(id+1) - sol(S_MASS,k+1)%data(d)%elts(id+1))
       
       sol(S_TEMP,kk)%data(d)%elts(id+1) = sol(S_TEMP,k+1)%data(d)%elts(id+1) + &
            dpressure * (sol(S_TEMP,k)%data(d)%elts(id+1) - sol(S_TEMP,k+1)%data(d)%elts(id+1))

       do e = 1, EDGE
          sol(S_VELO,kk)%data(d)%elts(EDGE*id+e) = sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e)  + &
               dpressure * (sol(S_VELO,k)%data(d)%elts(EDGE*id+e) - sol(S_VELO,k+1)%data(d)%elts(EDGE*id+e))
       end do

       ! For vertical velocity
       trend(S_TEMP,kk)%data(d)%elts(id+1) = trend(S_TEMP,k+1)%data(d)%elts(id+1) + &
            dpressure * (trend(S_TEMP,k)%data(d)%elts(id+1) - trend(S_TEMP,k+1)%data(d)%elts(id+1))
    end do
  end subroutine interp_save
end program save_vtk_data
