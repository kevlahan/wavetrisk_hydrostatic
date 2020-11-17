program spherical_harmonics
  ! Use SHTOOLS to compute global and local spherical harmonics power spectra from wavetrisk data interpolated to a uniform grid
  ! and projected to a uniform latitude longitude grid
  ! (see https://shtools.github.io/SHTOOLS/index.html for information about SHTOOLS)
  ! Wieczorek, M. A. and F. J. Simons, Minimum-variance multitaper spectral estimation on the sphere, J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692, 2007.
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  
  integer                                :: idata_loc, k, l, nmax
  integer, dimension(2)                  :: Nx, Ny
  integer, parameter                     :: nvar_save = 6, nvar_drake = 12, nvar_1layer = 5
  
  real(8), dimension(:),   allocatable   :: data, data_loc, lat, lat_loc, lon, lon_loc
  real(8), dimension(2)                  :: lon_lat_range
  real(4), dimension(:,:), allocatable   :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export
  
  character(2)                           :: var_file
  character(8)                           :: itype
  character(130)                         :: command

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read test case parameters
  call read_test_case_parameters

  if (trim (test_case) == 'DCMIP2012c4') then
     compressible   = .true.                      ! Compressible equations

     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
     omega          = 7.29212d-5                  ! Earth’s angular velocity in radians per second
     p_0            = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = p_0                         ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

     u_0            = 35.0_8                      ! maximum velocity of zonal wind
     eta_0          = 0.252_8                     ! value of eta at reference level (level of the jet)
  elseif (trim (test_case) == "DCMIP2008c5") then
     compressible   = .true.                      ! Compressible equations

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
  elseif (trim (test_case) == "Held_Suarez") then
     compressible   = .true.                      ! Compressible equations

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
  elseif (trim (test_case) == "drake") then
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
     density     = (/ ref_density, ref_density + drho/2 /)    ! densities in each layer
     height      = (/ abs(max_depth - halocline), abs(halocline) /) ! depths of each layer
     npts_penal  = 4

     mode_split     = .true.                             ! split barotropic mode if true
     compressible   = .false.                            ! always run with incompressible equations
     penalize       = .true.                             ! penalize land regions
  else
     write (6,'(A)') "Test case not supported"
     stop
  end if

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  resume = cp_idx

  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, run_id)

  resume = NONE
  call restart (set_thresholds, load, run_id)

  if (trim (spec_type) == 'sphere') then
     call spec_sphere
  elseif (trim (spec_type) == 'latlon') then
     if (zlevels == 2 .and. trim (test_case) == "drake") then
        call spec_latlon_2layer
     else
        call spec_latlon_1layer
     end if
  end if
  call finalize
contains
  subroutine spec_latlon_1layer
    ! Compute energy spectrum from 2d latitude-longitude projection
    use domain_mod
    implicit none
    integer :: d, i, j, l

    if (rank == 0) write (6,'(A,i6)') "Energy spectrum of checkpoint file = ", cp_idx

    call initialize_lonlat

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_fill
    call fill_up_grid_and_IWT (l)
    
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
    call project_vorticity_onto_plane (l, 1.0_8)
    if (rank == 0) call spectrum_lon_lat ('vort')
    deallocate (field2d)
  end subroutine spec_latlon_1layer

  subroutine spec_latlon_2layer
    ! Compute energy spectrum from 2d latitude-longitude projection for 2 layer case
    use domain_mod
    implicit none
    integer                              :: d, i, ibeg, ibeg_m, iend, iend_m, j, k, l
    real(8), dimension(:,:), allocatable :: dz
    character(4)                         :: var_file
    character(255)                       :: data_type

    if (rank == 0) write (6,'(A,i6)') "Energy spectrum of checkpoint file = ", cp_idx

    call initialize_lonlat_2layer

     ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_fill
    call fill_up_grid_and_IWT (l)
    trend = sol

    ! Set mean on filled grid
    do k = 1, zmax
       call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do

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
    field2d = 0.0_8
    call project_vorticity_onto_plane (l, 1.0_8)
    if (rank == 0) call spectrum_lon_lat ("barotropic")
    
    ! Baroclinic velocity and vorticity in each layer
    do k = 1, 2
      ! Baroclinic velocity
       do d = 1, size(grid)
          ibeg = (1+2*(POSIT(S_VELO)-1))*grid(d)%patch%elts(2+1)%elts_start + 1
          iend = sol(S_VELO,1)%data(d)%length
          trend(S_VELO,zlevels+1)%data(d)%elts(ibeg:iend) = &
               sol(S_VELO,k)%data(d)%elts(ibeg:iend) - sol(S_VELO,3)%data(d)%elts(ibeg:iend)
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
       
       field2d = 0.0_8
       call project_vorticity_onto_plane (l, 1.0_8)
       if (rank == 0) then
          write (data_type, '(a11,i1)') "baroclinic_", k
          call spectrum_lon_lat (data_type)
       end if
    end do
    deallocate (field2d)
  end subroutine spec_latlon_2layer
  
  subroutine spectrum_lon_lat (data_type)
    use SHTOOLS
    implicit none
!!$  input 
!!$  A 2D equally sampled (n,n) grid  (default) or equally spaced grid (n,2*n) that conforms to the sampling theorem of Driscoll and Healy (1994).
!!$  The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/n degrees.
!!$  The first longitudinal band is 0 E, the longitude band for 360 E is not included, and the longitudinal sampling interval is 360/n
!!$  for an equally sampled grid and 180/n for an equally spaced grid, respectively.
!!$  Input
!!$  n (integer): size of data griddh (n,2*n) or (n,n)
!!$  output 
!!$  cilm (real(8), dimension (2,n/2,n/2) or (2, lmax_calc+1, lmax_calc+1)):
!!$  The real spherical harmonic coefficients of the function. These will be exact if the function is bandlimited to degree lmax=n/2-1.
!!$  The coefficients c1lm and c2lm refer to the cosine (clm) and sine (slm) coefficients, respectively, with clm=cilm(1,l+1,m+1) and slm=cilm(2,l+1,m+1).
!!$
!!$  lmax (integer):     
!!$  The maximum spherical harmonic bandwidth of the input grid, which is n/2-1.
!!$  If the optional parameter lmax_calc is not specified, this corresponds to the maximum spherical harmonic degree of the output coefficients cilm.
    character(*)                            :: data_type
    
    integer                                 :: ierr, i, j, lmax, lwin
    real(8)                                 :: area
    real(8), dimension (:,:,:), allocatable :: cilm      ! spherical harmonic coefficients
    real(8), dimension (:),     allocatable :: pspectrum ! global power spectrum of the function
    character(4)                            :: var_file

    ! Save latlon data
    if (rank == 0) then
       write (var_file, '(i4.4)') cp_idx
       open (unit=10, file=trim(run_id)//'_'//var_file//'_'//trim(data_type), access="STREAM", &
            form="UNFORMATTED", status="REPLACE")
       do i = Ny(1), Ny(2)
          write (10) field2d(:,i)
       end do
       close (10)
    end if
    
    lmax = N/4 - 1
    allocate (cilm(2, lmax+1, lmax+1))

    ! Expand latitude-longitude data in spherical harmonics (reorder elements in latitude from N to S)
    call SHExpandDH (transpose(dble(field2d(Nx(1):Nx(2)-1,Ny(2):Ny(1)+1:-1))), N/2, cilm, lmax, norm=1, sampling=2, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHExpandDH"
       
    ! Calculate power spectrum
    allocate (pspectrum(lmax+1))
    call SHPowerSpectrum (cilm, lmax, pspectrum, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHPowerSpectrum"
    
    write (var_file, '(i4.4)') cp_idx
    open (unit=10, file=trim(run_id)//'_'//var_file//'_'//trim(data_type)//'_spec', form="FORMATTED", status="REPLACE")
    area = 4*MATH_PI*radius**2
    do j = 1, lmax + 1
       write (6, '(i4,1x,es10.4)') j, pspectrum(j) * area
       write (10,'(i4,1x,es10.4)') j, pspectrum(j) * area
    end do
    close (10)

    if (local_spec) then
       open (unit=10, file=trim(run_id)//'_'//var_file//'_'//trim(data_type)//'_local_spec', form="FORMATTED", status="REPLACE")
       
       ! Determine the spherical-harmonic bandwidth needed to achieve specified concentration factor
       lwin = SHFindLWin (theta0*DEG, m, concentration)

       ! Compute local spectrum and associated error
       call local_spectrum (cilm, lmax, lwin)
    end if
    deallocate (cilm, pspectrum)
  end subroutine spectrum_lon_lat
  
  subroutine local_spectrum (cilm, lmax, lwin)
    ! Computes a localized multitaper spectral analysis of an input function expressed in spherical harmonics. The maximum degree of the
    ! localized multitaper cross-power spectrum estimate is lmax-lmaxt. lat0 (deg), lon0 (deg), theta0 (radians) define local region
    use SHTOOLS
    implicit none
    integer                                :: lmax, lwin
    real(8), dimension (2, lmax+1, lmax+1) :: cilm
   
    integer                               :: ierr, j
    character(4)                          :: var_file

    real(8)                               :: area
    real(8), dimension (:), allocatable   :: eigenvalues ! concentration factors of the concentration windows
    real(8), dimension (:), allocatable   :: mtse        ! local multitaper power spectrum estimate
    real(8), dimension (:), allocatable   :: sd          ! errors of the localized multitaper power spectral estimates
    integer, dimension (:), allocatable   :: taper_order ! array containing the angular orders of the spherical harmonic coefficients in each column of the array tapers
                                                         ! (each window has non-zero coefficients for a single angular order that is specified in the array taper_order)
    real(8), dimension (:,:), allocatable :: tapers      ! array of the windowing functions, arranged in columns (obtained from SHReturnTapers)
                                                         ! (each window has non-zero coefficients for a single angular order that is specified in the array taper_order)

    allocate (eigenvalues((lwin+1)**2), mtse(lmax-lwin+1), sd(lmax-lwin+1), taper_order((lwin+1)**2), tapers(lwin+1,(lwin+1)**2))

    ! Compute the eigenfunctions of the spherical-cap concentration problem
    call SHReturnTapers (theta0*DEG, lwin, tapers, eigenvalues, taper_order, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHReturnTapers"

    ! Compute localized multitaper power spectrum
    call SHMultiTaperSE (mtse, sd, cilm, lmax, tapers, taper_order, lwin, ntaper, lat=lat0, lon=lon0+180, norm=1, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHMultiTaperSE"

    ! Save data
!!$    area = 2*MATH_PI*radius**2 * (1.0_8 - cos(theta0*DEG))
    area =  4*MATH_PI*radius**2
    do j = 1, lmax - lwin + 1
       write (6, '(i4,1x,2(es10.4,1x))') j, mtse(j) * area, sd(j) * area
       write (10,'(i4,1x,2(es10.4,1x))') j, mtse(j) * area, sd(j) * area
    end do
    close (10)
   
    ! Write out characteristics of local analysis
    write (6,'(a,i4)') "Spherical harmonic bandwidth of window = ", lwin
    do j = 1, ntaper
       write (6,'(a,i3,a,f6.4)') "Concentration factor of taper ", j, " = ", eigenvalues(j)
    end do
    do j = 1, ntaper
       write (6,'(2(a,i3))') "Angular order of taper ", j," = ", taper_order(j)
    end do

    deallocate (eigenvalues, mtse, sd, taper_order, tapers)
  end subroutine local_spectrum

  subroutine spec_sphere
    ! Combine vorticity and coordinate data on triangles into single vectors on rank 0 and compute spectrum
    use domain_mod
    implicit none
    integer                       :: d, j, l, n_loc
    integer, dimension(n_process) :: displs, n_glo
    logical, parameter            :: hex = .false.

    if (rank == 0) write (6,'(A,i6)') "Computing spherical harmonics spectrum of file = ", cp_idx

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_fill
    call fill_up_grid_and_IWT (l)

    ! Vorticity on triangles
    do d = 1, size(grid)
       velo  => sol(S_VELO,1)%data(d)%elts
       vort  => grid(d)%vort%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_vort, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
       end do
       call apply_to_penta_d (post_vort, grid(d), l, z_null)
       nullify (velo, vort)
    end do

    if (.not. hex) then
       nmax  = 2 * (4**l*10 + 1)                       ! total number of triangle points          
       n_loc = 2 * 4**(l - DOMAIN_LEVEL) * size(grid)  ! number of grid points on each rank (not including poles)
    end if

    ! Vorticity on hexagons
    if (hex) then
       do d = 1, size(grid)
          vort => grid(d)%press_lower%elts
          do j = 1, grid(d)%lev(l)%length
             call apply_onescale_to_patch (vort_triag_to_hex, grid(d), grid(d)%lev(l)%elts(j), z_null, -1, 1)
          end do
          nullify (vort)
       end do
       nmax  = 4**l*10 + 1                         ! total number of triangle points          
       n_loc = 4**(l - DOMAIN_LEVEL) * size(grid)  ! number of grid points on each rank (not including poles)
    end if
    if (rank == 0) n_loc = n_loc + 2                         ! store pole data on rank 0
    allocate (data_loc(n_loc), lat_loc(n_loc), lon_loc(n_loc)); data_loc = 0.0_8; lat_loc = 0.0_8; lon_loc = 0.0_8
    allocate (data(nmax), lat(nmax), lon(nmax));  data = 0.0_8; lat = 0.0_8; lon = 0.0_8

    ! Store data for spherical harmonics transform
    idata_loc = 0
    
    if (hex) then
       if (rank==0) call apply_to_pole (define_data_hex, min_level-1, z_null, 1, .false.)
       call apply_onescale (define_data_hex, l, z_null, 0, 0)
    else
       if (rank==0) call apply_to_pole (define_data_triag, min_level-1, z_null, 1, .false.)
       call apply_onescale (define_data_triag, l, z_null, 0, 0)
    end if

    ! Collect data lengths and compute displacements on rank 0
    call MPI_Gather (n_loc, 1, MPI_INTEGER, n_glo, 1, MPI_INTEGER, 0, MPI_COMM_WORLD, ierror)
    if (rank==0) then
       displs(1) = 0
       do j = 2, n_process
          displs(j) = displs(j-1) + n_glo(j-1)
       end do
    end if

    ! Gather data of different lengths onto rank 0
    call mpi_gatherv (data_loc, n_loc, MPI_DOUBLE_PRECISION, data, n_glo, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call mpi_gatherv (lat_loc,  n_loc, MPI_DOUBLE_PRECISION, lat,  n_glo, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call mpi_gatherv (lon_loc,  n_loc, MPI_DOUBLE_PRECISION, lon,  n_glo, displs, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

    ! Calculate spherical harmonics power spectrum on rank 0
    if (rank == 0) call spectrum
    deallocate (data, lat, lon, data_loc, lat_loc, lon_loc)
  end subroutine spec_sphere

  subroutine spectrum
!!$  Uses shtools routine SHExpandDH to expand data at irregularly spaced grid points on the sphere using a least squares inversion and then
!!$  the power spectrum is calculated using SHPowerSpectrum or SHPowerSpectrumDensity. For a given spherical harmonic degree l,
!!$
!!$ cilm : output, real(8), dimension (2, lmax+1, lmax+1)
!!$ The real spherical harmonic coefficients of the function. The coefficients C1lm and C2lm refer to the cosine (Clm) and sine (Slm) coefficients,
!!$ respectively, with Clm=cilm(1,l+1,m+1) and Slm=cilm(2,l+1,m+1).
!!$
!!$ lmax : input, integer
!!$ The maximum spherical harmonic bandwidth of the input grid, which is n/2-1.
!!$ If the optional parameter lmax_calc is not specified, this corresponds to the maximum spherical harmonic degree of the output coefficients cilm.
!!$ 
!!$ exitstatus: output, optional, integer
!!$ If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error.
!!$ 0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.
!!$
!!$ pspectrum : output, real(8), dimension (lmax+1)
!!$ pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2 
    use SHTOOLS
    implicit none
    integer                                 :: ierr, j, lmax
    real(8), dimension(:),      allocatable :: pspectrum
    real(8), dimension (:,:,:), allocatable :: cilm
    character(4)                            :: var_file

    ! Find spherical harmonics for irregularly spaced data
    lmax = floor (sqrt(dble(nmax)) - 1)
    write (6,'(a,i4)') "Maximum spherical harmonic degree for overdetermined least-squares inversion = ", lmax
    allocate (cilm(2,lmax+1,lmax+1), pspectrum(lmax+1))
    
    call SHExpandLSQ (cilm, data, lat, lon, nmax, lmax, norm=1, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHExpandLSQ"

    ! Calculate power spectrum
    call SHPowerSpectrum (cilm, lmax, pspectrum, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHPowerSpectrum"

    write (var_file, '(i4.4)') cp_idx
    open (unit=10, file=trim(run_id)//'_vort_spectrum_sphere_'//var_file, form="FORMATTED", status="REPLACE")
    do j = 1, lmax+1
       write (10,'(i4,1x,es10.4)') j, pspectrum(j) / (4*MATH_PI*radius**2)
       write (6,'(i4,1x,es10.4)')  j, pspectrum(j) / (4*MATH_PI*radius**2)
    end do
    close (10)
    deallocate (cilm, pspectrum)
  end subroutine spectrum

  subroutine define_data_triag (dom, i, j, zlev, offs, dims)
    ! Stores vorticity and associated latitude-longitude coordinates for spherical harmonic calculation
    ! (data on triangles)
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_t

    id = idx (i, j, offs, dims)

    idata_loc = idata_loc + 1
    id_t = TRIAG*id + LORT + 1
    data_loc(idata_loc) = dom%vort%elts(id_t)
    call cart2sph (dom%ccentre%elts(id_t), lon_loc(idata_loc), lat_loc(idata_loc))
    
    idata_loc = idata_loc + 1
    id_t = TRIAG*id + UPLT + 1
    data_loc(idata_loc) = dom%vort%elts(id_t)
    call cart2sph (dom%ccentre%elts(id_t), lon_loc(idata_loc), lat_loc(idata_loc))
  end subroutine define_data_triag

   subroutine define_data_hex (dom, i, j, zlev, offs, dims)
     ! Stores vorticity and associated latitude-longitude coordinates for spherical harmonic calculation
     ! (data on triangles)
    implicit none
    type (Domain)                  :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_t

    id = idx (i, j, offs, dims) + 1

    idata_loc = idata_loc + 1
    data_loc(idata_loc) = dom%press_lower%elts(id)
    call cart2sph (dom%node%elts(id), lon_loc(idata_loc), lat_loc(idata_loc))
  end subroutine define_data_hex

  subroutine initialize_lonlat
    implicit none

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export
    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
  end subroutine initialize_lonlat

   subroutine initialize_lonlat_2layer
    implicit none

    Nx = (/-N/2, N/2/)
    Ny = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export
    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
  end subroutine initialize_lonlat_2layer

  subroutine project_vorticity_onto_plane (l, default_val)
    ! Projects field from sphere at grid resolution l to longitude-latitude plane on grid defined by (Nx, Ny)
    use domain_mod
    use comm_mpi_mod
    integer :: l, itype
    real(8) :: default_val

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
end program spherical_harmonics
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




