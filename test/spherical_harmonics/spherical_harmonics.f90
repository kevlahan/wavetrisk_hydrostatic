program spherical_harmonics
  ! Use SHTOOLS to compute global and local spherical harmonics power spectra from wavetrisk data interpolated to a uniform grid
  ! and projected to a uniform latitude longitude grid
  ! (see https://shtools.github.io/SHTOOLS/index.html for information about SHTOOLS)
  ! Wieczorek, M. A. and F. J. Simons 2007 Minimum-variance multitaper spectral estimation on the sphere.
  ! J. Fourier Anal. Appl., 13, doi:10.1007/s00041-006-6904-1, 665-692.
  !
  ! Only need velocity field for spectra (not scalars)
  use main_mod
  use test_case_mod
  use io_mod
  use projection_mod
  use io_vtk_mod
  implicit none

  integer                              :: idata_loc, nmax
  real(8), dimension(:),   allocatable :: data, data_loc, lat_loc, lon_loc
  character(2)                         :: var_file
  character(130)                       :: command

  call init_arch_mod 
  call init_comm_mpi_mod
  call read_test_case_parameters

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !     Parameters (need radius for spectra)
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (trim (data_case) == "climate") then
     compressible            = .true.                    
     split_mean_perturbation = .true.           
     physics_model           = .true.
  elseif (trim (data_case) == "drake") then
     vert_diffuse            = .true.
     mode_split              = .true.                    
     penalize                = .true.                    
     compressible            = .false.                   
     split_mean_perturbation = .true.

     radius                  = radius/6                     
  elseif (trim (data_case) == "jet") then
     mode_split              = .true.                    
     compressible            = .false.                   
     penalize                = .true.                    
     tke_closure             = .true.
     
     radius                  = 10000 * KM               
  end if
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !    Initialization
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  call assign_functions
  call initialize (run_id)
  call print_test_case_parameters

  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Compute spectra for each checkpoint
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  do cp_idx = cp_beg, cp_end
     resume = NONE
     call restart
     
     if (trim (spec_type) == 'sphere') then
        call spec_sphere
     elseif (trim (spec_type) == 'latlon') then
        call spec_latlon_1layer
     end if
  end do
  
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !      Compute and save averages
  ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  if (cp_end /= cp_beg .and. rank == 0) then
     call avg_spec ('curlu')
     call avg_spec ('divu')
     call avg_spec ('u')
     if (local_spec) then
        call avg_spec ('curlu_local')
        call avg_spec ('divu_local')
        call avg_spec ('u_local')
     end if
  end if
  
  if (rank == 0) call cleanup
  
  call finalize
contains
  subroutine avg_spec (data_type)
    implicit none
    character(*) :: data_type

    integer                            :: cp, j, jj, k, lmax
    real(8), dimension(:), allocatable :: pspec, pspec_av
    character(4)                       :: var_file1, var_file2

    lmax = N/4 - 1

    allocate (pspec(lmax+1), pspec_av(lmax+1))

    do k = k_min, k_max
       pspec = 0d0; pspec_av = 0d0
       write (var_file2, '(i4.4)') k
       do cp = cp_beg, cp_end
          write (var_file1, '(i4.4)') cp
          open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//var_file2//'_'//trim(data_type)//'_spec', &
               form="FORMATTED", status="OLD")
          do j = 1, lmax + 1
             read (10,*) jj, pspec(j)
          end do
          close (10)
          pspec_av = pspec_av + pspec
       end do
       pspec_av = pspec_av / (cp_end - cp_beg + 1)

       open (unit=10, file=trim(run_id)//'_'//var_file2//'_'//trim(data_type)//'_spec', &
            form="FORMATTED", status="REPLACE")
       do j = 1, lmax + 1
          write (10,'(i4,1x,es10.4)') j, pspec_av(j)
       end do
       close (10)
    end do

    deallocate (pspec, pspec_av)
  end subroutine avg_spec

  subroutine avg_local_spec (data_type)
    implicit none
    character(*) :: data_type

    integer                            :: cp, j, jj, k, lmax
    real(8), dimension(:), allocatable :: mtse, mtse_av, sd, sd_av
    character(4)                       :: var_file1, var_file2

    lmax = N/4 - 1
    
    allocate (mtse(lmax-lwin+1),    sd(lmax-lwin+1))
    allocate (mtse_av(lmax-lwin+1), sd_av(lmax-lwin+1))
    
    do k = k_min, k_max
       write (var_file2, '(i4.4)') k
       mtse = 0d0; mtse_av = 0d0; sd = 0d0; sd_av = 0d0
       do cp = cp_beg, cp_end
          write (var_file1, '(i4.4)') cp
          open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//var_file2//'_'//trim(data_type)//'_spec', &
               form="FORMATTED", status="OLD")
          do j = 1, lmax - lwin + 1
             read (10,*) jj, mtse(j), sd(j)
          end do
          close (10)
          mtse_av = mtse_av + mtse
          sd_av = sd_av + sd
       end do
       mtse_av = mtse_av / (cp_end - cp_beg + 1)
       sd_av   = sd_av   / (cp_end - cp_beg + 1)

       open (unit=10, file=trim(run_id)//'_'//var_file2//'_'//trim(data_type)//'_spec', &
            form="FORMATTED", status="REPLACE")
       do j = 1, lmax - lwin + 1
          write (10,'(i4,1x,2(es10.4,1x))') j, mtse_av(j), sd_av(j)
       end do
       close (10)
    end do

    deallocate (mtse, mtse_av, sd, sd_av)
  end subroutine avg_local_spec

  subroutine spec_latlon_1layer
    ! Compute energy spectra of div-free and curl-free parts of the velocity field from 2d latitude-longitude projections
    ! of vorticity and div(u) respectively.
    use domain_mod
    use multi_level_mod
    implicit none
    integer :: d, i, j, k, l

    call initialize_projection (N)

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_fill
    call fill_up_grid_and_IWT (l)
    do k = 1, zlevels
       if (rank == 0) write (6,'(A,i3,A,i6,A,f10.2,A)') "Energy spectrum of vertical layer ", k, &
            " at checkpoint ", cp_idx, " at ", time/DAY, " days"
      
       ! Set mean on filled grid
       call apply_onescale (init_mean, l, k, -BDRY_THICKNESS, BDRY_THICKNESS)

       ! Vorticity at hexagon points
       do d = 1, size(grid)
          velo => sol(S_VELO,k)%data(d)%elts
          vort => grid(d)%vort%elts
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
       field2d = 0d0
       call project_array_onto_plane ("press_lower", l, 1d0)
       if (rank == 0) call spectrum_lon_lat ("curlu", k)

       ! Divergence at hexagon points
       call cal_divu_ml (sol(S_VELO,k))
       field2d = 0d0
       call project_array_onto_plane ("divu", l, 1d0)
       if (rank == 0) call spectrum_lon_lat ("divu", k)

       ! Velocity magnitude at hexagon points
       call umag (sol(S_VELO,k))
       field2d = 0d0
       call project_array_onto_plane ("ke", l, 1d0)
       if (rank == 0) call spectrum_lon_lat ("u", k)
    end do
    deallocate (field2d)
  end subroutine spec_latlon_1layer

  subroutine spectrum_lon_lat (data_type, k)
    use SHTOOLS
    implicit none
    !  Input: 
    !  A 2D equally sampled (n,n) grid  (default) or equally spaced grid (n,2*n) that conforms to the sampling theorem of Driscoll and Healy (1994).
    !  The first latitudinal band corresponds to 90 N, the latitudinal band for 90 S is not included, and the latitudinal sampling interval is 180/n degrees.
    !  The first longitudinal band is 0 E, the longitude band for 360 E is not included, and the longitudinal sampling interval is 360/n
    !  for an equally sampled grid and 180/n for an equally spaced grid, respectively.
    !    n (integer): size of data griddh (n,2*n) or (n,n)
    !
    !  Output: 
    !  cilm (real(8), dimension (2,n/2,n/2) or (2, lmax_calc+1, lmax_calc+1)):
    !  The real spherical harmonic coefficients of the function. These will be exact if the function is bandlimited to degree lmax=n/2-1.
    !  The coefficients c1lm and c2lm refer to the cosine (clm) and sine (slm) coefficients, respectively, with clm=cilm(1,l+1,m+1) and slm=cilm(2,l+1,m+1).
    !
    !  lmax (integer):     
    !  The maximum spherical harmonic bandwidth of the input grid, which is n/2-1.
    !  If the optional parameter lmax_calc is not specified, this corresponds to the maximum spherical harmonic degree of the output coefficients cilm.
    integer                                 :: k
    character(*)                            :: data_type

    integer                                 :: ierr, i, j, lmax
    real(8)                                 :: area
    real(8), dimension (:,:,:), allocatable :: cilm      ! spherical harmonic coefficients
    real(8), dimension (:),     allocatable :: pspectrum ! global power spectrum of the function
    character(4)                            :: var_file1, var_file2
    
    write (var_file1, '(i4.4)') cp_idx
    write (var_file2, '(i4.4)') k
   
    lmax = N/4 - 1
    allocate (cilm(2, lmax+1, lmax+1))

    ! Expand latitude-longitude data in spherical harmonics (reorder elements in latitude from N to S)
    call SHExpandDH (transpose(dble(field2d(Nx(1):Nx(2)-1,Ny(2):Ny(1)+1:-1))), N/2, cilm, lmax, norm=1, sampling=2, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHExpandDH"

    ! Calculate power spectrum
    allocate (pspectrum(lmax+1))
    call SHPowerSpectrum (cilm, lmax, pspectrum, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHPowerSpectrum"

    if (k /= 0) then
       open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//var_file2//'_'//trim(data_type)//'_spec', &
            form="FORMATTED", status="REPLACE")
    else
       open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//trim(data_type)//'_spec', &
            form="FORMATTED", status="REPLACE")
    end if
    area = 4d0*MATH_PI*radius**2
    do j = 1, lmax + 1
       write (10,'(i4,1x,es10.4)') j, pspectrum(j) * area
    end do
    close (10)

    if (local_spec) then
       if (k /= 0) then
          open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//var_file2//'_'//trim(data_type)//'_local_spec', &
               form="FORMATTED", status="REPLACE")
       else
          open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//trim(data_type)//'_local_spec', &
               form="FORMATTED", status="REPLACE")
       end if

       ! Determine the spherical-harmonic bandwidth needed to achieve specified concentration factor
       lwin = SHFindLWin (theta0*DEG, angular_order, concentration)

       ! Compute local spectrum and associated error
       call local_spectrum (cilm, lmax)
    end if
    deallocate (cilm, pspectrum)
  end subroutine spectrum_lon_lat
  
  subroutine local_spectrum (cilm, lmax)
    ! Computes a localized multitaper spectral analysis of an input function expressed in spherical harmonics. The maximum degree of the
    ! localized multitaper cross-power spectrum estimate is lmax-lmaxt. lat0 (deg), lon0 (deg), theta0 (radians) define local region
    use SHTOOLS
    implicit none
    integer                                :: lmax
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
    area = 2d0*MATH_PI*radius**2 * (1d0 - cos(theta0*DEG))
    do j = 1, lmax - lwin + 1
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
    integer                       :: d, j, k, l, n_loc
    integer, dimension(n_process) :: displs, n_glo
    logical, parameter            :: hex = .false.

    ! Fill up grid to level l and inverse wavelet transform onto the uniform grid at level l
    l = level_fill
    call fill_up_grid_and_IWT (l)

    do k = k_min, k_max
       if (rank == 0) write (6,'(A,i3,A,i6,A,f10.2,A)') "Energy spectrum of vertical layer ", k, &
            " at checkpoint ", cp_idx, " at ", time/DAY, " days"
       ! Vorticity on triangles
       do d = 1, size(grid)
          velo  => sol(S_VELO,k)%data(d)%elts
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
       allocate (data_loc(n_loc), lat_loc(n_loc), lon_loc(n_loc)); data_loc = 0d0; lat_loc = 0d0; lon_loc = 0d0
       allocate (data(nmax), lat(nmax), lon(nmax));  data = 0d0; lat = 0d0; lon = 0d0

       ! Store data for spherical harmonics transform
       idata_loc = 0

       if (hex) then
          if (rank==0) call apply_to_pole (define_data_hex_pole, min_level-1, z_null, 1, .false.)
          call apply_onescale (define_data_hex, l, z_null, 0, 0)
       else
          if (rank==0) call apply_to_pole (define_data_tri_pole, min_level-1, z_null, 1, .false.)
          call apply_onescale (define_data_tri, l, z_null, 0, 0)
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
       call mpi_gatherv (data_loc, n_loc, MPI_DOUBLE_PRECISION, data, n_glo, displs, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
       call mpi_gatherv (lat_loc,  n_loc, MPI_DOUBLE_PRECISION, lat,  n_glo, displs, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
       call mpi_gatherv (lon_loc,  n_loc, MPI_DOUBLE_PRECISION, lon,  n_glo, displs, &
            MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)

       ! Calculate spherical harmonics power spectrum on rank 0
       if (rank == 0) call spectrum_sphere ("vort_lsq", k)
    end do
    deallocate (data, lat, lon, data_loc, lat_loc, lon_loc)
  end subroutine spec_sphere

  subroutine spectrum_sphere (data_type, k)
    !  Uses shtools routine SHExpandDH to expand data at irregularly spaced grid points on the sphere using a least squares inversion and then
    !  the power spectrum is calculated using SHPowerSpectrum or SHPowerSpectrumDensity. For a given spherical harmonic degree l,
    !
    !  cilm : output, real(8), dimension (2, lmax+1, lmax+1)
    !  The real spherical harmonic coefficients of the function. The coefficients C1lm and C2lm refer to the cosine (Clm) and sine (Slm) coefficients,
    !  respectively, with Clm=cilm(1,l+1,m+1) and Slm=cilm(2,l+1,m+1).
    !
    !  lmax : input, integer
    !  The maximum spherical harmonic bandwidth of the input grid, which is n/2-1.
    !  If the optional parameter lmax_calc is not specified, this corresponds to the maximum spherical harmonic degree of the output coefficients cilm.
    ! 
    !  exitstatus: output, optional, integer
    !  If present, instead of executing a STOP when an error is encountered, the variable exitstatus will be returned describing the error.
    !  0 = No errors; 1 = Improper dimensions of input array; 2 = Improper bounds for input variable; 3 = Error allocating memory; 4 = File IO error.
    !
    !    pspectrum : output, real(8), dimension (lmax+1)
    !    pspectrum(l) = Sum_{i=1}^2 Sum_{m=0}^l cilm(i, l+1, m+1)**2
    use SHTOOLS
    implicit none
    integer      :: k
    character(*) :: data_type

    integer                                 :: ierr, j, lmax
    real(8), dimension(:),      allocatable :: pspectrum
    real(8), dimension (:,:,:), allocatable :: cilm
    character(4)                            :: var_file1, var_file2

    ! Find spherical harmonics for irregularly spaced data
    lmax = floor (sqrt(dble(nmax)) - 1)
    write (6,'(a,i4)') "Maximum spherical harmonic degree for overdetermined least-squares inversion = ", lmax
    allocate (cilm(2,lmax+1,lmax+1), pspectrum(lmax+1))
    
    call SHExpandLSQ (cilm, data, lat, lon, nmax, lmax, norm=1, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHExpandLSQ"

    ! Calculate power spectrum
    call SHPowerSpectrum (cilm, lmax, pspectrum, exitstatus=ierr)
    if (ierr /= 0) write(6,'(a,i1,a)') "Error status = ", ierr, " for routine SHPowerSpectrum"

    write (var_file1, '(i4.4)') cp_idx
    write (var_file1, '(i4.4)') k
    open (unit=10, file=trim(run_id)//'_vort_spectrum_sphere_'//var_file1//'_'//var_file2, &
         form="FORMATTED", status="REPLACE")
    do j = 1, lmax+1
       write (10,'(i4,1x,es10.4)') j, pspectrum(j) / (4*MATH_PI*radius**2)
       write (6,'(i4,1x,es10.4)')  j, pspectrum(j) / (4*MATH_PI*radius**2)
    end do
    close (10)

    if (local_spec) then
       open (unit=10, file=trim(run_id)//'_'//var_file1//'_'//var_file2//'_'//trim(data_type)//'_local_spec', &
            form="FORMATTED", status="REPLACE")

       ! Determine the spherical-harmonic bandwidth needed to achieve specified concentration factor
       lwin = SHFindLWin (theta0*DEG, angular_order, concentration)

       ! Compute local spectrum and associated error
       call local_spectrum (cilm, lmax)
    end if
    deallocate (cilm, pspectrum)
  end subroutine spectrum_sphere

  subroutine define_data_tri_pole (dom, p, i, j, zlev, offs, dims, ival)
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: ival, i, j, p, zlev
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
  end subroutine define_data_tri_pole

  subroutine define_data_tri (dom, i, j, zlev, offs, dims)
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
  end subroutine define_data_tri

  subroutine define_data_hex_pole (dom, p, i, j, zlev, offs, dims, ival)
    use domain_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: ival, i, j, p, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id, id_t

    id = idx (i, j, offs, dims) + 1

    idata_loc = idata_loc + 1
    data_loc(idata_loc) = dom%press_lower%elts(id)
    call cart2sph (dom%node%elts(id), lon_loc(idata_loc), lat_loc(idata_loc))
  end subroutine define_data_hex_pole

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

  subroutine cleanup
    implicit none

    ! Delete intermediate files
    command = '\rm -f '//trim(run_id)//'_[0-9][0-9][0-9][0-9]_[0-9][0-9][0-9][0-9]*spec'
    call system (trim(command))

    ! Compress output
    command = 'bash -c "ls -1 '//trim(run_id)//'*_spec > '//trim(run_id)//'_tmp1"'
    call system (trim(command))
    command = 'bash -c "gtar caf '//trim(run_id)//'_spec.tgz -T '//trim(run_id)//'_tmp1 --remove-files"'
    call system (trim(command))
    command = '\rm -f '//trim(run_id)//'_tmp1'
    call system (trim(command))
  end subroutine cleanup
end program





