program flat_projection_data
  ! Post-processing of checkpoint data to calculate flat projection
  use main_mod
  use test_case_mod
  use io_mod  
  implicit none
  
  integer                                :: nt, Ncumul
  integer, parameter                     :: nvar_save=6, nvar_zonal=8 ! Number of variables to save
  integer, dimension(2)                  :: Nx, Ny
  
  real(4), dimension(:,:),   allocatable :: field2d
  real(8)                                :: dx_export, dy_export, kx_export, ky_export
  real(8), dimension(2)                  :: lon_lat_range
  real(8), dimension(:,:,:), allocatable :: field2d_save, zonal_av
  
  character(2)                           :: var_file
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
     ref_press      = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = ref_press                   ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p

     u_0            = 35.0_8                      ! maximum velocity of zonal wind
     eta_0          = 0.252_8                     ! value of eta at reference level (level of the jet)
  elseif (trim (test_case) == "DCMIP2008c5") then
     compressible   = .true.                      ! Compressible equations

     radius         = 6.371229d6                  ! mean radius of the Earth in meters
     grav_accel     = 9.80616_8                   ! gravitational acceleration in meters per second squared
     omega          = 7.29211d-5                  ! Earth’s angular velocity in radians per second
     ref_press      = 100145.6_8                  ! reference pressure (mean surface pressure) in Pascals
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
     ref_press      = 1.0d5                       ! reference pressure (mean surface pressure) in Pascals
     ref_surf_press = ref_press                   ! reference surface pressure
     R_d            = 287.0_8                     ! ideal gas constant for dry air in joules per kilogram Kelvin
     gamma          = c_p/c_v                     ! heat capacity ratio
     kappa          = 2.0_8/7.0_8                 ! kappa=R_d/c_p
  else
     write (6,'(A)') "Test case not supported"
     stop
  end if
  resume = mean_beg
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
  ! Initialize variables
  call initialize (apply_initial_conditions, set_thresholds, dump, load, run_id)

  ! Initialize statistics
  call initialize_stat

  ! Calculate zonal average over all check points
  if (rank == 0) write (6,'(A)') "Calculating zonal averages over all checkpoints"
  Nt = 0
  do cp_idx = mean_beg, mean_end
     Nt = Nt + 1
     resume = NONE
     call restart (set_thresholds, load, run_id)
     call cal_zonal_av
     if (cp_idx == cp_2d) call latlon
  end do
  call barrier

  ! Finish covariance calculations
  zonal_av(:,:,2)   = zonal_av(:,:,2)   / (Ncumul-1)
  zonal_av(:,:,6:8) = zonal_av(:,:,6:8) / (Ncumul-1)
  
  if (rank==0) call write_out
  call finalize
contains
  subroutine cal_zonal_av
    ! Zonal average means and covariances over all checkpoints using stable online algorithm
    ! Uses Welford's stable onlne algorithm
    use domain_mod
    implicit none
    
    integer                                       :: d, ix, j, k
    real(8) ,dimension (Ny(1):Ny(2))              :: Tprime, Uprime, Vprime, Tprime_new, Uprime_new, Vprime_new
    real(8) ,dimension (Nx(1):Nx(2), Ny(1):Ny(2)) :: Tproj, Uproj, Vproj

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    ! Calculate temperature at all vertical levels (saved in exner_fun)
    call cal_surf_press (sol)
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE)

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
       
          ! Eddy kinetic energy (covariance)
          zonal_av(k,:,7) = zonal_av(k,:,7) + 0.5 * (Uprime * Uprime_new + Vprime * Vprime_new)

          ! Eddy heat flux (covariance)
          zonal_av(k,:,8) = zonal_av(k,:,8) + Vprime * Tprime_new
       end do
    end do
  end subroutine cal_zonal_av

  subroutine latlon
    ! Interpolate variables defined in valrange onto lon-lat grid of size (Nx(1):Nx(2), Ny(1):Ny(2), zlevels),
    use domain_mod
    integer :: d, j, k

    if (rank == 0) write (6,'(A,i6)') "Saving latitude-longitude projection of checkpoint file = ", cp_2d

    ! Fill up grid to level l and do inverse wavelet transform onto the uniform grid at level l
    call fill_up_grid_and_IWT (level_save)

    call cal_surf_press (sol)

    ! Remap to pressure_save vertical levels for saving data
    sol_save = sol(:,1:save_levels)
    call apply_onescale (interp_save, level_save, z_null, -1, 2)

    ! Calculate temperature at all vertical levels (saved in exner_fun) and temperature at interpolated saved vertical levels
    call apply_onescale (cal_temp, level_save, z_null, 0, 1)
    exner_fun%bdry_uptodate = .false.
    call update_vector_bdry (exner_fun, NONE)

    ! Latitude-longitude projections
    do k = 1, save_levels
       ! Temperature
       call project_onto_plane (horiz_flux(k), level_save, 0.0_8)
       field2d_save(:,:,1+k-1) = field2d

       ! Calculate zonal and meridional velocities and vorticity
       do d = 1, size(grid)
          velo => sol_save(S_VELO,k)%data(d)%elts
          vort => grid(d)%vort%elts
          do j = 1, grid(d)%lev(level_save)%length
             call apply_onescale_to_patch (interp_vel_hex, grid(d), grid(d)%lev(level_save)%elts(j), z_null,  0, 1)
             call apply_onescale_to_patch (cal_vort,       grid(d), grid(d)%lev(level_save)%elts(j), z_null, -1, 1)
          end do
          call apply_to_penta_d (post_vort, grid(d), level_save, z_null)
          nullify (velo, vort)
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
       call apply_onescale (vort_triag_to_hex, level_save, z_null, 0, 1)
       call project_vorticity_onto_plane (level_save, 1.0_8)
       field2d_save(:,:,5+k-1) = field2d

       ! Surface pressure
       call project_surf_press_onto_plane (level_save, 1.0_8)
       field2d_save(:,:,6+k-1) = field2d
    end do
  end subroutine latlon

  subroutine write_out
    ! Writes out results
    integer            :: i, k, v
    integer, parameter :: funit = 400

    ! 2d projections
    do v = 1, nvar_save*save_levels
       write (var_file, '(i1)') v
       open (unit=funit, file=trim(run_id)//'.3.0'//var_file)
       do i = Ny(1), Ny(2)
          write (funit,'(2047(E15.6, 1X))') field2d_save(:,i,v)
       end do
       close (funit)
    end do

    ! Zonal average of solution over all vertical levels
    do v = 1, nvar_zonal
       write (var_file, '(i2)') v+10
       open (unit=funit, file=trim(run_id)//'.3.'//var_file)
       do k = zlevels,1,-1
          write (funit,'(2047(E15.6, 1X))') zonal_av(k,:,v)
       end do
       close (funit)
    end do

    ! Coordinates

    ! Longitude values
    write (var_file, '(i2)') 20
    open (unit=funit, file=trim(run_id)//'.3.'//var_file) 
    write (funit,'(2047(E15.6, 1X))') (-180+dx_export*(i-1)/MATH_PI*180, i=1,Nx(2)-Nx(1)+1)
    close (funit)

    ! Latitude values
    write (var_file, '(i2)') 21
    open (unit=funit, file=trim(run_id)//'.3.'//var_file) 
    write (funit,'(2047(E15.6, 1X))') (-90+dy_export*(i-1)/MATH_PI*180, i=1,Ny(2)-Ny(1)+1)
    close (funit)

    ! Pressure vertical coordinates
    write (var_file, '(i2)') 22
    open (unit=funit, file=trim(run_id)//'.3.'//var_file) 
    write (funit,'(2047(E15.6, 1X))') (a_vert(k)*ref_press/ref_surf_press + b_vert(k), k=zlevels+1,2,-1)
    close (funit)

    ! Compress files
    command = 'ls -1 '//trim(run_id)//'.3.?? > tmp' 
    call system (command)
    command = 'tar czf '//trim(run_id)//'.3.tgz -T tmp --remove-files &'
    call system (command)
    deallocate (field2d, field2d_save, zonal_av)
  end subroutine write_out

  subroutine initialize_stat
    implicit none

    Nx     = (/-N/2, N/2/)
    Ny     = (/-N/4, N/4/)

    lon_lat_range = (/2*MATH_PI, MATH_PI/)
    dx_export = lon_lat_range(1)/(Nx(2)-Nx(1)+1); dy_export = lon_lat_range(2)/(Ny(2)-Ny(1)+1)
    kx_export = 1.0_8/dx_export; ky_export = 1.0_8/dy_export

    allocate (field2d(Nx(1):Nx(2),Ny(1):Ny(2)))
    allocate (zonal_av(1:zlevels,Ny(1):Ny(2),nvar_zonal))
    allocate (field2d_save(Nx(1):Nx(2),Ny(1):Ny(2),nvar_save*save_levels))
    zonal_av = 0.0_8
  end subroutine initialize_stat

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

                val   = grid(d)%adj_geopot%elts(id+1)
                valN  = grid(d)%adj_geopot%elts(idN+1)
                valE  = grid(d)%adj_geopot%elts(idE+1)
                valNE = grid(d)%adj_geopot%elts(idNE+1)

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

                val   = grid(d)%adj_mass%elts(id+1)
                valN  = grid(d)%adj_mass%elts(idN+1)
                valE  = grid(d)%adj_mass%elts(idE+1)
                valNE = grid(d)%adj_mass%elts(idNE+1)

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
    real(8) :: dpressure, pressure_lower, pressure_upper, spress

    d = dom%id + 1
    id = idx(i, j, offs, dims)

    spress = 0.0_8
    do k = 1, zlevels
       spress = spress + sol(S_MASS,k)%data(d)%elts(id+1)
    end do
    spress = press_infty + grav_accel*spress

    do kk = 1, save_levels
       ! Find pressure at current levels (not interfaces)
       pressure_lower = spress
       pressure_upper = 0.5*(a_vert(1)+a_vert(2))*ref_press + 0.5*(b_vert(1)+b_vert(2))*spress
       k = 1
       do while (pressure_upper > pressure_save(kk))
          k = k+1
          pressure_lower = pressure_upper
          pressure_upper = 0.5*(a_vert(k)+a_vert(k+1))*ref_press + 0.5*(b_vert(k)+b_vert(k+1))*spress
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

function grad_physics (scalar, dom, id, idE, idNE, idN, local_type)
  use domain_mod
  use test_case_mod
  implicit none

  real(8), dimension(1:EDGE)               :: grad_physics
  real(8), dimension(:), pointer           :: scalar
  type(Domain)                             :: dom
  integer                                  :: id, idE, idNE, idN
  logical                                  :: local_type

  if (.not.local_type) then ! Usual gradient at edges of hexagon E, NE, N
     grad_physics(RT+1) = (scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*id+RT+1) 
     grad_physics(DG+1) = (scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*id+DG+1) 
     grad_physics(UP+1) = (scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*id+UP+1) 
  else ! Gradient for southwest edges of hexagon W, SW, S
     grad_physics(RT+1) = -(scalar(idE+1) - scalar(id+1))  /dom%len%elts(EDGE*idE+RT+1) 
     grad_physics(DG+1) = -(scalar(id+1)  - scalar(idNE+1))/dom%len%elts(EDGE*idNE+DG+1)
     grad_physics(UP+1) = -(scalar(idN+1) - scalar(id+1))  /dom%len%elts(EDGE*idN+UP+1) 
  end if
end function grad_physics

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




