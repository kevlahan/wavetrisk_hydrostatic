Module test_case_mod
  ! Module file for DCMIP2008c5
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none

  ! Standard variables
  integer                              :: CP_EVERY, resume_init, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Local variables
  real(8)                              :: dH, mean_depth
contains
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer     :: d, id_i
    real (8)    :: depth, dz, z, z_bathy
    type(Coord) :: x_i
    
    d    = dom%id+1
    id_i = idx (i, j, offs, dims) + 1
    x_i  = dom%node%elts(id+1)

    ! Height of local bathymetry
    z_bathy = surf_geopot (x_i) / grav_accel 

    ! Local depth
    depth = mean_depth + eta (x_i) - z_bathy

    ! Vertical grid size
    dz = depth / zlevels

    ! Local z-coordinate (mid layer)
    z = z_bathy + (zlev - 0.5_8) * dz 

    ! Equally spaced sigma coordinates in z: sol(S_MASS) = ref_density * dz
    sol(S_MASS,zlev)%data(d)%elts(id_i) = ref_density * dz

    ! Density * dz
    sol(S_TEMP,zlev)%data(d)%elts(id_i) = sol(S_MASS,zlev)%data(d)%elts(id_i) * density (x_i, z)/ref_density

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  real(8) function surf_geopot (x_i)
    ! Surface geopotential
    implicit none
    type(Coord) :: x_i
    
    surf_geopot = 0.0_8
  end function surf_geopot

  real(8) function eta (x_i)
    ! Free surface perturbation
    implicit none
    type(Coord) :: x_i

    real(8) :: lon, lat, rgrc
    
    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    rgrc = radius * acos(sin(lat_c)*sin(lat) + cos(lat_c)*cos(lat)*cos(lon-lon_c))

    eta = dH * exp__flush (-rgrc**2/d2)
  end function eta

  real(8) function density (x_i, z)
    ! Density profile
    implicit none
    real(8)     :: z
    type(Coord) :: x_i
    
    density = ref_density
  end function density

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    u = 0.0_8 ! Zonal velocity component
    v = 0.0_8 ! Meridional velocity component
  end subroutine vel_fun
  
  subroutine set_thresholds
    ! Set thresholds dynamically (trend or sol must be known)
    use lnorms_mod
    use wavelet_mod
    implicit none
    integer                                     :: v
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: threshold_new
    character(3), parameter                     :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       if (adapt_trend) then
          call cal_lnorm_trend (trend, order)
       else
          call cal_lnorm_sol (sol,   order)
       end if
       threshold_new = tol*lnorm
    end if

    if (istep >= 10) then
       threshold = 0.01*threshold_new + 0.99*threshold
    else
       threshold = threshold_new
    end if
  end subroutine set_thresholds

  subroutine initialize_a_b_vert
    implicit none

    allocate (a_vert(1:zlevels+1), b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
  end subroutine initialize_a_b_vert

   subroutine read_test_case_parameters
    implicit none
    integer, parameter :: fid = 500
    real(8)            :: press_save
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
    read (fid,*) varname, compressible
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, uniform
    read (fid,*) varname, remap
    read (fid,*) varname, remapscalar_type
    read (fid,*) varname, remapvelo_type
    read (fid,*) varname, iremap
    read (fid,*) varname, adapt_trend
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, perfect
    read (fid,*) varname, tol
    read (fid,*) varname, optimize_grid
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, cfl_num
    read (fid,*) varname, timeint_type
    read (fid,*) varname, press_save
    read (fid,*) varname, Laplace_order_init
    read (fid,*) varname, n_diffuse
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, rebalance
    read (fid,*) varname, time_end
    read (fid,*) varname, resume_init
    close(fid)
    
    allocate (pressure_save(1))
    pressure_save(1) = 1.0d2*press_save
    dt_write = dt_write * DAY
    time_end = time_end * DAY
    resume   = resume_init
    Laplace_order = Laplace_order_init
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none
    
     if (rank==0) then
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A)')        "RUN PARAMETERS"
       write (6,'(A,A)')      "test_case           = ", trim (test_case)
       write (6,'(A,A)')      "run_id              = ", trim (run_id)
       write (6,'(A,L1)')     "compressible        = ", compressible
       write (6,'(A,i3)')     "min_level           = ", min_level
       write (6,'(A,i3)')     "max_level           = ", max_level
       write (6,'(A,i5)')     "number of domains   = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL        = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL         = ", PATCH_LEVEL
       write (6,'(A,i3)')     "zlevels             = ", zlevels
       write (6,'(A,L1)')     "uniform             = ", uniform
       write (6,'(A,L1)')     "remap               = ", remap
       write (6,'(a,a)')      "remapscalar_type    = ", trim (remapscalar_type)
       write (6,'(a,a)')      "remapvelo_type      = ", trim (remapvelo_type)
       write (6,'(a,i3)')     "iremap              = ", iremap
       write (6,'(A,L1)')     "adapt_trend         = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds  = ", default_thresholds
       write (6,'(A,L1)')     "perfect             = ", perfect
       write (6,'(A,es10.4)') "tolerance           = ", tol
       write (6,'(A,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt            = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(a,a)')      "timeint_type        = ", trim (timeint_type)
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", pressure_save(1)/100
       write (6,'(A,i1)')     "Laplace_order       = ", Laplace_order_init
       write (6,'(A,i4)')     "n_diffuse           = ", n_diffuse
       write (6,'(A,es10.4)') "dt_write (min)      = ", dt_write/DAY
       write (6,'(A,i6)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(a,l1)')     "rebalance           = ", rebalance
       write (6,'(A,es10.4)') "time_end (h)        = ", time_end/DAY
       write (6,'(A,i6)')     "resume              = ", resume_init
       
       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius              = ", radius
       write (6,'(A,es10.4)') "omega               = ", omega
       write (6,'(A,es10.4)') "p_top (hPa)         = ", p_top/100

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es10.4)') "d2                  = ", d2
       write (6,'(A,es10.4)') "dH                  = ", dH
       write (6,'(A,es10.4)') "mean_depth          = ", mean_depth
       write (6,'(A,es10.4)') "lon_c               = ", lon_c
       write (6,'(A,es10.4)') "lat_c               = ", lat_c

       write (6,'(A)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine print_log
    ! Prints out and saves logged data to a file
    implicit none

    integer :: min_load, max_load
    real(8) :: avg_load, rel_imbalance, timing

    timing = get_timing(); total_cpu_time = total_cpu_time + timing

    call cal_load_balance (min_load, avg_load, max_load, rel_imbalance)
    
    if (rank == 0) then
       write (6,'(a,es12.6,4(a,es8.2),a,i2,a,i9,4(a,es8.2,1x))') &
            'time [d] = ', time/DAY, &
            ' dt [s] = ', dt, &
            '  mass tol = ', sum (threshold(S_MASS,:))/zlevels, &
            ' temp tol = ', sum (threshold(S_TEMP,:))/zlevels, &
            ' velo tol = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax = ', level_end, &
            ' dof = ', sum (n_active), &
            ' min rel mass = ', min_mass, &
            ' mass error = ', mass_error, &
            ' balance = ', rel_imbalance, &
            ' cpu = ', timing

       write (12,'(5(es15.9,1x),i2,1x,i9,1x,4(es15.9,1x))')  &
            time/DAY, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    integer :: k
    
    allocate (threshold(S_MASS:S_VELO,1:zlevels));     threshold     = 0.0_8
    allocate (threshold_def(S_MASS:S_VELO,1:zlevels)); threshold_def = 0.0_8

    lnorm(S_MASS,:) = dH
    do k = 1, zlevels
       lnorm(S_TEMP,k) = (a_vert_mass(k) + b_vert_mass(k)*Pdim/grav_accel)*dTempdim
    end do
    lnorm(S_TEMP,:) = lnorm(S_TEMP,:) + Tempdim*lnorm(S_MASS,:) ! Add component due to tendency in mass 
    lnorm(S_VELO,:) = Udim
    if (adapt_trend) lnorm = lnorm/Tdim
    threshold_def = tol * lnorm
  end subroutine initialize_thresholds
  
  subroutine initialize_dt_viscosity 
    ! Initializes viscosity
    implicit none
    real(8) :: area, C_divu, C_sclr, C_rotu, tau_divu, tau_rotu, tau_sclr

    area = 4*MATH_PI*radius**2/(20*4**max_level) ! average area of a triangle
    dx_min = sqrt (4/sqrt(3.0_8) * area)         ! edge length of average triangle
      
    ! Diffusion constants
    C_sclr = 1.5d-2       ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
    C_divu = 1.5d-2    ! <= 1.75e-2 for hyperdiffusion (lower than exact limit 1/6^2 = 2.8e-2 due to non-uniform grid)
    C_rotu = C_sclr / 4**Laplace_order_init ! <= 1.09e-3 for hyperdiffusion (lower than exact limit 1/24^2 = 1.7e-3 due to non-uniform grid)
    
    ! CFL limit for time step
    dt_cfl = cfl_num*dx_min/(wave_speed+Udim) * 0.85 ! corrected for dynamic value
    dt_init = dt_cfl

    tau_sclr = dt_cfl / C_sclr
    tau_divu = dt_cfl / C_divu
    tau_rotu = dt_cfl / C_rotu

    if (Laplace_order_init == 0) then
       visc_sclr = 0.0_8
       visc_divu = 0.0_8
       visc_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       visc_sclr = dx_min**(2*Laplace_order_init) / tau_sclr
       visc_rotu = dx_min**(2*Laplace_order_init) / tau_rotu
       visc_divu = dx_min**(2*Laplace_order_init) / tau_divu
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if

    if (rank == 0) then
       write (6,'(/,3(a,es8.2),a,/)') "dx_min  = ", dx_min/1d3, " [km] dt_cfl = ", dt_cfl, " [s] tau_sclr = ", tau_sclr/HOUR, " [h]"
       write (6,'(3(a,es8.2),/)') "C_sclr = ", C_sclr, "  C_divu = ", C_divu, "  C_rotu = ", C_rotu
       write (6,'(4(a,es8.2))') "Viscosity_mass = ", visc_sclr(S_MASS)/n_diffuse, &
          " Viscosity_temp = ", visc_sclr(S_TEMP)/n_diffuse, &
          " Viscosity_divu = ", visc_divu/n_diffuse, " Viscosity_rotu = ", visc_rotu/n_diffuse
    end if
  end subroutine initialize_dt_viscosity

  subroutine apply_initial_conditions
    implicit none
    integer :: k, l

    do l = level_start, level_end
       do k = 1, zlevels
          call apply_onescale (init_sol, l, k, -1, 1)
       end do
    end do
  end subroutine apply_initial_conditions

  subroutine vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
    ! Sets the velocities on the computational grid given a function vel_fun that provides zonal and meridional velocities
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
    external                        :: vel_fun

    integer      :: d, id, idE, idN, idNE
    type (Coord) :: x_i, x_E, x_N, x_NE

    d = dom%id+1

    id   = idx (i,   j,   offs, dims)
    idN  = idx (i,   j+1, offs, dims)
    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, x_i,  x_E)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel(vel_fun, x_NE, x_i)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, x_i,  x_N)
  end subroutine vel2uvw

  subroutine set_save_level
    ! Determines closest vertical level to desired height z
    implicit none
    integer :: k
    real(8) :: dz, p, save_press

    dz = mean_depth / zlevels; save_zlev = 0
    do k = 1, zlevels
       z = (zlev - 0.5_8) * dz 
       if (abs(z-pressure_save(1)) < dz) then
          save_zlev = k
          save_press = z
       end if
    end do
    if (rank==0) write (6,'(/,A,i2,A,f5.1,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate height = ", save_press, " [m])"
  end subroutine set_save_level

  subroutine dump (fid)
    implicit none
    integer :: fid

    write (fid) itime
    write (fid) iwrite
    write (fid) threshold
  end subroutine dump

  subroutine load (fid)
    implicit none
    integer :: fid

    read (fid) itime
    read (fid) iwrite
    read (fid) threshold
  end subroutine load
end module test_case_mod
