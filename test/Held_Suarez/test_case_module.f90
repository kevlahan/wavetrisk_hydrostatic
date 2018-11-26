module test_case_mod
  ! Module file for Held & Suarez (1994) test case
  use shared_mod
  use domain_mod
  use comm_mpi_mod
  implicit none
  
  ! Standard variables
  integer                              :: CP_EVERY, save_zlev
  real(8)                              :: dt_cfl, initotalmass, mass_error, tau_diffusion, totalmass, total_cpu_time
  real(8)                              :: dPdim, Hdim, Ldim, Pdim, R_ddim, specvoldim, Tdim, Tempdim, dTempdim, Udim
  real(8), allocatable, dimension(:,:) :: threshold_def

  ! Test case variables
  real(8) :: delta_T, delta_theta, sigma_b, sigma_c, k_a, k_f, k_s, T_0, T_mean, T_tropo
  real(8) :: delta_T2, sigma_t, sigma_v, sigma_0, gamma_T, u_0
contains
  subroutine cal_theta_eq (p, p_s, lat, theta_equil, k_T)
    ! Returns equilibrium potential temperature theta_equil and Newton cooling constant k_T
    use domain_mod
    implicit none
    real(8) :: p, p_s, lat, theta_equil, k_T

    real(8) :: sn2, cs2, sigma, theta_force, theta_tropo

    sn2 = sin (lat)**2
    cs2 = cos (lat)**2

    sigma = (p - p_top) / (p_s - p_top)

    k_T = k_a + (k_s-k_a) * max (0.0_8, (sigma-sigma_b)/sigma_c) * cs2**2

    theta_tropo = T_tropo * (p/p_0)**(-kappa) ! Potential temperature at tropopause

    theta_force = T_mean - delta_T*sn2 - delta_theta*cs2 * log (p/p_0)

    theta_equil = max (theta_tropo, theta_force) ! Equilibrium temperature
  end subroutine cal_theta_eq

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! From Jablonowski and Williamson (2006) without perturbation
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    type(Coord) :: x_i, x_E, x_N, x_NE
    integer     :: id, d, idN, idE, idNE
    real(8)     :: p, p_s, pot_temp, sigma

    d = dom%id+1

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    ! Surface pressure
    dom%surf_press%elts(id+1) = surf_pressure (x_i)
    p_s = dom%surf_press%elts(id+1)

    ! Pressure at level zlev
    p = 0.5 * (a_vert(zlev)+a_vert(zlev+1) + (b_vert(zlev)+b_vert(zlev+1))*p_s)

    ! Normalized pressure
    sigma = (p - p_top) / (p_s - p_top)
    sigma_v = (sigma - sigma_0) * MATH_PI/2

    ! Mass/Area = rho*dz at level zlev
    sol(S_MASS,zlev)%data(d)%elts(id+1) = a_vert_mass(zlev) + b_vert_mass(zlev)*p_s/grav_accel
    
    ! Potential temperature
    pot_temp = set_temp (x_i, sigma) * (p/p_0)**(-kappa)
!     call cal_theta_eq (p, p_s, lat, theta_equil, k_T)

    ! Mass-weighted potential temperature
    sol(S_TEMP,zlev)%data(d)%elts(id+1) = sol(S_MASS,zlev)%data(d)%elts(id+1) * pot_temp

    ! Set initial velocity field
    call vel2uvw (dom, i, j, zlev, offs, dims, vel_fun)
  end subroutine init_sol

  real(8) function set_temp (x_i, sigma)
    ! From Jablonowski and Williamson (2006)
    implicit none
    type(Coord) :: x_i
    real(8)     :: sigma

    real(8) :: cs2, lon, lat, sn2, Tmean

    call cart2sph (x_i, lon, lat)
    sn2 = sin (lat)**2
    cs2 = cos (lat)**2

    if (sigma >= sigma_t) then
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel)
    else
       Tmean = T_0*sigma**(R_d*Gamma_T/grav_accel) + delta_T * (sigma_t - sigma)**5
    end if

    set_temp = Tmean + 0.75_8 * sigma*MATH_PI*u_0/R_d * sin(sigma_v) * sqrt(cos(sigma_v)) * &
         (2*u_0*cos(sigma_v)**1.5*(-2*sn2**3*(cs2+1/3.0_8) + 10/63.0_8) + radius*omega*(8/5.0_8*cs2**1.5*(sn2+2/3.0_8) - MATH_PI/4))
  end function set_temp

  real(8) function surf_geopot (x_i)
    ! Surface geopotential from Jablonowski and Williamson (2006)
    implicit none
    Type(Coord) :: x_i
    real(8)     :: c1, cs2, lon, lat, sn2

    ! Find latitude and longitude from Cartesian coordinates
    call cart2sph (x_i, lon, lat)
    cs2 = cos (lat)**2
    sn2 = sin (lat)**2

    c1 = u_0*cos((1.0_8-sigma_0)*MATH_PI/2)**1.5
    surf_geopot =  c1*(c1*(-2*sn2**3*(cs2 + 1/3.0_8) + 10/63.0_8)  + radius*omega*(8/5.0_8*cs2**1.5*(sn2 + 2/3.0_8) - MATH_PI/4))
!    surf_geopot = 0.0_8 ! Uniform
  end function surf_geopot

  real(8) function surf_pressure (x_i)
    ! Surface pressure
    implicit none
    type(Coord) :: x_i

    surf_pressure = p_0
  end function surf_pressure

  subroutine vel_fun (lon, lat, u, v)
    ! Zonal latitude-dependent wind
    implicit none
    real(8) :: lon, lat, u, v

    real(8) :: r

    call random_number (r)
    
    u = u_0 * cos (sigma_v)**1.5 * sin (2*lat)**2 + r ! Zonal velocity component
!    u = 0.0_8 ! Uniform
    v = 0.0_8                                 ! Meridional velocity component
  end subroutine vel_fun
 
  subroutine set_thresholds
    ! Set thresholds dynamically (trend or sol must be known)
    use wavelet_mod
    implicit none
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: lnorm, threshold_new

    character(3), parameter :: order = "inf"

    if (default_thresholds) then ! Initialize once
       threshold_new = threshold_def
    else
       if (adapt_trend) then
          call cal_lnorm (trend, order, lnorm)
       else
          call cal_lnorm (sol,   order, lnorm)
       end if
       threshold_new = max (tol*lnorm, threshold_def) ! Avoid very small thresholds before instability develops
    end if
    !threshold = 0.1*threshold_new + 0.9*threshold
    threshold = threshold_new
  end subroutine set_thresholds

   subroutine initialize_a_b_vert
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (uniform) then
       do k = 1, zlevels+1
          a_vert(k) = dble(k-1)/dble(zlevels) * p_top
          b_vert(k) = 1.0_8 - dble(k-1)/dble(zlevels)
       end do
    else
       ! LMDZ grid
       call cal_AB
    end if
    
    ! Set pressure at infinity
    p_top = a_vert(zlevels+1) ! note that b_vert at top level is 0, a_vert is small but non-zero

    ! Set mass coefficients
    a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
    b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
  end subroutine initialize_a_b_vert

  subroutine cal_AB
    ! Computes A and B coefficients for hybrid vertical grid as in LMDZ
    implicit none

    integer                         :: l
    real(8)                         :: snorm
    real(8), dimension(1:zlevels)   :: dsig
    real(8), dimension(1:zlevels+1) :: sig

    snorm  = 0.0_8
    do l = 1, zlevels
       dsig(l) = 1.0_8 + 7 * sin (MATH_PI*(l-0.5_8)/(zlevels+1))**2
       snorm = snorm + dsig(l)
    end do

    do l = 1, zlevels
       dsig(l) = dsig(l)/snorm
    end do

    sig(zlevels+1) = 0.0_8
    do l = zlevels, 1, -1
       sig(l) = sig(l+1) + dsig(l)
    end do

    b_vert(zlevels+1) = 0.0_8
    do  l = 1, zlevels
       b_vert(l) = exp (1.0_8 - 1/sig(l)**2)
       a_vert(l) = (sig(l) - b_vert(l)) * p_0
    end do
    b_vert(1) = 1.0_8
    a_vert(1) = 0.0_8
    a_vert(zlevels+1) = (sig(zlevels+1) - b_vert(zlevels+1)) * p_0
  end subroutine cal_AB

  subroutine read_test_case_parameters
    implicit none
    integer, parameter :: fid = 500
    real(8)            :: press_save
    character(255)     :: filename, varname

    ! Find input parameters file name
    if (iargc() >= 1) then
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
    read (fid,*) varname, min_allowed_mass
    read (fid,*) varname, adapt_trend
    read (fid,*) varname, default_thresholds
    read (fid,*) varname, perfect
    read (fid,*) varname, tol
    read (fid,*) varname, optimize_grid
    read (fid,*) varname, adapt_dt
    read (fid,*) varname, cfl_num
    read (fid,*) varname, press_save
    read (fid,*) varname, Laplace_order_init
    read (fid,*) varname, n_diffuse
    read (fid,*) varname, tau_diffusion
    read (fid,*) varname, dt_write
    read (fid,*) varname, CP_EVERY
    read (fid,*) varname, time_end
    read (fid,*) varname, resume
    close(fid)
    
    allocate (pressure_save(1))
    pressure_save(1) = 1.0d2*press_save
    tau_diffusion = tau_diffusion * HOUR
    dt_write = dt_write * MINUTE
    time_end = time_end * HOUR
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
       write (6,'(A,es10.4)') "min_allowed_mass    = ", min_allowed_mass
       write (6,'(A,L1)')     "adapt_trend         = ", adapt_trend
       write (6,'(A,L1)')     "default_thresholds  = ", default_thresholds
       write (6,'(A,L1)')     "perfect             = ", perfect
       write (6,'(A,es10.4)') "tolerance           = ", tol
       write (6,'(A,i1)')     "optimize_grid       = ", optimize_grid
       write (6,'(A,L1)')     "adapt_dt            = ", adapt_dt
       write (6,'(A,es10.4)') "cfl_num             = ", cfl_num
       write (6,'(A,es10.4)') "pressure_save (hPa) = ", pressure_save(1)/100
       write (6,'(A,i1)')     "Laplace_order       = ", Laplace_order_init
       write (6,'(A,i4)')     "n_diffuse           = ", n_diffuse
       write (6,'(A,es10.4)') "tau_diffusion (h)   = ", tau_diffusion/HOUR
       write (6,'(A,es10.4)') "dt_write (min)      = ", dt_write/MINUTE
       write (6,'(A,i6)')     "CP_EVERY            = ", CP_EVERY
       write (6,'(A,es10.4)') "time_end (h)        = ", time_end/HOUR
       write (6,'(A,i6)')     "resume              = ", resume
       
       write (6,'(/,A)')      "STANDARD PARAMETERS"
       write (6,'(A,es10.4)') "radius              = ", radius
       write (6,'(A,es10.4)') "omega               = ", omega
       write (6,'(A,es10.4)') "p_0   (hPa)         = ", p_0/100
       write (6,'(A,es10.4)') "p_top (hPa)         = ", p_top/100
       write (6,'(A,es10.4)') "R_d                 = ", R_d
       write (6,'(A,es10.4)') "c_p                 = ", c_p
       write (6,'(A,es10.4)') "c_v                 = ", c_v
       write (6,'(A,es10.4)') "gamma               = ", gamma
       write (6,'(A,es10.4)') "kappa               = ", kappa

       write (6,'(/,A)')      "TEST CASE PARAMETERS"
       write (6,'(A,es10.4)') "T_0                 = ", T_0
       write (6,'(A,es10.4)') "T_mean              = ", T_mean
       write (6,'(A,es10.4)') "T_tropo             = ", T_tropo
       write (6,'(A,es10.4)') "sigma_b             = ", sigma_b
       write (6,'(A,es10.4)') "k_a                 = ", k_a
       write (6,'(A,es10.4)') "k_f                 = ", k_f
       write (6,'(A,es10.4)') "k_s                 = ", k_s
       write (6,'(A,es10.4)') "delta_T             = ", delta_T
       write (6,'(A,es10.4)') "delta_theta         = ", delta_theta
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
            'time [h] = ', time/HOUR, &
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
            time/HOUR, dt, sum (threshold(S_MASS,:))/zlevels, sum (threshold(S_TEMP,:))/zlevels, &
            sum (threshold(S_VELO,:))/zlevels, level_end, sum (n_active), min_mass, mass_error, rel_imbalance, timing
    end if
  end subroutine print_log

  subroutine initialize_thresholds
    ! Set default thresholds based on dimensional scalings of norms
    implicit none

    integer :: k
    real(8), dimension(S_MASS:S_VELO,1:zlevels) :: lnorm
    
    allocate (threshold(S_MASS:S_VELO,1:zlevels));     threshold     = 0.0_8
    allocate (threshold_def(S_MASS:S_VELO,1:zlevels)); threshold_def = 0.0_8

    lnorm(S_MASS,:) = dPdim/grav_accel
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
    use wavelet_mod
    implicit none
    
    integer               :: k
    real(8)               :: k_max, visc
    real(8), dimension(3) :: L_scaled

    allocate (viscosity_divu(1:zlevels))
    
    ! Smallest edge length (scaled to account for non-uniform mesh)
    dx_min = 0.9 * sqrt (4*MATH_PI*radius**2/(sqrt(3.0_8)/2*10*4**max_level))

    ! Largest wavenumber on regular lozenge grid
    k_max = MATH_PI/(sqrt(3.0_8)*dx_min)

    ! CFL limit for time step
    dt_cfl = cfl_num*dx_min/(wave_speed+Udim)
    dt_init = dt_cfl

    ! Viscosity constant from eigenvalues of Laplacian
    if (Laplace_order_init == 0) then
       viscosity_mass = 0.0_8
       viscosity_temp = 0.0_8
       viscosity_divu = 0.0_8
       viscosity_rotu = 0.0_8
    elseif (Laplace_order_init == 1 .or. Laplace_order_init == 2) then
       L_scaled = L_diffusion / 2**(max_level-min_level) ! Correct length scales for finest grid

       viscosity_mass = L_scaled(1)**(2*Laplace_order_init) / tau_diffusion * n_diffuse
       viscosity_temp = L_scaled(1)**(2*Laplace_order_init) / tau_diffusion * n_diffuse
       viscosity_divu = L_scaled(2)**(2*Laplace_order_init) / tau_diffusion * n_diffuse
       viscosity_rotu = L_scaled(3)**(2*Laplace_order_init) / tau_diffusion * n_diffuse
    elseif (Laplace_order_init > 2) then
       if (rank == 0) write (6,'(A)') 'Unsupported iterated Laplacian (only 0, 1 or 2 supported)'
       stop
    end if
    visc = max (viscosity_mass, viscosity_temp, maxval (viscosity_divu), viscosity_rotu)
    
    if (rank == 0) then
       write (6,'(3(A,es8.2),/)') "dx_min  = ", dx_min, " k_max  = ", k_max, " dt_cfl = ", dt_cfl
       write (6,'(4(A,es8.2))') "Viscosity_mass = ", viscosity_mass/n_diffuse, " Viscosity_temp = ", viscosity_temp/n_diffuse, &
            " Viscosity_divu = ", sum (viscosity_divu)/zlevels/n_diffuse, " Viscosity_rotu = ", viscosity_rotu/n_diffuse
       write (6,'(A,es8.2,A)') "Diffusion stability constant = ", dt_cfl/dx_min**(2*Laplace_order_init) * visc, &
            " (should be < 0.25, or about < 0.5 for RK45 ssp if CFL <= 1.2)"
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

    id   = idx(i,   j,   offs, dims)
    idN  = idx(i,   j+1, offs, dims)
    idE  = idx(i+1, j,   offs, dims)
    idNE = idx(i+1, j+1, offs, dims)

    x_i  = dom%node%elts(id+1)
    x_E  = dom%node%elts(idE+1)
    x_N  = dom%node%elts(idN+1)
    x_NE = dom%node%elts(idNE+1)

    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+RT+1) = proj_vel(vel_fun, x_i,  x_E)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+DG+1) = proj_vel(vel_fun, x_NE, x_i)
    sol(S_VELO,zlev)%data(d)%elts(EDGE*id+UP+1) = proj_vel(vel_fun, x_i,  x_N)
  end subroutine vel2uvw

  subroutine set_save_level
    ! Determines closest vertical level to desired pressure
    implicit none
    integer :: k
    real(8) :: dpress, p, save_press

    dpress = 1d16; save_zlev = 0
    do k = 1, zlevels
       p = 0.5 * (a_vert(k)+a_vert(k+1) + (b_vert(k)+b_vert(k+1))*p_0)
       if (abs(p-pressure_save(1)) < dpress) then
          dpress = abs (p-pressure_save(1))
          save_zlev = k
          save_press = p
       end if
    end do
    if (rank==0) write (6,'(/,A,i2,A,f5.1,A,/)') "Saving vertical level ", save_zlev, &
         " (approximate pressure = ", save_press/100, " hPa)"
  end subroutine set_save_level

  subroutine initialize_seed
    implicit none
    
    integer                            :: k
    integer, dimension(1:8)            :: values
    integer, dimension(:), allocatable :: seed

    call date_and_time (values=values)
    call random_seed(size=k)
    allocate (seed(1:k))
    seed = values(8)
    call random_seed (put=seed)
  end subroutine initialize_seed

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
