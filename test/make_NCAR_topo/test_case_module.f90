module test_case_mod
  ! Module file for Held & Suarez (1994) test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use std_atm_profile_mod
  use io_mod
  implicit none
  integer         :: nsmth
  real(8)         :: Area_max, Area_min, dt_nu, smth_scl
  character(9999) :: topo_data
  logical         :: analytic_topo = .false.
contains
  subroutine assign_functions
    ! Assigns generic pointer functions to functions defined in test cases
    implicit none
    apply_initial_conditions => apply_initial_conditions_case
    initialize_a_b_vert      => initialize_a_b_vert_case
    initialize_dt_viscosity  => initialize_dt_viscosity_case
    set_save_level           => set_save_level_case
    initialize_thresholds    => initialize_thresholds_case
    set_thresholds           => set_thresholds_case
    update                   => update_case
    dump                     => dump_case
    load                     => load_case
  end subroutine assign_functions

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! From Jablonowski and Williamson (2006) without perturbation
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: d, id, k
    
    d = dom%id+1
    id = idx (i, j, offs, dims)

    do k = 1, zlevels
       sol(S_MASS,k)%data(d)%elts(id+1)                      = 0d0
       sol(S_TEMP,k)%data(d)%elts(id+1)                      = 0d0
       sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do
  end subroutine init_sol
  
  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! From Jablonowski and Williamson (2006) without perturbation
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: d, id, k

    d  = dom%id+1
    id = idx (i, j, offs, dims)

    do k = 1, zlevels
       sol_mean(S_MASS,k)%data(d)%elts(id+1)                      = 0d0
       sol_mean(S_TEMP,k)%data(d)%elts(id+1)                      = 0d0
       sol_mean(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1) = 0d0
    end do
  end subroutine init_mean

  subroutine init_topo (dom, i, j, zlev, offs, dims)
    ! Assigns analytic topography
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id
    real(8) :: lat, lon, width

    d  = dom%id + 1
    id = idx (i, j, offs, dims) + 1
    
    width = 8d0 * dx_min
    call cart2sph (grid(d)%node%elts(id), lon, lat)

    topography%data(d)%elts(id) = &
         tanh_profile (4d3*METRE,  65*DEG,   95*DEG,  25*DEG,  40*DEG) + & ! Himalayas
         tanh_profile (4d3*METRE, -60*DEG,  -50*DEG, -60*DEG, -10*DEG)     ! Andes
  contains
    real(8) function tanh_profile (height, lon_min, lon_max, lat_min, lat_max)
      implicit none
      real(8), intent(in) :: height, lon_min, lon_max, lat_min, lat_max

      tanh_profile = height * (profile1d (lat, lat_min, lat_max) * profile1d (lon, lon_min, lon_max))
    end function tanh_profile

    real(8) function profile1d (x, xmin, xmax)
      implicit none
      real(8) :: x, xmin, xmax

      profile1d = prof (x, xmax) - prof (x, xmin)
    end function profile1d

    real(8) function prof (x, x0)
      implicit none
      real(8) :: x, x0

      prof = 0.5d0 * (1d0 - tanh ((x - x0)/((width/5d0)/radius)))
    end function prof
    
    real(8) function ellipse_profile (lon_0, lat_0, e, height, sigma, theta)
      ! Elliptical smoothed mountain, ellipticity 0 < e <= 1
      ! Minimum resolution of gradient = N
      implicit none
      real(8), intent(in) :: lon_0, lat_0 ! centre of ellipse (in radians)
      real(8), intent(in) :: e            ! ellipticity
      real(8), intent(in) :: height       ! height of ellipse (in metres)
      real(8), intent(in) :: sigma        ! size of ellipse (in radians)
      real(8), intent(in) :: theta        ! orientation (in radians)

      real(8)            :: dtheta, p, lat_loc, lat_rot, lon_loc, lon_rot, rsq, sigma_x, sigma_y

      real(8), parameter :: npts_slope = 5d0 ! resolve slope with this many cells

      dtheta = dx_min / radius

      sigma_x = sigma
      sigma_y = sigma_x * sqrt (1d0 - e**2)

      ! Transform coordinates (shift and rotate)
      lon_loc = lon - lon_0; lat_loc = lat - lat_0

      lon_rot = lon_loc * cos (theta) - lat_loc * sin (theta)
      lat_rot = lon_loc * sin (theta) + lat_loc * cos (theta)

      rsq = (lon_rot/sigma_x)**2 + (lat_rot/sigma_y)**2

      ! Order of hyper Gaussian is 2 p
      p = (log (0.01d0) / log (1d0 - npts_slope * dtheta/(2d0*sigma_y))) / 2d0

      ellipse_profile = height * exp__flush (-rsq**p)
    end function ellipse_profile
  end subroutine init_topo

  subroutine set_thresholds_case
    use lnorms_mod
    use wavelet_mod
    implicit none

    threshold = tol
  end subroutine set_thresholds_case

  subroutine initialize_a_b_vert_case
    implicit none
    integer :: k

    ! Allocate vertical grid parameters
    allocate (a_vert(1:zlevels+1),    b_vert(1:zlevels+1))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    do k = 1, zlevels+1
       a_vert(k) = dble(k-1)/dble(zlevels) * p_top
       b_vert(k) = 1d0 - dble(k-1)/dble(zlevels)
    end do
    a_vert = a_vert(zlevels+1:1:-1) * p_0
    b_vert = b_vert(zlevels+1:1:-1)

    ! Set mass coefficients
    a_vert_mass = (a_vert(1:zlevels) - a_vert(2:zlevels+1))/grav_accel
    b_vert_mass =  b_vert(1:zlevels) - b_vert(2:zlevels+1)
  end subroutine initialize_a_b_vert_case

  subroutine read_test_case_parameters
    implicit none
    integer, parameter :: fid = 500
    character(255)     :: filename, varname
    character(2)       :: var_file
    logical            :: file_exists

    ! Find input parameters file name
    if (command_argument_count () >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'test_case.in'
    end if
    if (rank == 0) write (6,'(a,A)') "Input file = ", trim (filename)

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, max_level
    read (fid,*) varname, topo_data
    read (fid,*) varname, smth_scl
    close(fid)
  end subroutine read_test_case_parameters

  subroutine print_test_case_parameters
    implicit none
    if (rank==0) then
       write (6,'(a)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(a)')        "RUN PARAMETERS"
       write (6,'(a,i3)')     "min_level            = ", min_level
       write (6,'(a,i3)')     "max_level            = ", max_level
       write (6,'(a,i5)')     "number of domains    = ", N_GLO_DOMAIN
       write (6,'(a,i5)')     "number of processors = ", n_process
       write (6,'(a,i5)')     "DOMAIN_LEVEL         = ", DOMAIN_LEVEL
       write (6,'(a,i5)')     "PATCH_LEVEL          = ", PATCH_LEVEL
       write (6,'(a,es10.4)') "radius               = ", radius
       write (6,'(a,es10.4)') "grav_accel           = ", grav_accel
       write (6,'(a,a)')      "topo_data            = ", trim (topo_data)
       write (6,'(a,a)')      "topo_file            = ", trim (topo_file)
       write (6,'(a,es10.4)') "smth_scl             = ", smth_scl
       write (6,'(a)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine initialize_thresholds_case
    implicit none
    threshold_def = tol
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case
  end subroutine initialize_dt_viscosity_case

  subroutine apply_initial_conditions_case
    implicit none

    call apply_bdry (init_mean, z_null, 0, 1)
    call apply_bdry (init_sol,  z_null, 0, 1)
    sol%bdry_uptodate      = .false.
    sol_mean%bdry_uptodate = .false.

    call update_array_bdry (sol,        NONE)
    call update_array_bdry (sol_mean,   NONE)
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Not needed in this test case
    implicit none
   
  end subroutine update_case

  subroutine set_save_level_case
    implicit none
  end subroutine set_save_level_case

   subroutine dump_case (fid)
    implicit none
    integer :: fid
  end subroutine dump_case

  subroutine load_case (fid)
    implicit none
    integer :: fid
  end subroutine load_case

  subroutine smooth_topo (l)
    ! Smooths topography by taking an Euler time step with diffusive trend
    implicit none
    integer :: l

    real(8) :: Area

    Area = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**l) ! average hexagonal cell area
    dt_nu = Area/3d0 / dble (nsmth)                  ! amount of smoothing
    do istep = 1, nsmth
       call cal_trend_topo (l)

       call apply_onescale (Euler_step_topo, l, z_null, 0, 1)

       topography%bdry_uptodate = .false.
       call update_bdry (topography, l)
    end do
  end subroutine smooth_topo

  subroutine Euler_step_topo (dom, i, j, zlev, offs, dims)
    ! Euler step for topography smoothing
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id

    id = idx (i, j, offs, dims) + 1
    d = dom%id + 1
   
    topography%data(d)%elts(id) = topography%data(d)%elts(id) + trend(S_MASS,1)%data(d)%elts(id) 
  end subroutine Euler_step_topo

  subroutine cal_trend_topo (l)
    ! Computes Laplacian of topography, div (grad (topography) ), stored as trend(S_MASS,1)
    ! (used to smooth topography at scale l)
    implicit none
    integer :: l
    
    integer :: d, j

    call zero_float (trend(S_MASS,1))

    call update_bdry (topography, l)

    ! Compute scalar fluxes
    do d = 1, size(grid)
       scalar => topography%data(d)%elts
       h_flux => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call step1 (dom=grid(d), p=grid(d)%lev(l)%elts(j), itype=1)
       end do
       nullify (scalar, h_flux)
    end do
    horiz_flux%bdry_uptodate = .false.
    call update_bdry (horiz_flux(S_MASS), l)

    ! Compute divergence of fluxes
    do d = 1, size(grid)
       dscalar => trend(S_MASS,1)%data(d)%elts
       h_flux  => horiz_flux(S_MASS)%data(d)%elts
       do j = 1, grid(d)%lev(l)%length
          call apply_onescale_to_patch (cal_div_topo, grid(d), grid(d)%lev(l)%elts(j), z_null, 0, 1)
       end do
       nullify (dscalar, h_flux)
    end do
    trend(S_MASS,1)%bdry_uptodate = .false.
    call update_bdry (trend(S_MASS,1), l)
  end subroutine cal_trend_topo

  subroutine cal_div_topo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, id, idE, idNE, idN, idW, idS, idSW

    d = dom%id + 1

    id   = idx (i,   j,   offs, dims)

    idE  = idx (i+1, j,   offs, dims)
    idNE = idx (i+1, j+1, offs, dims)
    idN  = idx (i,   j+1, offs, dims)

    idW  = idx (i-1, j,   offs, dims)
    idSW = idx (i-1, j-1, offs, dims)
    idS  = idx (i,   j-1, offs, dims)

    dscalar(id+1) = ( nu_dt (id, idE)  * h_flux(EDGE*id+RT+1)   - nu_dt (id, idW)  * h_flux(EDGE*idW+RT+1)  &
         + nu_dt (id, idSW) * h_flux(EDGE*idSW+DG+1) - nu_dt (id, idNE) * h_flux(EDGE*id+DG+1)   &
         + nu_dt (id, idN)  * h_flux(EDGE*id+UP+1)   - nu_dt (id, idS)  * h_flux(EDGE*idS+UP+1) ) * dom%areas%elts(id+1)%hex_inv
  contains
    real(8) function nu_dt (id1, id2)
      implicit none
      integer :: id1, id2

      real(8)            :: Area, fac, p1, p2, rx_0, s
      real(8), parameter :: r_max = 0.15d0

      call std_surf_pres (topography%data(d)%elts(id1+1), p1)
      call std_surf_pres (topography%data(d)%elts(id2+1), p2)

      rx_0 = abs (p1 - p2) / (p1 + p2) ! rx_0 at edge

      Area = 1d0 / dom%areas%elts(id+1)%hex_inv

      s = r_max/10d0
      fac = ( 2d0 + tanh((rx_0-r_max/2d0)/s) - tanh( - (rx_0-r_max/2d0)/s) ) / 4d0 

      nu_dt = fac * Area/3d0
    end function nu_dt
  end subroutine cal_div_topo

  subroutine cal_grad_topo (dom, i, j, zlev, offs, dims)
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: id
    
    id = idx (i,   j,   offs, dims)
    velo(EDGE*id+RT+1:EDGE*id+UP+1) = gradi_e (scalar, dom, i, j, offs, dims)

  end subroutine cal_grad_topo
  
end module test_case_mod
