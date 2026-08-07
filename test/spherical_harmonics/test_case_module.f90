module test_case_mod

#ifdef MPI
  use mpi_f08
#endif

  use kind_mod
  use shared_mod
  use arch_mod
  use comm_mpi_mod
  use domain_mod
  use domain_ops_mod
  use geom_mod
  use init_mod
  use io_mod
  use utils_mod
  use vert_diffusion_mod

  implicit none

  integer        :: angular_order, cp_beg, cp_end, lwin, N, ntaper, k_min, k_max
  real(dp)        :: concentration, lat0, lon0, theta0, z_linear
  character(255) :: data_case, spec_type
  character(10)  :: stratification = "tanh"
  character(255) :: coords         = "uniform" ! not used if sigma_z = .true.
  logical        :: local_spec
  
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
    set_thresholds           => set_thresholds_case
    surf_geopot              => surf_geopot_case
    update                   => update_case
    z_coords                 => z_coords_case
  end subroutine assign_functions

  
  function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type) result(flux)
    ! Scalar diffusion flux
    !
    ! NOTE: call with arguments (d, id, idW, idSW, idS, type) if type = .true. to compute gradient at soutwest edges W, SW, S
    use domain_mod

    implicit none

    type(Float_Field), intent(inout)        :: q(1:N_VARIABLE,1:zlevels)
    type(domain),      intent(inout)        :: dom
    integer,           intent(in)           :: id, idE, idNE, idN, v, zlev
    logical,           intent(in), optional :: type
    real(dp)                                :: flux(EDGE)

    flux = 0.0_dp
  end function physics_scalar_flux_case
  

  function physics_velo_source_case (dom, i, j, zlev, offs, dims) result(source)
    ! Additional physics for the source term of the velocity trend
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)
    real(dp)                    :: source(EDGE)

    source = 0.0_dp
  end function physics_velo_source_case

  
  function surf_geopot_case (d, id) result(val)
    implicit none
    integer, intent(in) :: d, id
    real(dp) :: val

    val = 0.0_dp
  end function surf_geopot_case
  

  subroutine read_test_case_parameters
    implicit none
    integer        :: fid = 500
    character(255) :: filename, varname

    ! Find input parameters file name
    if (iargc() >= 1) then
       CALL getarg (1, filename)
    else
       filename = 'spherical_harmonics.in'
    end if
    if (rank == 0) write (6,'(/,a,a,/)') "Input file = ", trim (filename)

    test_case = "spherical_harmonics"

    open (unit=fid, file=filename, action='READ')
    read (fid,*) varname, spec_type
    read (fid,*) varname, data_case
    read (fid,*) varname, physics_type
    read (fid,*) varname, Nsoil
    read (fid,*) varname, run_id
    read (fid,*) varname, cp_beg
    read (fid,*) varname, cp_end
    read (fid,*) varname, max_level
    read (fid,*) varname, zlevels
    read (fid,*) varname, level_fill
    read (fid,*) varname, N
    read (fid,*) varname, local_spec
    read (fid,*) varname, lat0
    read (fid,*) varname, lon0
    read (fid,*) varname, theta0
    read (fid,*) varname, concentration
    read (fid,*) varname, ntaper
    read (fid,*) varname, angular_order
    close(fid)

    resume = cp_beg

    k_min = 1 ; k_max = zlevels
  end subroutine read_test_case_parameters

  
  subroutine print_test_case_parameters
    implicit none

    if (rank==0) then
       write (6,'(a)') &
            '***********************************************************************&
            &************************************************************'
       write (6,'(/,a,/)')        "SPHERICAL HARMONICS SPECTRUM PARAMETERS"
       write (6,'(a,a)')      "spec_type               = ", trim (spec_type)
       write (6,'(a,a)')      "run_id                  = ", trim (run_id)
       write (6,'(a,a)')      "data_case               = ", trim (data_case)
       write (6,'(a,i3)')     "min_level               = ", min_level
       write (6,'(a,i3)')     "max_level               = ", max_level
       write (6,'(a,i3)')     "zlevels                 = ", zlevels
       if (trim(data_case) == "climate") then
          write (6,'(a,a)')   "physics_type            = ", trim (physics_type)
          if (trim(physics_type) == "Simple") &
               write (6,'(a,i0)')  "Nsoil                  = ", Nsoil
          write (6,'(a,l1)')   "NCAR_topo              = ", NCAR_topo
          if (NCAR_topo) then
             write (6,'(a,a)') "topo_file              = ", trim (topo_file)
             write (6,'(a,i3)') "topo_min_level         = ", topo_min_level
             write (6,'(a,i3)') "topo_max_level         = ", topo_max_level
          end if
       end if
       write (6,'(a,i5)')     "number of domains       = ", N_GLO_DOMAIN
       write (6,'(a,i5)')     "number of processors    = ", n_process
       write (6,'(a,i5)')     "DOMAIN_LEVEL            = ", DOMAIN_LEVEL
       write (6,'(a,i5)')     "PATCH_LEVEL             = ", PATCH_LEVEL
       write (6,'(a,i5)')     "BDRY_THICKNESS          = ", BDRY_THICKNESS
       write (6,'(a,l1)')     "split_mean_perturbation = ", split_mean_perturbation
       write (6,'(a,l1)')     "mode_split              = ", mode_split
       write (6,'(a,i4)')     "First checkpoint        = ", cp_beg
       write (6,'(a,i4)')     "Last checkpoint         = ", cp_end
       write (6,'(a,l1)')     "vert_diffuse            = ", vert_diffuse
       write (6,'(a,i3)')     "level_fill              = ", level_fill
       write (6,'(a,i5)')     "N                       = ", N
       write (6,'(a,l1)')     "local_spec              = ", local_spec
       write (6,'(a,f4.1)')   "lat0                    = ", lat0
       write (6,'(a,f4.1)')   "lon0                    = ", lon0
       write (6,'(a,f4.1)')   "theta0                  = ", theta0
       write (6,'(a,f4.2)') "concentration           = ", concentration
       write (6,'(a,i3)')     "ntaper                  = ", ntaper
       write (6,'(a,i3)')     "angular_order           = ", angular_order
       write (6,*) ' '
    end if
  end subroutine print_test_case_parameters

  
  subroutine apply_initial_conditions_case
    implicit none

    ! Initialize variables
    call apply_bdry (init_mean, z_null, 0, 1)
    call apply_bdry (init_sol,  z_null, 0, 1)
    sol%bdry_uptodate      = .false.
    sol_mean%bdry_uptodate = .false.

    call update_bdry (sol,      NONE) 
    call update_bdry (sol_mean, NONE)

    if (NCAR_topo) call apply_bdry (assign_NCAR_topo, z_null, 0, 1)
    topography%bdry_uptodate = .false.
    call update_bdry (topography, NONE)
  end subroutine apply_initial_conditions_case

  
  subroutine update_case
    implicit none

    if (NCAR_topo) call apply_bdry (assign_NCAR_topo, z_null, 0, 1)
    topography%bdry_uptodate = .false.
    call update_bdry (topography, NONE)
  end subroutine update_case

  
  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Dummy routine
      
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer :: d, id, id_i, k

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1

    do k = 1, zlevels
       sol(S_MASS,k)%data(d)%elts(id_i)        = 0.0_dp
       sol(S_TEMP,k)%data(d)%elts(id_i)        = 0.0_dp
       sol(S_VELO,k)%data(d)%elts(id_edge(id)) = 0.0_dp
    end do

    if (mode_split) then
       sol(S_MASS,zlevels+1)%data(d)%elts(id_i)        = 0.0_dp
       sol(S_TEMP,zlevels+1)%data(d)%elts(id_i)        = 0.0_dp
       sol(S_VELO,zlevels+1)%data(d)%elts(id_edge(id)) = 0.0_dp
    end if
  end subroutine init_sol

  
  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values
    ! In split mean perturbation need to set sol_mean(S_MASS,:)
    
    implicit none
    
    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1)
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                        :: d, id, id_i, k
    real(dp)                       :: eta, rho, rho_dz, z_s
    real(dp), dimension(1:zlevels) :: dz
    real(dp), dimension(0:zlevels) :: z
   
    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
  
    if (trim (data_case) == "drake") then ! include mean component or rho_dz
       eta = 0.0_dp    ! no free surface perturbation
       z_s = max_depth ! no bathymetry

       if (sigma_z) then
          z = z_coords_case (eta, z_s)
       else
          z = a_vert * eta + b_vert * z_s
       end if
       dz = z(1:zlevels) - z(0:zlevels-1)

       do k = 1, zlevels
          rho = ref_density ! do not include penalization
          rho_dz = rho * dz(k)

          sol_mean(S_MASS,k)%data(d)%elts(id_i)        = rho_dz ! needed for Qperp in spectral flux computation
          sol_mean(S_TEMP,k)%data(d)%elts(id_i)        = 0.0_dp !! Assume not needed !!
          sol_mean(S_VELO,k)%data(d)%elts(id_edge(id)) = 0.0_dp
       end do
    else
       do k = 1, zlevels
          sol_mean(S_MASS,k)%data(d)%elts(id_i)        = 0.0_dp 
          sol_mean(S_TEMP,k)%data(d)%elts(id_i)        = 0.0_dp
          sol_mean(S_VELO,k)%data(d)%elts(id_edge(id)) = 0.0_dp
       end do
    end if
   
    if (mode_split) then
       sol_mean(S_MASS,zlevels+1)%data(d)%elts(id_i)        = 0.0_dp
       sol_mean(S_TEMP,zlevels+1)%data(d)%elts(id_i)        = 0.0_dp
       sol_mean(S_VELO,zlevels+1)%data(d)%elts(id_edge(id)) = 0.0_dp
    end if
  end subroutine init_mean

  
  function z_coords_case (eta_surf, z_s) result(val)
    ! Hybrid sigma-z vertical coordinates to minimize inclination of layers to geopotential
    ! near the free surface over strong bathymetry gradients.
    ! Reference: similar to Shchepetkin and McWilliams (JCP vol 228, 8985-9000, 2009)
    !
    ! Sets the a_vert parameter that depends on eta_surf (but not b_vert).
    
    implicit none
    
    real(dp), intent (in) :: eta_surf, z_s ! free surface and bathymetry
    real(dp)              :: val(0:zlevels)

    integer                        :: k
    real(dp)                       :: cff, cff1, cff2, hc, z_0
    real(dp)                       :: Cs(0:zlevels), sc(0:zlevels)

    real(dp), parameter            :: theta_b = 0.0_dp, theta_s = 7.0_dp
    real(dp), parameter            :: hc_min = -200 * METRE ! minimum depth of uniform layer region

    select case (stratification)
    case ("linear")
       hc = abs (min (z_mixed, hc_min))
    case ("tanh")
       hc = abs (min (z_linear, hc_min))
    case default
       hc = abs (min (z_linear, hc_min))
    end select
    hc = hc_min
    
    cff1 = 1.0_dp / sinh (theta_s)
    cff2 = 0.5 / tanh (0.5 * theta_s)
    
    sc(0) = -1.0_dp
    Cs(0) = -1.0_dp
    cff = 1.0_dp / dble(zlevels)
    do k = 1, zlevels
       sc(k) = cff * dble (k - zlevels)
       Cs(k) = (1.0_dp - theta_b) * cff1 * sinh (theta_s * sc(k)) + theta_b * (cff2 * tanh (theta_s * (sc(k) + 0.5_dp)) - 0.5_dp)
    end do

    val(0) = z_s
    do k = 1, zlevels
       cff = hc * (sc(k) - Cs(k))
       z_0 = cff - Cs(k) * z_s
       a_vert(k) = 1.0_dp - z_0 / z_s
      val(k) = eta_surf * a_vert(k) + z_0
    end do
  end function z_coords_case
  

  subroutine initialize_a_b_vert_case
    ! Initialize hybrid sigma-coordinate vertical grid
    implicit none
    integer :: k

    allocate (a_vert(0:zlevels), b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))

    if (zlevels == 1) then
       a_vert(0) = 0.0_dp; a_vert(1) = 1.0_dp
       b_vert(0) = 1.0_dp; b_vert(1) = 0.0_dp
    elseif (zlevels == 2) then 
       a_vert(0) = 0.0_dp; a_vert(1) = 0.0_dp;            a_vert(2) = 1.0_dp
       b_vert(0) = 1.0_dp; b_vert(1) = z_mixed/max_depth; b_vert(2) = 0.0_dp
    elseif (zlevels >= 3) then
       if (trim (coords) == "chebyshev") then
          do k = 0, zlevels
             b_vert(k) = (1.0_dp + cos (dble(k)/dble(zlevels) * MATH_PI)) / 2
          end do
       elseif (trim (coords) == "chebyshev_half") then
          do k = 0, zlevels
             b_vert(k) = 1.0_dp - sin (dble(k)/dble(zlevels) * MATH_PI/2)
          end do
       else ! default coordinates are uniform (not used if sigma_z = .true.)
          coords = "uniform"
          do k = 0, zlevels
             b_vert(k) = 1.0_dp - dble(k)/dble(zlevels)
          end do
       end if
       a_vert = 1.0_dp - b_vert
    end if

    ! Vertical grid spacing
    a_vert_mass = a_vert(1:zlevels) - a_vert(0:zlevels-1)
    b_vert_mass = b_vert(1:zlevels) - b_vert(0:zlevels-1)
  end subroutine initialize_a_b_vert_case

  subroutine set_thresholds_case
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    implicit none

  end subroutine initialize_dt_viscosity_case

  
  subroutine dump_case (fid)
    implicit none
    integer, intent(in) :: fid

  end subroutine dump_case

  
  subroutine load_case (fid)
    implicit none
    integer, intent(in) :: fid

    if (trim(data_case) /= "climate") then
       read (fid) iwrite
       read (fid) threshold
    end if
  end subroutine load_case

  
end module test_case_mod
