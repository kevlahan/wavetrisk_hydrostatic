module test_case_mod
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use ops_mod
  use io_mod
  implicit none
  integer        :: angular_order, cp_beg, cp_end, lwin, N, ntaper, k_min, k_max
  real(8)        :: concentration, lat0, lon0, theta0
  character(255) :: data_case, spec_type
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

  function physics_scalar_flux_case (q, dom, id, idE, idNE, idN, v, zlev, type)
    use domain_mod
    implicit none
    real(8), dimension(1:EDGE)                           :: physics_scalar_flux_case
    type(Float_Field), dimension(1:N_VARIABLE,1:zlevels) :: q
    type(domain)                                         :: dom
    integer                                              :: d, id, idE, idNE, idN, v, zlev
    logical, optional                                    :: type
  end function physics_scalar_flux_case

  function physics_velo_source_case (dom, i, j, zlev, offs, dims)
    use domain_mod
    implicit none
    real(8), dimension(1:EDGE)     :: physics_velo_source_case
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims
  end function physics_velo_source_case

  real(8) function surf_geopot_case (d, id)
    implicit none
    integer :: d, id
  end function surf_geopot_case

  subroutine initialize_a_b_vert_case
    implicit none
    allocate (a_vert(0:zlevels),      b_vert(0:zlevels))
    allocate (a_vert_mass(1:zlevels), b_vert_mass(1:zlevels))
  end subroutine initialize_a_b_vert_case

  function z_coords_case (eta_surf, z_s)
    implicit none
    real(8), intent(in)           :: eta_surf, z_s 
    real(8), dimension(0:zlevels) :: z_coords_case
  end function z_coords_case

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
    if (rank == 0) write (6,'(A,A)') "Input file = ", trim (filename)

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
       write (6,'(A)') &
            '********************************************************** Parameters &
            ************************************************************'
       write (6,'(A,A)')      "test_case              = ", trim (test_case)
       write (6,'(A,A)')      "data_case              = ", trim (data_case)
       write (6,'(A,A)')      "spec_type              = ", trim (spec_type)
       if (trim(data_case) == "climate") then
          write (6,'(A,A)')   "physics_type           = ", trim (physics_type)
          write (6,'(A,I0)')  "Nsoil                  = ", Nsoil
       end if
       write (6,'(A,A)')      "run_id                 = ", trim (run_id)
       write (6,'(A,i5)')     "number of domains      = ", N_GLO_DOMAIN
       write (6,'(A,i5)')     "number of processors   = ", n_process
       write (6,'(A,i5)')     "DOMAIN_LEVEL           = ", DOMAIN_LEVEL
       write (6,'(A,i5)')     "PATCH_LEVEL            = ", PATCH_LEVEL
       write (6,'(A,i4)')     "First checkpoint       = ", cp_beg
       write (6,'(A,i4)')     "Last checkpoint        = ", cp_end
       write (6,'(A,i3)')     "min_level              = ", min_level
       write (6,'(A,i3)')     "max_level              = ", max_level
       write (6,'(A,i3)')     "zlevels                = ", zlevels
       write (6,'(A,l1)')     "vert_diffuse           = ", vert_diffuse
       write (6,'(A,i3)')     "level_fill             = ", level_fill
       write (6,'(A,i5)')     "N                      = ", N
       write (6,'(a,l1)')     "local_spec             = ", local_spec
       write (6,'(A,es11.4)') "lat0                   = ", lat0
       write (6,'(A,es11.4)') "lon0                   = ", lon0
       write (6,'(A,es10.4)') "theta0                 = ", theta0
       write (6,'(A,es10.4)') "concentration          = ", concentration
       write (6,'(A,i3)')     "ntaper                 = ", ntaper
       write (6,'(A,i3)')     "angular_order          = ", angular_order
       write (6,*) ' '
    end if
  end subroutine print_test_case_parameters

  subroutine apply_initial_conditions_case
    implicit none
  end subroutine apply_initial_conditions_case

  subroutine update_case
    implicit none
  end subroutine update_case

  subroutine init_sol (dom, i, j, zlev, offs, dims)
    ! Dummy routine
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, k, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims
  end subroutine init_sol

  subroutine init_mean (dom, i, j, zlev, offs, dims)
    ! Initialize mean values
    use utils_mod
    implicit none
    type (Domain)                   :: dom
    integer                         :: i, j, zlev
    integer, dimension (N_BDRY+1)   :: offs
    integer, dimension (2,N_BDRY+1) :: dims

    integer :: d, id, k
   
    d    = dom%id+1
    id   = idx (i, j, offs, dims) 

    do k = 1, zlevels
       sol_mean(S_MASS,k)%data(d)%elts(id+1)        = 0d0
       sol_mean(S_TEMP,k)%data(d)%elts(id+1)        = 0d0
       sol_mean(S_VELO,k)%data(d)%elts(id_edge(id)) = 0d0
    end do
  end subroutine init_mean

  subroutine set_thresholds_case
  end subroutine set_thresholds_case

  subroutine initialize_thresholds_case
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case 
    implicit none

  end subroutine initialize_dt_viscosity_case

  subroutine dump_case (fid)
    implicit none
    integer :: fid

    write (fid) iwrite
    write (fid) threshold
  end subroutine dump_case

  subroutine load_case (fid)
    implicit none
    integer :: fid

    read (fid) iwrite
    read (fid) threshold
  end subroutine load_case
end module test_case_mod
