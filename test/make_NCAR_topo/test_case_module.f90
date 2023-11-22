module test_case_mod
  ! Module file for Held & Suarez (1994) test case
  use comm_mpi_mod
  use utils_mod
  use init_mod
  use std_atm_profile_mod
  use io_mod
  implicit none
  character(9999) :: topo_script
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
    read (fid,*) varname, topo_script
    read (fid,*) varname, topo_file
    read (fid,*) varname, topo_save_wav
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
       write (6,'(a,a)')      "topo_script          = ", trim (topo_script)
       write (6,'(a,a)')      "topo_file            = ", trim (topo_file)
       write (6,'(a,l1)')     "topo_save_wav        = ", topo_save_wav
       write (6,'(a)') &
            '*********************************************************************&
            ************************************************************'
    end if
  end subroutine print_test_case_parameters

  subroutine initialize_thresholds_case
    implicit none
    allocate (threshold(1:N_VARIABLE,zmin:zlevels));     threshold     = 0d0
    allocate (threshold_def(1:N_VARIABLE,zmin:zlevels)); threshold_def = 0d0
    threshold_def = tol
  end subroutine initialize_thresholds_case

  subroutine initialize_dt_viscosity_case
  end subroutine initialize_dt_viscosity_case

  subroutine apply_initial_conditions_case
    implicit none
    integer :: k, l

    do l = level_start, level_end
       call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       call apply_onescale (init_sol,  l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    end do
  end subroutine apply_initial_conditions_case

  subroutine update_case
    ! Update means, bathymetry and penalization mask
    ! not needed in this test case
    use wavelet_mod
    implicit none
    integer :: d, k, l, p

    if (istep /= 0) then
       do d = 1, size(grid)
          do p = n_patch_old(d)+1, grid(d)%patch%length
             call apply_onescale_to_patch (init_mean, grid(d), p-1, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
          end do
       end do
    else ! need to set values over entire grid on restart
       do l = level_start, level_end
          call apply_onescale (init_mean, l, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
       end do
   end if
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
end module test_case_mod
