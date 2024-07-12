program make_NCAR_topo
  ! Computes and saves NCAR topography data
  use main_mod
  use test_case_mod
  use io_mod
  use topo_grid_descriptor_mod
  implicit none
  integer         :: d, l
  real(8)         :: fine_mass
  character(9999) :: cmd, grid_name, jmin_txt, jmax_txt, smth_txt, topo_desc 

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod  
  
  if (n_process > 1) then
     if (rank == 0) write (6,'(/,a,/)') '!!! Must run make_NCAR_topo on a single core ... aborting !!!'
     call abort
  end if

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  test_case     = 'make_NCAR_topo'
  zlevels       = 1
  radius        = 6371d0    * KM              ! mean radius of the Earth
  grav_accel    = 9.80616d0 * METRE/SECOND**2 ! gravitational acceleration 
  tol           = 0d0                         ! non-adaptive grid
  resume        = -1                          ! fresh start
  analytic_topo = .false.                     ! use analytic topography

  Area_min  = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**max_level)
  Area_max  = 4d0*MATH_PI * radius**2 / (10d0 * 4d0**min_level)
  dx_min    = sqrt (2d0 / sqrt(3d0) * Area_min)              
  dx_max    = sqrt (2d0 / sqrt(3d0) * Area_max)

  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize ('topo')
  
  call print_test_case_parameters
 
  if (analytic_topo) then ! analytic topography (restricted from finest level)
     topo_file = "analytic_topo"
  else                    ! NCAR topography
     ! Write out grid coordinates for NCAR data
     write (jmin_txt,'(i2.2)') min_level
     write (jmax_txt,'(i2.2)') max_level
     grid_name = 'J'//trim(jmin_txt)//'J'//trim(jmax_txt)
     topo_desc = trim(grid_name)//'_topo_grid'

     call write_grid_coords (topo_desc)

     topo_desc = trim(grid_name)//'_topo_grid'

     write (smth_txt,'(f5.1)') smth_scl
     smth_txt = repeat ('0', 5-len_trim(adjustl(smth_txt))) // adjustl(smth_txt) ! add leading zeros

     ! Generate smoothed NCAR topography using cube_to_target
     write (cmd,'(a)') &
          "./cube_to_target &
          
          --no_ridges --smoothing_over_ocean &
          
          --grid_descriptor_file="//"'"//trim(topo_desc)//"' &
          
          --intermediate_cs_name="//"'"//trim(topo_data)//"' &
          
          --output_grid="//"'"//trim(grid_name)//"' &
          
          --smoothing_scale="//trim(smth_txt)// &
          
          " -u Nicholas Kevlahan kevlahan@mcmaster.ca -q"//" ." 

     call system (trim(cmd))

     ! Assign to grid topography at max_level
     topo_file = trim(grid_name)//"_"//trim(smth_txt)//"km"
     call assign_height (trim (topo_file))
  end if

  ! Restrict topography
  call topo_restriction (min_level, max_level)
  topography%bdry_uptodate = .false.
  call update_bdry (topography, NONE)

  ! Compute topography gradients at edges (stored in sol(S_VELO,1))
  do l = max_level, min_level, -1
     do d = 1, size(grid)
        scalar => topography%data(d)%elts
        velo   => sol(S_VELO,1)%data(d)%elts
        call apply_onescale_d (cal_grad_topo, grid(d), l, z_null, 0, 1)
        nullify (scalar, velo)
     end do
  end do
  sol(S_VELO,1)%bdry_uptodate = .false.
  call update_bdry (sol(S_VELO,1), NONE)

  ! Interpolate topography gradients to nodes
  do l = max_level, min_level, -1
     do d = 1, size(grid)
        velo  => sol(S_VELO,1)%data(d)%elts
        velo1 => grid(d)%u_zonal%elts
        velo2 => grid(d)%v_merid%elts
        call apply_onescale_d (interp_UVW_latlon, grid(d), l, z_null, 0, 1)
        nullify (velo1, velo2)
     end do
  end do
  
  ! Check mass conservation
  fine_mass = integrate_hex (topo, z_null, max_level)
  do l = max_level-1, min_level, -1
     write (6,'(a,i2,a,es10.4)') &
          "Relative mass error at level ", l, " = ", abs (integrate_hex (topo, z_null, l) - fine_mass)/fine_mass
  end do

  ! Save topography data on all non-adaptive levels
  call save_topo

  call finalize
end program
