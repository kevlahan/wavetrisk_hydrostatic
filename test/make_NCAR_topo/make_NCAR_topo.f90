program make_NCAR_topo
  ! Computes and saves NCAR topography data
  use main_mod
  use test_case_mod
  use io_mod
  use topo_grid_descriptor_mod
  implicit none
  integer         :: l, nsmth
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
  test_case  = 'make_NCAR_topo'
  zlevels    = 1
  radius     = 6371d0    * KM              ! mean radius of the Earth
  grav_accel = 9.80616d0 * METRE/SECOND**2 ! gravitational acceleration 
  tol        = 0d0                         ! non-adaptive grid
  resume     = -1                          ! fresh start
  
  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize ('topo')

  call print_test_case_parameters

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

  ! Smooth topography
  nsmth = 100
  !call smooth_topo (max_level, nsmth)
  call topo_restriction (max_level, max_level)
  do l = max_level-1, min_level, -1
     call smooth_topo (l, nsmth)
     call topo_restriction (l, l)
  end do
  
  ! Check mass conservation
  fine_mass = integrate_hex (topo, z_null, max_level)
  do l = max_level-1, min_level, -1
     write (6,'(a,i2,a,es10.4)') &
          "Relative mass error at level ", l, " = ", abs (integrate_hex (topo, z_null, l) - fine_mass)/fine_mass
  end do

  ! Save topography data (coarsest scaling function and wavelets) on all non-adaptive levels
  call save_topo

  call finalize
end program
