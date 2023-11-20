program make_NCAR_topo
  ! Computes and saves NCAR topography data
  use main_mod
  use ops_mod
  use test_case_mod
  use io_mod
  use topo_grid_descriptor_mod
  implicit none
  integer         :: l
  real(8)         :: fine_mass
  character(9999) :: cmd

  ! Initialize mpi, shared variables and domains
  call init_arch_mod 
  call init_comm_mpi_mod
  
  if (n_process > 1) then
     write (6,'(a)') 'Must run make_NCAR_topo on one core ... aborting'
     call abort
  end if

  ! Read test case parameters
  call read_test_case_parameters

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Standard (shared) parameter values for the simulation
  test_case   = 'make_NCAR_topo'
  zlevels     = 1
  radius      = 6371d0    * KM                   ! mean radius of the Earth
  grav_accel  = 9.80616d0 * METRE/SECOND**2      ! gravitational acceleration 
  tol         = 0d0                              ! non-adaptive grid
  resume      = -1                               ! fresh start
  
  ! Initialize functions
  call assign_functions

  ! Initialize variables
  call initialize ('topo')

  call print_test_case_parameters

  ! Write out grid coordinates for NCAR data
  call write_grid_coords

  ! Generate smoothed NCAR topography
  write (cmd,'(a,a)') 'source ', trim (topo_script)
  call system (cmd)

  ! Generate and save wavelet topography files
  call assign_height (trim (topo_file))
  call forward_topo_transform (topography, wav_topography)
  call inverse_topo_transform (wav_topography, topography)

  ! Check mass conservation
  fine_mass = integrate_hex (topo, z_null, level_end)
  do l = level_end-1, level_start, -1
     write (6,'(a,i2,a,es10.4)') &
          "Relative mass error at level ", l, " = ", abs (integrate_hex (topo, z_null, l) - fine_mass)/fine_mass
  end do

  ! Save compressed wavelet topography data
  call dump_topo

  call finalize
end program Make_NCAR_Topo
