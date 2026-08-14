module main_mod

  use, intrinsic :: iso_fortran_env, only : int8

  use mpi_f08, only : MPI_Allgather, MPI_Allgatherv, MPI_Allreduce, MPI_Exscan, &
       MPI_INTEGER, MPI_SUCCESS, MPI_SUM
  
  use, intrinsic :: ieee_arithmetic

  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, AT_NODE, AT_EDGE, BDRY_THICKNESS, CP_EVERY, DATA_GRID, DAY, EDGE, MATH_PI, MULT, N_BDRY, &
       N_CHDRN, N_GLO_DOMAIN, N_VARIABLE, NCAR_topo, NONE, &
       EAST, NORTH, NORTHEAST, NORTHWEST, SOUTH, SOUTHEAST, SOUTHWEST, WEST, &
       S_DIVU, S_ROTU, S_MASS, S_TEMP, S_VELO, TRIAG, XU_GRID, RT, DG, UP, LORT, UPLT, Laplace_divu, Laplace_rotu, Laplace_sclr, &
       adapt_dt, Area_avg, C_visc, a_vert, a_vert_mass, b_vert, b_vert_mass, cfl_safety, compressible, coord, cp_idx, &
       dt, dt_init, dt_write, dx_avg, gamma, grav_accel, iremap, iremap_max, istep, istep_cumul, itime, iwrite, &
       linear_solver, level_start, level_end, log_min_mass, log_total_mass, match_time, &
       min_mass, min_mass_remap, min_level, max_level, mode_split, n_active, n_domain, n_node_old, n_patch_old, optimize_grid, &
       P_top, penta_node, penta_node_std, physics_type, physics_model, r_adv, r_dif, radius, remap, &
       R_d, run_id, sigma_z, resume, scalars, sso, test_case, theta_grid, threshold, threshold_def, &
       time, time_end, timeint_type, vert_diffuse, vtk_grid, wave_speed, z_null, zlevels, zmin, zmax

  use adapt_mod,          only : adapt, WT_after_step
  use coarse_grid_mod,    only : read_optim_grid, smooth_Xu, update_geom_check_grid, zrotate 
  use comm_mod,           only : get_coord, set_coord
  use diagnostics_mod,    only : rho_dz_i, theta_i, theta2temp
  use domain_ops_mod,     only : apply_interscale, apply_no_bdry,  apply_onescale2
  use geom_mod,           only : number_hex
  use integrate_mod,      only : integrate_hex, integrate_tri
  use checkpoint_mod,     only : dump_adapt_mpi, load_adapt_mpi, read_checkpoint_directory
  use io_vtk_mod,         only : write_and_export
  use lin_solve_mod,      only : Full_Multigrid, Scheduled_Relaxation_Jacobi
  use lnorms_mod,         only : lnorm
  use mask_mod,           only : init_masks, mask_adj_child
  use multi_level_mod,    only : trend_ml
  use NCAR_topo_mod,      only : load_topo
  use ops_mod,            only : comp_offs3
  use patch_mod,          only : Patch, PATCH_SIZE
  use parallel_block_mod, only : Block_Bdry_Storage, Block_Data, &
       Block_Bdry_Link, Block_Ghost_Storage, Block_Stencil_Address, &
       STORE_NONE, STORE_PATCH, STORE_BDRY, STORE_GHOST, &
       NGB_INTERNAL, NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT, NGB_OTHER, &
       block_source, block_received, block_local, &
       block_source_catalog_index, block_retained_source_index, &
       block_migrating_source_index, &
       block_received_catalog_index, &
       block_local_catalog_index, &
       packed_block_nbyte, pack_block, unpack_block, &
       check_block_storage, clear_block_staging, &
       install_local_blocks
  use refine_patch_mod,   only : add_second_level
  use remap_mod,          only : remap_vertical_coordinates
  use utils_mod,          only : hex_len, hex_pedlen, interp, nu_scale, porous_density, tri_perim
  use vert_diffusion_mod, only : vertical_diffusion
  use wavelet_mod,        only : forward_wavelet_transform, init_wavelets, inverse_scalar_transform, inverse_wavelet_transform

  use arch_mod, only : abort_run, barrier, block_catalog, comm, distribute_grid, glo_id, init_arch_mod, loc_id, &
       n_glo_block, n_process, owner, Parallel_Block, rank

  use comm_mpi_mod, only : comm_nodes3_mpi, init_comm_mpi, recv_lengths, recv_offsets, req, send_lengths, send_offsets, &
       sum_int,  sync_max_int, sync_min_real, write_load_conn

  use domain_mod, only : Domain, bernoulli, divu, dscalar, grid, &
       dvelo, exner, exner_fun, h_flux, horiz_flux, ke, qe, scalar, mass, temp, velo, vort, wc_s, wc_u, &
       Kt, Kv, Laplacian_scalar, Laplacian_vector, penal_edge, penal_node, sso_param, &
       sol, sol_mean, tke, topography, topography_data, trend, wav_coeff, wav_tke, id_edge, idx, &
       subtree_weight_Domain, count_subtree_patches_Domain, extract_subtree_patches_Domain, subtree_depth_Domain, &
       compact_subtree_storage_Domain, copy_subtree_nodes_Domain, copy_subtree_field_Domain, renumber_subtree_neigh_Domain, &
       get_bdry_dims_domain

  use init_mod, only : apply_initial_conditions, elliptic_solver, init_geometry, init_grid, &
       initialize_a_b_vert, initialize_thresholds, initialize_dt_viscosity, &
       precompute_geometry, set_level, set_thresholds, z_coords

  use time_integr_mod, only : dt_step, dt_step_split, init_RK_mem, q1, q2, q3, q4, dq1, &
       Euler, Euler_split, RK2_split, RK3, RK3_split, RK33_opt, RK34_opt, RK4, RK4_split, RK45_opt

  use coord_arithmetic_mod

  use parallel_block_mpi_mod, only : Block_Migration_Manifest, &
       build_block_migration_manifest, &
       check_block_migration_manifest, &
       exchange_block_migration_sizes, &
       exchange_block_migration_payloads, &
       clear_block_migration_manifest

#ifdef PHYSICS
  use init_physics_mod,  only : init_physics, init_soil_grid
  use physics_trend_mod, only : physics_simple_step, trend_physics_Held_Suarez  
  use callkeys,          only : lverbose
#endif

  implicit none
  
  private
  public :: cpt_min_mass, initialize, restart, time_step
  
  
  type Initial_State
     integer              :: n_patch, n_bdry_patch, n_node, n_edge, n_tria
     integer, allocatable :: pack_len(:,:), unpk_len(:,:)
  end type Initial_State

integer, allocatable :: n_active_edges(:), n_active_nodes(:), node_level_start(:), edge_level_start(:)
real(dp)             :: dt_new, initial_total_mass, time_mult
real(dp)             :: dt_loc, min_mass_loc

type(Initial_State), allocatable :: ini_st(:)


type(Block_Migration_Manifest) :: block_migration

contains


  subroutine initialize (run_id) 
    ! Initialize from checkpoint or adapt to initialize conditions
    ! (solution is saved and restarted to balance load)
    implicit none
    character(*), intent(in) :: run_id

    if (max_level < min_level) then
       if (rank == 0) then
          write (6,'(//,a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a)')   "!!                                                                                 !!"
          write (6,'(2(a,i2),a)') "!!                max_level < min_level: ", max_level, " < ", min_level, &
               " ... aborting                      !!"
          write (6,'(a)')   "!!                                                                                 !!"
          write (6,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a,//)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       end if
       call abort_run
    end if

#ifdef PHYSICS
    if (physics_model .and. physics_type == "Simple") call init_soil_grid
#endif

    ! Linear solver
    select case (linear_solver)
    case ("FMG")
       elliptic_solver => Full_Multigrid
    case ("SRJ")
       elliptic_solver => Scheduled_Relaxation_Jacobi
    end select

    call set_time_integrator

    ! Check validity of various parameter choices
    if (compressible .and. mode_split) then
       if (rank == 0) write (6,'(a)') "Cannot use mode splitting with compressible dynamics ... aborting"
       call abort_run
    end if

    ! Ensure data and check point is saved at least once
    if (dt_write > time_end) then
       dt_write = time_end
       CP_EVERY = 1
    end if

    if (resume >= 0) then
       cp_idx = resume
       call read_checkpoint_directory (cp_idx)
       call restart
       call test_parallel_block_split
       call test_subtree_extraction
       call build_block_migration_manifest(block_migration)
       call check_block_migration_manifest(block_migration)
       call test_block_migration_sizes(block_migration)
       call test_block_migration_payloads(block_migration)
       call test_install_local_blocks
       call finalize_test_block_migration(block_migration)
    else
       call init_basic
       call init_structures

       if (rank == 0) write (6,'(/,a,/)') &
            '----------------------------------------------------- Adapting initial grid &
            &------------------------------------------------------'
       if (NCAR_topo) call load_topo

       ! Apply initial conditions 
       call count_active

       do while (level_end < max_level)
          if (rank == 0) write (6,'(a,i2,a,i2)') 'Initial refinement Level ', level_end, ' -> ', level_end+1

          ! Set current number of nodes
          node_level_start = grid%node%length+1; edge_level_start = grid%midpt%length+1
          n_patch_old      = grid%patch%length;  n_node_old       = grid%node%length

          call adapt (set_thresholds) ! extend grid
          call count_active           ! apply initial conditions on new grid and count new active wavelets

          if (n_active(AT_NODE) == 0 .and. n_active(AT_EDGE) == 0) exit ! no further active grid points
       end do

       if (rank == 0) write (6,'(a,/)') &
            '------------------------------------------------- Finished adapting initial grid &
            &-------------------------------------------------'
       if (rank==0) write (6,'(a,i12,/)') 'Initial number of active wavelets = ', sum (n_active)

       call adapt (set_thresholds)
       dt_new = cpt_dt ()
       call count_active

       if (trim (test_case) /= "make_NCAR_topo") call write_checkpoint
    end if
    if (trim (test_case) /= "make_NCAR_topo" .and. trim (test_case) /= "spherical_harmonics") call write_and_export (vtk_grid)
    call barrier

    call initialize_dt_viscosity
    call initialize_thresholds

#ifdef PHYSICS
    if (physics_model .and. physics_type == "Simple") call init_physics 
#endif
  end subroutine initialize

  subroutine count_active
    ! Apply initial conditions and count  number of active node and edge wavelets
    implicit none
    integer :: d, k, v

    call apply_initial_conditions
    call forward_wavelet_transform (sol, wav_coeff) 

    ! Initialize thresholds to default values for this grid (possibly based on mean values)
    call initialize_thresholds

    if (level_end /= level_start) then
       n_active = 0
       do k = zmin, zmax
          do d = 1, size(grid)
             do v = scalars(1), scalars(2)
                wc_s => wav_coeff(v,k)%data(d)%elts
                n_active(AT_NODE) = n_active(AT_NODE) &
                     + count (abs(wc_s(node_level_start(d):grid(d)%node%length)) >= threshold(v,k))
                nullify (wc_s)
             end do
             wc_u => wav_coeff(S_VELO,k)%data(d)%elts
             n_active(AT_EDGE) = n_active(AT_EDGE) &
                  + count (abs(wc_u(edge_level_start(d):grid(d)%midpt%length)) >= threshold(S_VELO,k))
             nullify (wc_u)
          end do
       end do

       ! Sum results over all ranks
       n_active(AT_NODE) = sum_int (n_active(AT_NODE)) ; n_active(AT_EDGE) = sum_int(n_active(AT_EDGE))

       if (rank == 0) write (6,'(a,i2,1x,2(a,i12,1x),/)') &
            'Level = ', level_end, 'number of active node wavelets = ', n_active(AT_NODE), &
            'number of active edge wavelets = ', n_active(AT_EDGE)
    end if
  end subroutine count_active

  subroutine record_init_state (init_state)
    implicit none
    type(Initial_State), allocatable, intent(out) :: init_state(:)

    integer :: d, i, v


    allocate (init_state(size(grid)))

    do d = 1, size(grid)
       allocate (init_state(d)%pack_len(AT_NODE:AT_EDGE,n_glo_block))
       allocate (init_state(d)%unpk_len(AT_NODE:AT_EDGE,n_glo_block))

       init_state(d)%n_patch      = grid(d)%patch%length
       init_state(d)%n_bdry_patch = grid(d)%bdry_patch%length
       init_state(d)%n_node       = grid(d)%node%length
       init_state(d)%n_edge       = grid(d)%midpt%length
       init_state(d)%n_tria       = grid(d)%ccentre%length

       do i = 1, n_glo_block
          do v = AT_NODE, AT_EDGE
             init_state(d)%pack_len(v,i) = grid(d)%pack(v,i)%length
             init_state(d)%unpk_len(v,i) = grid(d)%unpk(v,i)%length
          end do
       end do
    end do

  end subroutine record_init_state

  subroutine set_time_integrator
    ! Sets time integration scheme and associated stability factor as specified in test case
    implicit none

    if (mode_split) then 
       select case (timeint_type)
       case ("Euler")
          dt_step_split => Euler_split
          r_adv = 1.0_dp
          r_dif = 2.0_dp
       case ("RK2")
          dt_step_split => RK2_split
          r_adv = 1.0_dp
          r_dif = 2.0_dp
       case ("RK3")
          dt_step_split => RK3_split
          r_adv = sqrt (3.0_dp)
          r_dif = 2.51_dp
       case ("RK4")
          dt_step_split => RK4_split
          r_adv = 2 * sqrt (2.0_dp)
          r_dif = 2.77_dp
       case default
          dt_step_split => RK3_split
          r_adv = sqrt (3.0_dp)
          r_dif = 2.51_dp
       end select
    else
       select case (timeint_type)
       case ("Euler")
          dt_step => Euler
          r_adv = 1.0_dp
          r_dif = 2.0_dp
       case ("RK3")
          dt_step => RK3
          r_adv = sqrt (3.0_dp)
          r_dif = 2.5_dp
       case ("RK4")
          dt_step => RK4
          r_adv = 2 * sqrt (2.0_dp)
          r_dif = 2.77_dp
       case default
          dt_step => RK3
          r_adv = sqrt (3.0_dp)
          r_dif = 2.5_dp
       case ("RK33")
          dt_step => RK33_opt
          r_adv = sqrt (3.0_dp)
       case ("RK34")
          dt_step => RK34_opt
          r_adv = 2.0_dp
       case ("RK45")
          dt_step => RK45_opt
          r_adv = 3.28_dp
       end select
    end if
  end subroutine set_time_integrator

  subroutine write_checkpoint 
    implicit none

    cp_idx = cp_idx + 1
    if (rank == 0) then
       write (6,'(/,a,/)') &
            '************************************************************************&
            &**********************************************************'
       write (6,'(a,i4,a,es10.4,/)') 'Saving checkpoint ', cp_idx, ' at time [day] = ', time / DAY
    end if

    call dump_adapt_mpi (cp_idx)
    call restart

  end subroutine write_checkpoint

  subroutine restart 
    ! Fresh restart from checkpoint data (all structures reset)

    implicit none

    call deallocate_structures  ! deallocate all dynamic arrays and variables
    call init_basic

    if (rank == 0) then
       write (6,'(a,/)') &
            '********************************************************* Begin Restart &
            &**********************************************************'
       write (6,'(a,i4,/)') 'Restarting from checkpoint ', cp_idx
    end if
    call barrier

    call init_structures          ! initialize coarsest grid and distribute load
    call load_adapt_mpi (cp_idx)  ! load checkpoint data
    if (NCAR_topo) call load_topo ! load topography data

    ! Compute masks based on active wavelets in saved data
    call adapt (set_thresholds, .false.)  
    call inverse_wavelet_transform (wav_coeff, sol, jmin_in=level_start-1)
    if (vert_diffuse) call inverse_scalar_transform (wav_tke, tke, jmin_in=level_start-1)

    if (trim(test_case) /= "spherical_harmonics") then
       call initialize_dt_viscosity
       call initialize_thresholds

       if (log_total_mass) call cal_total_mass (.true.)

       itime = nint (time * time_mult, 8)

       dt_new = cpt_dt ()
       dt = dt_new

       if (rank == 0) then
          write (6,'(/,A,es12.6,3(A,es8.2),A,I2,A,I9,/)') &
               'time [d] = ', time/DAY, &
               '  mass threshold = ', sum (threshold(S_MASS,1:zlevels)) / dble (zlevels), &
               ' temp threshold = ', sum (threshold(S_TEMP,1:zlevels)) / dble (zlevels), &
               ' velo threshold = ', sum (threshold(S_VELO,1:zlevels)) / dble (zlevels), &
               ' Jmax = ', level_end, &
               '  dof = ', sum (n_active)
          write (6,'(A)') &
               '********************************************************** End Restart &
               &***********************************************************'
       end if
    end if
  end subroutine restart

  subroutine time_step 
    implicit none
    integer(8) :: idt, ialign
    logical    :: save_data 

    ! New time step
    istep       = istep       + 1
    istep_cumul = istep_cumul + 1
    dt          = dt_new
    idt         = nint (dt * time_mult, 8)

    ! Check whether to save data
    save_data = .false.
    if (dt_write <= time_end) then
       ialign = nint (dt_write * time_mult, 8)
       if (ialign > 0 .and. istep > 1) then
          save_data = (modulo (itime+idt, ialign) < modulo (itime, ialign))
       else
          save_data = .false.
       end if

       if (save_data .and. match_time) then ! adjust time step to match save time exactly
          idt = ialign - modulo (itime,ialign)
          dt = idt / time_mult
       end if
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Dynamics time step
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (mode_split) then
       call dt_step_split (dt)
    else
       call dt_step (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
    end if

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Physics split step
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Ocean (incompressible) 
    if (vert_diffuse) call vertical_diffusion ! ocean (incompressible) models

    ! Atmosphere (compressible) 
#ifdef PHYSICS
    if (physics_model) then
       select case (physics_type)
       case ("Held_Suarez")
          call Euler (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_physics_Held_Suarez, dt) 
       case ("Simple")
          call physics_simple_step (dt)
       end select
    end if
#endif 

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Grid adaptation
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (zmin < 1) call WT_after_step (sol(:,zmin:0), wav_coeff(:,zmin:0), level_start-1) ! compute wavelet coefficients in soil levels
    call adapt (set_thresholds)
    call inverse_wavelet_transform (wav_coeff, sol)

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Vertical remapping (after grid adaptation to ensure ADJCENT_ZONE and ZERO cells are remapped consistently)
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    min_mass = cpt_min_mass () 
    if (remap .and. min_mass < min_mass_remap .or. iremap == iremap_max) then
       if (rank == 0 .and. log_min_mass) write (6,'(a)') 'Remapping vertical coordinates ...'
       call remap_vertical_coordinates
       iremap = 1
    else
       iremap = iremap + 1
    end if
    if (log_total_mass) call cal_total_mass (.false.) ! change in total mass

    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !    Update time step and save data
    ! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    itime  = itime + idt
    time   = itime / time_mult
    dt_new = cpt_dt ()

    if (save_data) then
       iwrite = iwrite + 1
       if (modulo (iwrite, CP_EVERY) == 0) call write_checkpoint
       call write_and_export (vtk_grid)
    end if
  end subroutine time_step

  subroutine init_basic

    implicit none

    integer :: l

    call initialize_a_b_vert

    allocate (n_active_edges(min_level-1:max_level), n_active_nodes(min_level-1:max_level))
    n_active_edges = 0; n_active_nodes = 0

    allocate (dx_avg(min_level-1:max_level), Area_avg(min_level-1:max_level))
    do l = min_level-1, max_level
       Area_avg(l) = 4*MATH_PI * radius**2 / number_hex (l)
       dx_avg(l)   = sqrt (2 / sqrt(3.0_dp) * Area_avg(l))
    end do

    istep     = 0
    time_mult = 1.0_dp
  end subroutine init_basic

  subroutine init_structures
    ! Initialize dynamical arrays and structures
    implicit none
    integer :: ip

    level_start = min_level
    level_end   = level_start

    ! Distribute and balance grid over processors (necessary for correct restart!)
    call distribute_grid (cp_idx)

    call init_grid
    call init_comm_mpi
    call init_geometry
    call update_geom_check_grid

    ! Computational grid at min_level
    select case (optimize_grid)
    case (XU_GRID)
       call smooth_Xu
    case (DATA_GRID)
       call read_optim_grid
    end select

    ! Coordinates of pentagon nodes
    do ip = 1, 12
       call zrotate (penta_node_std(ip), penta_node(ip), theta_grid)
       penta_node(ip) = radius * penta_node(ip)
    end do

    call comm_nodes3_mpi (get_coord, set_coord)
    call precompute_geometry

    allocate (node_level_start(size(grid)), edge_level_start(size(grid)))
    edge_level_start = 0
    node_level_start = 0

    if (rank == 0) write (6,'(a,i2,a)') 'Make level J_min = ', min_level, ' ...'

    call init_wavelets
    call init_masks
    call add_second_level

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, BDRY_THICKNESS)
    call apply_interscale (mask_adj_child, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)

    if (time_end > 0.0_dp) then
       time_mult = 0.5_dp * real( huge(itime), kind=dp ) / time_end
    else
       time_mult = 0.0_dp
    end if

    allocate (n_patch_old(size(grid)), n_node_old(size(grid))); n_patch_old = 2

    call init_RK_mem

    allocate (threshold(1:N_VARIABLE,zmin:zmax), threshold_def(1:N_VARIABLE,zmin:zmax))
  end subroutine init_structures

  subroutine cal_total_mass (initialize_total_mass)
    ! Compute total mass over all vertical layers
    implicit none
    logical, intent(in) :: initialize_total_mass

    integer      :: k
    real(dp)     :: total_mass, mass_error
    character(3) :: int_type = "hex"

    total_mass = 0.0_dp
    do k = 1, zlevels
       select case (int_type)
       case ("hex") ! coarsest level only
          total_mass = total_mass + integrate_hex (rho_dz_i, k, level_start)
       case ("tri") ! all levels
          total_mass = total_mass + integrate_tri (rho_dz_i, k)
       end select
    end do
    if (initialize_total_mass) initial_total_mass = total_mass

    mass_error =  (total_mass - initial_total_mass) / initial_total_mass

    if (rank == 0 .and. .not. initialize_total_mass) then
       select case (int_type)
       case ("hex") 
          write (6,'(a,es11.4)') "Relative total mass error on coarsest level (hexagons) = ", mass_error
       case ("tri")
          write (6,'(a,es11.4)') "Relative total mass error on adaptive grid (triangles) = ", mass_error
       end select
    end if
  end subroutine cal_total_mass

  real(dp) function cpt_dt ()
    ! Calculates time step, minimum relative mass and active nodes and edges
    implicit none
    integer, dimension(2) :: n_active_loc

    if (adapt_dt) dt_loc = 1e16_dp
    n_active_nodes = 0
    n_active_edges = 0

    ! Calculate minimum time step, number of active nodes and edges
    call apply_no_bdry (cal_min_dt, z_null)

    ! Time step
    if (adapt_dt) then
       cpt_dt = sync_min_real (dt_loc)
    else
       cpt_dt = dt_init
    end if

    ! Active nodes and edges
    n_active_loc = (/ sum (n_active_nodes(level_start:level_end)), sum(n_active_edges(level_start:level_end)) /)
    n_active = sum_int (n_active_loc)
    level_end = sync_max_int (level_end)
  end function cpt_dt

  subroutine cal_min_dt (dom, i, j, zlev, offs, dims)
    ! Calculates time step and number of active nodes and edges
    ! (uses exact local CFL stability formula for hexagons/pentagons)

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                        :: d, e, id, id_i, k, l, lev
    integer,  dimension(1:EDGE)    :: ide
    real(dp)                       :: acoustic_speed, dt_adv, P_k, rho_dz, T, theta
    real(dp), dimension(0:zlevels) :: P
    real(dp), dimension(1:EDGE)    :: dx, speed

    d    = dom%id + 1
    id   = idx (i, j, offs, dims)
    id_i = id + 1
    ide  = id_edge (id)

    lev  = dom%level%elts(id_i)

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       n_active_nodes(lev) = n_active_nodes(lev) + 1

       dx = dom%len%elts(ide)

       if (adapt_dt) then
          if (compressible) then
             P(zlevels) = p_top
             do l = zlevels-1, 0, -1
                k = l + 1 
                rho_dz = sol_mean(S_MASS,k)%data(d)%elts(id_i) + sol(S_MASS,k)%data(d)%elts(id_i) 
                P(l)   = P(l+1) + grav_accel * rho_dz
                P_k    = interp (P(l), P(l+1))

                theta  = theta_i (dom, i, j, k, offs, dims)
                T      = theta2temp (theta, P_k)

                acoustic_speed = sqrt (gamma * R_d * T)
                speed          = abs (sol(S_VELO,k)%data(d)%elts(ide))

                dt_adv = cfl_safety * minval (r_adv * dx/4 / (speed + acoustic_speed))
                dt_loc = min (dt_init, dt_loc, dt_adv)
             end do
          else
             do k = 1, zlevels
                speed = abs (sol(S_VELO,k)%data(d)%elts(ide))

                dt_adv = cfl_safety * minval (r_adv * dx/4 / (speed + wave_speed))
                dt_loc = min (dt_init, dt_loc, dt_adv)
             end do
          end if
       end if
    end if

    do e = 1, EDGE
       if (dom%mask_e%elts(ide(e)) >= ADJZONE) n_active_edges(lev) = n_active_edges(lev) + 1
    end do
  contains
    ! Routines to compute exact amplification factors for diffusive stability on adaptive grid
    ! Example: dt_dif = r_dif / theta_max_sclr ()**Laplace_sclr / nu_scale (S_MASS,1)

    function dt_max_diffusive () 
      ! Maximum diffusive time step
      implicit none
      real(dp), dimension(3) :: dt_max_diffusive

      dt_max_diffusive = [ &
           r_dif / theta_max_sclr ()**Laplace_sclr / nu_scale (S_MASS,1), &
           r_dif / theta_max_divu ()**Laplace_divu / nu_scale (S_DIVU,1), &
           r_dif / theta_max_rotu ()**Laplace_rotu / nu_scale (S_ROTU,1) ]
    end function dt_max_diffusive

    real(dp) function theta_max_sclr ()
      ! Maximum amplification factor for scalar diffusion
      ! (conservative Gershgorin bounds to include irregular pentagons:
      ! about 33% too conservative compared to sharp bounds for regular hexagons)
      implicit none
      real(dp) :: sigma

      sigma = sum (hex_pedlen (dom, i, j, offs, dims) / hex_len (dom, i, j, offs, dims)) * dom%areas%elts(id+1)%hex_inv 

      theta_max_sclr = 2.0_dp * sigma
    end function theta_max_sclr

    real(dp) function theta_max_divu ()
      ! Maximum amplification factor for divergence diffusion
      ! (conservative Gershgorin bounds to include irregular pentagons:
      ! about 33% too conservative compared to sharp bounds for regular hexagons)
      implicit none
      integer                     :: idE, idNE, idN
      real(dp)                    :: chi, chiE, chiNE, chiN
      real(dp), dimension(1:EDGE) :: theta

      idE  = idx (i+1, j,   offs, dims)
      idNE = idx (i+1, j+1, offs, dims)
      idN  = idx (i,   j+1, offs, dims)

      chi   = sum (hex_pedlen (dom, i,   j,   offs, dims)) * dom%areas%elts(id  +1)%hex_inv
      chiE  = sum (hex_pedlen (dom, i+1, j,   offs, dims)) * dom%areas%elts(idE +1)%hex_inv
      chiNE = sum (hex_pedlen (dom, i+1, j+1, offs, dims)) * dom%areas%elts(idNE+1)%hex_inv
      chiN  = sum (hex_pedlen (dom, i,   j+1, offs, dims)) * dom%areas%elts(idN +1)%hex_inv

      theta(RT+1) = (chi + chiE ) / dom%len%elts(EDGE*id+RT+1)
      theta(DG+1) = (chi + chiNE) / dom%len%elts(EDGE*id+DG+1)
      theta(UP+1) = (chi + chiN ) / dom%len%elts(EDGE*id+UP+1)

      theta_max_divu = maxval (theta)
    end function theta_max_divu

    real(dp) function theta_max_rotu ()
      ! Maximum amplification factor for curl-curl diffusion
      ! (conservative Gershgorin bounds to include irregular pentagons:
      ! about 15% too conservative compared to sharp bounds for regular hexagons)
      implicit none
      integer                     :: id, idW, idS
      real(dp)                    :: chi_LORT, chi_LORT_W, chi_UPLT, chi_UPLT_S
      real(dp), dimension(1:EDGE) :: theta

      id  = idx (i,   j,   offs, dims)
      idW = idx (i-1, j,   offs, dims)
      idS = idx (i,   j-1, offs, dims)

      chi_UPLT_S = tri_perim (dom, i,   j-1, UPLT, offs, dims) / dom%triarea%elts(TRIAG*idS+UPLT+1)
      chi_LORT   = tri_perim (dom, i,   j,   LORT, offs, dims) / dom%triarea%elts(TRIAG*id +LORT+1)
      chi_UPLT   = tri_perim (dom, i,   j,   UPLT, offs, dims) / dom%triarea%elts(TRIAG*id +UPLT+1)
      chi_LORT_W = tri_perim (dom, i-1, j,   LORT, offs, dims) / dom%triarea%elts(TRIAG*idW+LORT+1)

      theta(RT+1) = (chi_UPLT_S + chi_LORT  ) / dom%pedlen%elts(EDGE*id+RT+1)
      theta(DG+1) = (chi_LORT   + chi_UPLT  ) / dom%pedlen%elts(EDGE*id+DG+1)
      theta(UP+1) = (chi_UPLT   + chi_LORT_W) / dom%pedlen%elts(EDGE*id+UP+1)

      theta_max_rotu = maxval (theta)
    end function theta_max_rotu
  end subroutine cal_min_dt

  real(dp) function cpt_min_mass ()
    ! Calculates minimum relative mass
    implicit none

    min_mass_loc = 1e16_dp
    call apply_no_bdry (cal_min_mass, z_null)

    cpt_min_mass = sync_min_real (min_mass_loc)

    if (log_min_mass .and. rank == 0) write (6,'(a,es11.4)') "Minimum relative mass = ", cpt_min_mass
  end function cpt_min_mass

  subroutine cal_min_mass (dom, i, j, zlev, offs, dims)
    ! Minimum mass compared to initial mass in a vertical column

    implicit none

    type(Domain), intent(inout) :: dom
    integer,      intent(in)    :: i, j, zlev
    integer,      intent(in)    :: offs(N_BDRY+1) 
    integer,      intent(in)    :: dims(2,N_BDRY+1)

    integer                        :: d, id, k
    real(dp)                       :: mass_ratio, P_s, z_s
    real(dp), dimension(0:zlevels) :: z
    real(dp), dimension(1:zlevels) :: dz, init_rho_dz, rho_dz

    id = idx (i, j, offs, dims) + 1
    d  = dom%id + 1

    do k = 1, zlevels
       rho_dz(k) = sol(S_MASS,k)%data(d)%elts(id) + sol_mean(S_MASS,k)%data(d)%elts(id)
    end do

    if (compressible) then
       P_s = grav_accel * sum (rho_dz) + P_top ! surface pressure
       init_rho_dz = a_vert_mass + b_vert_mass * P_s / grav_accel
    else
       z_s = topography%data(d)%elts(id)
       if (sigma_z) then
          z = z_coords (0.0_dp, z_s)
       else
          z = b_vert * z_s ! assumes zero free surface perturbation initial condition
       end if
       dz = z(1:zlevels) - z(0:zlevels-1)

       do k = 1, zlevels
          init_rho_dz(k) = porous_density (d, id, k) * dz(k)
       end do
    end if

    mass_ratio = minval (rho_dz / abs (init_rho_dz))

    if (mass_ratio <= 0.0_dp .or. ieee_is_nan (mass_ratio)) then
       write (6,'(a)') "A layer has collapsed ... aborting"
       call abort_run
    else
       min_mass_loc = min (min_mass_loc, mass_ratio)
    end if
  end subroutine cal_min_mass


subroutine test_subtree_extraction
  ! Build and verify every candidate block whose source Domain is
  ! currently local to this rank. Detailed diagnostics are printed for
  ! one representative block; all blocks contribute to rank totals.

  implicit none

  integer :: b
  integer :: b_verbose
  integer :: d
  integer :: i_block

  integer :: n_block_built
  integer :: n_block_owned
  integer :: n_block_migrating

  integer :: n_patch_block
  integer :: n_bdry_block
  integer :: n_ghost_block
  integer :: n_stencil_block
  integer :: n_remote_block
  integer :: n_value_block

  integer :: n_patch_total
  integer :: n_bdry_total
  integer :: n_ghost_total
  integer :: n_stencil_total
  integer :: n_remote_total
  integer :: n_value_total

  integer :: n_pack_block
  integer :: n_pack_byte_total
  integer :: n_pack_byte_max

  !
  ! Count the source-local candidate blocks before allocating the
  ! persistent block store and its retained/migrating index sets.
  !
  n_block_built     = 0
  n_block_owned     = 0
  n_block_migrating = 0

  do b = 1, size(block_catalog)

     if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

     n_block_built = n_block_built + 1

     if (block_catalog(b)%owner == rank) then
        n_block_owned = n_block_owned + 1
     else
        n_block_migrating = n_block_migrating + 1
     end if

  end do

  if (n_block_built < 1) then
     error stop "test_subtree_extraction: no source-local blocks"
  end if

  if (allocated(block_source)) then
     deallocate(block_source)
  end if

  if (allocated(block_source_catalog_index)) then
     deallocate(block_source_catalog_index)
  end if

  if (allocated(block_retained_source_index)) then
     deallocate(block_retained_source_index)
  end if

  if (allocated(block_migrating_source_index)) then
     deallocate(block_migrating_source_index)
  end if

  allocate(block_source(n_block_built))
  allocate(block_source_catalog_index(n_block_built))
  allocate(block_retained_source_index(n_block_owned))
  allocate(block_migrating_source_index(n_block_migrating))

  block_source_catalog_index = -1
  block_retained_source_index = -1
  block_migrating_source_index = -1

  b_verbose = -1

  !
  ! Prefer a representative local-source block of depth at least two.
  !
  if (rank == 0) then

  do b = 1, size(block_catalog)

     if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

     d = loc_id(block_catalog(b)%root_domain+1) + 1

     if (d < 1 .or. d > size(grid)) then
        error stop &
             "test_subtree_extraction: invalid local domain mapping"
     end if

     if (subtree_depth_Domain( &
          grid(d),block_catalog(b)%root_patch) >= 2) then

        b_verbose = b
        exit

     end if

  end do

  end if

  !
  ! Fallback to the first candidate with a local source Domain.
  !
  if (rank == 0 .and. b_verbose < 1) then

     do b = 1, size(block_catalog)

        if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

        d = loc_id(block_catalog(b)%root_domain+1) + 1

        if (d < 1 .or. d > size(grid)) then
           error stop &
                "test_subtree_extraction: invalid local domain mapping"
        end if

        b_verbose = b
        exit


     end do

  end if

  if (rank == 0 .and. b_verbose < 1) then
     error stop &
          "test_subtree_extraction: no local source-domain block"
  end if

  n_block_built = 0
  n_block_owned = 0
  n_block_migrating = 0
  n_patch_total  = 0
  n_bdry_total   = 0
  n_ghost_total  = 0
  n_stencil_total = 0
  n_remote_total = 0
  n_value_total  = 0

  i_block = 0

  do b = 1, size(block_catalog)

     if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

     d = loc_id(block_catalog(b)%root_domain+1) + 1

     if (d < 1 .or. d > size(grid)) then
        error stop &
             "test_subtree_extraction: invalid local domain mapping"
     end if

     i_block = i_block + 1

     block_source_catalog_index(i_block) = b

     call test_one_subtree_extraction( &
          b, block_source(i_block), b == b_verbose, &
          n_patch_block, n_bdry_block, n_ghost_block, &
          n_stencil_block, n_remote_block, n_value_block)

     n_block_built = n_block_built + 1

     if (block_catalog(b)%owner == rank) then
        n_block_owned = n_block_owned + 1
        block_retained_source_index(n_block_owned) = i_block
     else
        n_block_migrating = n_block_migrating + 1
        block_migrating_source_index(n_block_migrating) = i_block
     end if

     n_patch_total = n_patch_total + n_patch_block
     n_bdry_total = n_bdry_total + n_bdry_block
     n_ghost_total = n_ghost_total + n_ghost_block
     n_stencil_total = n_stencil_total + n_stencil_block
     n_remote_total = n_remote_total + n_remote_block
     n_value_total = n_value_total + n_value_block

  end do

  if (n_block_built < 1) then
     error stop "test_subtree_extraction: no blocks constructed"
  end if

  if (n_block_owned + n_block_migrating /= n_block_built) then
     error stop &
          "test_subtree_extraction: block migration count mismatch"
  end if

  if (i_block /= size(block_source)) then
     error stop &
          "test_subtree_extraction: persistent block count mismatch"
  end if

  call check_persistent_test_blocks


  call test_local_block_serialization( &
       n_pack_block, n_pack_byte_total, n_pack_byte_max)

  write(6,'(/,a,i0,a)') &
       "All-block extraction summary for rank ", rank, ":"

  write(6,'(a,i0)') &
       "  local-source candidate blocks built = ", &
       n_block_built

  write(6,'(a,i0)') &
       "  retained by this rank               = ", &
       n_block_owned

  write(6,'(a,i0)') &
       "  migrating to another rank           = ", &
       n_block_migrating

  write(6,'(a,i0)') &
       "  persistent block objects            = ", &
       size(block_source)

  write(6,'(a,i0)') &
       "  locally serialized migrating blocks = ", &
       n_pack_block

  write(6,'(a,i0)') &
       "  total packed migration bytes        = ", &
       n_pack_byte_total

  write(6,'(a,i0)') &
       "  maximum packed block bytes          = ", &
       n_pack_byte_max

  write(6,'(a,i0)') &
       "  extracted regular patches           = ", &
       n_patch_total

  write(6,'(a,i0)') &
       "  compact boundary patches            = ", &
       n_bdry_total

  write(6,'(a,i0)') &
       "  compact ghost patches               = ", &
       n_ghost_total

  write(6,'(a,i0)') &
       "  explicit stencil addresses          = ", &
       n_stencil_total

  write(6,'(a,i0)') &
       "  remote-owner inter-block links      = ", &
       n_remote_total

  write(6,'(a,i0)') &
       "  inter-block scalar values checked   = ", &
       n_value_total

  write(6,'(a)') &
       "  persistent block lifetime checks passed"

  write(6,'(a)') &
       "  local block pack/unpack checks passed"

  write(6,'(a,/)') &
       "  all local-source candidate block checks passed"

end subroutine test_subtree_extraction


subroutine check_persistent_test_blocks
  ! Verify that every constructed block and all of its allocatable
  ! components remain valid after the one-block builder has returned.
  ! Also verify that the retained and migrating index sets form an
  ! exact, non-overlapping partition of the persistent block store.

  implicit none

  integer :: b
  integer :: i
  integer :: ib

  integer, allocatable :: seen(:)

  if (.not. allocated(block_source)) then
     error stop &
          "check_persistent_test_blocks: block store is not allocated"
  end if

  if (.not. allocated(block_source_catalog_index)) then
     error stop &
          "check_persistent_test_blocks: catalog map is not allocated"
  end if

  if (.not. allocated(block_retained_source_index) .or. &
       .not. allocated(block_migrating_source_index)) then
     error stop &
          "check_persistent_test_blocks: ownership sets not allocated"
  end if

  if (size(block_source_catalog_index) /= size(block_source)) then
     error stop &
          "check_persistent_test_blocks: catalog map size mismatch"
  end if

  allocate(seen(size(block_source)))
  seen = 0

  do i = 1, size(block_source)

     b = block_source_catalog_index(i)

     if (b < 1 .or. b > size(block_catalog)) then
        error stop &
             "check_persistent_test_blocks: invalid catalog index"
     end if

     if (owner(block_catalog(b)%root_domain+1) /= rank) then
        error stop &
             "check_persistent_test_blocks: source owner mismatch"
     end if

     if (block_source(i)%id /= block_catalog(b)%id .or. &
          block_source(i)%root_domain /= &
          block_catalog(b)%root_domain .or. &
          block_source(i)%root_patch /= &
          block_catalog(b)%root_patch .or. &
          block_source(i)%level /= block_catalog(b)%level) then

        error stop &
             "check_persistent_test_blocks: block identity mismatch"

     end if

     call check_block_storage(block_source(i))

  end do

  do i = 1, size(block_retained_source_index)

     ib = block_retained_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "check_persistent_test_blocks: invalid retained index"
     end if

     if (seen(ib) /= 0) then
        error stop &
             "check_persistent_test_blocks: duplicate retained index"
     end if

     b = block_source_catalog_index(ib)

     if (block_catalog(b)%owner /= rank) then
        error stop &
             "check_persistent_test_blocks: retained owner mismatch"
     end if

     seen(ib) = 1

  end do

  do i = 1, size(block_migrating_source_index)

     ib = block_migrating_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "check_persistent_test_blocks: invalid migrating index"
     end if

     if (seen(ib) /= 0) then
        error stop &
             "check_persistent_test_blocks: duplicate migrating index"
     end if

     b = block_source_catalog_index(ib)

     if (block_catalog(b)%owner == rank) then
        error stop &
             "check_persistent_test_blocks: migrating owner mismatch"
     end if

     seen(ib) = 1

  end do

  if (any(seen /= 1)) then
     error stop &
          "check_persistent_test_blocks: ownership partition mismatch"
  end if

  deallocate(seen)

end subroutine check_persistent_test_blocks


subroutine test_local_block_serialization ( &
     n_block_out, n_byte_out, n_byte_max_out)
  ! Pack and unpack every block that will migrate away from this rank.
  ! The reconstructed block is packed again and the two byte streams
  ! must match exactly.  No MPI communication is performed here.

  implicit none

  integer, intent(out) :: n_block_out
  integer, intent(out) :: n_byte_out
  integer, intent(out) :: n_byte_max_out

  integer :: i
  integer :: ib
  integer :: nbyte

  integer(int8), allocatable :: buffer_source(:)
  integer(int8), allocatable :: buffer_copy(:)

  type(Block_Data) :: block_copy

  if (.not. allocated(block_source) .or. &
       .not. allocated(block_migrating_source_index)) then

     error stop &
          "test_local_block_serialization: block store not allocated"

  end if

  n_block_out    = 0
  n_byte_out     = 0
  n_byte_max_out = 0

  do i = 1, size(block_migrating_source_index)

     ib = block_migrating_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "test_local_block_serialization: invalid block index"
     end if

     call pack_block(block_source(ib),buffer_source)
     call unpack_block(buffer_source,block_copy)
     call pack_block(block_copy,buffer_copy)

     if (size(buffer_copy) /= size(buffer_source)) then
        error stop &
             "test_local_block_serialization: packed size mismatch"
     end if

     if (any(buffer_copy /= buffer_source)) then
        error stop &
             "test_local_block_serialization: round-trip mismatch"
     end if

     nbyte = size(buffer_source)

     n_block_out = n_block_out + 1
     n_byte_out = n_byte_out + nbyte
     n_byte_max_out = max(n_byte_max_out,nbyte)

  end do

  if (n_block_out /= size(block_migrating_source_index)) then
     error stop &
          "test_local_block_serialization: tested count mismatch"
  end if

end subroutine test_local_block_serialization


subroutine test_block_migration_sizes (manifest)
  ! Associate every outgoing manifest entry with its persistent source
  ! block, compute the exact serialized byte count, and exchange only
  ! those counts. No packed block payload is communicated here.

  implicit none

  type(Block_Migration_Manifest), intent(inout) :: manifest

  integer :: b
  integer :: i
  integer :: ib
  integer, allocatable :: seen(:)
  integer, allocatable :: send_nbyte(:)

  if (.not. manifest%validated) then
     error stop &
          "test_block_migration_sizes: manifest is not validated"
  end if

  if (.not. allocated(block_source) .or. &
       .not. allocated(block_source_catalog_index) .or. &
       .not. allocated(block_migrating_source_index)) then

     error stop &
          "test_block_migration_sizes: persistent block store missing"

  end if

  if (manifest%n_send /= size(block_migrating_source_index)) then
     error stop &
          "test_block_migration_sizes: outgoing block count mismatch"
  end if

  allocate(send_nbyte(manifest%n_send))
  allocate(seen(size(block_source)))

  send_nbyte = 0
  seen       = 0

  do i = 1, manifest%n_send

     b = manifest%send_block(i)

     if (b < 1 .or. b > size(block_catalog)) then
        error stop &
             "test_block_migration_sizes: invalid catalog index"
     end if

     ib = findloc(block_source_catalog_index,b,dim=1)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "test_block_migration_sizes: source block not found"
     end if

     if (findloc(block_migrating_source_index,ib,dim=1) < 1) then
        error stop &
             "test_block_migration_sizes: block is not migrating"
     end if

     if (seen(ib) /= 0) then
        error stop &
             "test_block_migration_sizes: duplicate source block"
     end if

     if (block_source(ib)%id /= block_catalog(b)%id) then
        error stop &
             "test_block_migration_sizes: block identity mismatch"
     end if

     send_nbyte(i) = &
          packed_block_nbyte(block_source(ib))

     if (send_nbyte(i) <= 0) then
        error stop &
             "test_block_migration_sizes: invalid packed byte count"
     end if

     seen(ib) = 1

  end do

  do i = 1, size(block_migrating_source_index)

     ib = block_migrating_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "test_block_migration_sizes: invalid migrating index"
     end if

     if (seen(ib) /= 1) then
        error stop &
             "test_block_migration_sizes: migrating block omitted"
     end if

  end do

  call exchange_block_migration_sizes(manifest,send_nbyte)

  deallocate(seen)
  deallocate(send_nbyte)

end subroutine test_block_migration_sizes


subroutine test_block_migration_payloads (manifest)
  ! Pack outgoing blocks in manifest order, exchange the byte streams,
  ! and unpack received blocks into a separate temporary persistent
  ! store. Source blocks and ownership maps remain unchanged.

  implicit none

  type(Block_Migration_Manifest), intent(inout) :: manifest

  integer :: b
  integer :: i
  integer :: ib
  integer :: nbyte
  integer :: n_recv_byte
  integer :: n_send_byte
  integer :: pos
  integer, allocatable :: seen_catalog(:)

  integer(int8), allocatable :: block_buffer(:)
  integer(int8), allocatable :: check_buffer(:)
  integer(int8), allocatable :: send_payload(:)

  if (.not. manifest%validated .or. &
       .not. manifest%sizes_validated) then
     error stop &
          "test_block_migration_payloads: sizes are not validated"
  end if

  if (.not. allocated(block_source) .or. &
       .not. allocated(block_source_catalog_index) .or. &
       .not. allocated(block_migrating_source_index)) then

     error stop &
          "test_block_migration_payloads: source block store missing"

  end if

  n_send_byte = int(manifest%total_send_nbyte)
  n_recv_byte = int(manifest%total_recv_nbyte)

  allocate(send_payload(max(1,n_send_byte)))
  send_payload = 0_int8

  pos = 0

  do i = 1, manifest%n_send

     b = manifest%send_block(i)
     ib = findloc(block_source_catalog_index,b,dim=1)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "test_block_migration_payloads: source block not found"
     end if

     if (findloc(block_migrating_source_index,ib,dim=1) < 1) then
        error stop &
             "test_block_migration_payloads: source is not migrating"
     end if

     call pack_block(block_source(ib),block_buffer)

     if (size(block_buffer) /= manifest%send_nbyte(i)) then
        error stop &
             "test_block_migration_payloads: outgoing size changed"
     end if

     nbyte = size(block_buffer)

     if (pos+nbyte > n_send_byte) then
        error stop &
             "test_block_migration_payloads: outgoing buffer overflow"
     end if

     send_payload(pos+1:pos+nbyte) = block_buffer
     pos = pos + nbyte

  end do

  if (pos /= n_send_byte) then
     error stop &
          "test_block_migration_payloads: outgoing extent mismatch"
  end if

  call exchange_block_migration_payloads(manifest,send_payload)

  if (.not. manifest%payload_validated) then
     error stop &
          "test_block_migration_payloads: transport is not validated"
  end if

  if (allocated(block_received)) then
     deallocate(block_received)
  end if

  if (allocated(block_received_catalog_index)) then
     deallocate(block_received_catalog_index)
  end if

  allocate(block_received(manifest%n_recv))
  allocate(block_received_catalog_index(manifest%n_recv))
  allocate(seen_catalog(size(block_catalog)))

  block_received_catalog_index = -1
  seen_catalog = 0
  pos = 0

  do i = 1, manifest%n_recv

     b     = manifest%recv_block(i)
     nbyte = manifest%recv_nbyte(i)

     if (b < 1 .or. b > size(block_catalog)) then
        error stop &
             "test_block_migration_payloads: invalid received catalog index"
     end if

     if (seen_catalog(b) /= 0) then
        error stop &
             "test_block_migration_payloads: duplicate received block"
     end if

     if (block_catalog(b)%owner /= rank) then
        error stop &
             "test_block_migration_payloads: wrong received owner"
     end if

     if (nbyte <= 0 .or. pos+nbyte > n_recv_byte) then
        error stop &
             "test_block_migration_payloads: invalid received extent"
     end if

     call unpack_block( &
          manifest%recv_payload(pos+1:pos+nbyte), &
          block_received(i))

     if (block_received(i)%id /= block_catalog(b)%id .or. &
          block_received(i)%root_domain /= &
          block_catalog(b)%root_domain .or. &
          block_received(i)%root_patch /= &
          block_catalog(b)%root_patch .or. &
          block_received(i)%level /= block_catalog(b)%level) then

        error stop &
             "test_block_migration_payloads: received identity mismatch"

     end if

     call pack_block(block_received(i),check_buffer)

     if (size(check_buffer) /= nbyte) then
        error stop &
             "test_block_migration_payloads: received size mismatch"
     end if

     if (any(check_buffer /= &
          manifest%recv_payload(pos+1:pos+nbyte))) then
        error stop &
             "test_block_migration_payloads: received round-trip mismatch"
     end if

     block_received_catalog_index(i) = b
     seen_catalog(b) = 1
     pos = pos + nbyte

  end do

  if (pos /= n_recv_byte) then
     error stop &
          "test_block_migration_payloads: incoming extent mismatch"
  end if

  if (count(seen_catalog /= 0) /= manifest%n_recv) then
     error stop &
          "test_block_migration_payloads: received inventory mismatch"
  end if

  call check_persistent_test_blocks

  write(6,'(/,a,i0,a)') &
       "Temporary received blocks for rank ", rank, ":"
  write(6,'(a,i0)') &
       "  received block objects = ", size(block_received)
  write(6,'(a,i0)') &
       "  received packed bytes  = ", n_recv_byte
  write(6,'(a)') &
       "  received identity and byte checks passed"
  write(6,'(a,/)') &
       "  source persistent block checks still passed"

  deallocate(seen_catalog)
  deallocate(send_payload)

end subroutine test_block_migration_payloads


subroutine test_install_local_blocks
  ! Install the final-owner local store through parallel_block_mod, then
  ! verify catalogue ownership and the global MPI inventory here.

  implicit none

  integer :: b
  integer :: expected_local
  integer :: global_count
  integer :: global_weight
  integer :: ierr
  integer :: local_weight
  integer :: n_local

  integer, allocatable :: global_seen(:)
  integer, allocatable :: local_seen(:)

  call install_local_blocks(size(block_catalog),local_seen)

  n_local = size(block_local)

  allocate(global_seen(size(block_catalog)))
  global_seen = 0

  expected_local = 0
  local_weight   = 0

  do b = 1, size(block_catalog)

     if (block_catalog(b)%owner == rank) then

        expected_local = expected_local + 1
        local_weight = local_weight + block_catalog(b)%weight

        if (local_seen(b) /= 1) then
           error stop &
                "test_install_local_blocks: owned block is missing"
        end if

     else

        if (local_seen(b) /= 0) then
           error stop &
                "test_install_local_blocks: nonlocal block was installed"
        end if

     end if

  end do

  if (n_local /= expected_local) then
     error stop &
          "test_install_local_blocks: local owner count mismatch"
  end if

  if (count(local_seen /= 0) /= n_local) then
     error stop &
          "test_install_local_blocks: local inventory mismatch"
  end if

  call MPI_Allreduce( &
       n_local, global_count, 1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "test_install_local_blocks: MPI_Allreduce count failed"
  end if

  call MPI_Allreduce( &
       local_weight, global_weight, 1, MPI_INTEGER, MPI_SUM, &
       comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "test_install_local_blocks: MPI_Allreduce weight failed"
  end if

  call MPI_Allreduce( &
       local_seen, global_seen, size(local_seen), MPI_INTEGER, &
       MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "test_install_local_blocks: MPI_Allreduce inventory failed"
  end if

  if (global_count /= size(block_catalog)) then
     error stop &
          "test_install_local_blocks: global block count mismatch"
  end if

  if (global_weight /= sum(block_catalog%weight)) then
     error stop &
          "test_install_local_blocks: global block weight mismatch"
  end if

  if (any(global_seen /= 1)) then
     error stop &
          "test_install_local_blocks: global ownership is not unique"
  end if

  call check_persistent_test_blocks

  write(6,'(/,a,i0,a)') &
       "Installed local block copies for rank ", rank, ":"
  write(6,'(a,i0)') &
       "  retained source copies = ", size(block_retained_source_index)
  write(6,'(a,i0)') &
       "  received block copies  = ", size(block_received)
  write(6,'(a,i0)') &
       "  installed local blocks = ", size(block_local)
  write(6,'(a,i0)') &
       "  installed block weight = ", local_weight
  write(6,'(a)') &
       "  local deep-copy checks passed"
  write(6,'(a,/)') &
       "  source and receive stores remain available"

  if (rank == 0) then
     write(6,'(/,a,i0)') &
          "Global installed block objects verified = ", global_count
     write(6,'(a,i0)') &
          "Global installed block weight verified  = ", global_weight
     write(6,'(a,/)') &
          "Unique final-owner block installation passed"
  end if

  deallocate(global_seen)
  deallocate(local_seen)

end subroutine test_install_local_blocks


subroutine check_local_test_blocks (verbose)
  ! Validate the final-owner local block store without referring to any
  ! source, receive or migration-manifest staging allocation.

  implicit none

  logical, optional, intent(in) :: verbose

  integer :: b
  integer :: expected_local
  integer :: global_count
  integer :: global_weight
  integer :: i
  integer :: ierr
  integer :: local_count
  integer :: local_weight

  integer, allocatable :: global_seen(:)
  integer, allocatable :: local_seen(:)

  logical :: print_summary

  print_summary = .true.
  if (present(verbose)) print_summary = verbose

  if (.not. allocated(block_local) .or. &
       .not. allocated(block_local_catalog_index)) then
     error stop &
          "check_local_test_blocks: local block store is not allocated"
  end if

  if (size(block_local) /= &
       size(block_local_catalog_index)) then
     error stop &
          "check_local_test_blocks: local catalog map size mismatch"
  end if

  allocate(local_seen(size(block_catalog)))
  allocate(global_seen(size(block_catalog)))

  local_seen  = 0
  global_seen = 0
  local_count = size(block_local)
  local_weight = 0

  do i = 1, size(block_local)

     b = block_local_catalog_index(i)

     if (b < 1 .or. b > size(block_catalog)) then
        error stop &
             "check_local_test_blocks: invalid catalog index"
     end if

     if (local_seen(b) /= 0) then
        error stop &
             "check_local_test_blocks: duplicate local block"
     end if

     if (block_catalog(b)%owner /= rank) then
        error stop &
             "check_local_test_blocks: local owner mismatch"
     end if

     if (block_local(i)%id /= block_catalog(b)%id .or. &
          block_local(i)%root_domain /= &
          block_catalog(b)%root_domain .or. &
          block_local(i)%root_patch /= &
          block_catalog(b)%root_patch .or. &
          block_local(i)%level /= block_catalog(b)%level) then

        error stop &
             "check_local_test_blocks: block identity mismatch"

     end if

     call check_block_storage(block_local(i),.true.)

     local_seen(b) = 1
     local_weight = local_weight + block_catalog(b)%weight

  end do

  expected_local = 0

  do b = 1, size(block_catalog)
     if (block_catalog(b)%owner == rank) then
        expected_local = expected_local + 1
     end if
  end do

  if (size(block_local) /= expected_local) then
     error stop &
          "check_local_test_blocks: expected local count mismatch"
  end if

  call MPI_Allreduce( &
       local_count, global_count, 1, &
       MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "check_local_test_blocks: MPI_Allreduce count failed"
  end if

  call MPI_Allreduce( &
       local_weight, global_weight, 1, MPI_INTEGER, MPI_SUM, &
       comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "check_local_test_blocks: MPI_Allreduce weight failed"
  end if

  call MPI_Allreduce( &
       local_seen, global_seen, size(local_seen), MPI_INTEGER, &
       MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop &
          "check_local_test_blocks: MPI_Allreduce inventory failed"
  end if

  if (global_count /= size(block_catalog)) then
     error stop &
          "check_local_test_blocks: global count mismatch"
  end if

  if (global_weight /= sum(block_catalog%weight)) then
     error stop &
          "check_local_test_blocks: global weight mismatch"
  end if

  if (any(global_seen /= 1)) then
     error stop &
          "check_local_test_blocks: nonunique global ownership"
  end if

  if (print_summary) then
     write(6,'(/,a,i0,a)') &
          "Standalone local block store for rank ", rank, ":"
     write(6,'(a,i0)') &
          "  final-owner blocks = ", size(block_local)
     write(6,'(a,i0)') &
          "  final-owner weight = ", local_weight
     write(6,'(a)') &
          "  component and serialization checks passed"
     write(6,'(a,/)') &
          "  unique global inventory check passed"
  end if

  if (print_summary .and. rank == 0) then
     write(6,'(/,a,i0)') &
          "Standalone global block objects verified = ", global_count
     write(6,'(a,i0)') &
          "Standalone global block weight verified  = ", global_weight
     write(6,'(a,/)') &
          "Final-owner block store is self-contained"
  end if

  deallocate(global_seen)
  deallocate(local_seen)

end subroutine check_local_test_blocks


subroutine finalize_test_block_migration (manifest)
  ! Verify the final-owner store, release every migration staging store,
  ! and verify the final-owner store again after cleanup.

  implicit none

  type(Block_Migration_Manifest), intent(inout) :: manifest

  call check_local_test_blocks(.false.)

  call clear_block_staging
  call clear_block_migration_manifest(manifest)

  if (allocated(manifest%send_count) .or. &
       allocated(manifest%recv_count) .or. &
       allocated(manifest%send_displ) .or. &
       allocated(manifest%recv_displ) .or. &
       allocated(manifest%send_block) .or. &
       allocated(manifest%recv_block) .or. &
       allocated(manifest%send_nbyte) .or. &
       allocated(manifest%recv_nbyte) .or. &
       allocated(manifest%send_byte_count) .or. &
       allocated(manifest%recv_byte_count) .or. &
       allocated(manifest%send_byte_displ) .or. &
       allocated(manifest%recv_byte_displ) .or. &
       allocated(manifest%recv_payload)) then

     error stop &
          "finalize_test_block_migration: manifest cleanup failed"

  end if

  call check_local_test_blocks(.true.)

  write(6,'(a,i0,a)') &
       "Migration staging cleanup passed on rank ", rank, "."

end subroutine finalize_test_block_migration




subroutine test_one_subtree_extraction ( &
     b_test, block_test, verbose, &
     n_patch_out, n_bdry_out, n_ghost_out, &
     n_stencil_out, n_remote_out, n_value_out)
  ! Extract one candidate subtree and verify:
  !
  !   - patch topology and compact interior storage;
  !   - copied node/scalar/vector data;
  !   - block-local neighbour classification;
  !   - explicit boundary-link catalogue;
  !   - complete compact boundary storage, including stencil closure;
  !   - unified effective-stencil and nominal inter-block ghost storage;
  !   - unique source block/owner and ghost ID for every NGB_BLOCK link;
  !   - block-local neighbour references;
  !   - explicit compact scalar stencil addressing.
  !
  ! Regular source patches required by effective stencil addresses but
  ! lying outside the extracted subtree are copied into compact ghost
  ! storage. NGB_BLOCK addressing is complete, but ghost field values
  ! are still copied directly until communication is constructed.

  implicit none

  integer, intent(in) :: b_test
  type(Block_Data), intent(out) :: block_test
  logical, intent(in) :: verbose

  integer, intent(out) :: n_patch_out
  integer, intent(out) :: n_bdry_out
  integer, intent(out) :: n_ghost_out
  integer, intent(out) :: n_stencil_out
  integer, intent(out) :: n_remote_out
  integer, intent(out) :: n_value_out

  integer :: b_src
  integer :: b_missing
  integer :: c, d, i, ib, is, j
  integer :: p_root

  integer :: n_old, n_new
  integer :: n_leaf_old, n_leaf_new
  integer :: depth_old, depth_new

  integer :: n_node_storage
  integer :: n_patch_field

  integer :: n_ngb_internal
  integer :: n_ngb_block
  integer :: n_ngb_domain
  integer :: n_ngb_adapt
  integer :: n_ngb_other

  integer :: n_bdry_local
  integer :: n_bdry_direct
  integer :: n_bdry_closure
  integer :: n_bdry_required
  integer :: n_bdry_node_total
  integer :: n_bdry_node_unique
  integer :: n_bdry_node_max

  integer :: n_block_source_local
  integer :: n_block_source_remote
  integer :: n_source_match
  integer :: source_block

  integer :: n_stencil_built
  integer :: n_stencil_patch
  integer :: n_stencil_bdry
  integer :: n_stencil_ghost
  integer :: n_stencil_block
  integer :: n_block_value_checked
  integer :: n_stencil_unresolved

  integer :: n_ghost
  integer :: n_ghost_node
  integer :: ghost_id
  integer :: source_patch

  integer :: n_unresolved_patch
  integer :: n_unresolved_bdry
  integer :: n_unresolved_none

  integer :: unresolved_patch
  integer :: unresolved_bdry

  integer :: target_patch
  integer :: target_bdry
  integer :: target_offset
  integer :: ghost_offset
  integer :: n_mapped
  integer :: q

  integer :: idx_src

  integer :: p_old
  integer :: p_chd_old
  integer :: p_chd_new
  integer :: p_ngb_old

  integer :: old_start
  integer :: new_start

  integer :: v_scalar
  integer :: v_vector
  integer :: k_test
  integer :: mult_scalar
  integer :: mult_vector

  integer :: offs_src(0:N_BDRY)
  integer :: dims_src(2,N_BDRY)

  integer, allocatable :: old_to_new(:)
  integer, allocatable :: old_elts_start(:)

  integer, allocatable :: bdry_required(:)
  integer, allocatable :: bdry_closure(:)
  integer, allocatable :: ghost_patch(:)

  real(dp) :: val_src
  real(dp) :: val_blk

  real(dp), allocatable :: scalar_copy(:)
  real(dp), allocatable :: vector_copy(:)

  type(Patch), allocatable :: patch_copy(:)
  type(Coord), allocatable :: node_copy(:)

  logical :: already_present

  if (owner(block_catalog(b_test)%root_domain+1) /= rank) then
     error stop &
          "test_subtree_extraction: source domain is not local"
  end if

  d = loc_id(block_catalog(b_test)%root_domain+1) + 1

  if (d < 1 .or. d > size(grid)) then
     error stop "test_subtree_extraction: invalid local domain"
  end if

  p_root   = block_catalog(b_test)%root_patch
  depth_old = subtree_depth_Domain(grid(d),p_root)

  !
  ! ===============================================================
  ! Extract and compact patch tree.
  ! ===============================================================
  !
  call extract_subtree_patches_Domain( &
       grid(d), p_root, patch_copy, old_to_new)

  n_old = count_subtree_patches_Domain(grid(d),p_root)
  n_new = size(patch_copy)

  if (n_new /= n_old) then
     error stop "test_subtree_extraction: patch count mismatch"
  end if

  if (old_to_new(p_root) /= 0) then
     error stop "test_subtree_extraction: extracted root is not patch zero"
  end if

  allocate(old_elts_start(size(patch_copy)))

  do i = 1, size(patch_copy)
     old_elts_start(i) = patch_copy(i)%elts_start
  end do

  do i = 1, size(patch_copy)

     if (old_elts_start(i) < 0) then
        error stop "test_subtree_extraction: invalid original elts_start"
     end if

     if (old_elts_start(i) + PATCH_SIZE**2 > grid(d)%node%length) then
        error stop &
             "test_subtree_extraction: original node storage out of bounds"
     end if

  end do

  call compact_subtree_storage_Domain(patch_copy)

  do i = 1, size(patch_copy)

     if (patch_copy(i)%elts_start /= (i-1)*PATCH_SIZE**2) then
        error stop &
             "test_subtree_extraction: incorrect compact elts_start"
     end if

  end do

  n_node_storage = size(patch_copy) * PATCH_SIZE**2

  !
  ! ===============================================================
  ! Copy and verify interior coordinates.
  ! ===============================================================
  !
  call copy_subtree_nodes_Domain( &
       grid(d), patch_copy, old_elts_start, node_copy)

  if (size(node_copy) /= n_node_storage) then
     error stop &
          "test_subtree_extraction: incorrect copied node storage size"
  end if

  do i = 1, size(patch_copy)

     old_start = old_elts_start(i)
     new_start = patch_copy(i)%elts_start

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%x - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%x)) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: node x-coordinate mismatch"

     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%y - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%y)) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: node y-coordinate mismatch"

     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%z - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%z)) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: node z-coordinate mismatch"

     end if

  end do

  !
  ! ===============================================================
  ! Copy and verify one scalar field.
  ! ===============================================================
  !
  v_scalar    = scalars(1)
  mult_scalar = MULT(v_scalar)
  k_test      = max(1,zmin)

  if (mult_scalar /= 1) then
     error stop &
          "test_subtree_extraction: unexpected scalar multiplier"
  end if

  call copy_subtree_field_Domain( &
       patch_copy, old_elts_start, mult_scalar, &
       sol(v_scalar,k_test)%data(d)%elts, scalar_copy)

  if (size(scalar_copy) /= mult_scalar*n_node_storage) then
     error stop &
          "test_subtree_extraction: incorrect scalar storage size"
  end if

  n_patch_field = mult_scalar * PATCH_SIZE**2

  do i = 1, size(patch_copy)

     old_start = mult_scalar * old_elts_start(i)
     new_start = mult_scalar * patch_copy(i)%elts_start

     if (maxval(abs( &
          scalar_copy(new_start+1:new_start+n_patch_field) - &
          sol(v_scalar,k_test)%data(d)%elts( &
          old_start+1:old_start+n_patch_field))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: scalar field copy mismatch"

     end if

  end do

  !
  ! ===============================================================
  ! Copy and verify one vector field.
  ! ===============================================================
  !
  v_vector    = S_VELO
  mult_vector = MULT(v_vector)

  if (mult_vector /= EDGE) then
     error stop &
          "test_subtree_extraction: unexpected vector multiplier"
  end if

  call copy_subtree_field_Domain( &
       patch_copy, old_elts_start, mult_vector, &
       sol(v_vector,k_test)%data(d)%elts, vector_copy)

  if (size(vector_copy) /= mult_vector*n_node_storage) then
     error stop &
          "test_subtree_extraction: incorrect vector storage size"
  end if

  n_patch_field = mult_vector * PATCH_SIZE**2

  do i = 1, size(patch_copy)

     old_start = mult_vector * old_elts_start(i)
     new_start = mult_vector * patch_copy(i)%elts_start

     if (maxval(abs( &
          vector_copy(new_start+1:new_start+n_patch_field) - &
          sol(v_vector,k_test)%data(d)%elts( &
          old_start+1:old_start+n_patch_field))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: vector field copy mismatch"

     end if

  end do

  !
  ! ===============================================================
  ! Verify child links, leaf count and depth.
  ! ===============================================================
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_CHDRN

        p_chd_old = grid(d)%patch%elts(p_old+1)%children(c)
        p_chd_new = patch_copy(old_to_new(p_old)+1)%children(c)

        if (p_chd_old <= 0) then

           if (p_chd_new /= 0) then
              error stop &
                   "test_subtree_extraction: unexpected child link"
           end if

        else if (grid(d)%patch%elts(p_chd_old+1)%deleted) then

           if (p_chd_new /= 0) then
              error stop &
                   "test_subtree_extraction: deleted child copied"
           end if

        else

           if (old_to_new(p_chd_old) < 0) then
              error stop &
                   "test_subtree_extraction: child missing from map"
           end if

           if (p_chd_new /= old_to_new(p_chd_old)) then
              error stop &
                   "test_subtree_extraction: incorrect child renumbering"
           end if

        end if

     end do

  end do

  n_leaf_old = 0

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     if (all(grid(d)%patch%elts(p_old+1)%children == 0)) then
        n_leaf_old = n_leaf_old + 1
     end if

  end do

  n_leaf_new = 0

  do i = 1, size(patch_copy)

     if (all(patch_copy(i)%children == 0)) then
        n_leaf_new = n_leaf_new + 1
     end if

  end do

  if (n_leaf_new /= n_leaf_old) then
     error stop "test_subtree_extraction: leaf count mismatch"
  end if

  depth_new = copied_depth(0)

  if (depth_new /= depth_old) then
     error stop "test_subtree_extraction: subtree depth mismatch"
  end if

  !
  ! ===============================================================
  ! Classify source neighbour links.
  ! ===============================================================
  !
  allocate(block_test%neigh_class(N_BDRY,size(patch_copy)))

  block_test%neigh_class = NGB_OTHER

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

        if (p_ngb_old > 0) then

           if (p_ngb_old >= grid(d)%patch%length) then
              error stop &
                   "test_subtree_extraction: invalid positive neighbour"
           end if

           if (old_to_new(p_ngb_old) >= 0) then

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_INTERNAL

           else

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_BLOCK

           end if

        else if (p_ngb_old < 0) then

           b_src = -p_ngb_old

           if (b_src >= grid(d)%bdry_patch%length) then
              error stop &
                   "test_subtree_extraction: invalid source boundary"
           end if

           if (grid(d)%bdry_patch%elts(b_src+1)%side > 0) then

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_DOMAIN

           else

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_ADAPT

           end if

        else

           block_test%neigh_class( &
                c,old_to_new(p_old)+1) = NGB_OTHER

        end if

     end do

  end do

  n_ngb_internal = count(block_test%neigh_class == NGB_INTERNAL)
  n_ngb_block    = count(block_test%neigh_class == NGB_BLOCK)
  n_ngb_domain   = count(block_test%neigh_class == NGB_DOMAIN)
  n_ngb_adapt    = count(block_test%neigh_class == NGB_ADAPT)
  n_ngb_other    = count(block_test%neigh_class == NGB_OTHER)

  if (n_ngb_internal + n_ngb_block + n_ngb_domain + &
       n_ngb_adapt + n_ngb_other /= &
       N_BDRY*size(patch_copy)) then

     error stop &
          "test_subtree_extraction: neighbour count mismatch"

  end if

  !
  ! ===============================================================
  ! Build explicit local boundary-link catalogue.
  ! ===============================================================
  !
  n_bdry_local = n_ngb_block + n_ngb_domain + n_ngb_adapt

  allocate(block_test%block_bdry(n_bdry_local))

  ib = 0

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL, NGB_OTHER)

           cycle

        case (NGB_BLOCK)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (p_ngb_old <= 0 .or. &
                p_ngb_old >= grid(d)%patch%length) then

              error stop &
                   "test_subtree_extraction: invalid inter-block neighbour"

           end if

           ib = ib + 1

           block_test%block_bdry(ib)%patch = old_to_new(p_old)
           block_test%block_bdry(ib)%side  = c
           block_test%block_bdry(ib)%class = NGB_BLOCK

           block_test%block_bdry(ib)%root_domain = &
                block_catalog(b_test)%root_domain

           block_test%block_bdry(ib)%neigh_patch = p_ngb_old

           block_test%block_bdry(ib)%source_block = -1
           block_test%block_bdry(ib)%source_block_id = -1
           block_test%block_bdry(ib)%source_owner = -1
           block_test%block_bdry(ib)%ghost_id     = -1
           block_test%block_bdry(ib)%storage_id = -1

        case (NGB_DOMAIN, NGB_ADAPT)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (p_ngb_old >= 0) then
              error stop &
                   "test_subtree_extraction: existing boundary is not negative"
           end if

           b_src = -p_ngb_old

           if (b_src >= grid(d)%bdry_patch%length) then
              error stop &
                   "test_subtree_extraction: invalid source bdry_patch"
           end if

           ib = ib + 1

           block_test%block_bdry(ib)%patch = old_to_new(p_old)
           block_test%block_bdry(ib)%side  = c

           block_test%block_bdry(ib)%class = &
                block_test%neigh_class(c,old_to_new(p_old)+1)

           block_test%block_bdry(ib)%root_domain = &
                block_catalog(b_test)%root_domain

           block_test%block_bdry(ib)%source_bdry = b_src

           block_test%block_bdry(ib)%elts_start = &
                grid(d)%bdry_patch%elts(b_src+1)%elts_start

           block_test%block_bdry(ib)%bdry_side = &
                grid(d)%bdry_patch%elts(b_src+1)%side

           block_test%block_bdry(ib)%bdry_neigh = &
                grid(d)%bdry_patch%elts(b_src+1)%neigh

           call get_bdry_dims_Domain( &
                grid(d), b_src, block_test%block_bdry(ib)%dims)

           block_test%block_bdry(ib)%n_node = &
                BDRY_THICKNESS * PATCH_SIZE

           block_test%block_bdry(ib)%storage_id = -1

        case default

           error stop &
                "test_subtree_extraction: unexpected neighbour class"

        end select

     end do

  end do

  if (ib /= n_bdry_local) then
     error stop &
          "test_subtree_extraction: local boundary count mismatch"
  end if

  !
  ! ===============================================================
  ! Map every inter-block link to its unique source patch, candidate
  ! source block, and prospective source owner.
  ! ===============================================================
  !
  n_block_source_local  = 0
  n_block_source_remote = 0

  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_block  = -1
     n_source_match = 0

     do i = 1, size(block_catalog)

        if (block_catalog(i)%root_domain /= &
             block_catalog(b_test)%root_domain) cycle

        if (.not. patch_in_subtree( &
             block_catalog(i)%root_patch, &
             block_test%block_bdry(ib)%neigh_patch)) cycle

        source_block   = i
        n_source_match = n_source_match + 1

     end do

     if (n_source_match /= 1) then
        error stop &
             "test_subtree_extraction: nonunique inter-block source"
     end if

     if (source_block == b_test) then
        error stop &
             "test_subtree_extraction: inter-block source is current block"
     end if

     block_test%block_bdry(ib)%source_block = source_block
     block_test%block_bdry(ib)%source_block_id = &
          block_catalog(source_block)%id
     block_test%block_bdry(ib)%source_owner = &
          block_catalog(source_block)%owner

     if (block_test%block_bdry(ib)%source_owner < 0 .or. &
          block_test%block_bdry(ib)%source_owner >= n_process) then

        error stop &
             "test_subtree_extraction: invalid inter-block source owner"

     end if

     if (block_test%block_bdry(ib)%source_owner == &
          block_catalog(b_test)%owner) then
        n_block_source_local = n_block_source_local + 1
     else
        n_block_source_remote = n_block_source_remote + 1
     end if

  end do

  if (n_block_source_local + n_block_source_remote /= &
       n_ngb_block) then

     error stop &
          "test_subtree_extraction: inter-block mapping count mismatch"

  end if

  !
  ! ===============================================================
  ! Build complete boundary-storage requirement.
  !
  ! First collect directly referenced boundaries, then add any source
  ! boundary records required by effective comp_offs3 stencil
  ! addresses.
  ! ===============================================================
  !
  allocate(bdry_required(max(1,grid(d)%bdry_patch%length)))
  allocate(bdry_closure(max(1,grid(d)%bdry_patch%length)))

  bdry_required = -1
  bdry_closure  = -1

  n_bdry_direct   = 0
  n_bdry_closure  = 0
  n_bdry_required = 0

  !
  ! Distinct directly referenced boundaries.
  !
  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_DOMAIN .and. &
          block_test%block_bdry(ib)%class /= NGB_ADAPT) cycle

     b_src = block_test%block_bdry(ib)%source_bdry

     already_present = .false.

     do j = 1, n_bdry_required

        if (bdry_required(j) == b_src) then
           already_present = .true.
           exit
        end if

     end do

     if (.not. already_present) then

        n_bdry_required = n_bdry_required + 1
        n_bdry_direct   = n_bdry_direct + 1

        bdry_required(n_bdry_required) = b_src

     end if

  end do

  !
  ! Add stencil-closure boundaries.
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     call comp_offs3( &
          grid(d), p_old, offs_src, dims_src)

     do c = 1, N_BDRY

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL, NGB_DOMAIN, NGB_ADAPT)

           idx_src = grid(d)%patch%elts(p_old+1)%elts_start + &
                offs_src(c)

        case default

           cycle

        end select

        !
        ! Already contained in an extracted patch?
        !
        target_patch = -1

        do i = 0, grid(d)%patch%length-1

           if (old_to_new(i) < 0) cycle

           old_start = grid(d)%patch%elts(i+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + PATCH_SIZE**2) then

              target_patch = old_to_new(i)
              exit

           end if

        end do

        if (target_patch >= 0) cycle

        !
        ! Already contained in a required boundary?
        !
        target_bdry = -1

        do j = 1, n_bdry_required

           b_src = bdry_required(j)

           old_start = &
                grid(d)%bdry_patch%elts(b_src+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + &
                BDRY_THICKNESS*PATCH_SIZE) then

              target_bdry = j
              exit

           end if

        end do

        if (target_bdry >= 0) cycle

        !
        ! Search complete source boundary catalogue.
        !
        b_missing = -1

        do is = 0, grid(d)%bdry_patch%length-1

           old_start = &
                grid(d)%bdry_patch%elts(is+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + &
                BDRY_THICKNESS*PATCH_SIZE) then

              b_missing = is
              exit

           end if

        end do

        !
        ! Remaining cases are outside-block regular patches or
        ! implicit source-domain ghost-layout addresses.
        !
        if (b_missing < 0) cycle

        already_present = .false.

        do j = 1, n_bdry_required

           if (bdry_required(j) == b_missing) then
              already_present = .true.
              exit
           end if

        end do

        if (already_present) cycle

        n_bdry_required = n_bdry_required + 1
        n_bdry_closure  = n_bdry_closure + 1

        bdry_required(n_bdry_required) = b_missing
        bdry_closure(n_bdry_closure)   = b_missing

     end do

  end do

  !
  ! Validate closure list generically.
  !
  do i = 1, n_bdry_closure

     b_src = bdry_closure(i)

     if (b_src < 0 .or. &
          b_src >= grid(d)%bdry_patch%length) then

        error stop &
             "test_subtree_extraction: invalid closure boundary"

     end if

     do j = 1, i-1

        if (bdry_closure(j) == b_src) then
           error stop &
                "test_subtree_extraction: duplicate closure boundary"
        end if

     end do

     do j = 1, n_bdry_direct

        if (bdry_required(j) == b_src) then
           error stop &
                "test_subtree_extraction: closure boundary already direct"
        end if

     end do

  end do

  if (n_bdry_required /= n_bdry_direct + n_bdry_closure) then
     error stop &
          "test_subtree_extraction: required boundary count mismatch"
  end if

  !
  ! ===============================================================
  ! Construct complete compact boundary-storage catalogue.
  ! ===============================================================
  !
  allocate(block_test%bdry_storage(n_bdry_required))

  n_bdry_node_unique = 0

  do is = 1, n_bdry_required

     b_src = bdry_required(is)

     block_test%bdry_storage(is)%source_bdry = b_src

     block_test%bdry_storage(is)%elts_start = &
          grid(d)%bdry_patch%elts(b_src+1)%elts_start

     call get_bdry_dims_Domain( &
          grid(d), b_src, block_test%bdry_storage(is)%dims)

     block_test%bdry_storage(is)%n_node = &
          BDRY_THICKNESS * PATCH_SIZE

     if (block_test%bdry_storage(is)%n_node <= 0) then
        error stop &
             "test_subtree_extraction: invalid boundary storage size"
     end if

     block_test%bdry_storage(is)%local_start = &
          n_bdry_node_unique

     n_bdry_node_unique = n_bdry_node_unique + &
          block_test%bdry_storage(is)%n_node

  end do

  !
  ! Assign storage IDs to directly referenced boundary links.
  !
  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class == NGB_BLOCK) cycle

     do is = 1, size(block_test%bdry_storage)

        if (block_test%bdry_storage(is)%source_bdry == &
             block_test%block_bdry(ib)%source_bdry) then

           block_test%block_bdry(ib)%storage_id = is
           exit

        end if

     end do

     if (block_test%block_bdry(ib)%storage_id < 1) then
        error stop &
             "test_subtree_extraction: missing boundary storage ID"
     end if

  end do

  !
  ! Verify compact storage layout and uniqueness.
  !
  do is = 1, size(block_test%bdry_storage)

     if (is == 1) then

        if (block_test%bdry_storage(is)%local_start /= 0) then
           error stop &
                "test_subtree_extraction: invalid first boundary offset"
        end if

     else

        if (block_test%bdry_storage(is)%local_start /= &
             block_test%bdry_storage(is-1)%local_start + &
             block_test%bdry_storage(is-1)%n_node) then

           error stop &
                "test_subtree_extraction: noncompact boundary storage"

        end if

     end if

     if (block_test%bdry_storage(is)%elts_start < 0 .or. &
          block_test%bdry_storage(is)%elts_start + &
          block_test%bdry_storage(is)%n_node > &
          grid(d)%node%length) then

        error stop &
             "test_subtree_extraction: boundary source range invalid"

     end if

     do j = is+1, size(block_test%bdry_storage)

        if (block_test%bdry_storage(is)%source_bdry == &
             block_test%bdry_storage(j)%source_bdry) then

           error stop &
                "test_subtree_extraction: duplicate boundary storage"

        end if

     end do

  end do

  !
  ! ===============================================================
  ! Copy complete compact boundary data.
  ! ===============================================================
  !
  allocate(block_test%bdry_node(n_bdry_node_unique))

  allocate(block_test%bdry_scalar( &
       mult_scalar*n_bdry_node_unique))

  allocate(block_test%bdry_vector( &
       mult_vector*n_bdry_node_unique))

  do is = 1, size(block_test%bdry_storage)

     old_start = block_test%bdry_storage(is)%elts_start
     new_start = block_test%bdry_storage(is)%local_start

     block_test%bdry_node( &
          new_start+1 : &
          new_start+block_test%bdry_storage(is)%n_node) = &
          grid(d)%node%elts( &
          old_start+1 : &
          old_start+block_test%bdry_storage(is)%n_node)

     block_test%bdry_scalar( &
          mult_scalar*new_start+1 : &
          mult_scalar*(new_start + &
          block_test%bdry_storage(is)%n_node)) = &
          sol(v_scalar,k_test)%data(d)%elts( &
          mult_scalar*old_start+1 : &
          mult_scalar*(old_start + &
          block_test%bdry_storage(is)%n_node))

     block_test%bdry_vector( &
          mult_vector*new_start+1 : &
          mult_vector*(new_start + &
          block_test%bdry_storage(is)%n_node)) = &
          sol(v_vector,k_test)%data(d)%elts( &
          mult_vector*old_start+1 : &
          mult_vector*(old_start + &
          block_test%bdry_storage(is)%n_node))

  end do

  !
  ! Verify copied boundary data.
  !
  do is = 1, size(block_test%bdry_storage)

     old_start = block_test%bdry_storage(is)%elts_start
     new_start = block_test%bdry_storage(is)%local_start

     if (maxval(abs( &
          block_test%bdry_node( &
          new_start+1 : &
          new_start+block_test%bdry_storage(is)%n_node)%x - &
          grid(d)%node%elts( &
          old_start+1 : &
          old_start+block_test%bdry_storage(is)%n_node)%x)) > &
          0.0_dp) then

        error stop &
             "test_subtree_extraction: boundary coordinate mismatch"

     end if

     if (maxval(abs( &
          block_test%bdry_scalar( &
          mult_scalar*new_start+1 : &
          mult_scalar*(new_start + &
          block_test%bdry_storage(is)%n_node)) - &
          sol(v_scalar,k_test)%data(d)%elts( &
          mult_scalar*old_start+1 : &
          mult_scalar*(old_start + &
          block_test%bdry_storage(is)%n_node)))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: boundary scalar mismatch"

     end if

     if (maxval(abs( &
          block_test%bdry_vector( &
          mult_vector*new_start+1 : &
          mult_vector*(new_start + &
          block_test%bdry_storage(is)%n_node)) - &
          sol(v_vector,k_test)%data(d)%elts( &
          mult_vector*old_start+1 : &
          mult_vector*(old_start + &
          block_test%bdry_storage(is)%n_node)))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: boundary vector mismatch"

     end if

  end do

  !
  ! Boundary-storage statistics for directly referenced links.
  !
  n_bdry_node_total = 0
  n_bdry_node_max   = 0

  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_DOMAIN .and. &
          block_test%block_bdry(ib)%class /= NGB_ADAPT) cycle

     n_bdry_node_total = n_bdry_node_total + &
          block_test%block_bdry(ib)%n_node

     n_bdry_node_max = max( &
          n_bdry_node_max, &
          block_test%block_bdry(ib)%n_node)

  end do

  !
  ! ===============================================================
  ! Identify distinct regular source patches required by effective
  ! stencil addresses but lying outside the extracted block.
  ! ===============================================================
  !
  allocate(ghost_patch(max(1,grid(d)%patch%length)))

  ghost_patch = -1
  n_ghost = 0

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     call comp_offs3( &
          grid(d), p_old, offs_src, dims_src)

     do c = 1, N_BDRY

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL, NGB_DOMAIN, NGB_ADAPT)

           idx_src = grid(d)%patch%elts(p_old+1)%elts_start + &
                offs_src(c)

        case default

           cycle

        end select

        target_patch = -1

        do i = 0, grid(d)%patch%length-1

           if (old_to_new(i) < 0) cycle

           old_start = grid(d)%patch%elts(i+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + PATCH_SIZE**2) then

              target_patch = old_to_new(i)
              exit

           end if

        end do

        if (target_patch >= 0) cycle

        target_bdry = -1

        do is = 1, size(block_test%bdry_storage)

           old_start = block_test%bdry_storage(is)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + &
                block_test%bdry_storage(is)%n_node) then

              target_bdry = is
              exit

           end if

        end do

        if (target_bdry >= 0) cycle

        source_patch = -1

        do i = 0, grid(d)%patch%length-1

           old_start = grid(d)%patch%elts(i+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + PATCH_SIZE**2) then

              source_patch = i
              exit

           end if

        end do

        if (source_patch < 0) cycle

        if (old_to_new(source_patch) >= 0) then
           error stop &
                "test_subtree_extraction: ghost source is inside block"
        end if

        already_present = .false.

        do j = 1, n_ghost

           if (ghost_patch(j) == source_patch) then
              already_present = .true.
              exit
           end if

        end do

        if (.not. already_present) then
           n_ghost = n_ghost + 1
           ghost_patch(n_ghost) = source_patch
        end if

     end do

  end do

  !
  ! Add all nominal NGB_BLOCK source patches to the same catalogue.
  ! Deduplicate them against effective-stencil ghost patches and
  ! against repeated inter-block links.
  !
  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_patch = block_test%block_bdry(ib)%neigh_patch
     already_present = .false.

     do j = 1, n_ghost

        if (ghost_patch(j) == source_patch) then
           already_present = .true.
           exit
        end if

     end do

     if (.not. already_present) then
        n_ghost = n_ghost + 1
        ghost_patch(n_ghost) = source_patch
     end if

  end do

  !
  ! Determine the unique source block and owner for every distinct
  ! ghost patch, including ghosts discovered only through effective
  ! stencil addressing.
  !
  allocate(block_test%ghost_storage(n_ghost))

  n_ghost_node = 0

  do ghost_id = 1, n_ghost

     source_block   = -1
     n_source_match = 0

     do i = 1, size(block_catalog)

        if (block_catalog(i)%root_domain /= &
             block_catalog(b_test)%root_domain) cycle

        if (.not. patch_in_subtree( &
             block_catalog(i)%root_patch, &
             ghost_patch(ghost_id))) cycle

        source_block   = i
        n_source_match = n_source_match + 1

     end do

     if (n_source_match /= 1) then
        error stop &
             "test_subtree_extraction: nonunique ghost source block"
     end if

     if (source_block == b_test) then
        error stop &
             "test_subtree_extraction: ghost source is current block"
     end if

     block_test%ghost_storage(ghost_id)%source_domain = &
          block_catalog(b_test)%root_domain

     block_test%ghost_storage(ghost_id)%source_patch = &
          ghost_patch(ghost_id)

     block_test%ghost_storage(ghost_id)%source_block = &
          source_block

     block_test%ghost_storage(ghost_id)%source_block_id = &
          block_catalog(source_block)%id

     block_test%ghost_storage(ghost_id)%source_owner = &
          block_catalog(source_block)%owner

     if (block_test%ghost_storage(ghost_id)%source_owner < 0 .or. &
          block_test%ghost_storage(ghost_id)%source_owner >= &
          n_process) then

        error stop &
             "test_subtree_extraction: invalid ghost source owner"

     end if

     block_test%ghost_storage(ghost_id)%elts_start = &
          grid(d)%patch%elts(ghost_patch(ghost_id)+1)%elts_start

     block_test%ghost_storage(ghost_id)%local_start = &
          n_ghost_node

     block_test%ghost_storage(ghost_id)%n_node = &
          PATCH_SIZE**2

     n_ghost_node = n_ghost_node + PATCH_SIZE**2

  end do

  !
  ! Assign every nominal inter-block link its compact ghost ID.
  !
  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     do ghost_id = 1, n_ghost

        if (block_test%ghost_storage(ghost_id)%source_patch == &
             block_test%block_bdry(ib)%neigh_patch) then

           block_test%block_bdry(ib)%ghost_id = ghost_id
           exit

        end if

     end do

     if (block_test%block_bdry(ib)%ghost_id < 1) then
        error stop &
             "test_subtree_extraction: missing inter-block ghost ID"
     end if

     ghost_id = block_test%block_bdry(ib)%ghost_id

     if (block_test%ghost_storage(ghost_id)%source_block /= &
          block_test%block_bdry(ib)%source_block .or. &
          block_test%ghost_storage(ghost_id)%source_block_id /= &
          block_test%block_bdry(ib)%source_block_id .or. &
          block_test%ghost_storage(ghost_id)%source_owner /= &
          block_test%block_bdry(ib)%source_owner) then

        error stop &
             "test_subtree_extraction: ghost source metadata mismatch"

     end if

  end do

  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_block = block_test%block_bdry(ib)%source_block

     if (source_block < 1 .or. &
          source_block > size(block_catalog)) then

        error stop &
             "test_subtree_extraction: invalid source catalog index"

     end if

     if (block_test%block_bdry(ib)%source_block_id /= &
          block_catalog(source_block)%id) then

        error stop &
             "test_subtree_extraction: source block ID mismatch"

     end if

  end do

  do ghost_id = 1, size(block_test%ghost_storage)

     source_block = &
          block_test%ghost_storage(ghost_id)%source_block

     if (source_block < 1 .or. &
          source_block > size(block_catalog)) then

        error stop &
             "test_subtree_extraction: invalid ghost catalog index"

     end if

     if (block_test%ghost_storage(ghost_id)%source_block_id /= &
          block_catalog(source_block)%id) then

        error stop &
             "test_subtree_extraction: ghost block ID mismatch"

     end if

  end do

  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     do j = ib+1, size(block_test%block_bdry)

        if (block_test%block_bdry(j)%class /= NGB_BLOCK) cycle

        if ((block_test%block_bdry(ib)%neigh_patch == &
             block_test%block_bdry(j)%neigh_patch) .neqv. &
             (block_test%block_bdry(ib)%ghost_id == &
             block_test%block_bdry(j)%ghost_id)) then

           error stop &
                "test_subtree_extraction: inconsistent ghost deduplication"

        end if

     end do

  end do

  allocate(block_test%ghost_node(n_ghost_node))

  allocate(block_test%ghost_scalar( &
       mult_scalar*n_ghost_node))

  allocate(block_test%ghost_vector( &
       mult_vector*n_ghost_node))

  !
  ! Copy temporary ghost data from the source domain. Eventually the
  ! field data will be supplied by inter-block communication.
  !
  do ghost_id = 1, n_ghost

     old_start = &
          block_test%ghost_storage(ghost_id)%elts_start

     new_start = &
          block_test%ghost_storage(ghost_id)%local_start

     block_test%ghost_node( &
          new_start+1 : new_start+PATCH_SIZE**2) = &
          grid(d)%node%elts( &
          old_start+1 : old_start+PATCH_SIZE**2)

     block_test%ghost_scalar( &
          mult_scalar*new_start+1 : &
          mult_scalar*(new_start+PATCH_SIZE**2)) = &
          sol(v_scalar,k_test)%data(d)%elts( &
          mult_scalar*old_start+1 : &
          mult_scalar*(old_start+PATCH_SIZE**2))

     block_test%ghost_vector( &
          mult_vector*new_start+1 : &
          mult_vector*(new_start+PATCH_SIZE**2)) = &
          sol(v_vector,k_test)%data(d)%elts( &
          mult_vector*old_start+1 : &
          mult_vector*(old_start+PATCH_SIZE**2))

  end do

  !
  ! Verify temporary ghost copies.
  !
  do ghost_id = 1, n_ghost

     old_start = &
          block_test%ghost_storage(ghost_id)%elts_start

     new_start = &
          block_test%ghost_storage(ghost_id)%local_start

     if (maxval(abs( &
          block_test%ghost_node( &
          new_start+1 : new_start+PATCH_SIZE**2)%x - &
          grid(d)%node%elts( &
          old_start+1 : old_start+PATCH_SIZE**2)%x)) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: ghost coordinate mismatch"

     end if

     if (maxval(abs( &
          block_test%ghost_scalar( &
          mult_scalar*new_start+1 : &
          mult_scalar*(new_start+PATCH_SIZE**2)) - &
          sol(v_scalar,k_test)%data(d)%elts( &
          mult_scalar*old_start+1 : &
          mult_scalar*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: ghost scalar mismatch"

     end if

     if (maxval(abs( &
          block_test%ghost_vector( &
          mult_vector*new_start+1 : &
          mult_vector*(new_start+PATCH_SIZE**2)) - &
          sol(v_vector,k_test)%data(d)%elts( &
          mult_vector*old_start+1 : &
          mult_vector*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

        error stop &
             "test_subtree_extraction: ghost vector mismatch"

     end if

  end do

  !
  ! ===============================================================
  ! Renumber neighbours and convert all external links to local
  ! block_bdry references.
  ! ===============================================================
  !
  call renumber_subtree_neigh_Domain( &
       grid(d), patch_copy, old_to_new)

  do ib = 1, size(block_test%block_bdry)

     patch_copy(block_test%block_bdry(ib)%patch+1)%neigh( &
          block_test%block_bdry(ib)%side) = -ib

  end do

  !
  ! Verify local neighbour representation.
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (patch_copy(old_to_new(p_old)+1)%neigh(c) /= &
                old_to_new(p_ngb_old)) then

              error stop &
                   "test_subtree_extraction: internal neighbour mismatch"

           end if

        case (NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT)

           ib = -patch_copy(old_to_new(p_old)+1)%neigh(c)

           if (ib < 1 .or. &
                ib > size(block_test%block_bdry)) then

              error stop &
                   "test_subtree_extraction: invalid boundary reference"

           end if

           if (block_test%block_bdry(ib)%patch /= &
                old_to_new(p_old) .or. &
                block_test%block_bdry(ib)%side /= c) then

              error stop &
                   "test_subtree_extraction: incorrect boundary reference"

           end if

        case (NGB_OTHER)

           if (patch_copy(old_to_new(p_old)+1)%neigh(c) /= 0) then
              error stop &
                   "test_subtree_extraction: zero neighbour changed"
           end if

        case default

           error stop &
                "test_subtree_extraction: unexpected neighbour class"

        end select

     end do

  end do

  !
  ! ===============================================================
  ! Package extracted interior arrays.
  ! ===============================================================
  !
  block_test%id          = block_catalog(b_test)%id
  block_test%root_domain = block_catalog(b_test)%root_domain
  block_test%root_patch  = block_catalog(b_test)%root_patch
  block_test%level       = block_catalog(b_test)%level

  call move_alloc(patch_copy,  block_test%patch)
  call move_alloc(node_copy,   block_test%node)
  call move_alloc(scalar_copy, block_test%scalar)
  call move_alloc(vector_copy, block_test%vector)

  !
  ! ===============================================================
  ! Build explicit compact stencil addresses.
  !
  ! NGB_BLOCK addresses use the unified compact ghost catalogue and
  ! the orientation-adjusted effective address from comp_offs3.
  ! ===============================================================
  !
  allocate(block_test%stencil( &
       N_BDRY,size(block_test%patch)))

  block_test%stencil%storage = STORE_NONE
  block_test%stencil%id      = -1
  block_test%stencil%offset  = 0

  do i = 1, size(block_test%stencil,2)
     do c = 1, size(block_test%stencil,1)
        block_test%stencil(c,i)%dims = 0
     end do
  end do

  n_stencil_built      = 0
  n_stencil_patch      = 0
  n_stencil_bdry       = 0
  n_stencil_ghost      = 0
  n_stencil_block      = 0
  n_block_value_checked = 0
  n_stencil_unresolved = 0

  n_unresolved_patch = 0
  n_unresolved_bdry  = 0
  n_unresolved_none  = 0

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     call comp_offs3( &
          grid(d), p_old, offs_src, dims_src)

     do c = 1, N_BDRY

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_OTHER)

           cycle

        case (NGB_INTERNAL, NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT)

           idx_src = grid(d)%patch%elts(p_old+1)%elts_start + &
                offs_src(c)

        case default

           error stop &
                "test_subtree_extraction: unexpected stencil class"

        end select

        target_patch  = -1
        target_bdry   = -1
        ghost_id      = -1
        target_offset = 0

        !
        ! Search extracted regular patches.
        !
        do i = 0, grid(d)%patch%length-1

           if (old_to_new(i) < 0) cycle

           old_start = grid(d)%patch%elts(i+1)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + PATCH_SIZE**2) then

              target_patch  = old_to_new(i)
              target_offset = idx_src - old_start
              exit

           end if

        end do

        !
        ! Search complete compact boundary storage.
        !
        if (target_patch < 0) then

           do is = 1, size(block_test%bdry_storage)

              old_start = block_test%bdry_storage(is)%elts_start

              if (idx_src >= old_start .and. &
                   idx_src < old_start + &
                   block_test%bdry_storage(is)%n_node) then

                 target_bdry   = is
                 target_offset = idx_src - old_start
                 exit

              end if

           end do

        end if

        !
        ! Search compact ghost-patch storage.
        !
        if (target_patch < 0 .and. target_bdry < 0) then

           do i = 1, size(block_test%ghost_storage)

              old_start = &
                   block_test%ghost_storage(i)%elts_start

              if (idx_src >= old_start .and. &
                   idx_src < old_start + &
                   block_test%ghost_storage(i)%n_node) then

                 ghost_id      = i
                 target_offset = idx_src - old_start
                 exit

              end if

           end do

        end if

        !
        ! For NGB_BLOCK, use the nominal neighbour ghost assigned to
        ! the explicit block-boundary record. The comp_offs3 address
        ! is an orientation-adjusted base and may lie outside the
        ! nominal ghost allocation. Therefore target_offset is a
        ! signed base displacement, not an immediately dereferenceable
        ! node index.
        !
        if (block_test%neigh_class( &
             c,old_to_new(p_old)+1) == NGB_BLOCK) then

           ib = -block_test%patch( &
                old_to_new(p_old)+1)%neigh(c)

           if (ib < 1 .or. &
                ib > size(block_test%block_bdry)) then

              error stop &
                   "test_subtree_extraction: invalid NGB_BLOCK record"

           end if

           if (block_test%block_bdry(ib)%class /= NGB_BLOCK) then
              error stop &
                   "test_subtree_extraction: incorrect NGB_BLOCK record"
           end if

           ghost_id = block_test%block_bdry(ib)%ghost_id

           if (ghost_id < 1 .or. &
                ghost_id > size(block_test%ghost_storage)) then

              error stop &
                   "test_subtree_extraction: invalid NGB_BLOCK ghost ID"

           end if

           if (block_test%ghost_storage(ghost_id)%source_patch /= &
                block_test%block_bdry(ib)%neigh_patch) then

              error stop &
                   "test_subtree_extraction: incorrect NGB_BLOCK ghost"

           end if

           target_patch = -1
           target_bdry  = -1

           old_start = &
                block_test%ghost_storage(ghost_id)%elts_start

           target_offset = idx_src - old_start

           if (target_offset > &
                block_test%ghost_storage(ghost_id)%n_node-1 .or. &
                target_offset + PATCH_SIZE**2-1 < 0) then

              error stop &
                   "test_subtree_extraction: empty NGB_BLOCK mapping"

           end if

        end if

        !
        ! Resolved interior target.
        !
        if (target_patch >= 0) then

           block_test%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_PATCH

           block_test%stencil( &
                c,old_to_new(p_old)+1)%id = target_patch

           block_test%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_test%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_patch = n_stencil_patch + 1
           n_stencil_built = n_stencil_built + 1

        !
        ! Resolved compact boundary target.
        !
        else if (target_bdry >= 0) then

           block_test%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_BDRY

           block_test%stencil( &
                c,old_to_new(p_old)+1)%id = target_bdry

           block_test%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_test%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_bdry  = n_stencil_bdry + 1
           n_stencil_built = n_stencil_built + 1

        !
        ! Resolved compact ghost-patch target.
        !
        else if (ghost_id >= 0) then

           block_test%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_GHOST

           block_test%stencil( &
                c,old_to_new(p_old)+1)%id = ghost_id

           block_test%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_test%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_ghost = n_stencil_ghost + 1
           n_stencil_built = n_stencil_built + 1

           if (block_test%neigh_class( &
                c,old_to_new(p_old)+1) == NGB_BLOCK) then

              n_stencil_block = n_stencil_block + 1

           end if

        !
        ! Not yet represented by the compact block.
        !
        else

           n_stencil_unresolved = n_stencil_unresolved + 1

           unresolved_patch = -1
           unresolved_bdry  = -1

           !
           ! Boundary closure should have captured every nominal
           ! boundary owner.
           !
           do is = 0, grid(d)%bdry_patch%length-1

              old_start = &
                   grid(d)%bdry_patch%elts(is+1)%elts_start

              if (idx_src >= old_start .and. &
                   idx_src < old_start + &
                   BDRY_THICKNESS*PATCH_SIZE) then

                 unresolved_bdry = is
                 exit

              end if

           end do

           !
           ! Otherwise search all source regular patches.
           !
           if (unresolved_bdry < 0) then

              do i = 0, grid(d)%patch%length-1

                 old_start = &
                      grid(d)%patch%elts(i+1)%elts_start

                 if (idx_src >= old_start .and. &
                      idx_src < old_start + PATCH_SIZE**2) then

                    unresolved_patch = i
                    exit

                 end if

              end do

           end if

           if (unresolved_bdry >= 0) then

              n_unresolved_bdry = n_unresolved_bdry + 1

           else if (unresolved_patch >= 0) then

              n_unresolved_patch = n_unresolved_patch + 1

           else

              n_unresolved_none = n_unresolved_none + 1

           end if

           cycle

        end if

        !
        ! Verify scalar value through explicit compact addressing.
        !
        select case (block_test%stencil( &
             c,old_to_new(p_old)+1)%storage)

        case (STORE_PATCH)

           target_patch = block_test%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_test%stencil( &
                c,old_to_new(p_old)+1)%offset

           val_src = sol(v_scalar,k_test)%data(d)%elts( &
                mult_scalar*idx_src + 1)

           val_blk = block_test%scalar( &
                mult_scalar * &
                (block_test%patch(target_patch+1)%elts_start + &
                target_offset) + 1)

        case (STORE_BDRY)

           target_bdry = block_test%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_test%stencil( &
                c,old_to_new(p_old)+1)%offset

           val_src = sol(v_scalar,k_test)%data(d)%elts( &
                mult_scalar*idx_src + 1)

           val_blk = block_test%bdry_scalar( &
                mult_scalar * &
                (block_test%bdry_storage(target_bdry)%local_start + &
                target_offset) + 1)

        case (STORE_GHOST)

           ghost_id = block_test%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_test%stencil( &
                c,old_to_new(p_old)+1)%offset

           if (block_test%neigh_class( &
                c,old_to_new(p_old)+1) == NGB_BLOCK) then

              n_mapped = 0

              do q = 0, PATCH_SIZE**2-1

                 ghost_offset = target_offset + q

                 if (ghost_offset < 0 .or. &
                      ghost_offset >= &
                      block_test%ghost_storage(ghost_id)%n_node) cycle

                 val_src = sol(v_scalar,k_test)%data(d)%elts( &
                      mult_scalar*(idx_src+q) + 1)

                 val_blk = block_test%ghost_scalar( &
                      mult_scalar * &
                      (block_test%ghost_storage( &
                      ghost_id)%local_start + ghost_offset) + 1)

                 if (abs(val_blk-val_src) > 0.0_dp) then
                    error stop &
                         "test_subtree_extraction: NGB_BLOCK value mismatch"
                 end if

                 n_mapped = n_mapped + 1

              end do

              if (n_mapped < 1) then
                 error stop &
                      "test_subtree_extraction: empty NGB_BLOCK validation"
              end if

              n_block_value_checked = &
                   n_block_value_checked + n_mapped

              ! Suppress the single-address comparison below; every
              ! physically mapped value has already been checked.
              val_src = 0.0_dp
              val_blk = 0.0_dp

           else

              val_src = sol(v_scalar,k_test)%data(d)%elts( &
                   mult_scalar*idx_src + 1)

              val_blk = block_test%ghost_scalar( &
                   mult_scalar * &
                   (block_test%ghost_storage( &
                   ghost_id)%local_start + target_offset) + 1)

           end if

        case default

           error stop &
                "test_subtree_extraction: invalid stencil storage"

        end select

        if (abs(val_blk-val_src) > 0.0_dp) then
           error stop &
                "test_subtree_extraction: explicit scalar stencil mismatch"
        end if

     end do

  end do

  !
  ! ===============================================================
  ! Final structural checks.
  ! ===============================================================
  !
  if (n_stencil_built + n_stencil_unresolved /= &
       n_ngb_internal + n_ngb_block + &
       n_ngb_domain + n_ngb_adapt) then

     error stop &
          "test_subtree_extraction: stencil count mismatch"

  end if

  if (n_stencil_patch + n_stencil_bdry + n_stencil_ghost /= &
       n_stencil_built) then

     error stop &
          "test_subtree_extraction: stencil class count mismatch"

  end if

  if (n_stencil_block /= n_ngb_block) then
     error stop &
          "test_subtree_extraction: NGB_BLOCK stencil count mismatch"
  end if

  if (n_block_value_checked < n_stencil_block) then
     error stop &
          "test_subtree_extraction: incomplete NGB_BLOCK value check"
  end if

  if (n_unresolved_patch + n_unresolved_bdry + &
       n_unresolved_none /= n_stencil_unresolved) then

     error stop &
          "test_subtree_extraction: unresolved count mismatch"

  end if

  !
  ! The boundary closure must now be complete.
  !
  if (n_unresolved_bdry /= 0) then
     error stop &
          "test_subtree_extraction: boundary-storage closure incomplete"
  end if

  !
  ! ===============================================================
  ! Diagnostic output.
  ! ===============================================================
  !
  if (verbose) then

  write(6,'(/,a,i0,a,i0)') &
       "Subtree extraction test: domain ", &
       block_catalog(b_test)%root_domain, &
       ", root patch ", p_root

  write(6,'(a,i0)') &
       "  block ID = ", block_catalog(b_test)%id

  write(6,'(a,i0)') &
       "  block weight = ", block_catalog(b_test)%weight

  write(6,'(a,i0)') &
       "  original subtree patches = ", n_old

  write(6,'(a,i0)') &
       "  extracted subtree patches = ", n_new

  write(6,'(a,i0)') &
       "  original leaves = ", n_leaf_old

  write(6,'(a,i0)') &
       "  extracted leaves = ", n_leaf_new

  write(6,'(a,i0)') &
       "  original subtree depth = ", depth_old

  write(6,'(a,i0)') &
       "  extracted subtree depth = ", depth_new

  write(6,'(a,i0)') &
       "  compact node storage size = ", n_node_storage

  write(6,'(a)') &
       "  interior coordinate/field copy checks passed"

  write(6,'(/,a)') &
       "  Block neighbour classification:"

  write(6,'(a,i0)') &
       "    internal neighbour links = ", n_ngb_internal

  write(6,'(a,i0)') &
       "    new inter-block boundary links = ", n_ngb_block

  write(6,'(a,i0)') &
       "    existing domain boundary links = ", n_ngb_domain

  write(6,'(a,i0)') &
       "    adaptive boundary links = ", n_ngb_adapt

  write(6,'(a,i0)') &
       "    other neighbour links = ", n_ngb_other

  write(6,'(/,a)') &
       "  Inter-block source mapping:"

  write(6,'(a,i0)') &
       "    mapped inter-block links = ", n_ngb_block

  write(6,'(a,i0)') &
       "    source owner matches block = ", &
       n_block_source_local

  write(6,'(a,i0)') &
       "    source owner differs       = ", &
       n_block_source_remote

  do ib = 1, size(block_test%block_bdry)

     if (block_test%block_bdry(ib)%class /= NGB_BLOCK) cycle

     write(6,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "    link ", ib, &
          ": source patch = ", &
          block_test%block_bdry(ib)%neigh_patch, &
          ", block index = ", &
          block_test%block_bdry(ib)%source_block, &
          ", block ID = ", &
          block_test%block_bdry(ib)%source_block_id, &
          ", owner = ", &
          block_test%block_bdry(ib)%source_owner, &
          ", ghost ID = ", &
          block_test%block_bdry(ib)%ghost_id

  end do

  write(6,'(a)') &
       "  unique inter-block source mapping check passed"

  write(6,'(/,a,i0)') &
       "  total local boundary records = ", &
       size(block_test%block_bdry)

  write(6,'(a,i0)') &
       "  directly referenced boundary patches = ", &
       n_bdry_direct

  write(6,'(a,i0)') &
       "  stencil-closure boundary patches = ", &
       n_bdry_closure

  write(6,'(a,i0)') &
       "  total compact boundary patches = ", &
       size(block_test%bdry_storage)

  write(6,'(a,i0)') &
       "  summed directly referenced boundary storage = ", &
       n_bdry_node_total

  write(6,'(a,i0)') &
       "  compact boundary node storage = ", &
       n_bdry_node_unique

  write(6,'(a,i0)') &
       "  maximum single boundary storage = ", &
       n_bdry_node_max

  write(6,'(a)') &
       "  boundary coordinate/scalar/vector copy checks passed"

  write(6,'(/,a,i0)') &
       "  compact ghost source patches = ", &
       size(block_test%ghost_storage)

  write(6,'(a,i0)') &
       "  compact ghost node storage = ", &
       size(block_test%ghost_node)

  do ghost_id = 1, size(block_test%ghost_storage)

     write(6,'(a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "    ghost ", ghost_id, &
          ": source patch = ", &
          block_test%ghost_storage(ghost_id)%source_patch, &
          ", block index = ", &
          block_test%ghost_storage(ghost_id)%source_block, &
          ", block ID = ", &
          block_test%ghost_storage(ghost_id)%source_block_id, &
          ", owner = ", &
          block_test%ghost_storage(ghost_id)%source_owner

  end do

  write(6,'(a)') &
       "  unified ghost catalogue and copy checks passed"

  write(6,'(/,a)') &
       "  Explicit compact stencil addressing:"

  write(6,'(a,i0)') &
       "    stencil addresses built       = ", &
       n_stencil_built

  write(6,'(a,i0)') &
       "    interior-patch targets        = ", &
       n_stencil_patch

  write(6,'(a,i0)') &
       "    boundary-storage targets      = ", &
       n_stencil_bdry

  write(6,'(a,i0)') &
       "    ghost-patch targets           = ", &
       n_stencil_ghost

  write(6,'(a,i0)') &
       "      nominal inter-block targets = ", &
       n_stencil_block

  write(6,'(a,i0)') &
       "      inter-block values checked  = ", &
       n_block_value_checked

  write(6,'(a,i0)') &
       "    unresolved targets            = ", &
       n_stencil_unresolved

  write(6,'(a,i0)') &
       "      outside-block patch targets = ", &
       n_unresolved_patch

  write(6,'(a,i0)') &
       "      uncopied boundary targets   = ", &
       n_unresolved_bdry

  write(6,'(a,i0)') &
       "      no nominal source owner     = ", &
       n_unresolved_none

  write(6,'(a)') &
       "  boundary-storage closure check passed"

  write(6,'(a)') &
       "  scalar explicit stencil addressing check passed"

  write(6,'(a,/)') &
       "  patch topology and storage layout checks passed"

  end if

  n_patch_out   = n_new
  n_bdry_out    = n_bdry_required
  n_ghost_out   = n_ghost
  n_stencil_out = n_stencil_built
  n_remote_out  = n_block_source_remote
  n_value_out   = n_block_value_checked

  deallocate(ghost_patch)
  deallocate(bdry_required)
  deallocate(bdry_closure)
  deallocate(old_to_new)
  deallocate(old_elts_start)

contains

  recursive logical function patch_in_subtree (p, p_target) &
       result(found_patch)

    implicit none

    integer, intent(in) :: p
    integer, intent(in) :: p_target

    integer :: c
    integer :: p_chd

    if (p < 0 .or. p >= grid(d)%patch%length) then
       error stop "patch_in_subtree: invalid source patch"
    end if

    found_patch = p == p_target

    if (found_patch) return

    do c = 1, N_CHDRN

       p_chd = grid(d)%patch%elts(p+1)%children(c)

       if (p_chd > 0) then

          if (patch_in_subtree(p_chd,p_target)) then
             found_patch = .true.
             return
          end if

       end if

    end do

  end function patch_in_subtree

  recursive integer function copied_depth (p) result(depth)

    implicit none

    integer, intent(in) :: p

    integer :: c
    integer :: p_chd

    depth = 0

    do c = 1, N_CHDRN

       p_chd = patch_copy(p+1)%children(c)

       if (p_chd > 0) then
          depth = max(depth,1+copied_depth(p_chd))
       end if

    end do

  end function copied_depth

end subroutine test_one_subtree_extraction



  subroutine test_parallel_block_split
  ! Diagnostic only: construct a complete prospective parallel-block
  ! decomposition of all locally owned geometric root domains, build
  ! a global block catalogue, and test prospective load balancing.
  !
  ! Rule:
  !   - if root patch 1 has all four children, use those four child
  !     subtrees as candidate blocks;
  !   - otherwise keep the whole root domain as one candidate block.
  !
  ! This routine does NOT modify the active decomposition.

  implicit none

  integer, parameter :: N_BLOCK_DATA = 6

  integer :: b, c, d, i, k, p, p_chd, r
  integer :: n_block_local, n_block_before, n_block_total
  integer :: n_data_local, n_data_total
  integer :: n_assigned, n_changed
  integer :: ierr

  integer, allocatable :: block_count(:)
  integer, allocatable :: recv_count(:)
  integer, allocatable :: recv_disp(:)
  integer, allocatable :: send_data(:)
  integer, allocatable :: recv_data(:)
  integer, allocatable :: domain_count(:)

  integer, allocatable :: block_owner_new(:)
  integer, allocatable :: load_current(:)
  integer, allocatable :: load_proposed(:)

  real(dp) :: balanced_weight
  real(dp) :: imbalance_goal

  real(dp), parameter :: init_goal = 0.05_dp
  real(dp), parameter :: incr_goal = 1.20_dp

  type(Parallel_Block), allocatable :: block_loc(:)

  !
  ! Count local candidate blocks.
  !
  n_block_local = 0

  do d = 1, size(grid)

     p = 1

     if (all(grid(d)%patch%elts(p+1)%children > 0)) then
        n_block_local = n_block_local + N_CHDRN
     else
        n_block_local = n_block_local + 1
     end if

  end do

  allocate(block_loc(n_block_local))

  !
  ! Construct local candidate block descriptors.
  !
  i = 0

  do d = 1, size(grid)

     p = 1

     if (all(grid(d)%patch%elts(p+1)%children > 0)) then

        do c = 1, N_CHDRN

           p_chd = grid(d)%patch%elts(p+1)%children(c)

           i = i + 1

           block_loc(i)%root_domain = glo_id(rank+1,d)
           block_loc(i)%root_patch  = p_chd
           block_loc(i)%level       = grid(d)%patch%elts(p_chd+1)%level
           block_loc(i)%owner       = rank
           block_loc(i)%weight      = &
                subtree_weight_Domain(grid(d), p_chd)

        end do

     else

        i = i + 1

        block_loc(i)%root_domain = glo_id(rank+1,d)
        block_loc(i)%root_patch  = p
        block_loc(i)%level       = grid(d)%patch%elts(p+1)%level
        block_loc(i)%owner       = rank
        block_loc(i)%weight      = &
             subtree_weight_Domain(grid(d), p)

     end if

  end do

  if (i /= n_block_local) then
     error stop &
          "test_parallel_block_split: local block count mismatch"
  end if

  !
  ! Determine the number of candidate blocks preceding those on
  ! this rank. This gives globally unique contiguous block IDs.
  !
  n_block_before = 0

  call MPI_Exscan( &
       n_block_local, n_block_before, &
       1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "test_parallel_block_split: MPI_Exscan failed"
  end if

  if (rank == 0) n_block_before = 0

  !
  ! Determine the total number of candidate blocks.
  !
  call MPI_Allreduce( &
       n_block_local, n_block_total, &
       1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "test_parallel_block_split: MPI_Allreduce failed"
  end if

  !
  ! Assign global candidate-block IDs.
  !
  do i = 1, n_block_local
     block_loc(i)%id = n_block_before + i - 1
  end do

  !
  ! Determine how many candidate blocks are contributed by every rank.
  !
  allocate(block_count(n_process))

  call MPI_Allgather( &
       n_block_local, 1, MPI_INTEGER, &
       block_count, 1, MPI_INTEGER, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "test_parallel_block_split: MPI_Allgather failed"
  end if

  if (sum(block_count) /= n_block_total) then
     error stop &
          "test_parallel_block_split: inconsistent global block count"
  end if

  !
  ! Pack each local candidate block as six integers:
  !
  !   id, root_domain, root_patch, level, owner, weight
  !
  n_data_local = N_BLOCK_DATA * n_block_local

  allocate(send_data(n_data_local))

  k = 0

  do i = 1, n_block_local

     send_data(k+1:k+N_BLOCK_DATA) = [ &
          block_loc(i)%id, &
          block_loc(i)%root_domain, &
          block_loc(i)%root_patch, &
          block_loc(i)%level, &
          block_loc(i)%owner, &
          block_loc(i)%weight ]

     k = k + N_BLOCK_DATA

  end do

  if (k /= n_data_local) then
     error stop &
          "test_parallel_block_split: local packing count mismatch"
  end if

  !
  ! Construct Allgatherv receive counts and displacements.
  !
  allocate(recv_count(n_process))
  allocate(recv_disp(n_process))

  recv_count = N_BLOCK_DATA * block_count

  recv_disp(1) = 0

  do i = 2, n_process
     recv_disp(i) = recv_disp(i-1) + recv_count(i-1)
  end do

  n_data_total = N_BLOCK_DATA * n_block_total

  allocate(recv_data(n_data_total))

  !
  ! Every rank receives the complete packed candidate-block catalogue.
  !
  call MPI_Allgatherv( &
       send_data, n_data_local, MPI_INTEGER, &
       recv_data, recv_count, recv_disp, MPI_INTEGER, &
       comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "test_parallel_block_split: MPI_Allgatherv failed"
  end if

  !
  ! Construct the global Parallel_Block catalogue on every rank.
  !
  if (allocated(block_catalog)) deallocate(block_catalog)

  allocate(block_catalog(n_block_total))

  do i = 1, n_block_total

     k = N_BLOCK_DATA * (i-1)

     block_catalog(i)%id          = recv_data(k+1)
     block_catalog(i)%root_domain = recv_data(k+2)
     block_catalog(i)%root_patch  = recv_data(k+3)
     block_catalog(i)%level       = recv_data(k+4)
     block_catalog(i)%owner       = recv_data(k+5)
     block_catalog(i)%weight      = recv_data(k+6)

  end do

  !
  ! Validate global block IDs and owners.
  !
  do i = 1, n_block_total

     if (block_catalog(i)%id /= i-1) then
        error stop &
             "test_parallel_block_split: invalid global block ID ordering"
     end if

     if (block_catalog(i)%owner < 0 .or. &
          block_catalog(i)%owner >= n_process) then
        error stop &
             "test_parallel_block_split: invalid block owner"
     end if

  end do

  !
  ! Verify that every geometric root domain occurs in the catalogue.
  ! With the current one-level split rule, each domain should occur
  ! either once or four times.
  !
  allocate(domain_count(N_GLO_DOMAIN))
  domain_count = 0

  do i = 1, n_block_total

     if (block_catalog(i)%root_domain < 0 .or. &
          block_catalog(i)%root_domain >= N_GLO_DOMAIN) then
        error stop &
             "test_parallel_block_split: invalid root domain"
     end if

     domain_count(block_catalog(i)%root_domain+1) = &
          domain_count(block_catalog(i)%root_domain+1) + 1

  end do

  if (any(domain_count == 0)) then
     error stop &
          "test_parallel_block_split: one or more root domains missing"
  end if

  if (any(domain_count /= 1 .and. domain_count /= N_CHDRN)) then
     error stop &
          "test_parallel_block_split: invalid number of blocks per root domain"
  end if

  !
  ! ---------------------------------------------------------------
  ! Test prospective load balancing of candidate parallel blocks.
  ! ---------------------------------------------------------------
  !
  ! This does not change actual ownership.
  !
  allocate(block_owner_new(n_block_total))
  allocate(load_current(n_process))
  allocate(load_proposed(n_process))

  block_owner_new = -1
  load_current     = 0
  load_proposed    = 0

  !
  ! Current load implied by the existing domain decomposition.
  !
  do b = 1, n_block_total

     r = block_catalog(b)%owner + 1

     if (r < 1 .or. r > n_process) then
        error stop &
             "test_parallel_block_split: invalid current block owner"
     end if

     load_current(r) = load_current(r) + block_catalog(b)%weight

  end do

  !
  ! Ideal average load per MPI rank.
  !
  balanced_weight = &
       real(sum(block_catalog%weight), kind=dp) / &
       real(n_process, kind=dp)

  !
  ! Use the same next-fit strategy as distribute_grid(), but apply it
  ! only to the prospective candidate blocks.
  !
  n_assigned     = 0
  imbalance_goal = init_goal

  do while (n_assigned < n_block_total)

     n_assigned      = 0
     load_proposed   = 0
     block_owner_new = -1

     do r = 1, n_process

        do while ( &
             real(load_proposed(r),kind=dp) < balanced_weight .and. &
             n_block_total - n_assigned > n_process - r)

           block_owner_new(n_assigned+1) = r - 1

           load_proposed(r) = load_proposed(r) + &
                block_catalog(n_assigned+1)%weight

           n_assigned = n_assigned + 1

        end do

        !
        ! If the final block makes this rank too heavy,
        ! move that block to the next rank.
        !
        if (n_assigned > 0) then

           if (real(load_proposed(r),kind=dp) > &
                balanced_weight * (1.0_dp + imbalance_goal)) then

              load_proposed(r) = load_proposed(r) - &
                   block_catalog(n_assigned)%weight

              block_owner_new(n_assigned) = -1

              n_assigned = n_assigned - 1

           end if

        end if

     end do

     !
     ! If all blocks could not be assigned, relax the allowed imbalance.
     !
     if (n_assigned < n_block_total) then
        imbalance_goal = imbalance_goal * incr_goal
     end if

  end do

  if (any(block_owner_new < 0)) then
     error stop &
          "test_parallel_block_split: one or more candidate blocks unassigned"
  end if

  !
  ! Independently reconstruct proposed rank loads from the owner map.
  !
  load_proposed = 0

  do b = 1, n_block_total

     r = block_owner_new(b) + 1

     if (r < 1 .or. r > n_process) then
        error stop &
             "test_parallel_block_split: invalid proposed block owner"
     end if

     load_proposed(r) = load_proposed(r) + block_catalog(b)%weight

  end do

  n_changed = count(block_owner_new /= block_catalog%owner)

  !
  ! Commit the prospective balanced ownership to the global block
  ! catalogue. The current source-domain owner remains available through
  ! owner(root_domain+1).
  !
  block_catalog%owner = block_owner_new

  !
  ! Validate the committed target owners.
  !
  if (any(block_catalog%owner < 0) .or. &
       any(block_catalog%owner >= n_process)) then

     error stop &
          "test_parallel_block_split: invalid committed block owner"

  end if
  
  !
  ! Print compact diagnostic summary.
  !
  if (rank == 0) then

     write(6,'(/,a,i0)') &
          "Total candidate parallel blocks = ", n_block_total

     write(6,'(a,i0)') &
          "Unsplit root domains = ", count(domain_count == 1)

     write(6,'(a,i0)') &
          "Split root domains   = ", count(domain_count == N_CHDRN)

     write(6,'(/,a)') &
          "Current candidate-block load:"

     write(6,'(a,i0)') &
          "  minimum = ", minval(load_current)

     write(6,'(a,f10.2)') &
          "  average = ", balanced_weight

     write(6,'(a,i0)') &
          "  maximum = ", maxval(load_current)

     write(6,'(a,f10.4)') &
          "  max/avg = ", &
          real(maxval(load_current),kind=dp) / balanced_weight

     write(6,'(/,a)') &
          "Prospective balanced block load:"

     write(6,'(a,i0)') &
          "  minimum = ", minval(load_proposed)

     write(6,'(a,f10.2)') &
          "  average = ", balanced_weight

     write(6,'(a,i0)') &
          "  maximum = ", maxval(load_proposed)

     write(6,'(a,f10.4)') &
          "  max/avg = ", &
          real(maxval(load_proposed),kind=dp) / balanced_weight

     write(6,'(a,i0,a,i0,/)') &
          "Blocks changing owner = ", n_changed, &
          " / ", n_block_total

  end if

  !
  ! Cleanup temporary arrays.
  !
  deallocate(block_loc)
  deallocate(block_count)
  deallocate(recv_count)
  deallocate(recv_disp)
  deallocate(send_data)
  deallocate(recv_data)
  deallocate(domain_count)

  deallocate(block_owner_new)
  deallocate(load_current)
  deallocate(load_proposed)

end subroutine test_parallel_block_split
 

  subroutine deallocate_structures
    ! Deallocate all dynamic arrays and structures for clean restart
    implicit none

    integer :: d, i, k, l, v, r
    integer :: istat
    character(len=256) :: emsg

    if (.not. allocated (grid)) return
    if (size(grid) <= 0) return

    if (allocated (Area_avg)) deallocate (Area_avg)
    if (allocated (C_visc))   deallocate (C_visc)
    if (allocated (dx_avg))   deallocate (dx_avg)

    ! Deallocate init_RK_mem allocations
    do k = 1, zmax
       do d = 1, n_domain(rank+1)
          do v = 1, N_VARIABLE
             if (allocated (q1(v,k)%data(d)%elts))  deallocate (q1(v,k)%data(d)%elts)
             if (allocated (q2(v,k)%data(d)%elts))  deallocate (q2(v,k)%data(d)%elts)
             if (allocated (q3(v,k)%data(d)%elts))  deallocate (q3(v,k)%data(d)%elts)
             if (allocated (q4(v,k)%data(d)%elts))  deallocate (q4(v,k)%data(d)%elts)
             if (allocated (dq1(v,k)%data(d)%elts)) deallocate (dq1(v,k)%data(d)%elts)
          end do
       end do
       do v = 1, N_VARIABLE
          if (allocated (q1(v,k)%data))  deallocate (q1(v,k)%data)
          if (allocated (q2(v,k)%data))  deallocate (q2(v,k)%data)
          if (allocated (q3(v,k)%data))  deallocate (q3(v,k)%data)
          if (allocated (q4(v,k)%data))  deallocate (q4(v,k)%data)
          if (allocated (dq1(v,k)%data)) deallocate (dq1(v,k)%data)
       end do
    end do

    if (allocated (q1))  deallocate (q1)
    if (allocated (q2))  deallocate (q2)
    if (allocated (q3))  deallocate (q3)
    if (allocated (q4))  deallocate (q4)
    if (allocated (dq1)) deallocate (dq1)

    ! Deallocate grid structure elements
    do d = 1, size(grid)
       if (allocated (grid(d)%mask_n%elts))             deallocate (grid(d)%mask_n%elts)
       if (allocated (grid(d)%mask_e%elts))             deallocate (grid(d)%mask_e%elts)

       if (allocated(grid(d)%level%elts)) then
          deallocate(grid(d)%level%elts, stat=istat, errmsg=emsg)
          if (istat /= 0) then
             write(*,*) "rank", rank, "d", d, "dealloc level%elts failed:", trim(emsg)
             error stop
          end if
       end if

       if (allocated (grid(d)%level%elts))              deallocate (grid(d)%level%elts)

       if (allocated (grid(d)%R_F_wgt%elts))            deallocate (grid(d)%R_F_wgt%elts)
       if (allocated (grid(d)%I_u_wgt%elts))            deallocate (grid(d)%I_u_wgt%elts)

       if (allocated (grid(d)%overl_areas%elts))        deallocate (grid(d)%overl_areas%elts)
       if (allocated (grid(d)%triarea%elts))            deallocate (grid(d)%triarea%elts)
       if (allocated (grid(d)%len%elts))                deallocate (grid(d)%len%elts)
       if (allocated (grid(d)%pedlen%elts))             deallocate (grid(d)%pedlen%elts)
       if (allocated (grid(d)%areas%elts))              deallocate (grid(d)%areas%elts)
       if (allocated (grid(d)%midpt%elts))              deallocate (grid(d)%midpt%elts)
       if (allocated (grid(d)%ccentre%elts))            deallocate (grid(d)%ccentre%elts)

       if (allocated (grid(d)%surf_press%elts))         deallocate (grid(d)%surf_press%elts)
       if (allocated (grid(d)%press%elts))              deallocate (grid(d)%press%elts)
       if (allocated (grid(d)%geopot%elts))             deallocate (grid(d)%geopot%elts)
       if (allocated (grid(d)%u_zonal%elts))            deallocate (grid(d)%u_zonal%elts)
       if (allocated (grid(d)%v_merid%elts))            deallocate (grid(d)%v_merid%elts)
       if (allocated (grid(d)%press_lower%elts))        deallocate (grid(d)%press_lower%elts)
       if (allocated (grid(d)%geopot_lower%elts))       deallocate (grid(d)%geopot_lower%elts)
       if (allocated (grid(d)%vort%elts))               deallocate (grid(d)%vort%elts)
       if (allocated (grid(d)%qe%elts))                 deallocate (grid(d)%qe%elts)
       if (allocated (grid(d)%bernoulli%elts))          deallocate (grid(d)%bernoulli%elts)
       if (allocated (grid(d)%ke%elts))                 deallocate (grid(d)%ke%elts)
       if (allocated (grid(d)%divu%elts))               deallocate (grid(d)%divu%elts)
       if (allocated (grid(d)%coriolis%elts))           deallocate (grid(d)%coriolis%elts)

       if (allocated (grid(d)%node%elts))               deallocate (grid(d)%node%elts) 
       if (allocated (grid(d)%bdry_patch%elts))         deallocate (grid(d)%bdry_patch%elts) 
       if (allocated (grid(d)%patch%elts))              deallocate (grid(d)%patch%elts) 
       if (allocated (grid(d)%neigh_pa_over_pole%elts)) deallocate (grid(d)%neigh_pa_over_pole%elts)
       if (allocated (grid(d)%send_pa_all%elts))        deallocate (grid(d)%send_pa_all%elts)

       do i = 1, n_glo_block
          if (allocated (grid(d)%recv_pa(i)%elts))   deallocate (grid(d)%recv_pa(i)%elts)
          if (allocated (grid(d)%send_conn(i)%elts)) deallocate (grid(d)%send_conn(i)%elts)
          do k = AT_NODE, AT_EDGE
             if (allocated(grid(d)%pack(k,i)%elts)) deallocate (grid(d)%pack(k,i)%elts)
             if (allocated(grid(d)%unpk(k,i)%elts)) deallocate (grid(d)%unpk(k,i)%elts)
          end do
       end do

       do i = lbound(grid(d)%lev,1), ubound(grid(d)%lev,1)
          if (allocated (grid(d)%lev(i)%elts)) deallocate (grid(d)%lev(i)%elts)
       end do
       deallocate (grid(d)%lev)

       do l = min_level, max_level
          do r = 1, n_process
             if (allocated (grid(d)%src_patch(r,l)%elts)) deallocate (grid(d)%src_patch(r,l)%elts) 
          end do
       end do
       deallocate (grid(d)%src_patch)

       if (allocated(Laplacian_vector(S_DIVU)%data(d)%elts)) deallocate (Laplacian_vector(S_DIVU)%data(d)%elts)
       if (allocated(Laplacian_vector(S_ROTU)%data(d)%elts)) deallocate (Laplacian_vector(S_ROTU)%data(d)%elts)

       do v = scalars(1), scalars(2)
          if (allocated(horiz_flux(v)%data(d)%elts))       deallocate (horiz_flux(v)%data(d)%elts)
          if (allocated(Laplacian_scalar(v)%data(d)%elts)) deallocate (Laplacian_scalar(v)%data(d)%elts)
       end do

       do k = zmin, zmax
          if (allocated(penal_node(k)%data(d)%elts))  deallocate (penal_node(k)%data(d)%elts)
          if (allocated (penal_edge(k)%data(d)%elts)) deallocate (penal_edge(k)%data(d)%elts)
          if (allocated (exner_fun(k)%data(d)%elts))  deallocate (exner_fun(k)%data(d)%elts)
       end do
       if (allocated (exner_fun(zmax+1)%data(d)%elts)) deallocate (exner_fun(zmax+1)%data(d)%elts)

       do v = 1, N_VARIABLE
          do k = zmin, zmax
             if (allocated(sol(v,k)%data(d)%elts))       deallocate (sol(v,k)%data(d)%elts)
             if (allocated(sol_mean(v,k)%data(d)%elts))  deallocate (sol_mean(v,k)%data(d)%elts)
             if (allocated(wav_coeff(v,k)%data(d)%elts)) deallocate (wav_coeff(v,k)%data(d)%elts)
          end do
          do k = 1, zlevels
             if (allocated(trend(v,k)%data(d)%elts)) deallocate (trend(v,k)%data(d)%elts)
          end do
       end do

       if (vert_diffuse) then
          if (allocated (Kt(0)%data(d)%elts)) deallocate (Kt(0)%data(d)%elts)
          if (allocated (Kv(0)%data(d)%elts)) deallocate (Kv(0)%data(d)%elts)
          do k = 1, zlevels
             if (allocated (Kt(k)%data(d)%elts))      deallocate (Kt(k)%data(d)%elts)
             if (allocated (Kv(k)%data(d)%elts))      deallocate (Kv(k)%data(d)%elts)
             if (allocated (tke(k)%data(d)%elts))     deallocate (tke(k)%data(d)%elts)
             if (allocated (wav_tke(k)%data(d)%elts)) deallocate (wav_tke(k)%data(d)%elts)
          end do
       end if

       if (NCAR_topo) then
          do l = min_level, max_level
             if (allocated (topography_data(l,d)%node)) deallocate (topography_data(l,d)%node)
             if (allocated (topography_data(l,d)%elts)) deallocate (topography_data(l,d)%elts)
          end do
       end if
    end do

    deallocate (topography%data)
    if (sso) then
       do k = 1, 4
          if (allocated (sso_param(k)%data)) deallocate (sso_param(k)%data)
       end do
    end if
    if (NCAR_topo .and. allocated (topography_data)) deallocate (topography_data)

    if (allocated (Laplacian_vector(S_DIVU)%data))  deallocate (Laplacian_vector(S_DIVU)%data)
    if (allocated (Laplacian_vector(S_ROTU)%data))  deallocate (Laplacian_vector(S_ROTU)%data)

    do k = zmin, zmax
       if (allocated(penal_node(k)%data)) deallocate (penal_node(k)%data)
       if (allocated(penal_edge(k)%data)) deallocate (penal_edge(k)%data)
       if (allocated(exner_fun(k)%data))  deallocate (exner_fun(k)%data)
    end do
    if (allocated(exner_fun(zmax+1)%data)) deallocate (exner_fun(zmax+1)%data)

    do v = scalars(1), scalars(2)
       if (allocated (horiz_flux(v)%data))       deallocate (horiz_flux(v)%data)
       if (allocated (Laplacian_scalar(v)%data)) deallocate (Laplacian_scalar(v)%data)
    end do

    do v = 1, N_VARIABLE
       do k = zmin, zmax
          if (allocated (sol(v,k)%data))       deallocate (sol(v,k)%data)
          if (allocated (sol_mean(v,k)%data))  deallocate (sol_mean(v,k)%data)
          if (allocated (wav_coeff(v,k)%data)) deallocate (wav_coeff(v,k)%data)
       end do
       do k = 1, zlevels
          if (allocated (trend(v,k)%data)) deallocate (trend(v,k)%data)
       end do
    end do

    if (vert_diffuse) then
       if (allocated (Kv(0)%data)) deallocate (Kv(0)%data)
       if (allocated (Kt(0)%data)) deallocate (Kt(0)%data)
       do k = 1, zlevels
          if (allocated (Kv(k)%data))      deallocate (Kv(k)%data)
          if (allocated (Kt(k)%data))      deallocate (Kt(k)%data)
          if (allocated (tke(k)%data))     deallocate (tke(k)%data)
          if (allocated (wav_tke(k)%data)) deallocate (wav_tke(k)%data)
       end do
    end if

    if (allocated (grid))             deallocate (grid)
    if (allocated (n_node_old))       deallocate (n_node_old)
    if (allocated (n_patch_old))      deallocate (n_patch_old)
    if (allocated (edge_level_start)) deallocate (edge_level_start)
    if (allocated (node_level_start)) deallocate (node_level_start)
    if (allocated (n_active_edges))   deallocate (n_active_edges)
    if (allocated (n_active_nodes))   deallocate (n_active_nodes)
    if (allocated (a_vert))           deallocate (a_vert)
    if (allocated (b_vert))           deallocate (b_vert)
    if (allocated (a_vert_mass))      deallocate (a_vert_mass)
    if (allocated (b_vert_mass))      deallocate (b_vert_mass)
    if (allocated (threshold))        deallocate (threshold)
    if (allocated (threshold_def))    deallocate (threshold_def)
    if (allocated (sol))              deallocate (sol)
    if (allocated (sol_mean))         deallocate (sol_mean)
    if (allocated (trend))            deallocate (trend)
    if (allocated (wav_coeff))        deallocate (wav_coeff)
    if (allocated (exner_fun))        deallocate (exner_fun)
    if (allocated (horiz_flux))       deallocate (horiz_flux)
    if (allocated (Laplacian_scalar)) deallocate (Laplacian_scalar)
    if (allocated (Laplacian_vector)) deallocate (Laplacian_vector)
    if (allocated (lnorm))            deallocate (lnorm)
    if (allocated (penal_edge))       deallocate (penal_edge)
    if (allocated (penal_node))       deallocate (penal_node)
    if (allocated (glo_id))           deallocate (glo_id)
    if (allocated (ini_st))           deallocate (ini_st)
    if (allocated (recv_lengths))     deallocate (recv_lengths)
    if (allocated (recv_offsets))     deallocate (recv_offsets)
    if (allocated (req))              deallocate (req)
    if (allocated (send_lengths))     deallocate (send_lengths)
    if (allocated (send_offsets))     deallocate (send_offsets)

    if (vert_diffuse) then
       if (allocated (Kv))      deallocate (Kv)
       if (allocated (Kt))      deallocate (Kt)
       if (allocated (tke))     deallocate (tke)
       if (allocated (wav_tke)) deallocate (wav_tke)
    end if

    nullify (mass, dscalar, h_flux, velo, dvelo, bernoulli, divu, exner, ke, qe, scalar, temp, vort, wc_s, wc_u)
  end subroutine deallocate_structures
end module main_mod
