module main_mod

  use mpi_f08, only : MPI_Allgather, MPI_Allgatherv, MPI_Allreduce, MPI_Exscan, &
       MPI_INTEGER, MPI_SUCCESS, MPI_SUM
  
  use, intrinsic :: ieee_arithmetic

  use kind_mod,   only : dp
  use shared_mod, only : ADJZONE, AT_NODE, AT_EDGE, BDRY_THICKNESS, CP_EVERY, DATA_GRID, DAY, EDGE, MATH_PI, MULT, N_BDRY, &
       N_CHDRN, N_GLO_DOMAIN, N_VARIABLE, NCAR_topo, NONE, &
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
  use patch_mod,          only : Patch, PATCH_SIZE
  use refine_patch_mod,   only : add_second_level
  use remap_mod,          only : remap_vertical_coordinates
  use utils_mod,          only : hex_len, hex_pedlen, interp, nu_scale, porous_density, tri_perim
  use vert_diffusion_mod, only : vertical_diffusion
  use wavelet_mod,        only : forward_wavelet_transform, init_wavelets, inverse_scalar_transform, inverse_wavelet_transform

  use arch_mod, only : abort_run, barrier, block_catalog, comm, distribute_grid, glo_id, init_arch_mod, loc_id, &
       n_glo_block, n_process, Parallel_Block, rank

  use comm_mpi_mod, only : comm_nodes3_mpi, init_comm_mpi, recv_lengths, recv_offsets, req, send_lengths, send_offsets, &
       sum_int,  sync_max_int, sync_min_real, write_load_conn

  use domain_mod, only : Domain, bernoulli, divu, dscalar, grid, &
       dvelo, exner, exner_fun, h_flux, horiz_flux, ke, qe, scalar, mass, temp, velo, vort, wc_s, wc_u, &
       Kt, Kv, Laplacian_scalar, Laplacian_vector, penal_edge, penal_node, sso_param, &
       sol, sol_mean, tke, topography, topography_data, trend, wav_coeff, wav_tke, id_edge, idx, &
       subtree_weight_Domain, count_subtree_patches_Domain, extract_subtree_patches_Domain, subtree_depth_Domain, &
       compact_subtree_storage_Domain, copy_subtree_nodes_Domain, copy_subtree_field_Domain

  use init_mod, only : apply_initial_conditions, elliptic_solver, init_geometry, init_grid, &
       initialize_a_b_vert, initialize_thresholds, initialize_dt_viscosity, &
       precompute_geometry, set_level, set_thresholds, z_coords

  use time_integr_mod, only : dt_step, dt_step_split, init_RK_mem, q1, q2, q3, q4, dq1, &
       Euler, Euler_split, RK2_split, RK3, RK3_split, RK33_opt, RK34_opt, RK4, RK4_split, RK45_opt

  use coord_arithmetic_mod

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

  type Test_Block
     integer :: id          = -1
     integer :: root_domain = -1
     integer :: root_patch  = -1
     integer :: level       = -1

     type(Patch), allocatable :: patch(:)
     type(Coord), allocatable :: node(:)

     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)

     integer, allocatable :: neigh_class(:,:)
  end type Test_Block

  integer, parameter :: NGB_INTERNAL = 1
  integer, parameter :: NGB_BLOCK    = 2
  integer, parameter :: NGB_DOMAIN   = 3
  integer, parameter :: NGB_ADAPT    = 4
  integer, parameter :: NGB_OTHER    = 5
 
  integer, allocatable :: n_active_edges(:), n_active_nodes(:), node_level_start(:), edge_level_start(:)
  real(dp)             :: dt_new, initial_total_mass, time_mult
  real(dp)             :: dt_loc, min_mass_loc

  type(Initial_State), allocatable :: ini_st(:)

  
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
  ! Extract one candidate subtree on rank zero and verify its topology,
  ! compact storage layout, copied geometry and fields, packaging into
  ! a temporary block, and classification of its patch neighbours.

  implicit none

  integer :: b_test
  integer :: c, d, i
  integer :: p_root
  integer :: n_old, n_new
  integer :: n_leaf_old, n_leaf_new
  integer :: depth_old, depth_new
  integer :: n_node_storage
  integer :: n_patch_field
  integer :: p_old, p_chd_old, p_chd_new
  integer :: p_ngb_old
  integer :: old_start, new_start
  integer :: v_scalar, v_vector, k_test
  integer :: mult_scalar, mult_vector

  integer :: n_ngb_internal
  integer :: n_ngb_block
  integer :: n_ngb_domain
  integer :: n_ngb_adapt
  integer :: n_ngb_other

  integer, allocatable :: old_to_new(:)
  integer, allocatable :: old_elts_start(:)

  real(dp), allocatable :: scalar_copy(:)
  real(dp), allocatable :: vector_copy(:)

  type(Patch), allocatable :: patch_copy(:)
  type(Coord), allocatable :: node_copy(:)

  type(Test_Block) :: block_test

  logical :: found

  !
  ! Keep the extraction test entirely local to rank zero.
  !
  if (rank /= 0) return

  found  = .false.
  b_test = -1

  !
  ! Prefer a candidate block containing a subtree of depth at least 2.
  !
  do i = 1, size(block_catalog)

     if (block_catalog(i)%owner /= rank) cycle

     d = loc_id(block_catalog(i)%root_domain+1) + 1

     if (d < 1 .or. d > size(grid)) cycle

     p_root = block_catalog(i)%root_patch

     if (subtree_depth_Domain(grid(d),p_root) >= 2) then
        b_test = i
        found  = .true.
        exit
     end if

  end do

  !
  ! If no sufficiently deep subtree exists, use the first local
  ! nontrivial candidate.
  !
  if (.not. found) then

     do i = 1, size(block_catalog)

        if (block_catalog(i)%owner /= rank) cycle

        d = loc_id(block_catalog(i)%root_domain+1) + 1

        if (d < 1 .or. d > size(grid)) cycle

        p_root = block_catalog(i)%root_patch

        if (subtree_weight_Domain(grid(d),p_root) > PATCH_SIZE**2) then
           b_test = i
           found  = .true.
           exit
        end if

     end do

  end if

  !
  ! Final fallback: use the first local candidate.
  !
  if (.not. found) then

     do i = 1, size(block_catalog)

        if (block_catalog(i)%owner /= rank) cycle

        d = loc_id(block_catalog(i)%root_domain+1) + 1

        if (d < 1 .or. d > size(grid)) cycle

        p_root = block_catalog(i)%root_patch

        b_test = i
        found  = .true.
        exit

     end do

  end if

  if (.not. found) then
     error stop "test_subtree_extraction: no local candidate block"
  end if

  !
  ! Recover local source domain and root patch.
  !
  d = loc_id(block_catalog(b_test)%root_domain+1) + 1

  if (d < 1 .or. d > size(grid)) then
     error stop "test_subtree_extraction: invalid local domain"
  end if

  p_root = block_catalog(b_test)%root_patch

  depth_old = subtree_depth_Domain(grid(d), p_root)

  !
  ! ---------------------------------------------------------------
  ! Extract and renumber the patch tree.
  ! ---------------------------------------------------------------
  !
  call extract_subtree_patches_Domain( &
       grid(d), p_root, patch_copy, old_to_new)

  n_old = count_subtree_patches_Domain(grid(d), p_root)
  n_new = size(patch_copy)

  if (n_new /= n_old) then
     error stop &
          "test_subtree_extraction: patch count mismatch"
  end if

  if (old_to_new(p_root) /= 0) then
     error stop &
          "test_subtree_extraction: extracted root is not patch zero"
  end if

  !
  ! Save original storage offsets before compacting them.
  !
  allocate(old_elts_start(size(patch_copy)))

  do i = 1, size(patch_copy)
     old_elts_start(i) = patch_copy(i)%elts_start
  end do

  !
  ! Check original node-storage ranges.
  !
  do i = 1, size(patch_copy)

     if (old_elts_start(i) < 0) then
        error stop &
             "test_subtree_extraction: invalid original elts_start"
     end if

     if (old_elts_start(i) + PATCH_SIZE**2 > grid(d)%node%length) then
        error stop &
             "test_subtree_extraction: original node storage range out of bounds"
     end if

  end do

  !
  ! Assign compact local storage offsets.
  !
  call compact_subtree_storage_Domain(patch_copy)

  do i = 1, size(patch_copy)

     if (patch_copy(i)%elts_start /= (i-1)*PATCH_SIZE**2) then
        error stop &
             "test_subtree_extraction: incorrect compact elts_start"
     end if

  end do

  n_node_storage = size(patch_copy) * PATCH_SIZE**2

  !
  ! ---------------------------------------------------------------
  ! Copy and verify node coordinates.
  ! ---------------------------------------------------------------
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
          grid(d)%node%elts(old_start+1:old_start+PATCH_SIZE**2)%x)) > &
          0.0_dp) then
        error stop &
             "test_subtree_extraction: node x-coordinate mismatch"
     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%y - &
          grid(d)%node%elts(old_start+1:old_start+PATCH_SIZE**2)%y)) > &
          0.0_dp) then
        error stop &
             "test_subtree_extraction: node y-coordinate mismatch"
     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%z - &
          grid(d)%node%elts(old_start+1:old_start+PATCH_SIZE**2)%z)) > &
          0.0_dp) then
        error stop &
             "test_subtree_extraction: node z-coordinate mismatch"
     end if

  end do

  !
  ! ---------------------------------------------------------------
  ! Copy and verify one node-based scalar field.
  ! ---------------------------------------------------------------
  !
  v_scalar    = scalars(1)
  mult_scalar = MULT(v_scalar)
  k_test      = max(1,zmin)

  if (mult_scalar /= 1) then
     error stop &
          "test_subtree_extraction: test scalar has unexpected multiplier"
  end if

  call copy_subtree_field_Domain( &
       patch_copy, old_elts_start, mult_scalar, &
       sol(v_scalar,k_test)%data(d)%elts, scalar_copy)

  if (size(scalar_copy) /= mult_scalar*n_node_storage) then
     error stop &
          "test_subtree_extraction: incorrect copied scalar storage size"
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
  ! ---------------------------------------------------------------
  ! Copy and verify one edge-based velocity field.
  ! ---------------------------------------------------------------
  !
  v_vector    = S_VELO
  mult_vector = MULT(v_vector)

  if (mult_vector /= EDGE) then
     error stop &
          "test_subtree_extraction: velocity has unexpected multiplier"
  end if

  call copy_subtree_field_Domain( &
       patch_copy, old_elts_start, mult_vector, &
       sol(v_vector,k_test)%data(d)%elts, vector_copy)

  if (size(vector_copy) /= mult_vector*n_node_storage) then
     error stop &
          "test_subtree_extraction: incorrect copied vector storage size"
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
  ! ---------------------------------------------------------------
  ! Verify copied parent-child links.
  ! ---------------------------------------------------------------
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_CHDRN

        p_chd_old = grid(d)%patch%elts(p_old+1)%children(c)

        p_chd_new = &
             patch_copy(old_to_new(p_old)+1)%children(c)

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

  !
  ! ---------------------------------------------------------------
  ! Compare leaf counts and subtree depth.
  ! ---------------------------------------------------------------
  !
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
     error stop &
          "test_subtree_extraction: leaf count mismatch"
  end if

  depth_new = copied_depth(0)

  if (depth_new /= depth_old) then
     error stop &
          "test_subtree_extraction: subtree depth mismatch"
  end if

  !
  ! ---------------------------------------------------------------
  ! Construct neighbour classification before MOVE_ALLOC.
  !
  ! NGB_INTERNAL : positive neighbour remains inside this block
  ! NGB_BLOCK    : positive neighbour lies outside this block
  ! NGB_DOMAIN   : negative neighbour representing Domain boundary
  ! NGB_ADAPT    : negative neighbour representing adaptive boundary
  ! NGB_OTHER    : zero/other neighbour entry
  ! ---------------------------------------------------------------
  !
  allocate(block_test%neigh_class(N_BDRY,size(patch_copy)))
  block_test%neigh_class = NGB_OTHER

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

        if (p_ngb_old > 0) then

           if (old_to_new(p_ngb_old) >= 0) then

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_INTERNAL

           else

              block_test%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_BLOCK

           end if

        else if (p_ngb_old < 0) then

           if (grid(d)%bdry_patch%elts(-p_ngb_old+1)%side > 0) then

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

  !
  ! Validate internal-neighbour classification.
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        if (block_test%neigh_class( &
             c,old_to_new(p_old)+1) /= NGB_INTERNAL) cycle

        p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

        if (p_ngb_old <= 0) then
           error stop &
                "test_subtree_extraction: invalid internal neighbour"
        end if

        if (old_to_new(p_ngb_old) < 0) then
           error stop &
                "test_subtree_extraction: internal neighbour missing"
        end if

     end do

  end do

  !
  ! Validate negative-neighbour classifications against bdry_patch.
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

        if (p_ngb_old >= 0) cycle

        select case (block_test%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_DOMAIN)

           if (grid(d)%bdry_patch%elts(-p_ngb_old+1)%side <= 0) then
              error stop &
                   "test_subtree_extraction: incorrect domain boundary classification"
           end if

        case (NGB_ADAPT)

           if (grid(d)%bdry_patch%elts(-p_ngb_old+1)%side > 0) then
              error stop &
                   "test_subtree_extraction: incorrect adaptive boundary classification"
           end if

        case default

           error stop &
                "test_subtree_extraction: invalid negative neighbour classification"

        end select

     end do

  end do

  n_ngb_internal = count(block_test%neigh_class == NGB_INTERNAL)
  n_ngb_block    = count(block_test%neigh_class == NGB_BLOCK)
  n_ngb_domain   = count(block_test%neigh_class == NGB_DOMAIN)
  n_ngb_adapt    = count(block_test%neigh_class == NGB_ADAPT)
  n_ngb_other    = count(block_test%neigh_class == NGB_OTHER)

  if (n_ngb_internal + n_ngb_block + n_ngb_domain + &
       n_ngb_adapt + n_ngb_other /= N_BDRY*size(patch_copy)) then

     error stop &
          "test_subtree_extraction: neighbour classification count mismatch"

  end if

  !
  ! ---------------------------------------------------------------
  ! Package validated data into the temporary block.
  ! ---------------------------------------------------------------
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
  ! Verify container dimensions.
  !
  if (size(block_test%patch) /= n_old) then
     error stop &
          "test_subtree_extraction: block patch size mismatch"
  end if

  if (size(block_test%node) /= n_node_storage) then
     error stop &
          "test_subtree_extraction: block node size mismatch"
  end if

  if (size(block_test%scalar) /= n_node_storage) then
     error stop &
          "test_subtree_extraction: block scalar size mismatch"
  end if

  if (size(block_test%vector) /= EDGE*n_node_storage) then
     error stop &
          "test_subtree_extraction: block vector size mismatch"
  end if

  if (size(block_test%neigh_class,1) /= N_BDRY .or. &
       size(block_test%neigh_class,2) /= size(block_test%patch)) then

     error stop &
          "test_subtree_extraction: block neighbour-class size mismatch"

  end if

  !
  ! Verify compact patch layout survived MOVE_ALLOC.
  !
  do i = 1, size(block_test%patch)

     if (block_test%patch(i)%elts_start /= &
          (i-1)*PATCH_SIZE**2) then

        error stop &
             "test_subtree_extraction: invalid block storage layout"

     end if

  end do

  !
  ! Source allocatables must no longer own storage.
  !
  if (allocated(patch_copy)) then
     error stop &
          "test_subtree_extraction: patch_copy still allocated"
  end if

  if (allocated(node_copy)) then
     error stop &
          "test_subtree_extraction: node_copy still allocated"
  end if

  if (allocated(scalar_copy)) then
     error stop &
          "test_subtree_extraction: scalar_copy still allocated"
  end if

  if (allocated(vector_copy)) then
     error stop &
          "test_subtree_extraction: vector_copy still allocated"
  end if

  !
  ! ---------------------------------------------------------------
  ! Diagnostic output.
  ! ---------------------------------------------------------------
  !
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
       "  node-coordinate copy check passed"

  write(6,'(a,i0,a,i0,a,i0)') &
       "  scalar field copy check passed: variable ", v_scalar, &
       ", level ", k_test, &
       ", multiplier ", mult_scalar

  write(6,'(a,i0,a,i0,a,i0)') &
       "  vector field copy check passed: variable ", v_vector, &
       ", level ", k_test, &
       ", multiplier ", mult_vector

  write(6,'(a,i0)') &
       "  temporary block ID = ", block_test%id

  write(6,'(a,i0)') &
       "  temporary block patches = ", size(block_test%patch)

  write(6,'(a,i0)') &
       "  temporary block nodes = ", size(block_test%node)

  write(6,'(a,i0)') &
       "  temporary block scalar values = ", size(block_test%scalar)

  write(6,'(a,i0)') &
       "  temporary block vector values = ", size(block_test%vector)

  write(6,'(a)') &
       "  temporary block container check passed"

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

  write(6,'(a)') &
       "  block boundary classification check passed"

  write(6,'(a,/)') &
       "  patch topology and storage layout checks passed"

  deallocate(old_to_new)
  deallocate(old_elts_start)

contains

  recursive integer function copied_depth (p) result(depth)
    implicit none

    integer, intent(in) :: p

    integer :: c, p_chd

    depth = 0

    do c = 1, N_CHDRN

       p_chd = patch_copy(p+1)%children(c)

       if (p_chd > 0) then
          depth = max(depth, 1 + copied_depth(p_chd))
       end if

    end do

  end function copied_depth

end subroutine test_subtree_extraction


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
