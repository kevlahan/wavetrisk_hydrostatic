module main_mod
  use io_mod
  use wavelet_mod
  use adapt_mod
  use remap_mod

  implicit none
  type Initial_State
     integer                                          :: n_patch, n_bdry_patch, n_node, n_edge, n_tria
     integer, dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack_len, unpk_len
  end type Initial_State
  integer                                        :: chkpt_info
  integer,             dimension(:), allocatable :: node_level_start, edge_level_start
  real(8)                                        :: dt_new, time_mult  
  type(Initial_State), dimension(:), allocatable :: ini_st
contains
  subroutine init_basic
    implicit none
    call init_comm_mod
    call init_init_mod
    call init_refine_patch_mod
    call init_time_integr_mod
    call init_io_mod
    call init_wavelet_mod
    call init_mask_mod
    call init_adapt_mod
    time_mult = 1d0
  end subroutine init_basic

  subroutine initialize (run_id)
    ! Initialize from checkpoint or adapt to initialize conditions
    ! Solution is saved and restarted to balance load
    implicit none
    character(*) :: run_id
    
    character(255) :: command
    integer        :: k, d, v

    ! Default elliptic solver (scheduled relaxation Jacobi method)
    elliptic_solver => SRJ
    
    if (max_level < min_level) then
       if (rank == 0) then
          write (6,'(//,a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a)')   "!!                                                                                 !!"
          write (6,'(2(a,i2),a)') "!!                max_level < max_level: ", max_level, " < ", min_level, &
               " ... aborting                      !!"
          write (6,'(a)')   "!!                                                                                 !!"
          write (6,'(a)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
          write (6,'(a,//)') "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
       end if
       call abort
    end if

    ! Check validity of various parameter choices
    if (compressible .and. mode_split) then
       if (rank == 0) write (6,'(a)') "Cannot use mode splitting with compressible dynamics ... aborting"
       call abort
    end if

#ifdef AMPI
    call MPI_Info_create (chkpt_info, ierror)
    call MPI_Info_set (chkpt_info, "ampi_checkpoint", "to_file=checkpoint", ierror)
#endif

    if (resume >= 0) then
       cp_idx = resume
       call restart (run_id)
       resume = NONE
    else
       ! Initialize basic structures
       call init_basic
       call init_structures (run_id)
       
       ! Initialize vertical grid
       call initialize_a_b_vert

       ! Determine vertical level to save
       call set_save_level

       ! Initialize time step and viscosities
       call initialize_dt_viscosity

        ! Initialize thresholds to default values 
       call initialize_thresholds

       ! Load topography
       if (NCAR_topo .and. trim (test_case) /= 'make_NCAR_topo') call load_topo
      
       if (rank == 0) write (6,'(/,A,/)') &
            '----------------------------------------------------- Adapting initial grid &
            ------------------------------------------------------'

       call apply_initial_conditions
       call forward_wavelet_transform (sol, wav_coeff)

       do while (level_end < max_level)
          if (rank == 0) write (6,'(A,i2,A,i2)') 'Initial refinement Level ', level_end, ' -> ', level_end+1
          node_level_start = grid(:)%node%length+1
          edge_level_start = grid(:)%midpt%length+1

          dt_new = cpt_dt ()
          call adapt (set_thresholds)

          call apply_initial_conditions
          call forward_wavelet_transform (sol, wav_coeff)

          ! Check whether there are any active nodes or edges at this scale
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
          if (rank == 0) write (6,'(A,i2,1x,2(A,i12,1x),/)') &
               'Level = ', level_end, 'number of active node wavelets = ', n_active(AT_NODE), &
               'number of active edge wavelets = ', n_active(AT_EDGE)
          if (n_active(AT_NODE) == 0 .and. n_active(AT_EDGE) == 0) exit ! No active nodes or edges at this scale
       end do
       if (rank == 0) write (6,'(A,/)') &
            '------------------------------------------------- Finished adapting initial grid &
            -------------------------------------------------'

       call adapt (set_thresholds) ; dt_new = cpt_dt ()
       if (rank==0) write (6,'(a,i8,/)') 'Initial number of dof = ', sum (n_active)

       if (trim (test_case) /= 'make_NCAR_topo') call write_checkpoint (run_id, .true.)
    end if
    call barrier
  end subroutine initialize

  subroutine record_init_state (init_state)
    implicit none
    type(Initial_State), dimension(:), allocatable :: init_state

    integer :: d, i, v

    allocate (init_state(size(grid)))

    do d = 1, size(grid)
       init_state(d)%n_patch      = grid(d)%patch%length
       init_state(d)%n_bdry_patch = grid(d)%bdry_patch%length
       init_state(d)%n_node       = grid(d)%node%length
       init_state(d)%n_edge       = grid(d)%midpt%length
       init_state(d)%n_tria       = grid(d)%ccentre%length
       do i = 1, N_GLO_DOMAIN
          do v = AT_NODE, AT_EDGE
             init_state(d)%pack_len(v,i) = grid(d)%pack(v,i)%length 
             init_state(d)%unpk_len(v,i) = grid(d)%unpk(v,i)%length 
          end do
       end do
    end do
  end subroutine record_init_state

  subroutine time_step (align_time, aligned)
    use vert_diffusion_mod
    implicit none
    real(8)              :: align_time
    logical, intent(out) :: aligned

    integer(8) :: idt, ialign
    real(8)    :: dx, dt_0
    
    istep       = istep+1
    istep_cumul = istep_cumul+1
    
    dt = dt_new

    n_patch_old = grid%patch%length
    n_node_old  = grid%node%length

    idt    = nint (dt*time_mult, 8)
    ialign = nint (align_time*time_mult, 8)
    if (ialign > 0 .and. istep > 20) then
       aligned = (modulo (itime+idt,ialign) < modulo (itime,ialign))
    else
       aligned = .false.
    end if

    ! Modify time step
    if (aligned .and. match_time) then
       idt = ialign - modulo (itime,ialign)
       dt = idt/time_mult
    end if

    ! Take time step
    if (mode_split) then ! 2D barotropic mode splitting (implicit Euler)
       if (istep <= nstep_init) then  ! small time steps to start
          dx = sqrt (4d0/sqrt(3d0) * 4d0*MATH_PI*radius**2/(20d0*4**level_end)) 
          dt_0 = 0.8d0 * dx / wave_speed

          dt = dt_0 + (dt - dt_0) * sin (MATH_PI/2d0 * dble(istep-1)/dble(nstep_init-1))
       end if

       select case (timeint_type)
       case ("RK4")
          call RK4_split (dt)
       case ("RK3")
          call RK3_split (dt)
       case ("RK2")
          call RK2_split (dt)
       case ("Euler")
          call Euler_split (dt)
       case default
          call RK3_split (dt)
       end select
    else
       select case (timeint_type)
       case ("Euler")
          call Euler (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case ("RK33")
          call RK33_opt (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case ("RK34")
          call RK34_opt (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case ("RK45")
          call RK45_opt (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case ("RK3")
          call RK3 (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case ("RK4")
          call RK4 (sol(1:N_VARIABLE,1:zlevels), wav_coeff(1:N_VARIABLE,1:zlevels), trend_ml, dt)
       case default
          if (rank == 0) write (6,'(a)') "Invalid timestepping choice ... aborting"
          call abort
       end select
    end if

    ! Split step routines
    if (vert_diffuse) call vertical_diffusion 

    ! If necessary, remap vertical coordinates
    if (remap .and. modulo (istep, iremap) == 0) call remap_vertical_coordinates

    ! Compute change in vertical layer depths
    if (log_mass) min_mass = cpt_min_mass ()
    
    ! Add diffusion
    if (modulo (istep_cumul, n_diffuse) == 0 .and. Laplace_order_init /= 0) then
       Laplace_order = Laplace_order_init
    else
       Laplace_order = 0
    end if
    
    ! Adapt grid
    if (vert_diffuse .or. (remap .and. modulo (istep, iremap) == 0)) &
         call WT_after_step (sol(:,1:zlevels), wav_coeff(:,1:zlevels), level_start-1)
    if (zmin < 1) &
         call WT_after_step (sol(:,zmin:0), wav_coeff(:,zmin:0), level_start-1)
    call adapt (set_thresholds)
    call inverse_wavelet_transform (wav_coeff, sol, jmin_in=level_start)

    call update
    
    if (log_mass) call sum_total_mass (.false.)

    itime = itime + idt
    
    if (match_time) then
       time = itime/time_mult
    else
       time = time + dt
    end if
    
    ! Set new time step and count active nodes
    if (modulo (istep, iadapt) == 0) dt_new = cpt_dt ()

    ! Check load balance
#ifdef AMPI
    if (modulo (istep, irebalance) == 0) then
       if (rank == 0) write (6,'(a)') "Checking load balance and rebalancing if necessary using AMPI ..."
       call MPI_Barrier (MPI_COMM_WORLD, ierror)
       call AMPI_Migrate (AMPI_INFO_LB_SYNC, ierror)
    end if
#endif
  end subroutine time_step

  subroutine restart (run_id)
    ! Fresh restart from checkpoint data (all structures reset)
    implicit none
    character(*) :: run_id
    
    integer         :: l
    character(9999) :: bash_cmd, cmd_archive, cmd_files, command
    
    if (maxval (C_visc(S_MASS:S_TEMP)) > (1d0/6d0)**Laplace_order .or. C_visc(S_VELO) > (1d0/24d0)**Laplace_order) then
       if (rank == 0) write (6,*) "Dimensional viscosity too large ... aborting"
       call abort
    end if

    if (rank == 0) then
       write (6,'(A,/)') &
            '********************************************************* Begin Restart &
            **********************************************************'
       write (6,'(A,i4,/)') 'Reloading from checkpoint ', cp_idx
    end if

    ! Deallocate all dynamic arrays and variables
    if (resume == NONE) call deallocate_structures

    ! Initialize basic structures
    call init_basic

    ! Initialize vertical grid
    call initialize_a_b_vert
    
    ! Determine vertical level to save
    call set_save_level

    ! Uncompress checkpoint data (needed for init_structures and load_adapt_mpi)
    if (rank == 0) then
       write (cmd_archive, '(a,i4.4,a)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (6,              '(a,a,/)') 'Loading file ', trim (cmd_archive)
       write (command,          '(a,a)') 'gtar xzf ', trim (cmd_archive)
       write (bash_cmd,       '(a,a,a)') 'bash -c "', trim (command), '"'
       call system (trim(bash_cmd))
    end if
    call barrier ! make sure all archive files have been uncompressed

    ! Rebalance adaptive grid and re-initialize structures
    call init_structures (run_id)
    
    ! Load NCAR topography data (defined on non-adaptive grid from min_level to max_level)
    if (NCAR_topo) call load_topo 

    ! Initialize thresholds to default values 
    call initialize_thresholds

    ! Load checkpoint data
    call load_adapt_mpi (cp_idx, run_id)

    ! Compute masks based on active wavelets
    ! (do not re-calculate thresholds)
    call adapt (set_thresholds, .false.) 
    call inverse_wavelet_transform (wav_coeff, sol, jmin_in=level_start-1)
    if (vert_diffuse) call inverse_scalar_transform (wav_tke, tke, jmin_in=level_start-1)

    ! Initialize time step and viscosities
    call initialize_dt_viscosity

    ! Initialize mean values and other test case defined variables
    call update

    ! Initialize total mass value
    if (log_mass) call sum_total_mass (.true.)
    
    ! Initialize time step and counters
    dt_new = cpt_dt ()
    itime = nint (time*time_mult, 8)
    istep = 0

    if (rank == 0) then
       write (6,'(/,A,es12.6,3(A,es8.2),A,I2,A,I9,/)') &
            'time [d] = ', time/DAY, &
            '  mass threshold = ', sum (threshold(S_MASS,:))/size(threshold,2), &
            ' temp threshold = ', sum (threshold(S_TEMP,:))/size(threshold,2), &
            ' velo threshold = ', sum (threshold(S_VELO,:))/size(threshold,2), &
            ' Jmax = ', level_end, &
            '  dof = ', sum (n_active)
       write (6,'(A)') &
            '********************************************************** End Restart &
            ***********************************************************'
    end if

#ifdef AMPI
    if (rank == 0) write (6,'(/,a)') "Rebalancing using AMPI ..."
    call MPI_Barrier (MPI_COMM_WORLD, ierror)
    call AMPI_Migrate (AMPI_INFO_LB_SYNC, ierror)
#endif
  end subroutine restart

  subroutine write_checkpoint (run_id, rebal)
    implicit none
    character(*) :: run_id
    logical      :: rebal
    
    cp_idx = cp_idx + 1 

    if (rank == 0) then
       write (6,'(/,a,/)') &
            '************************************************************************&
            **********************************************************'
       write (6,'(a,i4,a,es10.4,/)') 'Saving checkpoint ', cp_idx, ' at time [day] = ', time/DAY
    end if
    
! #ifdef AMPI
!     if (rank == 0) write (6,'(a)') "Checkpointing using AMPI ..."
!     call MPI_Info_set (chkpt_info, "ampi_checkpoint", "to_file=checkpoint", ierror)
!     call MPI_Barrier (MPI_COMM_WORLD, ierror)
!     call AMPI_Migrate (chkpt_info, ierror)
!     if (log_mass) call sum_total_mass (.true.) 
! #else
    call write_load_conn (cp_idx, run_id)
    call dump_adapt_mpi  (cp_idx, run_id)
    
    call restart (run_id)
!#endif
  end subroutine write_checkpoint

  subroutine init_structures (run_id)
    ! Initialize dynamical arrays and structures
    implicit none
    character(*) :: run_id

    level_start = min_level
    level_end = level_start
    
    ! Distribute and balance grid over processors (necessary for correct restart!)
    call distribute_grid (cp_idx, run_id)
    
    call init_grid
    call init_comm_mpi
    call init_geometry

    if (optimize_grid == XU_GRID) call smooth_Xu (1.0d6*eps())
    if (optimize_grid == HR_GRID) call read_HR_optim_grid

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call precompute_geometry

    allocate (node_level_start(size(grid)), edge_level_start(size(grid)))
    edge_level_start = 0
    node_level_start = 0

    if (rank == 0) write (6,'(A,i2,A)') 'Make level J_min = ', min_level, ' ...'
    
    call init_wavelets
    call init_masks
    call add_second_level

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, +BDRY_THICKNESS)
    call apply_interscale (mask_adj_scale, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)
    if (time_end > 0d0) time_mult = huge (itime)/2/time_end

    allocate (n_patch_old(size(grid)), n_node_old(size(grid))); n_patch_old = 2

    call init_RK_mem
  end subroutine init_structures

  real(8) function cpt_dt ()
    ! Calculates time step, minimum relative mass and active nodes and edges
    implicit none
    integer               :: l, ierror, level_end_glo
    integer, dimension(2) :: n_active_loc
    
    if (adapt_dt) dt_loc = 1d16
    n_active_nodes = 0
    n_active_edges = 0

    ! Calculate minimum time step, number of active nodes and edges
    do l = level_start, level_end
       call apply_onescale (cal_min_dt, l, z_null, 0, 0)
    end do

    ! Time step
    if (adapt_dt) then
       cpt_dt = sync_min_real (dt_loc)
    else
       cpt_dt = dt_init
    end if

    ! Active nodes and edges
    n_active_loc = (/ sum (n_active_nodes(level_start:level_end)), sum(n_active_edges(level_start:level_end)) /)
    n_active = sum_int_vector (n_active_loc, 2)
    level_end = sync_max_int (level_end)
  end function cpt_dt

  real(8) function cpt_min_mass ()
    ! Calculates minimum relative mass
    implicit none
    integer :: ierror, l

    min_mass_loc = 1d16
    do l = level_start, level_end
       call apply_onescale (cal_min_mass, l, z_null, 0, 0)
    end do

    cpt_min_mass = sync_min_real (min_mass_loc)
  end function cpt_min_mass

  subroutine cal_min_dt (dom, i, j, zlev, offs, dims)
    ! Calculates time step and number of active nodes and edges
    ! time step is smallest of barotropic time step, advective time step and internal wave time step for mode split case
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer :: d, e, id, id_e, id_i, k, l
    real(8) :: dx, v_mag

    id = idx (i, j, offs, dims)
    id_i = id + 1
    d  = dom%id + 1
    l  = dom%level%elts(id_i)
        
    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       n_active_nodes(l) = n_active_nodes(l) + 1 
       if (adapt_dt) then
          dx = minval (dom%len%elts(EDGE*id+RT+1:EDGE*id+UP+1))
          do k = 1, zlevels
             v_mag = maxval (abs(sol(S_VELO,k)%data(d)%elts(EDGE*id+RT+1:EDGE*id+UP+1)))
             if (mode_split) then
                dt_loc = min (dt_loc, cfl_num * dx / wave_speed, cfl_adv * dx / v_mag, cfl_bar * dx / c1)
             else
                dt_loc = min (dt_loc, cfl_num * dx / (v_mag + wave_speed))
             end if
          end do
       end if
    end if

    do e = 1, EDGE
       id_e = EDGE * id + e 
       if (dom%mask_e%elts(id_e) >= ADJZONE) n_active_edges(l) = n_active_edges(l) + 1
    end do
  end subroutine cal_min_dt

  subroutine cal_min_mass (dom, i, j, zlev, offs, dims)
    ! Calculates minimum relative mass and checks diffusion stability limits
    use init_mod
    implicit none
    type(Domain)                   :: dom
    integer                        :: i, j, zlev
    integer, dimension(N_BDRY+1)   :: offs
    integer, dimension(2,N_BDRY+1) :: dims

    integer                       :: d, e, id, id_e, id_i, k, l
    real(8)                       :: col_mass, d_e, fac, full_mass, init_mass, rho, z_s
    real(8)                       :: beta_sclr, beta_divu, beta_rotu
    real(8), dimension(1:zlevels) :: dz
    real(8), dimension(0:zlevels) :: z

    id   = idx (i, j, offs, dims)
    id_i = id + 1
    d    = dom%id + 1

    if (dom%mask_n%elts(id_i) >= ADJZONE) then
       col_mass = 0d0
       do k = 1, zlevels
          full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
          if (full_mass < 0d0 .or. full_mass /= full_mass) then
             write (6,'(A,i8,A,3(es9.2,1x),A,i2,A)') "Mass negative at id = ", id_i, &
                  " with position ", dom%node%elts(id_i)%x,  dom%node%elts(id_i)%y, dom%node%elts(id_i)%z, &
                  " vertical level k = ", k, " ... aborting"
             call abort
          end if
          col_mass = col_mass + full_mass
       end do

       ! Measure relative change in mass
       if (compressible) then
          do k = 1, zlevels
             init_mass = a_vert_mass(k) + b_vert_mass(k)*col_mass
             min_mass_loc = min (min_mass_loc, sol(S_MASS,k)%data(d)%elts(id_i)/init_mass)
          end do
       else
          z_s = dom%topo%elts(id_i)
          if (sigma_z) then
             z = z_coords (0d0, z_s)
          else
             z = b_vert * z_s ! assumes zero free surface perturbation initial condition
          end if
          dz = z(1:zlevels) - z(0:zlevels-1)
          do k = 1, zlevels
             rho = porous_density (d, id_i, k)
             init_mass = rho * dz(k)
             full_mass = sol(S_MASS,k)%data(d)%elts(id_i) + sol_mean(S_MASS,k)%data(d)%elts(id_i)
             min_mass_loc = min (min_mass_loc, full_mass/init_mass)
          end do
       end if
    end if
  end subroutine cal_min_mass
  
  integer function write_active_per_level ()
    ! Write out distribution of active nodes over levels
    implicit none
    integer                                         :: l, n_full, fillin, n_lev_cur, recommended_level_start
    integer, dimension(2*(level_end-level_start+1)) :: n_active_all_loc, n_active_all
    integer, dimension(level_start:level_end)       :: n_active_per_lev
    real(8)                                         :: dt

    dt = cpt_dt () ! to set n_active_*

    n_lev_cur = level_end - level_start + 1
    
    n_active_all_loc = (/n_active_nodes(level_start:level_end), n_active_edges(level_start:level_end)/)
    
    ! Sum n_active_all_loc up across all processes and distribute result n_active_all_glo among all processes
    n_active_all = sum_int_vector (n_active_all_loc, n_lev_cur*2)
    
    n_active_nodes(level_start:level_end) = n_active_all(1:n_lev_cur)
    n_active_edges(level_start:level_end) = n_active_all(n_lev_cur+1:n_lev_cur*2)
    n_active_per_lev = n_active_edges(level_start:level_end) + n_active_nodes(level_start:level_end)

    if (rank == 0) write (6,'(6X,A,A,3(1X,A))') '   N_p   ','   N_u   ','of all active', 'of full level', 'fill-in'

    recommended_level_start = level_start

    do l = level_start, level_end
       n_full = max_nodes_per_level(l) + max_nodes_per_level(l,EDGE)

       ! Fill-in: additional nodes on level `l` if it'd become lowest level 
       ! minus the nodes on lower levels which would be removed
       fillin = n_full-n_active_per_lev(l)-sum(n_active_per_lev(level_start:l-1))

       if (rank == 0) then
          write (6,'(A,I2,I9,I9,2(1X,F9.1,A),1X,I9,1X,F9.1,A)') &
               'lev', l, n_active_nodes(l), n_active_edges(l), &
               float(n_active_per_lev(l))/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%', &
               float(n_active_per_lev(l))/float(n_full)*100.0, '%', &
               fillin, float(fillin)/float(sum(n_active(AT_NODE:AT_EDGE)))*100.0, '%'
       end if

       if (fillin <= 0) recommended_level_start = l
    end do

    if (rank == 0) then
       write (6,'(A,I9,I9,2(1X,F9.1,A),9X,I9)') 'total', n_active(AT_NODE:AT_EDGE), 100.0, '%', &
            float(sum(n_active(AT_NODE:AT_EDGE)))/float(n_full)*100.0, '%', &
            n_full/sum(n_active(AT_NODE:AT_EDGE))
    end if

    write_active_per_level = recommended_level_start
  end function write_active_per_level

  subroutine deallocate_structures
    ! Deallocate all dynamic arrays and structures for clean restart
    implicit none

    integer :: d, i, k, l, v, r

    ! Deallocate init_RK_mem allocations
    do k = 1, zmax
       do d = 1, n_domain(rank+1)
          do v = 1, N_VARIABLE
             deallocate (q1(v,k)%data(d)%elts)
             deallocate (q2(v,k)%data(d)%elts)
             deallocate (q3(v,k)%data(d)%elts)
             deallocate (q4(v,k)%data(d)%elts)
             deallocate (dq1(v,k)%data(d)%elts)
          end do
       end do
       do v = 1, N_VARIABLE
          deallocate (q1(v,k)%data)
          deallocate (q2(v,k)%data)
          deallocate (q3(v,k)%data)
          deallocate (q4(v,k)%data)
          deallocate (dq1(v,k)%data)
       end do
    end do
    deallocate (q1, q2, q3, q4, dq1)

    ! Deallocate grid structure elements
    do d = 1, size(grid)
       deallocate (grid(d)%mask_n%elts)
       deallocate (grid(d)%mask_e%elts)
       
       deallocate (grid(d)%level%elts)
       
       deallocate (grid(d)%R_F_wgt%elts)
       deallocate (grid(d)%I_u_wgt%elts)
       
       deallocate (grid(d)%overl_areas%elts)
       deallocate (grid(d)%triarea%elts)
       deallocate (grid(d)%len%elts)
       deallocate (grid(d)%pedlen%elts)
       deallocate (grid(d)%areas%elts)
       deallocate (grid(d)%midpt%elts)
       deallocate (grid(d)%ccentre%elts)

       deallocate (grid(d)%surf_press%elts)
       deallocate (grid(d)%press%elts)
       deallocate (grid(d)%geopot%elts)
       deallocate (grid(d)%u_zonal%elts)
       deallocate (grid(d)%v_merid%elts)
       deallocate (grid(d)%press_lower%elts)
       deallocate (grid(d)%geopot_lower%elts)
       deallocate (grid(d)%vort%elts)
       deallocate (grid(d)%qe%elts)
       deallocate (grid(d)%bernoulli%elts)
       deallocate (grid(d)%ke%elts)
       deallocate (grid(d)%divu%elts)
       deallocate (grid(d)%topo%elts)
       deallocate (grid(d)%coriolis%elts)
       
       deallocate (grid(d)%node%elts) 
       deallocate (grid(d)%bdry_patch%elts) 
       deallocate (grid(d)%patch%elts) 
       deallocate (grid(d)%neigh_pa_over_pole%elts)
       deallocate (grid(d)%send_pa_all%elts)
       
       do i = 1, N_GLO_DOMAIN
          deallocate (grid(d)%recv_pa(i)%elts)
          deallocate (grid(d)%send_conn(i)%elts)
          do k = AT_NODE, AT_EDGE
             deallocate (grid(d)%pack(k,i)%elts)
             deallocate (grid(d)%unpk(k,i)%elts)
          end do
       end do

       do i = lbound(grid(d)%lev,1), ubound(grid(d)%lev,1)
          deallocate (grid(d)%lev(i)%elts)
       end do
       deallocate (grid(d)%lev)

       do l = min_level, max_level
          do r = 1, n_process
             deallocate (grid(d)%src_patch(r,l)%elts) 
          end do
       end do
       deallocate (grid(d)%src_patch)

       deallocate (Laplacian_vector(S_DIVU)%data(d)%elts)
       deallocate (Laplacian_vector(S_ROTU)%data(d)%elts)
       
       do v = scalars(1), scalars(2)
          deallocate (horiz_flux(v)%data(d)%elts)
          deallocate (Laplacian_scalar(v)%data(d)%elts)
       end do

       do k = zmin, zmax
          deallocate (penal_node(k)%data(d)%elts)
          deallocate (penal_edge(k)%data(d)%elts)
          deallocate (exner_fun(k)%data(d)%elts)
       end do
       deallocate (exner_fun(zmax+1)%data(d)%elts)

       do v = 1, N_VARIABLE
          do k = zmin, zmax
             deallocate (sol(v,k)%data(d)%elts)
             deallocate (sol_mean(v,k)%data(d)%elts)
             if (k > 0) deallocate (trend(v,k)%data(d)%elts)
             deallocate (wav_coeff(v,k)%data(d)%elts)
          end do
         
          do k = 1, save_levels
             deallocate (sol_save(v,k)%data(d)%elts) 
          end do
       end do
       
       if (vert_diffuse) then
          deallocate (Kt(0)%data(d)%elts)
          deallocate (Kv(0)%data(d)%elts)
          do k = 1, zlevels
             deallocate (Kt(k)%data(d)%elts)
             deallocate (Kv(k)%data(d)%elts)
             deallocate (tke(k)%data(d)%elts)
             deallocate (wav_tke(k)%data(d)%elts)
          end do
       end if
    end do

    deallocate (topography%data)
    deallocate (wav_topography%data)

    deallocate (Laplacian_vector(S_DIVU)%data)
    deallocate (Laplacian_vector(S_ROTU)%data)

    do k = zmin, zmax
       deallocate (penal_node(k)%data)
       deallocate (penal_edge(k)%data)
       deallocate (exner_fun(k)%data)
    end do
    deallocate (exner_fun(zmax+1)%data)
    
    do v = scalars(1), scalars(2)
       deallocate (horiz_flux(v)%data)
       deallocate (Laplacian_scalar(v)%data)
    end do

    do v = 1, N_VARIABLE
       do k = zmin, zmax
          deallocate (sol(v,k)%data)
          deallocate (sol_mean(v,k)%data)
          if (k > 0) deallocate (trend(v,k)%data)
          deallocate (wav_coeff(v,k)%data)
       end do
       do k = 1, save_levels
          deallocate (sol_save(v,k)%data)
       end do
    end do
    
    if (vert_diffuse) then
       deallocate (Kv(0)%data)
       deallocate (Kt(0)%data)
       do k = 1, zlevels
          deallocate (Kv(k)%data)
          deallocate (Kt(k)%data)
          deallocate (tke(k)%data)
          deallocate (wav_tke(k)%data)
       end do
    end if

    deallocate (grid, n_patch_old, n_node_old)
    deallocate (edge_level_start, node_level_start, n_active_edges, n_active_nodes)
    deallocate (a_vert, b_vert, a_vert_mass, b_vert_mass)
    deallocate (threshold, threshold_def)
    deallocate (sol, sol_mean, sol_save, trend, wav_coeff)       
    deallocate (exner_fun, horiz_flux, Laplacian_scalar, Laplacian_vector, lnorm, penal_node, penal_edge)
    deallocate (glo_id, ini_st, recv_lengths, recv_offsets, req, send_lengths, send_offsets)
    if (vert_diffuse) deallocate (Kv, Kt, tke, wav_tke)

    nullify (mass, dscalar, h_flux, velo, dvelo, bernoulli, divu, exner, ke, qe, scalar, temp, vort, wc_s, wc_u)
  end subroutine deallocate_structures
end module main_mod
