module main_mod
  use time_integr_mod
  use io_mod
  use wavelet_mod
  use adapt_mod
  use remap_mod

  implicit none
  type Initial_State
     integer                                          :: n_patch, n_bdry_patch, n_node, n_edge, n_tria
     integer, dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack_len, unpk_len
  end type Initial_State
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
    time_mult = 1.0_8
  end subroutine init_basic

  subroutine initialize (run_id)
    ! Initialize from checkpoint or adapt to initialize conditions
    ! Solution is saved and restarted to balance load
    implicit none
    character(*) :: run_id
    
    character(255) :: command
    integer        :: k, d, v

    ! Check validity of various parameter choices
    if (compressible .and. mode_split) then
       write (6,'(a)') "Cannot use mode splitting with compressible dynamics ... aborting"
       call abort
    end if

    if (max_level < min_level) then
       if (rank == 0) write (6,'(a)') "Max_level not >= min_level ... aborting"
       call abort
    end if

    if (resume >= 0) then
       cp_idx = resume
       call restart (run_id)
       resume = NONE
    else
       ! Initialize basic structures
       call init_basic

       ! Initialize vertical grid
       call initialize_a_b_vert

       ! Determine vertical level to save
       call set_save_level

       ! Initialize time step and viscosities
       call initialize_dt_viscosity

       call init_structures (run_id)
       call apply_initial_conditions
      
       ! Initialize thresholds to default values 
       call initialize_thresholds

       if (rank == 0) write (6,'(/,A,/)') &
            '----------------------------------------------------- Adapting initial grid &
            -----------------------------------------------------'

       call forward_wavelet_transform (sol, wav_coeff)
       if (adapt_trend) then
          call trend_ml (sol, trend)
          call forward_wavelet_transform (trend, trend_wav_coeff)
       end if

       do while (level_end < max_level)
          if (rank == 0) write (6,'(A,i2,A,i2)') 'Initial refinement Level ', level_end, ' -> ', level_end+1
          node_level_start = grid(:)%node%length+1
          edge_level_start = grid(:)%midpt%length+1
          
          dt_new = cpt_dt()
          call adapt (set_thresholds)

          call apply_initial_conditions
         
          call forward_wavelet_transform (sol, wav_coeff)
          if (adapt_trend) then
             call trend_ml (sol, trend)
             call forward_wavelet_transform (trend, trend_wav_coeff)
          end if

          ! Check whether there are any active nodes or edges at this scale
          n_active = 0
          do k = 1, zmax
             do d = 1, size(grid)
                do v = scalars(1), scalars(2)
                   if (adapt_trend) then
                      wc_s => trend_wav_coeff(v,k)%data(d)%elts
                   else
                      wc_s => wav_coeff(v,k)%data(d)%elts
                   end if
                   n_active(AT_NODE) = n_active(AT_NODE) &
                        + count (abs(wc_s(node_level_start(d):grid(d)%node%length)) >= threshold(v,k))
                   nullify (wc_s)
                end do
                if (adapt_trend) then
                   wc_u => trend_wav_coeff(S_VELO,k)%data(d)%elts
                else
                   wc_u => wav_coeff(S_VELO,k)%data(d)%elts
                end if
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
            ------------------------------------------------'

       call adapt (set_thresholds) ; dt_new = cpt_dt()
       if (rank==0) write (6,'(A,i8,/)') 'Initial number of dof = ', sum (n_active)
       
       call write_checkpoint (run_id, .true.)
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

  subroutine time_step (align_time, aligned, bottom_friction, wind_d, source_b, source_t)
    use lateral_diffusion_mod
    use vert_diffusion_mod
    implicit none
    real(8)              :: align_time
    logical, intent(out) :: aligned
    real(8), optional    :: bottom_friction
    optional             :: wind_d, source_b, source_t

    integer(8) :: idt, ialign
    real(8)    :: dx

    interface
     real(8) function source_b (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function source_b
     real(8) function source_t (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
     end function source_t
     function wind_d (dom, i, j, z_lev, offs, dims)
       use domain_mod
       implicit none
       type(Domain)                   :: dom
       integer                        :: i, j, z_lev
       integer, dimension(N_BDRY+1)   :: offs
       integer, dimension(2,N_BDRY+1) :: dims
       real(8), dimension(1:EDGE)     :: wind_d
     end function wind_d
  end interface
    
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
       if (istep <= 5) then  ! small time steps to start
          dx = sqrt (4/sqrt(3.0_8) * 4*MATH_PI*radius**2/(20*4**level_end)) 
          dt = 0.1*dx/wave_speed
       end if
       
       select case (timeint_type)
       case ("Euler")
          call Euler_split (dt)
       case ("RK4")
          call RK4_split (dt)
       case default
          if (rank == 0) &
               write (6,'(a)') "Invalid timestepping choice ... only Euler and RK4 supported for free surface ... aborting"
          call abort
       end select
    else
       select case (timeint_type)
       case ("Euler")
          call Euler (sol, wav_coeff, trend_ml, dt)
       case ("RK33")
          call RK33_opt (sol, wav_coeff, trend_ml, dt)
       case ("RK34")
          call RK34_opt (sol, wav_coeff, trend_ml, dt)
       case ("RK45")
          call RK45_opt (sol, wav_coeff, trend_ml, dt)
       case ("RK4")
          call RK4 (sol, wav_coeff, trend_ml, dt)
       case default
          if (rank == 0) write (6,'(a)') "Invalid timestepping choice ... aborting"
          call abort
       end select
    end if

    ! Split step routines
    if (implicit_diff_sclr .or. implicit_diff_divu) call implicit_lateral_diffusion
    if (vert_diffuse) call implicit_vertical_diffusion (bottom_friction, wind_d, source_b, source_t)

    ! If necessary, remap vertical coordinates
    if (remap .and. modulo (istep, iremap) == 0) call remap_vertical_coordinates
    
    min_mass = cpt_min_mass ()
    
    ! Add diffusion
    if (modulo (istep_cumul, n_diffuse) == 0 .and. Laplace_order_init /= 0) then
       Laplace_order = Laplace_order_init
    else
       Laplace_order = 0
    end if

    ! Adapt grid
    if (min_level /= max_level) call adapt_grid (set_thresholds)

    call update

    ! Apply velocity penalization (no slip boundary condition)
    if (penalize) call apply_penal (sol)
    
    call sum_total_mass (.false.)

    itime = itime + idt
    
    if (match_time) then
       time  = itime/time_mult
    else
       time = time + dt
    end if
    
    ! Set new time step, find change in vertical levels and count active nodes
    dt_new = cpt_dt ()
  end subroutine time_step

  subroutine restart (run_id)
    ! Fresh restart from checkpoint data (all structures reset)
    implicit none
    character(*) :: run_id

    integer        :: ierror, l    
    character(255) :: cmd_archive, cmd_files, command

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

    ! Uncompress checkpoint data
    if (rank == 0) then
       write (cmd_archive, '(A,I4.4,A)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (6,'(A,A,/)') 'Loading file ', trim (cmd_archive)
       write (command, '(A,A)') 'tar xzf ', trim (cmd_archive)
       call system (command)
    end if
    call barrier ! Make sure all archive files have been uncompressed

    ! Rebalance adaptive grid and re-initialize structures
    call init_structures (run_id)

    ! Initialize thresholds to default values 
    call initialize_thresholds

    ! Load checkpoint data
    call load_adapt_mpi (cp_idx, load, run_id)

    ! Delete temporary files
    call barrier ! Do not delete files before everyone has read them
    if (rank == 0) then
       write (cmd_files, '(A,A,I4.4,A,A,A,I4.4)') &
            trim (run_id), '{_grid,_coef}.', cp_idx , '_????? ', trim (run_id), '_conn.', cp_idx
       write (command, '(A,A)') '\rm ', trim (cmd_files)
       call system (command)
    end if
    call barrier

    call adapt (set_thresholds, .false.) ! Do not re-calculate thresholds, compute masks based on active wavelets
    call inverse_wavelet_transform (wav_coeff, sol, level_start-1)

    ! Apply velocity penalization
    if (penalize) call apply_penal (sol)

    ! Initialize time step and viscosities
    call initialize_dt_viscosity

    ! Set penalization and depth
    call update

    ! Initialize total mass value
    call sum_total_mass (.true.)
    
    ! Initialize time step and counters
    dt_new = cpt_dt ()
    itime = nint (time*time_mult, 8)
    istep = 0

    if (rank == 0) then
       write (6,'(/,A,es12.6,3(A,es8.2),A,I2,A,I9,/)') &
            'time [d] = ', time/DAY, &
            '  mass threshold = ', sum (threshold(S_MASS,:))/zlevels, &
             ' temp threshold = ', sum (threshold(S_TEMP,:))/zlevels, &
             ' velo threshold = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax = ', level_end, &
            '  dof = ', sum (n_active)
       write (6,'(A)') &
            '********************************************************** End Restart &
            ***********************************************************'
    end if
  end subroutine restart

  subroutine write_checkpoint (run_id, rebal)
    implicit none
    character(*) :: run_id
    logical      :: rebal

    character(255) :: cmd_archive, cmd_files, command
    
    cp_idx = cp_idx + 1

    if (rank == 0) then
       write (6,'(A,/)') &
            '************************************************************************&
            **********************************************************'
       write (6,'(a,i4,a,es10.4,/)') 'Saving checkpoint ', cp_idx, ' at time [day] = ', time/DAY
    end if
    
    call write_load_conn (cp_idx, run_id)
    call dump_adapt_mpi  (cp_idx, dump, run_id)
    
    ! Archive checkpoint (overwriting existing checkpoint if present)
    call barrier ! Make sure all processors have written data
    if (rank == 0) then
       write (cmd_files, '(A,A,I4.4,A,A,A,I4.4)') &
            trim (run_id), '{_grid,_coef}.', cp_idx , '_????? ', trim (run_id), '_conn.', cp_idx
       write (cmd_archive, '(A,I4.4,A)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (command, '(A,A,A,A,A)') 'tar cfz ', trim (cmd_archive), ' ', trim (cmd_files), ' --remove-files'
       call system (command)
    end if
    call barrier ! Make sure data is archived before restarting
    
    ! Must restart if want to load balance (compiled with mpi-lb)
    if (rebal) call restart (run_id)
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

    if (rank == 0) write (6,'(A,i2,A)') 'Make level J_min = ', min_level, ' ...'
    
    call init_wavelets
    call init_masks
    call add_second_level

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, +BDRY_THICKNESS)
    call apply_interscale (mask_adj_children, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)
    if (time_end > 0.0_8) time_mult = huge (itime)/2/time_end

    allocate (n_patch_old(size(grid)), n_node_old(size(grid))); n_patch_old = 2

    call init_RK_mem
  end subroutine init_structures

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

       do k = 1, zmax
          deallocate (penal_node(k)%data(d)%elts)
          deallocate (penal_edge(k)%data(d)%elts)
          deallocate (exner_fun(k)%data(d)%elts)
       end do
       deallocate (exner_fun(zmax+1)%data(d)%elts)

       do v = 1, N_VARIABLE
          do k = 1, zmax
             deallocate (sol(v,k)%data(d)%elts)
             deallocate (sol_mean(v,k)%data(d)%elts)
             deallocate (trend(v,k)%data(d)%elts)
             deallocate (wav_coeff(v,k)%data(d)%elts)
             deallocate (trend_wav_coeff(v,k)%data(d)%elts)
          end do
          do k = 1, save_levels
             deallocate (sol_save(v,k)%data(d)%elts) 
          end do
       end do
    end do
    
    deallocate (Laplacian_vector(S_DIVU)%data)
    deallocate (Laplacian_vector(S_ROTU)%data)

    do k = 1, zmax
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
       do k = 1, zmax
          deallocate (sol(v,k)%data)
          deallocate (sol_mean(v,k)%data)
          deallocate (trend(v,k)%data)
          deallocate (wav_coeff(v,k)%data)
          deallocate (trend_wav_coeff(v,k)%data)
       end do
       do k = 1, save_levels
          deallocate (sol_save(v,k)%data)
       end do
    end do

    deallocate (grid, n_patch_old, n_node_old)
    deallocate (edge_level_start, node_level_start, n_active_edges, n_active_nodes)
    deallocate (a_vert, b_vert, a_vert_mass, b_vert_mass)
    deallocate (threshold, threshold_def)
    deallocate (sol, sol_mean, sol_save, trend, trend_wav_coeff, wav_coeff)       
    deallocate (exner_fun, horiz_flux, Laplacian_scalar, Laplacian_vector, lnorm, penal_node, penal_edge)
    deallocate (glo_id, ini_st, recv_lengths, recv_offsets, req, send_lengths, send_offsets)

    nullify (mass, dscalar, h_flux, velo, dvelo, bernoulli, divu, exner, ke, qe, scalar, temp, vort, wc_s, wc_u)
  end subroutine deallocate_structures
end module main_mod
