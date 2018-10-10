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

  subroutine initialize (apply_init_cond, set_thresholds, custom_dump, custom_load, run_id)
    ! Initialize from checkpoint or adapt to initialize conditions
    ! Solution is saved and restarted to balance load
    implicit none
    external     :: apply_init_cond, set_thresholds, custom_dump, custom_load
    character(*) :: run_id
    
    character(255)                 :: command
    integer                        :: k, d
    real(8), dimension(:), pointer :: wc_m, wc_t, wc_u

    if (resume >= 0) then
       cp_idx = resume
       call restart (set_thresholds, custom_load, run_id)
       resume = NONE
    else
       ! Initialize basic structures
       call init_basic

       ! Initialize vertical grid
       call initialize_a_b_vert

       ! Determine vertical level to save
       call set_save_level

       ! Initialize thresholds to default values 
       call initialize_thresholds

       call init_structures
       call apply_init_cond

       ! Calculate diffusion length scales
       if (Laplace_order /= 0) call evals_diffusion

       ! Initialize time step and viscosities
       call initialize_dt_viscosity

       if (rank == 0) write (6,'(/,A,/)') &
            '----------------------------------------------------- Adapting initial grid &
            -----------------------------------------------------'

       call forward_wavelet_transform (sol, wav_coeff)
       call trend_ml (sol, trend)
       call forward_wavelet_transform (trend, trend_wav_coeff)

       do while (level_end < max_level)
          if (rank == 0) write (6,'(A,i2,A,i2)') 'Initial refinement Level', level_end, ' -> ', level_end+1
          node_level_start = grid(:)%node%length+1
          edge_level_start = grid(:)%midpt%length+1
          
          dt_new = cpt_dt_mpi()
          call adapt (set_thresholds)

          call apply_init_cond
         
          call forward_wavelet_transform (sol, wav_coeff)
          call trend_ml (sol, trend)
          call forward_wavelet_transform (trend, trend_wav_coeff)

          ! Check whether there are any active nodes or edges at this scale
          n_active = 0
          do k = 1, zlevels
             do d = 1, size(grid)
                if (adapt_trend) then
                   wc_m => trend_wav_coeff(S_MASS,k)%data(d)%elts
                   wc_t => trend_wav_coeff(S_TEMP,k)%data(d)%elts
                   wc_u => trend_wav_coeff(S_VELO,k)%data(d)%elts
                else
                   wc_m => wav_coeff(S_MASS,k)%data(d)%elts
                   wc_t => wav_coeff(S_TEMP,k)%data(d)%elts
                   wc_u => wav_coeff(S_VELO,k)%data(d)%elts
                end if
                n_active = n_active + (/ count(abs(wc_m(node_level_start(d):grid(d)%node%length))  >= threshold(S_MASS,k) .or. &
                                               abs(wc_t(node_level_start(d):grid(d)%node%length))  >= threshold(S_TEMP,k)), &
                                         count(abs(wc_u(edge_level_start(d):grid(d)%midpt%length)) >= threshold(S_VELO,k)) /)
                nullify (wc_m, wc_t, wc_u)
             end do
          end do
          
          ! Sum results over all ranks
          n_active(AT_NODE) = sum_int (n_active(AT_NODE)) ; n_active(AT_EDGE) = sum_int(n_active(AT_EDGE))          
          if (rank == 0) write (6,'(A,i2,1x,2(A,i8,1x),/)') &
               'Level = ', level_end, 'number of active node wavelets = ', n_active(AT_NODE), &
               'number of active edge wavelets = ', n_active(AT_EDGE)
          if (n_active(AT_NODE) == 0 .and. n_active(AT_EDGE) == 0) exit ! No active nodes or edges at this scale
       end do
       if (rank == 0) write (6,'(A,/)') &
            '------------------------------------------------- Finished adapting initial grid &
            ------------------------------------------------'

       call adapt (set_thresholds) ; dt_new = cpt_dt_mpi()
       if (rank==0) write (6,'(A,i8,/)') 'Initial number of dof = ', sum (n_active)
       
       call write_checkpoint (custom_dump, custom_load, run_id)
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

  subroutine time_step (align_time, aligned, set_thresholds)
    implicit none
    real(8)              :: align_time
    logical, intent(out) :: aligned
    external             :: set_thresholds

    integer(8) :: idt, ialign
    
    istep = istep+1
    dt = dt_new

    ! Match certain times exactly
    idt    = nint (dt*time_mult, 8)
    ialign = nint (align_time*time_mult, 8)
    if (ialign > 0 .and. istep /= 1) then
       aligned = (modulo (itime+idt,ialign) < modulo (itime,ialign))
    else
       aligned = .false.
    end if
    if (aligned) idt = ialign - modulo (itime,ialign)
    dt = idt/time_mult ! Modify time step

    ! Add diffusion
    if (modulo (nint(time/dt_cfl), n_diffuse) == 0 .and. Laplace_order_init /= 0) then
       Laplace_order = Laplace_order_init
    else
       Laplace_order = 0
    end if
        
    ! Take time step
    sol%bdry_uptodate = .false.
    call update_array_bdry (sol, NONE)
    call RK4 (trend_ml, dt)
    
    ! If necessary, remap vertical coordinates
    if (remap .and. min_mass < min_allowed_mass) call remap_vertical_coordinates (set_thresholds)

    ! Adapt grid
    if (min_level /= max_level) call adapt_grid (set_thresholds)
    
    ! Set new time step, find change in vertical levels and count active nodes
    dt_new = cpt_dt_mpi() 

    itime = itime + idt
    time  = itime/time_mult
  end subroutine time_step

  subroutine restart (set_thresholds, custom_load, run_id)
    ! Fresh restart from checkpoint data (all structures reset)
    implicit none
    external     :: set_thresholds, custom_load
    character(*) :: run_id

    integer :: ierror    
    character(255) :: cmd_archive, cmd_files, command

    if (rank == 0) then
       write (6,'(A,/)') &
            '********************************************************* Begin Restart &
            *********************************************************'
       write (6,'(A,i4,/)') 'Reloading from checkpoint ', cp_idx
    end if

    ! Deallocate all dynamic arrays and variables
    if (resume == NONE) then
       call deallocate_structures
       allocate (n_domain(n_process))
       call init (send_buf_i, 0)
       call init (recv_buf_i, 0) 
       call init (send_buf,   0)
       call init (recv_buf,   0)
    end if

    ! Initialize basic structures
    call init_basic

    ! Initialize vertical grid
    call initialize_a_b_vert

    ! Determine vertical level to save
    call set_save_level

    ! Initialize thresholds to default values 
    call initialize_thresholds

    ! Uncompress checkpoint data
    if (rank == 0) then
       write (cmd_archive, '(A,I4.4,A)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (6,'(A,A,/)') 'Loading file ', trim (cmd_archive)
       write (command, '(A,A)') 'tar xzf ', trim (cmd_archive)
       call system (command)
    end if
    call barrier ! Make sure all archive files have been uncompressed

    ! Rebalance adaptive grid and re-initialize structures
    call init_structures
    
    ! Load checkpoint data
    call load_adapt_mpi (cp_idx, custom_load)

    ! Delete temporary files
    call barrier ! Do not delete files before everyone has read them
    if (rank == 0) then
       write (cmd_files, '(A,I4.4,A,I4.4)') '{grid,coef}.', cp_idx , '_????? conn.', cp_idx
       write (command, '(A,A)') '\rm ', trim (cmd_files)
       call system (command)
    end if
    call barrier

    call adapt (set_thresholds, .false.) ! Do not re-calculate thresholds, compute masks based on active wavelets
    call inverse_wavelet_transform (wav_coeff, sol, level_start-1)
      
    ! Initialize total mass value
    call sum_total_mass (.true.)

    ! Calculate diffusion length scales
    if (Laplace_order /= 0 .and. resume /= NONE) call evals_diffusion

    ! Initialize time step and viscosities
    call initialize_dt_viscosity

    ! Initialize time step and counters
    dt_new = cpt_dt_mpi()
    itime = nint (time*time_mult, 8)
    istep = 0

    if (rank == 0) then
       write (6,'(/,A,es12.6,3(A,es8.2),A,I2,A,I9,/)') &
            'time [h] = ', time/HOUR, &
            '  mass threshold = ', sum (threshold(S_MASS,:))/zlevels, &
             ' temp threshold = ', sum (threshold(S_TEMP,:))/zlevels, &
             ' velo threshold = ', sum (threshold(S_VELO,:))/zlevels, &
            ' Jmax =', level_end, &
            '  dof = ', sum (n_active)
       write (6,'(A)') &
            '********************************************************** End Restart &
            **********************************************************'
    end if
  end subroutine restart

  subroutine write_checkpoint (custom_dump, custom_load, run_id)
    implicit none
    external :: custom_dump, custom_load
    character(*) :: run_id

    character(255) :: cmd_archive, cmd_files, command
    
    cp_idx = cp_idx + 1

    if (rank == 0) then
       write (6,'(A,/)') &
            '***********************************************************************&
            **********************************************************'
       write (6,'(A,i4,A,es10.4,/)') 'Saving checkpoint ', cp_idx, ' at time [h] = ', time/HOUR
    end if
    
    call write_load_conn (cp_idx)
    call dump_adapt_mpi (cp_idx, custom_dump)
    
    ! Archive checkpoint (overwriting existing checkpoint if present)
    call barrier ! Make sure all processors have written data
    if (rank == 0) then
       write (cmd_files, '(A,I4.4,A,I4.4)') '{grid,coef}.', cp_idx , '_????? conn.', cp_idx
       write (cmd_archive, '(A,I4.4,A)') trim (run_id)//'_checkpoint_' , cp_idx, ".tgz"
       write (command, '(A,A,A,A)') 'tar c --remove-files -z -f ', trim (cmd_archive), ' ', trim (cmd_files)
       call system (command)
    end if
    call barrier ! Make sure data is archived before restarting
    
    ! Must restart after checkpoint and load balance (if compiled with mpi-lb)
    call restart (set_thresholds, custom_load, run_id)
  end subroutine write_checkpoint

  subroutine init_structures
    ! Initialize dynamical arrays and structures
    implicit none

    level_start = min_level
    level_end = level_start
    
    ! Distribute and balance grid over processors (necessary for correct restart!)
    call distribute_grid (cp_idx)
    
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

    call init_RK_mem
  end subroutine init_structures

  subroutine deallocate_structures
    ! Deallocate all dynamic arrays and structures for clean restart
    implicit none

    integer :: d, i, k, l, v, r

    ! Deallocate init_RK_mem allocations
    do k = 1, zlevels
       do d = 1, n_domain(rank+1)
          do v = S_MASS, S_VELO
             deallocate (q1(v,k)%data(d)%elts)
             deallocate (q2(v,k)%data(d)%elts)
             deallocate (q3(v,k)%data(d)%elts)
             deallocate (q4(v,k)%data(d)%elts)
             deallocate (dq1(v,k)%data(d)%elts)
          end do
       end do
       do v = S_MASS, S_VELO
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
       deallocate (grid(d)%adj_mass%elts)
       deallocate (grid(d)%adj_temp%elts)
       deallocate (grid(d)%adj_geopot%elts)
       deallocate (grid(d)%vort%elts)
       deallocate (grid(d)%qe%elts)
       deallocate (grid(d)%bernoulli%elts)
       deallocate (grid(d)%divu%elts)
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

       deallocate (Laplacian_divu%data(d)%elts)
       deallocate (Laplacian_rotu%data(d)%elts)
       
       do v = S_MASS, S_TEMP
          deallocate (horiz_flux(v)%data(d)%elts)
          deallocate (Laplacian_scalar(v)%data(d)%elts)
       end do

       do k = 1, zlevels+1
          deallocate (exner_fun(k)%data(d)%elts)
       end do
       
       do v = S_MASS, S_VELO
          do k = 1, zlevels
             deallocate (sol(v,k)%data(d)%elts)
             deallocate (trend(v,k)%data(d)%elts)
             deallocate (wav_coeff(v,k)%data(d)%elts)
             deallocate (trend_wav_coeff(v,k)%data(d)%elts)
          end do
          do k = 1, save_levels
             deallocate (sol_save(v,k)%data(d)%elts) 
          end do
       end do
    end do
    
    deallocate (Laplacian_divu%data)
    deallocate (Laplacian_rotu%data)

    do k = 1, zlevels+1
       deallocate (exner_fun(k)%data)
    end do
    
    do v = S_MASS, S_TEMP
       deallocate (horiz_flux(v)%data)
       deallocate (Laplacian_scalar(v)%data)
    end do
    
    do v = S_MASS, S_VELO
       do k = 1, zlevels
          deallocate (sol(v,k)%data)
          deallocate (trend(v,k)%data)
          deallocate (wav_coeff(v,k)%data)
          deallocate (trend_wav_coeff(v,k)%data)
       end do
       do k = 1, save_levels
          deallocate (sol_save(v,k)%data)
       end do
    end do

    deallocate (recv_buf%elts, send_buf%elts, recv_buf_i%elts, send_buf_i%elts)

    deallocate (grid)
    deallocate (edge_level_start, node_level_start, n_active_edges, n_active_nodes)
    deallocate (a_vert, b_vert, a_vert_mass, b_vert_mass)
    deallocate (viscosity_divu, threshold, threshold_def)
    deallocate (sol, sol_save, trend, trend_wav_coeff, wav_coeff)       
    deallocate (exner_fun, horiz_flux, Laplacian_scalar)
    deallocate (glo_id, ini_st, n_domain, recv_lengths, recv_offsets, req, send_lengths, send_offsets, stat_ray)

    nullify (mass, dmass, h_mflux, temp, dtemp, h_tflux, velo, dvelo, wc_u, wc_m, wc_t, bernoulli, divu, exner, &
         qe, vort, wc_u, wc_m, wc_t)
  end subroutine deallocate_structures
end module main_mod
