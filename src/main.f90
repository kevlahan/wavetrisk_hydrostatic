module main_mod
  use param_mod
  use domain_mod
  use comm_mod
  use comm_mpi_mod
  use init_mod
  use refine_patch_mod
  use arch_mod
  use time_integr_mod
  use io_mod
  use wavelet_mod
  use mask_mod
  use adapt_mod
  use smooth_mod

  implicit none
  integer(8)                         :: itime
  integer                            :: cp_idx
  integer, dimension(:), allocatable :: node_level_start, edge_level_start
  real(8)                            :: dt, dt_new, time_mult  
  type Initial_State
     integer                                          :: n_patch, n_bdry_patch, n_node, n_edge, n_tria
     integer, dimension(AT_NODE:AT_EDGE,N_GLO_DOMAIN) :: pack_len, unpk_len
  end type Initial_State
  
  type(Initial_State), dimension(:), allocatable :: ini_st
contains
  subroutine init_main_mod
    implicit none
    call init_arch_mod
    call init_domain_mod
    call init_comm_mod
    call init_comm_mpi_mod
    call init_init_mod
    call init_refine_patch_mod
    call init_time_integr_mod
    call init_io_mod
    call init_wavelet_mod
    call init_mask_mod
    call init_adapt_mod
    time_mult = 1.0_8
  end subroutine init_main_mod

  subroutine initialize (apply_init_cond, set_thresholds, custom_dump, custom_load, test_case)
    implicit none
    real(8), dimension(:), pointer :: wc_m, wc_t, wc_u
    external     :: apply_init_cond, set_thresholds, custom_dump, custom_load
    character(*) :: test_case
    
    character(255) :: command
    integer        :: k, d, ierr

    if (min_level > max_level) then
       if (rank == 0) write(6,'(A,I4,1X,A,I4,A,I4)') 'ERROR: max_level < min_level:', max_level, &
            '<', min_level, '. Setting max_level to', min_level
       max_level = min_level
    end if

    if (resume >= 0) then
       if (rank == 0) write(6,'(A,i6)') 'Resuming from checkpoint ', resume
       cp_idx = resume
       write(command, '(A,I4.4,A)')  'tar xzf '//trim(test_case)//'_checkpoint_' , cp_idx , ".tgz"
       if (rank == 0) call system (command)
       call barrier ! make sure all files are extracted before everyone starts reading them
    end if

    call init_structures

    if (resume >= 0) then
       if (rank == 0) write(6,'(A,i6)') 'Resuming from checkpoint ', resume
    else
       if (rank == 0) write(6,'(/,A,/)') '------------- Adapting initial grid -------------'

       call apply_init_cond
       call forward_wavelet_transform (sol, wav_coeff)
       call trend_ml (sol, trend)
       call forward_wavelet_transform (trend, trend_wav_coeff)

       do while (level_end < max_level)
          if (rank == 0) write(6,'(A,i2,A,i2)') 'Initial refinement Level', level_end, ' -> ', level_end+1
          node_level_start = grid(:)%node%length+1
          edge_level_start = grid(:)%midpt%length+1
          
          if (rank == 0) write(6,'(A,i2)') 'Initialize solution on level ', level_end

          dt_new = cpt_dt_mpi()
          call adapt (set_thresholds)

          call apply_init_cond
         
          call forward_wavelet_transform (sol, wav_coeff)
          call trend_ml (sol, trend)
          call forward_wavelet_transform (trend, trend_wav_coeff)

          !--Check whether there are any active nodes at this scale
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
                n_active = n_active + (/ count( abs(wc_m(node_level_start(d):grid(d)%node%length))  >= tol_mass(k) .or. &
                                                abs(wc_t(node_level_start(d):grid(d)%node%length))  >= tol_temp(k)), &
                                         count( abs(wc_u(edge_level_start(d):grid(d)%midpt%length)) >= tol_velo(k)) /)
                nullify (wc_m, wc_t, wc_u)
             end do
          end do
          ! Sum results over all ranks
          n_active(AT_NODE) = sum_int (n_active(AT_NODE)) ; n_active(AT_EDGE) = sum_int(n_active(AT_EDGE))
          
          if (rank == 0) write(6,'(A,i2,1x,2(A,i8,1x),/)') &
               'Level = ', level_end, 'number of active nodes = ',n_active(AT_NODE), 'number of active edges = ', n_active(AT_EDGE)

          if (n_active(AT_NODE) == 0 .and. n_active(AT_EDGE) == 0) exit !--No active nodes at this scale
       end do
       if (rank == 0) write(6,'(/,A,/)') '--------- Finished adapting initial grid ---------'

       dt_new = cpt_dt_mpi()
       call adapt (set_thresholds)
       dt_new = cpt_dt_mpi()
       if (rank==0) write(6,'(A,i8)') 'Initial dof = ', sum(n_active)

       cp_idx = -1
       ierr = write_checkpoint (custom_dump)
    end if

    if (rank == 0) write(6,'(A,i6,/)') 'Restarting from initial checkpoint ', cp_idx
    call restart_full (set_thresholds, custom_load, test_case)
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

    integer(8)           :: idt, ialign
    
    istep = istep+1
    dt = dt_new

    ! Match certain times exactly
    idt    = nint(dt*time_mult, 8)
    ialign = nint(align_time*time_mult, 8)
    if (ialign>0 .and. cp_idx /= resume .and. istep /= 1) then
       aligned = (modulo(itime+idt,ialign) < modulo(itime,ialign))
    else
       resume = NONE ! Set unequal cp_idx => only first step after resume is protected from alignment
       aligned = .False.
    end if
    if (aligned) idt = ialign - modulo(itime,ialign)
    dt = idt/time_mult ! Modify time step

    call RK34_opt (trend_ml, dt)

    if (min_level < max_level) call adapt_grid (set_thresholds)
    dt_new = cpt_dt_mpi() ! Set new time step and count active nodes
    
    itime = itime + idt
    time  = itime/time_mult
  end subroutine time_step

  subroutine reset (init_state)
    implicit none
    type(Initial_State), dimension (:), allocatable :: init_state
    integer                                         :: k, l, d, v, i
    integer, dimension (AT_NODE:AT_EDGE)            :: num

    do d = 1, size(grid)
       grid(d)%lev(min_level+1:max_level)%length = 0
    end do

    level_start = min_level
    level_end = level_start

    do d = 1, size(grid)
       num = (/init_state(d)%n_node, init_state(d)%n_edge/)  
       grid(d)%patch%length       = init_state(d)%n_patch 
       grid(d)%bdry_patch%length  = init_state(d)%n_bdry_patch 
       grid(d)%node%length        = init_state(d)%n_node 
       grid(d)%midpt%length       = init_state(d)%n_edge
       grid(d)%ccentre%length     = init_state(d)%n_tria
       grid(d)%areas%length       = init_state(d)%n_node
       grid(d)%triarea%length     = init_state(d)%n_tria
       grid(d)%pedlen%length      = init_state(d)%n_edge
       grid(d)%len%length         = init_state(d)%n_edge
       grid(d)%coriolis%length    = init_state(d)%n_tria
       grid(d)%overl_areas%length = init_state(d)%n_node
       grid(d)%I_u_wgt%length     = init_state(d)%n_node
       grid(d)%R_F_wgt%length     = init_state(d)%n_node
       grid(d)%mask_n%length      = init_state(d)%n_node
       grid(d)%mask_e%length      = init_state(d)%n_edge
       grid(d)%level%length       = init_state(d)%n_node

       grid(d)%surf_press%length  = init_state(d)%n_node
       grid(d)%press%length       = init_state(d)%n_node
       grid(d)%surf_geopot%length = init_state(d)%n_node
       grid(d)%geopot%length      = init_state(d)%n_node
       grid(d)%u_zonal%length     = init_state(d)%n_node
       grid(d)%v_merid%length     = init_state(d)%n_node
       grid(d)%adj_mass%length    = init_state(d)%n_node
       grid(d)%adj_temp%length    = init_state(d)%n_node
       grid(d)%adj_geopot%length  = init_state(d)%n_node
       grid(d)%bernoulli%length   = init_state(d)%n_node
       grid(d)%divu%length        = init_state(d)%n_node
       grid(d)%vort%length        = init_state(d)%n_tria
       grid(d)%qe%length          = init_state(d)%n_edge
       grid(d)%Laplacian_u%length = init_state(d)%n_edge

       do k = 1, zlevels
          exner_fun(k)%data(d)%length = num(AT_NODE)
          do v = S_MASS, S_VELO
             wav_coeff(v,k)%data(d)%length = num(POSIT(v))
             trend_wav_coeff(v,k)%data(d)%length = num(POSIT(v))
             trend(v,k)%data(d)%length = num(POSIT(v))
             sol(v,k)%data(d)%length = num(POSIT(v))
             dq1(v,k)%data(d)%length = num(POSIT(v))
             q1(v,k)%data(d)%length = num(POSIT(v))
             q2(v,k)%data(d)%length = num(POSIT(v))
             q3(v,k)%data(d)%length = num(POSIT(v))
             q4(v,k)%data(d)%length = num(POSIT(v))
          end do
       end do
       exner_fun(zlevels+1)%data(d)%length = num(AT_NODE)

       do k = 1, save_levels
          do v = S_MASS, S_VELO
             sol_save(v,k)%data(d)%length = num(POSIT(v))
          end do
       end do

       do v = S_MASS, S_TEMP
          horiz_flux(v)%data(d)%length = num(AT_EDGE)
          Laplacian_scalar(v)%data(d)%length = num(AT_NODE)
       end do

       do i = 1, N_GLO_DOMAIN
          grid(d)%send_conn(i)%length = 0
          grid(d)%recv_pa(i)%length = 0
          do v = AT_NODE, AT_EDGE
             grid(d)%pack(v,i)%length = init_state(d)%pack_len(v,i)
             grid(d)%unpk(v,i)%length = init_state(d)%unpk_len(v,i)
          end do
       end do

       grid(d)%send_pa_all%length = 0
       grid(d)%neigh_pa_over_pole%length = level_end*2 + 2

       ! level 2: set patch children to zero
       do i = 2, 5
          grid(d)%patch%elts(i+1)%children = 0
       end do
    end do
  end subroutine reset

  subroutine restart_full (set_thresholds, custom_load, test_case)
    implicit none
    external     :: set_thresholds, custom_load
    character(*) :: test_case
    
    integer, parameter       :: len_cmd_files = 12 + 4 + 12 + 4
    character(len_cmd_files) :: cmd_files
    character(255)           :: cmd_archive, command

    ! Deallocate all dynamic arrays and variables
    call deallocate_structures

    call init_structures

    if (rank == 0) write (6,*) 'Reloading from checkpoint', cp_idx

    call load_adapt_mpi (cp_idx, custom_load)
        
    itime = nint(time*time_mult, 8)
    resume = cp_idx ! to disable alignment for next step

    ! Do not override existing checkpoint archive
    write(cmd_files, '(A,I4.4,A,I4.4)') '{grid,coef}.', cp_idx , '_????? conn.', cp_idx
    write(cmd_archive, '(A,I4.4,A)') trim(test_case)//'_checkpoint_' , cp_idx, ".tgz"
    write(command, '(A,A,A,A,A,A,A,A,A)') 'if [ -e ', trim(cmd_archive), ' ]; then rm -f ', cmd_files, &
         '; else tar c --remove-files -z -f ', trim(cmd_archive), ' ', cmd_files, '; fi'

    call barrier ! do not delete files before everyone has read them

    if (rank == 0) call system (command)
    istep = 0
    call adapt (set_thresholds, .false.) ! Do not re-calculate thresholds, compute masks based on active wavelets
    call inverse_wavelet_transform (wav_coeff,         sol, level_start-1)
    call inverse_wavelet_transform (trend_wav_coeff, trend, level_start-1)
    dt_new = cpt_dt_mpi()
  end subroutine restart_full

  function write_checkpoint (custom_dump)
    implicit none
    integer  :: write_checkpoint
    external :: custom_dump
    
    external :: custom_load

    character (38+4+22+4+6) :: command
    
    cp_idx = cp_idx + 1
    call write_load_conn (cp_idx)
    write_checkpoint = dump_adapt_mpi (cp_idx, custom_dump)
  end function write_checkpoint

  subroutine compress_files (iwrite, test_case)
    implicit none
    integer      :: iwrite
    character(*) :: test_case

    character(4)   :: s_time
    character(130) :: command

    write (s_time, '(i4.4)') iwrite

    command = 'ls -1 '//trim(test_case)//'.1'//s_time//'?? > tmp1'
    
    call system (command)

    command = 'tar czf '//trim(test_case)//'.1'//s_time//'.tgz -T tmp1 --remove-files &'
    call system (command)

    command = 'ls -1 '//trim(test_case)//'.2'//s_time //'?? > tmp2' 
    call system (command)

    command = 'tar czf '//trim(test_case)//'.2'//s_time//'.tgz -T tmp2 --remove-files &'
    call system (command)
  end subroutine compress_files

  subroutine init_structures
    ! Initialize dynamical arrays and structures
    implicit none

    level_start = min_level
    level_end = level_start
    
    call distribute_grid (resume)
    call init_grid
    call init_comm_mpi
    call init_geometry

    if (optimize_grid == XU_GRID) call smooth_Xu (1.0d6*eps())
    if (optimize_grid == HR_GRID) call read_HR_optim_grid

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call precompute_geometry

    allocate (node_level_start(size(grid)), edge_level_start(size(grid)))

    if (rank == 0) write(6,'(A,i2,A,/)') 'Make level J_min = ', min_level, ' ...'

    call init_wavelets
    call init_masks
    call add_second_level

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, +BDRY_THICKNESS)
    call apply_interscale (mask_adj_children, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)
    if (time_end > 0.0_8) time_mult = huge(itime)/2/time_end

    call init_RK_mem
  end subroutine init_structures

  subroutine deallocate_structures
    ! Deallocate dynamic arrays and structures for clean restart
    implicit none

    integer :: d, i, k, l, v, r

     ! deallocate init_RK_mem allocations
    do k = 1, zlevels
       do d = 1, n_domain(rank+1)
          do v = S_MASS, S_VELO
             deallocate(q1(v,k)%data(d)%elts)
             deallocate(q2(v,k)%data(d)%elts)
             deallocate(q3(v,k)%data(d)%elts)
             deallocate(q4(v,k)%data(d)%elts)
             deallocate(dq1(v,k)%data(d)%elts)
          end do
       end do
       do v = S_MASS, S_VELO
          deallocate(q1(v,k)%data)
          deallocate(q2(v,k)%data)
          deallocate(q3(v,k)%data)
          deallocate(q4(v,k)%data)
          deallocate(dq1(v,k)%data)
       end do
    end do
    deallocate(q1, q2, q3, q4, dq1)

    deallocate(ini_st)

    ! deallocate mask and geometry allocations
    do d = 1, size(grid)
       deallocate (grid(d)%mask_n%elts)
       deallocate (grid(d)%mask_e%elts)
       deallocate (grid(d)%level%elts)
       deallocate (grid(d)%R_F_wgt%elts)
       deallocate (grid(d)%I_u_wgt%elts)
       deallocate (grid(d)%overl_areas%elts)
       deallocate (grid(d)%surf_press%elts)
       deallocate (grid(d)%press%elts)
       deallocate (grid(d)%surf_geopot%elts)
       deallocate (grid(d)%geopot%elts)
       deallocate (grid(d)%u_zonal%elts)
       deallocate (grid(d)%v_merid%elts)
       deallocate (grid(d)%adj_mass%elts)
       deallocate (grid(d)%adj_temp%elts)
       deallocate (grid(d)%adj_geopot%elts)
       deallocate (grid(d)%vort%elts)
       deallocate (grid(d)%qe%elts)
       deallocate (grid(d)%Laplacian_u%elts)
       deallocate (grid(d)%bernoulli%elts)
       deallocate (grid(d)%divu%elts)
       deallocate (grid(d)%coriolis%elts)
       deallocate (grid(d)%triarea%elts)
       deallocate (grid(d)%len%elts)
       deallocate (grid(d)%pedlen%elts)
       deallocate (grid(d)%areas%elts)
       deallocate (grid(d)%midpt%elts)
       deallocate (grid(d)%ccentre%elts)
    end do

    ! deallocate wavelet allocations
    do k = 1, zlevels
       do d = 1, size(grid)
          do v = S_MASS, S_VELO
             deallocate (wav_coeff(v,k)%data(d)%elts)
             deallocate (trend_wav_coeff(v,k)%data(d)%elts)
          end do
       end do
       do v = S_MASS, S_VELO
          deallocate (wav_coeff(v,k)%data)
          deallocate (trend_wav_coeff(v,k)%data)
       end do
    end do
    deallocate (wav_coeff, trend_wav_coeff)

    deallocate (node_level_start, edge_level_start)

    ! deallocate precompute_geometry allocations
    do k = 1, zlevels
       do d = 1, size(grid)
          deallocate (exner_fun(k)%data(d)%elts)
          do v = S_MASS, S_VELO
             deallocate (trend(v,k)%data(d)%elts)
          end do
       end do
    end do
    do d = 1, size(grid)
       deallocate (exner_fun(zlevels+1)%data(d)%elts)
    end do

    do d = 1, size(grid)
       do v = S_MASS, S_TEMP
          deallocate (horiz_flux(v)%data(d)%elts)
          deallocate (Laplacian_scalar(v)%data(d)%elts)
       end do
    end do

    ! deallocate solution arrays
    do k = 1, zlevels
       do d = 1, size(grid)
          do v = S_MASS, S_VELO
             deallocate (sol(v,k)%data(d)%elts) 
          end do
       end do
    end do

    do k = 1, save_levels
       do d = 1, size(grid)
          do v = S_MASS, S_VELO
             deallocate (sol_save(v,k)%data(d)%elts) 
          end do
       end do
    end do

    ! deallocate init_grid allocations
    do d = 1, size(grid)
       deallocate (grid(d)%neigh_pa_over_pole%elts)

       do k = AT_NODE, AT_EDGE
          do i = 1, N_GLO_DOMAIN
             deallocate (grid(d)%pack(k,i)%elts)
             deallocate (grid(d)%unpk(k,i)%elts)
          end do
       end do

       do i = 1, N_GLO_DOMAIN
          deallocate (grid(d)%recv_pa(i)%elts)
       end do

       deallocate (grid(d)%send_pa_all%elts)

       do i = 1, N_GLO_DOMAIN
          deallocate (grid(d)%send_conn(i)%elts)
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

       deallocate(grid(d)%src_patch)
       deallocate(grid(d)%node%elts) 
       deallocate(grid(d)%bdry_patch%elts) 
       deallocate(grid(d)%patch%elts) 
    end do

    do k = 1, zlevels
       deallocate (exner_fun(k)%data)
       do v = S_MASS, S_VELO
          deallocate (sol(v,k)%data)
          deallocate (trend(v,k)%data)
       end do
    end do
    deallocate (exner_fun(zlevels+1)%data)

    do k = 1, save_levels
       do v = S_MASS, S_VELO
          deallocate (sol_save(v,k)%data)
       end do
    end do

    do v = S_MASS, S_TEMP
       deallocate (horiz_flux(v)%data)
       deallocate (Laplacian_scalar(v)%data)
    end do

    deallocate (grid, sol, sol_save, trend, exner_fun, horiz_flux, Laplacian_scalar)
    deallocate (n_active_edges, n_active_nodes)
    deallocate (send_lengths, send_offsets, recv_lengths, recv_offsets, req, stat_ray)
  end subroutine deallocate_structures
end module main_mod
