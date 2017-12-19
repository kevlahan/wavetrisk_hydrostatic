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
  use remap_mod

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

  subroutine initialize (apply_init_cond, stage, set_thresholds, custom_dump, custom_load)
    external :: apply_init_cond, set_thresholds, custom_dump, custom_load
    
    character(20+4+4) :: command
    integer           :: k, d, stage, ierr

    if (min_level .gt. max_level) then
       if (rank .eq. 0) write(*,'(A,I4,1X,A,I4,A,I4)') 'ERROR: max_level < min_level:', max_level, &
            '<', min_level, '. Setting max_level to', min_level
       max_level = min_level
    end if

    if (resume .ge. 0) then
       cp_idx = resume
       write(command, '(A,I4.4,A)')  "tar -xzf checkpoint_" , cp_idx , ".tgz"
       if (rank .eq. 0) call system (command)
       call barrier ! make sure all files are extracted before everyone starts reading them
    end if

    call distribute_grid (resume)
    call init_grid
    call init_comm_mpi
    call init_geometry

    if (optimize_grid .eq. XU_GRID) call smooth_Xu (1d6*eps())
    if (optimize_grid .eq. HR_GRID) call read_HR_optim_grid

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call precompute_geometry

    allocate (node_level_start(size(grid)), edge_level_start(size(grid)))

    if (rank .eq. 0) write(6,*) 'Make level J_min =', min_level, '...'

    call init_wavelets
    call init_masks
    call add_second_level

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, +BDRY_THICKNESS)
    call apply_interscale (mask_adj_scale, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)
    if (time_end .gt. 0.0_8) time_mult = huge(itime)/2/time_end
    if (stage .eq. 0) return

    call init_RK_mem

    if (resume .ge. 0) then
       if (rank .eq. 0) write(6,*) 'Resuming from checkpoint', resume
    else
       call apply_init_cond
       call forward_wavelet_transform (sol, wav_coeff)
       if (adapt_trend) then
          call trend_ml (sol, trend, 0)
          call forward_wavelet_transform (trend, trend_wav_coeff)
       end if

       dt_init = cpt_dt_mpi()
       do while (level_end .lt. max_level)
          if (rank .eq. 0) write(*,*) 'Initial refine. Level', level_end, ' -> ', level_end+1
          node_level_start = grid(:)%node%length+1
          edge_level_start = grid(:)%midpt%length+1
          
          call adapt (set_thresholds)

          if (rank .eq. 0) write(6,*) 'Initialize solution on level', level_end

          call apply_init_cond
          
          call forward_wavelet_transform (sol, wav_coeff)

          if (adapt_trend) then
             call trend_ml (sol, trend, 0)
             call forward_wavelet_transform (trend, trend_wav_coeff)
          end if
                    
          !--Check whether there are any active nodes at this scale
          n_active = 0
          do k = 1, zlevels
             n_active = n_active + (/ &
                  sum((/(count( &
                  abs(wav_coeff(S_MASS,k)%data(d)%elts(node_level_start(d):grid(d)%node%length)) .gt. tol_mass &
                  .or. &
                  abs(wav_coeff(S_TEMP,k)%data(d)%elts(node_level_start(d):grid(d)%node%length)) .gt. tol_temp ), &
                  d = 1, n_domain(rank+1)) /)), &
                  
                  sum((/(count(abs(wav_coeff(S_VELO,k)%data(d)%elts(edge_level_start(d):grid(d)%midpt%length)) .gt. tol_velo), &
                  d = 1, n_domain(rank+1)) /)) &
                  /)
          end do
          n_active(AT_NODE) = sync_max(n_active(AT_NODE))
          n_active(AT_EDGE) = sync_max(n_active(AT_EDGE))
          if (rank.eq.0) write(6,'(A,i2,1x,2(A,i8,1x))') &
               'Level = ', level_end, 'number of active nodes = ',n_active(AT_NODE), 'number of active edges = ', n_active(AT_EDGE)

          if (n_active(AT_NODE) .eq. 0 .and. n_active(AT_EDGE) .eq. 0) exit !--No active nodes at this scale
       end do

       cp_idx = 0
       call write_load_conn(0)
       ierr = dump_adapt_mpi(write_mt_wc, write_u_wc, cp_idx, custom_dump)
    end if
    dt_init = cpt_dt_mpi()
    call restart_full (set_thresholds, custom_load)
  end subroutine initialize

  subroutine record_init_state (init_state)
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
    real(8)              :: align_time
    logical, intent(out) :: aligned
    external             :: set_thresholds

    integer(8)           :: idt, ialign
    
    istep = istep+1
    dt = dt_new

    ! Match certain times exactly
    idt    = nint(dt*time_mult, 8)
    ialign = nint(align_time*time_mult, 8)
    if (ialign.gt.0 .and. cp_idx.ne.resume .and. istep.ne.1) then
       aligned = (modulo(itime+idt,ialign) .lt. modulo(itime,ialign))
    else
       resume = NONE ! Set unequal cp_idx => only first step after resume is protected from alignment
       aligned = .False.
    end if
    if (aligned) idt = ialign - modulo(itime,ialign)
    dt = idt/time_mult ! Modify time step

    call RK34_opt (trend_ml, dt)
    !call RK45_opt (trend_ml, dt)
    !call euler (trend_ml, dt)

    if (min_level .lt. max_level) call adapt_grid (set_thresholds)
    dt_new = cpt_dt_mpi() ! Set new time step and count active nodes
    
    itime = itime + idt
    time  = itime/time_mult
  end subroutine time_step

  subroutine time_step_diffuse (dt)
    ! Euler time step to diffuse solution
    real(8) :: dt
    call euler (trend_ml, dt)
  end subroutine time_step_diffuse

  subroutine reset (init_state)
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
       grid(d)%bern_slow%length   = init_state(d)%n_node
       grid(d)%bern_fast%length   = init_state(d)%n_node
       grid(d)%divu%length        = init_state(d)%n_node
       grid(d)%vort%length        = init_state(d)%n_tria
       grid(d)%qe%length          = init_state(d)%n_edge

       do k = 1, zlevels
          bernoulli_fast(k)%data(d)%length = num(AT_NODE)
          exner_fun(k)%data(d)%length      = num(AT_NODE)
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

       do v = S_MASS, S_TEMP
          horiz_flux(v)%data(d)%length = num(AT_EDGE)
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

  subroutine restart_full (set_thresholds, custom_load)
    external           :: set_thresholds, custom_load
    
    integer            :: i, d, k, r, l, v
    integer, parameter :: len_cmd_files = 12 + 4 + 12 + 4
    integer, parameter :: len_cmd_archive = 11 + 4 + 4
    character(len_cmd_files) :: cmd_files
    character(len_cmd_archive) :: cmd_archive
    character(8+len_cmd_archive+15+len_cmd_files+34+len_cmd_archive+1+len_cmd_files+4) :: command

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
       deallocate (grid(d)%bernoulli%elts)
       deallocate (grid(d)%bern_slow%elts)
       deallocate (grid(d)%bern_fast%elts)
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
          deallocate (bernoulli_fast(k)%data(d)%elts)
          deallocate (exner_fun(k)%data(d)%elts)
          do v = S_MASS, S_VELO
             deallocate (trend(v,k)%data(d)%elts)
          end do
       end do
    end do

    do d = 1, size(grid)
       do v = S_MASS, S_TEMP
          deallocate (horiz_flux(v)%data(d)%elts)
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
       deallocate (bernoulli_fast(k)%data)
       deallocate (exner_fun(k)%data)
       do v = S_MASS, S_VELO
          deallocate (sol(v,k)%data)
          deallocate (trend(v,k)%data)
       end do
    end do

    do v = S_MASS, S_TEMP
       deallocate (horiz_flux(v)%data)
    end do

    deallocate (grid, sol, trend, bernoulli_fast, exner_fun, horiz_flux)

    ! init_shared_mod
    level_start = min_level
    level_end = level_start

    call distribute_grid (cp_idx)
    deallocate (n_active_edges, n_active_nodes)

    call init_grid 
    call init_comm 
    call comm_communication_mpi 
    call init_geometry 

    if (optimize_grid .eq. XU_GRID) call smooth_Xu (1d6*eps())
    if (optimize_grid .eq. HR_GRID) call read_HR_optim_grid 

    call comm_nodes3_mpi (get_coord, set_coord, NONE)
    call precompute_geometry 

    allocate (node_level_start(size(grid)), edge_level_start(size(grid)))

    if (rank .eq. 0) write(*,*) 'Make level J_min =', min_level, '...'
    call init_wavelets 
    call init_masks 
    call add_second_level 

    call apply_onescale2 (set_level, level_start, z_null, -BDRY_THICKNESS, +BDRY_THICKNESS)
    call apply_interscale (mask_adj_scale, level_start-1, z_null, 0, 1) ! level 0 = TOLRNZ => level 1 = ADJZONE

    call record_init_state (ini_st)

    call init_RK_mem

    if (rank .eq. 0) write (6,*) 'Reloading from checkpoint', cp_idx

    call load_adapt_mpi (read_mt_wc_and_mask, read_u_wc_and_mask, cp_idx, custom_load)
        
    itime = nint(time*time_mult, 8)

    ! Do not override existing checkpoint archive
    write(cmd_files, '(A,I4.4,A,I4.4)') "{grid,coef}.", cp_idx , "_????? conn.", cp_idx
    write(cmd_archive, '(A,I4.4,A)') "checkpoint_" , cp_idx, ".tgz"
    write(command, '(A,A,A,A,A,A,A,A,A)') "if [ -e ", cmd_archive, " ]; then rm -f ", cmd_files, &
         "; else tar c --remove-files -z -f ", cmd_archive, " ", cmd_files, "; fi"

    call barrier ! do not delete files before everyone has read them

    if (rank .eq. 0) call system (command)

    call inverse_wavelet_transform (wav_coeff, sol, level_start-1)
    istep = 0
    dt_init = cpt_dt_mpi()
    if (adapt_trend) call trend_ml (sol, trend, 0)
    call adapt (set_thresholds)
    call inverse_wavelet_transform (wav_coeff, sol, level_start)
    dt_new = cpt_dt_mpi()
  end subroutine restart_full

  function write_checkpoint (custom_dump)
    integer  :: write_checkpoint
    external :: custom_dump
    
    external :: custom_load

    character (38+4+22+4+6) :: command
    
    cp_idx = cp_idx + 1
    call write_load_conn (cp_idx)
    write_checkpoint = dump_adapt_mpi (write_mt_wc, write_u_wc, cp_idx, custom_dump)
  end function write_checkpoint

  subroutine compress_files (iwrite)
    integer :: iwrite

    character(3)   :: s_time
    character(130) :: command

    write (s_time, '(i3.3)') iwrite

    command = 'ls -1 fort.1'//s_time // '* > tmp' 
    call system (command)

    command = 'tar czf fort.1'//s_time //'.tgz -T tmp --remove-files'
    call system (command)

    command = '\rm tmp; ls -1 fort.2' //s_time // '* > tmp' 
    call system (command)

    command = 'tar czf fort.2'//s_time //'.tgz -T tmp --remove-files;\rm tmp'
    call system (command)
  end subroutine compress_files
end module main_mod
