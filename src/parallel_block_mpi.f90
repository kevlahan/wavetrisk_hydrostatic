module parallel_block_mpi_mod

  use iso_fortran_env, only : error_unit, int8, int64
  use mpi_f08,        only : MPI_Allgather, MPI_Allgatherv, MPI_Allreduce, &
       MPI_Alltoall, MPI_Alltoallv, MPI_Exscan, MPI_BYTE, MPI_INTEGER, &
       MPI_INTEGER8, MPI_DOUBLE_PRECISION, MPI_MAX, MPI_MIN, MPI_SUCCESS, &
       MPI_SUM

  use kind_mod,   only : dp
  use shared_mod, only : EDGE, N_BDRY, N_CHDRN, N_GLO_DOMAIN, &
       S_MASS, S_TEMP, &
       c_p, compressible, grav_accel, kappa, p_0, p_top, zlevels

  use domain_mod, only : grid, sol, sol_mean, tke, topography, &
       wav_coeff, wav_tke, &
       subtree_weight_Domain

  use patch_mod, only : PATCH_SIZE

  use arch_mod, only : abort_run, block_catalog, comm, glo_id, loc_id, &
       n_process, owner, Parallel_Block, rank

  use parallel_block_mod, only : Block_Data, block_source, block_received, &
       block_source_catalog_index, &
       block_migrating_source_index, block_received_catalog_index, &
       packed_block_nbyte, pack_block, unpack_block, &
       check_block_storage, install_local_blocks, clear_block_staging, &
       clear_local_blocks, local_block_store_ready, n_local_blocks, &
       local_block_catalog, catalog_local_block, &
       get_local_block_identity, check_local_block_storage, &
       get_block_field_layout, get_block_turbulence_layout, &
       get_local_block_field_layout, get_local_block_turbulence_layout, &
       local_block_field_statistics, local_block_wavelet_statistics, &
       local_block_mean_field_statistics, &
       local_block_turbulence_statistics, &
       local_block_topography_statistics, &
       source_block_scalar_stencil_statistics, &
       local_block_scalar_stencil_statistics, &
       source_block_vector_stencil_statistics, &
       local_block_vector_stencil_statistics, &
       source_block_boundary_route_statistics, &
       local_block_boundary_route_statistics, &
       source_block_ghost_source_statistics, &
       local_block_ghost_source_statistics, &
       validate_local_block_ghost_sources, &
       get_local_block_ghost_requests, local_block_patch_count, &
       local_block_boundary_count, &
       local_block_ghost_count, local_block_scalar_patch_nvalue, &
       local_block_scalar_family_patch_nvalue, &
       local_block_scalar_boundary_nvalue, &
       local_block_scalar_family_boundary_nvalue, &
       get_local_block_scalar_patch_values, &
       get_local_block_scalar_patch_family_values, &
       set_local_block_scalar_patch_values, &
       set_local_block_scalar_patch_family_values, &
       fill_local_block_scalar_patch_values, &
       fill_local_block_scalar_patch_family_values, &
       get_local_block_scalar_boundary_values, &
       get_local_block_scalar_boundary_family_values, &
       set_local_block_scalar_boundary_family_values, &
       fill_local_block_scalar_boundary_values, &
       fill_local_block_scalar_boundary_family_values, &
       get_local_block_scalar_ghost_values, &
       get_local_block_scalar_ghost_family_values, &
       set_local_block_scalar_ghost_family_values, &
       fill_local_block_scalar_ghost_family_values, &
       local_block_vector_patch_nvalue, &
       local_block_vector_family_patch_nvalue, &
       local_block_vector_boundary_nvalue, &
       local_block_vector_family_boundary_nvalue, &
       get_local_block_vector_patch_values, &
       get_local_block_vector_patch_family_values, &
       set_local_block_vector_patch_values, &
       set_local_block_vector_patch_family_values, &
       fill_local_block_vector_patch_values, &
       fill_local_block_vector_patch_family_values, &
       get_local_block_vector_boundary_values, &
       get_local_block_vector_boundary_family_values, &
       set_local_block_vector_boundary_family_values, &
       fill_local_block_vector_boundary_values, &
       fill_local_block_vector_boundary_family_values, &
       get_local_block_vector_ghost_values, &
       get_local_block_vector_ghost_family_values, &
       set_local_block_vector_ghost_family_values, &
       fill_local_block_vector_ghost_family_values, &
       ensure_local_block_hydrostatic_state, &
       local_block_hydrostatic_state_ready, &
       local_block_hydrostatic_refresh_count, &
       local_block_hydrostatic_block_refresh_count, &
       local_block_hydrostatic_surface_nvalue, &
       local_block_hydrostatic_column_nvalue, &
       get_local_block_hydrostatic_patch_values, &
       get_local_block_hydrostatic_values, &
       apply_local_block_field_consumer, &
       Local_Block_Tendency_Kernel, &
       apply_local_block_tendency_kernel, &
       apply_local_block_tendency_consumer, &
       local_block_tendency_state_ready, &
       local_block_tendency_execution_count, &
       local_block_tendency_allocation_count, &
       local_block_tendency_statistics, &
       reset_local_block_tendency_accumulator, &
       accumulate_local_block_tendency, &
       begin_local_block_accumulated_tendency_trial, &
       local_block_tendency_accumulator_state_ready, &
       local_block_tendency_accumulator_allocation_count, &
       local_block_tendency_accumulator_statistics, &
       begin_local_block_tendency_trial, &
       commit_local_block_tendency_trial, &
       finalize_local_block_tendency_commit, &
       restore_local_block_tendency_commit, &
       local_block_tendency_commit_checkpoint_is_ready, &
       local_block_tendency_commit_checkpoint_statistics, &
       rollback_local_block_tendency_trial, &
       local_block_tendency_trial_is_active, &
       apply_local_block_hydrostatic_consumer, &
       local_block_hydrostatic_statistics, &
       BLOCK_PAYLOAD_SOL, BLOCK_PAYLOAD_WAV_COEFF, &
       STORE_PATCH, STORE_BDRY, STORE_GHOST

  implicit none

  private

  real(dp), parameter :: BLOCK_GHOST_POISON = &
       -0.25_dp*huge(0.0_dp)
  real(dp), parameter :: BLOCK_BOUNDARY_POISON = &
       -0.125_dp*huge(0.0_dp)
  real(dp), parameter :: BLOCK_PATCH_POISON = &
       -0.0625_dp*huge(0.0_dp)

  type, public :: Block_Migration_Manifest
     integer :: n_send = 0
     integer :: n_recv = 0
     integer, allocatable :: send_count(:)
     integer, allocatable :: recv_count(:)
     integer, allocatable :: send_displ(:)
     integer, allocatable :: recv_displ(:)
     integer, allocatable :: send_block(:)
     integer, allocatable :: recv_block(:)
     integer, allocatable :: send_nbyte(:)
     integer, allocatable :: recv_nbyte(:)
     integer, allocatable :: send_byte_count(:)
     integer, allocatable :: recv_byte_count(:)
     integer, allocatable :: send_byte_displ(:)
     integer, allocatable :: recv_byte_displ(:)
     integer(int8), allocatable :: recv_payload(:)
     integer(int64) :: total_send_nbyte = 0_int64
     integer(int64) :: total_recv_nbyte = 0_int64
     integer :: max_send_nbyte = 0
     integer :: max_recv_nbyte = 0
     logical :: validated = .false.
     logical :: sizes_validated = .false.
     logical :: payload_validated = .false.
  end type Block_Migration_Manifest

  type :: Block_Ghost_Exchange_Plan
     integer :: n_request = 0
     integer :: n_local_request = 0
     integer :: n_remote_send = 0
     integer :: n_remote_recv = 0
     integer, allocatable :: source_block(:)
     integer, allocatable :: source_local_patch(:)
     integer, allocatable :: source_owner(:)
     integer, allocatable :: destination_block(:)
     integer, allocatable :: destination_ghost(:)
     integer, allocatable :: send_record_count(:)
     integer, allocatable :: recv_record_count(:)
     integer, allocatable :: send_record_displ(:)
     integer, allocatable :: recv_record_displ(:)
     integer, allocatable :: send_data(:)
     integer, allocatable :: recv_data(:)
     integer, allocatable :: request_index(:)
     integer :: scalar_n_value = 0
     integer :: vector_n_value = 0
     integer, allocatable :: scalar_send_count(:)
     integer, allocatable :: scalar_recv_count(:)
     integer, allocatable :: scalar_send_displ(:)
     integer, allocatable :: scalar_recv_displ(:)
     integer, allocatable :: vector_send_count(:)
     integer, allocatable :: vector_recv_count(:)
     integer, allocatable :: vector_send_displ(:)
     integer, allocatable :: vector_recv_displ(:)
     real(dp), allocatable :: scalar_send_buffer(:)
     real(dp), allocatable :: scalar_recv_buffer(:)
     real(dp), allocatable :: scalar_patch_buffer(:)
     real(dp), allocatable :: vector_send_buffer(:)
     real(dp), allocatable :: vector_recv_buffer(:)
     real(dp), allocatable :: vector_patch_buffer(:)
     logical :: ready = .false.
  end type Block_Ghost_Exchange_Plan

  type(Block_Ghost_Exchange_Plan), save :: ghost_exchange_plan

  type :: Block_Writeback_Plan_Type
     integer :: catalog_size = 0
     integer :: domain_count = 0
     integer :: installed_block_count = 0
     integer :: n_send = 0
     integer :: n_recv = 0
     integer :: n_retained = 0
     integer, allocatable :: send_count(:)
     integer, allocatable :: recv_count(:)
     integer, allocatable :: send_displ(:)
     integer, allocatable :: recv_displ(:)
     integer, allocatable :: send_block(:)
     integer, allocatable :: recv_block(:)
     integer, allocatable :: send_patch_count(:)
     integer, allocatable :: recv_patch_count(:)
     integer, allocatable :: send_scalar_nvalue(:)
     integer, allocatable :: recv_scalar_nvalue(:)
     integer, allocatable :: send_vector_nvalue(:)
     integer, allocatable :: recv_vector_nvalue(:)
     integer, allocatable :: scalar_send_count(:)
     integer, allocatable :: scalar_recv_count(:)
     integer, allocatable :: scalar_send_displ(:)
     integer, allocatable :: scalar_recv_displ(:)
     integer, allocatable :: vector_send_count(:)
     integer, allocatable :: vector_recv_count(:)
     integer, allocatable :: vector_send_displ(:)
     integer, allocatable :: vector_recv_displ(:)
     real(dp), allocatable :: scalar_send_buffer(:)
     real(dp), allocatable :: scalar_recv_buffer(:)
     real(dp), allocatable :: vector_send_buffer(:)
     real(dp), allocatable :: vector_recv_buffer(:)
     integer, allocatable :: domain_patch_displ(:)
     real(dp), allocatable :: scalar_domain_stage(:)
     real(dp), allocatable :: vector_domain_stage(:)
     logical, allocatable :: domain_patch_covered(:)
     integer :: scalar_patch_nvalue = 0
     integer :: vector_patch_nvalue = 0
     integer :: reconstructed_patch_count = 0
     integer :: preserved_patch_count = 0
     integer(int64) :: buffer_allocations = 0_int64
     integer(int64) :: stage_allocations = 0_int64
     integer(int64) :: production_writeback_count = 0_int64
     logical :: ready = .false.
  end type Block_Writeback_Plan_Type

  type(Block_Writeback_Plan_Type), save :: block_writeback_plan

  type :: Block_Boundary_Snapshot
     integer :: catalog_index = 0
     integer :: boundary_index = 0
     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
  end type Block_Boundary_Snapshot

  type :: Block_Patch_Snapshot
     integer :: catalog_index = 0
     integer :: local_patch = -1
     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
  end type Block_Patch_Snapshot

  type :: Block_Hydrostatic_Traversal_Context
     integer(int64) :: block_count = 0_int64
     integer(int64) :: patch_count = 0_int64
     integer(int64) :: surface_count = 0_int64
     integer(int64) :: column_count = 0_int64
     real(dp) :: surface_moment(3) = 0.0_dp
     real(dp) :: exner_moment(3) = 0.0_dp
     real(dp) :: temperature_moment(3) = 0.0_dp
  end type Block_Hydrostatic_Traversal_Context

  type :: Block_Field_Traversal_Context
     integer(int64) :: block_count = 0_int64
     integer(int64) :: patch_count = 0_int64
     integer(int64) :: boundary_count = 0_int64
     integer(int64) :: ghost_count = 0_int64
     integer(int64) :: node_count = 0_int64
     integer(int64) :: boundary_node_count = 0_int64
     integer(int64) :: ghost_node_count = 0_int64
     integer(int64) :: scalar_count(3) = 0_int64
     integer(int64) :: vector_count(3) = 0_int64
     integer(int64) :: surface_count = 0_int64
     real(dp) :: scalar_moment(3,3) = 0.0_dp
     real(dp) :: vector_moment(3,3) = 0.0_dp
     real(dp) :: surface_moment(3) = 0.0_dp
  end type Block_Field_Traversal_Context

  type :: Block_Stencil_Kernel_Context
     integer(int64) :: block_count = 0_int64
     integer(int64) :: address_count(3) = 0_int64
     integer(int64) :: scalar_count = 0_int64
     integer(int64) :: vector_count = 0_int64
     real(dp) :: scalar_moment(3) = 0.0_dp
     real(dp) :: vector_moment(3) = 0.0_dp
     real(dp) :: scalar_difference_moment(3) = 0.0_dp
     real(dp) :: vector_difference_moment(3) = 0.0_dp
  end type Block_Stencil_Kernel_Context

  type :: Block_Tendency_Traversal_Context
     integer(int64) :: block_count = 0_int64
     integer(int64) :: scalar_count = 0_int64
     integer(int64) :: vector_count = 0_int64
     integer(int64) :: scalar_changed_block_count = 0_int64
     real(dp) :: scalar_moment(3) = 0.0_dp
     real(dp) :: vector_moment(3) = 0.0_dp
  end type Block_Tendency_Traversal_Context

  type, public :: Block_Two_Stage_Step_Result
     integer(int64) :: scalar_count = 0_int64
     integer(int64) :: vector_count = 0_int64
     integer(int64) :: scalar_changed_block_count = 0_int64
     integer(int64) :: vector_changed_block_count = 0_int64
     integer(int64) :: stage_count = 0_int64
     real(dp) :: scalar_moment(3) = 0.0_dp
     real(dp) :: vector_moment(3) = 0.0_dp
     real(dp) :: scalar_max_update = 0.0_dp
     real(dp) :: vector_max_update = 0.0_dp
  end type Block_Two_Stage_Step_Result

  public :: build_block_migration_manifest
  public :: check_block_migration_manifest
  public :: exchange_block_migration_sizes
  public :: exchange_block_migration_payloads
  public :: clear_block_migration_manifest
  public :: build_parallel_block_catalog
  public :: clear_parallel_block_state
  public :: parallel_block_state_is_ready
  public :: invalidate_parallel_block_domain_shadow
  public :: synchronize_parallel_block_checkpoint
  public :: prepare_parallel_block_grid_change
  public :: migrate_blocks
  public :: check_local_blocks
  public :: check_block_field_inventory
  public :: check_block_scalar_stencil_consumer
  public :: check_block_vector_stencil_consumer
  public :: check_block_boundary_routes
  public :: check_block_ghost_source_addresses
  public :: check_block_ghost_request_manifest
  public :: build_block_ghost_exchange_plan
  public :: clear_block_ghost_exchange_plan
  public :: check_block_ghost_exchange_plan
  public :: build_block_writeback_plan
  public :: clear_block_writeback_plan
  public :: check_block_writeback_plan
  public :: exchange_block_writeback_payloads
  public :: write_block_field_family_to_domains
  public :: block_domain_production_writeback_count
  public :: check_block_writeback_payload_exchange
  public :: check_block_writeback_domain_reconstruction
  public :: block_writeback_plan_is_ready
  public :: block_writeback_plan_allocation_count
  public :: check_block_scalar_ghost_payload_exchange
  public :: check_block_vector_ghost_payload_exchange
  public :: refresh_block_sol_ghosts
  public :: apply_refreshed_block_tendency_kernel
  public :: begin_block_two_stage_tendency_step
  public :: complete_block_two_stage_tendency_step
  public :: refresh_block_wav_coeff_ghosts
  public :: refresh_block_sol_wav_coeff_ghosts
  public :: check_production_block_ghost_refresh
  public :: check_block_field_family_accessors
  public :: check_block_patch_writable_storage
  public :: check_block_boundary_family_mutators
  public :: check_block_boundary_family_bulk_fill
  public :: check_refreshed_block_stencil_consumers
  public :: check_block_hydrostatic_state_accessors
  public :: check_block_hydrostatic_consumer
  public :: check_block_field_consumer
  public :: check_block_stencil_kernel
  public :: check_block_tendency_kernel
  public :: check_block_tendency_trial_update
  public :: check_block_tendency_commit
  public :: check_block_tendency_step_driver
  public :: check_block_tendency_accepted_step
  public :: check_block_multistage_tendency_accumulator
  public :: check_block_multistage_tendency_commit
  public :: check_block_two_stage_step_driver
  public :: check_block_two_stage_step_completion
  public :: check_parallel_block_lifecycle
  public :: check_parallel_block_scaling
  public :: check_block_hydrostatic_reconstruction

contains


  subroutine clear_parallel_block_state
    ! Release persistent and staging block data before a checkpoint
    ! restart invalidates the geometric-domain state from which those
    ! blocks were constructed. Safe to call when nothing is allocated.

    implicit none

    call clear_block_ghost_exchange_plan
    call clear_block_writeback_plan
    call clear_local_blocks
    call clear_block_staging

    if (allocated(block_catalog)) deallocate(block_catalog)

    if (local_block_store_ready()) then
       call fail("local block store remained ready after reset")
    end if

    if (allocated(block_catalog)) then
       call fail("block catalogue remained allocated after reset")
    end if
    if (ghost_exchange_plan%ready) then
       call fail("ghost exchange plan remained ready after reset")
    end if
    if (block_writeback_plan%ready) then
       call fail("writeback plan remained ready after reset")
    end if

  end subroutine clear_parallel_block_state


  logical function parallel_block_state_is_ready () result(ready)
    ! Report whether the installed store and both persistent communication
    ! plans form a complete state for stepping or lifecycle synchronization.

    implicit none

    logical :: local_store_ready
    logical :: writeback_ready

    local_store_ready = local_block_store_ready()
    writeback_ready = block_writeback_plan_is_ready()
    ready = allocated(block_catalog) .and. local_store_ready .and. &
         ghost_exchange_plan%ready .and. writeback_ready

  end function parallel_block_state_is_ready


  subroutine invalidate_parallel_block_domain_shadow
    ! The current production time integrator updates authoritative Domain
    ! fields directly. Discard any installed block copy before that update;
    ! otherwise a later block writeback could restore stale field values.
    ! This transition is intentionally idempotent.

    implicit none

    logical :: block_state_present
    logical :: local_store_ready

    local_store_ready = local_block_store_ready()
    block_state_present = allocated(block_catalog)
    if (local_store_ready) block_state_present = .true.
    if (ghost_exchange_plan%ready) block_state_present = .true.
    if (block_writeback_plan%ready) block_state_present = .true.
    if (.not. block_state_present) return

    call clear_parallel_block_state

  end subroutine invalidate_parallel_block_domain_shadow


  subroutine synchronize_parallel_block_checkpoint
    ! In a block-authoritative integration, materialize both prognostic field
    ! families in Domain storage while retaining every persistent object.

    implicit none

    logical :: checkpoint_ready
    logical :: state_ready
    logical :: trial_active

    trial_active = local_block_tendency_trial_is_active()
    checkpoint_ready = &
         local_block_tendency_commit_checkpoint_is_ready()
    if (trial_active) then
       call fail("checkpoint synchronization during active block trial")
    end if
    if (checkpoint_ready) then
       call fail("checkpoint synchronization before step completion")
    end if

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("checkpoint synchronization before block state is ready")
    end if

    call write_block_field_family_to_domains(BLOCK_PAYLOAD_SOL)
    call write_block_field_family_to_domains(BLOCK_PAYLOAD_WAV_COEFF)

  end subroutine synchronize_parallel_block_checkpoint


  subroutine prepare_parallel_block_grid_change
    ! In a block-authoritative integration, synchronize the current fields,
    ! then invalidate objects whose topology changes during adaptation.

    implicit none

    logical :: state_ready

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("grid-change preparation before block state is ready")
    end if

    call synchronize_parallel_block_checkpoint
    call clear_parallel_block_state

    state_ready = parallel_block_state_is_ready()
    if (state_ready) then
       call fail("grid-change preparation retained stale block state")
    end if

  end subroutine prepare_parallel_block_grid_change


  subroutine build_parallel_block_catalog
  ! Construct the parallel-block decomposition of all locally owned
  ! geometric root domains, build the replicated global block catalogue,
  ! and assign the target block owners using load balancing.
  !
  ! Rule:
  !   - if root patch 1 has all four children, use those four child
  !     subtrees as candidate blocks;
  !   - otherwise keep the whole root domain as one candidate block.
  !
  ! This routine does not modify the active geometric-domain storage.

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

  if (local_block_store_ready()) then
     call fail("existing local block store was not reset")
  end if

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
          "build_parallel_block_catalog: local block count mismatch"
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
     error stop "build_parallel_block_catalog: MPI_Exscan failed"
  end if

  if (rank == 0) n_block_before = 0

  !
  ! Determine the total number of candidate blocks.
  !
  call MPI_Allreduce( &
       n_block_local, n_block_total, &
       1, MPI_INTEGER, MPI_SUM, comm, ierr)

  if (ierr /= MPI_SUCCESS) then
     error stop "build_parallel_block_catalog: MPI_Allreduce failed"
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
     error stop "build_parallel_block_catalog: MPI_Allgather failed"
  end if

  if (sum(block_count) /= n_block_total) then
     error stop &
          "build_parallel_block_catalog: inconsistent global block count"
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
          "build_parallel_block_catalog: local packing count mismatch"
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
     error stop "build_parallel_block_catalog: MPI_Allgatherv failed"
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
             "build_parallel_block_catalog: invalid global block ID ordering"
     end if

     if (block_catalog(i)%owner < 0 .or. &
          block_catalog(i)%owner >= n_process) then
        error stop &
             "build_parallel_block_catalog: invalid block owner"
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
             "build_parallel_block_catalog: invalid root domain"
     end if

     domain_count(block_catalog(i)%root_domain+1) = &
          domain_count(block_catalog(i)%root_domain+1) + 1

  end do

  if (any(domain_count == 0)) then
     error stop &
          "build_parallel_block_catalog: one or more root domains missing"
  end if

  if (any(domain_count /= 1 .and. domain_count /= N_CHDRN)) then
     error stop &
          "build_parallel_block_catalog: invalid number of blocks per root domain"
  end if

  !
  ! ---------------------------------------------------------------
  ! Compute the target load balancing of candidate parallel blocks.
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
             "build_parallel_block_catalog: invalid current block owner"
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
          "build_parallel_block_catalog: one or more candidate blocks unassigned"
  end if

  !
  ! Independently reconstruct proposed rank loads from the owner map.
  !
  load_proposed = 0

  do b = 1, n_block_total

     r = block_owner_new(b) + 1

     if (r < 1 .or. r > n_process) then
        error stop &
             "build_parallel_block_catalog: invalid proposed block owner"
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
          "build_parallel_block_catalog: invalid committed block owner"

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

end subroutine build_parallel_block_catalog


  subroutine migrate_blocks (verbose)
    ! Migrate the source-local blocks to their final catalogue owners,
    ! install a self-contained local block store, and release all
    ! source, receive and MPI migration staging storage.

    implicit none

    logical, optional, intent(in) :: verbose

    type(Block_Migration_Manifest) :: manifest

    integer :: b
    integer :: expected_local
    integer :: i
    integer :: ib
    integer :: n_recv_byte
    integer :: n_received
    integer :: n_send_byte
    integer :: n_sent
    integer :: nbyte
    integer :: pos

    integer, allocatable :: local_seen(:)
    integer, allocatable :: seen_catalog(:)
    integer, allocatable :: seen_source(:)
    integer, allocatable :: send_nbyte(:)

    integer(int8), allocatable :: block_buffer(:)
    integer(int8), allocatable :: send_payload(:)

    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. allocated(block_source) .or. &
         .not. allocated(block_source_catalog_index) .or. &
         .not. allocated(block_migrating_source_index)) then
       call fail("source block store is incomplete")
    end if

    if (size(block_source_catalog_index) /= size(block_source)) then
       call fail("source block catalogue map has the wrong extent")
    end if

    call build_block_migration_manifest(manifest)
    call check_block_migration_manifest(manifest,print_local)

    if (manifest%n_send /= size(block_migrating_source_index)) then
       call fail("manifest and source migration counts differ")
    end if

    allocate(send_nbyte(manifest%n_send))
    allocate(seen_source(size(block_source)))

    send_nbyte = 0
    seen_source = 0

    do i = 1, manifest%n_send

       b = manifest%send_block(i)
       ib = findloc(block_source_catalog_index,b,dim=1)

       if (ib < 1 .or. ib > size(block_source)) then
          call fail("manifest source block was not constructed locally")
       end if

       if (findloc(block_migrating_source_index,ib,dim=1) < 1) then
          call fail("manifest source block is not marked for migration")
       end if

       if (seen_source(ib) /= 0) then
          call fail("source block occurs more than once in manifest")
       end if

       send_nbyte(i) = packed_block_nbyte(block_source(ib))

       if (send_nbyte(i) <= 0) then
          call fail("source block has an invalid packed byte count")
       end if

       seen_source(ib) = 1

    end do

    do i = 1, size(block_migrating_source_index)

       ib = block_migrating_source_index(i)

       if (ib < 1 .or. ib > size(block_source)) then
          call fail("migrating source index is invalid")
       end if

       if (seen_source(ib) /= 1) then
          call fail("migrating source block is absent from manifest")
       end if

    end do

    call exchange_block_migration_sizes( &
         manifest,send_nbyte,print_local)

    if (manifest%total_send_nbyte > int(huge(0),int64) .or. &
         manifest%total_recv_nbyte > int(huge(0),int64)) then
       call fail("packed migration size exceeds default integer range")
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
          call fail("source block disappeared before packing")
       end if

       call pack_block(block_source(ib),block_buffer)

       if (size(block_buffer) /= manifest%send_nbyte(i)) then
          call fail("source block packed size changed")
       end if

       nbyte = size(block_buffer)

       if (pos+nbyte > n_send_byte) then
          call fail("outgoing packed payload exceeds its buffer")
       end if

       send_payload(pos+1:pos+nbyte) = block_buffer
       pos = pos + nbyte

    end do

    if (pos /= n_send_byte) then
       call fail("outgoing packed payload has the wrong extent")
    end if

    call exchange_block_migration_payloads( &
         manifest,send_payload,print_local)

    if (allocated(block_received)) deallocate(block_received)
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

       b = manifest%recv_block(i)
       nbyte = manifest%recv_nbyte(i)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("received block catalogue index is invalid")
       end if

       if (seen_catalog(b) /= 0) then
          call fail("received block occurs more than once")
       end if

       if (block_catalog(b)%owner /= rank) then
          call fail("received block has the wrong final owner")
       end if

       if (nbyte <= 0 .or. pos+nbyte > n_recv_byte) then
          call fail("received packed block extent is invalid")
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
          call fail("received block identity does not match catalogue")
       end if

       call check_block_storage(block_received(i),.true.)

       block_received_catalog_index(i) = b
       seen_catalog(b) = 1
       pos = pos + nbyte

    end do

    if (pos /= n_recv_byte) then
       call fail("received packed payload has the wrong extent")
    end if

    if (count(seen_catalog /= 0) /= manifest%n_recv) then
       call fail("received block inventory is incomplete")
    end if

    call install_local_blocks(size(block_catalog),local_seen)

    expected_local = 0

    do b = 1, size(block_catalog)
       if (block_catalog(b)%owner == rank) then
          expected_local = expected_local + 1
          if (local_seen(b) /= 1) then
             call fail("owned block is absent from installed store")
          end if
       else
          if (local_seen(b) /= 0) then
             call fail("nonlocal block appears in installed store")
          end if
       end if
    end do

    if (n_local_blocks() /= expected_local) then
       call fail("installed block count does not match final ownership")
    end if

    call check_block_scalar_stencil_consumer(print_local)
    call check_block_vector_stencil_consumer(print_local)
    call check_block_boundary_routes(print_local)
    call check_block_ghost_source_addresses(print_local)
    call check_block_ghost_request_manifest(print_local)
    call build_block_ghost_exchange_plan
    call check_block_ghost_exchange_plan(print_local)
    call build_block_writeback_plan
    call check_block_writeback_plan(print_local)
    call check_block_writeback_payload_exchange(print_local)
    call check_block_writeback_domain_reconstruction(print_local)
    call check_block_scalar_ghost_payload_exchange(print_local)
    call check_block_vector_ghost_payload_exchange(print_local)
    call check_production_block_ghost_refresh(print_local)
    call check_block_field_family_accessors(print_local)
    call check_block_patch_writable_storage(print_local)
    call check_block_boundary_family_mutators(print_local)
    call check_block_boundary_family_bulk_fill(print_local)

    n_sent     = manifest%n_send
    n_received = manifest%n_recv

    call clear_block_staging
    call clear_block_migration_manifest(manifest)

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Completed block migration for rank ", rank, ":"
       write(6,'(a,i0)') "  blocks sent      = ", n_sent
       write(6,'(a,i0)') "  blocks received  = ", n_received
       write(6,'(a,i0)') &
            "  blocks installed = ", n_local_blocks()
       write(6,'(a,/)') &
            "  migration staging storage released"
    end if

    deallocate(local_seen)
    deallocate(seen_catalog)
    deallocate(seen_source)
    deallocate(send_nbyte)
    deallocate(send_payload)

    call check_local_blocks(print_local)
    call check_block_field_inventory(print_local)
    if (compressible) call ensure_local_block_hydrostatic_state
    call check_block_hydrostatic_state_accessors(print_local)
    call check_block_hydrostatic_consumer(print_local)
    call check_block_field_consumer(print_local)
    call check_block_stencil_kernel(print_local)
    call check_block_tendency_kernel(print_local)
    call check_block_tendency_trial_update(print_local)
    call check_block_tendency_commit(print_local)
    call check_block_tendency_step_driver(print_local)
    call check_block_tendency_accepted_step(print_local)
    call check_block_multistage_tendency_accumulator(print_local)
    call check_block_multistage_tendency_commit(print_local)
    call check_block_two_stage_step_driver(print_local)
    call check_block_two_stage_step_completion(print_local)
    call check_parallel_block_lifecycle(print_local)
    call check_parallel_block_scaling(print_local)
    call check_block_hydrostatic_reconstruction(print_local)

  end subroutine migrate_blocks


  subroutine check_block_scalar_stencil_consumer (verbose)
    ! Compare a read-only traversal of every scalar compact-stencil window
    ! before migration with the installed final-owner block store. This is
    ! the first persistent-store consumer of patch, boundary and ghost
    ! addressing across every scalar variable and represented level.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: field_kind
    integer :: ierr

    integer(int64) :: source_address_local(3)
    integer(int64) :: source_address_global(3)
    integer(int64) :: local_address_local(3)
    integer(int64) :: local_address_global(3)
    integer(int64) :: source_value_local
    integer(int64) :: source_value_global
    integer(int64) :: local_value_local
    integer(int64) :: local_value_global

    real(dp) :: source_moment_local(3,3)
    real(dp) :: source_moment_global(3,3)
    real(dp) :: local_moment_local(3,3)
    real(dp) :: local_moment_global(3,3)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call source_block_scalar_stencil_statistics( &
         source_address_local,source_value_local,source_moment_local)

    call local_block_scalar_stencil_statistics( &
         local_address_local,local_value_local,local_moment_local)

    call MPI_Allreduce( &
         source_address_local,source_address_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source stencil addresses")

    call MPI_Allreduce( &
         local_address_local,local_address_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final stencil addresses")

    call MPI_Allreduce( &
         source_value_local,source_value_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source stencil values")

    call MPI_Allreduce( &
         local_value_local,local_value_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final stencil values")

    call MPI_Allreduce( &
         source_moment_local,source_moment_global,9, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source stencil moments")

    call MPI_Allreduce( &
         local_moment_local,local_moment_global,9, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final stencil moments")

    if (any(source_address_global /= local_address_global)) then
       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Source/final scalar-stencil-address mismatch:"
          write(error_unit,'(a,3(i0,1x))') &
               "  source patch/boundary/ghost = ", &
               source_address_global
          write(error_unit,'(a,3(i0,1x))') &
               "  final  patch/boundary/ghost = ", &
               local_address_global
       end if
       call fail("source and final scalar stencil addresses differ")
    end if

    if (source_value_global /= local_value_global) then
       call fail("source and final scalar stencil value counts differ")
    end if

    do field_kind = 1, 3
       if (.not. field_moments_match( &
            source_moment_global(:,field_kind), &
            local_moment_global(:,field_kind), &
            source_value_global)) then

          if (rank == 0) then
             write(error_unit,'(/,a,i0,a)') &
                  "Source/final scalar-stencil moments for field kind ", &
                  field_kind,":"
             write(error_unit,'(a,3(es24.16,1x))') &
                  "  source moments = ", &
                  source_moment_global(:,field_kind)
             write(error_unit,'(a,3(es24.16,1x))') &
                  "  final moments  = ", &
                  local_moment_global(:,field_kind)
          end if

          call fail("source and final scalar stencil moments differ")
       end if
    end do

    if (any(local_address_global <= 0_int64) .or. &
         local_value_global <= 0_int64) then
       call fail("incomplete final-owner scalar stencil inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Final-owner scalar stencil consumer for rank ",rank,":"
       write(6,'(a,3(i0,1x))') &
            "  local patch/boundary/ghost addresses = ", &
            local_address_local
       write(6,'(a,i0)') &
            "  local scalar field samples = ",local_value_local
       write(6,'(a,/)') &
            "  source/final scalar stencil checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,3(i0,1x))') &
            "Global patch/boundary/ghost stencil addresses = ", &
            local_address_global
       write(6,'(a,i0)') &
            "Global scalar stencil field samples verified = ", &
            local_value_global
       write(6,'(a,3(es24.16,1x))') &
            "Global sol stencil moments verified = ", &
            local_moment_global(:,1)
       write(6,'(a,3(es24.16,1x))') &
            "Global sol_mean stencil moments verified = ", &
            local_moment_global(:,2)
       write(6,'(a,3(es24.16,1x))') &
            "Global wav_coeff stencil moments verified = ", &
            local_moment_global(:,3)
       write(6,'(a,/)') &
            "Final-owner scalar stencil consumer matches source blocks"
    end if

  end subroutine check_block_scalar_stencil_consumer


  subroutine check_block_vector_stencil_consumer (verbose)
    ! Compare a read-only traversal of every vector compact-stencil window
    ! before migration with the installed final-owner block store. Every
    ! stored edge component and represented level is included.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: field_kind
    integer :: ierr

    integer(int64) :: source_address_local(3)
    integer(int64) :: source_address_global(3)
    integer(int64) :: local_address_local(3)
    integer(int64) :: local_address_global(3)
    integer(int64) :: source_value_local
    integer(int64) :: source_value_global
    integer(int64) :: local_value_local
    integer(int64) :: local_value_global

    real(dp) :: source_moment_local(3,3)
    real(dp) :: source_moment_global(3,3)
    real(dp) :: local_moment_local(3,3)
    real(dp) :: local_moment_global(3,3)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call source_block_vector_stencil_statistics( &
         source_address_local,source_value_local,source_moment_local)

    call local_block_vector_stencil_statistics( &
         local_address_local,local_value_local,local_moment_local)

    call MPI_Allreduce( &
         source_address_local,source_address_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source vector addresses")

    call MPI_Allreduce( &
         local_address_local,local_address_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final vector addresses")

    call MPI_Allreduce( &
         source_value_local,source_value_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source vector values")

    call MPI_Allreduce( &
         local_value_local,local_value_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final vector values")

    call MPI_Allreduce( &
         source_moment_local,source_moment_global,9, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source vector moments")

    call MPI_Allreduce( &
         local_moment_local,local_moment_global,9, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final vector moments")

    if (any(source_address_global /= local_address_global)) then
       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Source/final vector-stencil-address mismatch:"
          write(error_unit,'(a,3(i0,1x))') &
               "  source patch/boundary/ghost = ", &
               source_address_global
          write(error_unit,'(a,3(i0,1x))') &
               "  final  patch/boundary/ghost = ", &
               local_address_global
       end if
       call fail("source and final vector stencil addresses differ")
    end if

    if (source_value_global /= local_value_global) then
       call fail("source and final vector stencil value counts differ")
    end if

    do field_kind = 1, 3
       if (.not. field_moments_match( &
            source_moment_global(:,field_kind), &
            local_moment_global(:,field_kind), &
            source_value_global)) then

          if (rank == 0) then
             write(error_unit,'(/,a,i0,a)') &
                  "Source/final vector-stencil moments for field kind ", &
                  field_kind,":"
             write(error_unit,'(a,3(es24.16,1x))') &
                  "  source moments = ", &
                  source_moment_global(:,field_kind)
             write(error_unit,'(a,3(es24.16,1x))') &
                  "  final moments  = ", &
                  local_moment_global(:,field_kind)
          end if

          call fail("source and final vector stencil moments differ")
       end if
    end do

    if (any(local_address_global <= 0_int64) .or. &
         local_value_global <= 0_int64) then
       call fail("incomplete final-owner vector stencil inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Final-owner vector stencil consumer for rank ",rank,":"
       write(6,'(a,3(i0,1x))') &
            "  local patch/boundary/ghost addresses = ", &
            local_address_local
       write(6,'(a,i0)') &
            "  local vector field samples = ",local_value_local
       write(6,'(a,/)') &
            "  source/final vector stencil checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,3(i0,1x))') &
            "Global vector patch/boundary/ghost addresses = ", &
            local_address_global
       write(6,'(a,i0)') &
            "Global vector stencil field samples verified = ", &
            local_value_global
       write(6,'(a,3(es24.16,1x))') &
            "Global velocity stencil moments verified = ", &
            local_moment_global(:,1)
       write(6,'(a,3(es24.16,1x))') &
            "Global mean-velocity stencil moments verified = ", &
            local_moment_global(:,2)
       write(6,'(a,3(es24.16,1x))') &
            "Global wavelet-velocity stencil moments verified = ", &
            local_moment_global(:,3)
       write(6,'(a,/)') &
            "Final-owner vector stencil consumer matches source blocks"
    end if

  end subroutine check_block_vector_stencil_consumer


  subroutine check_refreshed_block_stencil_consumers (verbose)
    ! Both ghost-payload exchanges deliberately poison their destination
    ! storage before installing the requested sol and wav_coeff values.
    ! Re-run the scalar and vector compact-stencil consumers only after
    ! both installations have completed. This verifies that the refreshed
    ! ghost payloads are immediately usable through the same explicit
    ! addressing paths used by block-local numerical kernels.

    implicit none

    logical, optional, intent(in) :: verbose

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call check_block_scalar_stencil_consumer(.false.)
    call check_block_vector_stencil_consumer(.false.)

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Refreshed ghost stencil consumers for rank ",rank,":"
       write(6,'(a)') &
            "  scalar compact-stencil reads passed"
       write(6,'(a,/)') &
            "  vector compact-stencil reads passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Installed sol/wav_coeff ghost stencil-read checks passed"
    end if

  end subroutine check_refreshed_block_stencil_consumers


  subroutine check_block_boundary_routes (verbose)
    ! Validate the common route catalogue used by update_bdry for scalar,
    ! rank-one and rank-two Float_Field arguments. AT_NODE uses one value
    ! per stored node; AT_EDGE uses EDGE values and the signed receive rule.
    ! This stage validates topology and payload extents but performs no
    ! production boundary exchange and changes no field values.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: source_link_local(3)
    integer(int64) :: source_link_global(3)
    integer(int64) :: local_link_local(3)
    integer(int64) :: local_link_global(3)
    integer(int64) :: source_storage_local(2)
    integer(int64) :: source_storage_global(2)
    integer(int64) :: local_storage_local(2)
    integer(int64) :: local_storage_global(2)
    integer(int64) :: source_node_local(2)
    integer(int64) :: source_node_global(2)
    integer(int64) :: local_node_local(2)
    integer(int64) :: local_node_global(2)
    integer(int64) :: node_payload_local
    integer(int64) :: node_payload_global
    integer(int64) :: edge_payload_local
    integer(int64) :: edge_payload_global

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call source_block_boundary_route_statistics( &
         source_link_local,source_storage_local,source_node_local)
    call local_block_boundary_route_statistics( &
         local_link_local,local_storage_local,local_node_local)

    call MPI_Allreduce( &
         source_link_local,source_link_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source boundary links")

    call MPI_Allreduce( &
         local_link_local,local_link_global,3, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final boundary links")

    call MPI_Allreduce( &
         source_storage_local,source_storage_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source boundary storage")

    call MPI_Allreduce( &
         local_storage_local,local_storage_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final boundary storage")

    call MPI_Allreduce( &
         source_node_local,source_node_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source boundary nodes")

    call MPI_Allreduce( &
         local_node_local,local_node_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final boundary nodes")

    if (any(source_link_global /= local_link_global)) then
       call fail("source and final boundary link counts differ")
    end if
    if (any(source_storage_global /= local_storage_global)) then
       call fail("source and final boundary storage counts differ")
    end if
    if (any(source_node_global /= local_node_global)) then
       call fail("source and final boundary node counts differ")
    end if

    if (sum(local_link_global) <= 0_int64 .or. &
         local_link_global(1) <= 0_int64 .or. &
         sum(local_storage_global) <= 0_int64 .or. &
         sum(local_node_global) <= 0_int64) then
       call fail("incomplete final-owner boundary route inventory")
    end if

    node_payload_local  = sum(local_node_local)
    node_payload_global = sum(local_node_global)
    edge_payload_local  = int(EDGE,int64)*node_payload_local
    edge_payload_global = int(EDGE,int64)*node_payload_global

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Generic Float_Field boundary routes for rank ",rank,":"
       write(6,'(a,3(i0,1x))') &
            "  local block/domain/adaptive links = ",local_link_local
       write(6,'(a,2(i0,1x))') &
            "  local ghost/boundary records      = ",local_storage_local
       write(6,'(a,i0)') &
            "  local AT_NODE values per field    = ",node_payload_local
       write(6,'(a,i0)') &
            "  local AT_EDGE values per field    = ",edge_payload_local
       write(6,'(a,/)') &
            "  field-independent boundary route checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,3(i0,1x))') &
            "Global block/domain/adaptive boundary links = ", &
            local_link_global
       write(6,'(a,2(i0,1x))') &
            "Global ghost/boundary storage records       = ", &
            local_storage_global
       write(6,'(a,i0)') &
            "Global AT_NODE values per Float_Field       = ", &
            node_payload_global
       write(6,'(a,i0)') &
            "Global AT_EDGE values per Float_Field       = ", &
            edge_payload_global
       write(6,'(a,/)') &
            "Generic Float_Field boundary catalogue matches source blocks"
    end if

  end subroutine check_block_boundary_routes


  subroutine check_block_ghost_source_addresses (verbose)
    ! Validate the compact source patch used to pack every NGB_BLOCK
    ! ghost dynamically. The same address applies to any Float_Field;
    ! field rank controls repetition and field position controls the
    ! one-versus-EDGE value multiplier.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: i
    integer :: ierr

    integer(int64) :: source_ghost_local
    integer(int64) :: source_ghost_global
    integer(int64) :: local_ghost_local
    integer(int64) :: local_ghost_global
    integer(int64) :: source_value_local(2)
    integer(int64) :: source_value_global(2)
    integer(int64) :: local_value_local(2)
    integer(int64) :: local_value_global(2)
    integer(int64) :: source_sum_local(5)
    integer(int64) :: source_sum_global(5)
    integer(int64) :: local_sum_local(5)
    integer(int64) :: local_sum_global(5)

    integer, allocatable :: catalog_domain(:)
    integer, allocatable :: catalog_id(:)
    integer, allocatable :: catalog_owner(:)
    integer, allocatable :: patch_count_global(:)
    integer, allocatable :: patch_count_local(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    allocate(catalog_domain(size(block_catalog)))
    allocate(catalog_id(size(block_catalog)))
    allocate(catalog_owner(size(block_catalog)))
    allocate(patch_count_global(size(block_catalog)))
    allocate(patch_count_local(size(block_catalog)))

    patch_count_local = 0

    do b = 1, size(block_catalog)
       catalog_domain(b) = block_catalog(b)%root_domain
       catalog_id(b)     = block_catalog(b)%id
       catalog_owner(b)  = block_catalog(b)%owner
    end do

    do i = 1, size(block_source)
       b = block_source_catalog_index(i)
       if (b < 1 .or. b > size(block_catalog)) then
          call fail("invalid source block catalogue index")
       end if
       if (patch_count_local(b) /= 0) then
          call fail("duplicate source block patch count")
       end if
       patch_count_local(b) = size(block_source(i)%patch)
    end do

    call MPI_Allreduce( &
         patch_count_local,patch_count_global,size(block_catalog), &
         MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source block patch counts")

    if (any(patch_count_global <= 0)) then
       call fail("incomplete global source block patch counts")
    end if

    call validate_local_block_ghost_sources( &
         patch_count_global,catalog_owner,catalog_id,catalog_domain)

    call source_block_ghost_source_statistics( &
         source_ghost_local,source_value_local,source_sum_local)
    call local_block_ghost_source_statistics( &
         local_ghost_local,local_value_local,local_sum_local)

    call MPI_Allreduce( &
         source_ghost_local,source_ghost_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source ghost count")

    call MPI_Allreduce( &
         local_ghost_local,local_ghost_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final ghost count")

    call MPI_Allreduce( &
         source_value_local,source_value_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source ghost values")

    call MPI_Allreduce( &
         local_value_local,local_value_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final ghost values")

    call MPI_Allreduce( &
         source_sum_local,source_sum_global,5, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce source ghost addresses")

    call MPI_Allreduce( &
         local_sum_local,local_sum_global,5, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce final ghost addresses")

    if (source_ghost_global /= local_ghost_global .or. &
         any(source_value_global /= local_value_global) .or. &
         any(source_sum_global /= local_sum_global)) then
       call fail("source and final ghost source mappings differ")
    end if

    if (local_ghost_global <= 0_int64 .or. &
         local_value_global(1) /= &
         int(PATCH_SIZE**2,int64)*local_ghost_global .or. &
         local_value_global(2) /= &
         int(EDGE,int64)*local_value_global(1)) then
       call fail("invalid dynamic ghost source payload extent")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Dynamic Float_Field ghost sources for rank ",rank,":"
       write(6,'(a,i0)') &
            "  local compact source patches = ",local_ghost_local
       write(6,'(a,i0)') &
            "  local AT_NODE values per field = ",local_value_local(1)
       write(6,'(a,i0)') &
            "  local AT_EDGE values per field = ",local_value_local(2)
       write(6,'(a,/)') &
            "  source catalogue and compact-patch checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global dynamically addressable ghost patches = ", &
            local_ghost_global
       write(6,'(a,i0)') &
            "Global AT_NODE inter-block values per field   = ", &
            local_value_global(1)
       write(6,'(a,i0)') &
            "Global AT_EDGE inter-block values per field   = ", &
            local_value_global(2)
       write(6,'(a,/)') &
            "Dynamic Float_Field ghost source mappings match source blocks"
    end if

    deallocate(catalog_domain)
    deallocate(catalog_id)
    deallocate(catalog_owner)
    deallocate(patch_count_global)
    deallocate(patch_count_local)

  end subroutine check_block_ghost_source_addresses


  subroutine clear_block_ghost_exchange_plan
    ! Release field-independent block-ghost routing metadata. The plan is
    ! invalid whenever a checkpoint restart replaces the installed blocks.

    implicit none

    if (allocated(ghost_exchange_plan%source_block)) then
       deallocate(ghost_exchange_plan%source_block)
    end if
    if (allocated(ghost_exchange_plan%source_local_patch)) then
       deallocate(ghost_exchange_plan%source_local_patch)
    end if
    if (allocated(ghost_exchange_plan%source_owner)) then
       deallocate(ghost_exchange_plan%source_owner)
    end if
    if (allocated(ghost_exchange_plan%destination_block)) then
       deallocate(ghost_exchange_plan%destination_block)
    end if
    if (allocated(ghost_exchange_plan%destination_ghost)) then
       deallocate(ghost_exchange_plan%destination_ghost)
    end if
    if (allocated(ghost_exchange_plan%send_record_count)) then
       deallocate(ghost_exchange_plan%send_record_count)
    end if
    if (allocated(ghost_exchange_plan%recv_record_count)) then
       deallocate(ghost_exchange_plan%recv_record_count)
    end if
    if (allocated(ghost_exchange_plan%send_record_displ)) then
       deallocate(ghost_exchange_plan%send_record_displ)
    end if
    if (allocated(ghost_exchange_plan%recv_record_displ)) then
       deallocate(ghost_exchange_plan%recv_record_displ)
    end if
    if (allocated(ghost_exchange_plan%send_data)) then
       deallocate(ghost_exchange_plan%send_data)
    end if
    if (allocated(ghost_exchange_plan%recv_data)) then
       deallocate(ghost_exchange_plan%recv_data)
    end if
    if (allocated(ghost_exchange_plan%request_index)) then
       deallocate(ghost_exchange_plan%request_index)
    end if
    if (allocated(ghost_exchange_plan%scalar_send_count)) then
       deallocate(ghost_exchange_plan%scalar_send_count)
    end if
    if (allocated(ghost_exchange_plan%scalar_recv_count)) then
       deallocate(ghost_exchange_plan%scalar_recv_count)
    end if
    if (allocated(ghost_exchange_plan%scalar_send_displ)) then
       deallocate(ghost_exchange_plan%scalar_send_displ)
    end if
    if (allocated(ghost_exchange_plan%scalar_recv_displ)) then
       deallocate(ghost_exchange_plan%scalar_recv_displ)
    end if
    if (allocated(ghost_exchange_plan%vector_send_count)) then
       deallocate(ghost_exchange_plan%vector_send_count)
    end if
    if (allocated(ghost_exchange_plan%vector_recv_count)) then
       deallocate(ghost_exchange_plan%vector_recv_count)
    end if
    if (allocated(ghost_exchange_plan%vector_send_displ)) then
       deallocate(ghost_exchange_plan%vector_send_displ)
    end if
    if (allocated(ghost_exchange_plan%vector_recv_displ)) then
       deallocate(ghost_exchange_plan%vector_recv_displ)
    end if
    if (allocated(ghost_exchange_plan%scalar_send_buffer)) then
       deallocate(ghost_exchange_plan%scalar_send_buffer)
    end if
    if (allocated(ghost_exchange_plan%scalar_recv_buffer)) then
       deallocate(ghost_exchange_plan%scalar_recv_buffer)
    end if
    if (allocated(ghost_exchange_plan%scalar_patch_buffer)) then
       deallocate(ghost_exchange_plan%scalar_patch_buffer)
    end if
    if (allocated(ghost_exchange_plan%vector_send_buffer)) then
       deallocate(ghost_exchange_plan%vector_send_buffer)
    end if
    if (allocated(ghost_exchange_plan%vector_recv_buffer)) then
       deallocate(ghost_exchange_plan%vector_recv_buffer)
    end if
    if (allocated(ghost_exchange_plan%vector_patch_buffer)) then
       deallocate(ghost_exchange_plan%vector_patch_buffer)
    end if

    ghost_exchange_plan%n_request = 0
    ghost_exchange_plan%n_local_request = 0
    ghost_exchange_plan%n_remote_send = 0
    ghost_exchange_plan%n_remote_recv = 0
    ghost_exchange_plan%scalar_n_value = 0
    ghost_exchange_plan%vector_n_value = 0
    ghost_exchange_plan%ready = .false.

  end subroutine clear_block_ghost_exchange_plan


  subroutine build_block_ghost_exchange_plan
    ! Build the request routing once for the installed final-owner block
    ! store. Subsequent scalar and vector field-family refreshes reuse this
    ! plan, its MPI payload buffers and one-patch packing workspaces.

    implicit none

    integer, parameter :: REQUEST_SIZE = 4

    integer :: destination
    integer :: destination_nghost
    integer :: field_level
    integer :: fill_record
    integer :: i
    integer :: ierr
    integer :: level_count
    integer :: pos
    integer :: r
    integer :: scalar_count
    integer :: scalar_mult
    integer :: scalar_variable
    integer :: source
    integer :: source_npatch
    integer :: vector_mult
    integer :: vector_variable

    integer, allocatable :: fill(:)
    integer, allocatable :: recv_count(:)
    integer, allocatable :: recv_displ(:)
    integer, allocatable :: send_count(:)
    integer, allocatable :: send_displ(:)

    call clear_block_ghost_exchange_plan

    if (.not. local_block_store_ready()) then
       call fail("block ghost exchange plan before block installation")
    end if
    if (.not. allocated(block_catalog)) then
       call fail("block ghost exchange plan without block catalogue")
    end if

    call get_local_block_ghost_requests( &
         ghost_exchange_plan%source_block, &
         ghost_exchange_plan%source_local_patch, &
         ghost_exchange_plan%source_owner, &
         ghost_exchange_plan%destination_block, &
         ghost_exchange_plan%destination_ghost)

    ghost_exchange_plan%n_request = &
         size(ghost_exchange_plan%source_block)

    if (size(ghost_exchange_plan%source_local_patch) /= &
         ghost_exchange_plan%n_request .or. &
         size(ghost_exchange_plan%source_owner) /= &
         ghost_exchange_plan%n_request .or. &
         size(ghost_exchange_plan%destination_block) /= &
         ghost_exchange_plan%n_request .or. &
         size(ghost_exchange_plan%destination_ghost) /= &
         ghost_exchange_plan%n_request) then
       call fail("inconsistent block ghost exchange plan arrays")
    end if

    allocate(ghost_exchange_plan%send_record_count(n_process))
    allocate(ghost_exchange_plan%recv_record_count(n_process))
    allocate(ghost_exchange_plan%send_record_displ(n_process))
    allocate(ghost_exchange_plan%recv_record_displ(n_process))
    allocate(fill(n_process))

    ghost_exchange_plan%send_record_count = 0
    ghost_exchange_plan%n_local_request = 0

    do i = 1, ghost_exchange_plan%n_request
       source = ghost_exchange_plan%source_block(i)
       destination = ghost_exchange_plan%destination_block(i)

       if (source < 1 .or. source > size(block_catalog)) then
          call fail("invalid source block in ghost exchange plan")
       end if
       if (destination < 1 .or. &
            destination > size(block_catalog)) then
          call fail("invalid destination block in ghost exchange plan")
       end if
       if (ghost_exchange_plan%source_owner(i) /= &
            block_catalog(source)%owner) then
          call fail("stale source owner in ghost exchange plan")
       end if
       if (block_catalog(destination)%owner /= rank) then
          call fail("nonlocal destination in ghost exchange plan")
       end if

       destination_nghost = local_block_ghost_count(destination)
       if (ghost_exchange_plan%destination_ghost(i) < 1 .or. &
            ghost_exchange_plan%destination_ghost(i) > &
            destination_nghost) then
          call fail("invalid destination ghost in exchange plan")
       end if

       if (ghost_exchange_plan%source_owner(i) == rank) then
          if (catalog_local_block(source) < 1) then
             call fail("local source block absent from ghost plan")
          end if
          source_npatch = local_block_patch_count(source)
          if (ghost_exchange_plan%source_local_patch(i) < 0 .or. &
               ghost_exchange_plan%source_local_patch(i) >= &
               source_npatch) then
             call fail("invalid local source patch in ghost plan")
          end if
          ghost_exchange_plan%n_local_request = &
               ghost_exchange_plan%n_local_request + 1
       else
          r = ghost_exchange_plan%source_owner(i) + 1
          ghost_exchange_plan%send_record_count(r) = &
               ghost_exchange_plan%send_record_count(r) + 1
       end if
    end do

    call MPI_Alltoall( &
         ghost_exchange_plan%send_record_count,1,MPI_INTEGER, &
         ghost_exchange_plan%recv_record_count,1,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoall persistent ghost request counts")

    ghost_exchange_plan%n_remote_send = &
         sum(ghost_exchange_plan%send_record_count)
    ghost_exchange_plan%n_remote_recv = &
         sum(ghost_exchange_plan%recv_record_count)

    ghost_exchange_plan%send_record_displ(1) = 0
    ghost_exchange_plan%recv_record_displ(1) = 0

    do r = 2, n_process
       ghost_exchange_plan%send_record_displ(r) = &
            ghost_exchange_plan%send_record_displ(r-1) + &
            ghost_exchange_plan%send_record_count(r-1)
       ghost_exchange_plan%recv_record_displ(r) = &
            ghost_exchange_plan%recv_record_displ(r-1) + &
            ghost_exchange_plan%recv_record_count(r-1)
    end do

    allocate(ghost_exchange_plan%send_data( &
         REQUEST_SIZE*ghost_exchange_plan%n_remote_send))
    allocate(ghost_exchange_plan%recv_data( &
         REQUEST_SIZE*ghost_exchange_plan%n_remote_recv))
    allocate(ghost_exchange_plan%request_index( &
         ghost_exchange_plan%n_remote_send))

    fill = 0

    do i = 1, ghost_exchange_plan%n_request
       if (ghost_exchange_plan%source_owner(i) == rank) cycle

       r = ghost_exchange_plan%source_owner(i) + 1
       pos = REQUEST_SIZE*( &
            ghost_exchange_plan%send_record_displ(r)+fill(r))
       ghost_exchange_plan%send_data(pos+1:pos+REQUEST_SIZE) = [ &
            ghost_exchange_plan%source_block(i), &
            ghost_exchange_plan%source_local_patch(i), &
            ghost_exchange_plan%destination_block(i), &
            ghost_exchange_plan%destination_ghost(i) ]
       fill_record = ghost_exchange_plan%send_record_displ(r) + &
            fill(r) + 1
       ghost_exchange_plan%request_index(fill_record) = i
       fill(r) = fill(r) + 1
    end do

    if (any(fill /= ghost_exchange_plan%send_record_count)) then
       call fail("persistent ghost request packing count mismatch")
    end if

    allocate(send_count(n_process))
    allocate(recv_count(n_process))
    allocate(send_displ(n_process))
    allocate(recv_displ(n_process))

    send_count = REQUEST_SIZE*ghost_exchange_plan%send_record_count
    recv_count = REQUEST_SIZE*ghost_exchange_plan%recv_record_count
    send_displ = REQUEST_SIZE*ghost_exchange_plan%send_record_displ
    recv_displ = REQUEST_SIZE*ghost_exchange_plan%recv_record_displ

    call MPI_Alltoallv( &
         ghost_exchange_plan%send_data,send_count,send_displ,MPI_INTEGER, &
         ghost_exchange_plan%recv_data,recv_count,recv_displ,MPI_INTEGER, &
         comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv persistent ghost requests")

    do r = 1, n_process
       do i = 0, ghost_exchange_plan%recv_record_count(r)-1
          pos = REQUEST_SIZE*( &
               ghost_exchange_plan%recv_record_displ(r)+i)
          source = ghost_exchange_plan%recv_data(pos+1)

          if (source < 1 .or. source > size(block_catalog)) then
             call fail("received persistent ghost source is invalid")
          end if
          if (block_catalog(source)%owner /= rank) then
             call fail("received persistent ghost source is not local")
          end if
          if (catalog_local_block(source) < 1) then
             call fail("received persistent ghost source is absent")
          end if

          source_npatch = local_block_patch_count(source)
          if (ghost_exchange_plan%recv_data(pos+2) < 0 .or. &
               ghost_exchange_plan%recv_data(pos+2) >= &
               source_npatch) then
             call fail("received persistent source patch is invalid")
          end if

          destination = ghost_exchange_plan%recv_data(pos+3)
          if (destination < 1 .or. &
               destination > size(block_catalog)) then
             call fail("received persistent destination is invalid")
          end if
          if (block_catalog(destination)%owner /= r-1) then
             call fail("received persistent destination owner mismatch")
          end if
          if (ghost_exchange_plan%recv_data(pos+4) < 1) then
             call fail("received persistent ghost index is invalid")
          end if
       end do
    end do

    call get_block_field_layout( &
         scalar_variable,scalar_count,vector_variable,field_level, &
         level_count,scalar_mult,vector_mult)

    if (scalar_count < 1 .or. level_count < 1 .or. &
         scalar_mult < 1 .or. vector_mult < 1) then
       call fail("invalid persistent ghost payload layout")
    end if

    ghost_exchange_plan%scalar_n_value = &
         scalar_count*level_count*scalar_mult*PATCH_SIZE**2
    ghost_exchange_plan%vector_n_value = &
         level_count*vector_mult*PATCH_SIZE**2

    allocate(ghost_exchange_plan%scalar_send_count(n_process))
    allocate(ghost_exchange_plan%scalar_recv_count(n_process))
    allocate(ghost_exchange_plan%scalar_send_displ(n_process))
    allocate(ghost_exchange_plan%scalar_recv_displ(n_process))
    allocate(ghost_exchange_plan%vector_send_count(n_process))
    allocate(ghost_exchange_plan%vector_recv_count(n_process))
    allocate(ghost_exchange_plan%vector_send_displ(n_process))
    allocate(ghost_exchange_plan%vector_recv_displ(n_process))

    ghost_exchange_plan%scalar_send_count = &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%recv_record_count
    ghost_exchange_plan%scalar_recv_count = &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%send_record_count
    ghost_exchange_plan%vector_send_count = &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%recv_record_count
    ghost_exchange_plan%vector_recv_count = &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%send_record_count

    ghost_exchange_plan%scalar_send_displ(1) = 0
    ghost_exchange_plan%scalar_recv_displ(1) = 0
    ghost_exchange_plan%vector_send_displ(1) = 0
    ghost_exchange_plan%vector_recv_displ(1) = 0

    do r = 2, n_process
       ghost_exchange_plan%scalar_send_displ(r) = &
            ghost_exchange_plan%scalar_send_displ(r-1) + &
            ghost_exchange_plan%scalar_send_count(r-1)
       ghost_exchange_plan%scalar_recv_displ(r) = &
            ghost_exchange_plan%scalar_recv_displ(r-1) + &
            ghost_exchange_plan%scalar_recv_count(r-1)
       ghost_exchange_plan%vector_send_displ(r) = &
            ghost_exchange_plan%vector_send_displ(r-1) + &
            ghost_exchange_plan%vector_send_count(r-1)
       ghost_exchange_plan%vector_recv_displ(r) = &
            ghost_exchange_plan%vector_recv_displ(r-1) + &
            ghost_exchange_plan%vector_recv_count(r-1)
    end do

    allocate(ghost_exchange_plan%scalar_send_buffer( &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%n_remote_recv))
    allocate(ghost_exchange_plan%scalar_recv_buffer( &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%n_remote_send))
    allocate(ghost_exchange_plan%scalar_patch_buffer( &
         ghost_exchange_plan%scalar_n_value))
    allocate(ghost_exchange_plan%vector_send_buffer( &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%n_remote_recv))
    allocate(ghost_exchange_plan%vector_recv_buffer( &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%n_remote_send))
    allocate(ghost_exchange_plan%vector_patch_buffer( &
         ghost_exchange_plan%vector_n_value))

    ghost_exchange_plan%scalar_send_buffer = 0.0_dp
    ghost_exchange_plan%scalar_recv_buffer = 0.0_dp
    ghost_exchange_plan%scalar_patch_buffer = 0.0_dp
    ghost_exchange_plan%vector_send_buffer = 0.0_dp
    ghost_exchange_plan%vector_recv_buffer = 0.0_dp
    ghost_exchange_plan%vector_patch_buffer = 0.0_dp

    ghost_exchange_plan%ready = .true.

    deallocate(fill)
    deallocate(recv_count)
    deallocate(recv_displ)
    deallocate(send_count)
    deallocate(send_displ)

  end subroutine build_block_ghost_exchange_plan


  subroutine check_block_ghost_exchange_plan (verbose)
    ! Compare the cached plan with a fresh local inventory and verify the
    ! global request counts and request-record checksums.

    implicit none

    integer, parameter :: REQUEST_SIZE = 4

    logical, optional, intent(in) :: verbose

    integer :: ierr
    integer :: n_request
    integer :: r

    integer(int64) :: count_global(3)
    integer(int64) :: count_local(3)
    integer(int64) :: buffer_count_global(4)
    integer(int64) :: buffer_count_local(4)
    integer(int64) :: recv_sum_global(REQUEST_SIZE)
    integer(int64) :: recv_sum_local(REQUEST_SIZE)
    integer(int64) :: send_sum_global(REQUEST_SIZE)
    integer(int64) :: send_sum_local(REQUEST_SIZE)

    integer, allocatable :: destination_block(:)
    integer, allocatable :: destination_ghost(:)
    integer, allocatable :: source_block(:)
    integer, allocatable :: source_local_patch(:)
    integer, allocatable :: source_owner(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. ghost_exchange_plan%ready) then
       call fail("persistent block ghost exchange plan is not ready")
    end if

    if (.not. allocated(ghost_exchange_plan%scalar_send_count)) then
       call fail("persistent scalar send counts are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_recv_count)) then
       call fail("persistent scalar receive counts are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_send_displ)) then
       call fail("persistent scalar send displacements are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_recv_displ)) then
       call fail("persistent scalar receive displacements are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_send_count)) then
       call fail("persistent vector send counts are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_recv_count)) then
       call fail("persistent vector receive counts are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_send_displ)) then
       call fail("persistent vector send displacements are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_recv_displ)) then
       call fail("persistent vector receive displacements are absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_send_buffer)) then
       call fail("persistent scalar send buffer is absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_recv_buffer)) then
       call fail("persistent scalar receive buffer is absent")
    end if
    if (.not. allocated(ghost_exchange_plan%scalar_patch_buffer)) then
       call fail("persistent scalar patch buffer is absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_send_buffer)) then
       call fail("persistent vector send buffer is absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_recv_buffer)) then
       call fail("persistent vector receive buffer is absent")
    end if
    if (.not. allocated(ghost_exchange_plan%vector_patch_buffer)) then
       call fail("persistent vector patch buffer is absent")
    end if

    if (ghost_exchange_plan%scalar_n_value <= 0 .or. &
         ghost_exchange_plan%vector_n_value <= 0) then
       call fail("persistent ghost payload size is invalid")
    end if

    if (size(ghost_exchange_plan%scalar_send_count) /= n_process .or. &
         size(ghost_exchange_plan%scalar_recv_count) /= n_process .or. &
         size(ghost_exchange_plan%scalar_send_displ) /= n_process .or. &
         size(ghost_exchange_plan%scalar_recv_displ) /= n_process .or. &
         size(ghost_exchange_plan%vector_send_count) /= n_process .or. &
         size(ghost_exchange_plan%vector_recv_count) /= n_process .or. &
         size(ghost_exchange_plan%vector_send_displ) /= n_process .or. &
         size(ghost_exchange_plan%vector_recv_displ) /= n_process) then
       call fail("persistent ghost payload process layout mismatch")
    end if

    if (any(ghost_exchange_plan%scalar_send_count /= &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%recv_record_count) .or. &
         any(ghost_exchange_plan%scalar_recv_count /= &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%send_record_count) .or. &
         any(ghost_exchange_plan%vector_send_count /= &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%recv_record_count) .or. &
         any(ghost_exchange_plan%vector_recv_count /= &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%send_record_count)) then
       call fail("persistent ghost payload counts are invalid")
    end if

    if (ghost_exchange_plan%scalar_send_displ(1) /= 0 .or. &
         ghost_exchange_plan%scalar_recv_displ(1) /= 0 .or. &
         ghost_exchange_plan%vector_send_displ(1) /= 0 .or. &
         ghost_exchange_plan%vector_recv_displ(1) /= 0) then
       call fail("persistent ghost payload displacement origin mismatch")
    end if

    do r = 2, n_process
       if (ghost_exchange_plan%scalar_send_displ(r) /= &
            ghost_exchange_plan%scalar_send_displ(r-1) + &
            ghost_exchange_plan%scalar_send_count(r-1)) then
          call fail("persistent scalar send displacement mismatch")
       end if
       if (ghost_exchange_plan%scalar_recv_displ(r) /= &
            ghost_exchange_plan%scalar_recv_displ(r-1) + &
            ghost_exchange_plan%scalar_recv_count(r-1)) then
          call fail("persistent scalar receive displacement mismatch")
       end if
       if (ghost_exchange_plan%vector_send_displ(r) /= &
            ghost_exchange_plan%vector_send_displ(r-1) + &
            ghost_exchange_plan%vector_send_count(r-1)) then
          call fail("persistent vector send displacement mismatch")
       end if
       if (ghost_exchange_plan%vector_recv_displ(r) /= &
            ghost_exchange_plan%vector_recv_displ(r-1) + &
            ghost_exchange_plan%vector_recv_count(r-1)) then
          call fail("persistent vector receive displacement mismatch")
       end if
    end do

    if (size(ghost_exchange_plan%scalar_send_buffer) /= &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%n_remote_recv .or. &
         size(ghost_exchange_plan%scalar_recv_buffer) /= &
         ghost_exchange_plan%scalar_n_value* &
         ghost_exchange_plan%n_remote_send .or. &
         size(ghost_exchange_plan%vector_send_buffer) /= &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%n_remote_recv .or. &
         size(ghost_exchange_plan%vector_recv_buffer) /= &
         ghost_exchange_plan%vector_n_value* &
         ghost_exchange_plan%n_remote_send .or. &
         size(ghost_exchange_plan%scalar_patch_buffer) /= &
         ghost_exchange_plan%scalar_n_value .or. &
         size(ghost_exchange_plan%vector_patch_buffer) /= &
         ghost_exchange_plan%vector_n_value) then
       call fail("persistent ghost payload buffer extent mismatch")
    end if

    call get_local_block_ghost_requests( &
         source_block,source_local_patch,source_owner, &
         destination_block,destination_ghost)

    n_request = size(source_block)

    if (n_request /= ghost_exchange_plan%n_request) then
       call fail("persistent ghost plan request count changed")
    end if
    if (any(source_block /= ghost_exchange_plan%source_block) .or. &
         any(source_local_patch /= &
         ghost_exchange_plan%source_local_patch) .or. &
         any(source_owner /= ghost_exchange_plan%source_owner) .or. &
         any(destination_block /= &
         ghost_exchange_plan%destination_block) .or. &
         any(destination_ghost /= &
         ghost_exchange_plan%destination_ghost)) then
       call fail("persistent ghost plan differs from block inventory")
    end if

    count_local = int([ &
         ghost_exchange_plan%n_local_request, &
         ghost_exchange_plan%n_remote_send, &
         ghost_exchange_plan%n_remote_recv ],int64)

    send_sum_local = 0_int64
    recv_sum_local = 0_int64

    if (ghost_exchange_plan%n_remote_send > 0) then
       send_sum_local = sum(reshape(int( &
            ghost_exchange_plan%send_data,int64), &
            [REQUEST_SIZE,ghost_exchange_plan%n_remote_send]),dim=2)
    end if
    if (ghost_exchange_plan%n_remote_recv > 0) then
       recv_sum_local = sum(reshape(int( &
            ghost_exchange_plan%recv_data,int64), &
            [REQUEST_SIZE,ghost_exchange_plan%n_remote_recv]),dim=2)
    end if

    call MPI_Allreduce( &
         count_local,count_global,3,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce persistent ghost request totals")

    buffer_count_local = int([ &
         size(ghost_exchange_plan%scalar_send_buffer), &
         size(ghost_exchange_plan%scalar_recv_buffer), &
         size(ghost_exchange_plan%vector_send_buffer), &
         size(ghost_exchange_plan%vector_recv_buffer) ],int64)

    call MPI_Allreduce( &
         buffer_count_local,buffer_count_global,4, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce persistent payload buffer totals")

    call MPI_Allreduce( &
         send_sum_local,send_sum_global,REQUEST_SIZE, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce persistent ghost send checksum")

    call MPI_Allreduce( &
         recv_sum_local,recv_sum_global,REQUEST_SIZE, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce persistent ghost recv checksum")

    if (count_global(2) /= count_global(3) .or. &
         any(send_sum_global /= recv_sum_global)) then
       call fail("persistent ghost plan global routing mismatch")
    end if

    if (buffer_count_global(1) /= buffer_count_global(2) .or. &
         buffer_count_global(3) /= buffer_count_global(4)) then
       call fail("persistent ghost payload global buffer mismatch")
    end if
    if (buffer_count_global(1) /= &
         int(ghost_exchange_plan%scalar_n_value,int64)* &
         count_global(2) .or. &
         buffer_count_global(3) /= &
         int(ghost_exchange_plan%vector_n_value,int64)* &
         count_global(2)) then
       call fail("persistent ghost payload global extent mismatch")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Persistent block ghost exchange plan for rank ",rank,":"
       write(6,'(a,i0)') &
            "  local-source requests    = ", &
            ghost_exchange_plan%n_local_request
       write(6,'(a,i0)') &
            "  remote requests sent     = ", &
            ghost_exchange_plan%n_remote_send
       write(6,'(a,i0)') &
            "  remote requests received = ", &
            ghost_exchange_plan%n_remote_recv
       write(6,'(a,i0)') &
            "  reusable scalar send values = ", &
            size(ghost_exchange_plan%scalar_send_buffer)
       write(6,'(a,i0)') &
            "  reusable scalar receive values = ", &
            size(ghost_exchange_plan%scalar_recv_buffer)
       write(6,'(a,i0)') &
            "  reusable scalar patch values = ", &
            size(ghost_exchange_plan%scalar_patch_buffer)
       write(6,'(a,i0)') &
            "  reusable vector send values = ", &
            size(ghost_exchange_plan%vector_send_buffer)
       write(6,'(a,i0)') &
            "  reusable vector receive values = ", &
            size(ghost_exchange_plan%vector_recv_buffer)
       write(6,'(a,i0)') &
            "  reusable vector patch values = ", &
            size(ghost_exchange_plan%vector_patch_buffer)
       write(6,'(a,/)') &
            "  reusable request and payload buffer checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global persistent block ghost requests = ", &
            count_global(1)+count_global(2)
       write(6,'(a,i0)') &
            "Global reusable scalar send values = ", &
            buffer_count_global(1)
       write(6,'(a,i0)') &
            "Global reusable scalar receive values = ", &
            buffer_count_global(2)
       write(6,'(a,i0)') &
            "Global reusable vector send values = ", &
            buffer_count_global(3)
       write(6,'(a,i0)') &
            "Global reusable vector receive values = ", &
            buffer_count_global(4)
       write(6,'(a,/)') &
            "Persistent block ghost exchange buffers passed"
    end if

    deallocate(destination_block)
    deallocate(destination_ghost)
    deallocate(source_block)
    deallocate(source_local_patch)
    deallocate(source_owner)

  end subroutine check_block_ghost_exchange_plan


  subroutine clear_block_writeback_plan
    ! Release reverse-routing metadata and payload buffers before the block
    ! catalogue or installed store is invalidated.

    implicit none

    if (allocated(block_writeback_plan%send_count)) then
       deallocate(block_writeback_plan%send_count)
    end if
    if (allocated(block_writeback_plan%recv_count)) then
       deallocate(block_writeback_plan%recv_count)
    end if
    if (allocated(block_writeback_plan%send_displ)) then
       deallocate(block_writeback_plan%send_displ)
    end if
    if (allocated(block_writeback_plan%recv_displ)) then
       deallocate(block_writeback_plan%recv_displ)
    end if
    if (allocated(block_writeback_plan%send_block)) then
       deallocate(block_writeback_plan%send_block)
    end if
    if (allocated(block_writeback_plan%recv_block)) then
       deallocate(block_writeback_plan%recv_block)
    end if
    if (allocated(block_writeback_plan%send_patch_count)) then
       deallocate(block_writeback_plan%send_patch_count)
    end if
    if (allocated(block_writeback_plan%recv_patch_count)) then
       deallocate(block_writeback_plan%recv_patch_count)
    end if
    if (allocated(block_writeback_plan%send_scalar_nvalue)) then
       deallocate(block_writeback_plan%send_scalar_nvalue)
    end if
    if (allocated(block_writeback_plan%recv_scalar_nvalue)) then
       deallocate(block_writeback_plan%recv_scalar_nvalue)
    end if
    if (allocated(block_writeback_plan%send_vector_nvalue)) then
       deallocate(block_writeback_plan%send_vector_nvalue)
    end if
    if (allocated(block_writeback_plan%recv_vector_nvalue)) then
       deallocate(block_writeback_plan%recv_vector_nvalue)
    end if
    if (allocated(block_writeback_plan%scalar_send_count)) then
       deallocate(block_writeback_plan%scalar_send_count)
    end if
    if (allocated(block_writeback_plan%scalar_recv_count)) then
       deallocate(block_writeback_plan%scalar_recv_count)
    end if
    if (allocated(block_writeback_plan%scalar_send_displ)) then
       deallocate(block_writeback_plan%scalar_send_displ)
    end if
    if (allocated(block_writeback_plan%scalar_recv_displ)) then
       deallocate(block_writeback_plan%scalar_recv_displ)
    end if
    if (allocated(block_writeback_plan%vector_send_count)) then
       deallocate(block_writeback_plan%vector_send_count)
    end if
    if (allocated(block_writeback_plan%vector_recv_count)) then
       deallocate(block_writeback_plan%vector_recv_count)
    end if
    if (allocated(block_writeback_plan%vector_send_displ)) then
       deallocate(block_writeback_plan%vector_send_displ)
    end if
    if (allocated(block_writeback_plan%vector_recv_displ)) then
       deallocate(block_writeback_plan%vector_recv_displ)
    end if
    if (allocated(block_writeback_plan%scalar_send_buffer)) then
       deallocate(block_writeback_plan%scalar_send_buffer)
    end if
    if (allocated(block_writeback_plan%scalar_recv_buffer)) then
       deallocate(block_writeback_plan%scalar_recv_buffer)
    end if
    if (allocated(block_writeback_plan%vector_send_buffer)) then
       deallocate(block_writeback_plan%vector_send_buffer)
    end if
    if (allocated(block_writeback_plan%vector_recv_buffer)) then
       deallocate(block_writeback_plan%vector_recv_buffer)
    end if
    if (allocated(block_writeback_plan%domain_patch_displ)) then
       deallocate(block_writeback_plan%domain_patch_displ)
    end if
    if (allocated(block_writeback_plan%scalar_domain_stage)) then
       deallocate(block_writeback_plan%scalar_domain_stage)
    end if
    if (allocated(block_writeback_plan%vector_domain_stage)) then
       deallocate(block_writeback_plan%vector_domain_stage)
    end if
    if (allocated(block_writeback_plan%domain_patch_covered)) then
       deallocate(block_writeback_plan%domain_patch_covered)
    end if

    block_writeback_plan%catalog_size = 0
    block_writeback_plan%domain_count = 0
    block_writeback_plan%installed_block_count = 0
    block_writeback_plan%n_send = 0
    block_writeback_plan%n_recv = 0
    block_writeback_plan%n_retained = 0
    block_writeback_plan%scalar_patch_nvalue = 0
    block_writeback_plan%vector_patch_nvalue = 0
    block_writeback_plan%reconstructed_patch_count = 0
    block_writeback_plan%preserved_patch_count = 0
    block_writeback_plan%production_writeback_count = 0_int64
    block_writeback_plan%ready = .false.

  end subroutine clear_block_writeback_plan


  subroutine build_block_writeback_plan
    ! Build persistent reverse routes from final block owners to the ranks
    ! that own the corresponding legacy Domains. No field is overwritten.

    implicit none

    integer :: b
    integer :: current_block_count
    integer :: d
    integer :: destination
    integer :: first_field_level
    integer :: first
    integer :: ierr
    integer :: last
    integer :: local_index
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_patch
    integer :: n_scalar_variable
    integer :: r
    integer :: slot
    integer :: v_scalar
    integer :: v_vector

    integer, allocatable :: cursor(:)

    integer(int64) :: rank_nvalue
    integer(int64) :: total_patch

    if (.not. allocated(block_catalog) .or. &
         .not. allocated(owner)) then
       call fail("writeback plan before catalogue ownership is ready")
    end if
    if (.not. local_block_store_ready()) then
       call fail("writeback plan before local block installation")
    end if

    current_block_count = n_local_blocks()

    if (block_writeback_plan%ready) then
       if (block_writeback_plan%catalog_size /= size(block_catalog) .or. &
            block_writeback_plan%installed_block_count /= &
            current_block_count) then
          call fail("ready writeback plan has stale block coverage")
       end if
       return
    end if

    call clear_block_writeback_plan

    allocate(block_writeback_plan%send_count(n_process))
    allocate(block_writeback_plan%recv_count(n_process))
    allocate(block_writeback_plan%send_displ(n_process))
    allocate(block_writeback_plan%recv_displ(n_process))

    block_writeback_plan%send_count = 0
    block_writeback_plan%recv_count = 0
    block_writeback_plan%send_displ = 0
    block_writeback_plan%recv_displ = 0
    block_writeback_plan%n_retained = 0

    do local_index = 1,n_local_blocks()
       b = local_block_catalog(local_index)
       if (block_catalog(b)%owner /= rank) then
          call fail("writeback source block has wrong final owner")
       end if
       destination = source_rank(b)
       if (destination == rank) then
          block_writeback_plan%n_retained = &
               block_writeback_plan%n_retained + 1
       else
          block_writeback_plan%send_count(destination+1) = &
               block_writeback_plan%send_count(destination+1) + 1
       end if
    end do

    call MPI_Alltoall( &
         block_writeback_plan%send_count,1,MPI_INTEGER, &
         block_writeback_plan%recv_count,1,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoall writeback block counts")

    do r = 2,n_process
       block_writeback_plan%send_displ(r) = &
            block_writeback_plan%send_displ(r-1) + &
            block_writeback_plan%send_count(r-1)
       block_writeback_plan%recv_displ(r) = &
            block_writeback_plan%recv_displ(r-1) + &
            block_writeback_plan%recv_count(r-1)
    end do

    block_writeback_plan%n_send = &
         sum(block_writeback_plan%send_count)
    block_writeback_plan%n_recv = &
         sum(block_writeback_plan%recv_count)

    allocate(block_writeback_plan%send_block( &
         max(1,block_writeback_plan%n_send)))
    allocate(block_writeback_plan%recv_block( &
         max(1,block_writeback_plan%n_recv)))
    block_writeback_plan%send_block = 0
    block_writeback_plan%recv_block = 0

    allocate(cursor(n_process))
    cursor = block_writeback_plan%send_displ
    do local_index = 1,n_local_blocks()
       b = local_block_catalog(local_index)
       destination = source_rank(b)
       if (destination == rank) cycle
       slot = cursor(destination+1) + 1
       if (slot < 1 .or. slot > block_writeback_plan%n_send) then
          call fail("writeback send-block position is invalid")
       end if
       block_writeback_plan%send_block(slot) = b
       cursor(destination+1) = cursor(destination+1) + 1
    end do
    do r = 1,n_process
       if (cursor(r) /= block_writeback_plan%send_displ(r) + &
            block_writeback_plan%send_count(r)) then
          call fail("writeback send-block count mismatch")
       end if
    end do
    deallocate(cursor)

    call MPI_Alltoallv( &
         block_writeback_plan%send_block, &
         block_writeback_plan%send_count, &
         block_writeback_plan%send_displ,MPI_INTEGER, &
         block_writeback_plan%recv_block, &
         block_writeback_plan%recv_count, &
         block_writeback_plan%recv_displ,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback block manifest")

    allocate(block_writeback_plan%send_patch_count( &
         max(1,block_writeback_plan%n_send)))
    allocate(block_writeback_plan%recv_patch_count( &
         max(1,block_writeback_plan%n_recv)))
    allocate(block_writeback_plan%send_scalar_nvalue( &
         max(1,block_writeback_plan%n_send)))
    allocate(block_writeback_plan%recv_scalar_nvalue( &
         max(1,block_writeback_plan%n_recv)))
    allocate(block_writeback_plan%send_vector_nvalue( &
         max(1,block_writeback_plan%n_send)))
    allocate(block_writeback_plan%recv_vector_nvalue( &
         max(1,block_writeback_plan%n_recv)))

    block_writeback_plan%send_patch_count = 0
    block_writeback_plan%recv_patch_count = 0
    block_writeback_plan%send_scalar_nvalue = 0
    block_writeback_plan%recv_scalar_nvalue = 0
    block_writeback_plan%send_vector_nvalue = 0
    block_writeback_plan%recv_vector_nvalue = 0

    do slot = 1,block_writeback_plan%n_send
       b = block_writeback_plan%send_block(slot)
       n_patch = local_block_patch_count(b)
       block_writeback_plan%send_patch_count(slot) = n_patch
       block_writeback_plan%send_scalar_nvalue(slot) = n_patch * &
            local_block_scalar_family_patch_nvalue(b)
       block_writeback_plan%send_vector_nvalue(slot) = n_patch * &
            local_block_vector_family_patch_nvalue(b)
       if (n_patch <= 0 .or. &
            block_writeback_plan%send_scalar_nvalue(slot) <= 0 .or. &
            block_writeback_plan%send_vector_nvalue(slot) <= 0) then
          call fail("writeback send payload extent is invalid")
       end if
    end do

    call MPI_Alltoallv( &
         block_writeback_plan%send_patch_count, &
         block_writeback_plan%send_count, &
         block_writeback_plan%send_displ,MPI_INTEGER, &
         block_writeback_plan%recv_patch_count, &
         block_writeback_plan%recv_count, &
         block_writeback_plan%recv_displ,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback patch counts")
    call MPI_Alltoallv( &
         block_writeback_plan%send_scalar_nvalue, &
         block_writeback_plan%send_count, &
         block_writeback_plan%send_displ,MPI_INTEGER, &
         block_writeback_plan%recv_scalar_nvalue, &
         block_writeback_plan%recv_count, &
         block_writeback_plan%recv_displ,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback scalar extents")
    call MPI_Alltoallv( &
         block_writeback_plan%send_vector_nvalue, &
         block_writeback_plan%send_count, &
         block_writeback_plan%send_displ,MPI_INTEGER, &
         block_writeback_plan%recv_vector_nvalue, &
         block_writeback_plan%recv_count, &
         block_writeback_plan%recv_displ,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback vector extents")

    allocate(block_writeback_plan%scalar_send_count(n_process))
    allocate(block_writeback_plan%scalar_recv_count(n_process))
    allocate(block_writeback_plan%scalar_send_displ(n_process))
    allocate(block_writeback_plan%scalar_recv_displ(n_process))
    allocate(block_writeback_plan%vector_send_count(n_process))
    allocate(block_writeback_plan%vector_recv_count(n_process))
    allocate(block_writeback_plan%vector_send_displ(n_process))
    allocate(block_writeback_plan%vector_recv_displ(n_process))

    block_writeback_plan%scalar_send_count = 0
    block_writeback_plan%scalar_recv_count = 0
    block_writeback_plan%scalar_send_displ = 0
    block_writeback_plan%scalar_recv_displ = 0
    block_writeback_plan%vector_send_count = 0
    block_writeback_plan%vector_recv_count = 0
    block_writeback_plan%vector_send_displ = 0
    block_writeback_plan%vector_recv_displ = 0

    do r = 1,n_process
       first = block_writeback_plan%send_displ(r) + 1
       last = block_writeback_plan%send_displ(r) + &
            block_writeback_plan%send_count(r)
       rank_nvalue = 0_int64
       if (last >= first) then
          rank_nvalue = sum(int( &
               block_writeback_plan%send_scalar_nvalue(first:last), &
               int64))
       end if
       if (rank_nvalue > int(huge(0),int64)) then
          call fail("writeback scalar send count exceeds MPI range")
       end if
       block_writeback_plan%scalar_send_count(r) = int(rank_nvalue)

       rank_nvalue = 0_int64
       if (last >= first) then
          rank_nvalue = sum(int( &
               block_writeback_plan%send_vector_nvalue(first:last), &
               int64))
       end if
       if (rank_nvalue > int(huge(0),int64)) then
          call fail("writeback vector send count exceeds MPI range")
       end if
       block_writeback_plan%vector_send_count(r) = int(rank_nvalue)

       first = block_writeback_plan%recv_displ(r) + 1
       last = block_writeback_plan%recv_displ(r) + &
            block_writeback_plan%recv_count(r)
       rank_nvalue = 0_int64
       if (last >= first) then
          rank_nvalue = sum(int( &
               block_writeback_plan%recv_scalar_nvalue(first:last), &
               int64))
       end if
       if (rank_nvalue > int(huge(0),int64)) then
          call fail("writeback scalar receive count exceeds MPI range")
       end if
       block_writeback_plan%scalar_recv_count(r) = int(rank_nvalue)

       rank_nvalue = 0_int64
       if (last >= first) then
          rank_nvalue = sum(int( &
               block_writeback_plan%recv_vector_nvalue(first:last), &
               int64))
       end if
       if (rank_nvalue > int(huge(0),int64)) then
          call fail("writeback vector receive count exceeds MPI range")
       end if
       block_writeback_plan%vector_recv_count(r) = int(rank_nvalue)
    end do

    do r = 2,n_process
       block_writeback_plan%scalar_send_displ(r) = &
            block_writeback_plan%scalar_send_displ(r-1) + &
            block_writeback_plan%scalar_send_count(r-1)
       block_writeback_plan%scalar_recv_displ(r) = &
            block_writeback_plan%scalar_recv_displ(r-1) + &
            block_writeback_plan%scalar_recv_count(r-1)
       block_writeback_plan%vector_send_displ(r) = &
            block_writeback_plan%vector_send_displ(r-1) + &
            block_writeback_plan%vector_send_count(r-1)
       block_writeback_plan%vector_recv_displ(r) = &
            block_writeback_plan%vector_recv_displ(r-1) + &
            block_writeback_plan%vector_recv_count(r-1)
    end do

    allocate(block_writeback_plan%scalar_send_buffer(max(1, &
         sum(block_writeback_plan%scalar_send_count))))
    allocate(block_writeback_plan%scalar_recv_buffer(max(1, &
         sum(block_writeback_plan%scalar_recv_count))))
    allocate(block_writeback_plan%vector_send_buffer(max(1, &
         sum(block_writeback_plan%vector_send_count))))
    allocate(block_writeback_plan%vector_recv_buffer(max(1, &
         sum(block_writeback_plan%vector_recv_count))))
    block_writeback_plan%buffer_allocations = &
         block_writeback_plan%buffer_allocations + 4_int64

    block_writeback_plan%scalar_send_buffer = 0.0_dp
    block_writeback_plan%scalar_recv_buffer = 0.0_dp
    block_writeback_plan%vector_send_buffer = 0.0_dp
    block_writeback_plan%vector_recv_buffer = 0.0_dp

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,first_field_level, &
         n_field_level,mult_scalar,mult_vector)
    if (v_scalar < 1 .or. v_vector < 1 .or. &
         n_scalar_variable < 1 .or. n_field_level < 1 .or. &
         mult_scalar < 1 .or. mult_vector < 1) then
       call fail("writeback Domain stage field layout is invalid")
    end if
    block_writeback_plan%scalar_patch_nvalue = &
         n_scalar_variable*n_field_level*mult_scalar*PATCH_SIZE**2
    block_writeback_plan%vector_patch_nvalue = &
         n_field_level*mult_vector*PATCH_SIZE**2

    allocate(block_writeback_plan%domain_patch_displ(size(grid)))
    block_writeback_plan%domain_patch_displ = 0
    do d = 2,size(grid)
       block_writeback_plan%domain_patch_displ(d) = &
            block_writeback_plan%domain_patch_displ(d-1) + &
            grid(d-1)%patch%length
    end do

    total_patch = 0_int64
    if (size(grid) > 0) then
       total_patch = int( &
            block_writeback_plan%domain_patch_displ(size(grid)),int64) + &
            int(grid(size(grid))%patch%length,int64)
    end if
    if (total_patch > int(huge(0),int64)) then
       call fail("writeback Domain stage patch count exceeds range")
    end if
    if (total_patch*int( &
         block_writeback_plan%scalar_patch_nvalue,int64) > &
         int(huge(0),int64) .or. &
         total_patch*int( &
         block_writeback_plan%vector_patch_nvalue,int64) > &
         int(huge(0),int64)) then
       call fail("writeback Domain stage value count exceeds range")
    end if

    allocate(block_writeback_plan%scalar_domain_stage(max(1, &
         int(total_patch)*block_writeback_plan%scalar_patch_nvalue)))
    allocate(block_writeback_plan%vector_domain_stage(max(1, &
         int(total_patch)*block_writeback_plan%vector_patch_nvalue)))
    allocate(block_writeback_plan%domain_patch_covered(max(1, &
         int(total_patch))))
    block_writeback_plan%stage_allocations = &
         block_writeback_plan%stage_allocations + 3_int64
    block_writeback_plan%scalar_domain_stage = 0.0_dp
    block_writeback_plan%vector_domain_stage = 0.0_dp
    block_writeback_plan%domain_patch_covered = .false.
    block_writeback_plan%catalog_size = size(block_catalog)
    block_writeback_plan%domain_count = size(grid)
    block_writeback_plan%installed_block_count = n_local_blocks()
    block_writeback_plan%ready = .true.

  end subroutine build_block_writeback_plan


  logical function block_writeback_plan_is_ready () result(ready)
    ! Report whether the reverse Domain-owner routes match this block store.

    implicit none

    ready = .false.
    if (.not. block_writeback_plan%ready) return
    if (.not. allocated(block_catalog)) return
    if (.not. local_block_store_ready()) return
    if (block_writeback_plan%catalog_size /= size(block_catalog)) return
    if (block_writeback_plan%domain_count /= size(grid)) return
    if (block_writeback_plan%installed_block_count /= &
         n_local_blocks()) return
    ready = .true.

  end function block_writeback_plan_is_ready


  integer(int64) function block_writeback_plan_allocation_count () &
       result(n_allocation)
    ! Count persistent scalar/vector communication-buffer allocations.

    implicit none

    n_allocation = block_writeback_plan%buffer_allocations

  end function block_writeback_plan_allocation_count


  integer(int64) function block_domain_production_writeback_count () &
       result(n_writeback)
    ! Count successful production-facing field-family synchronizations.

    implicit none

    n_writeback = block_writeback_plan%production_writeback_count

  end function block_domain_production_writeback_count


  subroutine check_block_writeback_plan (verbose)
    ! Prove that every final-owner block has exactly one reverse route to
    ! its Domain owner and that all future payload extents balance globally.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: destination
    integer :: expected_remote
    integer :: expected_retained
    integer :: ierr
    integer :: pos
    integer :: r
    integer :: stage_patch_count

    integer(int64) :: allocation_before
    integer(int64) :: expected_checksum
    integer(int64) :: global_recv
    integer(int64) :: global_recv_checksum
    integer(int64) :: global_recv_patch
    integer(int64) :: global_recv_scalar
    integer(int64) :: global_recv_vector
    integer(int64) :: global_retained
    integer(int64) :: global_send
    integer(int64) :: global_send_checksum
    integer(int64) :: global_send_patch
    integer(int64) :: global_send_scalar
    integer(int64) :: global_send_vector
    integer(int64) :: local_recv
    integer(int64) :: local_recv_checksum
    integer(int64) :: local_recv_patch
    integer(int64) :: local_recv_scalar
    integer(int64) :: local_recv_vector
    integer(int64) :: local_retained
    integer(int64) :: local_send
    integer(int64) :: local_send_checksum
    integer(int64) :: local_send_patch
    integer(int64) :: local_send_scalar
    integer(int64) :: local_send_vector

    logical :: print_summary
    logical, allocatable :: seen(:)

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. block_writeback_plan_is_ready()) then
       call fail("writeback plan check before plan is ready")
    end if
    if (.not. allocated(block_writeback_plan%send_block) .or. &
         .not. allocated(block_writeback_plan%recv_block) .or. &
         .not. allocated(block_writeback_plan%send_patch_count) .or. &
         .not. allocated(block_writeback_plan%recv_patch_count) .or. &
         .not. allocated(block_writeback_plan%send_scalar_nvalue) .or. &
         .not. allocated(block_writeback_plan%recv_scalar_nvalue) .or. &
         .not. allocated(block_writeback_plan%send_vector_nvalue) .or. &
         .not. allocated(block_writeback_plan%recv_vector_nvalue)) then
       call fail("writeback plan record storage is incomplete")
    end if
    if (.not. allocated(block_writeback_plan%scalar_send_buffer) .or. &
         .not. allocated(block_writeback_plan%scalar_recv_buffer) .or. &
         .not. allocated(block_writeback_plan%vector_send_buffer) .or. &
         .not. allocated(block_writeback_plan%vector_recv_buffer)) then
       call fail("writeback plan payload storage is incomplete")
    end if
    if (.not. allocated(block_writeback_plan%domain_patch_displ) .or. &
         .not. allocated(block_writeback_plan%scalar_domain_stage) .or. &
         .not. allocated(block_writeback_plan%vector_domain_stage) .or. &
         .not. allocated(block_writeback_plan%domain_patch_covered)) then
       call fail("writeback Domain stage storage is incomplete")
    end if

    if (block_writeback_plan%n_send /= &
         sum(block_writeback_plan%send_count) .or. &
         block_writeback_plan%n_recv /= &
         sum(block_writeback_plan%recv_count)) then
       call fail("writeback plan block totals are inconsistent")
    end if
    if (block_writeback_plan%n_send + &
         block_writeback_plan%n_retained /= n_local_blocks()) then
       call fail("writeback plan does not cover every local block")
    end if
    if (size(block_writeback_plan%scalar_send_buffer) /= max(1, &
         sum(block_writeback_plan%scalar_send_count)) .or. &
         size(block_writeback_plan%scalar_recv_buffer) /= max(1, &
         sum(block_writeback_plan%scalar_recv_count)) .or. &
         size(block_writeback_plan%vector_send_buffer) /= max(1, &
         sum(block_writeback_plan%vector_send_count)) .or. &
         size(block_writeback_plan%vector_recv_buffer) /= max(1, &
         sum(block_writeback_plan%vector_recv_count))) then
       call fail("writeback plan payload buffer extent mismatch")
    end if

    stage_patch_count = 0
    if (size(grid) > 0) then
       stage_patch_count = &
            block_writeback_plan%domain_patch_displ(size(grid)) + &
            grid(size(grid))%patch%length
    end if
    if (size(block_writeback_plan%domain_patch_covered) /= &
         max(1,stage_patch_count) .or. &
         size(block_writeback_plan%scalar_domain_stage) /= max(1, &
         stage_patch_count* &
         block_writeback_plan%scalar_patch_nvalue) .or. &
         size(block_writeback_plan%vector_domain_stage) /= max(1, &
         stage_patch_count* &
         block_writeback_plan%vector_patch_nvalue)) then
       call fail("writeback Domain stage extent mismatch")
    end if

    allocate(seen(size(block_catalog)))
    seen = .false.

    do r = 0,n_process-1
       do pos = block_writeback_plan%send_displ(r+1)+1, &
            block_writeback_plan%send_displ(r+1) + &
            block_writeback_plan%send_count(r+1)
          b = block_writeback_plan%send_block(pos)
          if (b < 1 .or. b > size(block_catalog)) then
             call fail("writeback send catalogue index is invalid")
          end if
          if (seen(b)) then
             call fail("duplicate block in writeback send manifest")
          end if
          seen(b) = .true.
          destination = source_rank(b)
          if (block_catalog(b)%owner /= rank .or. &
               destination /= r .or. r == rank) then
             call fail("writeback send route is invalid")
          end if
          if (catalog_local_block(b) < 1) then
             call fail("writeback send block is not installed locally")
          end if
          if (block_writeback_plan%send_patch_count(pos) <= 0 .or. &
               block_writeback_plan%send_scalar_nvalue(pos) <= 0 .or. &
               block_writeback_plan%send_vector_nvalue(pos) <= 0) then
             call fail("writeback send record has invalid extent")
          end if
       end do
    end do

    seen = .false.
    do r = 0,n_process-1
       do pos = block_writeback_plan%recv_displ(r+1)+1, &
            block_writeback_plan%recv_displ(r+1) + &
            block_writeback_plan%recv_count(r+1)
          b = block_writeback_plan%recv_block(pos)
          if (b < 1 .or. b > size(block_catalog)) then
             call fail("writeback receive catalogue index is invalid")
          end if
          if (seen(b)) then
             call fail("duplicate block in writeback receive manifest")
          end if
          seen(b) = .true.
          destination = source_rank(b)
          if (block_catalog(b)%owner /= r .or. &
               destination /= rank .or. r == rank) then
             call fail("writeback receive route is invalid")
          end if
          if (block_writeback_plan%recv_patch_count(pos) <= 0 .or. &
               block_writeback_plan%recv_scalar_nvalue(pos) <= 0 .or. &
               block_writeback_plan%recv_vector_nvalue(pos) <= 0) then
             call fail("writeback receive record has invalid extent")
          end if
       end do
    end do

    expected_retained = 0
    do b = 1,size(block_catalog)
       if (source_rank(b) /= rank) cycle
       if (block_catalog(b)%owner == rank) then
          expected_retained = expected_retained + 1
          if (catalog_local_block(b) < 1) then
             call fail("retained writeback block is not installed")
          end if
       else
          if (.not. seen(b)) then
             call fail("Domain-owner writeback block is missing")
          end if
       end if
    end do
    if (block_writeback_plan%n_retained /= expected_retained) then
       call fail("writeback retained-block count mismatch")
    end if
    deallocate(seen)

    allocation_before = block_writeback_plan_allocation_count()
    call build_block_writeback_plan
    if (block_writeback_plan_allocation_count() /= allocation_before) then
       call fail("ready writeback plan reallocated payload buffers")
    end if

    local_send_checksum = 0_int64
    local_recv_checksum = 0_int64
    if (block_writeback_plan%n_send > 0) then
       local_send_checksum = sum(int( &
            block_writeback_plan%send_block( &
            1:block_writeback_plan%n_send),int64))
    end if
    if (block_writeback_plan%n_recv > 0) then
       local_recv_checksum = sum(int( &
            block_writeback_plan%recv_block( &
            1:block_writeback_plan%n_recv),int64))
    end if

    local_send = int(block_writeback_plan%n_send,int64)
    local_recv = int(block_writeback_plan%n_recv,int64)
    local_retained = int(block_writeback_plan%n_retained,int64)
    local_send_patch = 0_int64
    local_recv_patch = 0_int64
    if (block_writeback_plan%n_send > 0) then
       local_send_patch = sum(int( &
            block_writeback_plan%send_patch_count( &
            1:block_writeback_plan%n_send),int64))
    end if
    if (block_writeback_plan%n_recv > 0) then
       local_recv_patch = sum(int( &
            block_writeback_plan%recv_patch_count( &
            1:block_writeback_plan%n_recv),int64))
    end if
    local_send_scalar = sum(int( &
         block_writeback_plan%scalar_send_count,int64))
    local_recv_scalar = sum(int( &
         block_writeback_plan%scalar_recv_count,int64))
    local_send_vector = sum(int( &
         block_writeback_plan%vector_send_count,int64))
    local_recv_vector = sum(int( &
         block_writeback_plan%vector_recv_count,int64))

    call MPI_Allreduce(local_send,global_send,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback send blocks")
    call MPI_Allreduce(local_recv,global_recv,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback receive blocks")
    call MPI_Allreduce(local_retained,global_retained,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback retained blocks")
    call MPI_Allreduce(local_send_checksum,global_send_checksum,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback send checksum")
    call MPI_Allreduce(local_recv_checksum,global_recv_checksum,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback receive checksum")
    call MPI_Allreduce(local_send_patch,global_send_patch,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback send patches")
    call MPI_Allreduce(local_recv_patch,global_recv_patch,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback receive patches")
    call MPI_Allreduce(local_send_scalar,global_send_scalar,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback send scalar values")
    call MPI_Allreduce(local_recv_scalar,global_recv_scalar,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback receive scalar values")
    call MPI_Allreduce(local_send_vector,global_send_vector,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback send vector values")
    call MPI_Allreduce(local_recv_vector,global_recv_vector,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writeback receive vector values")

    expected_remote = 0
    expected_checksum = 0_int64
    do b = 1,size(block_catalog)
       if (block_catalog(b)%owner == source_rank(b)) cycle
       expected_remote = expected_remote + 1
       expected_checksum = expected_checksum + int(b,int64)
    end do

    if (global_send /= int(expected_remote,int64) .or. &
         global_recv /= int(expected_remote,int64) .or. &
         global_retained+global_send /= &
         int(size(block_catalog),int64)) then
       call fail("global writeback block coverage mismatch")
    end if
    if (global_send_checksum /= expected_checksum .or. &
         global_recv_checksum /= expected_checksum) then
       call fail("global writeback catalogue checksum mismatch")
    end if
    if (global_send_patch /= global_recv_patch .or. &
         global_send_scalar /= global_recv_scalar .or. &
         global_send_vector /= global_recv_vector) then
       call fail("global writeback payload extent mismatch")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Persistent block writeback plan for rank ",rank,":"
       write(6,'(a,i0)') "  remote blocks sent     = ", &
            block_writeback_plan%n_send
       write(6,'(a,i0)') "  remote blocks received = ", &
            block_writeback_plan%n_recv
       write(6,'(a,i0)') "  locally retained blocks = ", &
            block_writeback_plan%n_retained
       write(6,'(a,i0)') "  scalar send values      = ", &
            sum(block_writeback_plan%scalar_send_count)
       write(6,'(a,i0)') "  vector send values      = ", &
            sum(block_writeback_plan%vector_send_count)
       write(6,'(a)') "  final-owner to Domain-owner routes passed"
       write(6,'(a)') "  patch and payload extent exchange passed"
       write(6,'(a,/)') "  persistent buffer reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global remote block writeback routes = ",global_send
       write(6,'(a,i0)') &
            "Global retained block writeback routes = ",global_retained
       write(6,'(a,i0)') &
            "Global writeback patches = ",global_send_patch
       write(6,'(a,i0)') &
            "Global writeback scalar values = ",global_send_scalar
       write(6,'(a,i0)') &
            "Global writeback vector values = ",global_send_vector
       write(6,'(a,/)') &
            "Persistent block-to-Domain writeback plan passed"
    end if

  end subroutine check_block_writeback_plan


  subroutine exchange_block_writeback_payloads (payload_family)
    ! Pack one scalar/vector prognostic-field family owned by final block
    ! ranks and transport it back to the legacy Domain owners.
    ! The persistent plan buffers are reused and Domain fields are unchanged.

    implicit none

    integer, optional, intent(in) :: payload_family

    integer :: b
    integer :: family
    integer :: ierr
    integer :: local_patch
    integer :: n_patch
    integer :: n_scalar_patch
    integer :: n_vector_patch
    integer :: pos_scalar
    integer :: pos_vector
    integer :: r
    integer :: slot

    if (.not. block_writeback_plan_is_ready()) then
       call fail("writeback payload exchange before plan is ready")
    end if

    family = BLOCK_PAYLOAD_SOL
    if (present(payload_family)) family = payload_family
    if (family /= BLOCK_PAYLOAD_SOL .and. &
         family /= BLOCK_PAYLOAD_WAV_COEFF) then
       call fail("invalid writeback payload family")
    end if

    block_writeback_plan%scalar_send_buffer = 0.0_dp
    block_writeback_plan%scalar_recv_buffer = 0.0_dp
    block_writeback_plan%vector_send_buffer = 0.0_dp
    block_writeback_plan%vector_recv_buffer = 0.0_dp

    do r = 1,n_process
       pos_scalar = block_writeback_plan%scalar_send_displ(r) + 1
       pos_vector = block_writeback_plan%vector_send_displ(r) + 1

       do slot = block_writeback_plan%send_displ(r)+1, &
            block_writeback_plan%send_displ(r) + &
            block_writeback_plan%send_count(r)
          b = block_writeback_plan%send_block(slot)
          n_patch = local_block_patch_count(b)
          n_scalar_patch = local_block_scalar_family_patch_nvalue(b)
          n_vector_patch = local_block_vector_family_patch_nvalue(b)

          if (n_patch /= block_writeback_plan%send_patch_count(slot) .or. &
               n_patch*n_scalar_patch /= &
               block_writeback_plan%send_scalar_nvalue(slot) .or. &
               n_patch*n_vector_patch /= &
               block_writeback_plan%send_vector_nvalue(slot)) then
             call fail("writeback payload source extent changed")
          end if

          do local_patch = 0,n_patch-1
             call get_local_block_scalar_patch_family_values( &
                  b,local_patch,family, &
                  block_writeback_plan%scalar_send_buffer( &
                  pos_scalar:pos_scalar+n_scalar_patch-1))
             call get_local_block_vector_patch_family_values( &
                  b,local_patch,family, &
                  block_writeback_plan%vector_send_buffer( &
                  pos_vector:pos_vector+n_vector_patch-1))
             pos_scalar = pos_scalar + n_scalar_patch
             pos_vector = pos_vector + n_vector_patch
          end do
       end do

       if (pos_scalar /= block_writeback_plan%scalar_send_displ(r) + &
            block_writeback_plan%scalar_send_count(r) + 1 .or. &
            pos_vector /= block_writeback_plan%vector_send_displ(r) + &
            block_writeback_plan%vector_send_count(r) + 1) then
          call fail("writeback packed send extent mismatch")
       end if
    end do

    call MPI_Alltoallv( &
         block_writeback_plan%scalar_send_buffer, &
         block_writeback_plan%scalar_send_count, &
         block_writeback_plan%scalar_send_displ, &
         MPI_DOUBLE_PRECISION, &
         block_writeback_plan%scalar_recv_buffer, &
         block_writeback_plan%scalar_recv_count, &
         block_writeback_plan%scalar_recv_displ, &
         MPI_DOUBLE_PRECISION,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback scalar-family payload")

    call MPI_Alltoallv( &
         block_writeback_plan%vector_send_buffer, &
         block_writeback_plan%vector_send_count, &
         block_writeback_plan%vector_send_displ, &
         MPI_DOUBLE_PRECISION, &
         block_writeback_plan%vector_recv_buffer, &
         block_writeback_plan%vector_recv_count, &
         block_writeback_plan%vector_recv_displ, &
         MPI_DOUBLE_PRECISION,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv writeback vector-family payload")

  end subroutine exchange_block_writeback_payloads


  subroutine check_block_writeback_payload_exchange (verbose)
    ! Compare every remotely transported sol and wav_coeff value with the
    ! exact values still held by its destination Domain owner. Repeat each
    ! exchange to prove that the persistent communication storage is reused
    ! without changing the authoritative Domain representation.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: d
    integer :: ierr
    integer :: n_patch
    integer :: pos_scalar
    integer :: pos_vector
    integer :: r
    integer :: scalar_start
    integer :: slot
    integer :: vector_start

    integer(int64) :: allocation_before
    integer(int64) :: global_scalar
    integer(int64) :: global_vector
    integer(int64) :: local_scalar
    integer(int64) :: local_vector

    logical :: print_summary

    real(dp), allocatable :: expected_scalar(:)
    real(dp), allocatable :: expected_vector(:)

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. block_writeback_plan_is_ready()) then
       call fail("writeback payload check before plan is ready")
    end if

    allocate(expected_scalar(size( &
         block_writeback_plan%scalar_recv_buffer)))
    allocate(expected_vector(size( &
         block_writeback_plan%vector_recv_buffer)))

    allocation_before = block_writeback_plan_allocation_count()

    call check_writeback_family(BLOCK_PAYLOAD_SOL)
    call check_writeback_family(BLOCK_PAYLOAD_WAV_COEFF)

    if (block_writeback_plan_allocation_count() /= allocation_before) then
       call fail("writeback payload exchange reallocated plan buffers")
    end if

    local_scalar = int(sum( &
         block_writeback_plan%scalar_recv_count),int64)
    local_vector = int(sum( &
         block_writeback_plan%vector_recv_count),int64)

    call MPI_Allreduce(local_scalar,global_scalar,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce transported scalar values")
    call MPI_Allreduce(local_vector,global_vector,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce transported vector values")

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Block writeback payload exchange for rank ",rank,":"
       write(6,'(a,i0)') "  scalar values received = ",local_scalar
       write(6,'(a,i0)') "  vector values received = ",local_vector
       write(6,'(a)') "  exact sol Domain payload comparison passed"
       write(6,'(a)') "  exact wav_coeff Domain payload comparison passed"
       write(6,'(a)') "  repeated sol exchange passed"
       write(6,'(a,/)') "  repeated wav_coeff exchange passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global transported writeback scalar values = ", &
            global_scalar
       write(6,'(a,i0)') &
            "Global transported writeback vector values = ", &
            global_vector
       write(6,'(a,/)') &
            "Exact sol/wav_coeff block-to-Domain transport passed"
    end if

    deallocate(expected_vector)
    deallocate(expected_scalar)

  contains

    subroutine check_writeback_family (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      call exchange_block_writeback_payloads(payload_family)
      call build_expected_writeback_payloads( &
           payload_family,expected_scalar,expected_vector)
      call compare_writeback_payloads( &
           payload_family,expected_scalar,expected_vector)

      call exchange_block_writeback_payloads(payload_family)
      call build_expected_writeback_payloads( &
           payload_family,expected_scalar,expected_vector)
      call compare_writeback_payloads( &
           payload_family,expected_scalar,expected_vector)

    end subroutine check_writeback_family

    subroutine build_expected_writeback_payloads ( &
         payload_family,scalar_payload,vector_payload)

      implicit none

      integer, intent(in) :: payload_family
      real(dp), intent(out) :: scalar_payload(:)
      real(dp), intent(out) :: vector_payload(:)

      scalar_payload = 0.0_dp
      vector_payload = 0.0_dp

      do r = 1,n_process
         pos_scalar = block_writeback_plan%scalar_recv_displ(r) + 1
         pos_vector = block_writeback_plan%vector_recv_displ(r) + 1

         do slot = block_writeback_plan%recv_displ(r)+1, &
              block_writeback_plan%recv_displ(r) + &
              block_writeback_plan%recv_count(r)
            b = block_writeback_plan%recv_block(slot)
            if (source_rank(b) /= rank) then
               call fail("writeback payload arrived at wrong Domain owner")
            end if
            if (block_catalog(b)%owner /= r-1) then
               call fail("writeback payload came from wrong block owner")
            end if

            d = loc_id(block_catalog(b)%root_domain+1) + 1
            if (d < 1 .or. d > size(grid)) then
               call fail("writeback payload has invalid local Domain")
            end if

            scalar_start = pos_scalar
            vector_start = pos_vector
            n_patch = 0
            call pack_domain_subtree_prognostic( &
                 d,block_catalog(b)%root_patch,payload_family, &
                 scalar_payload,pos_scalar, &
                 vector_payload,pos_vector,n_patch)

            if (n_patch /= &
                 block_writeback_plan%recv_patch_count(slot) .or. &
                 pos_scalar-scalar_start /= &
                 block_writeback_plan%recv_scalar_nvalue(slot) .or. &
                 pos_vector-vector_start /= &
                 block_writeback_plan%recv_vector_nvalue(slot)) then
               call fail("writeback expected payload extent mismatch")
            end if
         end do

         if (pos_scalar /= block_writeback_plan%scalar_recv_displ(r) + &
              block_writeback_plan%scalar_recv_count(r) + 1 .or. &
              pos_vector /= block_writeback_plan%vector_recv_displ(r) + &
              block_writeback_plan%vector_recv_count(r) + 1) then
            call fail("writeback expected receive extent mismatch")
         end if
      end do

    end subroutine build_expected_writeback_payloads


    subroutine compare_writeback_payloads ( &
         payload_family,scalar_payload,vector_payload)

      implicit none

      integer, intent(in) :: payload_family
      real(dp), intent(in) :: scalar_payload(:)
      real(dp), intent(in) :: vector_payload(:)

      integer :: n_scalar
      integer :: n_vector

      n_scalar = sum(block_writeback_plan%scalar_recv_count)
      n_vector = sum(block_writeback_plan%vector_recv_count)

      if (n_scalar > 0) then
         if (any(abs(block_writeback_plan%scalar_recv_buffer( &
              1:n_scalar)-scalar_payload(1:n_scalar)) > 0.0_dp)) then
            if (payload_family == BLOCK_PAYLOAD_SOL) then
               call fail("transported scalar sol writeback mismatch")
            else
               call fail("transported scalar wav_coeff writeback mismatch")
            end if
         end if
      end if
      if (n_vector > 0) then
         if (any(abs(block_writeback_plan%vector_recv_buffer( &
              1:n_vector)-vector_payload(1:n_vector)) > 0.0_dp)) then
            if (payload_family == BLOCK_PAYLOAD_SOL) then
               call fail("transported vector sol writeback mismatch")
            else
               call fail("transported vector wav_coeff writeback mismatch")
            end if
         end if
      end if

    end subroutine compare_writeback_payloads

  end subroutine check_block_writeback_payload_exchange


  recursive subroutine pack_domain_subtree_prognostic ( &
       d,p,payload_family,scalar_payload,scalar_pos, &
       vector_payload,vector_pos,n_patch)
    ! Pack one legacy Domain subtree in the same preorder and field order
    ! used by the compact block field-family patch accessors.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: payload_family
    real(dp), intent(inout) :: scalar_payload(:)
    integer, intent(inout) :: scalar_pos
    real(dp), intent(inout) :: vector_payload(:)
    integer, intent(inout) :: vector_pos
    integer, intent(inout) :: n_patch

    integer :: c
    integer :: p_child

    if (p < 0 .or. p >= grid(d)%patch%length) then
       call fail("invalid patch in writeback Domain payload")
    end if
    if (grid(d)%patch%elts(p+1)%deleted) return

    call pack_domain_patch_prognostic( &
         d,p,payload_family,scalar_payload,scalar_pos, &
         vector_payload,vector_pos)

    n_patch = n_patch + 1

    do c = 1,N_CHDRN
       p_child = grid(d)%patch%elts(p+1)%children(c)
       if (p_child <= 0) cycle
       call pack_domain_subtree_prognostic( &
            d,p_child,payload_family,scalar_payload,scalar_pos, &
            vector_payload,vector_pos,n_patch)
    end do

  end subroutine pack_domain_subtree_prognostic


  subroutine pack_domain_patch_prognostic ( &
       d,p,payload_family,scalar_payload,scalar_pos, &
       vector_payload,vector_pos)
    ! Pack one active legacy Domain patch in compact field-family order.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: payload_family
    real(dp), intent(inout) :: scalar_payload(:)
    integer, intent(inout) :: scalar_pos
    real(dp), intent(inout) :: vector_payload(:)
    integer, intent(inout) :: vector_pos

    integer :: field_level
    integer :: first_field_level
    integer :: level_slot
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_scalar_patch
    integer :: n_scalar_variable
    integer :: n_value
    integer :: n_vector_patch
    integer :: scalar_id
    integer :: scalar_slot
    integer :: start
    integer :: v_scalar
    integer :: v_vector

    if (d < 1 .or. d > size(grid)) then
       call fail("invalid Domain in writeback record")
    end if
    if (p < 0 .or. p >= grid(d)%patch%length) then
       call fail("invalid patch in Domain writeback record")
    end if
    if (grid(d)%patch%elts(p+1)%deleted) then
       call fail("deleted patch in Domain writeback record")
    end if
    if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
         payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
       call fail("invalid Domain writeback payload family")
    end if

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,first_field_level, &
         n_field_level,mult_scalar,mult_vector)

    n_scalar_patch = &
         n_scalar_variable*n_field_level*mult_scalar*PATCH_SIZE**2
    n_vector_patch = n_field_level*mult_vector*PATCH_SIZE**2
    if (scalar_pos < 1 .or. &
         scalar_pos+n_scalar_patch-1 > size(scalar_payload) .or. &
         vector_pos < 1 .or. &
         vector_pos+n_vector_patch-1 > size(vector_payload)) then
       call fail("Domain writeback record buffer extent is invalid")
    end if

    start = mult_scalar*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_scalar*PATCH_SIZE**2

    do scalar_slot = 1,n_scalar_variable
       scalar_id = v_scalar + scalar_slot - 1
       do level_slot = 1,n_field_level
          field_level = first_field_level + level_slot - 1
          select case (payload_family)
          case (BLOCK_PAYLOAD_SOL)
             if (start < 0 .or. start+n_value > &
                  size(sol(scalar_id,field_level)%data(d)%elts)) then
                call fail("legacy scalar sol writeback extent is invalid")
             end if
             scalar_payload(scalar_pos:scalar_pos+n_value-1) = &
                  sol(scalar_id,field_level)%data(d)%elts( &
                  start+1:start+n_value)
          case (BLOCK_PAYLOAD_WAV_COEFF)
             if (start < 0 .or. start+n_value > &
                  size(wav_coeff(scalar_id,field_level)%data(d)%elts)) then
                call fail( &
                     "legacy scalar wav_coeff writeback extent is invalid")
             end if
             scalar_payload(scalar_pos:scalar_pos+n_value-1) = &
                  wav_coeff(scalar_id,field_level)%data(d)%elts( &
                  start+1:start+n_value)
          end select
          scalar_pos = scalar_pos + n_value
       end do
    end do

    start = mult_vector*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_vector*PATCH_SIZE**2

    do level_slot = 1,n_field_level
       field_level = first_field_level + level_slot - 1
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (start < 0 .or. start+n_value > &
               size(sol(v_vector,field_level)%data(d)%elts)) then
             call fail("legacy vector sol writeback extent is invalid")
          end if
          vector_payload(vector_pos:vector_pos+n_value-1) = &
               sol(v_vector,field_level)%data(d)%elts( &
               start+1:start+n_value)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (start < 0 .or. start+n_value > &
               size(wav_coeff(v_vector,field_level)%data(d)%elts)) then
             call fail( &
                  "legacy vector wav_coeff writeback extent is invalid")
          end if
          vector_payload(vector_pos:vector_pos+n_value-1) = &
               wav_coeff(v_vector,field_level)%data(d)%elts( &
               start+1:start+n_value)
       end select
       vector_pos = vector_pos + n_value
    end do

  end subroutine pack_domain_patch_prognostic


  subroutine reconstruct_block_writeback_domain_stage (payload_family)
    ! Reconstruct every active patch of each locally owned legacy Domain in
    ! persistent non-authoritative staging arrays. Retained blocks are read
    ! directly, remote blocks are consumed from the writeback buffers and
    ! fixed coarse-scaffold patches outside the catalogue are preserved.

    implicit none

    integer, intent(in) :: payload_family

    integer :: active_patch_count
    integer :: b
    integer :: d
    integer :: destination
    integer :: expected_patch_count
    integer :: local_index
    integer :: local_patch
    integer :: n_patch
    integer :: p
    integer :: patch_slot
    integer :: pos_scalar
    integer :: pos_vector
    integer :: preserved_patch_count
    integer :: r
    integer :: reconstructed_patch_count
    integer :: scalar_limit
    integer :: scalar_start
    integer :: slot
    integer :: vector_limit
    integer :: vector_stage_start
    integer :: vector_start
    integer :: scalar_stage_start

    if (.not. block_writeback_plan_is_ready()) then
       call fail("Domain reconstruction before writeback plan is ready")
    end if
    if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
         payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
       call fail("invalid Domain reconstruction payload family")
    end if

    call exchange_block_writeback_payloads(payload_family)

    block_writeback_plan%scalar_domain_stage = 0.0_dp
    block_writeback_plan%vector_domain_stage = 0.0_dp
    block_writeback_plan%domain_patch_covered = .false.
    reconstructed_patch_count = 0
    preserved_patch_count = 0

    do local_index = 1,n_local_blocks()
       b = local_block_catalog(local_index)
       destination = source_rank(b)
       if (destination /= rank) cycle
       d = loc_id(block_catalog(b)%root_domain+1) + 1
       if (d < 1 .or. d > size(grid)) then
          call fail("retained writeback block has invalid local Domain")
       end if

       local_patch = 0
       n_patch = 0
       call stage_local_subtree( &
            d,block_catalog(b)%root_patch,b,payload_family, &
            local_patch,n_patch)
       expected_patch_count = local_block_patch_count(b)
       if (local_patch /= expected_patch_count .or. &
            n_patch /= expected_patch_count) then
          call fail("retained writeback block reconstruction is incomplete")
       end if
       reconstructed_patch_count = reconstructed_patch_count + n_patch
    end do

    do r = 1,n_process
       pos_scalar = block_writeback_plan%scalar_recv_displ(r) + 1
       pos_vector = block_writeback_plan%vector_recv_displ(r) + 1

       do slot = block_writeback_plan%recv_displ(r)+1, &
            block_writeback_plan%recv_displ(r) + &
            block_writeback_plan%recv_count(r)
          b = block_writeback_plan%recv_block(slot)
          destination = source_rank(b)
          if (destination /= rank) then
             call fail("received reconstruction block has wrong Domain owner")
          end if
          if (block_catalog(b)%owner /= r-1) then
             call fail("received reconstruction block has wrong source owner")
          end if

          d = loc_id(block_catalog(b)%root_domain+1) + 1
          if (d < 1 .or. d > size(grid)) then
             call fail("received writeback block has invalid local Domain")
          end if

          scalar_start = pos_scalar
          vector_start = pos_vector
          scalar_limit = scalar_start + &
               block_writeback_plan%recv_scalar_nvalue(slot) - 1
          vector_limit = vector_start + &
               block_writeback_plan%recv_vector_nvalue(slot) - 1
          n_patch = 0
          call stage_buffer_subtree( &
               d,block_catalog(b)%root_patch, &
               block_writeback_plan%scalar_recv_buffer,pos_scalar, &
               scalar_limit,block_writeback_plan%vector_recv_buffer, &
               pos_vector,vector_limit,n_patch)

          if (n_patch /= &
               block_writeback_plan%recv_patch_count(slot) .or. &
               pos_scalar-scalar_start /= &
               block_writeback_plan%recv_scalar_nvalue(slot) .or. &
               pos_vector-vector_start /= &
               block_writeback_plan%recv_vector_nvalue(slot)) then
             call fail("received Domain reconstruction extent mismatch")
          end if
          reconstructed_patch_count = reconstructed_patch_count + n_patch
       end do

       if (pos_scalar /= block_writeback_plan%scalar_recv_displ(r) + &
            block_writeback_plan%scalar_recv_count(r) + 1 .or. &
            pos_vector /= block_writeback_plan%vector_recv_displ(r) + &
            block_writeback_plan%vector_recv_count(r) + 1) then
          call fail("Domain reconstruction receive extent mismatch")
       end if
    end do

    if (count(block_writeback_plan%domain_patch_covered) /= &
         reconstructed_patch_count) then
       call fail("block-derived Domain patch coverage is inconsistent")
    end if

    active_patch_count = 0
    do d = 1,size(grid)
       do p = 0,grid(d)%patch%length-1
          patch_slot = block_writeback_plan%domain_patch_displ(d) + p + 1
          if (grid(d)%patch%elts(p+1)%deleted) then
             if (block_writeback_plan%domain_patch_covered(patch_slot)) then
                call fail("deleted Domain patch was reconstructed")
             end if
          else
             active_patch_count = active_patch_count + 1
             if (.not. &
                  block_writeback_plan%domain_patch_covered(patch_slot)) then
                scalar_stage_start = (patch_slot-1)* &
                     block_writeback_plan%scalar_patch_nvalue + 1
                vector_stage_start = (patch_slot-1)* &
                     block_writeback_plan%vector_patch_nvalue + 1
                call pack_domain_patch_prognostic( &
                     d,p,payload_family, &
                     block_writeback_plan%scalar_domain_stage, &
                     scalar_stage_start, &
                     block_writeback_plan%vector_domain_stage, &
                     vector_stage_start)
                block_writeback_plan%domain_patch_covered(patch_slot) = &
                     .true.
                preserved_patch_count = preserved_patch_count + 1
             end if
          end if
       end do
    end do

    if (count(block_writeback_plan%domain_patch_covered) /= &
         active_patch_count) then
       call fail("Domain reconstruction patch coverage is not exact")
    end if
    if (active_patch_count /= reconstructed_patch_count + &
         preserved_patch_count) then
       call fail("complete Domain reconstruction count is inconsistent")
    end if

    block_writeback_plan%reconstructed_patch_count = &
         reconstructed_patch_count
    block_writeback_plan%preserved_patch_count = preserved_patch_count

  contains

    subroutine claim_domain_patch ( &
         d,p,scalar_stage_start,vector_stage_start)

      implicit none

      integer, intent(in) :: d
      integer, intent(in) :: p
      integer, intent(out) :: scalar_stage_start
      integer, intent(out) :: vector_stage_start

      integer :: claimed_slot

      if (p < 0 .or. p >= grid(d)%patch%length) then
         call fail("invalid patch in Domain reconstruction")
      end if
      if (grid(d)%patch%elts(p+1)%deleted) then
         call fail("deleted patch in Domain reconstruction")
      end if

      claimed_slot = &
           block_writeback_plan%domain_patch_displ(d) + p + 1
      if (block_writeback_plan%domain_patch_covered(claimed_slot)) then
         call fail("Domain patch was reconstructed more than once")
      end if
      block_writeback_plan%domain_patch_covered(claimed_slot) = .true.

      scalar_stage_start = (claimed_slot-1)* &
           block_writeback_plan%scalar_patch_nvalue + 1
      vector_stage_start = (claimed_slot-1)* &
           block_writeback_plan%vector_patch_nvalue + 1

    end subroutine claim_domain_patch


    recursive subroutine stage_local_subtree ( &
         d,p,b,payload_family,local_patch,n_patch)

      implicit none

      integer, intent(in) :: d
      integer, intent(in) :: p
      integer, intent(in) :: b
      integer, intent(in) :: payload_family
      integer, intent(inout) :: local_patch
      integer, intent(inout) :: n_patch

      integer :: c
      integer :: p_child
      integer :: scalar_stage_start
      integer :: vector_stage_start

      if (p < 0 .or. p >= grid(d)%patch%length) then
         call fail("invalid retained patch in Domain reconstruction")
      end if
      if (grid(d)%patch%elts(p+1)%deleted) return
      if (local_patch >= local_block_patch_count(b)) then
         call fail("retained block patch traversal exceeded extent")
      end if

      call claim_domain_patch( &
           d,p,scalar_stage_start,vector_stage_start)
      call get_local_block_scalar_patch_family_values( &
           b,local_patch,payload_family, &
           block_writeback_plan%scalar_domain_stage( &
           scalar_stage_start:scalar_stage_start+ &
           block_writeback_plan%scalar_patch_nvalue-1))
      call get_local_block_vector_patch_family_values( &
           b,local_patch,payload_family, &
           block_writeback_plan%vector_domain_stage( &
           vector_stage_start:vector_stage_start+ &
           block_writeback_plan%vector_patch_nvalue-1))
      local_patch = local_patch + 1
      n_patch = n_patch + 1

      do c = 1,N_CHDRN
         p_child = grid(d)%patch%elts(p+1)%children(c)
         if (p_child <= 0) cycle
         call stage_local_subtree( &
              d,p_child,b,payload_family,local_patch,n_patch)
      end do

    end subroutine stage_local_subtree


    recursive subroutine stage_buffer_subtree ( &
         d,p,scalar_payload,scalar_pos,scalar_limit, &
         vector_payload,vector_pos,vector_limit,n_patch)

      implicit none

      integer, intent(in) :: d
      integer, intent(in) :: p
      real(dp), intent(in) :: scalar_payload(:)
      integer, intent(inout) :: scalar_pos
      integer, intent(in) :: scalar_limit
      real(dp), intent(in) :: vector_payload(:)
      integer, intent(inout) :: vector_pos
      integer, intent(in) :: vector_limit
      integer, intent(inout) :: n_patch

      integer :: c
      integer :: p_child
      integer :: scalar_stage_start
      integer :: vector_stage_start

      if (p < 0 .or. p >= grid(d)%patch%length) then
         call fail("invalid received patch in Domain reconstruction")
      end if
      if (grid(d)%patch%elts(p+1)%deleted) return
      if (scalar_pos+block_writeback_plan%scalar_patch_nvalue-1 > &
           scalar_limit .or. &
           vector_pos+block_writeback_plan%vector_patch_nvalue-1 > &
           vector_limit) then
         call fail("received block payload ended during reconstruction")
      end if

      call claim_domain_patch( &
           d,p,scalar_stage_start,vector_stage_start)
      block_writeback_plan%scalar_domain_stage( &
           scalar_stage_start:scalar_stage_start+ &
           block_writeback_plan%scalar_patch_nvalue-1) = &
           scalar_payload(scalar_pos:scalar_pos+ &
           block_writeback_plan%scalar_patch_nvalue-1)
      block_writeback_plan%vector_domain_stage( &
           vector_stage_start:vector_stage_start+ &
           block_writeback_plan%vector_patch_nvalue-1) = &
           vector_payload(vector_pos:vector_pos+ &
           block_writeback_plan%vector_patch_nvalue-1)
      scalar_pos = scalar_pos + &
           block_writeback_plan%scalar_patch_nvalue
      vector_pos = vector_pos + &
           block_writeback_plan%vector_patch_nvalue
      n_patch = n_patch + 1

      do c = 1,N_CHDRN
         p_child = grid(d)%patch%elts(p+1)%children(c)
         if (p_child <= 0) cycle
         call stage_buffer_subtree( &
              d,p_child,scalar_payload,scalar_pos,scalar_limit, &
              vector_payload,vector_pos,vector_limit,n_patch)
      end do

    end subroutine stage_buffer_subtree

  end subroutine reconstruct_block_writeback_domain_stage


  logical function domain_patch_prognostic_extent_is_valid ( &
       d,p,payload_family) result(valid)
    ! Preflight every authoritative array section touched by one patch.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: payload_family

    integer :: field_level
    integer :: first_field_level
    integer :: level_slot
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: n_value
    integer :: scalar_id
    integer :: scalar_slot
    integer :: start
    integer :: v_scalar
    integer :: v_vector

    valid = .false.
    if (d < 1 .or. d > size(grid)) return
    if (p < 0 .or. p >= grid(d)%patch%length) return
    if (grid(d)%patch%elts(p+1)%deleted) return
    if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
         payload_family /= BLOCK_PAYLOAD_WAV_COEFF) return

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,first_field_level, &
         n_field_level,mult_scalar,mult_vector)
    if (v_scalar < 1 .or. v_vector < 1 .or. &
         n_scalar_variable < 1 .or. n_field_level < 1 .or. &
         mult_scalar < 1 .or. mult_vector < 1) return

    start = mult_scalar*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_scalar*PATCH_SIZE**2
    if (start < 0) return
    do scalar_slot = 1,n_scalar_variable
       scalar_id = v_scalar + scalar_slot - 1
       do level_slot = 1,n_field_level
          field_level = first_field_level + level_slot - 1
          select case (payload_family)
          case (BLOCK_PAYLOAD_SOL)
             if (start+n_value > &
                  size(sol(scalar_id,field_level)%data(d)%elts)) return
          case (BLOCK_PAYLOAD_WAV_COEFF)
             if (start+n_value > &
                  size(wav_coeff(scalar_id,field_level)%data(d)%elts)) &
                  return
          end select
       end do
    end do

    start = mult_vector*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_vector*PATCH_SIZE**2
    if (start < 0) return
    do level_slot = 1,n_field_level
       field_level = first_field_level + level_slot - 1
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (start+n_value > &
               size(sol(v_vector,field_level)%data(d)%elts)) return
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (start+n_value > &
               size(wav_coeff(v_vector,field_level)%data(d)%elts)) return
       end select
    end do

    valid = .true.

  end function domain_patch_prognostic_extent_is_valid


  logical function block_writeback_domain_stage_is_valid ( &
       payload_family) result(valid)
    ! Validate complete coverage and every destination extent before the
    ! first authoritative value is changed.

    implicit none

    integer, intent(in) :: payload_family

    integer :: active_patch_count
    integer :: d
    integer :: p
    integer :: patch_slot
    integer :: scalar_start
    integer :: vector_start

    logical :: plan_ready
    logical :: patch_extent_valid

    valid = .false.
    plan_ready = block_writeback_plan_is_ready()
    if (.not. plan_ready) return
    if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
         payload_family /= BLOCK_PAYLOAD_WAV_COEFF) return
    if (.not. allocated(block_writeback_plan%domain_patch_displ)) return
    if (.not. allocated(block_writeback_plan%domain_patch_covered)) return
    if (.not. allocated(block_writeback_plan%scalar_domain_stage)) return
    if (.not. allocated(block_writeback_plan%vector_domain_stage)) return
    if (block_writeback_plan%scalar_patch_nvalue <= 0 .or. &
         block_writeback_plan%vector_patch_nvalue <= 0) return

    active_patch_count = 0
    do d = 1,size(grid)
       do p = 0,grid(d)%patch%length-1
          patch_slot = &
               block_writeback_plan%domain_patch_displ(d) + p + 1
          if (patch_slot < 1 .or. patch_slot > &
               size(block_writeback_plan%domain_patch_covered)) return
          if (grid(d)%patch%elts(p+1)%deleted) then
             if (block_writeback_plan%domain_patch_covered(patch_slot)) &
                  return
             cycle
          end if

          active_patch_count = active_patch_count + 1
          if (.not. &
               block_writeback_plan%domain_patch_covered(patch_slot)) return

          scalar_start = (patch_slot-1)* &
               block_writeback_plan%scalar_patch_nvalue + 1
          vector_start = (patch_slot-1)* &
               block_writeback_plan%vector_patch_nvalue + 1
          if (scalar_start < 1 .or. scalar_start + &
               block_writeback_plan%scalar_patch_nvalue - 1 > &
               size(block_writeback_plan%scalar_domain_stage)) return
          if (vector_start < 1 .or. vector_start + &
               block_writeback_plan%vector_patch_nvalue - 1 > &
               size(block_writeback_plan%vector_domain_stage)) return

          patch_extent_valid = &
               domain_patch_prognostic_extent_is_valid( &
               d,p,payload_family)
          if (.not. patch_extent_valid) return
       end do
    end do

    if (active_patch_count /= &
         block_writeback_plan%reconstructed_patch_count + &
         block_writeback_plan%preserved_patch_count) return
    if (count(block_writeback_plan%domain_patch_covered) /= &
         active_patch_count) return

    valid = .true.

  end function block_writeback_domain_stage_is_valid


  subroutine write_domain_patch_prognostic ( &
       d,p,payload_family,scalar_payload,scalar_pos, &
       vector_payload,vector_pos)
    ! Copy one already validated staged patch into authoritative fields.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: payload_family
    real(dp), intent(in) :: scalar_payload(:)
    integer, intent(inout) :: scalar_pos
    real(dp), intent(in) :: vector_payload(:)
    integer, intent(inout) :: vector_pos

    integer :: field_level
    integer :: first_field_level
    integer :: level_slot
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: n_value
    integer :: scalar_id
    integer :: scalar_slot
    integer :: start
    integer :: v_scalar
    integer :: v_vector

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,first_field_level, &
         n_field_level,mult_scalar,mult_vector)

    start = mult_scalar*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_scalar*PATCH_SIZE**2
    do scalar_slot = 1,n_scalar_variable
       scalar_id = v_scalar + scalar_slot - 1
       do level_slot = 1,n_field_level
          field_level = first_field_level + level_slot - 1
          select case (payload_family)
          case (BLOCK_PAYLOAD_SOL)
             sol(scalar_id,field_level)%data(d)%elts( &
                  start+1:start+n_value) = &
                  scalar_payload(scalar_pos:scalar_pos+n_value-1)
          case (BLOCK_PAYLOAD_WAV_COEFF)
             wav_coeff(scalar_id,field_level)%data(d)%elts( &
                  start+1:start+n_value) = &
                  scalar_payload(scalar_pos:scalar_pos+n_value-1)
          end select
          scalar_pos = scalar_pos + n_value
       end do
    end do

    start = mult_vector*grid(d)%patch%elts(p+1)%elts_start
    n_value = mult_vector*PATCH_SIZE**2
    do level_slot = 1,n_field_level
       field_level = first_field_level + level_slot - 1
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          sol(v_vector,field_level)%data(d)%elts( &
               start+1:start+n_value) = &
               vector_payload(vector_pos:vector_pos+n_value-1)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          wav_coeff(v_vector,field_level)%data(d)%elts( &
               start+1:start+n_value) = &
               vector_payload(vector_pos:vector_pos+n_value-1)
       end select
       vector_pos = vector_pos + n_value
    end do

  end subroutine write_domain_patch_prognostic


  logical function try_commit_block_writeback_domain_stage ( &
       payload_family) result(committed)
    ! Commit only after a complete preflight. A rejected transaction leaves
    ! every authoritative Domain value untouched.

    implicit none

    integer, intent(in) :: payload_family

    integer :: d
    integer :: p
    integer :: patch_slot
    integer :: scalar_pos
    integer :: vector_pos

    logical :: stage_valid

    committed = .false.
    stage_valid = &
         block_writeback_domain_stage_is_valid(payload_family)
    if (.not. stage_valid) return

    do d = 1,size(grid)
       do p = 0,grid(d)%patch%length-1
          if (grid(d)%patch%elts(p+1)%deleted) cycle
          patch_slot = &
               block_writeback_plan%domain_patch_displ(d) + p + 1
          scalar_pos = (patch_slot-1)* &
               block_writeback_plan%scalar_patch_nvalue + 1
          vector_pos = (patch_slot-1)* &
               block_writeback_plan%vector_patch_nvalue + 1
          call write_domain_patch_prognostic( &
               d,p,payload_family, &
               block_writeback_plan%scalar_domain_stage,scalar_pos, &
               block_writeback_plan%vector_domain_stage,vector_pos)
       end do
    end do

    committed = .true.

  end function try_commit_block_writeback_domain_stage


  subroutine write_block_field_family_to_domains (payload_family)
    ! Production-facing transaction for one prognostic field family.
    ! Reconstruction and complete validation precede authoritative mutation.

    implicit none

    integer, intent(in) :: payload_family

    logical :: committed

    call reconstruct_block_writeback_domain_stage(payload_family)
    committed = &
         try_commit_block_writeback_domain_stage(payload_family)
    if (.not. committed) then
       call fail("complete block-to-Domain transaction was rejected")
    end if
    block_writeback_plan%production_writeback_count = &
         block_writeback_plan%production_writeback_count + 1_int64

  end subroutine write_block_field_family_to_domains


  subroutine assert_block_domain_field_family_match (payload_family)
    ! Reconstruct one family from final-owner blocks and compare every
    ! active authoritative Domain value exactly. Domain fields are read only.

    implicit none

    integer, intent(in) :: payload_family

    integer :: d
    integer :: p
    integer :: patch_slot
    integer :: scalar_pos
    integer :: scalar_start
    integer :: vector_pos
    integer :: vector_start

    real(dp) :: current_scalar( &
         block_writeback_plan%scalar_patch_nvalue)
    real(dp) :: current_vector( &
         block_writeback_plan%vector_patch_nvalue)

    call reconstruct_block_writeback_domain_stage(payload_family)

    do d = 1,size(grid)
       do p = 0,grid(d)%patch%length-1
          if (grid(d)%patch%elts(p+1)%deleted) cycle
          patch_slot = &
               block_writeback_plan%domain_patch_displ(d) + p + 1
          scalar_start = (patch_slot-1)* &
               block_writeback_plan%scalar_patch_nvalue + 1
          vector_start = (patch_slot-1)* &
               block_writeback_plan%vector_patch_nvalue + 1
          scalar_pos = 1
          vector_pos = 1
          call pack_domain_patch_prognostic( &
               d,p,payload_family,current_scalar,scalar_pos, &
               current_vector,vector_pos)
          if (any(abs(current_scalar- &
               block_writeback_plan%scalar_domain_stage( &
               scalar_start:scalar_start+ &
               block_writeback_plan%scalar_patch_nvalue-1)) > &
               0.0_dp)) then
             call fail("production scalar block/Domain mismatch")
          end if
          if (any(abs(current_vector- &
               block_writeback_plan%vector_domain_stage( &
               vector_start:vector_start+ &
               block_writeback_plan%vector_patch_nvalue-1)) > &
               0.0_dp)) then
             call fail("production vector block/Domain mismatch")
          end if
       end do
    end do

  end subroutine assert_block_domain_field_family_match


  subroutine check_block_writeback_domain_reconstruction (verbose)
    ! Validate complete sol and wav_coeff Domain staging and transactional
    ! writeback. Rejected commits must not mutate authoritative fields;
    ! accepted commits are verified exactly and restored to their references.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: active_patch_count
    integer :: d
    integer :: ierr
    integer :: p
    integer :: patch_slot
    integer :: scalar_start
    integer :: sol_preserved_patch_count
    integer :: sol_reconstructed_patch_count
    integer :: vector_start

    integer(int64) :: allocation_before
    integer(int64) :: global_patch_count
    integer(int64) :: global_preserved_patch_count
    integer(int64) :: global_reconstructed_patch_count
    integer(int64) :: global_scalar_count
    integer(int64) :: global_vector_count
    integer(int64) :: local_patch_count
    integer(int64) :: local_preserved_patch_count
    integer(int64) :: local_reconstructed_patch_count
    integer(int64) :: local_scalar_count
    integer(int64) :: local_vector_count

    logical :: print_summary

    real(dp), allocatable :: current_scalar(:)
    real(dp), allocatable :: current_vector(:)
    real(dp), allocatable :: reference_scalar(:)
    real(dp), allocatable :: reference_vector(:)

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. block_writeback_plan_is_ready()) then
       call fail("Domain reconstruction check before plan is ready")
    end if

    allocate(reference_scalar( &
         size(block_writeback_plan%scalar_domain_stage)))
    allocate(reference_vector( &
         size(block_writeback_plan%vector_domain_stage)))
    allocate(current_scalar( &
         block_writeback_plan%scalar_patch_nvalue))
    allocate(current_vector( &
         block_writeback_plan%vector_patch_nvalue))

    allocation_before = block_writeback_plan%stage_allocations

    call check_family(BLOCK_PAYLOAD_SOL)
    sol_reconstructed_patch_count = &
         block_writeback_plan%reconstructed_patch_count
    sol_preserved_patch_count = &
         block_writeback_plan%preserved_patch_count
    call check_family(BLOCK_PAYLOAD_WAV_COEFF)
    if (block_writeback_plan%reconstructed_patch_count /= &
         sol_reconstructed_patch_count .or. &
         block_writeback_plan%preserved_patch_count /= &
         sol_preserved_patch_count) then
       call fail("sol and wav_coeff Domain patch coverage differs")
    end if

    if (block_writeback_plan%stage_allocations /= allocation_before) then
       call fail("Domain reconstruction reallocated persistent staging")
    end if

    active_patch_count = 0
    do d = 1,size(grid)
       do p = 0,grid(d)%patch%length-1
          if (.not. grid(d)%patch%elts(p+1)%deleted) then
             active_patch_count = active_patch_count + 1
          end if
       end do
    end do
    local_patch_count = int(active_patch_count,int64)
    local_reconstructed_patch_count = int( &
         block_writeback_plan%reconstructed_patch_count,int64)
    local_preserved_patch_count = int( &
         block_writeback_plan%preserved_patch_count,int64)
    local_scalar_count = local_patch_count*int( &
         block_writeback_plan%scalar_patch_nvalue,int64)
    local_vector_count = local_patch_count*int( &
         block_writeback_plan%vector_patch_nvalue,int64)

    call MPI_Allreduce(local_patch_count,global_patch_count,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce reconstructed Domain patches")
    call MPI_Allreduce( &
         local_reconstructed_patch_count,global_reconstructed_patch_count, &
         1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block-derived Domain patches")
    call MPI_Allreduce( &
         local_preserved_patch_count,global_preserved_patch_count, &
         1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce preserved Domain patches")
    call MPI_Allreduce(local_scalar_count,global_scalar_count,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce reconstructed scalar values")
    call MPI_Allreduce(local_vector_count,global_vector_count,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce reconstructed vector values")

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Complete Domain reconstruction for rank ",rank,":"
       write(6,'(a,i0)') "  active Domain patches = ",local_patch_count
       write(6,'(a,i0)') "  block-derived patches = ", &
            local_reconstructed_patch_count
       write(6,'(a,i0)') "  preserved coarse patches = ", &
            local_preserved_patch_count
       write(6,'(a,i0)') "  scalar values staged = ",local_scalar_count
       write(6,'(a,i0)') "  vector values staged = ",local_vector_count
       write(6,'(a)') "  exact complete sol reconstruction passed"
       write(6,'(a)') "  exact complete wav_coeff reconstruction passed"
       write(6,'(a)') "  repeated persistent Domain staging passed"
       write(6,'(a)') "  incomplete transaction rejected without mutation"
       write(6,'(a)') "  exact transactional sol writeback passed"
       write(6,'(a)') "  exact transactional wav_coeff writeback passed"
       write(6,'(a)') "  repeated transactional Domain writeback passed"
       write(6,'(a,/)') "  authoritative Domain fields restored exactly"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global reconstructed Domain patches = ",global_patch_count
       write(6,'(a,i0)') &
            "Global block-derived Domain patches = ", &
            global_reconstructed_patch_count
       write(6,'(a,i0)') &
            "Global preserved coarse Domain patches = ", &
            global_preserved_patch_count
       write(6,'(a,i0)') &
            "Global reconstructed scalar values = ",global_scalar_count
       write(6,'(a,i0)') &
            "Global reconstructed vector values = ",global_vector_count
       write(6,'(a,/)') &
            "Transactional block-to-Domain writeback passed"
    end if

    deallocate(reference_vector)
    deallocate(reference_scalar)
    deallocate(current_vector)
    deallocate(current_scalar)

  contains

    subroutine build_reference (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      reference_scalar = 0.0_dp
      reference_vector = 0.0_dp
      do d = 1,size(grid)
         do p = 0,grid(d)%patch%length-1
            if (grid(d)%patch%elts(p+1)%deleted) cycle
            patch_slot = &
                 block_writeback_plan%domain_patch_displ(d) + p + 1
            scalar_start = (patch_slot-1)* &
                 block_writeback_plan%scalar_patch_nvalue + 1
            vector_start = (patch_slot-1)* &
                 block_writeback_plan%vector_patch_nvalue + 1
            call pack_domain_patch_prognostic( &
                 d,p,payload_family,reference_scalar,scalar_start, &
                 reference_vector,vector_start)
         end do
      end do

    end subroutine build_reference


    subroutine compare_stage_and_authority (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      integer :: current_scalar_pos
      integer :: current_vector_pos

      if (any(abs(block_writeback_plan%scalar_domain_stage- &
           reference_scalar) > 0.0_dp)) then
         if (payload_family == BLOCK_PAYLOAD_SOL) then
            call fail("complete scalar sol reconstruction mismatch")
         else
            call fail("complete scalar wav_coeff reconstruction mismatch")
         end if
      end if
      if (any(abs(block_writeback_plan%vector_domain_stage- &
           reference_vector) > 0.0_dp)) then
         if (payload_family == BLOCK_PAYLOAD_SOL) then
            call fail("complete vector sol reconstruction mismatch")
         else
            call fail("complete vector wav_coeff reconstruction mismatch")
         end if
      end if

      do d = 1,size(grid)
         do p = 0,grid(d)%patch%length-1
            if (grid(d)%patch%elts(p+1)%deleted) cycle
            patch_slot = &
                 block_writeback_plan%domain_patch_displ(d) + p + 1
            scalar_start = (patch_slot-1)* &
                 block_writeback_plan%scalar_patch_nvalue + 1
            vector_start = (patch_slot-1)* &
                 block_writeback_plan%vector_patch_nvalue + 1
            current_scalar_pos = 1
            current_vector_pos = 1
            call pack_domain_patch_prognostic( &
                 d,p,payload_family,current_scalar,current_scalar_pos, &
                 current_vector,current_vector_pos)
            if (any(abs(current_scalar-reference_scalar( &
                 scalar_start:scalar_start+ &
                 block_writeback_plan%scalar_patch_nvalue-1)) > &
                 0.0_dp)) then
               call fail("authoritative scalar Domain field changed")
            end if
            if (any(abs(current_vector-reference_vector( &
                 vector_start:vector_start+ &
                 block_writeback_plan%vector_patch_nvalue-1)) > &
                 0.0_dp)) then
               call fail("authoritative vector Domain field changed")
            end if
         end do
      end do

    end subroutine compare_stage_and_authority


    subroutine compare_authority_to_stage (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      integer :: current_scalar_pos
      integer :: current_vector_pos

      do d = 1,size(grid)
         do p = 0,grid(d)%patch%length-1
            if (grid(d)%patch%elts(p+1)%deleted) cycle
            patch_slot = &
                 block_writeback_plan%domain_patch_displ(d) + p + 1
            scalar_start = (patch_slot-1)* &
                 block_writeback_plan%scalar_patch_nvalue + 1
            vector_start = (patch_slot-1)* &
                 block_writeback_plan%vector_patch_nvalue + 1
            current_scalar_pos = 1
            current_vector_pos = 1
            call pack_domain_patch_prognostic( &
                 d,p,payload_family,current_scalar,current_scalar_pos, &
                 current_vector,current_vector_pos)
            if (any(abs(current_scalar- &
                 block_writeback_plan%scalar_domain_stage( &
                 scalar_start:scalar_start+ &
                 block_writeback_plan%scalar_patch_nvalue-1)) > &
                 0.0_dp)) then
               call fail("transactional scalar Domain writeback mismatch")
            end if
            if (any(abs(current_vector- &
                 block_writeback_plan%vector_domain_stage( &
                 vector_start:vector_start+ &
                 block_writeback_plan%vector_patch_nvalue-1)) > &
                 0.0_dp)) then
               call fail("transactional vector Domain writeback mismatch")
            end if
         end do
      end do

    end subroutine compare_authority_to_stage


    subroutine check_transaction (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      integer :: first_active_slot

      logical :: committed

      first_active_slot = 0
      do d = 1,size(grid)
         do p = 0,grid(d)%patch%length-1
            if (grid(d)%patch%elts(p+1)%deleted) cycle
            first_active_slot = &
                 block_writeback_plan%domain_patch_displ(d) + p + 1
            exit
         end do
         if (first_active_slot > 0) exit
      end do
      if (first_active_slot <= 0) then
         call fail("transactional writeback has no active Domain patch")
      end if

      block_writeback_plan%domain_patch_covered(first_active_slot) = &
           .false.
      committed = &
           try_commit_block_writeback_domain_stage(payload_family)
      block_writeback_plan%domain_patch_covered(first_active_slot) = &
           .true.
      if (committed) then
         call fail("incomplete Domain writeback transaction was accepted")
      end if
      call compare_stage_and_authority(payload_family)

      block_writeback_plan%scalar_domain_stage = 0.125_dp
      block_writeback_plan%vector_domain_stage = -0.125_dp
      committed = &
           try_commit_block_writeback_domain_stage(payload_family)
      if (.not. committed) then
         call fail("complete Domain writeback transaction was rejected")
      end if
      call compare_authority_to_stage(payload_family)

      block_writeback_plan%scalar_domain_stage = reference_scalar
      block_writeback_plan%vector_domain_stage = reference_vector
      committed = &
           try_commit_block_writeback_domain_stage(payload_family)
      if (.not. committed) then
         call fail("Domain writeback restoration was rejected")
      end if
      call compare_stage_and_authority(payload_family)

    end subroutine check_transaction


    subroutine check_family (payload_family)

      implicit none

      integer, intent(in) :: payload_family

      call build_reference(payload_family)
      call reconstruct_block_writeback_domain_stage(payload_family)
      call compare_stage_and_authority(payload_family)
      call reconstruct_block_writeback_domain_stage(payload_family)
      call compare_stage_and_authority(payload_family)
      call check_transaction(payload_family)

    end subroutine check_family

  end subroutine check_block_writeback_domain_reconstruction


  subroutine check_block_ghost_request_manifest (verbose)
    ! Exchange field-independent NGB_BLOCK ghost requests from each final
    ! owner to the current owner of the compact source block. This checks
    ! the routing layer only; Float_Field payload values are not moved yet.

    implicit none

    integer, parameter :: REQUEST_SIZE = 4

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: i
    integer :: ierr
    integer :: n_local_request
    integer :: n_remote_recv
    integer :: n_remote_send
    integer :: n_request
    integer :: pos
    integer :: r
    integer :: source

    integer(int64) :: count_local(3)
    integer(int64) :: count_global(3)
    integer(int64) :: recv_sum_local(REQUEST_SIZE)
    integer(int64) :: recv_sum_global(REQUEST_SIZE)
    integer(int64) :: send_sum_local(REQUEST_SIZE)
    integer(int64) :: send_sum_global(REQUEST_SIZE)

    integer, allocatable :: destination_block(:)
    integer, allocatable :: destination_ghost(:)
    integer, allocatable :: fill(:)
    integer, allocatable :: ghost_count_global(:)
    integer, allocatable :: ghost_count_local(:)
    integer, allocatable :: patch_count_global(:)
    integer, allocatable :: patch_count_local(:)
    integer, allocatable :: recv_count(:)
    integer, allocatable :: recv_data(:)
    integer, allocatable :: recv_displ(:)
    integer, allocatable :: recv_record_count(:)
    integer, allocatable :: send_count(:)
    integer, allocatable :: send_data(:)
    integer, allocatable :: send_displ(:)
    integer, allocatable :: send_record_count(:)
    integer, allocatable :: source_block(:)
    integer, allocatable :: source_local_patch(:)
    integer, allocatable :: source_owner(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    allocate(ghost_count_global(size(block_catalog)))
    allocate(ghost_count_local(size(block_catalog)))
    allocate(patch_count_global(size(block_catalog)))
    allocate(patch_count_local(size(block_catalog)))

    ghost_count_local = 0
    patch_count_local = 0

    do b = 1, size(block_catalog)
       if (block_catalog(b)%owner /= rank) cycle
       if (catalog_local_block(b) < 1) then
          call fail("owned block missing from request manifest store")
       end if
       ghost_count_local(b) = local_block_ghost_count(b)
       patch_count_local(b) = local_block_patch_count(b)
    end do

    call MPI_Allreduce( &
         ghost_count_local,ghost_count_global,size(block_catalog), &
         MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block ghost counts")

    call MPI_Allreduce( &
         patch_count_local,patch_count_global,size(block_catalog), &
         MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block patch counts")

    if (any(patch_count_global <= 0) .or. &
         any(ghost_count_global < 0)) then
       call fail("invalid global request-manifest block extents")
    end if

    call get_local_block_ghost_requests( &
         source_block,source_local_patch,source_owner, &
         destination_block,destination_ghost)

    n_request = size(source_block)

    if (size(source_local_patch) /= n_request .or. &
         size(source_owner) /= n_request .or. &
         size(destination_block) /= n_request .or. &
         size(destination_ghost) /= n_request) then
       call fail("inconsistent local ghost request arrays")
    end if

    allocate(send_record_count(n_process))
    allocate(recv_record_count(n_process))
    allocate(send_count(n_process))
    allocate(recv_count(n_process))
    allocate(send_displ(n_process))
    allocate(recv_displ(n_process))
    allocate(fill(n_process))

    send_record_count = 0
    n_local_request = 0

    do i = 1, n_request
       source = source_block(i)

       if (source < 1 .or. source > size(block_catalog)) then
          call fail("ghost request has invalid source block")
       end if
       if (source_owner(i) /= block_catalog(source)%owner) then
          call fail("ghost request has stale source owner")
       end if
       if (source_local_patch(i) < 0 .or. &
            source_local_patch(i) >= patch_count_global(source)) then
          call fail("ghost request has invalid compact source patch")
       end if
       if (destination_block(i) < 1 .or. &
            destination_block(i) > size(block_catalog)) then
          call fail("ghost request has invalid destination block")
       end if
       if (block_catalog(destination_block(i))%owner /= rank .or. &
            destination_ghost(i) < 1 .or. &
            destination_ghost(i) > &
            ghost_count_global(destination_block(i))) then
          call fail("ghost request has invalid destination record")
       end if

       if (source_owner(i) == rank) then
          if (catalog_local_block(source) < 1) then
             call fail("local ghost source block is unavailable")
          end if
          n_local_request = n_local_request + 1
       else
          send_record_count(source_owner(i)+1) = &
               send_record_count(source_owner(i)+1) + 1
       end if
    end do

    call MPI_Alltoall( &
         send_record_count,1,MPI_INTEGER, &
         recv_record_count,1,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoall ghost request counts")

    n_remote_send = sum(send_record_count)
    n_remote_recv = sum(recv_record_count)

    send_count = REQUEST_SIZE*send_record_count
    recv_count = REQUEST_SIZE*recv_record_count
    send_displ(1) = 0
    recv_displ(1) = 0

    do r = 2, n_process
       send_displ(r) = send_displ(r-1) + send_count(r-1)
       recv_displ(r) = recv_displ(r-1) + recv_count(r-1)
    end do

    allocate(send_data(REQUEST_SIZE*n_remote_send))
    allocate(recv_data(REQUEST_SIZE*n_remote_recv))

    fill = 0
    send_sum_local = 0_int64

    do i = 1, n_request
       if (source_owner(i) == rank) cycle
       r = source_owner(i) + 1
       pos = send_displ(r) + REQUEST_SIZE*fill(r)
       send_data(pos+1:pos+REQUEST_SIZE) = [ &
            source_block(i),source_local_patch(i), &
            destination_block(i),destination_ghost(i) ]
       send_sum_local = send_sum_local + &
            int(send_data(pos+1:pos+REQUEST_SIZE),int64)
       fill(r) = fill(r) + 1
    end do

    if (any(fill /= send_record_count)) then
       call fail("ghost request send packing count mismatch")
    end if

    call MPI_Alltoallv( &
         send_data,send_count,send_displ,MPI_INTEGER, &
         recv_data,recv_count,recv_displ,MPI_INTEGER,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv ghost requests")

    recv_sum_local = 0_int64

    do r = 1, n_process
       do i = 0, recv_record_count(r)-1
          pos = recv_displ(r) + REQUEST_SIZE*i
          source = recv_data(pos+1)

          if (source < 1 .or. source > size(block_catalog)) then
             call fail("received ghost request has invalid source block")
          end if
          if (block_catalog(source)%owner /= rank) then
             call fail("received ghost request source is not local")
          end if
          if (catalog_local_block(source) < 1) then
             call fail("received ghost request source block is absent")
          end if
          if (recv_data(pos+2) < 0 .or. &
               recv_data(pos+2) >= patch_count_global(source)) then
             call fail("received ghost request has invalid source patch")
          end if
          b = recv_data(pos+3)
          if (b < 1 .or. b > size(block_catalog)) then
             call fail("received ghost request has invalid destination")
          end if
          if (block_catalog(b)%owner /= r-1 .or. &
               recv_data(pos+4) < 1 .or. &
               recv_data(pos+4) > ghost_count_global(b)) then
             call fail("received ghost request has invalid ghost record")
          end if

          recv_sum_local = recv_sum_local + &
               int(recv_data(pos+1:pos+REQUEST_SIZE),int64)
       end do
    end do

    count_local = int([ &
         n_local_request,n_remote_send,n_remote_recv ],int64)

    call MPI_Allreduce( &
         count_local,count_global,3,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce ghost request totals")

    call MPI_Allreduce( &
         send_sum_local,send_sum_global,REQUEST_SIZE, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce sent ghost request checksum")

    call MPI_Allreduce( &
         recv_sum_local,recv_sum_global,REQUEST_SIZE, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce received ghost request checksum")

    if (count_global(2) /= count_global(3) .or. &
         any(send_sum_global /= recv_sum_global)) then
       call fail("remote ghost request exchange mismatch")
    end if

    if (count_global(1)+count_global(2) /= &
         int(sum(ghost_count_global),int64)) then
       call fail("global ghost request inventory is incomplete")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Float_Field ghost request manifest for rank ",rank,":"
       write(6,'(a,i0)') &
            "  local-source requests  = ",n_local_request
       write(6,'(a,i0)') &
            "  remote requests sent   = ",n_remote_send
       write(6,'(a,i0)') &
            "  remote requests received = ",n_remote_recv
       write(6,'(a,/)') &
            "  request source/destination checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global local-source ghost requests  = ",count_global(1)
       write(6,'(a,i0)') &
            "Global remote ghost requests        = ",count_global(2)
       write(6,'(a,i0)') &
            "Global Float_Field ghost requests   = ", &
            count_global(1)+count_global(2)
       write(6,'(a,/)') &
            "Float_Field ghost request manifest exchange passed"
    end if

    deallocate(destination_block)
    deallocate(destination_ghost)
    deallocate(fill)
    deallocate(ghost_count_global)
    deallocate(ghost_count_local)
    deallocate(patch_count_global)
    deallocate(patch_count_local)
    deallocate(recv_count)
    deallocate(recv_data)
    deallocate(recv_displ)
    deallocate(recv_record_count)
    deallocate(send_count)
    deallocate(send_data)
    deallocate(send_displ)
    deallocate(send_record_count)
    deallocate(source_block)
    deallocate(source_local_patch)
    deallocate(source_owner)

  end subroutine check_block_ghost_request_manifest


  subroutine check_block_scalar_ghost_payload_exchange (verbose)
    ! Exercise the scalar exchange with poisoning and exact installation
    ! checks enabled.

    implicit none

    logical, optional, intent(in) :: verbose

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call exchange_block_scalar_ghost_payloads( &
         BLOCK_PAYLOAD_SOL,print_summary,.true.)
    call exchange_block_scalar_ghost_payloads( &
         BLOCK_PAYLOAD_WAV_COEFF,print_summary,.true.)

  end subroutine check_block_scalar_ghost_payload_exchange


  subroutine exchange_block_scalar_ghost_payloads ( &
       payload_family,print_summary,verify_installation)
    ! Use the persistent field-independent request plan to return the
    ! complete scalar sol/wav_coeff bundle for every ghost patch. Quiet
    ! production calls reuse request metadata and communication buffers,
    ! pack directly into retained storage, and perform no allocation or
    ! diagnostic reductions.

    implicit none

    integer, parameter :: REQUEST_SIZE = 4
    logical, intent(in) :: print_summary
    logical, intent(in) :: verify_installation
    integer, intent(in) :: payload_family

    integer :: destination
    integer :: field_level
    integer :: fill_record
    integer :: i
    integer :: ierr
    integer :: level_count
    integer :: n_local_request
    integer :: n_remote_send
    integer :: n_request
    integer :: n_value
    integer :: pos
    integer :: r
    integer :: scalar_count
    integer :: scalar_mult
    integer :: scalar_variable
    integer :: source
    integer :: vector_mult
    integer :: vector_variable

    integer(int64) :: count_global(3)
    integer(int64) :: count_local(3)
    integer(int64) :: local_value_count

    real(dp), allocatable :: expected(:)
    character(len=9) :: payload_name

    payload_name = ""

    select case (payload_family)
    case (BLOCK_PAYLOAD_SOL)
       payload_name = "sol"
    case (BLOCK_PAYLOAD_WAV_COEFF)
       payload_name = "wav_coeff"
    case default
       call fail("invalid scalar ghost payload family")
    end select

    call get_block_field_layout( &
         scalar_variable,scalar_count,vector_variable,field_level, &
         level_count,scalar_mult,vector_mult)

    if (scalar_count < 1 .or. level_count < 1 .or. &
         scalar_mult < 1) then
       call fail("invalid scalar payload layout")
    end if

    n_value = &
         scalar_count*level_count*scalar_mult*PATCH_SIZE**2

    if (.not. ghost_exchange_plan%ready) then
       call fail("scalar ghost payload exchange without persistent plan")
    end if

    n_request = ghost_exchange_plan%n_request
    n_local_request = ghost_exchange_plan%n_local_request
    n_remote_send = ghost_exchange_plan%n_remote_send

    if (n_value /= ghost_exchange_plan%scalar_n_value) then
       call fail("persistent scalar ghost payload size changed")
    end if

    do i = 1, n_request
       destination = ghost_exchange_plan%destination_block(i)
       if (local_block_scalar_family_patch_nvalue(destination) /= &
            n_value) then
          call fail("destination scalar ghost payload layout mismatch")
       end if
    end do

    if (verify_installation) allocate(expected(n_value))

    do r = 1, n_process
       do i = 0, ghost_exchange_plan%recv_record_count(r)-1
          pos = REQUEST_SIZE*( &
               ghost_exchange_plan%recv_record_displ(r)+i)
          source = ghost_exchange_plan%recv_data(pos+1)

          if (source < 1 .or. source > size(block_catalog)) then
             call fail("received scalar ghost source is invalid")
          end if

          if (block_catalog(source)%owner /= rank) then
             call fail("received scalar ghost source owner mismatch")
          end if

          if (catalog_local_block(source) < 1) then
             call fail("received scalar ghost source is not local")
          end if

          call get_local_block_scalar_patch_family_values( &
               source,ghost_exchange_plan%recv_data(pos+2), &
               payload_family, &
               ghost_exchange_plan%scalar_send_buffer( &
               ghost_exchange_plan%scalar_send_displ(r) + &
               n_value*i + 1: &
               ghost_exchange_plan%scalar_send_displ(r) + &
               n_value*(i+1)))
       end do
    end do

    call MPI_Alltoallv( &
         ghost_exchange_plan%scalar_send_buffer,ghost_exchange_plan%scalar_send_count,ghost_exchange_plan%scalar_send_displ, &
         MPI_DOUBLE_PRECISION,ghost_exchange_plan%scalar_recv_buffer,ghost_exchange_plan%scalar_recv_count, &
         ghost_exchange_plan%scalar_recv_displ,MPI_DOUBLE_PRECISION,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv scalar ghost payloads")

    if (verify_installation) then
       do i = 1, n_request
          if (ghost_exchange_plan%source_owner(i) /= rank) cycle

          call get_local_block_scalar_patch_family_values( &
               ghost_exchange_plan%source_block(i), &
               ghost_exchange_plan%source_local_patch(i), &
               payload_family, &
               ghost_exchange_plan%scalar_patch_buffer)
          call get_local_block_scalar_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs( &
               ghost_exchange_plan%scalar_patch_buffer-expected)) > &
               0.0_dp) then
             call fail("local scalar ghost payload values do not match")
          end if
       end do

       do r = 1, n_process
          do i = 0, ghost_exchange_plan%send_record_count(r)-1
             fill_record = &
                  ghost_exchange_plan%send_record_displ(r)+i+1
             destination = ghost_exchange_plan%destination_block( &
                  ghost_exchange_plan%request_index(fill_record))

             call get_local_block_scalar_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)

             pos = ghost_exchange_plan%scalar_recv_displ(r) + n_value*i
             if (maxval(abs( &
                  ghost_exchange_plan%scalar_recv_buffer(pos+1:pos+n_value)-expected)) > 0.0_dp) then
                call fail("remote scalar ghost payload values do not match")
             end if
          end do
       end do

       call fill_local_block_scalar_ghost_family_values( &
            payload_family,BLOCK_GHOST_POISON)
    end if

    do i = 1, n_request
       if (ghost_exchange_plan%source_owner(i) /= rank) cycle

       if (verify_installation) then
          call get_local_block_scalar_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs(expected-BLOCK_GHOST_POISON)) > 0.0_dp) then
             call fail("scalar ghost storage was not invalidated")
          end if
       end if

       call get_local_block_scalar_patch_family_values( &
            ghost_exchange_plan%source_block(i), &
            ghost_exchange_plan%source_local_patch(i), &
            payload_family,ghost_exchange_plan%scalar_patch_buffer)
       call set_local_block_scalar_ghost_family_values( &
            ghost_exchange_plan%destination_block(i), &
            ghost_exchange_plan%destination_ghost(i), &
            payload_family,ghost_exchange_plan%scalar_patch_buffer)
       if (verify_installation) then
          call get_local_block_scalar_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs( &
               ghost_exchange_plan%scalar_patch_buffer-expected)) > &
               0.0_dp) then
             call fail("local scalar ghost payload installation failed")
          end if
       end if
    end do

    do r = 1, n_process
       do i = 0, ghost_exchange_plan%send_record_count(r)-1
          fill_record = &
               ghost_exchange_plan%send_record_displ(r)+i+1
          destination = ghost_exchange_plan%destination_block( &
               ghost_exchange_plan%request_index(fill_record))

          if (verify_installation) then
             call get_local_block_scalar_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)
             if (maxval(abs(expected-BLOCK_GHOST_POISON)) > 0.0_dp) then
                call fail("scalar ghost storage was not invalidated")
             end if
          end if

          pos = ghost_exchange_plan%scalar_recv_displ(r) + n_value*i
          call set_local_block_scalar_ghost_family_values( &
               destination, &
               ghost_exchange_plan%destination_ghost( &
               ghost_exchange_plan%request_index(fill_record)), &
               payload_family,ghost_exchange_plan%scalar_recv_buffer(pos+1:pos+n_value))
          if (verify_installation) then
             call get_local_block_scalar_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)
             if (maxval(abs( &
                  ghost_exchange_plan%scalar_recv_buffer(pos+1:pos+n_value)-expected)) > 0.0_dp) then
                call fail("remote scalar ghost payload installation failed")
             end if
          end if
       end do
    end do

    if (print_summary) then
       local_value_count = int(n_value,int64)*int(n_request,int64)
       count_local(1) = int(n_local_request,int64)
       count_local(2) = int(n_remote_send,int64)
       count_local(3) = local_value_count

       call MPI_Allreduce( &
            count_local,count_global,3,MPI_INTEGER8,MPI_SUM,comm,ierr)
       call check_mpi(ierr,"MPI_Allreduce scalar ghost payload totals")
    end if

    if (print_summary) then
       write(6,'(/,a,a,a,i0,a)') &
            "Scalar ",trim(payload_name), &
            " ghost payloads for rank ",rank,":"
       write(6,'(a,i0)') &
            "  values per ghost patch = ",n_value
       write(6,'(a,i0)') &
            "  local payload patches  = ",n_local_request
       write(6,'(a,i0)') &
            "  remote payload patches = ",n_remote_send
       write(6,'(a,/)') &
            "  selected scalar ghost family restored from payloads"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,a,a,i0)') &
            "Global scalar ",trim(payload_name), &
            " ghost patches = ", &
            count_global(1)+count_global(2)
       write(6,'(a,a,a,i0)') &
            "Global scalar ",trim(payload_name), &
            " values = ",count_global(3)
       write(6,'(a,a,a,/)') &
            "Scalar ",trim(payload_name), &
            " ghost payload installation passed"
    end if

    if (allocated(expected)) deallocate(expected)

  end subroutine exchange_block_scalar_ghost_payloads


  subroutine check_block_vector_ghost_payload_exchange (verbose)
    ! Exercise the vector exchange with poisoning and exact installation
    ! checks enabled.

    implicit none

    logical, optional, intent(in) :: verbose

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call exchange_block_vector_ghost_payloads( &
         BLOCK_PAYLOAD_SOL,print_summary,.true.)
    call exchange_block_vector_ghost_payloads( &
         BLOCK_PAYLOAD_WAV_COEFF,print_summary,.true.)

  end subroutine check_block_vector_ghost_payload_exchange


  subroutine exchange_block_vector_ghost_payloads ( &
       payload_family,print_summary,verify_installation)
    ! Use the persistent field-independent request plan to return the
    ! complete vector sol/wav_coeff bundle for every ghost patch. Quiet
    ! production calls reuse request metadata and communication buffers,
    ! pack directly into retained storage, and perform no allocation or
    ! diagnostic reductions.

    implicit none

    integer, parameter :: REQUEST_SIZE = 4
    logical, intent(in) :: print_summary
    logical, intent(in) :: verify_installation
    integer, intent(in) :: payload_family

    integer :: destination
    integer :: field_level
    integer :: fill_record
    integer :: i
    integer :: ierr
    integer :: level_count
    integer :: n_local_request
    integer :: n_remote_send
    integer :: n_request
    integer :: n_value
    integer :: pos
    integer :: r
    integer :: scalar_count
    integer :: scalar_mult
    integer :: scalar_variable
    integer :: source
    integer :: vector_mult
    integer :: vector_variable

    integer(int64) :: count_global(3)
    integer(int64) :: count_local(3)
    integer(int64) :: local_value_count

    real(dp), allocatable :: expected(:)
    character(len=9) :: payload_name

    payload_name = ""

    select case (payload_family)
    case (BLOCK_PAYLOAD_SOL)
       payload_name = "sol"
    case (BLOCK_PAYLOAD_WAV_COEFF)
       payload_name = "wav_coeff"
    case default
       call fail("invalid vector ghost payload family")
    end select

    call get_block_field_layout( &
         scalar_variable,scalar_count,vector_variable,field_level, &
         level_count,scalar_mult,vector_mult)

    if (level_count < 1 .or. vector_mult < 1) then
       call fail("invalid vector payload layout")
    end if

    n_value = level_count*vector_mult*PATCH_SIZE**2

    if (.not. ghost_exchange_plan%ready) then
       call fail("vector ghost payload exchange without persistent plan")
    end if

    n_request = ghost_exchange_plan%n_request
    n_local_request = ghost_exchange_plan%n_local_request
    n_remote_send = ghost_exchange_plan%n_remote_send

    if (n_value /= ghost_exchange_plan%vector_n_value) then
       call fail("persistent vector ghost payload size changed")
    end if

    do i = 1, n_request
       destination = ghost_exchange_plan%destination_block(i)
       if (local_block_vector_family_patch_nvalue(destination) /= &
            n_value) then
          call fail("destination vector ghost payload layout mismatch")
       end if
    end do

    if (verify_installation) allocate(expected(n_value))

    do r = 1, n_process
       do i = 0, ghost_exchange_plan%recv_record_count(r)-1
          pos = REQUEST_SIZE*( &
               ghost_exchange_plan%recv_record_displ(r)+i)
          source = ghost_exchange_plan%recv_data(pos+1)

          if (source < 1 .or. source > size(block_catalog)) then
             call fail("received vector ghost source is invalid")
          end if

          if (block_catalog(source)%owner /= rank) then
             call fail("received vector ghost source owner mismatch")
          end if

          if (catalog_local_block(source) < 1) then
             call fail("received vector ghost source is not local")
          end if

          call get_local_block_vector_patch_family_values( &
               source,ghost_exchange_plan%recv_data(pos+2), &
               payload_family, &
               ghost_exchange_plan%vector_send_buffer( &
               ghost_exchange_plan%vector_send_displ(r) + &
               n_value*i + 1: &
               ghost_exchange_plan%vector_send_displ(r) + &
               n_value*(i+1)))
       end do
    end do

    call MPI_Alltoallv( &
         ghost_exchange_plan%vector_send_buffer,ghost_exchange_plan%vector_send_count,ghost_exchange_plan%vector_send_displ, &
         MPI_DOUBLE_PRECISION,ghost_exchange_plan%vector_recv_buffer,ghost_exchange_plan%vector_recv_count, &
         ghost_exchange_plan%vector_recv_displ,MPI_DOUBLE_PRECISION,comm,ierr)
    call check_mpi(ierr,"MPI_Alltoallv vector ghost payloads")

    if (verify_installation) then
       do i = 1, n_request
          if (ghost_exchange_plan%source_owner(i) /= rank) cycle

          call get_local_block_vector_patch_family_values( &
               ghost_exchange_plan%source_block(i), &
               ghost_exchange_plan%source_local_patch(i), &
               payload_family, &
               ghost_exchange_plan%vector_patch_buffer)
          call get_local_block_vector_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs( &
               ghost_exchange_plan%vector_patch_buffer-expected)) > &
               0.0_dp) then
             call fail("local vector ghost payload values do not match")
          end if
       end do

       do r = 1, n_process
          do i = 0, ghost_exchange_plan%send_record_count(r)-1
             fill_record = &
                  ghost_exchange_plan%send_record_displ(r)+i+1
             destination = ghost_exchange_plan%destination_block( &
                  ghost_exchange_plan%request_index(fill_record))

             call get_local_block_vector_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)

             pos = ghost_exchange_plan%vector_recv_displ(r) + n_value*i
             if (maxval(abs( &
                  ghost_exchange_plan%vector_recv_buffer(pos+1:pos+n_value)-expected)) > 0.0_dp) then
                call fail("remote vector ghost payload values do not match")
             end if
          end do
       end do

       call fill_local_block_vector_ghost_family_values( &
            payload_family,BLOCK_GHOST_POISON)
    end if

    do i = 1, n_request
       if (ghost_exchange_plan%source_owner(i) /= rank) cycle

       if (verify_installation) then
          call get_local_block_vector_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs(expected-BLOCK_GHOST_POISON)) > 0.0_dp) then
             call fail("vector ghost storage was not invalidated")
          end if
       end if

       call get_local_block_vector_patch_family_values( &
            ghost_exchange_plan%source_block(i), &
            ghost_exchange_plan%source_local_patch(i), &
            payload_family,ghost_exchange_plan%vector_patch_buffer)
       call set_local_block_vector_ghost_family_values( &
            ghost_exchange_plan%destination_block(i), &
            ghost_exchange_plan%destination_ghost(i), &
            payload_family,ghost_exchange_plan%vector_patch_buffer)
       if (verify_installation) then
          call get_local_block_vector_ghost_family_values( &
               ghost_exchange_plan%destination_block(i), &
               ghost_exchange_plan%destination_ghost(i), &
               payload_family,expected)
          if (maxval(abs( &
               ghost_exchange_plan%vector_patch_buffer-expected)) > &
               0.0_dp) then
             call fail("local vector ghost payload installation failed")
          end if
       end if
    end do

    do r = 1, n_process
       do i = 0, ghost_exchange_plan%send_record_count(r)-1
          fill_record = &
               ghost_exchange_plan%send_record_displ(r)+i+1
          destination = ghost_exchange_plan%destination_block( &
               ghost_exchange_plan%request_index(fill_record))

          if (verify_installation) then
             call get_local_block_vector_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)
             if (maxval(abs(expected-BLOCK_GHOST_POISON)) > 0.0_dp) then
                call fail("vector ghost storage was not invalidated")
             end if
          end if

          pos = ghost_exchange_plan%vector_recv_displ(r) + n_value*i
          call set_local_block_vector_ghost_family_values( &
               destination, &
               ghost_exchange_plan%destination_ghost( &
               ghost_exchange_plan%request_index(fill_record)), &
               payload_family,ghost_exchange_plan%vector_recv_buffer(pos+1:pos+n_value))
          if (verify_installation) then
             call get_local_block_vector_ghost_family_values( &
                  destination, &
                  ghost_exchange_plan%destination_ghost( &
                  ghost_exchange_plan%request_index(fill_record)), &
                  payload_family,expected)
             if (maxval(abs( &
                  ghost_exchange_plan%vector_recv_buffer(pos+1:pos+n_value)-expected)) > 0.0_dp) then
                call fail("remote vector ghost payload installation failed")
             end if
          end if
       end do
    end do

    if (print_summary) then
       local_value_count = int(n_value,int64)*int(n_request,int64)
       count_local(1) = int(n_local_request,int64)
       count_local(2) = int(n_remote_send,int64)
       count_local(3) = local_value_count

       call MPI_Allreduce( &
            count_local,count_global,3,MPI_INTEGER8,MPI_SUM,comm,ierr)
       call check_mpi(ierr,"MPI_Allreduce vector ghost payload totals")
    end if

    if (print_summary) then
       write(6,'(/,a,a,a,i0,a)') &
            "Vector ",trim(payload_name), &
            " ghost payloads for rank ",rank,":"
       write(6,'(a,i0)') &
            "  values per ghost patch = ",n_value
       write(6,'(a,i0)') &
            "  local payload patches  = ",n_local_request
       write(6,'(a,i0)') &
            "  remote payload patches = ",n_remote_send
       write(6,'(a,/)') &
            "  selected vector ghost family restored from payloads"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,a,a,i0)') &
            "Global vector ",trim(payload_name), &
            " ghost patches = ", &
            count_global(1)+count_global(2)
       write(6,'(a,a,a,i0)') &
            "Global vector ",trim(payload_name), &
            " values = ",count_global(3)
       write(6,'(a,a,a,/)') &
            "Vector ",trim(payload_name), &
            " ghost payload installation passed"
    end if

    if (allocated(expected)) deallocate(expected)

  end subroutine exchange_block_vector_ghost_payloads


  subroutine refresh_block_sol_ghosts
    ! Refresh only scalar and vector sol ghost payloads.

    implicit none

    if (.not. local_block_store_ready()) then
       call fail("block sol ghost refresh before block installation")
    end if

    call exchange_block_scalar_ghost_payloads( &
         BLOCK_PAYLOAD_SOL,.false.,.false.)
    call exchange_block_vector_ghost_payloads( &
         BLOCK_PAYLOAD_SOL,.false.,.false.)

  end subroutine refresh_block_sol_ghosts


  subroutine apply_refreshed_block_tendency_kernel (kernel,context)
    ! Enforce the production ordering required by every stencil-dependent
    ! block tendency: refresh sol ghosts before executing the local kernel.

    implicit none

    procedure(Local_Block_Tendency_Kernel) :: kernel
    class(*), intent(inout) :: context

    if (.not. local_block_store_ready()) then
       call fail("refreshed tendency kernel before block installation")
    end if
    if (local_block_tendency_trial_is_active()) then
       call fail("refreshed tendency kernel during active trial")
    end if

    call refresh_block_sol_ghosts
    call apply_local_block_tendency_kernel(kernel,context)

    if (.not. local_block_tendency_state_ready()) then
       call fail("refreshed tendency kernel output is not ready")
    end if

  end subroutine apply_refreshed_block_tendency_kernel


  subroutine begin_block_two_stage_tendency_step ( &
       kernel,context,scale,weight,result)
    ! Execute a guarded two-stage block update and retain its one-level
    ! checkpoint. The caller must subsequently finalize or restore it.

    implicit none

    procedure(Local_Block_Tendency_Kernel) :: kernel
    class(*), intent(inout) :: context
    real(dp), intent(in) :: scale
    real(dp), intent(in) :: weight(2)
    type(Block_Two_Stage_Step_Result), intent(out) :: result

    integer(int64) :: accumulator_changed_block_count(2)

    result = Block_Two_Stage_Step_Result()

    if (.not. local_block_store_ready()) then
       call fail("two-stage block step before block installation")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage block step found pending transaction")
    end if

    if (compressible) call ensure_local_block_hydrostatic_state
    call reset_local_block_tendency_accumulator

    call apply_refreshed_block_tendency_kernel(kernel,context)
    call accumulate_local_block_tendency(weight(1))

    call begin_local_block_tendency_trial(scale)
    call commit_local_block_tendency_trial

    if (compressible) call ensure_local_block_hydrostatic_state
    call apply_refreshed_block_tendency_kernel(kernel,context)
    call accumulate_local_block_tendency(weight(2))

    call local_block_tendency_accumulator_statistics( &
         result%scalar_count,result%vector_count, &
         accumulator_changed_block_count(1), &
         accumulator_changed_block_count(2),result%stage_count, &
         result%scalar_moment,result%vector_moment)

    call restore_local_block_tendency_commit
    if (local_block_tendency_commit_checkpoint_is_ready() .or. &
         local_block_tendency_state_ready()) then
       call fail("two-stage intermediate state was not restored")
    end if

    call begin_local_block_accumulated_tendency_trial(scale)
    call commit_local_block_tendency_trial
    call local_block_tendency_commit_checkpoint_statistics( &
         result%scalar_changed_block_count, &
         result%vector_changed_block_count, &
         result%scalar_max_update,result%vector_max_update)

    if (abs(scale) > 0.0_dp) then
       if (result%scalar_changed_block_count /= &
            accumulator_changed_block_count(1) .or. &
            result%vector_changed_block_count /= &
            accumulator_changed_block_count(2)) then
          call fail("two-stage committed update coverage mismatch")
       end if
    end if
    if (local_block_tendency_trial_is_active() .or. &
         .not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage block step did not retain its checkpoint")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("two-stage block step retained stale tendencies")
    end if

  end subroutine begin_block_two_stage_tendency_step


  subroutine complete_block_two_stage_tendency_step (accept)
    ! Resolve the checkpoint retained by begin_block_two_stage_tendency_step.
    ! Accepted sol fields are transactionally synchronized to their Domain
    ! owners before the block checkpoint is finalized. Rejected fields are
    ! restored without modifying the authoritative Domain representation.

    implicit none

    logical, intent(in) :: accept

    if (local_block_tendency_trial_is_active()) then
       call fail("two-stage completion found an active trial")
    end if
    if (.not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage completion found no pending checkpoint")
    end if

    if (accept) then
       call write_block_field_family_to_domains(BLOCK_PAYLOAD_SOL)
       call finalize_local_block_tendency_commit
    else
       call restore_local_block_tendency_commit
    end if

    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage completion left pending transaction state")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("two-stage completion retained stale tendencies")
    end if

  end subroutine complete_block_two_stage_tendency_step


  subroutine refresh_block_wav_coeff_ghosts
    ! Refresh only scalar and vector wav_coeff ghost payloads.

    implicit none

    if (.not. local_block_store_ready()) then
       call fail("block wav_coeff ghost refresh before block installation")
    end if

    call exchange_block_scalar_ghost_payloads( &
         BLOCK_PAYLOAD_WAV_COEFF,.false.,.false.)
    call exchange_block_vector_ghost_payloads( &
         BLOCK_PAYLOAD_WAV_COEFF,.false.,.false.)

  end subroutine refresh_block_wav_coeff_ghosts


  subroutine refresh_block_sol_wav_coeff_ghosts
    ! Refresh both payload families using the two selective production
    ! entry points while retaining the combined convenience API.

    implicit none

    call refresh_block_sol_ghosts
    call refresh_block_wav_coeff_ghosts

  end subroutine refresh_block_sol_wav_coeff_ghosts


  subroutine check_production_block_ghost_refresh (verbose)
    ! Poison and restore each family independently through its production
    ! entry point. After each selective refresh, verify the complete scalar
    ! and vector state through the compact stencil consumers.

    implicit none

    logical, optional, intent(in) :: verbose

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call fill_local_block_scalar_ghost_family_values( &
         BLOCK_PAYLOAD_SOL,BLOCK_GHOST_POISON)
    call fill_local_block_vector_ghost_family_values( &
         BLOCK_PAYLOAD_SOL,BLOCK_GHOST_POISON)

    call refresh_block_sol_ghosts
    call check_refreshed_block_stencil_consumers(.false.)

    call fill_local_block_scalar_ghost_family_values( &
         BLOCK_PAYLOAD_WAV_COEFF,BLOCK_GHOST_POISON)
    call fill_local_block_vector_ghost_family_values( &
         BLOCK_PAYLOAD_WAV_COEFF,BLOCK_GHOST_POISON)

    call refresh_block_wav_coeff_ghosts
    call check_refreshed_block_stencil_consumers(.false.)
    call check_block_ghost_exchange_plan(.false.)

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Production block ghost refresh for rank ",rank,":"
       write(6,'(a)') &
            "  selective sol ghost refresh passed"
       write(6,'(a)') &
            "  selective wav_coeff ghost refresh passed"
       write(6,'(a,/)') &
            "  persistent request plan and buffers retained"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Reusable-buffer production block ghost refreshes passed"
    end if

  end subroutine check_production_block_ghost_refresh


  subroutine check_block_field_family_accessors (verbose)
    ! Verify that independent sol and wav_coeff views reconstruct the
    ! existing combined payload exactly for every local patch, boundary
    ! record and ghost. Boundary access is read-only in this stage.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: boundary_index
    integer :: catalog_index
    integer :: ghost_index
    integer :: ierr
    integer :: n_scalar_boundary_combined
    integer :: n_scalar_boundary_family
    integer :: n_scalar_combined
    integer :: n_scalar_family
    integer :: n_vector_boundary_combined
    integer :: n_vector_boundary_family
    integer :: n_vector_combined
    integer :: n_vector_family
    integer :: patch_index

    integer(int64) :: boundary_count_global(2)
    integer(int64) :: boundary_count_local(2)
    integer(int64) :: count_global(4)
    integer(int64) :: count_local(4)

    real(dp), allocatable :: scalar_boundary_combined(:)
    real(dp), allocatable :: scalar_boundary_family(:)
    real(dp), allocatable :: scalar_combined(:)
    real(dp), allocatable :: scalar_family(:)
    real(dp), allocatable :: vector_boundary_combined(:)
    real(dp), allocatable :: vector_boundary_family(:)
    real(dp), allocatable :: vector_combined(:)
    real(dp), allocatable :: vector_family(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("field-family accessor check before block installation")
    end if

    boundary_count_local = 0_int64
    count_local = 0_int64

    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)

       n_scalar_family = &
            local_block_scalar_family_patch_nvalue(catalog_index)
       n_vector_family = &
            local_block_vector_family_patch_nvalue(catalog_index)
       n_scalar_combined = &
            local_block_scalar_patch_nvalue(catalog_index)
       n_vector_combined = &
            local_block_vector_patch_nvalue(catalog_index)

       if (n_scalar_combined /= 2*n_scalar_family) then
          call fail("combined and scalar-family payload extents differ")
       end if

       if (n_vector_combined /= 2*n_vector_family) then
          call fail("combined and vector-family payload extents differ")
       end if

       allocate(scalar_combined(2*n_scalar_family))
       allocate(scalar_family(n_scalar_family))
       allocate(vector_combined(2*n_vector_family))
       allocate(vector_family(n_vector_family))

       do patch_index = 0, local_block_patch_count(catalog_index)-1
          call get_local_block_scalar_patch_values( &
               catalog_index,patch_index,scalar_combined)
          call get_local_block_scalar_patch_family_values( &
               catalog_index,patch_index,BLOCK_PAYLOAD_SOL, &
               scalar_family)
          if (maxval(abs( &
               scalar_family-scalar_combined(1:n_scalar_family))) > &
               0.0_dp) then
             call fail("selective scalar sol patch accessor mismatch")
          end if
          call get_local_block_scalar_patch_family_values( &
               catalog_index,patch_index,BLOCK_PAYLOAD_WAV_COEFF, &
               scalar_family)
          if (maxval(abs( &
               scalar_family- &
               scalar_combined(n_scalar_family+1: &
               2*n_scalar_family))) > 0.0_dp) then
             call fail("selective scalar wav_coeff patch accessor mismatch")
          end if
          count_local(1) = count_local(1) + &
               int(2*n_scalar_family,int64)

          call get_local_block_vector_patch_values( &
               catalog_index,patch_index,vector_combined)
          call get_local_block_vector_patch_family_values( &
               catalog_index,patch_index,BLOCK_PAYLOAD_SOL, &
               vector_family)
          if (maxval(abs( &
               vector_family-vector_combined(1:n_vector_family))) > &
               0.0_dp) then
             call fail("selective vector sol patch accessor mismatch")
          end if
          call get_local_block_vector_patch_family_values( &
               catalog_index,patch_index,BLOCK_PAYLOAD_WAV_COEFF, &
               vector_family)
          if (maxval(abs( &
               vector_family- &
               vector_combined(n_vector_family+1: &
               2*n_vector_family))) > 0.0_dp) then
             call fail("selective vector wav_coeff patch accessor mismatch")
          end if
          count_local(3) = count_local(3) + &
               int(2*n_vector_family,int64)
       end do

       do boundary_index = 1, &
            local_block_boundary_count(catalog_index)
          n_scalar_boundary_family = &
               local_block_scalar_family_boundary_nvalue( &
               catalog_index,boundary_index)
          n_scalar_boundary_combined = &
               local_block_scalar_boundary_nvalue( &
               catalog_index,boundary_index)
          n_vector_boundary_family = &
               local_block_vector_family_boundary_nvalue( &
               catalog_index,boundary_index)
          n_vector_boundary_combined = &
               local_block_vector_boundary_nvalue( &
               catalog_index,boundary_index)

          if (n_scalar_boundary_combined /= &
               2*n_scalar_boundary_family) then
             call fail( &
                  "combined and scalar boundary-family extents differ")
          end if
          if (n_vector_boundary_combined /= &
               2*n_vector_boundary_family) then
             call fail( &
                  "combined and vector boundary-family extents differ")
          end if

          allocate(scalar_boundary_combined( &
               n_scalar_boundary_combined))
          allocate(scalar_boundary_family( &
               n_scalar_boundary_family))
          allocate(vector_boundary_combined( &
               n_vector_boundary_combined))
          allocate(vector_boundary_family( &
               n_vector_boundary_family))

          call get_local_block_scalar_boundary_values( &
               catalog_index,boundary_index, &
               scalar_boundary_combined)
          call get_local_block_scalar_boundary_family_values( &
               catalog_index,boundary_index,BLOCK_PAYLOAD_SOL, &
               scalar_boundary_family)
          if (maxval(abs(scalar_boundary_family - &
               scalar_boundary_combined( &
               1:n_scalar_boundary_family))) > 0.0_dp) then
             call fail( &
                  "selective scalar sol boundary accessor mismatch")
          end if
          call get_local_block_scalar_boundary_family_values( &
               catalog_index,boundary_index, &
               BLOCK_PAYLOAD_WAV_COEFF,scalar_boundary_family)
          if (maxval(abs(scalar_boundary_family - &
               scalar_boundary_combined( &
               n_scalar_boundary_family+1: &
               n_scalar_boundary_combined))) > 0.0_dp) then
             call fail( &
                  "selective scalar wav_coeff boundary accessor mismatch")
          end if
          boundary_count_local(1) = boundary_count_local(1) + &
               int(n_scalar_boundary_combined,int64)

          call get_local_block_vector_boundary_values( &
               catalog_index,boundary_index, &
               vector_boundary_combined)
          call get_local_block_vector_boundary_family_values( &
               catalog_index,boundary_index,BLOCK_PAYLOAD_SOL, &
               vector_boundary_family)
          if (maxval(abs(vector_boundary_family - &
               vector_boundary_combined( &
               1:n_vector_boundary_family))) > 0.0_dp) then
             call fail( &
                  "selective vector sol boundary accessor mismatch")
          end if
          call get_local_block_vector_boundary_family_values( &
               catalog_index,boundary_index, &
               BLOCK_PAYLOAD_WAV_COEFF,vector_boundary_family)
          if (maxval(abs(vector_boundary_family - &
               vector_boundary_combined( &
               n_vector_boundary_family+1: &
               n_vector_boundary_combined))) > 0.0_dp) then
             call fail( &
                  "selective vector wav_coeff boundary accessor mismatch")
          end if
          boundary_count_local(2) = boundary_count_local(2) + &
               int(n_vector_boundary_combined,int64)

          deallocate(scalar_boundary_combined)
          deallocate(scalar_boundary_family)
          deallocate(vector_boundary_combined)
          deallocate(vector_boundary_family)
       end do

       do ghost_index = 1, local_block_ghost_count(catalog_index)
          call get_local_block_scalar_ghost_values( &
               catalog_index,ghost_index,scalar_combined)
          call get_local_block_scalar_ghost_family_values( &
               catalog_index,ghost_index,BLOCK_PAYLOAD_SOL, &
               scalar_family)
          if (maxval(abs( &
               scalar_family-scalar_combined(1:n_scalar_family))) > &
               0.0_dp) then
             call fail("selective scalar sol ghost accessor mismatch")
          end if
          call get_local_block_scalar_ghost_family_values( &
               catalog_index,ghost_index,BLOCK_PAYLOAD_WAV_COEFF, &
               scalar_family)
          if (maxval(abs( &
               scalar_family- &
               scalar_combined(n_scalar_family+1: &
               2*n_scalar_family))) > 0.0_dp) then
             call fail("selective scalar wav_coeff ghost accessor mismatch")
          end if
          count_local(2) = count_local(2) + &
               int(2*n_scalar_family,int64)

          call get_local_block_vector_ghost_values( &
               catalog_index,ghost_index,vector_combined)
          call get_local_block_vector_ghost_family_values( &
               catalog_index,ghost_index,BLOCK_PAYLOAD_SOL, &
               vector_family)
          if (maxval(abs( &
               vector_family-vector_combined(1:n_vector_family))) > &
               0.0_dp) then
             call fail("selective vector sol ghost accessor mismatch")
          end if
          call get_local_block_vector_ghost_family_values( &
               catalog_index,ghost_index,BLOCK_PAYLOAD_WAV_COEFF, &
               vector_family)
          if (maxval(abs( &
               vector_family- &
               vector_combined(n_vector_family+1: &
               2*n_vector_family))) > 0.0_dp) then
             call fail("selective vector wav_coeff ghost accessor mismatch")
          end if
          count_local(4) = count_local(4) + &
               int(2*n_vector_family,int64)
       end do

       deallocate(scalar_combined)
       deallocate(scalar_family)
       deallocate(vector_combined)
       deallocate(vector_family)
    end do

    call MPI_Allreduce( &
         count_local,count_global,4,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce field-family accessor values")

    call MPI_Allreduce( &
         boundary_count_local,boundary_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi( &
         ierr,"MPI_Allreduce boundary field-family accessor values")

    if (any(count_global <= 0_int64)) then
       call fail("incomplete field-family accessor inventory")
    end if
    if (any(boundary_count_global <= 0_int64)) then
       call fail("incomplete boundary field-family accessor inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Selective field-family accessors for rank ",rank,":"
       write(6,'(a,i0)') &
            "  scalar patch values verified = ",count_local(1)
       write(6,'(a,i0)') &
            "  scalar boundary values verified = ", &
            boundary_count_local(1)
       write(6,'(a,i0)') &
            "  scalar ghost values verified = ",count_local(2)
       write(6,'(a,i0)') &
            "  vector patch values verified = ",count_local(3)
       write(6,'(a,i0)') &
            "  vector boundary values verified = ", &
            boundary_count_local(2)
       write(6,'(a,i0)') &
            "  vector ghost values verified = ",count_local(4)
       write(6,'(a,/)') &
            "  sol/wav_coeff family reconstruction checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,4(i0,1x))') &
            "Global scalar-patch/scalar-ghost/vector-patch/" // &
            "vector-ghost values = ",count_global
       write(6,'(a,2(i0,1x))') &
            "Global scalar/vector boundary values = ", &
            boundary_count_global
       write(6,'(a,/)') &
            "Selective sol/wav_coeff field-family accessors passed"
    end if

  end subroutine check_block_field_family_accessors


  subroutine check_block_patch_writable_storage (verbose)
    ! Snapshot every compact patch, exercise family-selective and combined
    ! whole-store writes, and prove exact isolation and restoration.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: catalog_index
    integer :: ierr
    integer :: local_patch
    integer :: n_patch_record
    integer :: record_index

    integer(int64) :: count_global(2)
    integer(int64) :: count_local(2)

    logical :: print_summary

    type(Block_Patch_Snapshot), allocatable :: snapshot(:)

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("writable patch check before block installation")
    end if

    n_patch_record = 0
    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)
       n_patch_record = n_patch_record + &
            local_block_patch_count(catalog_index)
    end do

    allocate(snapshot(n_patch_record))
    count_local = 0_int64
    record_index = 0

    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)
       do local_patch = 0,local_block_patch_count(catalog_index)-1
          record_index = record_index + 1
          snapshot(record_index)%catalog_index = catalog_index
          snapshot(record_index)%local_patch = local_patch

          allocate(snapshot(record_index)%scalar( &
               local_block_scalar_patch_nvalue(catalog_index)))
          allocate(snapshot(record_index)%vector( &
               local_block_vector_patch_nvalue(catalog_index)))

          call get_local_block_scalar_patch_values( &
               catalog_index,local_patch,snapshot(record_index)%scalar)
          call get_local_block_vector_patch_values( &
               catalog_index,local_patch,snapshot(record_index)%vector)

          count_local(1) = count_local(1) + &
               int(size(snapshot(record_index)%scalar),int64)
          count_local(2) = count_local(2) + &
               int(size(snapshot(record_index)%vector),int64)
       end do
    end do

    if (record_index /= n_patch_record) then
       call fail("writable patch snapshot count mismatch")
    end if

    call check_scalar_patch_family_bulk_write( &
         snapshot,BLOCK_PAYLOAD_SOL)
    call check_scalar_patch_family_bulk_write( &
         snapshot,BLOCK_PAYLOAD_WAV_COEFF)
    call check_vector_patch_family_bulk_write( &
         snapshot,BLOCK_PAYLOAD_SOL)
    call check_vector_patch_family_bulk_write( &
         snapshot,BLOCK_PAYLOAD_WAV_COEFF)
    call check_scalar_patch_combined_bulk_write(snapshot)
    call check_vector_patch_combined_bulk_write(snapshot)

    call check_refreshed_block_stencil_consumers(.false.)

    call MPI_Allreduce( &
         count_local,count_global,2,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce writable patch values")

    if (any(count_global <= 0_int64)) then
       call fail("incomplete writable patch inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Writable compact patch storage for rank ",rank,":"
       write(6,'(a,i0)') &
            "  scalar patch values bulk-tested = ",count_local(1)
       write(6,'(a,i0)') &
            "  vector patch values bulk-tested = ",count_local(2)
       write(6,'(a,/)') &
            "  family and combined poison, isolation and restore passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global scalar/vector patch values bulk-tested = ", &
            count_global
       write(6,'(a,/)') &
            "Writable sol/wav_coeff compact patch storage passed"
    end if

    deallocate(snapshot)

  end subroutine check_block_patch_writable_storage


  subroutine check_scalar_patch_family_bulk_write ( &
       snapshot,payload_family)
    ! Fill and restore one scalar family across the complete patch store.

    implicit none

    type(Block_Patch_Snapshot), intent(in) :: snapshot(:)
    integer, intent(in) :: payload_family

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)

    call fill_local_block_scalar_patch_family_values( &
         payload_family,BLOCK_PATCH_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%scalar),2) /= 0 .or. &
            size(snapshot(i)%scalar) <= 0) then
          call fail("invalid scalar writable patch snapshot")
       end if
       n_family = size(snapshot(i)%scalar)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))

       call get_local_block_scalar_patch_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            payload_family,family)
       if (maxval(abs(family-BLOCK_PATCH_POISON)) > 0.0_dp) then
          call fail("scalar patch family bulk poison failed")
       end if

       call get_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (maxval(abs(combined(n_family+1:2*n_family) - &
               snapshot(i)%scalar(n_family+1:2*n_family))) > &
               0.0_dp) then
             call fail("scalar patch wav_coeff changed during family fill")
          end if
          family = snapshot(i)%scalar(1:n_family)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (maxval(abs(combined(1:n_family) - &
               snapshot(i)%scalar(1:n_family))) > 0.0_dp) then
             call fail("scalar patch sol changed during family fill")
          end if
          family = snapshot(i)%scalar(n_family+1:2*n_family)
       case default
          call fail("invalid scalar writable patch family")
       end select

       call set_local_block_scalar_patch_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            payload_family,family)
       call get_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-snapshot(i)%scalar)) > 0.0_dp) then
          call fail("scalar patch family bulk restore failed")
       end if

       deallocate(combined)
       deallocate(family)
    end do

  end subroutine check_scalar_patch_family_bulk_write


  subroutine check_vector_patch_family_bulk_write ( &
       snapshot,payload_family)
    ! Fill and restore one vector family across the complete patch store.

    implicit none

    type(Block_Patch_Snapshot), intent(in) :: snapshot(:)
    integer, intent(in) :: payload_family

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)

    call fill_local_block_vector_patch_family_values( &
         payload_family,BLOCK_PATCH_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%vector),2) /= 0 .or. &
            size(snapshot(i)%vector) <= 0) then
          call fail("invalid vector writable patch snapshot")
       end if
       n_family = size(snapshot(i)%vector)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))

       call get_local_block_vector_patch_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            payload_family,family)
       if (maxval(abs(family-BLOCK_PATCH_POISON)) > 0.0_dp) then
          call fail("vector patch family bulk poison failed")
       end if

       call get_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (maxval(abs(combined(n_family+1:2*n_family) - &
               snapshot(i)%vector(n_family+1:2*n_family))) > &
               0.0_dp) then
             call fail("vector patch wav_coeff changed during family fill")
          end if
          family = snapshot(i)%vector(1:n_family)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (maxval(abs(combined(1:n_family) - &
               snapshot(i)%vector(1:n_family))) > 0.0_dp) then
             call fail("vector patch sol changed during family fill")
          end if
          family = snapshot(i)%vector(n_family+1:2*n_family)
       case default
          call fail("invalid vector writable patch family")
       end select

       call set_local_block_vector_patch_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            payload_family,family)
       call get_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-snapshot(i)%vector)) > 0.0_dp) then
          call fail("vector patch family bulk restore failed")
       end if

       deallocate(combined)
       deallocate(family)
    end do

  end subroutine check_vector_patch_family_bulk_write


  subroutine check_scalar_patch_combined_bulk_write (snapshot)
    ! Fill both scalar families, verify vector isolation and restore.

    implicit none

    type(Block_Patch_Snapshot), intent(in) :: snapshot(:)

    integer :: i

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: vector_observed(:)

    call fill_local_block_scalar_patch_values(BLOCK_PATCH_POISON)

    do i = 1, size(snapshot)
       if (size(snapshot(i)%scalar) <= 0 .or. &
            size(snapshot(i)%vector) <= 0) then
          call fail("invalid scalar combined patch snapshot")
       end if
       allocate(combined(size(snapshot(i)%scalar)))
       allocate(vector_observed(size(snapshot(i)%vector)))

       call get_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-BLOCK_PATCH_POISON)) > 0.0_dp) then
          call fail("scalar combined patch poison failed")
       end if

       call get_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            vector_observed)
       if (maxval(abs(vector_observed-snapshot(i)%vector)) > &
            0.0_dp) then
          call fail("vector patch changed during scalar combined fill")
       end if

       call set_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            snapshot(i)%scalar)
       call get_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-snapshot(i)%scalar)) > 0.0_dp) then
          call fail("scalar combined patch restore failed")
       end if

       deallocate(combined)
       deallocate(vector_observed)
    end do

  end subroutine check_scalar_patch_combined_bulk_write


  subroutine check_vector_patch_combined_bulk_write (snapshot)
    ! Fill both vector families, verify scalar isolation and restore.

    implicit none

    type(Block_Patch_Snapshot), intent(in) :: snapshot(:)

    integer :: i

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: scalar_observed(:)

    call fill_local_block_vector_patch_values(BLOCK_PATCH_POISON)

    do i = 1, size(snapshot)
       if (size(snapshot(i)%vector) <= 0 .or. &
            size(snapshot(i)%scalar) <= 0) then
          call fail("invalid vector combined patch snapshot")
       end if
       allocate(combined(size(snapshot(i)%vector)))
       allocate(scalar_observed(size(snapshot(i)%scalar)))

       call get_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-BLOCK_PATCH_POISON)) > 0.0_dp) then
          call fail("vector combined patch poison failed")
       end if

       call get_local_block_scalar_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            scalar_observed)
       if (maxval(abs(scalar_observed-snapshot(i)%scalar)) > &
            0.0_dp) then
          call fail("scalar patch changed during vector combined fill")
       end if

       call set_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            snapshot(i)%vector)
       call get_local_block_vector_patch_values( &
            snapshot(i)%catalog_index,snapshot(i)%local_patch, &
            combined)
       if (maxval(abs(combined-snapshot(i)%vector)) > 0.0_dp) then
          call fail("vector combined patch restore failed")
       end if

       deallocate(combined)
       deallocate(scalar_observed)
    end do

  end subroutine check_vector_patch_combined_bulk_write


  subroutine check_block_boundary_family_mutators (verbose)
    ! Poison, read back and restore each scalar and vector boundary family
    ! independently. The opposite family must remain unchanged, and the
    ! complete combined payload must be recovered exactly after each test.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: boundary_index
    integer :: catalog_index
    integer :: ierr
    integer :: n_scalar_family
    integer :: n_vector_family
    integer :: payload_family

    integer(int64) :: count_global(2)
    integer(int64) :: count_local(2)

    real(dp), allocatable :: scalar_combined_observed(:)
    real(dp), allocatable :: scalar_combined_original(:)
    real(dp), allocatable :: scalar_family_observed(:)
    real(dp), allocatable :: scalar_family_original(:)
    real(dp), allocatable :: vector_combined_observed(:)
    real(dp), allocatable :: vector_combined_original(:)
    real(dp), allocatable :: vector_family_observed(:)
    real(dp), allocatable :: vector_family_original(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("boundary-family mutation before block installation")
    end if

    count_local = 0_int64

    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)

       do boundary_index = 1, &
            local_block_boundary_count(catalog_index)
          n_scalar_family = &
               local_block_scalar_family_boundary_nvalue( &
               catalog_index,boundary_index)
          n_vector_family = &
               local_block_vector_family_boundary_nvalue( &
               catalog_index,boundary_index)

          allocate(scalar_combined_original(2*n_scalar_family))
          allocate(scalar_combined_observed(2*n_scalar_family))
          allocate(scalar_family_original(n_scalar_family))
          allocate(scalar_family_observed(n_scalar_family))
          allocate(vector_combined_original(2*n_vector_family))
          allocate(vector_combined_observed(2*n_vector_family))
          allocate(vector_family_original(n_vector_family))
          allocate(vector_family_observed(n_vector_family))

          call get_local_block_scalar_boundary_values( &
               catalog_index,boundary_index, &
               scalar_combined_original)
          call get_local_block_vector_boundary_values( &
               catalog_index,boundary_index, &
               vector_combined_original)

          do payload_family = &
               BLOCK_PAYLOAD_SOL,BLOCK_PAYLOAD_WAV_COEFF
             select case (payload_family)
             case (BLOCK_PAYLOAD_SOL)
                scalar_family_original = &
                     scalar_combined_original(1:n_scalar_family)
                vector_family_original = &
                     vector_combined_original(1:n_vector_family)
             case (BLOCK_PAYLOAD_WAV_COEFF)
                scalar_family_original = &
                     scalar_combined_original(n_scalar_family+1: &
                     2*n_scalar_family)
                vector_family_original = &
                     vector_combined_original(n_vector_family+1: &
                     2*n_vector_family)
             end select

             scalar_family_observed = BLOCK_BOUNDARY_POISON
             call set_local_block_scalar_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  scalar_family_observed)
             call get_local_block_scalar_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  scalar_family_observed)
             if (maxval(abs( &
                  scalar_family_observed-BLOCK_BOUNDARY_POISON)) > &
                  0.0_dp) then
                call fail("scalar boundary family poison failed")
             end if
             call get_local_block_scalar_boundary_values( &
                  catalog_index,boundary_index, &
                  scalar_combined_observed)
             select case (payload_family)
             case (BLOCK_PAYLOAD_SOL)
                if (maxval(abs( &
                     scalar_combined_observed(n_scalar_family+1: &
                     2*n_scalar_family) - &
                     scalar_combined_original(n_scalar_family+1: &
                     2*n_scalar_family))) > 0.0_dp) then
                   call fail("scalar wav_coeff boundary family changed")
                end if
             case (BLOCK_PAYLOAD_WAV_COEFF)
                if (maxval(abs( &
                     scalar_combined_observed(1:n_scalar_family) - &
                     scalar_combined_original(1:n_scalar_family))) > &
                     0.0_dp) then
                   call fail("scalar sol boundary family changed")
                end if
             end select
             call set_local_block_scalar_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  scalar_family_original)
             call get_local_block_scalar_boundary_values( &
                  catalog_index,boundary_index, &
                  scalar_combined_observed)
             if (maxval(abs(scalar_combined_observed - &
                  scalar_combined_original)) > 0.0_dp) then
                call fail("scalar boundary family restore failed")
             end if

             vector_family_observed = BLOCK_BOUNDARY_POISON
             call set_local_block_vector_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  vector_family_observed)
             call get_local_block_vector_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  vector_family_observed)
             if (maxval(abs( &
                  vector_family_observed-BLOCK_BOUNDARY_POISON)) > &
                  0.0_dp) then
                call fail("vector boundary family poison failed")
             end if
             call get_local_block_vector_boundary_values( &
                  catalog_index,boundary_index, &
                  vector_combined_observed)
             select case (payload_family)
             case (BLOCK_PAYLOAD_SOL)
                if (maxval(abs( &
                     vector_combined_observed(n_vector_family+1: &
                     2*n_vector_family) - &
                     vector_combined_original(n_vector_family+1: &
                     2*n_vector_family))) > 0.0_dp) then
                   call fail("vector wav_coeff boundary family changed")
                end if
             case (BLOCK_PAYLOAD_WAV_COEFF)
                if (maxval(abs( &
                     vector_combined_observed(1:n_vector_family) - &
                     vector_combined_original(1:n_vector_family))) > &
                     0.0_dp) then
                   call fail("vector sol boundary family changed")
                end if
             end select
             call set_local_block_vector_boundary_family_values( &
                  catalog_index,boundary_index,payload_family, &
                  vector_family_original)
             call get_local_block_vector_boundary_values( &
                  catalog_index,boundary_index, &
                  vector_combined_observed)
             if (maxval(abs(vector_combined_observed - &
                  vector_combined_original)) > 0.0_dp) then
                call fail("vector boundary family restore failed")
             end if
          end do

          count_local(1) = count_local(1) + &
               int(2*n_scalar_family,int64)
          count_local(2) = count_local(2) + &
               int(2*n_vector_family,int64)

          deallocate(scalar_combined_observed)
          deallocate(scalar_combined_original)
          deallocate(scalar_family_observed)
          deallocate(scalar_family_original)
          deallocate(vector_combined_observed)
          deallocate(vector_combined_original)
          deallocate(vector_family_observed)
          deallocate(vector_family_original)
       end do
    end do

    call check_refreshed_block_stencil_consumers(.false.)

    call MPI_Allreduce( &
         count_local,count_global,2,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce boundary-family mutation values")

    if (any(count_global <= 0_int64)) then
       call fail("incomplete boundary-family mutation inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Writable boundary-family accessors for rank ",rank,":"
       write(6,'(a,i0)') &
            "  scalar boundary values round-tripped = ", &
            count_local(1)
       write(6,'(a,i0)') &
            "  vector boundary values round-tripped = ", &
            count_local(2)
       write(6,'(a,/)') &
            "  independent sol/wav_coeff poison and restore passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global scalar/vector boundary values round-tripped = ", &
            count_global
       write(6,'(a,/)') &
            "Writable sol/wav_coeff boundary-family accessors passed"
    end if

  end subroutine check_block_boundary_family_mutators


  subroutine check_block_boundary_family_bulk_fill (verbose)
    ! Snapshot every compact boundary, fill one complete payload family,
    ! and prove exact family isolation and restoration over the store.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: boundary_index
    integer :: catalog_index
    integer :: ierr
    integer :: n_boundary_record
    integer :: record_index

    integer(int64) :: count_global(2)
    integer(int64) :: count_local(2)

    logical :: print_summary

    type(Block_Boundary_Snapshot), allocatable :: snapshot(:)

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("boundary-family bulk fill before block installation")
    end if

    n_boundary_record = 0
    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)
       n_boundary_record = n_boundary_record + &
            local_block_boundary_count(catalog_index)
    end do

    allocate(snapshot(n_boundary_record))
    count_local = 0_int64
    record_index = 0

    do b = 1, n_local_blocks()
       catalog_index = local_block_catalog(b)
       do boundary_index = 1, &
            local_block_boundary_count(catalog_index)
          record_index = record_index + 1
          snapshot(record_index)%catalog_index = catalog_index
          snapshot(record_index)%boundary_index = boundary_index

          allocate(snapshot(record_index)%scalar( &
               local_block_scalar_boundary_nvalue( &
               catalog_index,boundary_index)))
          allocate(snapshot(record_index)%vector( &
               local_block_vector_boundary_nvalue( &
               catalog_index,boundary_index)))

          call get_local_block_scalar_boundary_values( &
               catalog_index,boundary_index, &
               snapshot(record_index)%scalar)
          call get_local_block_vector_boundary_values( &
               catalog_index,boundary_index, &
               snapshot(record_index)%vector)

          count_local(1) = count_local(1) + &
               int(size(snapshot(record_index)%scalar),int64)
          count_local(2) = count_local(2) + &
               int(size(snapshot(record_index)%vector),int64)
       end do
    end do

    if (record_index /= n_boundary_record) then
       call fail("boundary-family bulk snapshot count mismatch")
    end if

    call check_scalar_boundary_family_bulk_fill( &
         snapshot,BLOCK_PAYLOAD_SOL)
    call check_scalar_boundary_family_bulk_fill( &
         snapshot,BLOCK_PAYLOAD_WAV_COEFF)
    call check_vector_boundary_family_bulk_fill( &
         snapshot,BLOCK_PAYLOAD_SOL)
    call check_vector_boundary_family_bulk_fill( &
         snapshot,BLOCK_PAYLOAD_WAV_COEFF)
    call check_scalar_boundary_combined_bulk_fill(snapshot)
    call check_vector_boundary_combined_bulk_fill(snapshot)

    call check_refreshed_block_stencil_consumers(.false.)

    call MPI_Allreduce( &
         count_local,count_global,2,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce boundary-family bulk values")

    if (any(count_global <= 0_int64)) then
       call fail("incomplete boundary-family bulk-fill inventory")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Bulk boundary-family fill for rank ",rank,":"
       write(6,'(a,i0)') &
            "  scalar boundary values preserved = ",count_local(1)
       write(6,'(a,i0)') &
            "  vector boundary values preserved = ",count_local(2)
       write(6,'(a)') &
            "  whole-store sol/wav_coeff poison and restore passed"
       write(6,'(a,/)') &
            "  combined-family poison, isolation and restore passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global scalar/vector boundary values bulk-tested = ", &
            count_global
       write(6,'(a)') &
            "Bulk sol/wav_coeff boundary-family fills passed"
       write(6,'(a,/)') &
            "Bulk combined scalar/vector boundary fills passed"
    end if

    deallocate(snapshot)

  end subroutine check_block_boundary_family_bulk_fill


  subroutine check_scalar_boundary_family_bulk_fill ( &
       snapshot,payload_family)
    ! Check one scalar whole-store fill and restore its saved family.

    implicit none

    type(Block_Boundary_Snapshot), intent(in) :: snapshot(:)
    integer, intent(in) :: payload_family

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)

    call fill_local_block_scalar_boundary_family_values( &
         payload_family,BLOCK_BOUNDARY_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%scalar),2) /= 0 .or. &
            size(snapshot(i)%scalar) <= 0) then
          call fail("invalid scalar boundary bulk snapshot")
       end if
       n_family = size(snapshot(i)%scalar)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))

       call get_local_block_scalar_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            payload_family,family)
       if (maxval(abs(family-BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("scalar boundary bulk poison failed")
       end if

       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (maxval(abs(combined(n_family+1:2*n_family) - &
               snapshot(i)%scalar(n_family+1:2*n_family))) > &
               0.0_dp) then
             call fail("scalar wav_coeff changed during bulk fill")
          end if
          family = snapshot(i)%scalar(1:n_family)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (maxval(abs(combined(1:n_family) - &
               snapshot(i)%scalar(1:n_family))) > 0.0_dp) then
             call fail("scalar sol changed during bulk fill")
          end if
          family = snapshot(i)%scalar(n_family+1:2*n_family)
       case default
          call fail("invalid scalar boundary bulk family")
       end select

       call set_local_block_scalar_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            payload_family,family)
       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-snapshot(i)%scalar)) > 0.0_dp) then
          call fail("scalar boundary bulk restore failed")
       end if

       deallocate(combined)
       deallocate(family)
    end do

  end subroutine check_scalar_boundary_family_bulk_fill


  subroutine check_vector_boundary_family_bulk_fill ( &
       snapshot,payload_family)
    ! Check one vector whole-store fill and restore its saved family.

    implicit none

    type(Block_Boundary_Snapshot), intent(in) :: snapshot(:)
    integer, intent(in) :: payload_family

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)

    call fill_local_block_vector_boundary_family_values( &
         payload_family,BLOCK_BOUNDARY_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%vector),2) /= 0 .or. &
            size(snapshot(i)%vector) <= 0) then
          call fail("invalid vector boundary bulk snapshot")
       end if
       n_family = size(snapshot(i)%vector)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))

       call get_local_block_vector_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            payload_family,family)
       if (maxval(abs(family-BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("vector boundary bulk poison failed")
       end if

       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       select case (payload_family)
       case (BLOCK_PAYLOAD_SOL)
          if (maxval(abs(combined(n_family+1:2*n_family) - &
               snapshot(i)%vector(n_family+1:2*n_family))) > &
               0.0_dp) then
             call fail("vector wav_coeff changed during bulk fill")
          end if
          family = snapshot(i)%vector(1:n_family)
       case (BLOCK_PAYLOAD_WAV_COEFF)
          if (maxval(abs(combined(1:n_family) - &
               snapshot(i)%vector(1:n_family))) > 0.0_dp) then
             call fail("vector sol changed during bulk fill")
          end if
          family = snapshot(i)%vector(n_family+1:2*n_family)
       case default
          call fail("invalid vector boundary bulk family")
       end select

       call set_local_block_vector_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            payload_family,family)
       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-snapshot(i)%vector)) > 0.0_dp) then
          call fail("vector boundary bulk restore failed")
       end if

       deallocate(combined)
       deallocate(family)
    end do

  end subroutine check_vector_boundary_family_bulk_fill


  subroutine check_scalar_boundary_combined_bulk_fill (snapshot)
    ! Check a combined scalar fill, vector isolation and family restore.

    implicit none

    type(Block_Boundary_Snapshot), intent(in) :: snapshot(:)

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)
    real(dp), allocatable :: vector_observed(:)

    call fill_local_block_scalar_boundary_values(BLOCK_BOUNDARY_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%scalar),2) /= 0 .or. &
            size(snapshot(i)%scalar) <= 0 .or. &
            size(snapshot(i)%vector) <= 0) then
          call fail("invalid scalar combined boundary snapshot")
       end if
       n_family = size(snapshot(i)%scalar)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))
       allocate(vector_observed(size(snapshot(i)%vector)))

       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("scalar combined boundary poison failed")
       end if

       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            vector_observed)
       if (maxval(abs(vector_observed-snapshot(i)%vector)) > &
            0.0_dp) then
          call fail("vector boundary changed during scalar fill")
       end if

       family = snapshot(i)%scalar(1:n_family)
       call set_local_block_scalar_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            BLOCK_PAYLOAD_SOL,family)
       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined(1:n_family) - &
            snapshot(i)%scalar(1:n_family))) > 0.0_dp) then
          call fail("scalar sol combined boundary restore failed")
       end if
       if (maxval(abs(combined(n_family+1:2*n_family) - &
            BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("scalar wav_coeff restored before installation")
       end if

       family = snapshot(i)%scalar(n_family+1:2*n_family)
       call set_local_block_scalar_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            BLOCK_PAYLOAD_WAV_COEFF,family)
       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-snapshot(i)%scalar)) > 0.0_dp) then
          call fail("scalar combined boundary restore failed")
       end if

       deallocate(combined)
       deallocate(family)
       deallocate(vector_observed)
    end do

  end subroutine check_scalar_boundary_combined_bulk_fill


  subroutine check_vector_boundary_combined_bulk_fill (snapshot)
    ! Check a combined vector fill, scalar isolation and family restore.

    implicit none

    type(Block_Boundary_Snapshot), intent(in) :: snapshot(:)

    integer :: i
    integer :: n_family

    real(dp), allocatable :: combined(:)
    real(dp), allocatable :: family(:)
    real(dp), allocatable :: scalar_observed(:)

    call fill_local_block_vector_boundary_values(BLOCK_BOUNDARY_POISON)

    do i = 1, size(snapshot)
       if (mod(size(snapshot(i)%vector),2) /= 0 .or. &
            size(snapshot(i)%vector) <= 0 .or. &
            size(snapshot(i)%scalar) <= 0) then
          call fail("invalid vector combined boundary snapshot")
       end if
       n_family = size(snapshot(i)%vector)/2
       allocate(combined(2*n_family))
       allocate(family(n_family))
       allocate(scalar_observed(size(snapshot(i)%scalar)))

       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("vector combined boundary poison failed")
       end if

       call get_local_block_scalar_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            scalar_observed)
       if (maxval(abs(scalar_observed-snapshot(i)%scalar)) > &
            0.0_dp) then
          call fail("scalar boundary changed during vector fill")
       end if

       family = snapshot(i)%vector(1:n_family)
       call set_local_block_vector_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            BLOCK_PAYLOAD_SOL,family)
       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined(1:n_family) - &
            snapshot(i)%vector(1:n_family))) > 0.0_dp) then
          call fail("vector sol combined boundary restore failed")
       end if
       if (maxval(abs(combined(n_family+1:2*n_family) - &
            BLOCK_BOUNDARY_POISON)) > 0.0_dp) then
          call fail("vector wav_coeff restored before installation")
       end if

       family = snapshot(i)%vector(n_family+1:2*n_family)
       call set_local_block_vector_boundary_family_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            BLOCK_PAYLOAD_WAV_COEFF,family)
       call get_local_block_vector_boundary_values( &
            snapshot(i)%catalog_index,snapshot(i)%boundary_index, &
            combined)
       if (maxval(abs(combined-snapshot(i)%vector)) > 0.0_dp) then
          call fail("vector combined boundary restore failed")
       end if

       deallocate(combined)
       deallocate(family)
       deallocate(scalar_observed)
    end do

  end subroutine check_vector_boundary_combined_bulk_fill


  subroutine check_local_blocks (verbose)
    ! Validate the final-owner local block store without referring to
    ! source, receive or migration-manifest staging allocation.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: expected_local
    integer :: expected_field_level
    integer :: expected_n_field_level
    integer :: expected_n_scalar_variable
    integer :: expected_scalar_mult
    integer :: expected_scalar_variable
    integer :: expected_vector_mult
    integer :: expected_vector_variable
    integer :: expected_tke_level
    integer :: expected_n_tke_level
    integer :: field_level
    integer :: global_count
    integer :: global_weight
    integer :: i
    integer :: id
    integer :: ierr
    integer :: level
    integer :: local_count
    integer :: local_weight
    integer :: root_domain
    integer :: root_patch
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: scalar_mult
    integer :: scalar_variable
    integer :: vector_mult
    integer :: vector_variable
    integer :: tke_level
    integer :: n_tke_level

    integer, allocatable :: global_seen(:)
    integer, allocatable :: local_seen(:)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("local block store is not ready")
    end if

    allocate(local_seen(size(block_catalog)))
    allocate(global_seen(size(block_catalog)))

    local_seen   = 0
    global_seen  = 0
    local_count  = n_local_blocks()
    local_weight = 0

    call get_block_field_layout( &
         expected_scalar_variable,expected_n_scalar_variable, &
         expected_vector_variable, &
         expected_field_level,expected_n_field_level, &
         expected_scalar_mult, &
         expected_vector_mult)

    call get_block_turbulence_layout( &
         expected_tke_level,expected_n_tke_level)

    do i = 1, local_count

       b = local_block_catalog(i)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("local block has an invalid catalogue index")
       end if

       if (local_seen(b) /= 0) then
          call fail("local block occurs more than once")
       end if

       if (catalog_local_block(b) /= i) then
          call fail("local/catalogue inverse mapping mismatch")
       end if

       if (block_catalog(b)%owner /= rank) then
          call fail("local block has the wrong final owner")
       end if

       call get_local_block_identity( &
            i,id,root_domain,root_patch,level)

       call get_local_block_field_layout( &
            i,scalar_variable,n_scalar_variable, &
            vector_variable,field_level, &
            n_field_level,scalar_mult,vector_mult)

       call get_local_block_turbulence_layout( &
            i,tke_level,n_tke_level)

       if (id /= block_catalog(b)%id .or. &
            root_domain /= block_catalog(b)%root_domain .or. &
            root_patch /= block_catalog(b)%root_patch .or. &
            level /= block_catalog(b)%level) then
          call fail("local block identity does not match catalogue")
       end if

       if (scalar_variable /= expected_scalar_variable .or. &
            n_scalar_variable /= expected_n_scalar_variable .or. &
            vector_variable /= expected_vector_variable .or. &
            field_level /= expected_field_level .or. &
            n_field_level /= expected_n_field_level .or. &
            scalar_mult /= expected_scalar_mult .or. &
            vector_mult /= expected_vector_mult) then
          call fail("local block field layout does not match format")
       end if

       if (tke_level /= expected_tke_level .or. &
            n_tke_level /= expected_n_tke_level) then
          call fail("local block turbulence layout does not match format")
       end if

       call check_local_block_storage(i,.true.)

       local_seen(b) = 1
       local_weight = local_weight + block_catalog(b)%weight

    end do

    expected_local = count(block_catalog%owner == rank)

    if (local_count /= expected_local) then
       call fail("local block count does not match final ownership")
    end if

    call MPI_Allreduce( &
         local_count,global_count,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block count")

    call MPI_Allreduce( &
         local_weight,global_weight,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block weight")

    call MPI_Allreduce( &
         local_seen,global_seen,size(local_seen),MPI_INTEGER, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce local block inventory")

    if (global_count /= size(block_catalog)) then
       call fail("global block count mismatch")
    end if

    if (global_weight /= sum(block_catalog%weight)) then
       call fail("global block weight mismatch")
    end if

    if (any(global_seen /= 1)) then
       call fail("global block ownership is not unique")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Standalone local block store for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  final-owner blocks = ", local_count
       write(6,'(a,i0)') &
            "  final-owner weight = ", local_weight
       write(6,'(a)') &
            "  component and serialization checks passed"
       write(6,'(a)') &
            "  persistent store readiness check passed"
       write(6,'(a)') &
            "  bidirectional local/catalogue mapping check passed"
       write(6,'(a)') &
            "  self-describing field metadata check passed"
       write(6,'(a,i0)') &
            "  scalar variables represented = ", &
            expected_n_scalar_variable
       write(6,'(a,i0)') &
            "  first sol field level = ", expected_field_level
       write(6,'(a,i0)') &
            "  sol field levels represented = ", &
            expected_n_field_level
       write(6,'(a,i0)') &
            "  tke field levels represented = ", &
            expected_n_tke_level
       write(6,'(a,/)') &
            "  unique global inventory check passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Standalone global block objects verified = ", &
            global_count
       write(6,'(a,i0)') &
            "Standalone global block weight verified  = ", &
            global_weight
       write(6,'(a,/)') &
            "Final-owner block store is self-contained"
    end if

    deallocate(global_seen)
    deallocate(local_seen)

  end subroutine check_local_blocks


  subroutine check_block_field_inventory (verbose)
    ! Compare sol, sol_mean, wav_coeff, tke, wav_tke and topography in the
    ! migrated final-owner block store with the still-authoritative
    ! legacy Domain fields covered by the catalogue-rooted subtrees. The
    ! fixed coarse scaffold above those roots is not block storage.
    ! Only global counts and order-independent moments are compared,
    ! because block and Domain ownership differ after migration.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: d
    integer :: ierr
    integer :: mult_scalar
    integer :: mult_vector
    integer :: n_field_level
    integer :: n_scalar_variable
    integer :: v_scalar
    integer :: v_vector
    integer :: k_field
    integer :: tke_level
    integer :: n_tke_level

    integer(int64) :: block_count_local(2)
    integer(int64) :: block_count_global(2)
    integer(int64) :: block_mean_count_local(2)
    integer(int64) :: block_mean_count_global(2)
    integer(int64) :: block_wavelet_count_local(2)
    integer(int64) :: block_wavelet_count_global(2)
    integer(int64) :: domain_count_local(2)
    integer(int64) :: domain_count_global(2)
    integer(int64) :: domain_mean_count_local(2)
    integer(int64) :: domain_mean_count_global(2)
    integer(int64) :: domain_wavelet_count_local(2)
    integer(int64) :: domain_wavelet_count_global(2)
    integer(int64) :: block_tke_count_local
    integer(int64) :: block_tke_count_global
    integer(int64) :: block_wav_tke_count_local
    integer(int64) :: block_wav_tke_count_global
    integer(int64) :: domain_tke_count_local
    integer(int64) :: domain_tke_count_global
    integer(int64) :: domain_wav_tke_count_local
    integer(int64) :: domain_wav_tke_count_global
    integer(int64) :: block_topography_count_local
    integer(int64) :: block_topography_count_global
    integer(int64) :: domain_topography_count_local
    integer(int64) :: domain_topography_count_global

    real(dp) :: block_moment_local(3,2)
    real(dp) :: block_moment_global(3,2)
    real(dp) :: block_mean_moment_local(3,2)
    real(dp) :: block_mean_moment_global(3,2)
    real(dp) :: block_wavelet_moment_local(3,2)
    real(dp) :: block_wavelet_moment_global(3,2)
    real(dp) :: domain_moment_local(3,2)
    real(dp) :: domain_moment_global(3,2)
    real(dp) :: domain_mean_moment_local(3,2)
    real(dp) :: domain_mean_moment_global(3,2)
    real(dp) :: domain_wavelet_moment_local(3,2)
    real(dp) :: domain_wavelet_moment_global(3,2)
    real(dp) :: block_tke_moment_local(3)
    real(dp) :: block_tke_moment_global(3)
    real(dp) :: block_wav_tke_moment_local(3)
    real(dp) :: block_wav_tke_moment_global(3)
    real(dp) :: domain_tke_moment_local(3)
    real(dp) :: domain_tke_moment_global(3)
    real(dp) :: domain_wav_tke_moment_local(3)
    real(dp) :: domain_wav_tke_moment_global(3)
    real(dp) :: block_topography_moment_local(3)
    real(dp) :: block_topography_moment_global(3)
    real(dp) :: domain_topography_moment_local(3)
    real(dp) :: domain_topography_moment_global(3)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_store_ready()) then
       call fail("field inventory requested before block installation")
    end if

    call local_block_field_statistics( &
         block_count_local(1),block_count_local(2), &
         block_moment_local(:,1),block_moment_local(:,2))

    call local_block_wavelet_statistics( &
         block_wavelet_count_local(1),block_wavelet_count_local(2), &
         block_wavelet_moment_local(:,1), &
         block_wavelet_moment_local(:,2))

    call local_block_mean_field_statistics( &
         block_mean_count_local(1),block_mean_count_local(2), &
         block_mean_moment_local(:,1),block_mean_moment_local(:,2))

    call local_block_turbulence_statistics( &
         block_tke_count_local,block_wav_tke_count_local, &
         block_tke_moment_local,block_wav_tke_moment_local)

    call local_block_topography_statistics( &
         block_topography_count_local,block_topography_moment_local)

    call get_block_field_layout( &
         v_scalar,n_scalar_variable,v_vector,k_field, &
         n_field_level,mult_scalar,mult_vector)

    call get_block_turbulence_layout(tke_level,n_tke_level)

    domain_count_local  = 0_int64
    domain_moment_local = 0.0_dp
    domain_mean_count_local  = 0_int64
    domain_mean_moment_local = 0.0_dp
    domain_wavelet_count_local  = 0_int64
    domain_wavelet_moment_local = 0.0_dp
    domain_tke_count_local      = 0_int64
    domain_wav_tke_count_local  = 0_int64
    domain_tke_moment_local     = 0.0_dp
    domain_wav_tke_moment_local = 0.0_dp
    domain_topography_count_local  = 0_int64
    domain_topography_moment_local = 0.0_dp

    do b = 1, size(block_catalog)

       if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

       d = loc_id(block_catalog(b)%root_domain+1) + 1

       if (d < 1 .or. d > size(grid)) then
          call fail("invalid source Domain in field inventory")
       end if

       call accumulate_domain_subtree_fields( &
            d,block_catalog(b)%root_patch, &
            v_scalar,n_scalar_variable,v_vector,k_field,n_field_level, &
            mult_scalar,mult_vector, &
            domain_count_local,domain_moment_local, &
            domain_mean_count_local,domain_mean_moment_local, &
            domain_wavelet_count_local,domain_wavelet_moment_local, &
            tke_level,n_tke_level, &
            domain_tke_count_local,domain_tke_moment_local, &
            domain_wav_tke_count_local,domain_wav_tke_moment_local, &
            domain_topography_count_local, &
            domain_topography_moment_local)

    end do

    call MPI_Allreduce( &
         block_count_local,block_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block field counts")

    call MPI_Allreduce( &
         block_mean_count_local,block_mean_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block mean-field counts")

    call MPI_Allreduce( &
         block_wavelet_count_local,block_wavelet_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wavelet counts")

    call MPI_Allreduce( &
         domain_count_local,domain_count_global,2,MPI_INTEGER8, &
         MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain field counts")

    call MPI_Allreduce( &
         domain_mean_count_local,domain_mean_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain mean-field counts")

    call MPI_Allreduce( &
         domain_wavelet_count_local,domain_wavelet_count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wavelet counts")

    call MPI_Allreduce( &
         block_tke_count_local,block_tke_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block tke count")

    call MPI_Allreduce( &
         block_wav_tke_count_local,block_wav_tke_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wav_tke count")

    call MPI_Allreduce( &
         domain_tke_count_local,domain_tke_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain tke count")

    call MPI_Allreduce( &
         domain_wav_tke_count_local,domain_wav_tke_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wav_tke count")

    call MPI_Allreduce( &
         block_topography_count_local,block_topography_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block topography count")

    call MPI_Allreduce( &
         domain_topography_count_local,domain_topography_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain topography count")

    call MPI_Allreduce( &
         block_moment_local,block_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block field moments")

    call MPI_Allreduce( &
         block_mean_moment_local,block_mean_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block mean-field moments")

    call MPI_Allreduce( &
         block_wavelet_moment_local,block_wavelet_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wavelet moments")

    call MPI_Allreduce( &
         domain_moment_local,domain_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain field moments")

    call MPI_Allreduce( &
         domain_mean_moment_local,domain_mean_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain mean-field moments")

    call MPI_Allreduce( &
         domain_wavelet_moment_local,domain_wavelet_moment_global,6, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wavelet moments")

    call MPI_Allreduce( &
         block_tke_moment_local,block_tke_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block tke moments")

    call MPI_Allreduce( &
         block_wav_tke_moment_local,block_wav_tke_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block wav_tke moments")

    call MPI_Allreduce( &
         domain_tke_moment_local,domain_tke_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain tke moments")

    call MPI_Allreduce( &
         domain_wav_tke_moment_local,domain_wav_tke_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain wav_tke moments")

    call MPI_Allreduce( &
         block_topography_moment_local,block_topography_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block topography moments")

    call MPI_Allreduce( &
         domain_topography_moment_local,domain_topography_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain topography moments")

    if (any(block_count_global /= domain_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain field-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block counts  = ", block_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain counts = ", domain_count_global
       end if

       call fail("block and Domain field counts differ")
    end if

    if (any(block_mean_count_global /= domain_mean_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean-field-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block mean counts  = ", block_mean_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain mean counts = ", domain_mean_count_global
       end if

       call fail("block and Domain mean-field counts differ")
    end if

    if (any(block_wavelet_count_global /= &
         domain_wavelet_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain wavelet-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block wavelet counts  = ", block_wavelet_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain wavelet counts = ", domain_wavelet_count_global
       end if

       call fail("block and Domain wavelet counts differ")
    end if

    if (block_tke_count_global /= domain_tke_count_global .or. &
         block_wav_tke_count_global /= domain_wav_tke_count_global) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain turbulence-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block tke/wav_tke counts  = ", &
               block_tke_count_global,block_wav_tke_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain tke/wav_tke counts = ", &
               domain_tke_count_global,domain_wav_tke_count_global
       end if

       call fail("block and Domain turbulence counts differ")
    end if

    if (block_topography_count_global /= &
         domain_topography_count_global) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain topography-count mismatch:"
          write(error_unit,'(a,i0)') &
               "  block topography count  = ", &
               block_topography_count_global
          write(error_unit,'(a,i0)') &
               "  Domain topography count = ", &
               domain_topography_count_global
       end if

       call fail("block and Domain topography counts differ")
    end if

    if (.not. field_moments_match( &
         block_moment_global(:,1),domain_moment_global(:,1), &
         block_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain scalar-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_moment_global(:,1)
       end if

       call fail("block and Domain scalar moments differ")
    end if

    if (.not. field_moments_match( &
         block_moment_global(:,2),domain_moment_global(:,2), &
         block_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain vector-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_moment_global(:,2)
       end if

       call fail("block and Domain vector moments differ")
    end if

    if (.not. field_moments_match( &
         block_mean_moment_global(:,1), &
         domain_mean_moment_global(:,1), &
         block_mean_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean scalar-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_mean_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_mean_moment_global(:,1)
       end if

       call fail("block and Domain mean scalar moments differ")
    end if

    if (.not. field_moments_match( &
         block_mean_moment_global(:,2), &
         domain_mean_moment_global(:,2), &
         block_mean_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain mean vector-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_mean_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_mean_moment_global(:,2)
       end if

       call fail("block and Domain mean vector moments differ")
    end if

    if (.not. field_moments_match( &
         block_wavelet_moment_global(:,1), &
         domain_wavelet_moment_global(:,1), &
         block_wavelet_count_global(1))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain scalar wavelet-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_wavelet_moment_global(:,1)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_wavelet_moment_global(:,1)
       end if

       call fail("block and Domain scalar wavelet moments differ")
    end if

    if (.not. field_moments_match( &
         block_wavelet_moment_global(:,2), &
         domain_wavelet_moment_global(:,2), &
         block_wavelet_count_global(2))) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain vector wavelet-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_wavelet_moment_global(:,2)
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_wavelet_moment_global(:,2)
       end if

       call fail("block and Domain vector wavelet moments differ")
    end if

    if (.not. field_moments_match( &
         block_tke_moment_global,domain_tke_moment_global, &
         block_tke_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') "Block/Domain tke-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_tke_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_tke_moment_global
       end if

       call fail("block and Domain tke moments differ")
    end if

    if (.not. field_moments_match( &
         block_wav_tke_moment_global,domain_wav_tke_moment_global, &
         block_wav_tke_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain wav_tke-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_wav_tke_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_wav_tke_moment_global
       end if

       call fail("block and Domain wav_tke moments differ")
    end if

    if (.not. field_moments_match( &
         block_topography_moment_global, &
         domain_topography_moment_global, &
         block_topography_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain topography-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_topography_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_topography_moment_global
       end if

       call fail("block and Domain topography moments differ")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Read-only block field consumer for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  local block scalar values = ", block_count_local(1)
       write(6,'(a,i0)') &
            "  local block vector values = ", block_count_local(2)
       write(6,'(a,i0)') &
            "  local block wavelet scalar values = ", &
            block_wavelet_count_local(1)
       write(6,'(a,i0)') &
            "  local block wavelet vector values = ", &
            block_wavelet_count_local(2)
       write(6,'(a,i0)') &
            "  local block mean scalar values = ", &
            block_mean_count_local(1)
       write(6,'(a,i0)') &
            "  local block mean vector values = ", &
            block_mean_count_local(2)
       write(6,'(a,i0)') &
            "  local block topography values = ", &
            block_topography_count_local
       if (n_tke_level > 0) then
          write(6,'(a,i0)') &
               "  local block tke values = ", block_tke_count_local
          write(6,'(a,i0)') &
               "  local block wav_tke values = ", &
               block_wav_tke_count_local
       end if
       write(6,'(a,/)') &
            "  global sol/sol_mean/wav_coeff inventory checks passed"
       write(6,'(a,/)') &
            "  global topography inventory check passed"
       if (n_tke_level > 0) then
          write(6,'(a,/)') &
               "  global tke/wav_tke inventory checks passed"
       end if
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global scalar interior values verified = ", &
            block_count_global(1)
       write(6,'(a,i0)') &
            "Global vector interior values verified = ", &
            block_count_global(2)
       write(6,'(a,i0)') &
            "Global wavelet scalar values verified = ", &
            block_wavelet_count_global(1)
       write(6,'(a,i0)') &
            "Global wavelet vector values verified = ", &
            block_wavelet_count_global(2)
       write(6,'(a,i0)') &
            "Global mean scalar interior values verified = ", &
            block_mean_count_global(1)
       write(6,'(a,i0)') &
            "Global mean vector interior values verified = ", &
            block_mean_count_global(2)
       write(6,'(a,i0)') &
            "Global topography values verified = ", &
            block_topography_count_global
       write(6,'(a,3(es24.16,1x))') &
            "Global topography moments verified = ", &
            block_topography_moment_global
       if (n_tke_level > 0) then
          write(6,'(a,i0)') &
               "Global tke interior values verified = ", &
               block_tke_count_global
          write(6,'(a,i0)') &
               "Global wav_tke values verified = ", &
               block_wav_tke_count_global
       end if
       write(6,'(a,/)') &
            "Block sol, sol_mean and wav_coeff match legacy Domain data"
       write(6,'(a,/)') &
            "Block topography matches legacy Domain data"
       if (n_tke_level > 0) then
          write(6,'(a,/)') &
               "Block tke and wav_tke match legacy Domain data"
       end if
    end if

  end subroutine check_block_field_inventory


  subroutine check_block_hydrostatic_state_accessors (verbose)
    ! Verify patch and whole-block views of persistent thermodynamic
    ! fields, dependency-selective invalidation and lazy cache refresh.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: catalog_index
    integer :: clean_catalog_index
    integer :: column_base
    integer :: column_nvalue
    integer :: local_index
    integer :: local_patch
    integer :: n_patch
    integer :: scalar_nvalue
    integer :: surface_base
    integer :: surface_nvalue

    integer(int64) :: column_count_after
    integer(int64) :: column_count_before
    integer(int64) :: clean_refresh_before
    integer(int64) :: dirty_refresh_before
    integer(int64) :: refresh_count_before
    integer(int64) :: refresh_count_expected
    integer(int64) :: surface_count_after
    integer(int64) :: surface_count_before

    real(dp), allocatable :: air_temperature(:)
    real(dp), allocatable :: dynamic_exner(:)
    real(dp), allocatable :: scalar_snapshot(:)
    real(dp), allocatable :: surface_pressure(:)

    real(dp) :: air_temperature_patch(zlevels*PATCH_SIZE**2)
    real(dp) :: dynamic_exner_patch(zlevels*PATCH_SIZE**2)
    real(dp) :: surface_pressure_patch(PATCH_SIZE**2)
    real(dp) :: exner_moment_after(3)
    real(dp) :: exner_moment_before(3)
    real(dp) :: surface_moment_after(3)
    real(dp) :: surface_moment_before(3)
    real(dp) :: temperature_moment_after(3)
    real(dp) :: temperature_moment_before(3)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. compressible) then
       if (print_summary) then
          write(6,'(/,a,i0,a)') &
               "Persistent hydrostatic accessors for rank ", rank, ":"
          write(6,'(a,/)') &
               "  skipped for incompressible configuration"
       end if
       return
    end if

    if (.not. local_block_hydrostatic_state_ready()) then
       call fail("hydrostatic accessor check before state refresh")
    end if

    call local_block_hydrostatic_statistics( &
         surface_count_before,column_count_before, &
         surface_moment_before,exner_moment_before, &
         temperature_moment_before)
    refresh_count_before = local_block_hydrostatic_refresh_count()

    do local_index = 1,n_local_blocks()
       catalog_index = local_block_catalog(local_index)
       n_patch = local_block_patch_count(catalog_index)
       surface_nvalue = &
            local_block_hydrostatic_surface_nvalue(catalog_index)
       column_nvalue = &
            local_block_hydrostatic_column_nvalue(catalog_index)

       if (surface_nvalue /= n_patch*PATCH_SIZE**2 .or. &
            column_nvalue /= n_patch*zlevels*PATCH_SIZE**2) then
          call fail("whole-block hydrostatic extent mismatch")
       end if

       allocate(surface_pressure(surface_nvalue))
       allocate(dynamic_exner(column_nvalue))
       allocate(air_temperature(column_nvalue))

       call get_local_block_hydrostatic_values( &
            catalog_index,surface_pressure,dynamic_exner, &
            air_temperature)

       do local_patch = 0,n_patch-1
          call get_local_block_hydrostatic_patch_values( &
               catalog_index,local_patch,surface_pressure_patch, &
               dynamic_exner_patch,air_temperature_patch)

          surface_base = local_patch*PATCH_SIZE**2
          column_base = local_patch*zlevels*PATCH_SIZE**2

          if (any(abs(surface_pressure( &
               surface_base+1:surface_base+PATCH_SIZE**2) - &
               surface_pressure_patch) > 0.0_dp)) then
             call fail("surface-pressure patch/block view mismatch")
          end if
          if (any(abs(dynamic_exner( &
               column_base+1:column_base+zlevels*PATCH_SIZE**2) - &
               dynamic_exner_patch) > 0.0_dp)) then
             call fail("dynamic-Exner patch/block view mismatch")
          end if
          if (any(abs(air_temperature( &
               column_base+1:column_base+zlevels*PATCH_SIZE**2) - &
               air_temperature_patch) > 0.0_dp)) then
             call fail("temperature patch/block view mismatch")
          end if
       end do

       deallocate(surface_pressure)
       deallocate(dynamic_exner)
       deallocate(air_temperature)
    end do

    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_before) then
       call fail("read access caused redundant hydrostatic refresh")
    end if

    if (n_local_blocks() > 0) then
       catalog_index = local_block_catalog(1)
       surface_nvalue = &
            local_block_hydrostatic_surface_nvalue(catalog_index)
       column_nvalue = &
            local_block_hydrostatic_column_nvalue(catalog_index)
       dirty_refresh_before = &
            local_block_hydrostatic_block_refresh_count(catalog_index)

       allocate(surface_pressure(surface_nvalue))
       allocate(dynamic_exner(column_nvalue))
       allocate(air_temperature(column_nvalue))

       clean_refresh_before = 0_int64
       if (n_local_blocks() > 1) then
          clean_catalog_index = local_block_catalog(2)
          clean_refresh_before = &
               local_block_hydrostatic_block_refresh_count( &
               clean_catalog_index)
       end if

       scalar_nvalue = &
            local_block_scalar_family_patch_nvalue(catalog_index)
       allocate(scalar_snapshot(scalar_nvalue))

       call get_local_block_scalar_patch_family_values( &
            catalog_index,0,BLOCK_PAYLOAD_WAV_COEFF,scalar_snapshot)
       call set_local_block_scalar_patch_family_values( &
            catalog_index,0,BLOCK_PAYLOAD_WAV_COEFF,scalar_snapshot)

       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("wavelet write invalidated hydrostatic state")
       end if
       if (local_block_hydrostatic_refresh_count() /= &
            refresh_count_before) then
          call fail("wavelet write refreshed hydrostatic state")
       end if

       call get_local_block_scalar_patch_family_values( &
            catalog_index,0,BLOCK_PAYLOAD_SOL,scalar_snapshot)
       call set_local_block_scalar_patch_family_values( &
            catalog_index,0,BLOCK_PAYLOAD_SOL,scalar_snapshot)

       if (local_block_hydrostatic_state_ready()) then
          call fail("scalar sol write did not invalidate hydrostatic state")
       end if
       if (local_block_hydrostatic_refresh_count() /= &
            refresh_count_before) then
          call fail("scalar sol write eagerly refreshed hydrostatic state")
       end if

       if (n_local_blocks() > 1) then
          surface_nvalue = &
               local_block_hydrostatic_surface_nvalue( &
               clean_catalog_index)
          if (local_block_hydrostatic_state_ready()) then
             call fail("clean-block access refreshed a dirty block")
          end if
          if (local_block_hydrostatic_refresh_count() /= &
               refresh_count_before) then
             call fail("clean-block access changed global refresh count")
          end if
          if (local_block_hydrostatic_block_refresh_count( &
               clean_catalog_index) /= clean_refresh_before) then
             call fail("clean block was redundantly refreshed")
          end if
       end if

       call get_local_block_hydrostatic_values( &
            catalog_index,surface_pressure,dynamic_exner, &
            air_temperature)

       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("dirty-block access did not restore cache readiness")
       end if
       if (local_block_hydrostatic_block_refresh_count(catalog_index) /= &
            dirty_refresh_before+1_int64) then
          call fail("dirty block refresh count mismatch")
       end if
       if (n_local_blocks() > 1) then
          if (local_block_hydrostatic_block_refresh_count( &
               clean_catalog_index) /= clean_refresh_before) then
             call fail("dirty-block access refreshed a clean block")
          end if
       end if

       deallocate(scalar_snapshot)
       deallocate(surface_pressure)
       deallocate(dynamic_exner)
       deallocate(air_temperature)
    end if

    ! A whole-store consumer must not refresh blocks already made current.
    call local_block_hydrostatic_statistics( &
         surface_count_after,column_count_after, &
         surface_moment_after,exner_moment_after, &
         temperature_moment_after)

    refresh_count_expected = refresh_count_before
    if (n_local_blocks() > 0) then
       refresh_count_expected = refresh_count_expected + 1_int64
    end if

    if (.not. local_block_hydrostatic_state_ready()) then
       call fail("hydrostatic state not ready after lazy refresh")
    end if
    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_expected) then
       call fail("lazy hydrostatic refresh count mismatch")
    end if

    call ensure_local_block_hydrostatic_state
    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_expected) then
       call fail("redundant hydrostatic refresh was not suppressed")
    end if

    if (surface_count_after /= surface_count_before .or. &
         column_count_after /= column_count_before) then
       call fail("hydrostatic refresh changed local value counts")
    end if
    if (.not. field_moments_match( &
         surface_moment_after,surface_moment_before, &
         surface_count_before)) then
       call fail("hydrostatic refresh changed surface-pressure values")
    end if
    if (.not. field_moments_match( &
         exner_moment_after,exner_moment_before,column_count_before)) then
       call fail("hydrostatic refresh changed dynamic Exner values")
    end if
    if (.not. field_moments_match( &
         temperature_moment_after,temperature_moment_before, &
         column_count_before)) then
       call fail("hydrostatic refresh changed temperature values")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Persistent hydrostatic accessors for rank ", rank, ":"
       write(6,'(a)') &
            "  exact patch/whole-block layout checks passed"
       write(6,'(a)') &
            "  wavelet-independent cache retention checks passed"
       write(6,'(a)') &
            "  selective dirty-block refresh checks passed"
       write(6,'(a,/)') &
            "  per-block lazy cache coherence checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Selective per-block thermodynamic caches passed"
    end if

  end subroutine check_block_hydrostatic_state_accessors


  subroutine accumulate_block_field_consumer (catalog_index,block,context)
    ! Production prognostic-field kernel. Consume the complete read-only
    ! block view, validate its topology/storage contract and accumulate
    ! interior field inventories without accessor copies. For compressible
    ! cases, also diagnose surface pressure directly from mass fields.

    implicit none

    integer, intent(in) :: catalog_index
    type(Block_Data), intent(in) :: block
    class(*), intent(inout) :: context

    integer :: expected_scalar
    integer :: expected_vector
    integer :: i
    integer :: k
    integer :: mass_index
    integer :: mass_slot
    integer :: n_node
    integer :: scalar_variable_size
    integer :: temperature_slot

    real(dp) :: pressure
    real(dp) :: rho_dz

    if (catalog_index < 1) then
       call fail("field consumer received invalid catalogue index")
    end if
    if (.not. allocated(block%patch) .or. &
         .not. allocated(block%node) .or. &
         .not. allocated(block%bdry_node) .or. &
         .not. allocated(block%ghost_node) .or. &
         .not. allocated(block%stencil) .or. &
         .not. allocated(block%block_bdry) .or. &
         .not. allocated(block%bdry_storage) .or. &
         .not. allocated(block%ghost_storage)) then
       call fail("field consumer received incomplete topology")
    end if

    n_node = size(block%node)
    if (n_node /= size(block%patch)*PATCH_SIZE**2 .or. &
         size(block%stencil,1) /= N_BDRY .or. &
         size(block%stencil,2) /= size(block%patch)) then
       call fail("field consumer received invalid patch topology")
    end if
    if (sum(block%bdry_storage%n_node) /= size(block%bdry_node) .or. &
         sum(block%ghost_storage%n_node) /= size(block%ghost_node)) then
       call fail("field consumer received invalid compact topology")
    end if

    expected_scalar = block%n_scalar_variable*block%n_field_level* &
         block%scalar_mult*n_node
    expected_vector = block%n_field_level*block%vector_mult*n_node

    if (size(block%scalar) /= expected_scalar .or. &
         size(block%scalar_mean) /= expected_scalar .or. &
         size(block%wavelet_scalar) /= expected_scalar .or. &
         size(block%vector) /= expected_vector .or. &
         size(block%vector_mean) /= expected_vector .or. &
         size(block%wavelet_vector) /= expected_vector) then
       call fail("field consumer received invalid interior fields")
    end if

    expected_scalar = block%n_scalar_variable*block%n_field_level* &
         block%scalar_mult*size(block%bdry_node)
    expected_vector = block%n_field_level*block%vector_mult* &
         size(block%bdry_node)
    if (size(block%bdry_scalar) /= expected_scalar .or. &
         size(block%bdry_scalar_mean) /= expected_scalar .or. &
         size(block%bdry_wavelet_scalar) /= expected_scalar .or. &
         size(block%bdry_vector) /= expected_vector .or. &
         size(block%bdry_vector_mean) /= expected_vector .or. &
         size(block%bdry_wavelet_vector) /= expected_vector) then
       call fail("field consumer received invalid boundary fields")
    end if

    expected_scalar = block%n_scalar_variable*block%n_field_level* &
         block%scalar_mult*size(block%ghost_node)
    expected_vector = block%n_field_level*block%vector_mult* &
         size(block%ghost_node)
    if (size(block%ghost_scalar) /= expected_scalar .or. &
         size(block%ghost_scalar_mean) /= expected_scalar .or. &
         size(block%ghost_wavelet_scalar) /= expected_scalar .or. &
         size(block%ghost_vector) /= expected_vector .or. &
         size(block%ghost_vector_mean) /= expected_vector .or. &
         size(block%ghost_wavelet_vector) /= expected_vector) then
       call fail("field consumer received invalid ghost fields")
    end if

    select type (statistics => context)
    type is (Block_Field_Traversal_Context)
       statistics%block_count = statistics%block_count + 1_int64
       statistics%patch_count = statistics%patch_count + &
            int(size(block%patch),int64)
       statistics%boundary_count = statistics%boundary_count + &
            int(size(block%bdry_storage),int64)
       statistics%ghost_count = statistics%ghost_count + &
            int(size(block%ghost_storage),int64)
       statistics%node_count = statistics%node_count + int(n_node,int64)
       statistics%boundary_node_count = &
            statistics%boundary_node_count + &
            int(size(block%bdry_node),int64)
       statistics%ghost_node_count = statistics%ghost_node_count + &
            int(size(block%ghost_node),int64)

       statistics%scalar_count = statistics%scalar_count + &
            int(size(block%scalar),int64)
       statistics%vector_count = statistics%vector_count + &
            int(size(block%vector),int64)

       statistics%scalar_moment(:,1) = &
            statistics%scalar_moment(:,1) + [ &
            sum(block%scalar),sum(abs(block%scalar)), &
            sum(block%scalar**2) ]
       statistics%scalar_moment(:,2) = &
            statistics%scalar_moment(:,2) + [ &
            sum(block%scalar_mean),sum(abs(block%scalar_mean)), &
            sum(block%scalar_mean**2) ]
       statistics%scalar_moment(:,3) = &
            statistics%scalar_moment(:,3) + [ &
            sum(block%wavelet_scalar), &
            sum(abs(block%wavelet_scalar)), &
            sum(block%wavelet_scalar**2) ]

       statistics%vector_moment(:,1) = &
            statistics%vector_moment(:,1) + [ &
            sum(block%vector),sum(abs(block%vector)), &
            sum(block%vector**2) ]
       statistics%vector_moment(:,2) = &
            statistics%vector_moment(:,2) + [ &
            sum(block%vector_mean),sum(abs(block%vector_mean)), &
            sum(block%vector_mean**2) ]
       statistics%vector_moment(:,3) = &
            statistics%vector_moment(:,3) + [ &
            sum(block%wavelet_vector), &
            sum(abs(block%wavelet_vector)), &
            sum(block%wavelet_vector**2) ]

       if (compressible) then
          if (block%scalar_mult /= 1 .or. &
               block%field_level > 1 .or. &
               block%field_level+block%n_field_level-1 < zlevels) then
             call fail("field consumer cannot diagnose surface pressure")
          end if

          mass_slot = S_MASS - block%scalar_variable
          temperature_slot = S_TEMP - block%scalar_variable
          if (mass_slot < 0 .or. &
               mass_slot >= block%n_scalar_variable .or. &
               temperature_slot < 0 .or. &
               temperature_slot >= block%n_scalar_variable) then
             call fail("field consumer lacks thermodynamic variables")
          end if

          scalar_variable_size = block%n_field_level*n_node
          do i = 1,n_node
             pressure = p_top
             do k = 1,zlevels
                mass_index = mass_slot*scalar_variable_size + &
                     (k-block%field_level)*n_node + i
                rho_dz = block%scalar(mass_index) + &
                     block%scalar_mean(mass_index)
                if (rho_dz <= 0.0_dp) then
                   call fail("field consumer diagnosed nonpositive mass")
                end if
                pressure = pressure + grav_accel*rho_dz
             end do

             statistics%surface_count = &
                  statistics%surface_count + 1_int64
             statistics%surface_moment(1) = &
                  statistics%surface_moment(1) + pressure
             statistics%surface_moment(2) = &
                  statistics%surface_moment(2) + abs(pressure)
             statistics%surface_moment(3) = &
                  statistics%surface_moment(3) + pressure**2
          end do
       end if
    class default
       call fail("field consumer received invalid context")
    end select

  end subroutine accumulate_block_field_consumer


  subroutine check_block_field_consumer (verbose)
    ! Exercise the generalized production traversal, compare direct field
    ! inventories with established consumers and validate the first compact
    ! prognostic-to-surface-pressure diagnostic in shadow mode.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: catalog_index
    integer :: local_index

    integer(int64) :: expected_boundary_count
    integer(int64) :: expected_ghost_count
    integer(int64) :: expected_patch_count
    integer(int64) :: column_count
    integer(int64) :: field_count(2)
    integer(int64) :: mean_count(2)
    integer(int64) :: refresh_count_before
    integer(int64) :: surface_count
    integer(int64) :: wavelet_count(2)

    real(dp) :: column_moment(3)
    real(dp) :: field_moment(3,2)
    real(dp) :: mean_moment(3,2)
    real(dp) :: surface_moment(3)
    real(dp) :: temperature_moment(3)
    real(dp) :: wavelet_moment(3,2)

    logical :: print_summary

    type(Block_Field_Traversal_Context) :: statistics

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_mean_field_statistics( &
         mean_count(1),mean_count(2), &
         mean_moment(:,1),mean_moment(:,2))
    call local_block_wavelet_statistics( &
         wavelet_count(1),wavelet_count(2), &
         wavelet_moment(:,1),wavelet_moment(:,2))

    surface_count = 0_int64
    surface_moment = 0.0_dp
    if (compressible) then
       call local_block_hydrostatic_statistics( &
            surface_count,column_count,surface_moment, &
            column_moment,temperature_moment)
    end if

    refresh_count_before = local_block_hydrostatic_refresh_count()
    statistics = Block_Field_Traversal_Context()
    call apply_local_block_field_consumer( &
         accumulate_block_field_consumer,statistics)

    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_before) then
       call fail("production field consumer refreshed hydrostatic cache")
    end if

    expected_patch_count = 0_int64
    expected_boundary_count = 0_int64
    expected_ghost_count = 0_int64
    do local_index = 1,n_local_blocks()
       catalog_index = local_block_catalog(local_index)
       expected_patch_count = expected_patch_count + &
            int(local_block_patch_count(catalog_index),int64)
       expected_boundary_count = expected_boundary_count + &
            int(local_block_boundary_count(catalog_index),int64)
       expected_ghost_count = expected_ghost_count + &
            int(local_block_ghost_count(catalog_index),int64)
    end do

    if (statistics%block_count /= int(n_local_blocks(),int64) .or. &
         statistics%patch_count /= expected_patch_count .or. &
         statistics%boundary_count /= expected_boundary_count .or. &
         statistics%ghost_count /= expected_ghost_count .or. &
         statistics%node_count /= expected_patch_count*PATCH_SIZE**2) then
       call fail("production field consumer topology mismatch")
    end if

    if (statistics%scalar_count(1) /= field_count(1) .or. &
         statistics%scalar_count(2) /= mean_count(1) .or. &
         statistics%scalar_count(3) /= wavelet_count(1) .or. &
         statistics%vector_count(1) /= field_count(2) .or. &
         statistics%vector_count(2) /= mean_count(2) .or. &
         statistics%vector_count(3) /= wavelet_count(2)) then
       call fail("production field consumer inventory count mismatch")
    end if

    if (.not. field_moments_match( &
         statistics%scalar_moment(:,1),field_moment(:,1), &
         field_count(1)) .or. &
         .not. field_moments_match( &
         statistics%scalar_moment(:,2),mean_moment(:,1), &
         mean_count(1)) .or. &
         .not. field_moments_match( &
         statistics%scalar_moment(:,3),wavelet_moment(:,1), &
         wavelet_count(1)) .or. &
         .not. field_moments_match( &
         statistics%vector_moment(:,1),field_moment(:,2), &
         field_count(2)) .or. &
         .not. field_moments_match( &
         statistics%vector_moment(:,2),mean_moment(:,2), &
         mean_count(2)) .or. &
         .not. field_moments_match( &
         statistics%vector_moment(:,3),wavelet_moment(:,2), &
         wavelet_count(2))) then
       call fail("production field consumer inventory moment mismatch")
    end if

    if (compressible) then
       if (statistics%surface_count /= surface_count) then
          call fail("production surface-pressure count mismatch")
       end if
       if (column_count /= int(zlevels,int64)*surface_count) then
          call fail("production hydrostatic column count mismatch")
       end if
       if (.not. field_moments_match( &
            statistics%surface_moment,surface_moment,surface_count)) then
          call fail("production surface-pressure diagnostic mismatch")
       end if
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Block prognostic traversal for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  local blocks consumed       = ", statistics%block_count
       write(6,'(a,i0)') &
            "  local patches consumed      = ", statistics%patch_count
       write(6,'(a,i0)') &
            "  local boundary records seen = ", statistics%boundary_count
       write(6,'(a,i0)') &
            "  local ghost records seen    = ", statistics%ghost_count
       write(6,'(a)') &
            "  direct sol/sol_mean/wav_coeff traversal passed"
       write(6,'(a)') &
            "  topology and boundary/ghost views passed"
       if (compressible) then
          write(6,'(a,/)') &
               "  prognostic surface-pressure shadow kernel passed"
       end if
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Production block prognostic-field traversal passed"
    end if

  end subroutine check_block_field_consumer


  subroutine accumulate_block_stencil_kernel (catalog_index,block,context)
    ! Read the sol family through every compact stencil address and form
    ! neighbour-value and neighbour-minus-centre diagnostics. Boundary and
    ! ghost addresses use the same topology that production dynamics sees.

    implicit none

    integer, intent(in) :: catalog_index
    type(Block_Data), intent(in) :: block
    class(*), intent(inout) :: context

    integer :: address_offset
    integer :: center_base
    integer :: center_index
    integer :: center_node
    integer :: component_slot
    integer :: field_base
    integer :: field_index
    integer :: level_slot
    integer :: n_storage_node
    integer :: node_index
    integer :: p
    integer :: q
    integer :: record
    integer :: scalar_slot
    integer :: side
    integer :: storage_class
    integer :: storage_start

    real(dp) :: center_value
    real(dp) :: difference
    real(dp) :: value

    if (catalog_index < 1) then
       call fail("stencil kernel received invalid catalogue index")
    end if
    if (block%scalar_mult /= 1 .or. block%vector_mult < 1) then
       call fail("stencil kernel received invalid field multipliers")
    end if
    if (size(block%stencil,1) /= N_BDRY .or. &
         size(block%stencil,2) /= size(block%patch)) then
       call fail("stencil kernel received invalid topology")
    end if

    select type (statistics => context)
    type is (Block_Stencil_Kernel_Context)
       statistics%block_count = statistics%block_count + 1_int64

       do p = 1,size(block%patch)
          do side = 1,N_BDRY
             storage_class = block%stencil(side,p)%storage
             record = block%stencil(side,p)%id
             address_offset = block%stencil(side,p)%offset
             storage_start = 0
             n_storage_node = 0

             select case (storage_class)
             case (STORE_PATCH)
                if (record < 0 .or. record >= size(block%patch)) then
                   call fail("stencil kernel received invalid patch address")
                end if
                storage_start = block%patch(record+1)%elts_start
                n_storage_node = PATCH_SIZE**2
             case (STORE_BDRY)
                if (record < 1 .or. &
                     record > size(block%bdry_storage)) then
                   call fail( &
                        "stencil kernel received invalid boundary address")
                end if
                storage_start = block%bdry_storage(record)%local_start
                n_storage_node = block%bdry_storage(record)%n_node
             case (STORE_GHOST)
                if (record < 1 .or. &
                     record > size(block%ghost_storage)) then
                   call fail( &
                        "stencil kernel received invalid ghost address")
                end if
                storage_start = block%ghost_storage(record)%local_start
                n_storage_node = block%ghost_storage(record)%n_node
             case default
                call fail("stencil kernel received invalid storage class")
             end select

             do q = 0,PATCH_SIZE**2-1
                if (address_offset+q < 0 .or. &
                     address_offset+q >= n_storage_node) cycle

                node_index = storage_start + address_offset + q
                center_node = block%patch(p)%elts_start + q
                if (center_node < 0 .or. &
                     center_node >= size(block%node)) then
                   call fail("stencil kernel received invalid centre node")
                end if

                select case (storage_class)
                case (STORE_PATCH)
                   if (node_index < 0 .or. &
                        node_index >= size(block%node)) then
                      call fail("stencil kernel patch node is invalid")
                   end if
                case (STORE_BDRY)
                   if (node_index < 0 .or. &
                        node_index >= size(block%bdry_node)) then
                      call fail("stencil kernel boundary node is invalid")
                   end if
                case (STORE_GHOST)
                   if (node_index < 0 .or. &
                        node_index >= size(block%ghost_node)) then
                      call fail("stencil kernel ghost node is invalid")
                   end if
                end select

                statistics%address_count(storage_class) = &
                     statistics%address_count(storage_class) + 1_int64

                do scalar_slot = 1,block%n_scalar_variable
                   do level_slot = 1,block%n_field_level
                      center_base = &
                           ((scalar_slot-1)*block%n_field_level + &
                           level_slot-1)*size(block%node)
                      center_index = center_base + center_node + 1
                      center_value = block%scalar(center_index)
                      value = 0.0_dp

                      select case (storage_class)
                      case (STORE_PATCH)
                         field_base = center_base
                         field_index = field_base + node_index + 1
                         value = block%scalar(field_index)
                      case (STORE_BDRY)
                         field_base = &
                              ((scalar_slot-1)*block%n_field_level + &
                              level_slot-1)*size(block%bdry_node)
                         field_index = field_base + node_index + 1
                         value = block%bdry_scalar(field_index)
                      case (STORE_GHOST)
                         field_base = &
                              ((scalar_slot-1)*block%n_field_level + &
                              level_slot-1)*size(block%ghost_node)
                         field_index = field_base + node_index + 1
                         value = block%ghost_scalar(field_index)
                      end select

                      difference = value-center_value
                      statistics%scalar_count = &
                           statistics%scalar_count + 1_int64
                      statistics%scalar_moment = &
                           statistics%scalar_moment + &
                           [value,abs(value),value**2]
                      statistics%scalar_difference_moment = &
                           statistics%scalar_difference_moment + &
                           [difference,abs(difference),difference**2]
                   end do
                end do

                do level_slot = 1,block%n_field_level
                   do component_slot = 1,block%vector_mult
                      center_base = (level_slot-1)* &
                           block%vector_mult*size(block%node)
                      center_index = center_base + &
                           block%vector_mult*center_node + component_slot
                      center_value = block%vector(center_index)
                      value = 0.0_dp

                      select case (storage_class)
                      case (STORE_PATCH)
                         field_base = center_base
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%vector(field_index)
                      case (STORE_BDRY)
                         field_base = (level_slot-1)* &
                              block%vector_mult*size(block%bdry_node)
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%bdry_vector(field_index)
                      case (STORE_GHOST)
                         field_base = (level_slot-1)* &
                              block%vector_mult*size(block%ghost_node)
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%ghost_vector(field_index)
                      end select

                      difference = value-center_value
                      statistics%vector_count = &
                           statistics%vector_count + 1_int64
                      statistics%vector_moment = &
                           statistics%vector_moment + &
                           [value,abs(value),value**2]
                      statistics%vector_difference_moment = &
                           statistics%vector_difference_moment + &
                           [difference,abs(difference),difference**2]
                   end do
                end do
             end do
          end do
       end do
    class default
       call fail("stencil kernel received invalid context")
    end select

  end subroutine accumulate_block_stencil_kernel


  subroutine check_block_stencil_kernel (verbose)
    ! Refresh production sol ghosts and validate a stencil-dependent block
    ! kernel after migration staging has been released. A second refresh
    ! must reproduce the complete neighbour-difference diagnostic.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: count_global(5)
    integer(int64) :: count_local(5)
    integer(int64) :: refresh_count_before
    integer(int64) :: scalar_address_count(3)
    integer(int64) :: scalar_value_count
    integer(int64) :: vector_address_count(3)
    integer(int64) :: vector_value_count

    real(dp) :: scalar_value_moment(3,3)
    real(dp) :: vector_value_moment(3,3)

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: first
    type(Block_Stencil_Kernel_Context) :: second

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    refresh_count_before = local_block_hydrostatic_refresh_count()

    call refresh_block_sol_ghosts
    first = Block_Stencil_Kernel_Context()
    call apply_local_block_field_consumer( &
         accumulate_block_stencil_kernel,first)

    call local_block_scalar_stencil_statistics( &
         scalar_address_count,scalar_value_count,scalar_value_moment)
    call local_block_vector_stencil_statistics( &
         vector_address_count,vector_value_count,vector_value_moment)

    if (first%block_count /= int(n_local_blocks(),int64)) then
       call fail("production stencil kernel block-count mismatch")
    end if
    if (any(first%address_count /= scalar_address_count) .or. &
         any(first%address_count /= vector_address_count)) then
       call fail("production stencil kernel address mismatch")
    end if
    if (first%scalar_count /= scalar_value_count .or. &
         first%vector_count /= vector_value_count) then
       call fail("production stencil kernel value-count mismatch")
    end if
    if (.not. field_moments_match( &
         first%scalar_moment,scalar_value_moment(:,1), &
         scalar_value_count) .or. &
         .not. field_moments_match( &
         first%vector_moment,vector_value_moment(:,1), &
         vector_value_count)) then
       call fail("production stencil kernel sol moment mismatch")
    end if

    call refresh_block_sol_ghosts
    second = Block_Stencil_Kernel_Context()
    call apply_local_block_field_consumer( &
         accumulate_block_stencil_kernel,second)

    if (second%block_count /= first%block_count .or. &
         any(second%address_count /= first%address_count) .or. &
         second%scalar_count /= first%scalar_count .or. &
         second%vector_count /= first%vector_count) then
       call fail("repeated production stencil traversal changed counts")
    end if
    if (.not. field_moments_match( &
         second%scalar_moment,first%scalar_moment, &
         first%scalar_count) .or. &
         .not. field_moments_match( &
         second%vector_moment,first%vector_moment, &
         first%vector_count) .or. &
         .not. field_moments_match( &
         second%scalar_difference_moment, &
         first%scalar_difference_moment,first%scalar_count) .or. &
         .not. field_moments_match( &
         second%vector_difference_moment, &
         first%vector_difference_moment,first%vector_count)) then
       call fail("repeated production stencil diagnostic changed")
    end if

    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_before) then
       call fail("production stencil kernel refreshed hydrostatic cache")
    end if

    count_local(1:3) = first%address_count
    count_local(4) = first%scalar_count
    count_local(5) = first%vector_count
    call MPI_Allreduce( &
         count_local,count_global,5,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce production stencil counts")

    if (any(count_global <= 0_int64)) then
       call fail("production stencil kernel global inventory incomplete")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Production sol stencil kernel for rank ",rank,":"
       write(6,'(a,3(i0,1x))') &
            "  patch/boundary/ghost addresses = ",first%address_count
       write(6,'(a,i0)') &
            "  scalar neighbour samples       = ",first%scalar_count
       write(6,'(a,i0)') &
            "  vector neighbour samples       = ",first%vector_count
       write(6,'(a)') &
            "  direct neighbour-value moments passed"
       write(6,'(a)') &
            "  neighbour-minus-centre diagnostics passed"
       write(6,'(a,/)') &
            "  repeated sol ghost refresh is stencil-stable"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,3(i0,1x))') &
            "Global production patch/boundary/ghost addresses = ", &
            count_global(1:3)
       write(6,'(a,i0)') &
            "Global production scalar neighbour samples = ", &
            count_global(4)
       write(6,'(a,i0)') &
            "Global production vector neighbour samples = ", &
            count_global(5)
       write(6,'(a,/)') &
            "Production block sol stencil kernel passed"
    end if

  end subroutine check_block_stencil_kernel


  subroutine accumulate_block_tendency_kernel ( &
       catalog_index,block,scalar_tendency,vector_tendency,context)
    ! Write neighbour-minus-centre sol differences into persistent interior
    ! scalar/vector tendency arrays using compact patch, boundary and ghost
    ! stencil addresses.

    implicit none

    integer, intent(in) :: catalog_index
    type(Block_Data), intent(in) :: block
    real(dp), intent(inout) :: scalar_tendency(:)
    real(dp), intent(inout) :: vector_tendency(:)
    class(*), intent(inout) :: context

    integer :: address_offset
    integer :: center_base
    integer :: center_index
    integer :: center_node
    integer :: component_slot
    integer :: field_base
    integer :: field_index
    integer :: level_slot
    integer :: n_storage_node
    integer :: node_index
    integer :: p
    integer :: q
    integer :: record
    integer :: scalar_slot
    integer :: side
    integer :: storage_class
    integer :: storage_start

    real(dp) :: center_value
    real(dp) :: difference
    real(dp) :: value

    if (catalog_index < 1) then
       call fail("tendency kernel received invalid catalogue index")
    end if
    if (size(scalar_tendency) /= size(block%scalar) .or. &
         size(vector_tendency) /= size(block%vector)) then
       call fail("tendency kernel received invalid output extents")
    end if
    if (block%scalar_mult /= 1 .or. block%vector_mult < 1) then
       call fail("tendency kernel received invalid field multipliers")
    end if
    if (size(block%stencil,1) /= N_BDRY .or. &
         size(block%stencil,2) /= size(block%patch)) then
       call fail("tendency kernel received invalid topology")
    end if

    select type (statistics => context)
    type is (Block_Stencil_Kernel_Context)
       statistics%block_count = statistics%block_count + 1_int64

       do p = 1,size(block%patch)
          do side = 1,N_BDRY
             storage_class = block%stencil(side,p)%storage
             record = block%stencil(side,p)%id
             address_offset = block%stencil(side,p)%offset
             storage_start = 0
             n_storage_node = 0

             select case (storage_class)
             case (STORE_PATCH)
                if (record < 0 .or. record >= size(block%patch)) then
                   call fail("tendency kernel received invalid patch address")
                end if
                storage_start = block%patch(record+1)%elts_start
                n_storage_node = PATCH_SIZE**2
             case (STORE_BDRY)
                if (record < 1 .or. &
                     record > size(block%bdry_storage)) then
                   call fail( &
                        "tendency kernel received invalid boundary address")
                end if
                storage_start = block%bdry_storage(record)%local_start
                n_storage_node = block%bdry_storage(record)%n_node
             case (STORE_GHOST)
                if (record < 1 .or. &
                     record > size(block%ghost_storage)) then
                   call fail( &
                        "tendency kernel received invalid ghost address")
                end if
                storage_start = block%ghost_storage(record)%local_start
                n_storage_node = block%ghost_storage(record)%n_node
             case default
                call fail("tendency kernel received invalid storage class")
             end select

             do q = 0,PATCH_SIZE**2-1
                if (address_offset+q < 0 .or. &
                     address_offset+q >= n_storage_node) cycle

                node_index = storage_start + address_offset + q
                center_node = block%patch(p)%elts_start + q

                if (center_node < 0 .or. &
                     center_node >= size(block%node)) then
                   call fail("tendency kernel received invalid centre node")
                end if
                select case (storage_class)
                case (STORE_PATCH)
                   if (node_index < 0 .or. &
                        node_index >= size(block%node)) then
                      call fail("tendency kernel patch node is invalid")
                   end if
                case (STORE_BDRY)
                   if (node_index < 0 .or. &
                        node_index >= size(block%bdry_node)) then
                      call fail("tendency kernel boundary node is invalid")
                   end if
                case (STORE_GHOST)
                   if (node_index < 0 .or. &
                        node_index >= size(block%ghost_node)) then
                      call fail("tendency kernel ghost node is invalid")
                   end if
                end select

                statistics%address_count(storage_class) = &
                     statistics%address_count(storage_class) + 1_int64

                do scalar_slot = 1,block%n_scalar_variable
                   do level_slot = 1,block%n_field_level
                      center_base = &
                           ((scalar_slot-1)*block%n_field_level + &
                           level_slot-1)*size(block%node)
                      center_index = center_base + center_node + 1
                      center_value = block%scalar(center_index)
                      value = 0.0_dp

                      select case (storage_class)
                      case (STORE_PATCH)
                         field_base = center_base
                         field_index = field_base + node_index + 1
                         value = block%scalar(field_index)
                      case (STORE_BDRY)
                         field_base = &
                              ((scalar_slot-1)*block%n_field_level + &
                              level_slot-1)*size(block%bdry_node)
                         field_index = field_base + node_index + 1
                         value = block%bdry_scalar(field_index)
                      case (STORE_GHOST)
                         field_base = &
                              ((scalar_slot-1)*block%n_field_level + &
                              level_slot-1)*size(block%ghost_node)
                         field_index = field_base + node_index + 1
                         value = block%ghost_scalar(field_index)
                      end select

                      difference = value-center_value
                      scalar_tendency(center_index) = &
                           scalar_tendency(center_index) + difference
                      statistics%scalar_count = &
                           statistics%scalar_count + 1_int64
                      statistics%scalar_difference_moment = &
                           statistics%scalar_difference_moment + &
                           [difference,abs(difference),difference**2]
                   end do
                end do

                do level_slot = 1,block%n_field_level
                   do component_slot = 1,block%vector_mult
                      center_base = (level_slot-1)* &
                           block%vector_mult*size(block%node)
                      center_index = center_base + &
                           block%vector_mult*center_node + component_slot
                      center_value = block%vector(center_index)
                      value = 0.0_dp

                      select case (storage_class)
                      case (STORE_PATCH)
                         field_base = center_base
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%vector(field_index)
                      case (STORE_BDRY)
                         field_base = (level_slot-1)* &
                              block%vector_mult*size(block%bdry_node)
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%bdry_vector(field_index)
                      case (STORE_GHOST)
                         field_base = (level_slot-1)* &
                              block%vector_mult*size(block%ghost_node)
                         field_index = field_base + &
                              block%vector_mult*node_index + &
                              component_slot
                         value = block%ghost_vector(field_index)
                      end select

                      difference = value-center_value
                      vector_tendency(center_index) = &
                           vector_tendency(center_index) + difference
                      statistics%vector_count = &
                           statistics%vector_count + 1_int64
                      statistics%vector_difference_moment = &
                           statistics%vector_difference_moment + &
                           [difference,abs(difference),difference**2]
                   end do
                end do
             end do
          end do
       end do
    class default
       call fail("tendency kernel received invalid context")
    end select

  end subroutine accumulate_block_tendency_kernel


  subroutine check_block_tendency_kernel (verbose)
    ! Validate reusable writable tendency storage and a complete stencil
    ! kernel while leaving authoritative prognostic fields unchanged.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: allocation_after_first
    integer(int64) :: allocation_before
    integer(int64) :: count_global(2)
    integer(int64) :: count_local(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: refresh_count_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_second(2)

    real(dp) :: factor
    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: scale
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_second(3,2)

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: reference
    type(Block_Stencil_Kernel_Context) :: writable_first
    type(Block_Stencil_Kernel_Context) :: writable_second

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))

    refresh_count_before = local_block_hydrostatic_refresh_count()
    allocation_before = local_block_tendency_allocation_count()
    execution_before = local_block_tendency_execution_count()

    call refresh_block_sol_ghosts
    reference = Block_Stencil_Kernel_Context()
    call apply_local_block_field_consumer( &
         accumulate_block_stencil_kernel,reference)

    writable_first = Block_Stencil_Kernel_Context()
    call apply_local_block_tendency_kernel( &
         accumulate_block_tendency_kernel,writable_first)

    if (.not. local_block_tendency_state_ready()) then
       call fail("production tendency output state is not ready")
    end if
    if (local_block_tendency_execution_count() /= execution_before+1_int64) then
       call fail("production tendency execution count mismatch")
    end if

    allocation_after_first = local_block_tendency_allocation_count()
    if (allocation_after_first < allocation_before) then
       call fail("production tendency allocation count regressed")
    end if

    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (any(tendency_count /= field_count)) then
       call fail("production tendency output coverage mismatch")
    end if
    if (writable_first%block_count /= int(n_local_blocks(),int64) .or. &
         any(writable_first%address_count /= reference%address_count) .or. &
         writable_first%scalar_count /= reference%scalar_count .or. &
         writable_first%vector_count /= reference%vector_count) then
       call fail("production tendency stencil traversal mismatch")
    end if
    if (.not. field_moments_match( &
         writable_first%scalar_difference_moment, &
         reference%scalar_difference_moment,reference%scalar_count) .or. &
         .not. field_moments_match( &
         writable_first%vector_difference_moment, &
         reference%vector_difference_moment,reference%vector_count)) then
       call fail("production tendency difference inventory mismatch")
    end if

    factor = 256.0_dp*epsilon(1.0_dp)* &
         real(max(1_int64,reference%scalar_count),dp)
    scale = max(1.0_dp,reference%scalar_difference_moment(2), &
         tendency_moment(2,1))
    if (abs(tendency_moment(1,1)- &
         reference%scalar_difference_moment(1)) > factor*scale) then
       call fail("scalar tendency accumulation is not conservative")
    end if

    factor = 256.0_dp*epsilon(1.0_dp)* &
         real(max(1_int64,reference%vector_count),dp)
    scale = max(1.0_dp,reference%vector_difference_moment(2), &
         tendency_moment(2,2))
    if (abs(tendency_moment(1,2)- &
         reference%vector_difference_moment(1)) > factor*scale) then
       call fail("vector tendency accumulation is not conservative")
    end if

    call refresh_block_sol_ghosts
    writable_second = Block_Stencil_Kernel_Context()
    call apply_local_block_tendency_kernel( &
         accumulate_block_tendency_kernel,writable_second)
    call local_block_tendency_statistics( &
         tendency_count_second(1),tendency_count_second(2), &
         tendency_moment_second(:,1),tendency_moment_second(:,2))

    if (local_block_tendency_allocation_count() /= &
         allocation_after_first) then
       call fail("production tendency workspace was reallocated")
    end if
    if (local_block_tendency_execution_count() /= execution_before+2_int64) then
       call fail("repeated tendency execution count mismatch")
    end if
    if (any(tendency_count_second /= tendency_count)) then
       call fail("repeated tendency output coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_second(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_second(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("repeated tendency output changed")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))

    if (any(field_count_after /= field_count)) then
       call fail("tendency kernel changed prognostic field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("tendency kernel changed prognostic fields")
    end if
    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_before) then
       call fail("tendency kernel refreshed hydrostatic cache")
    end if

    count_local = tendency_count
    call MPI_Allreduce( &
         count_local,count_global,2,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce production tendency counts")

    if (any(count_global <= 0_int64)) then
       call fail("production tendency global output is incomplete")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Writable block tendency kernel for rank ",rank,":"
       write(6,'(a,i0)') &
            "  scalar tendency values = ",tendency_count(1)
       write(6,'(a,i0)') &
            "  vector tendency values = ",tendency_count(2)
       write(6,'(a)') &
            "  complete writable stencil output passed"
       write(6,'(a)') &
            "  conservative difference accumulation passed"
       write(6,'(a)') &
            "  persistent workspace reuse passed"
       write(6,'(a,/)') &
            "  prognostic fields and hydrostatic cache unchanged"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global writable scalar tendency values = ",count_global(1)
       write(6,'(a,i0)') &
            "Global writable vector tendency values = ",count_global(2)
       write(6,'(a,/)') &
            "Persistent writable block stencil tendencies passed"
    end if

  end subroutine check_block_tendency_kernel


  subroutine accumulate_block_tendency_consumer ( &
       catalog_index,scalar_tendency,vector_tendency,context)
    ! Read the persistent tendency arrays supplied directly by the
    ! read-only production traversal.

    implicit none

    integer, intent(in) :: catalog_index
    real(dp), intent(in) :: scalar_tendency(:)
    real(dp), intent(in) :: vector_tendency(:)
    class(*), intent(inout) :: context

    if (catalog_index < 1) then
       call fail("tendency consumer received invalid catalogue index")
    end if

    select type (statistics => context)
    type is (Block_Tendency_Traversal_Context)
       statistics%block_count = statistics%block_count + 1_int64
       statistics%scalar_count = statistics%scalar_count + &
            int(size(scalar_tendency),int64)
       statistics%vector_count = statistics%vector_count + &
            int(size(vector_tendency),int64)
       if (maxval(abs(scalar_tendency)) > 0.0_dp) then
          statistics%scalar_changed_block_count = &
               statistics%scalar_changed_block_count + 1_int64
       end if

       statistics%scalar_moment(1) = &
            statistics%scalar_moment(1) + sum(scalar_tendency)
       statistics%scalar_moment(2) = &
            statistics%scalar_moment(2) + sum(abs(scalar_tendency))
       statistics%scalar_moment(3) = &
            statistics%scalar_moment(3) + sum(scalar_tendency**2)
       statistics%vector_moment(1) = &
            statistics%vector_moment(1) + sum(vector_tendency)
       statistics%vector_moment(2) = &
            statistics%vector_moment(2) + sum(abs(vector_tendency))
       statistics%vector_moment(3) = &
            statistics%vector_moment(3) + sum(vector_tendency**2)
    class default
       call fail("tendency consumer received invalid context")
    end select

  end subroutine accumulate_block_tendency_consumer


  subroutine check_block_tendency_trial_update (verbose)
    ! Validate direct read-only tendency traversal and a reversible shadow
    ! update without changing the active Domain timestep.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: changed_block_global
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: refresh_after
    integer(int64) :: refresh_before
    integer(int64) :: tendency_count(2)

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: trial_scale

    logical :: print_summary

    type(Block_Tendency_Traversal_Context) :: traversal

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("tendency trial requested before output is ready")
    end if
    if (local_block_tendency_trial_is_active()) then
       call fail("tendency trial unexpectedly active")
    end if

    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    traversal = Block_Tendency_Traversal_Context()
    call apply_local_block_tendency_consumer( &
         accumulate_block_tendency_consumer,traversal)

    if (traversal%block_count /= int(n_local_blocks(),int64) .or. &
         traversal%scalar_count /= tendency_count(1) .or. &
         traversal%vector_count /= tendency_count(2)) then
       call fail("read-only tendency traversal coverage mismatch")
    end if
    if (.not. field_moments_match( &
         traversal%scalar_moment,tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         traversal%vector_moment,tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("read-only tendency traversal moment mismatch")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    refresh_before = local_block_hydrostatic_refresh_count()

    trial_scale = epsilon(1.0_dp)**0.25_dp
    call begin_local_block_tendency_trial(trial_scale)
    if (.not. local_block_tendency_trial_is_active()) then
       call fail("reversible tendency trial did not become active")
    end if

    if (compressible .and. &
         traversal%scalar_changed_block_count > 0_int64) then
       if (local_block_hydrostatic_state_ready()) then
          call fail("scalar tendency trial did not invalidate cache")
       end if
       if (local_block_hydrostatic_refresh_count() /= refresh_before) then
          call fail("scalar tendency trial refreshed cache eagerly")
       end if
    end if

    call rollback_local_block_tendency_trial
    if (local_block_tendency_trial_is_active()) then
       call fail("reversible tendency trial remained active")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("tendency rollback changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("tendency rollback did not recover prognostic fields")
    end if

    if (.not. local_block_tendency_state_ready()) then
       call fail("tendency output was lost during trial rollback")
    end if

    if (compressible) then
       call ensure_local_block_hydrostatic_state
       refresh_after = local_block_hydrostatic_refresh_count()
       if (refresh_after-refresh_before /= &
            traversal%scalar_changed_block_count) then
          call fail("selective hydrostatic refresh count mismatch")
       end if
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("hydrostatic cache not ready after rollback")
       end if
    end if

    call MPI_Allreduce( &
         traversal%scalar_changed_block_count,changed_block_global, &
         1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce tendency trial blocks")

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Reversible block tendency trial for rank ",rank,":"
       write(6,'(a,i0)') &
            "  directly traversed scalar values = ", &
            traversal%scalar_count
       write(6,'(a,i0)') &
            "  directly traversed vector values = ", &
            traversal%vector_count
       write(6,'(a,i0)') &
            "  scalar-modified blocks            = ", &
            traversal%scalar_changed_block_count
       write(6,'(a)') "  direct persistent tendency traversal passed"
       write(6,'(a)') "  reversible scalar/vector trial update passed"
       write(6,'(a)') "  selective hydrostatic invalidation passed"
       write(6,'(a)') "  exact prognostic rollback passed"
       write(6,'(a,/)') "  persistent tendency outputs retained"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global scalar-modified trial blocks = ", &
            changed_block_global
       write(6,'(a,/)') &
            "Reversible block tendency shadow update passed"
    end if

  end subroutine check_block_tendency_trial_update


  subroutine check_block_tendency_commit (verbose)
    ! Validate the commit transition without advancing the physical state.
    ! A zero increment keeps the Domain/block shadow comparison exact while
    ! exercising stale-output invalidation and allocation reuse.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: allocation_before
    integer(int64) :: count_global(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: refresh_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: regenerated

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("tendency commit requested before output is ready")
    end if
    if (local_block_tendency_trial_is_active()) then
       call fail("tendency commit found an active trial")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    refresh_before = local_block_hydrostatic_refresh_count()
    allocation_before = local_block_tendency_allocation_count()
    execution_before = local_block_tendency_execution_count()

    call begin_local_block_tendency_trial(0.0_dp)
    if (.not. local_block_tendency_trial_is_active()) then
       call fail("zero-increment tendency trial did not become active")
    end if

    call commit_local_block_tendency_trial
    if (local_block_tendency_trial_is_active()) then
       call fail("committed tendency trial remained active")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("commit retained stale tendency outputs")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("zero-increment commit changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("zero-increment commit changed prognostic fields")
    end if

    if (compressible) then
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("zero-increment commit invalidated hydrostatic cache")
       end if
       if (local_block_hydrostatic_refresh_count() /= refresh_before) then
          call fail("zero-increment commit refreshed hydrostatic cache")
       end if
    end if

    if (.not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("zero-increment commit checkpoint is not ready")
    end if
    call finalize_local_block_tendency_commit
    if (local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("zero-increment commit checkpoint was not finalized")
    end if

    regenerated = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,regenerated)

    if (.not. local_block_tendency_state_ready()) then
       call fail("committed tendency output was not regenerated")
    end if
    if (local_block_tendency_allocation_count() /= allocation_before) then
       call fail("committed tendency regeneration reallocated workspace")
    end if
    if (local_block_tendency_execution_count() /= execution_before+1_int64) then
       call fail("committed tendency regeneration count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("regenerated tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("regenerated tendency output changed")
    end if

    call MPI_Allreduce( &
         tendency_count_after,count_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce committed tendency counts")
    if (any(count_global <= 0_int64)) then
       call fail("committed tendency regeneration is incomplete")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Committed block tendency lifecycle for rank ",rank,":"
       write(6,'(a)') "  zero-increment commit preserved fields and cache"
       write(6,'(a)') "  pre-commit tendency outputs invalidated"
       write(6,'(a)') "  persistent tendency allocation reuse passed"
       write(6,'(a,/)') "  regenerated tendency inventory passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global regenerated scalar tendencies = ",count_global(1)
       write(6,'(a,i0)') &
            "Global regenerated vector tendencies = ",count_global(2)
       write(6,'(a,/)') &
            "Committed block tendency lifecycle passed"
    end if

  end subroutine check_block_tendency_commit


  subroutine check_block_tendency_step_driver (verbose)
    ! Exercise two complete refresh, tendency, nonzero trial, rollback and
    ! hydrostatic-recovery cycles through the guarded production ordering.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: cycle
    integer :: ierr

    integer(int64) :: allocation_before
    integer(int64) :: changed_block_global
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: refresh_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: trial_scale

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: kernel_context
    type(Block_Tendency_Traversal_Context) :: traversal

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("block step driver requested before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active()) then
       call fail("block step driver found an active trial")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    traversal = Block_Tendency_Traversal_Context()
    call apply_local_block_tendency_consumer( &
         accumulate_block_tendency_consumer,traversal)

    if (compressible) call ensure_local_block_hydrostatic_state
    refresh_before = local_block_hydrostatic_refresh_count()
    allocation_before = local_block_tendency_allocation_count()
    execution_before = local_block_tendency_execution_count()
    trial_scale = epsilon(1.0_dp)**0.25_dp

    do cycle = 1,2
       kernel_context = Block_Stencil_Kernel_Context()
       call apply_refreshed_block_tendency_kernel( &
            accumulate_block_tendency_kernel,kernel_context)

       if (kernel_context%block_count /= int(n_local_blocks(),int64)) then
          call fail("guarded block step kernel coverage mismatch")
       end if

       call begin_local_block_tendency_trial(trial_scale)
       if (.not. local_block_tendency_trial_is_active()) then
          call fail("guarded block step trial did not become active")
       end if

       if (compressible .and. &
            traversal%scalar_changed_block_count > 0_int64) then
          if (local_block_hydrostatic_state_ready()) then
             call fail("guarded block step retained stale hydrostatic cache")
          end if
       end if

       call rollback_local_block_tendency_trial
       if (local_block_tendency_trial_is_active()) then
          call fail("guarded block step rollback remained active")
       end if

       if (compressible) then
          call ensure_local_block_hydrostatic_state
          if (local_block_hydrostatic_refresh_count()-refresh_before /= &
               int(cycle,int64)* &
               traversal%scalar_changed_block_count) then
             call fail("guarded block step hydrostatic refresh mismatch")
          end if
       end if

       call local_block_field_statistics( &
            field_count_after(1),field_count_after(2), &
            field_moment_after(:,1),field_moment_after(:,2))
       if (any(field_count_after /= field_count)) then
          call fail("guarded block step changed field coverage")
       end if
       if (.not. field_moments_match( &
            field_moment_after(:,1),field_moment(:,1), &
            field_count(1)) .or. &
            .not. field_moments_match( &
            field_moment_after(:,2),field_moment(:,2), &
            field_count(2))) then
          call fail("guarded block step rollback changed fields")
       end if
    end do

    if (local_block_tendency_allocation_count() /= allocation_before) then
       call fail("guarded block step reallocated tendency workspace")
    end if
    if (local_block_tendency_execution_count() /= execution_before+2_int64) then
       call fail("guarded block step execution count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("guarded block step tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("guarded block step tendency inventory changed")
    end if

    call MPI_Allreduce( &
         traversal%scalar_changed_block_count,changed_block_global, &
         1,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce guarded block step blocks")

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Guarded block tendency step driver for rank ",rank,":"
       write(6,'(a,i0)') &
            "  complete shadow cycles = ",2
       write(6,'(a,i0)') &
            "  scalar-modified blocks = ", &
            traversal%scalar_changed_block_count
       write(6,'(a)') "  sol ghost refresh before each kernel passed"
       write(6,'(a)') "  nonzero update and exact rollback passed"
       write(6,'(a)') "  selective hydrostatic recovery passed"
       write(6,'(a,/)') "  persistent tendency workspace reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global guarded-step scalar-modified blocks = ", &
            changed_block_global
       write(6,'(a,/)') &
            "Guarded block tendency step sequence passed"
    end if

  end subroutine check_block_tendency_step_driver


  subroutine check_block_tendency_accepted_step (verbose)
    ! Commit one nonzero block update, consume its derived state, then use
    ! the retained one-level checkpoint to recover the exact Domain shadow.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: allocation_before
    integer(int64) :: changed_block_count(2)
    integer(int64) :: changed_block_global(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: refresh_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: max_update(2)
    real(dp) :: max_update_global(2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: trial_scale

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: committed_kernel
    type(Block_Stencil_Kernel_Context) :: restored_kernel
    type(Block_Tendency_Traversal_Context) :: traversal

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("accepted block step requested before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("accepted block step found pending transaction state")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    traversal = Block_Tendency_Traversal_Context()
    call apply_local_block_tendency_consumer( &
         accumulate_block_tendency_consumer,traversal)

    if (compressible) call ensure_local_block_hydrostatic_state
    refresh_before = local_block_hydrostatic_refresh_count()
    allocation_before = local_block_tendency_allocation_count()
    execution_before = local_block_tendency_execution_count()
    trial_scale = epsilon(1.0_dp)**0.25_dp

    call begin_local_block_tendency_trial(trial_scale)
    call commit_local_block_tendency_trial

    if (local_block_tendency_trial_is_active() .or. &
         .not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("nonzero block update was not committed")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("accepted block update retained stale tendencies")
    end if

    call local_block_tendency_commit_checkpoint_statistics( &
         changed_block_count(1),changed_block_count(2), &
         max_update(1),max_update(2))
    if (changed_block_count(1) /= &
         traversal%scalar_changed_block_count) then
       call fail("accepted scalar update coverage mismatch")
    end if

    committed_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,committed_kernel)
    if (committed_kernel%block_count /= int(n_local_blocks(),int64)) then
       call fail("accepted-state tendency traversal incomplete")
    end if

    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (local_block_hydrostatic_refresh_count()-refresh_before /= &
            changed_block_count(1)) then
          call fail("accepted-state hydrostatic refresh mismatch")
       end if
    end if

    call restore_local_block_tendency_commit
    if (local_block_tendency_commit_checkpoint_is_ready() .or. &
         local_block_tendency_state_ready()) then
       call fail("accepted block checkpoint restore left stale state ready")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("accepted block checkpoint changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("accepted block checkpoint did not restore fields")
    end if

    restored_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,restored_kernel)
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (local_block_hydrostatic_refresh_count()-refresh_before /= &
            2_int64*changed_block_count(1)) then
          call fail("restored-state hydrostatic refresh mismatch")
       end if
    end if

    if (local_block_tendency_allocation_count() /= allocation_before) then
       call fail("accepted block step reallocated tendency workspace")
    end if
    if (local_block_tendency_execution_count() /= execution_before+2_int64) then
       call fail("accepted block step execution count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("restored tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("restored tendency inventory changed")
    end if

    call MPI_Allreduce( &
         changed_block_count,changed_block_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce accepted block update counts")
    call MPI_Allreduce( &
         max_update,max_update_global,2, &
         MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce accepted block update maxima")
    if (changed_block_global(1) <= 0_int64 .or. &
         maxval(max_update_global) <= 0.0_dp) then
       call fail("accepted block update produced no change")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Accepted nonzero block step for rank ",rank,":"
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector modified blocks = ",changed_block_count
       write(6,'(a,2(es14.6,1x))') &
            "  maximum scalar/vector updates = ",max_update
       write(6,'(a)') "  committed-state ghost refresh passed"
       write(6,'(a)') "  committed-state tendency regeneration passed"
       write(6,'(a)') "  committed-state hydrostatic rebuild passed"
       write(6,'(a,/)') "  exact checkpoint recovery passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global accepted scalar/vector modified blocks = ", &
            changed_block_global
       write(6,'(a,2(es14.6,1x))') &
            "Global maximum scalar/vector updates = ", &
            max_update_global
       write(6,'(a,/)') &
            "Accepted nonzero block tendency step passed"
    end if

  end subroutine check_block_tendency_accepted_step


  subroutine check_block_multistage_tendency_accumulator (verbose)
    ! Combine tendencies from the original and one accepted intermediate
    ! state, apply the weighted register reversibly, and recover the exact
    ! Domain-shadow fields and original tendency inventory.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: accumulator_allocation_after_reset
    integer(int64) :: accumulator_allocation_before
    integer(int64) :: accumulator_changed_count(2)
    integer(int64) :: accumulator_changed_global(2)
    integer(int64) :: accumulator_count(2)
    integer(int64) :: accumulator_stage_count
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: tendency_allocation_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)

    real(dp) :: accumulator_moment(3,2)
    real(dp) :: accumulator_abs_local(2)
    real(dp) :: accumulator_abs_global(2)
    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: trial_scale

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: intermediate_kernel
    type(Block_Stencil_Kernel_Context) :: restored_kernel

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("multi-stage accumulator before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("multi-stage accumulator found pending transaction")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    tendency_allocation_before = local_block_tendency_allocation_count()
    accumulator_allocation_before = &
         local_block_tendency_accumulator_allocation_count()
    execution_before = local_block_tendency_execution_count()
    trial_scale = epsilon(1.0_dp)**0.25_dp

    call reset_local_block_tendency_accumulator
    if (.not. local_block_tendency_accumulator_state_ready()) then
       call fail("multi-stage accumulator did not become ready")
    end if
    accumulator_allocation_after_reset = &
         local_block_tendency_accumulator_allocation_count()
    if (accumulator_allocation_after_reset < &
         accumulator_allocation_before) then
       call fail("multi-stage accumulator allocation count regressed")
    end if

    call accumulate_local_block_tendency(0.5_dp)

    call begin_local_block_tendency_trial(trial_scale)
    call commit_local_block_tendency_trial
    intermediate_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,intermediate_kernel)
    if (intermediate_kernel%block_count /= &
         int(n_local_blocks(),int64)) then
       call fail("multi-stage intermediate tendency traversal incomplete")
    end if

    call accumulate_local_block_tendency(0.5_dp)
    call local_block_tendency_accumulator_statistics( &
         accumulator_count(1),accumulator_count(2), &
         accumulator_changed_count(1),accumulator_changed_count(2), &
         accumulator_stage_count, &
         accumulator_moment(:,1),accumulator_moment(:,2))

    if (any(accumulator_count /= field_count)) then
       call fail("multi-stage accumulator field coverage mismatch")
    end if
    if (accumulator_stage_count /= 2_int64) then
       call fail("multi-stage accumulator stage count mismatch")
    end if

    if (compressible) call ensure_local_block_hydrostatic_state
    call restore_local_block_tendency_commit

    call begin_local_block_accumulated_tendency_trial(trial_scale)
    if (.not. local_block_tendency_trial_is_active()) then
       call fail("multi-stage accumulated trial did not become active")
    end if
    if (compressible .and. &
         accumulator_changed_count(1) > 0_int64) then
       if (local_block_hydrostatic_state_ready()) then
          call fail("multi-stage accumulated trial retained stale cache")
       end if
    end if

    call rollback_local_block_tendency_trial
    if (compressible) call ensure_local_block_hydrostatic_state

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("multi-stage recovery changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("multi-stage recovery did not restore fields")
    end if

    restored_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,restored_kernel)

    if (local_block_tendency_allocation_count() /= &
         tendency_allocation_before) then
       call fail("multi-stage cycle reallocated tendency workspace")
    end if
    if (local_block_tendency_execution_count() /= &
         execution_before+2_int64) then
       call fail("multi-stage tendency execution count mismatch")
    end if
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_after_reset) then
       call fail("multi-stage cycle reallocated accumulator storage")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("multi-stage restored tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("multi-stage restored tendency inventory changed")
    end if

    call MPI_Allreduce( &
         accumulator_changed_count,accumulator_changed_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce multi-stage accumulator blocks")
    accumulator_abs_local = accumulator_moment(2,:)
    call MPI_Allreduce( &
         accumulator_abs_local,accumulator_abs_global,2, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce multi-stage accumulator moments")
    if (any(accumulator_changed_global <= 0_int64) .or. &
         maxval(accumulator_abs_global) <= 0.0_dp) then
       call fail("multi-stage accumulator produced no change")
    end if

    call reset_local_block_tendency_accumulator
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_after_reset) then
       call fail("multi-stage accumulator reset reallocated storage")
    end if

    if (compressible) then
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("multi-stage hydrostatic recovery incomplete")
       end if
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Persistent multi-stage tendency register for rank ",rank,":"
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector register values = ",accumulator_count
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector modified blocks = ", &
            accumulator_changed_count
       write(6,'(a)') "  two weighted tendency stages accumulated"
       write(6,'(a)') "  accepted intermediate-state evaluation passed"
       write(6,'(a)') "  accumulated nonzero update and rollback passed"
       write(6,'(a,/)') "  tendency and accumulator workspace reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global multi-stage scalar/vector modified blocks = ", &
            accumulator_changed_global
       write(6,'(a,/)') &
            "Persistent multi-stage block tendency accumulation passed"
    end if

  end subroutine check_block_multistage_tendency_accumulator


  subroutine check_block_multistage_tendency_commit (verbose)
    ! Commit a weighted two-stage tendency, consume all derived state from
    ! the accepted result, then recover the exact retained Domain shadow.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: accumulator_allocation_before
    integer(int64) :: accumulator_changed_count(2)
    integer(int64) :: accumulator_changed_global(2)
    integer(int64) :: accumulator_count(2)
    integer(int64) :: accumulator_stage_count
    integer(int64) :: changed_block_count(2)
    integer(int64) :: changed_block_global(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: tendency_allocation_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)

    real(dp) :: accumulator_moment(3,2)
    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: max_update(2)
    real(dp) :: max_update_global(2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: trial_scale

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: accepted_kernel
    type(Block_Stencil_Kernel_Context) :: intermediate_kernel
    type(Block_Stencil_Kernel_Context) :: restored_kernel

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("multi-stage commit before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("multi-stage commit found pending transaction")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    tendency_allocation_before = local_block_tendency_allocation_count()
    accumulator_allocation_before = &
         local_block_tendency_accumulator_allocation_count()
    execution_before = local_block_tendency_execution_count()
    trial_scale = epsilon(1.0_dp)**0.25_dp

    call reset_local_block_tendency_accumulator
    if (.not. local_block_tendency_accumulator_state_ready()) then
       call fail("multi-stage commit accumulator did not become ready")
    end if
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("multi-stage commit reset reallocated accumulator")
    end if

    call accumulate_local_block_tendency(0.5_dp)

    call begin_local_block_tendency_trial(trial_scale)
    call commit_local_block_tendency_trial
    intermediate_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,intermediate_kernel)
    if (intermediate_kernel%block_count /= &
         int(n_local_blocks(),int64)) then
       call fail("multi-stage commit intermediate traversal incomplete")
    end if

    call accumulate_local_block_tendency(0.5_dp)
    call local_block_tendency_accumulator_statistics( &
         accumulator_count(1),accumulator_count(2), &
         accumulator_changed_count(1),accumulator_changed_count(2), &
         accumulator_stage_count, &
         accumulator_moment(:,1),accumulator_moment(:,2))
    if (any(accumulator_count /= field_count)) then
       call fail("multi-stage commit accumulator coverage mismatch")
    end if
    if (accumulator_stage_count /= 2_int64) then
       call fail("multi-stage commit accumulator stage mismatch")
    end if

    if (compressible) call ensure_local_block_hydrostatic_state
    call restore_local_block_tendency_commit
    if (local_block_tendency_commit_checkpoint_is_ready() .or. &
         local_block_tendency_state_ready()) then
       call fail("multi-stage intermediate restore left stale state ready")
    end if

    call begin_local_block_accumulated_tendency_trial(trial_scale)
    call commit_local_block_tendency_trial
    if (local_block_tendency_trial_is_active() .or. &
         .not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("multi-stage accumulated update was not committed")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("multi-stage committed update retained stale tendencies")
    end if

    call local_block_tendency_commit_checkpoint_statistics( &
         changed_block_count(1),changed_block_count(2), &
         max_update(1),max_update(2))
    if (any(changed_block_count /= accumulator_changed_count)) then
       call fail("multi-stage committed update coverage mismatch")
    end if

    accepted_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,accepted_kernel)
    if (accepted_kernel%block_count /= int(n_local_blocks(),int64)) then
       call fail("multi-stage accepted-state traversal incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("multi-stage accepted hydrostatic rebuild incomplete")
       end if
    end if

    call restore_local_block_tendency_commit
    if (local_block_tendency_commit_checkpoint_is_ready() .or. &
         local_block_tendency_state_ready()) then
       call fail("multi-stage checkpoint restore left stale state ready")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("multi-stage checkpoint changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("multi-stage checkpoint did not restore fields")
    end if

    restored_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,restored_kernel)
    if (restored_kernel%block_count /= int(n_local_blocks(),int64)) then
       call fail("multi-stage restored-state traversal incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("multi-stage restored hydrostatic rebuild incomplete")
       end if
    end if

    if (local_block_tendency_allocation_count() /= &
         tendency_allocation_before) then
       call fail("multi-stage commit reallocated tendency workspace")
    end if
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("multi-stage commit reallocated accumulator storage")
    end if
    if (local_block_tendency_execution_count() /= &
         execution_before+3_int64) then
       call fail("multi-stage commit execution count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("multi-stage commit restored tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("multi-stage commit restored tendency inventory changed")
    end if

    call MPI_Allreduce( &
         changed_block_count,changed_block_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce multi-stage committed blocks")
    call MPI_Allreduce( &
         accumulator_changed_count,accumulator_changed_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce multi-stage register blocks")
    call MPI_Allreduce( &
         max_update,max_update_global,2, &
         MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce multi-stage committed maxima")
    if (any(changed_block_global /= accumulator_changed_global)) then
       call fail("global multi-stage committed coverage mismatch")
    end if
    if (any(changed_block_global <= 0_int64) .or. &
         any(max_update_global <= 0.0_dp)) then
       call fail("multi-stage committed update produced no change")
    end if

    call reset_local_block_tendency_accumulator
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("multi-stage post-commit reset reallocated storage")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Accepted multi-stage block step for rank ",rank,":"
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector committed blocks = ",changed_block_count
       write(6,'(a,2(es14.6,1x))') &
            "  maximum scalar/vector updates = ",max_update
       write(6,'(a)') "  weighted two-stage update committed"
       write(6,'(a)') "  accepted-state ghost refresh passed"
       write(6,'(a)') "  accepted-state tendency regeneration passed"
       write(6,'(a)') "  accepted-state hydrostatic rebuild passed"
       write(6,'(a)') "  exact multi-stage checkpoint recovery passed"
       write(6,'(a,/)') &
            "  tendency and accumulator workspace reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global committed multi-stage scalar/vector blocks = ", &
            changed_block_global
       write(6,'(a,2(es14.6,1x))') &
            "Global committed multi-stage maximum updates = ", &
            max_update_global
       write(6,'(a,/)') &
            "Accepted multi-stage block tendency step passed"
    end if

  end subroutine check_block_multistage_tendency_commit


  subroutine check_block_two_stage_step_driver (verbose)
    ! Validate the reusable production-facing two-stage transaction through
    ! accepted-state consumers and exact checkpoint recovery.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: accumulator_allocation_before
    integer(int64) :: changed_block_count(2)
    integer(int64) :: changed_block_global(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: tendency_allocation_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)
    integer(int64) :: writeback_before

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: max_update(2)
    real(dp) :: max_update_global(2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: trial_scale
    real(dp) :: weight(2)

    logical :: print_summary

    type(Block_Stencil_Kernel_Context) :: accepted_kernel
    type(Block_Stencil_Kernel_Context) :: driver_kernel
    type(Block_Stencil_Kernel_Context) :: restored_kernel
    type(Block_Two_Stage_Step_Result) :: result

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("two-stage driver check before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage driver check found pending transaction")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    tendency_allocation_before = local_block_tendency_allocation_count()
    accumulator_allocation_before = &
         local_block_tendency_accumulator_allocation_count()
    execution_before = local_block_tendency_execution_count()
    writeback_before = block_domain_production_writeback_count()
    trial_scale = epsilon(1.0_dp)**0.25_dp
    weight = 0.5_dp

    driver_kernel = Block_Stencil_Kernel_Context()
    call begin_block_two_stage_tendency_step( &
         accumulate_block_tendency_kernel,driver_kernel, &
         trial_scale,weight,result)

    if (driver_kernel%block_count /= &
         2_int64*int(n_local_blocks(),int64)) then
       call fail("two-stage production driver traversal incomplete")
    end if
    if (result%scalar_count /= field_count(1) .or. &
         result%vector_count /= field_count(2)) then
       call fail("two-stage production driver coverage mismatch")
    end if
    if (result%stage_count /= 2_int64) then
       call fail("two-stage production driver stage count mismatch")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         .not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage production driver checkpoint is not ready")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("two-stage production driver tendencies are not stale")
    end if
    if (compressible .and. &
         result%scalar_changed_block_count > 0_int64) then
       if (local_block_hydrostatic_state_ready()) then
          call fail("two-stage production driver retained stale cache")
       end if
    end if

    accepted_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,accepted_kernel)
    if (accepted_kernel%block_count /= int(n_local_blocks(),int64)) then
       call fail("two-stage driver accepted traversal incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("two-stage driver accepted hydrostatic state not ready")
       end if
    end if

    call complete_block_two_stage_tendency_step(.false.)
    if (block_domain_production_writeback_count() /= writeback_before) then
       call fail("rejected two-stage step performed Domain writeback")
    end if
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
    if (local_block_tendency_commit_checkpoint_is_ready() .or. &
         local_block_tendency_state_ready()) then
       call fail("two-stage driver recovery left stale state ready")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("two-stage driver recovery changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("two-stage driver recovery did not restore fields")
    end if

    restored_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,restored_kernel)
    if (restored_kernel%block_count /= int(n_local_blocks(),int64)) then
       call fail("two-stage driver restored traversal incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("two-stage driver restored hydrostatic state not ready")
       end if
    end if

    if (local_block_tendency_allocation_count() /= &
         tendency_allocation_before) then
       call fail("two-stage driver reallocated tendency workspace")
    end if
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("two-stage driver reallocated accumulator workspace")
    end if
    if (local_block_tendency_execution_count() /= &
         execution_before+4_int64) then
       call fail("two-stage production driver execution count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("two-stage driver restored tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("two-stage driver restored tendency inventory changed")
    end if

    changed_block_count = [result%scalar_changed_block_count, &
         result%vector_changed_block_count]
    max_update = [result%scalar_max_update,result%vector_max_update]
    call MPI_Allreduce( &
         changed_block_count,changed_block_global,2, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce two-stage driver blocks")
    call MPI_Allreduce( &
         max_update,max_update_global,2, &
         MPI_DOUBLE_PRECISION,MPI_MAX,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce two-stage driver maxima")
    if (any(changed_block_global <= 0_int64) .or. &
         any(max_update_global <= 0.0_dp)) then
       call fail("two-stage production driver produced no change")
    end if

    call reset_local_block_tendency_accumulator
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("two-stage driver reset reallocated accumulator")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Reusable two-stage block driver for rank ",rank,":"
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector committed blocks = ",changed_block_count
       write(6,'(a,2(es14.6,1x))') &
            "  maximum scalar/vector updates = ",max_update
       write(6,'(a)') "  guarded stage refresh and evaluation passed"
       write(6,'(a)') "  weighted accumulated commit passed"
       write(6,'(a)') "  caller-visible recovery checkpoint passed"
       write(6,'(a)') "  rejected step performed no Domain writeback"
       write(6,'(a)') "  accepted-state derived consumers passed"
       write(6,'(a,/)') "  exact recovery and workspace reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global two-stage driver scalar/vector blocks = ", &
            changed_block_global
       write(6,'(a,2(es14.6,1x))') &
            "Global two-stage driver maximum updates = ", &
            max_update_global
       write(6,'(a,/)') &
            "Reusable production two-stage block step driver passed"
    end if

  end subroutine check_block_two_stage_step_driver


  subroutine check_block_two_stage_step_completion (verbose)
    ! Exercise permanent acceptance through the production completion API
    ! twice with a zero increment. Each accepted step must synchronize sol
    ! exactly while preserving the physical state and wav_coeff shadow.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: accumulator_allocation_before
    integer(int64) :: changed_block_count(2)
    integer(int64) :: count_global(2)
    integer(int64) :: execution_before
    integer(int64) :: field_count(2)
    integer(int64) :: field_count_after(2)
    integer(int64) :: tendency_allocation_before
    integer(int64) :: tendency_count(2)
    integer(int64) :: tendency_count_after(2)
    integer(int64) :: writeback_allocation_before
    integer(int64) :: writeback_before

    real(dp) :: field_moment(3,2)
    real(dp) :: field_moment_after(3,2)
    real(dp) :: max_update(2)
    real(dp) :: tendency_moment(3,2)
    real(dp) :: tendency_moment_after(3,2)
    real(dp) :: weight(2)

    logical :: print_summary
    logical :: scalar_moments_match
    logical :: vector_moments_match

    type(Block_Stencil_Kernel_Context) :: driver_kernel
    type(Block_Stencil_Kernel_Context) :: regenerated_kernel
    type(Block_Stencil_Kernel_Context) :: repeated_kernel
    type(Block_Two_Stage_Step_Result) :: result

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. local_block_tendency_state_ready()) then
       call fail("two-stage completion check before tendency is ready")
    end if
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("two-stage completion check found pending transaction")
    end if

    call local_block_field_statistics( &
         field_count(1),field_count(2), &
         field_moment(:,1),field_moment(:,2))
    call local_block_tendency_statistics( &
         tendency_count(1),tendency_count(2), &
         tendency_moment(:,1),tendency_moment(:,2))

    if (compressible) call ensure_local_block_hydrostatic_state
    tendency_allocation_before = local_block_tendency_allocation_count()
    accumulator_allocation_before = &
         local_block_tendency_accumulator_allocation_count()
    execution_before = local_block_tendency_execution_count()
    writeback_allocation_before = &
         block_writeback_plan_allocation_count()
    writeback_before = block_domain_production_writeback_count()
    weight = 0.5_dp

    driver_kernel = Block_Stencil_Kernel_Context()
    call begin_block_two_stage_tendency_step( &
         accumulate_block_tendency_kernel,driver_kernel, &
         0.0_dp,weight,result)

    if (driver_kernel%block_count /= &
         2_int64*int(n_local_blocks(),int64)) then
       call fail("accepted two-stage completion traversal incomplete")
    end if
    if (result%scalar_count /= field_count(1) .or. &
         result%vector_count /= field_count(2) .or. &
         result%stage_count /= 2_int64) then
       call fail("accepted two-stage completion coverage mismatch")
    end if

    changed_block_count = [result%scalar_changed_block_count, &
         result%vector_changed_block_count]
    max_update = [result%scalar_max_update,result%vector_max_update]
    if (any(changed_block_count /= 0_int64) .or. &
         maxval(abs(max_update)) > 0.0_dp) then
       call fail("zero-increment two-stage completion changed fields")
    end if
    if (.not. local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("accepted two-stage completion checkpoint is not ready")
    end if
    if (compressible) then
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("zero-increment two-stage step invalidated cache")
       end if
    end if

    call complete_block_two_stage_tendency_step(.true.)
    if (block_domain_production_writeback_count() /= &
         writeback_before+1_int64) then
       call fail("accepted step did not perform one Domain writeback")
    end if
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)
    if (local_block_tendency_trial_is_active() .or. &
         local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("accepted two-stage completion remained pending")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("accepted two-stage completion retained stale tendency")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("accepted two-stage completion changed field coverage")
    end if
    if (.not. field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1)) .or. &
         .not. field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))) then
       call fail("accepted zero-increment step changed fields")
    end if

    regenerated_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,regenerated_kernel)
    if (regenerated_kernel%block_count /= &
         int(n_local_blocks(),int64)) then
       call fail("accepted two-stage regeneration incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("accepted two-stage hydrostatic state not ready")
       end if
    end if

    repeated_kernel = Block_Stencil_Kernel_Context()
    call begin_block_two_stage_tendency_step( &
         accumulate_block_tendency_kernel,repeated_kernel, &
         0.0_dp,weight,result)
    if (repeated_kernel%block_count /= &
         2_int64*int(n_local_blocks(),int64)) then
       call fail("repeated accepted completion traversal incomplete")
    end if
    changed_block_count = [result%scalar_changed_block_count, &
         result%vector_changed_block_count]
    max_update = [result%scalar_max_update,result%vector_max_update]
    if (any(changed_block_count /= 0_int64) .or. &
         maxval(abs(max_update)) > 0.0_dp) then
       call fail("repeated zero-increment completion changed fields")
    end if

    call complete_block_two_stage_tendency_step(.true.)
    if (block_domain_production_writeback_count() /= &
         writeback_before+2_int64) then
       call fail("repeated accepted step Domain writeback count mismatch")
    end if
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)
    if (local_block_tendency_trial_is_active()) then
       call fail("repeated accepted completion left active trial")
    end if
    if (local_block_tendency_commit_checkpoint_is_ready()) then
       call fail("repeated accepted completion left checkpoint ready")
    end if
    if (local_block_tendency_state_ready()) then
       call fail("repeated accepted completion left stale tendency")
    end if

    call local_block_field_statistics( &
         field_count_after(1),field_count_after(2), &
         field_moment_after(:,1),field_moment_after(:,2))
    if (any(field_count_after /= field_count)) then
       call fail("repeated accepted completion changed field coverage")
    end if
    scalar_moments_match = field_moments_match( &
         field_moment_after(:,1),field_moment(:,1),field_count(1))
    vector_moments_match = field_moments_match( &
         field_moment_after(:,2),field_moment(:,2),field_count(2))
    if (.not. scalar_moments_match .or. &
         .not. vector_moments_match) then
       call fail("repeated accepted zero-increment step changed fields")
    end if

    regenerated_kernel = Block_Stencil_Kernel_Context()
    call apply_refreshed_block_tendency_kernel( &
         accumulate_block_tendency_kernel,regenerated_kernel)
    if (regenerated_kernel%block_count /= &
         int(n_local_blocks(),int64)) then
       call fail("repeated accepted two-stage regeneration incomplete")
    end if
    if (compressible) then
       call ensure_local_block_hydrostatic_state
       if (.not. local_block_hydrostatic_state_ready()) then
          call fail("repeated accepted hydrostatic state not ready")
       end if
    end if

    if (local_block_tendency_allocation_count() /= &
         tendency_allocation_before) then
       call fail("two-stage completion reallocated tendency workspace")
    end if
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("two-stage completion reallocated accumulator workspace")
    end if
    if (block_writeback_plan_allocation_count() /= &
         writeback_allocation_before) then
       call fail("two-stage completion reallocated writeback buffers")
    end if
    if (local_block_tendency_execution_count() /= &
         execution_before+6_int64) then
       call fail("two-stage completion execution count mismatch")
    end if

    call local_block_tendency_statistics( &
         tendency_count_after(1),tendency_count_after(2), &
         tendency_moment_after(:,1),tendency_moment_after(:,2))
    if (any(tendency_count_after /= tendency_count)) then
       call fail("two-stage completion tendency coverage changed")
    end if
    if (.not. field_moments_match( &
         tendency_moment_after(:,1),tendency_moment(:,1), &
         tendency_count(1)) .or. &
         .not. field_moments_match( &
         tendency_moment_after(:,2),tendency_moment(:,2), &
         tendency_count(2))) then
       call fail("two-stage completion tendency inventory changed")
    end if

    call MPI_Allreduce( &
         field_count,count_global,2,MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce two-stage completion values")
    if (any(count_global <= 0_int64)) then
       call fail("two-stage completion field inventory is empty")
    end if

    call reset_local_block_tendency_accumulator
    if (local_block_tendency_accumulator_allocation_count() /= &
         accumulator_allocation_before) then
       call fail("two-stage completion reset reallocated accumulator")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Completed two-stage block transaction for rank ",rank,":"
       write(6,'(a,2(i0,1x))') &
            "  scalar/vector accepted values = ",field_count
       write(6,'(a)') "  caller completion decision cleared checkpoint"
       write(6,'(a)') "  explicit permanent acceptance path passed"
       write(6,'(a)') "  zero-increment acceptance preserved fields"
       write(6,'(a)') "  accepted sol synchronized to Domain owners"
       write(6,'(a)') "  wav_coeff Domain shadow remained unchanged"
       write(6,'(a)') "  repeated accepted Domain writeback passed"
       write(6,'(a)') "  accepted-state tendency readiness passed"
       write(6,'(a,/)') "  persistent writeback workspace reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,2(i0,1x))') &
            "Global completed two-stage scalar/vector values = ", &
            count_global
       write(6,'(a,/)') &
            "Production two-stage Domain synchronization passed"
    end if

  end subroutine check_block_two_stage_step_completion


  subroutine check_parallel_block_lifecycle (verbose)
    ! Validate checkpoint synchronization and non-destructive teardown and
    ! reconstruction of topology-dependent persistent communication plans.

    implicit none

    logical, optional, intent(in) :: verbose

    integer(int64) :: allocation_after_rebuild
    integer(int64) :: allocation_before
    integer(int64) :: stage_allocation_after_rebuild
    integer(int64) :: stage_allocation_before
    integer(int64) :: writeback_before

    logical :: local_store_ready
    logical :: print_summary
    logical :: state_ready

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("lifecycle check before parallel block state is ready")
    end if

    allocation_before = block_writeback_plan_allocation_count()
    stage_allocation_before = block_writeback_plan%stage_allocations
    writeback_before = block_domain_production_writeback_count()

    call synchronize_parallel_block_checkpoint
    if (block_domain_production_writeback_count() /= &
         writeback_before+2_int64) then
       call fail("checkpoint did not synchronize both field families")
    end if
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("checkpoint synchronization invalidated block state")
    end if
    if (block_writeback_plan_allocation_count() /= allocation_before) then
       call fail("checkpoint synchronization reallocated writeback buffers")
    end if
    if (block_writeback_plan%stage_allocations /= &
         stage_allocation_before) then
       call fail("checkpoint synchronization reallocated Domain staging")
    end if

    call clear_block_ghost_exchange_plan
    call clear_block_writeback_plan

    state_ready = parallel_block_state_is_ready()
    if (state_ready) then
       call fail("cleared communication plans remained ready")
    end if
    local_store_ready = local_block_store_ready()
    if (.not. local_store_ready) then
       call fail("communication-plan teardown cleared local blocks")
    end if
    if (.not. allocated(block_catalog)) then
       call fail("communication-plan teardown cleared block catalogue")
    end if

    call build_block_ghost_exchange_plan
    call check_block_ghost_exchange_plan(.false.)
    call build_block_writeback_plan
    call check_block_writeback_plan(.false.)

    allocation_after_rebuild = &
         block_writeback_plan_allocation_count()
    stage_allocation_after_rebuild = &
         block_writeback_plan%stage_allocations
    if (allocation_after_rebuild /= allocation_before+4_int64) then
       call fail("writeback plan lifecycle allocation count mismatch")
    end if
    if (stage_allocation_after_rebuild /= &
         stage_allocation_before+3_int64) then
       call fail("Domain stage lifecycle allocation count mismatch")
    end if

    call build_block_ghost_exchange_plan
    call build_block_writeback_plan
    if (block_writeback_plan_allocation_count() /= &
         allocation_after_rebuild) then
       call fail("ready lifecycle plan reallocated persistent buffers")
    end if
    if (block_writeback_plan%stage_allocations /= &
         stage_allocation_after_rebuild) then
       call fail("ready lifecycle plan reallocated Domain staging")
    end if

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("rebuilt parallel block state is not ready")
    end if
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_SOL)
    call assert_block_domain_field_family_match(BLOCK_PAYLOAD_WAV_COEFF)

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Parallel block lifecycle for rank ",rank,":"
       write(6,'(a)') "  checkpoint sol/wav_coeff synchronization passed"
       write(6,'(a)') "  checkpoint retained installed block state"
       write(6,'(a)') "  topology-dependent plan invalidation passed"
       write(6,'(a)') "  persistent communication-plan rebuild passed"
       write(6,'(a)') "  rebuilt sol/wav_coeff Domain comparison passed"
       write(6,'(a,/)') "  ready-plan allocation reuse passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Checkpoint/adaptation block lifecycle interface passed"
    end if

  end subroutine check_parallel_block_lifecycle


  subroutine check_parallel_block_scaling (verbose)
    ! Collect repeatable min/average/max work and communication metrics.
    ! This is a read-only validation: it does not rebuild plans, exchange
    ! prognostic values or alter either block or Domain field storage.

    implicit none

    integer, parameter :: SCALE_BLOCK = 1
    integer, parameter :: SCALE_WEIGHT = 2
    integer, parameter :: SCALE_PATCH = 3
    integer, parameter :: SCALE_GHOST_SEND_PEER = 4
    integer, parameter :: SCALE_GHOST_RECV_PEER = 5
    integer, parameter :: SCALE_GHOST_SCALAR_SEND = 6
    integer, parameter :: SCALE_GHOST_SCALAR_RECV = 7
    integer, parameter :: SCALE_GHOST_VECTOR_SEND = 8
    integer, parameter :: SCALE_GHOST_VECTOR_RECV = 9
    integer, parameter :: SCALE_WRITEBACK_SEND_PEER = 10
    integer, parameter :: SCALE_WRITEBACK_RECV_PEER = 11
    integer, parameter :: SCALE_WRITEBACK_SCALAR_SEND = 12
    integer, parameter :: SCALE_WRITEBACK_SCALAR_RECV = 13
    integer, parameter :: SCALE_WRITEBACK_VECTOR_SEND = 14
    integer, parameter :: SCALE_WRITEBACK_VECTOR_RECV = 15
    integer, parameter :: SCALE_PERSISTENT_REAL = 16
    integer, parameter :: SCALE_DOMAIN_STAGE_REAL = 17
    integer, parameter :: SCALE_METRIC_COUNT = 17

    logical, optional, intent(in) :: verbose

    integer :: ierr

    integer(int64) :: accumulator_allocation_after
    integer(int64) :: accumulator_allocation_before
    integer(int64) :: global_max(SCALE_METRIC_COUNT)
    integer(int64) :: global_max_repeat(SCALE_METRIC_COUNT)
    integer(int64) :: global_min(SCALE_METRIC_COUNT)
    integer(int64) :: global_min_repeat(SCALE_METRIC_COUNT)
    integer(int64) :: global_sum(SCALE_METRIC_COUNT)
    integer(int64) :: global_sum_repeat(SCALE_METRIC_COUNT)
    integer(int64) :: local_metric(SCALE_METRIC_COUNT)
    integer(int64) :: local_metric_repeat(SCALE_METRIC_COUNT)
    integer(int64) :: tendency_allocation_after
    integer(int64) :: tendency_allocation_before
    integer(int64) :: writeback_allocation_after
    integer(int64) :: writeback_allocation_before

    logical :: print_summary
    logical :: state_ready

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    state_ready = parallel_block_state_is_ready()
    if (.not. state_ready) then
       call fail("scaling check before parallel block state is ready")
    end if

    writeback_allocation_before = &
         block_writeback_plan_allocation_count()
    tendency_allocation_before = &
         local_block_tendency_allocation_count()
    accumulator_allocation_before = &
         local_block_tendency_accumulator_allocation_count()

    call collect_scaling_snapshot( &
         local_metric,global_min,global_max,global_sum)
    call collect_scaling_snapshot( &
         local_metric_repeat,global_min_repeat,global_max_repeat, &
         global_sum_repeat)

    writeback_allocation_after = &
         block_writeback_plan_allocation_count()
    tendency_allocation_after = &
         local_block_tendency_allocation_count()
    accumulator_allocation_after = &
         local_block_tendency_accumulator_allocation_count()

    if (any(local_metric_repeat /= local_metric) .or. &
         any(global_min_repeat /= global_min) .or. &
         any(global_max_repeat /= global_max) .or. &
         any(global_sum_repeat /= global_sum)) then
       call fail("repeated scaling snapshot changed")
    end if
    if (writeback_allocation_after /= writeback_allocation_before) then
       call fail("scaling snapshot reallocated writeback buffers")
    end if
    if (tendency_allocation_after /= tendency_allocation_before) then
       call fail("scaling snapshot reallocated tendency storage")
    end if
    if (accumulator_allocation_after /= &
         accumulator_allocation_before) then
       call fail("scaling snapshot reallocated accumulator storage")
    end if

    if (global_sum(SCALE_BLOCK) /= int(size(block_catalog),int64)) then
       call fail("scaling snapshot global block count mismatch")
    end if
    if (global_sum(SCALE_WEIGHT) /= &
         sum(int(block_catalog%weight,int64))) then
       call fail("scaling snapshot global block weight mismatch")
    end if
    if (global_sum(SCALE_GHOST_SEND_PEER) /= &
         global_sum(SCALE_GHOST_RECV_PEER)) then
       call fail("scaling snapshot ghost peer count mismatch")
    end if
    if (global_sum(SCALE_GHOST_SCALAR_SEND) /= &
         global_sum(SCALE_GHOST_SCALAR_RECV) .or. &
         global_sum(SCALE_GHOST_VECTOR_SEND) /= &
         global_sum(SCALE_GHOST_VECTOR_RECV)) then
       call fail("scaling snapshot ghost payload mismatch")
    end if
    if (global_sum(SCALE_WRITEBACK_SEND_PEER) /= &
         global_sum(SCALE_WRITEBACK_RECV_PEER)) then
       call fail("scaling snapshot writeback peer count mismatch")
    end if
    if (global_sum(SCALE_WRITEBACK_SCALAR_SEND) /= &
         global_sum(SCALE_WRITEBACK_SCALAR_RECV) .or. &
         global_sum(SCALE_WRITEBACK_VECTOR_SEND) /= &
         global_sum(SCALE_WRITEBACK_VECTOR_RECV)) then
       call fail("scaling snapshot writeback payload mismatch")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Parallel block scaling metrics for rank ",rank,":"
       write(6,'(a,i0)') "  final-owner blocks = ", &
            local_metric(SCALE_BLOCK)
       write(6,'(a,i0)') "  final-owner weight = ", &
            local_metric(SCALE_WEIGHT)
       write(6,'(a,i0)') "  installed patches = ", &
            local_metric(SCALE_PATCH)
       write(6,'(a,i0,a,i0)') "  ghost send/receive peers = ", &
            local_metric(SCALE_GHOST_SEND_PEER)," / ", &
            local_metric(SCALE_GHOST_RECV_PEER)
       write(6,'(a,i0,a,i0)') "  ghost scalar send/receive values = ", &
            local_metric(SCALE_GHOST_SCALAR_SEND)," / ", &
            local_metric(SCALE_GHOST_SCALAR_RECV)
       write(6,'(a,i0,a,i0)') "  ghost vector send/receive values = ", &
            local_metric(SCALE_GHOST_VECTOR_SEND)," / ", &
            local_metric(SCALE_GHOST_VECTOR_RECV)
       write(6,'(a,i0,a,i0)') "  writeback send/receive peers = ", &
            local_metric(SCALE_WRITEBACK_SEND_PEER)," / ", &
            local_metric(SCALE_WRITEBACK_RECV_PEER)
       write(6,'(a,i0,a,i0)') &
            "  writeback scalar send/receive values = ", &
            local_metric(SCALE_WRITEBACK_SCALAR_SEND)," / ", &
            local_metric(SCALE_WRITEBACK_SCALAR_RECV)
       write(6,'(a,i0,a,i0)') &
            "  writeback vector send/receive values = ", &
            local_metric(SCALE_WRITEBACK_VECTOR_SEND)," / ", &
            local_metric(SCALE_WRITEBACK_VECTOR_RECV)
       write(6,'(a,i0)') "  persistent real-value capacity = ", &
            local_metric(SCALE_PERSISTENT_REAL)
       write(6,'(a,/)') "  repeated read-only snapshot passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a)') "Parallel block global scaling summary:"
       write(6,'(a,i0)') "  MPI ranks = ",n_process
       call print_min_average_max("final-owner blocks",SCALE_BLOCK)
       call print_min_average_max("final-owner weight",SCALE_WEIGHT)
       call print_min_average_max("installed patches",SCALE_PATCH)
       call print_min_average_max( &
            "ghost send peers",SCALE_GHOST_SEND_PEER)
       call print_min_average_max( &
            "ghost receive peers",SCALE_GHOST_RECV_PEER)
       call print_min_average_max( &
            "writeback send peers",SCALE_WRITEBACK_SEND_PEER)
       call print_min_average_max( &
            "writeback receive peers",SCALE_WRITEBACK_RECV_PEER)
       call print_min_average_max( &
            "persistent real values",SCALE_PERSISTENT_REAL)
       write(6,'(a,i0)') "  global ghost scalar values = ", &
            global_sum(SCALE_GHOST_SCALAR_SEND)
       write(6,'(a,i0)') "  global ghost vector values = ", &
            global_sum(SCALE_GHOST_VECTOR_SEND)
       write(6,'(a,i0)') "  global writeback scalar values = ", &
            global_sum(SCALE_WRITEBACK_SCALAR_SEND)
       write(6,'(a,i0)') "  global writeback vector values = ", &
            global_sum(SCALE_WRITEBACK_VECTOR_SEND)
       write(6,'(a,i0)') "  global Domain stage real values = ", &
            global_sum(SCALE_DOMAIN_STAGE_REAL)
       write(6,'(a)') "  global route and payload balance passed"
       write(6,'(a)') "  repeated allocation-free snapshot passed"
       write(6,'(a,/)') "Parallel block scaling validation passed"
    end if

  contains

    subroutine collect_scaling_snapshot ( &
         local_value,minimum_value,maximum_value,total_value)

      implicit none

      integer(int64), intent(out) :: local_value(SCALE_METRIC_COUNT)
      integer(int64), intent(out) :: maximum_value(SCALE_METRIC_COUNT)
      integer(int64), intent(out) :: minimum_value(SCALE_METRIC_COUNT)
      integer(int64), intent(out) :: total_value(SCALE_METRIC_COUNT)

      integer :: b
      integer :: i
      integer :: local_block_count

      local_value = 0_int64
      local_block_count = n_local_blocks()
      local_value(SCALE_BLOCK) = int(local_block_count,int64)

      do i = 1,local_block_count
         b = local_block_catalog(i)
         local_value(SCALE_WEIGHT) = local_value(SCALE_WEIGHT) + &
              int(block_catalog(b)%weight,int64)
         local_value(SCALE_PATCH) = local_value(SCALE_PATCH) + &
              int(local_block_patch_count(b),int64)
      end do

      local_value(SCALE_GHOST_SEND_PEER) = int(count( &
           ghost_exchange_plan%scalar_send_count > 0),int64)
      local_value(SCALE_GHOST_RECV_PEER) = int(count( &
           ghost_exchange_plan%scalar_recv_count > 0),int64)
      local_value(SCALE_GHOST_SCALAR_SEND) = sum(int( &
           ghost_exchange_plan%scalar_send_count,int64))
      local_value(SCALE_GHOST_SCALAR_RECV) = sum(int( &
           ghost_exchange_plan%scalar_recv_count,int64))
      local_value(SCALE_GHOST_VECTOR_SEND) = sum(int( &
           ghost_exchange_plan%vector_send_count,int64))
      local_value(SCALE_GHOST_VECTOR_RECV) = sum(int( &
           ghost_exchange_plan%vector_recv_count,int64))

      local_value(SCALE_WRITEBACK_SEND_PEER) = int(count( &
           block_writeback_plan%scalar_send_count > 0),int64)
      local_value(SCALE_WRITEBACK_RECV_PEER) = int(count( &
           block_writeback_plan%scalar_recv_count > 0),int64)
      local_value(SCALE_WRITEBACK_SCALAR_SEND) = sum(int( &
           block_writeback_plan%scalar_send_count,int64))
      local_value(SCALE_WRITEBACK_SCALAR_RECV) = sum(int( &
           block_writeback_plan%scalar_recv_count,int64))
      local_value(SCALE_WRITEBACK_VECTOR_SEND) = sum(int( &
           block_writeback_plan%vector_send_count,int64))
      local_value(SCALE_WRITEBACK_VECTOR_RECV) = sum(int( &
           block_writeback_plan%vector_recv_count,int64))

      local_value(SCALE_PERSISTENT_REAL) = int(size( &
           ghost_exchange_plan%scalar_send_buffer),int64) + int(size( &
           ghost_exchange_plan%scalar_recv_buffer),int64) + int(size( &
           ghost_exchange_plan%scalar_patch_buffer),int64) + int(size( &
           ghost_exchange_plan%vector_send_buffer),int64) + int(size( &
           ghost_exchange_plan%vector_recv_buffer),int64) + int(size( &
           ghost_exchange_plan%vector_patch_buffer),int64) + int(size( &
           block_writeback_plan%scalar_send_buffer),int64) + int(size( &
           block_writeback_plan%scalar_recv_buffer),int64) + int(size( &
           block_writeback_plan%vector_send_buffer),int64) + int(size( &
           block_writeback_plan%vector_recv_buffer),int64) + int(size( &
           block_writeback_plan%scalar_domain_stage),int64) + int(size( &
           block_writeback_plan%vector_domain_stage),int64)
      local_value(SCALE_DOMAIN_STAGE_REAL) = int(size( &
           block_writeback_plan%scalar_domain_stage),int64) + int(size( &
           block_writeback_plan%vector_domain_stage),int64)

      if (any(local_value < 0_int64)) then
         call fail("scaling snapshot contains a negative metric")
      end if

      call MPI_Allreduce(local_value,minimum_value,SCALE_METRIC_COUNT, &
           MPI_INTEGER8,MPI_MIN,comm,ierr)
      call check_mpi(ierr,"MPI_Allreduce scaling minima")
      call MPI_Allreduce(local_value,maximum_value,SCALE_METRIC_COUNT, &
           MPI_INTEGER8,MPI_MAX,comm,ierr)
      call check_mpi(ierr,"MPI_Allreduce scaling maxima")
      call MPI_Allreduce(local_value,total_value,SCALE_METRIC_COUNT, &
           MPI_INTEGER8,MPI_SUM,comm,ierr)
      call check_mpi(ierr,"MPI_Allreduce scaling totals")

    end subroutine collect_scaling_snapshot


    subroutine print_min_average_max (description,metric_index)

      implicit none

      character(*), intent(in) :: description
      integer, intent(in) :: metric_index

      real(dp) :: average
      real(dp) :: maximum_over_average

      average = real(global_sum(metric_index),dp) / real(n_process,dp)
      maximum_over_average = 0.0_dp
      if (abs(average) > 0.0_dp) then
         maximum_over_average = &
              real(global_max(metric_index),dp) / average
      end if

      write(6,'(2a,i0,a,f12.2,a,i0,a,f10.4)') "  ", &
           trim(description)//": min/avg/max = ", &
           global_min(metric_index)," / ",average," / ", &
           global_max(metric_index),"  max/avg = ", &
           maximum_over_average

    end subroutine print_min_average_max

  end subroutine check_parallel_block_scaling


  subroutine accumulate_block_hydrostatic_consumer ( &
       catalog_index,n_patch,surface_pressure,dynamic_exner, &
       air_temperature,context)
    ! Production-driver test kernel. Accumulate a complete inventory
    ! directly from the read-only block arrays supplied by the driver.

    implicit none

    integer, intent(in) :: catalog_index
    integer, intent(in) :: n_patch
    real(dp), intent(in) :: surface_pressure(:)
    real(dp), intent(in) :: dynamic_exner(:)
    real(dp), intent(in) :: air_temperature(:)
    class(*), intent(inout) :: context

    if (catalog_index < 1) then
       call fail("hydrostatic consumer received invalid catalogue index")
    end if
    if (n_patch < 1) then
       call fail("hydrostatic consumer received empty block")
    end if
    if (size(surface_pressure) /= n_patch*PATCH_SIZE**2 .or. &
         size(dynamic_exner) /= n_patch*zlevels*PATCH_SIZE**2 .or. &
         size(air_temperature) /= n_patch*zlevels*PATCH_SIZE**2) then
       call fail("hydrostatic consumer received invalid field extents")
    end if

    select type (statistics => context)
    type is (Block_Hydrostatic_Traversal_Context)
       statistics%block_count = statistics%block_count + 1_int64
       statistics%patch_count = statistics%patch_count + &
            int(n_patch,int64)
       statistics%surface_count = statistics%surface_count + &
            int(size(surface_pressure),int64)
       statistics%column_count = statistics%column_count + &
            int(size(dynamic_exner),int64)

       statistics%surface_moment(1) = &
            statistics%surface_moment(1) + sum(surface_pressure)
       statistics%surface_moment(2) = &
            statistics%surface_moment(2) + &
            sum(abs(surface_pressure))
       statistics%surface_moment(3) = &
            statistics%surface_moment(3) + sum(surface_pressure**2)

       statistics%exner_moment(1) = &
            statistics%exner_moment(1) + sum(dynamic_exner)
       statistics%exner_moment(2) = &
            statistics%exner_moment(2) + sum(abs(dynamic_exner))
       statistics%exner_moment(3) = &
            statistics%exner_moment(3) + sum(dynamic_exner**2)

       statistics%temperature_moment(1) = &
            statistics%temperature_moment(1) + sum(air_temperature)
       statistics%temperature_moment(2) = &
            statistics%temperature_moment(2) + &
            sum(abs(air_temperature))
       statistics%temperature_moment(3) = &
            statistics%temperature_moment(3) + &
            sum(air_temperature**2)
    class default
       call fail("hydrostatic consumer received invalid context")
    end select

  end subroutine accumulate_block_hydrostatic_consumer


  subroutine check_block_hydrostatic_consumer (verbose)
    ! Exercise the production traversal interface and compare its direct
    ! block-array inventory with the established cache statistics.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: catalog_index
    integer :: local_index

    integer(int64) :: column_count
    integer(int64) :: expected_patch_count
    integer(int64) :: refresh_count_before
    integer(int64) :: surface_count

    real(dp) :: exner_moment(3)
    real(dp) :: surface_moment(3)
    real(dp) :: temperature_moment(3)

    logical :: print_summary

    type(Block_Hydrostatic_Traversal_Context) :: statistics

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. compressible) then
       if (print_summary) then
          write(6,'(/,a,i0,a)') &
               "Block thermodynamic traversal for rank ", rank, ":"
          write(6,'(a,/)') &
               "  skipped for incompressible configuration"
       end if
       return
    end if

    call local_block_hydrostatic_statistics( &
         surface_count,column_count,surface_moment,exner_moment, &
         temperature_moment)

    refresh_count_before = local_block_hydrostatic_refresh_count()
    statistics = Block_Hydrostatic_Traversal_Context()

    call apply_local_block_hydrostatic_consumer( &
         accumulate_block_hydrostatic_consumer,statistics)

    if (local_block_hydrostatic_refresh_count() /= &
         refresh_count_before) then
       call fail("production hydrostatic consumer refreshed cache")
    end if

    expected_patch_count = 0_int64
    do local_index = 1,n_local_blocks()
       catalog_index = local_block_catalog(local_index)
       expected_patch_count = expected_patch_count + &
            int(local_block_patch_count(catalog_index),int64)
    end do

    if (statistics%block_count /= int(n_local_blocks(),int64) .or. &
         statistics%patch_count /= expected_patch_count) then
       call fail("production hydrostatic consumer traversal mismatch")
    end if
    if (statistics%surface_count /= surface_count .or. &
         statistics%column_count /= column_count) then
       call fail("production hydrostatic consumer count mismatch")
    end if
    if (.not. field_moments_match( &
         statistics%surface_moment,surface_moment,surface_count)) then
       call fail("production hydrostatic consumer surface mismatch")
    end if
    if (.not. field_moments_match( &
         statistics%exner_moment,exner_moment,column_count)) then
       call fail("production hydrostatic consumer Exner mismatch")
    end if
    if (.not. field_moments_match( &
         statistics%temperature_moment,temperature_moment, &
         column_count)) then
       call fail("production hydrostatic consumer temperature mismatch")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Block thermodynamic traversal for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  local blocks consumed  = ", statistics%block_count
       write(6,'(a,i0)') &
            "  local patches consumed = ", statistics%patch_count
       write(6,'(a,/)') &
            "  direct production traversal checks passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,/)') &
            "Production block thermodynamic consumer traversal passed"
    end if

  end subroutine check_block_hydrostatic_consumer


  subroutine check_block_hydrostatic_reconstruction (verbose)
    ! Validate persistent production block hydrostatic state in shadow
    ! mode against an independent reconstruction from legacy fields.

    implicit none

    logical, optional, intent(in) :: verbose

    integer :: b
    integer :: d
    integer :: ierr

    integer(int64) :: block_surface_count_local
    integer(int64) :: block_surface_count_global
    integer(int64) :: block_column_count_local
    integer(int64) :: block_column_count_global
    integer(int64) :: domain_surface_count_local
    integer(int64) :: domain_surface_count_global
    integer(int64) :: domain_column_count_local
    integer(int64) :: domain_column_count_global

    real(dp) :: block_surface_moment_local(3)
    real(dp) :: block_surface_moment_global(3)
    real(dp) :: block_exner_moment_local(3)
    real(dp) :: block_exner_moment_global(3)
    real(dp) :: block_temperature_moment_local(3)
    real(dp) :: block_temperature_moment_global(3)
    real(dp) :: domain_surface_moment_local(3)
    real(dp) :: domain_surface_moment_global(3)
    real(dp) :: domain_exner_moment_local(3)
    real(dp) :: domain_exner_moment_global(3)
    real(dp) :: domain_temperature_moment_local(3)
    real(dp) :: domain_temperature_moment_global(3)

    logical :: print_summary

    print_summary = .true.
    if (present(verbose)) print_summary = verbose

    if (.not. compressible) then
       if (print_summary) then
          write(6,'(/,a,i0,a)') &
               "Hydrostatic block reconstruction for rank ", rank, ":"
          write(6,'(a,/)') &
               "  skipped for incompressible configuration"
       end if
       return
    end if

    if (.not. local_block_store_ready()) then
       call fail("hydrostatic reconstruction before block installation")
    end if
    if (.not. local_block_hydrostatic_state_ready()) then
       call fail("persistent block hydrostatic state is not ready")
    end if

    call local_block_hydrostatic_statistics( &
         block_surface_count_local,block_column_count_local, &
         block_surface_moment_local,block_exner_moment_local, &
         block_temperature_moment_local)

    domain_surface_count_local      = 0_int64
    domain_column_count_local       = 0_int64
    domain_surface_moment_local     = 0.0_dp
    domain_exner_moment_local       = 0.0_dp
    domain_temperature_moment_local = 0.0_dp

    do b = 1, size(block_catalog)

       if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

       d = loc_id(block_catalog(b)%root_domain+1) + 1

       if (d < 1 .or. d > size(grid)) then
          call fail("invalid source Domain in hydrostatic reconstruction")
       end if

       call accumulate_domain_subtree_hydrostatic( &
            d,block_catalog(b)%root_patch, &
            domain_surface_count_local,domain_column_count_local, &
            domain_surface_moment_local,domain_exner_moment_local, &
            domain_temperature_moment_local)

    end do

    call MPI_Allreduce( &
         block_surface_count_local,block_surface_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block surface-pressure count")

    call MPI_Allreduce( &
         block_column_count_local,block_column_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block hydrostatic-column count")

    call MPI_Allreduce( &
         domain_surface_count_local,domain_surface_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain surface-pressure count")

    call MPI_Allreduce( &
         domain_column_count_local,domain_column_count_global,1, &
         MPI_INTEGER8,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain hydrostatic-column count")

    call MPI_Allreduce( &
         block_surface_moment_local,block_surface_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block surface-pressure moments")

    call MPI_Allreduce( &
         block_exner_moment_local,block_exner_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block Exner moments")

    call MPI_Allreduce( &
         block_temperature_moment_local, &
         block_temperature_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce block temperature moments")

    call MPI_Allreduce( &
         domain_surface_moment_local,domain_surface_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain surface-pressure moments")

    call MPI_Allreduce( &
         domain_exner_moment_local,domain_exner_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain Exner moments")

    call MPI_Allreduce( &
         domain_temperature_moment_local, &
         domain_temperature_moment_global,3, &
         MPI_DOUBLE_PRECISION,MPI_SUM,comm,ierr)
    call check_mpi(ierr,"MPI_Allreduce Domain temperature moments")

    if (block_surface_count_global /= domain_surface_count_global .or. &
         block_column_count_global /= domain_column_count_global) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain hydrostatic-count mismatch:"
          write(error_unit,'(a,2(i0,1x))') &
               "  block surface/column counts  = ", &
               block_surface_count_global,block_column_count_global
          write(error_unit,'(a,2(i0,1x))') &
               "  Domain surface/column counts = ", &
               domain_surface_count_global,domain_column_count_global
       end if

       call fail("block and Domain hydrostatic counts differ")
    end if

    if (.not. field_moments_match( &
         block_surface_moment_global,domain_surface_moment_global, &
         block_surface_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain surface-pressure-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_surface_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_surface_moment_global
       end if

       call fail("block and Domain surface-pressure moments differ")
    end if

    if (.not. field_moments_match( &
         block_exner_moment_global,domain_exner_moment_global, &
         block_column_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain dynamic-Exner-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", block_exner_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", domain_exner_moment_global
       end if

       call fail("block and Domain dynamic Exner moments differ")
    end if

    if (.not. field_moments_match( &
         block_temperature_moment_global, &
         domain_temperature_moment_global, &
         block_column_count_global)) then

       if (rank == 0) then
          write(error_unit,'(/,a)') &
               "Block/Domain air-temperature-moment mismatch:"
          write(error_unit,'(a,3(es24.16,1x))') &
               "  block moments  = ", &
               block_temperature_moment_global
          write(error_unit,'(a,3(es24.16,1x))') &
               "  Domain moments = ", &
               domain_temperature_moment_global
       end if

       call fail("block and Domain air-temperature moments differ")
    end if

    if (print_summary) then
       write(6,'(/,a,i0,a)') &
            "Hydrostatic block reconstruction for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  local surface-pressure values = ", &
            block_surface_count_local
       write(6,'(a,i0)') &
            "  local column diagnostic values = ", &
            block_column_count_local
       write(6,'(a)') &
            "  global hydrostatic reconstruction checks passed"
       write(6,'(a,/)') &
            "  persistent production state shadow comparison passed"
    end if

    if (print_summary .and. rank == 0) then
       write(6,'(/,a,i0)') &
            "Global surface-pressure values reconstructed = ", &
            block_surface_count_global
       write(6,'(a,i0)') &
            "Global dynamic Exner/temperature values reconstructed = ", &
            block_column_count_global
       write(6,'(a,3(es24.16,1x))') &
            "Global surface-pressure moments verified = ", &
            block_surface_moment_global
       write(6,'(a,3(es24.16,1x))') &
            "Global dynamic Exner moments verified = ", &
            block_exner_moment_global
       write(6,'(a,3(es24.16,1x))') &
            "Global air-temperature moments verified = ", &
            block_temperature_moment_global
       write(6,'(a,/)') &
            "Persistent production block hydrostatic state matches legacy Domain state"
    end if

  end subroutine check_block_hydrostatic_reconstruction


  logical function field_moments_match (first,second,n_value) &
       result(match)
    ! Allow only the floating-point reassociation error introduced by
    ! summing an identical value inventory in different block/rank
    ! orders.

    implicit none

    real(dp), intent(in) :: first(3)
    real(dp), intent(in) :: second(3)
    integer(int64), intent(in) :: n_value

    real(dp) :: factor
    real(dp) :: scale

    factor = 256.0_dp * epsilon(1.0_dp) * &
         real(max(1_int64,n_value),dp)

    scale = max(1.0_dp,first(2),second(2))

    match = abs(first(1)-second(1)) <= factor*scale
    if (.not. match) return

    match = abs(first(2)-second(2)) <= factor*scale
    if (.not. match) return

    scale = max(1.0_dp,first(3),second(3))
    match = abs(first(3)-second(3)) <= factor*scale

  end function field_moments_match


  recursive subroutine accumulate_domain_subtree_fields ( &
       d,p,v_scalar,n_scalar_variable,v_vector,k_field,n_field_level, &
       mult_scalar,mult_vector,field_count,field_moment, &
       mean_field_count,mean_field_moment, &
       wavelet_field_count,wavelet_field_moment, &
       tke_level,n_tke_level,tke_count,tke_moment, &
       wav_tke_count,wav_tke_moment, &
       topography_count,topography_moment)
    ! Accumulate one catalogue-rooted subtree from the authoritative
    ! Domain representation using the same patch coverage copied into
    ! Block_Data.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p
    integer, intent(in) :: v_scalar
    integer, intent(in) :: n_scalar_variable
    integer, intent(in) :: v_vector
    integer, intent(in) :: k_field
    integer, intent(in) :: n_field_level
    integer, intent(in) :: mult_scalar
    integer, intent(in) :: mult_vector
    integer, intent(in) :: tke_level
    integer, intent(in) :: n_tke_level

    integer(int64), intent(inout) :: field_count(2)
    real(dp), intent(inout) :: field_moment(3,2)
    integer(int64), intent(inout) :: mean_field_count(2)
    real(dp), intent(inout) :: mean_field_moment(3,2)
    integer(int64), intent(inout) :: wavelet_field_count(2)
    real(dp), intent(inout) :: wavelet_field_moment(3,2)
    integer(int64), intent(inout) :: tke_count
    real(dp), intent(inout) :: tke_moment(3)
    integer(int64), intent(inout) :: wav_tke_count
    real(dp), intent(inout) :: wav_tke_moment(3)
    integer(int64), intent(inout) :: topography_count
    real(dp), intent(inout) :: topography_moment(3)

    integer :: c
    integer :: field_level
    integer :: level_slot
    integer :: n_value
    integer :: p_child
    integer :: scalar_id
    integer :: scalar_slot
    integer :: start

    if (p < 0 .or. p >= grid(d)%patch%length) then
       call fail("invalid patch in Domain subtree field inventory")
    end if

    if (grid(d)%patch%elts(p+1)%deleted) return

    do scalar_slot = 1, n_scalar_variable

       scalar_id = v_scalar + scalar_slot - 1

       do level_slot = 1, n_field_level

          field_level = k_field + level_slot - 1
          start = mult_scalar * grid(d)%patch%elts(p+1)%elts_start
          n_value = mult_scalar * PATCH_SIZE**2

          if (start < 0 .or. &
               start+n_value > &
               size(sol(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy scalar patch extent is invalid")
          end if

          field_count(1) = field_count(1) + int(n_value,int64)
          field_moment(1,1) = field_moment(1,1) + &
               sum(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          field_moment(2,1) = field_moment(2,1) + &
               sum(abs(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          field_moment(3,1) = field_moment(3,1) + &
               sum(sol(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

          if (start+n_value > &
               size(wav_coeff(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy scalar wavelet patch extent is invalid")
          end if

          wavelet_field_count(1) = &
               wavelet_field_count(1) + int(n_value,int64)
          wavelet_field_moment(1,1) = wavelet_field_moment(1,1) + &
               sum(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          wavelet_field_moment(2,1) = wavelet_field_moment(2,1) + &
               sum(abs(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          wavelet_field_moment(3,1) = wavelet_field_moment(3,1) + &
               sum(wav_coeff(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

          if (start+n_value > &
               size(sol_mean(scalar_id,field_level)%data(d)%elts)) then
             call fail("legacy mean scalar patch extent is invalid")
          end if

          mean_field_count(1) = &
               mean_field_count(1) + int(n_value,int64)
          mean_field_moment(1,1) = mean_field_moment(1,1) + &
               sum(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value))
          mean_field_moment(2,1) = mean_field_moment(2,1) + &
               sum(abs(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)))
          mean_field_moment(3,1) = mean_field_moment(3,1) + &
               sum(sol_mean(scalar_id,field_level)%data(d)%elts( &
               start+1:start+n_value)**2)

       end do

    end do

    do level_slot = 1, n_field_level

       field_level = k_field + level_slot - 1
       start = mult_vector * grid(d)%patch%elts(p+1)%elts_start
       n_value = mult_vector * PATCH_SIZE**2

       if (start < 0 .or. &
            start+n_value > &
            size(sol(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy vector patch extent is invalid")
       end if

       field_count(2) = field_count(2) + int(n_value,int64)
       field_moment(1,2) = field_moment(1,2) + &
            sum(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       field_moment(2,2) = field_moment(2,2) + &
            sum(abs(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       field_moment(3,2) = field_moment(3,2) + &
            sum(sol(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

       if (start+n_value > &
            size(wav_coeff(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy vector wavelet patch extent is invalid")
       end if

       wavelet_field_count(2) = &
            wavelet_field_count(2) + int(n_value,int64)
       wavelet_field_moment(1,2) = wavelet_field_moment(1,2) + &
            sum(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       wavelet_field_moment(2,2) = wavelet_field_moment(2,2) + &
            sum(abs(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       wavelet_field_moment(3,2) = wavelet_field_moment(3,2) + &
            sum(wav_coeff(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

       if (start+n_value > &
            size(sol_mean(v_vector,field_level)%data(d)%elts)) then
          call fail("legacy mean vector patch extent is invalid")
       end if

       mean_field_count(2) = &
            mean_field_count(2) + int(n_value,int64)
       mean_field_moment(1,2) = mean_field_moment(1,2) + &
            sum(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value))
       mean_field_moment(2,2) = mean_field_moment(2,2) + &
            sum(abs(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)))
       mean_field_moment(3,2) = mean_field_moment(3,2) + &
            sum(sol_mean(v_vector,field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

    end do

    do level_slot = 1, n_tke_level

       field_level = tke_level + level_slot - 1
       start = grid(d)%patch%elts(p+1)%elts_start
       n_value = PATCH_SIZE**2

       if (start < 0 .or. &
            start+n_value > size(tke(field_level)%data(d)%elts)) then
          call fail("legacy tke patch extent is invalid")
       end if

       if (start+n_value > &
            size(wav_tke(field_level)%data(d)%elts)) then
          call fail("legacy wav_tke patch extent is invalid")
       end if

       tke_count = tke_count + int(n_value,int64)
       tke_moment(1) = tke_moment(1) + &
            sum(tke(field_level)%data(d)%elts(start+1:start+n_value))
       tke_moment(2) = tke_moment(2) + &
            sum(abs(tke(field_level)%data(d)%elts( &
            start+1:start+n_value)))
       tke_moment(3) = tke_moment(3) + &
            sum(tke(field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

       wav_tke_count = wav_tke_count + int(n_value,int64)
       wav_tke_moment(1) = wav_tke_moment(1) + &
            sum(wav_tke(field_level)%data(d)%elts( &
            start+1:start+n_value))
       wav_tke_moment(2) = wav_tke_moment(2) + &
            sum(abs(wav_tke(field_level)%data(d)%elts( &
            start+1:start+n_value)))
       wav_tke_moment(3) = wav_tke_moment(3) + &
            sum(wav_tke(field_level)%data(d)%elts( &
            start+1:start+n_value)**2)

    end do

    start = grid(d)%patch%elts(p+1)%elts_start
    n_value = PATCH_SIZE**2

    if (start < 0 .or. &
         start+n_value > size(topography%data(d)%elts)) then
       call fail("legacy topography patch extent is invalid")
    end if

    topography_count = topography_count + int(n_value,int64)
    topography_moment(1) = topography_moment(1) + &
         sum(topography%data(d)%elts(start+1:start+n_value))
    topography_moment(2) = topography_moment(2) + &
         sum(abs(topography%data(d)%elts(start+1:start+n_value)))
    topography_moment(3) = topography_moment(3) + &
         sum(topography%data(d)%elts(start+1:start+n_value)**2)

    do c = 1, N_CHDRN

       p_child = grid(d)%patch%elts(p+1)%children(c)
       if (p_child == 0) cycle

       call accumulate_domain_subtree_fields( &
            d,p_child,v_scalar,n_scalar_variable,v_vector,k_field, &
            n_field_level,mult_scalar,mult_vector, &
            field_count,field_moment, &
            mean_field_count,mean_field_moment, &
            wavelet_field_count,wavelet_field_moment, &
            tke_level,n_tke_level,tke_count,tke_moment, &
            wav_tke_count,wav_tke_moment, &
            topography_count,topography_moment)

    end do

  end subroutine accumulate_domain_subtree_fields


  recursive subroutine accumulate_domain_subtree_hydrostatic ( &
       d,p,surface_count,column_count,surface_moment,exner_moment, &
       temperature_moment)
    ! Independently reconstruct compressible hydrostatic diagnostics
    ! over one catalogue-rooted legacy Domain subtree.

    implicit none

    integer, intent(in) :: d
    integer, intent(in) :: p

    integer(int64), intent(inout) :: surface_count
    integer(int64), intent(inout) :: column_count

    real(dp), intent(inout) :: surface_moment(3)
    real(dp), intent(inout) :: exner_moment(3)
    real(dp), intent(inout) :: temperature_moment(3)

    integer :: c
    integer :: i
    integer :: k
    integer :: node_index
    integer :: p_child
    integer :: start

    real(dp) :: exner_value
    real(dp) :: layer_pressure
    real(dp) :: pressure_factor
    real(dp) :: pressure_lower
    real(dp) :: pressure_upper
    real(dp) :: rho_dz
    real(dp) :: rho_dz_theta
    real(dp) :: surface_pressure
    real(dp) :: temperature_value

    if (p < 0 .or. p >= grid(d)%patch%length) then
       call fail("invalid patch in Domain hydrostatic reconstruction")
    end if

    if (grid(d)%patch%elts(p+1)%deleted) return

    start = grid(d)%patch%elts(p+1)%elts_start

    if (start < 0 .or. &
         start+PATCH_SIZE**2 > &
         size(sol(S_MASS,1)%data(d)%elts)) then
       call fail("legacy hydrostatic patch extent is invalid")
    end if

    do i = 1, PATCH_SIZE**2

       node_index = start + i
       pressure_lower = p_top

       do k = 1, zlevels
          rho_dz = sol(S_MASS,k)%data(d)%elts(node_index) + &
               sol_mean(S_MASS,k)%data(d)%elts(node_index)

          if (rho_dz <= 0.0_dp) then
             call fail("nonpositive legacy mass in reconstruction")
          end if

          pressure_lower = pressure_lower + grav_accel*rho_dz
       end do

       surface_pressure = pressure_lower
       surface_count = surface_count + 1_int64
       surface_moment(1) = surface_moment(1) + surface_pressure
       surface_moment(2) = surface_moment(2) + abs(surface_pressure)
       surface_moment(3) = surface_moment(3) + surface_pressure**2

       do k = 1, zlevels
          rho_dz = sol(S_MASS,k)%data(d)%elts(node_index) + &
               sol_mean(S_MASS,k)%data(d)%elts(node_index)
          rho_dz_theta = &
               sol(S_TEMP,k)%data(d)%elts(node_index) + &
               sol_mean(S_TEMP,k)%data(d)%elts(node_index)

          pressure_upper = pressure_lower - grav_accel*rho_dz
          layer_pressure = 0.5_dp*(pressure_lower+pressure_upper)

          if (layer_pressure <= 0.0_dp) then
             call fail("nonpositive legacy pressure in reconstruction")
          end if

          pressure_factor = (layer_pressure/p_0)**kappa
          exner_value = c_p*pressure_factor
          temperature_value = rho_dz_theta/rho_dz*pressure_factor

          column_count = column_count + 1_int64

          exner_moment(1) = exner_moment(1) + exner_value
          exner_moment(2) = exner_moment(2) + abs(exner_value)
          exner_moment(3) = exner_moment(3) + exner_value**2

          temperature_moment(1) = &
               temperature_moment(1) + temperature_value
          temperature_moment(2) = &
               temperature_moment(2) + abs(temperature_value)
          temperature_moment(3) = &
               temperature_moment(3) + temperature_value**2

          pressure_lower = pressure_upper
       end do

    end do

    do c = 1, N_CHDRN

       p_child = grid(d)%patch%elts(p+1)%children(c)
       if (p_child == 0) cycle

       call accumulate_domain_subtree_hydrostatic( &
            d,p_child,surface_count,column_count,surface_moment, &
            exner_moment,temperature_moment)

    end do

  end subroutine accumulate_domain_subtree_hydrostatic


  subroutine build_block_migration_manifest (manifest)
    ! Build and exchange a migration manifest containing only indices
    ! into the replicated block_catalog. No block payload is moved and
    ! neither grid nor any persistent block object is modified.

    implicit none

    type(Block_Migration_Manifest), intent(out) :: manifest

    integer :: b
    integer :: destination
    integer :: ierr
    integer :: r
    integer :: slot
    integer :: source
    integer, allocatable :: cursor(:)

    if (.not. allocated(block_catalog)) then
       call fail("block_catalog is not allocated")
    end if

    if (.not. allocated(owner)) then
       call fail("domain owner map is not allocated")
    end if

    allocate(manifest%send_count(n_process))
    allocate(manifest%recv_count(n_process))
    allocate(manifest%send_displ(n_process))
    allocate(manifest%recv_displ(n_process))

    manifest%send_count = 0
    manifest%recv_count = 0
    manifest%send_displ = 0
    manifest%recv_displ = 0
    manifest%validated  = .false.

    do b = 1, size(block_catalog)
       source      = source_rank(b)
       destination = block_catalog(b)%owner

       call check_rank(destination,"destination owner")

       if (source == rank .and. destination /= rank) then
          manifest%send_count(destination+1) = &
               manifest%send_count(destination+1) + 1
       end if
    end do

    call MPI_Alltoall( &
         manifest%send_count, 1, MPI_INTEGER, &
         manifest%recv_count, 1, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoall migration counts")

    do r = 2, n_process
       manifest%send_displ(r) = manifest%send_displ(r-1) + &
            manifest%send_count(r-1)
       manifest%recv_displ(r) = manifest%recv_displ(r-1) + &
            manifest%recv_count(r-1)
    end do

    manifest%n_send = sum(manifest%send_count)
    manifest%n_recv = sum(manifest%recv_count)

    ! A one-element allocation avoids passing an unallocated or
    ! zero-sized buffer to MPI implementations with older Fortran
    ! choice-buffer handling. n_send/n_recv remain the true lengths.
    allocate(manifest%send_block(max(1,manifest%n_send)))
    allocate(manifest%recv_block(max(1,manifest%n_recv)))

    manifest%send_block = 0
    manifest%recv_block = 0

    allocate(cursor(n_process))
    cursor = manifest%send_displ

    do b = 1, size(block_catalog)
       source      = source_rank(b)
       destination = block_catalog(b)%owner

       if (source /= rank .or. destination == rank) cycle

       slot = cursor(destination+1) + 1

       if (slot < 1 .or. slot > manifest%n_send) then
          call fail("send-manifest fill position is invalid")
       end if

       manifest%send_block(slot) = b
       cursor(destination+1) = cursor(destination+1) + 1
    end do

    do r = 1, n_process
       if (cursor(r) /= manifest%send_displ(r) + &
            manifest%send_count(r)) then
          call fail("send-manifest count mismatch")
       end if
    end do

    deallocate(cursor)

    call MPI_Alltoallv( &
         manifest%send_block, manifest%send_count, &
         manifest%send_displ, MPI_INTEGER, &
         manifest%recv_block, manifest%recv_count, &
         manifest%recv_displ, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv migration manifest")

  end subroutine build_block_migration_manifest


  subroutine check_block_migration_manifest (manifest, verbose)
    ! Verify that every manifest entry has the expected source rank and
    ! this rank as its prospective destination. Global counts and a
    ! catalogue-index checksum must agree with the replicated catalogue.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    logical, optional, intent(in)                   :: verbose

    integer :: b
    integer :: destination
    integer :: expected_changed
    integer :: expected_checksum
    integer :: global_recv
    integer :: global_recv_checksum
    integer :: global_send
    integer :: global_send_checksum
    integer :: ierr
    integer :: local_recv_checksum
    integer :: local_send_checksum
    integer :: pos
    integer :: r
    integer :: source
    logical :: print_local
    logical, allocatable :: seen(:)

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. allocated(manifest%send_count) .or. &
         .not. allocated(manifest%recv_count) .or. &
         .not. allocated(manifest%send_displ) .or. &
         .not. allocated(manifest%recv_displ) .or. &
         .not. allocated(manifest%send_block) .or. &
         .not. allocated(manifest%recv_block)) then
       call fail("migration manifest is incomplete")
    end if

    if (size(manifest%send_count) /= n_process .or. &
         size(manifest%recv_count) /= n_process .or. &
         size(manifest%send_displ) /= n_process .or. &
         size(manifest%recv_displ) /= n_process) then
       call fail("migration manifest has an invalid rank dimension")
    end if

    if (manifest%n_send /= sum(manifest%send_count) .or. &
         manifest%n_recv /= sum(manifest%recv_count)) then
       call fail("migration manifest has an invalid local total")
    end if

    allocate(seen(size(block_catalog)))
    seen = .false.

    do r = 0, n_process-1
       do pos = manifest%recv_displ(r+1) + 1, &
            manifest%recv_displ(r+1) + manifest%recv_count(r+1)

          b = manifest%recv_block(pos)

          if (b < 1 .or. b > size(block_catalog)) then
             call fail("received block-catalog index is invalid")
          end if

          if (seen(b)) then
             call fail("duplicate block in received migration manifest")
          end if
          seen(b) = .true.

          source      = source_rank(b)
          destination = block_catalog(b)%owner

          if (source /= r) then
             call fail("received block has the wrong source rank")
          end if

          if (destination /= rank) then
             call fail("received block has the wrong destination rank")
          end if

          if (source == destination) then
             call fail("retained block appears in migration manifest")
          end if
       end do
    end do

    deallocate(seen)

    do pos = 1, manifest%n_send
       b = manifest%send_block(pos)

       if (b < 1 .or. b > size(block_catalog)) then
          call fail("sent block-catalog index is invalid")
       end if

       if (source_rank(b) /= rank) then
          call fail("sent block does not originate on this rank")
       end if

       if (block_catalog(b)%owner == rank) then
          call fail("retained block appears in send manifest")
       end if
    end do

    local_send_checksum = 0
    local_recv_checksum = 0

    if (manifest%n_send > 0) then
       local_send_checksum = &
            sum(manifest%send_block(1:manifest%n_send))
    end if

    if (manifest%n_recv > 0) then
       local_recv_checksum = &
            sum(manifest%recv_block(1:manifest%n_recv))
    end if

    call MPI_Allreduce( &
         manifest%n_send, global_send, 1, MPI_INTEGER, MPI_SUM, &
         comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce sent migration count")

    call MPI_Allreduce( &
         manifest%n_recv, global_recv, 1, MPI_INTEGER, MPI_SUM, &
         comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce received migration count")

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce sent manifest checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce received manifest checksum")

    expected_changed  = 0
    expected_checksum = 0

    do b = 1, size(block_catalog)
       if (source_rank(b) == block_catalog(b)%owner) cycle
       expected_changed  = expected_changed + 1
       expected_checksum = expected_checksum + b
    end do

    if (global_send /= expected_changed .or. &
         global_recv /= expected_changed) then
       call fail("global migration count does not match block_catalog")
    end if

    if (global_send_checksum /= expected_checksum .or. &
         global_recv_checksum /= expected_checksum) then
       call fail("global migration checksum does not match block_catalog")
    end if

    manifest%validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration manifest for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  blocks sent     = ", manifest%n_send
       write(6,'(a,i0)') &
            "  blocks received = ", manifest%n_recv
       write(6,'(a)') &
            "  source/destination checks passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global migrating blocks verified = ", expected_changed
       write(6,'(a,/)') &
            "Block migration manifest exchange passed"
    end if

  end subroutine check_block_migration_manifest


  subroutine exchange_block_migration_sizes (manifest, send_nbyte, verbose)
    ! Exchange only the packed byte count associated with each entry in
    ! the already validated migration manifest. Block payloads are not
    ! communicated and persistent block storage is not modified.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    integer, intent(in)                            :: send_nbyte(:)
    logical, optional, intent(in)                  :: verbose

    integer :: global_max_recv
    integer :: global_max_send
    integer :: ierr
    integer :: pos
    integer(int64) :: global_recv_checksum
    integer(int64) :: global_recv_nbyte
    integer(int64) :: global_send_checksum
    integer(int64) :: global_send_nbyte
    integer(int64) :: local_recv_checksum
    integer(int64) :: local_send_checksum
    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. manifest%validated) then
       call fail("migration manifest must be validated before sizes")
    end if

    if (size(send_nbyte) /= manifest%n_send) then
       call fail("outgoing byte-count array has the wrong extent")
    end if

    if (manifest%n_send > 0) then
       if (any(send_nbyte <= 0)) then
          call fail("outgoing packed block has a nonpositive byte count")
       end if
    end if

    if (allocated(manifest%send_nbyte)) then
       deallocate(manifest%send_nbyte)
    end if

    if (allocated(manifest%recv_nbyte)) then
       deallocate(manifest%recv_nbyte)
    end if

    allocate(manifest%send_nbyte(max(1,manifest%n_send)))
    allocate(manifest%recv_nbyte(max(1,manifest%n_recv)))

    manifest%send_nbyte = 0
    manifest%recv_nbyte = 0

    if (manifest%n_send > 0) then
       manifest%send_nbyte(1:manifest%n_send) = send_nbyte
    end if

    call MPI_Alltoallv( &
         manifest%send_nbyte, manifest%send_count, &
         manifest%send_displ, MPI_INTEGER, &
         manifest%recv_nbyte, manifest%recv_count, &
         manifest%recv_displ, MPI_INTEGER, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv packed block sizes")

    if (manifest%n_recv > 0) then
       if (any(manifest%recv_nbyte(1:manifest%n_recv) <= 0)) then
          call fail("received packed block has a nonpositive byte count")
       end if
    end if

    manifest%total_send_nbyte = 0_int64
    manifest%total_recv_nbyte = 0_int64
    manifest%max_send_nbyte   = 0
    manifest%max_recv_nbyte   = 0

    if (manifest%n_send > 0) then
       manifest%total_send_nbyte = &
            sum(int(manifest%send_nbyte(1:manifest%n_send),int64))
       manifest%max_send_nbyte = &
            maxval(manifest%send_nbyte(1:manifest%n_send))
    end if

    if (manifest%n_recv > 0) then
       manifest%total_recv_nbyte = &
            sum(int(manifest%recv_nbyte(1:manifest%n_recv),int64))
       manifest%max_recv_nbyte = &
            maxval(manifest%recv_nbyte(1:manifest%n_recv))
    end if

    local_send_checksum = 0_int64
    local_recv_checksum = 0_int64

    do pos = 1, manifest%n_send
       local_send_checksum = local_send_checksum + &
            int(manifest%send_block(pos),int64) * &
            int(manifest%send_nbyte(pos),int64)
    end do

    do pos = 1, manifest%n_recv
       local_recv_checksum = local_recv_checksum + &
            int(manifest%recv_block(pos),int64) * &
            int(manifest%recv_nbyte(pos),int64)
    end do

    call MPI_Allreduce( &
         manifest%total_send_nbyte, global_send_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing packed bytes")

    call MPI_Allreduce( &
         manifest%total_recv_nbyte, global_recv_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming packed bytes")

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing size checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming size checksum")

    call MPI_Allreduce( &
         manifest%max_send_nbyte, global_max_send, 1, &
         MPI_INTEGER, MPI_MAX, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce maximum outgoing block size")

    call MPI_Allreduce( &
         manifest%max_recv_nbyte, global_max_recv, 1, &
         MPI_INTEGER, MPI_MAX, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce maximum incoming block size")

    if (global_send_nbyte /= global_recv_nbyte) then
       call fail("global outgoing and incoming byte totals differ")
    end if

    if (global_send_checksum /= global_recv_checksum) then
       call fail("global block/size checksums differ")
    end if

    if (global_max_send /= global_max_recv) then
       call fail("global maximum outgoing and incoming sizes differ")
    end if

    manifest%sizes_validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration sizes for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  outgoing blocks       = ", manifest%n_send
       write(6,'(a,i0)') &
            "  outgoing packed bytes = ", manifest%total_send_nbyte
       write(6,'(a,i0)') &
            "  incoming blocks       = ", manifest%n_recv
       write(6,'(a,i0)') &
            "  incoming packed bytes = ", manifest%total_recv_nbyte
       write(6,'(a)') &
            "  packed-size routing checks passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global packed migration bytes verified = ", &
            global_send_nbyte
       write(6,'(a,i0)') &
            "Maximum packed migrating block bytes   = ", &
            global_max_send
       write(6,'(a,/)') &
            "Block migration size exchange passed"
    end if

  end subroutine exchange_block_migration_sizes


  subroutine exchange_block_migration_payloads ( &
       manifest, send_payload, verbose)
    ! Exchange the serialized payloads using the already validated
    ! per-block byte counts. The received byte stream remains separate
    ! from all source and installed block storage.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest
    integer(int8), intent(in)                      :: send_payload(:)
    logical, optional, intent(in)                  :: verbose

    integer :: first
    integer :: ierr
    integer :: last
    integer :: n_recv_byte
    integer :: n_send_byte
    integer :: r
    integer(int64) :: global_recv_checksum
    integer(int64) :: global_send_checksum
    integer(int64) :: global_send_nbyte
    integer(int64) :: local_recv_checksum
    integer(int64) :: local_send_checksum
    integer(int64) :: rank_nbyte
    logical :: print_local

    print_local = .true.
    if (present(verbose)) print_local = verbose

    if (.not. manifest%validated .or. &
         .not. manifest%sizes_validated) then
       call fail("manifest and sizes must precede payload exchange")
    end if

    if (manifest%total_send_nbyte > int(huge(0),int64) .or. &
         manifest%total_recv_nbyte > int(huge(0),int64)) then
       call fail("packed byte total exceeds MPI default-count range")
    end if

    n_send_byte = int(manifest%total_send_nbyte)
    n_recv_byte = int(manifest%total_recv_nbyte)

    if (n_send_byte > 0) then
       if (size(send_payload) /= n_send_byte) then
          call fail("outgoing payload buffer has the wrong extent")
       end if
    else
       if (size(send_payload) < 1) then
          call fail("zero-byte sender requires a one-element MPI buffer")
       end if
    end if

    if (allocated(manifest%send_byte_count)) then
       deallocate(manifest%send_byte_count)
    end if
    if (allocated(manifest%recv_byte_count)) then
       deallocate(manifest%recv_byte_count)
    end if
    if (allocated(manifest%send_byte_displ)) then
       deallocate(manifest%send_byte_displ)
    end if
    if (allocated(manifest%recv_byte_displ)) then
       deallocate(manifest%recv_byte_displ)
    end if
    if (allocated(manifest%recv_payload)) then
       deallocate(manifest%recv_payload)
    end if

    allocate(manifest%send_byte_count(n_process))
    allocate(manifest%recv_byte_count(n_process))
    allocate(manifest%send_byte_displ(n_process))
    allocate(manifest%recv_byte_displ(n_process))

    manifest%send_byte_count = 0
    manifest%recv_byte_count = 0
    manifest%send_byte_displ = 0
    manifest%recv_byte_displ = 0
    manifest%payload_validated = .false.

    do r = 1, n_process

       first = manifest%send_displ(r) + 1
       last  = manifest%send_displ(r) + manifest%send_count(r)

       rank_nbyte = 0_int64
       if (last >= first) then
          rank_nbyte = &
               sum(int(manifest%send_nbyte(first:last),int64))
       end if

       if (rank_nbyte > int(huge(0),int64)) then
          call fail("per-rank outgoing payload exceeds MPI count range")
       end if

       manifest%send_byte_count(r) = int(rank_nbyte)

       first = manifest%recv_displ(r) + 1
       last  = manifest%recv_displ(r) + manifest%recv_count(r)

       rank_nbyte = 0_int64
       if (last >= first) then
          rank_nbyte = &
               sum(int(manifest%recv_nbyte(first:last),int64))
       end if

       if (rank_nbyte > int(huge(0),int64)) then
          call fail("per-rank incoming payload exceeds MPI count range")
       end if

       manifest%recv_byte_count(r) = int(rank_nbyte)

    end do

    do r = 2, n_process
       manifest%send_byte_displ(r) = &
            manifest%send_byte_displ(r-1) + &
            manifest%send_byte_count(r-1)
       manifest%recv_byte_displ(r) = &
            manifest%recv_byte_displ(r-1) + &
            manifest%recv_byte_count(r-1)
    end do

    if (sum(int(manifest%send_byte_count,int64)) /= &
         manifest%total_send_nbyte) then
       call fail("outgoing per-rank byte counts do not sum correctly")
    end if

    if (sum(int(manifest%recv_byte_count,int64)) /= &
         manifest%total_recv_nbyte) then
       call fail("incoming per-rank byte counts do not sum correctly")
    end if

    allocate(manifest%recv_payload(max(1,n_recv_byte)))
    manifest%recv_payload = 0_int8

    call MPI_Alltoallv( &
         send_payload, manifest%send_byte_count, &
         manifest%send_byte_displ, MPI_BYTE, &
         manifest%recv_payload, manifest%recv_byte_count, &
         manifest%recv_byte_displ, MPI_BYTE, comm, ierr)
    call check_mpi(ierr,"MPI_Alltoallv packed block payloads")

    local_send_checksum = 0_int64
    local_recv_checksum = 0_int64

    if (n_send_byte > 0) then
       local_send_checksum = &
            sum(int(send_payload(1:n_send_byte),int64))
    end if

    if (n_recv_byte > 0) then
       local_recv_checksum = &
            sum(int(manifest%recv_payload(1:n_recv_byte),int64))
    end if

    call MPI_Allreduce( &
         local_send_checksum, global_send_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce outgoing payload checksum")

    call MPI_Allreduce( &
         local_recv_checksum, global_recv_checksum, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce incoming payload checksum")

    call MPI_Allreduce( &
         manifest%total_send_nbyte, global_send_nbyte, 1, &
         MPI_INTEGER8, MPI_SUM, comm, ierr)
    call check_mpi(ierr,"MPI_Allreduce transferred payload bytes")

    if (global_send_checksum /= global_recv_checksum) then
       call fail("global payload byte checksums differ")
    end if

    manifest%payload_validated = .true.

    if (print_local) then
       write(6,'(/,a,i0,a)') &
            "Block migration payload transport for rank ", rank, ":"
       write(6,'(a,i0)') &
            "  outgoing payload bytes = ", &
            manifest%total_send_nbyte
       write(6,'(a,i0)') &
            "  incoming payload bytes = ", &
            manifest%total_recv_nbyte
       write(6,'(a)') &
            "  payload transport checksum passed"
    end if

    if (rank == 0) then
       write(6,'(/,a,i0)') &
            "Global packed payload bytes transferred = ", &
            global_send_nbyte
       write(6,'(a,/)') &
            "Block migration payload transport passed"
    end if

  end subroutine exchange_block_migration_payloads


  subroutine clear_block_migration_manifest (manifest)
    ! Release all count, routing, size and payload staging storage after
    ! the received blocks have been installed and independently checked.

    implicit none

    type(Block_Migration_Manifest), intent(inout) :: manifest

    if (allocated(manifest%send_count)) then
       deallocate(manifest%send_count)
    end if
    if (allocated(manifest%recv_count)) then
       deallocate(manifest%recv_count)
    end if
    if (allocated(manifest%send_displ)) then
       deallocate(manifest%send_displ)
    end if
    if (allocated(manifest%recv_displ)) then
       deallocate(manifest%recv_displ)
    end if
    if (allocated(manifest%send_block)) then
       deallocate(manifest%send_block)
    end if
    if (allocated(manifest%recv_block)) then
       deallocate(manifest%recv_block)
    end if
    if (allocated(manifest%send_nbyte)) then
       deallocate(manifest%send_nbyte)
    end if
    if (allocated(manifest%recv_nbyte)) then
       deallocate(manifest%recv_nbyte)
    end if
    if (allocated(manifest%send_byte_count)) then
       deallocate(manifest%send_byte_count)
    end if
    if (allocated(manifest%recv_byte_count)) then
       deallocate(manifest%recv_byte_count)
    end if
    if (allocated(manifest%send_byte_displ)) then
       deallocate(manifest%send_byte_displ)
    end if
    if (allocated(manifest%recv_byte_displ)) then
       deallocate(manifest%recv_byte_displ)
    end if
    if (allocated(manifest%recv_payload)) then
       deallocate(manifest%recv_payload)
    end if

    manifest%n_send = 0
    manifest%n_recv = 0
    manifest%total_send_nbyte = 0_int64
    manifest%total_recv_nbyte = 0_int64
    manifest%max_send_nbyte = 0
    manifest%max_recv_nbyte = 0
    manifest%validated = .false.
    manifest%sizes_validated = .false.
    manifest%payload_validated = .false.

  end subroutine clear_block_migration_manifest


  integer function source_rank (b) result(source)

    implicit none

    integer, intent(in) :: b
    integer             :: root_domain

    if (b < 1 .or. b > size(block_catalog)) then
       call fail("block-catalog index is invalid")
    end if

    root_domain = block_catalog(b)%root_domain

    if (root_domain < 0 .or. root_domain >= size(owner)) then
       call fail("block root-domain index is invalid")
    end if

    source = owner(root_domain+1)
    call check_rank(source,"source domain owner")

  end function source_rank


  subroutine check_rank (r, description)

    implicit none

    integer, intent(in)          :: r
    character(*), intent(in)     :: description

    if (r < 0 .or. r >= n_process) then
       call fail(trim(description)//" is outside the MPI rank range")
    end if

  end subroutine check_rank


  subroutine check_mpi (ierr, operation)

    implicit none

    integer, intent(in)      :: ierr
    character(*), intent(in) :: operation

    if (ierr /= MPI_SUCCESS) then
       write(error_unit,'(a,i0,2a,i0)') &
            "Rank ", rank, ": ", trim(operation)// &
            " failed with MPI error ", ierr
       call abort_run
    end if

  end subroutine check_mpi


  subroutine fail (message)

    implicit none

    character(*), intent(in) :: message

    write(error_unit,'(a,i0,2a)') &
         "Rank ", rank, ": parallel_block_mpi_mod: ", trim(message)
    call abort_run

  end subroutine fail


end module parallel_block_mpi_mod
