module parallel_block_mod

  use, intrinsic :: iso_fortran_env, only : int8, int64

  use kind_mod,   only : dp
  use shared_mod, only : Coord, EDGE, MULT, N_BDRY, S_MASS, S_TEMP, &
       S_VELO, c_p, compressible, grav_accel, kappa, p_0, p_top, &
       scalars, vert_diffuse, zlevels, zmin, zmax
  use patch_mod,  only : Patch, PATCH_SIZE

  implicit none

  private

  integer, parameter, public :: STORE_NONE  = 0
  integer, parameter, public :: STORE_PATCH = 1
  integer, parameter, public :: STORE_BDRY  = 2
  integer, parameter, public :: STORE_GHOST = 3

  integer, parameter, public :: NGB_INTERNAL = 1
  integer, parameter, public :: NGB_BLOCK    = 2
  integer, parameter, public :: NGB_DOMAIN   = 3
  integer, parameter, public :: NGB_ADAPT    = 4
  integer, parameter, public :: NGB_OTHER    = 5

  integer, parameter, public :: BLOCK_PAYLOAD_SOL = 1
  integer, parameter, public :: BLOCK_PAYLOAD_WAV_COEFF = 2

  integer, parameter :: BLOCK_PACK_MAGIC = &
       int(z'54424C4B')
  integer, parameter :: BLOCK_PACK_VERSION = 9
  integer, parameter :: BLOCK_PACK_HEADER_SIZE = 44


  type, public :: Block_Bdry_Storage
     integer :: source_bdry = -1
     integer :: elts_start  = -1
     integer :: dims(2)     = 0
     integer :: n_node      = 0
     integer :: local_start = -1
  end type Block_Bdry_Storage


  type, public :: Block_Bdry_Link
     integer :: patch       = -1
     integer :: side        = -1
     integer :: class       = -1

     integer :: root_domain = -1

     integer :: neigh_patch = -1
     integer :: source_bdry = -1

     integer :: elts_start  = -1
     integer :: bdry_side   = 0
     integer :: bdry_neigh  = 0

     integer :: dims(2)     = 0
     integer :: n_node      = 0

     integer :: storage_id  = -1

     integer :: source_block    = -1
     integer :: source_block_id = -1
     integer :: source_owner    = -1
     integer :: ghost_id        = -1
  end type Block_Bdry_Link


  type, public :: Block_Ghost_Storage
     integer :: source_patch = -1
     integer :: elts_start   = -1
     integer :: local_start  = -1
     integer :: n_node       = 0

     integer :: source_domain   = -1
     integer :: source_block    = -1
     integer :: source_block_id = -1
     integer :: source_owner    = -1
     integer :: source_local_patch = -1
  end type Block_Ghost_Storage


  type, public :: Block_Stencil_Address
     integer :: storage = STORE_NONE
     integer :: id      = -1
     integer :: offset  = 0
     integer :: dims(2) = 0
  end type Block_Stencil_Address


  type, public :: Block_Data
     integer :: id          = -1
     integer :: root_domain = -1
     integer :: root_patch  = -1
     integer :: level       = -1

     integer :: scalar_variable = -1
     integer :: n_scalar_variable = 0
     integer :: vector_variable = -1
     integer :: field_level     = -1
     integer :: n_field_level   = 0
     integer :: scalar_mult     = 0
     integer :: vector_mult     = 0
     integer :: tke_level       = 1
     integer :: n_tke_level     = 0

     type(Patch), allocatable :: patch(:)
     type(Coord), allocatable :: node(:)
     type(Coord), allocatable :: bdry_node(:)

     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
     real(dp), allocatable :: wavelet_scalar(:)
     real(dp), allocatable :: wavelet_vector(:)
     real(dp), allocatable :: scalar_mean(:)
     real(dp), allocatable :: vector_mean(:)
     real(dp), allocatable :: tke(:)
     real(dp), allocatable :: wavelet_tke(:)
     real(dp), allocatable :: topography(:)
     real(dp), allocatable :: bdry_scalar(:)
     real(dp), allocatable :: bdry_vector(:)
     real(dp), allocatable :: bdry_wavelet_scalar(:)
     real(dp), allocatable :: bdry_wavelet_vector(:)
     real(dp), allocatable :: bdry_scalar_mean(:)
     real(dp), allocatable :: bdry_vector_mean(:)
     real(dp), allocatable :: bdry_tke(:)
     real(dp), allocatable :: bdry_wavelet_tke(:)
     real(dp), allocatable :: bdry_topography(:)

     integer, allocatable :: neigh_class(:,:)

     type(Block_Bdry_Link), allocatable :: block_bdry(:)
     type(Block_Bdry_Storage), allocatable :: bdry_storage(:)
     type(Block_Stencil_Address), allocatable :: stencil(:,:)
     type(Block_Ghost_Storage), allocatable :: ghost_storage(:)

     type(Coord), allocatable :: ghost_node(:)

     real(dp), allocatable :: ghost_scalar(:)
     real(dp), allocatable :: ghost_vector(:)
     real(dp), allocatable :: ghost_wavelet_scalar(:)
     real(dp), allocatable :: ghost_wavelet_vector(:)
     real(dp), allocatable :: ghost_scalar_mean(:)
     real(dp), allocatable :: ghost_vector_mean(:)
     real(dp), allocatable :: ghost_tke(:)
     real(dp), allocatable :: ghost_wavelet_tke(:)
     real(dp), allocatable :: ghost_topography(:)
  end type Block_Data


  type :: Block_Hydrostatic_Storage
     integer :: catalog_index = 0
     logical :: ready = .false.
     integer(int64) :: refreshes = 0_int64
     real(dp), allocatable :: surface_pressure(:)
     real(dp), allocatable :: dynamic_exner(:)
     real(dp), allocatable :: air_temperature(:)
  end type Block_Hydrostatic_Storage


  type :: Block_Tendency_Storage
     integer :: catalog_index = 0
     logical :: ready = .false.
     integer(int64) :: generation = 0_int64
     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
  end type Block_Tendency_Storage


  type :: Block_Tendency_Trial_Storage
     integer :: catalog_index = 0
     logical :: active = .false.
     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
  end type Block_Tendency_Trial_Storage


  type :: Block_Tendency_Accumulator_Storage
     integer :: catalog_index = 0
     real(dp), allocatable :: scalar(:)
     real(dp), allocatable :: vector(:)
  end type Block_Tendency_Accumulator_Storage


  abstract interface
     subroutine Local_Block_Field_Consumer (catalog_index,block,context)
       import :: Block_Data

       integer, intent(in) :: catalog_index
       type(Block_Data), intent(in) :: block
       class(*), intent(inout) :: context
     end subroutine Local_Block_Field_Consumer

     subroutine Local_Block_Tendency_Kernel ( &
          catalog_index,block,scalar_tendency,vector_tendency,context)
       import :: Block_Data, dp

       integer, intent(in) :: catalog_index
       type(Block_Data), intent(in) :: block
       real(dp), intent(inout) :: scalar_tendency(:)
       real(dp), intent(inout) :: vector_tendency(:)
       class(*), intent(inout) :: context
     end subroutine Local_Block_Tendency_Kernel

     subroutine Local_Block_Tendency_Consumer ( &
          catalog_index,scalar_tendency,vector_tendency,context)
       import :: dp

       integer, intent(in) :: catalog_index
       real(dp), intent(in) :: scalar_tendency(:)
       real(dp), intent(in) :: vector_tendency(:)
       class(*), intent(inout) :: context
     end subroutine Local_Block_Tendency_Consumer

     subroutine Local_Block_Hydrostatic_Consumer ( &
          catalog_index,n_patch,surface_pressure,dynamic_exner, &
          air_temperature,context)
       import :: dp

       integer, intent(in) :: catalog_index
       integer, intent(in) :: n_patch
       real(dp), intent(in) :: surface_pressure(:)
       real(dp), intent(in) :: dynamic_exner(:)
       real(dp), intent(in) :: air_temperature(:)
       class(*), intent(inout) :: context
     end subroutine Local_Block_Hydrostatic_Consumer
  end interface


  type(Block_Data), allocatable, public :: block_source(:)
  type(Block_Data), allocatable, public :: block_received(:)
  type(Block_Data), allocatable :: block_local(:)
  type(Block_Hydrostatic_Storage), allocatable :: block_hydrostatic(:)
  type(Block_Tendency_Storage), allocatable :: block_tendency(:)
  type(Block_Tendency_Trial_Storage), allocatable :: block_tendency_trial(:)
  type(Block_Tendency_Accumulator_Storage), allocatable :: &
       block_tendency_accumulator(:)

  integer, allocatable, public :: block_source_catalog_index(:)
  integer, allocatable, public :: block_retained_source_index(:)
  integer, allocatable, public :: block_migrating_source_index(:)
  integer, allocatable, public :: block_received_catalog_index(:)
  integer, allocatable :: block_local_catalog_index(:)
  integer, allocatable :: block_catalog_local_index(:)
  integer, allocatable :: block_tendency_import_patch_count(:)

  logical :: block_store_ready = .false.
  logical :: block_hydrostatic_ready = .false.
  integer(int64) :: block_hydrostatic_refreshes = 0_int64
  logical :: block_tendency_ready = .false.
  logical :: block_tendency_import_active = .false.
  logical :: block_tendency_trial_active = .false.
  logical :: block_tendency_commit_checkpoint_ready = .false.
  logical :: block_tendency_accumulator_ready = .false.
  integer(int64) :: block_tendency_executions = 0_int64
  integer(int64) :: block_tendency_allocations = 0_int64
  integer(int64) :: block_tendency_import_allocations = 0_int64
  integer(int64) :: block_tendency_accumulator_allocations = 0_int64
  integer(int64) :: block_tendency_accumulator_stages = 0_int64

  public :: packed_block_nbyte
  public :: pack_block
  public :: unpack_block
  public :: check_block_storage
  public :: clear_block_staging
  public :: clear_local_blocks
  public :: local_block_store_ready
  public :: n_local_blocks
  public :: local_block_catalog
  public :: catalog_local_block
  public :: get_local_block_identity
  public :: get_block_field_layout
  public :: get_block_turbulence_layout
  public :: get_local_block_field_layout
  public :: get_local_block_turbulence_layout
  public :: check_local_block_storage
  public :: local_block_field_statistics
  public :: local_block_wavelet_statistics
  public :: local_block_mean_field_statistics
  public :: local_block_turbulence_statistics
  public :: local_block_topography_statistics
  public :: source_block_scalar_stencil_statistics
  public :: local_block_scalar_stencil_statistics
  public :: source_block_vector_stencil_statistics
  public :: local_block_vector_stencil_statistics
  public :: source_block_boundary_route_statistics
  public :: local_block_boundary_route_statistics
  public :: source_block_ghost_source_statistics
  public :: local_block_ghost_source_statistics
  public :: validate_local_block_ghost_sources
  public :: get_local_block_ghost_requests
  public :: local_block_patch_count
  public :: local_block_boundary_count
  public :: get_local_block_boundary_source
  public :: local_block_ghost_count
  public :: local_block_scalar_patch_nvalue
  public :: local_block_scalar_family_patch_nvalue
  public :: local_block_scalar_boundary_nvalue
  public :: local_block_scalar_family_boundary_nvalue
  public :: get_local_block_scalar_patch_values
  public :: get_local_block_scalar_patch_family_values
  public :: set_local_block_scalar_patch_values
  public :: set_local_block_scalar_patch_family_values
  public :: fill_local_block_scalar_patch_values
  public :: fill_local_block_scalar_patch_family_values
  public :: get_local_block_scalar_boundary_values
  public :: get_local_block_scalar_boundary_family_values
  public :: set_local_block_scalar_boundary_family_values
  public :: fill_local_block_scalar_boundary_values
  public :: fill_local_block_scalar_boundary_family_values
  public :: get_local_block_scalar_ghost_values
  public :: get_local_block_scalar_ghost_family_values
  public :: set_local_block_scalar_ghost_values
  public :: set_local_block_scalar_ghost_family_values
  public :: fill_local_block_scalar_ghost_values
  public :: fill_local_block_scalar_ghost_family_values
  public :: local_block_vector_patch_nvalue
  public :: local_block_vector_family_patch_nvalue
  public :: local_block_vector_boundary_nvalue
  public :: local_block_vector_family_boundary_nvalue
  public :: get_local_block_vector_patch_values
  public :: get_local_block_vector_patch_family_values
  public :: set_local_block_vector_patch_values
  public :: set_local_block_vector_patch_family_values
  public :: fill_local_block_vector_patch_values
  public :: fill_local_block_vector_patch_family_values
  public :: get_local_block_vector_boundary_values
  public :: get_local_block_vector_boundary_family_values
  public :: set_local_block_vector_boundary_family_values
  public :: fill_local_block_vector_boundary_values
  public :: fill_local_block_vector_boundary_family_values
  public :: get_local_block_vector_ghost_values
  public :: get_local_block_vector_ghost_family_values
  public :: set_local_block_vector_ghost_values
  public :: set_local_block_vector_ghost_family_values
  public :: fill_local_block_vector_ghost_values
  public :: fill_local_block_vector_ghost_family_values
  public :: compute_local_block_hydrostatic_patch
  public :: refresh_local_block_hydrostatic_state
  public :: ensure_local_block_hydrostatic_state
  public :: local_block_hydrostatic_state_ready
  public :: local_block_hydrostatic_refresh_count
  public :: local_block_hydrostatic_block_refresh_count
  public :: local_block_hydrostatic_surface_nvalue
  public :: local_block_hydrostatic_column_nvalue
  public :: get_local_block_hydrostatic_patch_values
  public :: get_local_block_hydrostatic_values
  public :: Local_Block_Field_Consumer
  public :: apply_local_block_field_consumer
  public :: Local_Block_Tendency_Kernel
  public :: apply_local_block_tendency_kernel
  public :: Local_Block_Tendency_Consumer
  public :: apply_local_block_tendency_consumer
  public :: local_block_tendency_state_ready
  public :: discard_local_block_tendency_output
  public :: invalidate_local_block_tendency_products
  public :: prepare_local_block_tendency_workspace
  public :: local_block_tendency_execution_count
  public :: local_block_tendency_allocation_count
  public :: local_block_tendency_statistics
  public :: begin_local_block_tendency_import
  public :: set_local_block_tendency_patch_values
  public :: finish_local_block_tendency_import
  public :: get_local_block_tendency_patch_values
  public :: assert_local_block_tendency_patch_values
  public :: local_block_tendency_import_is_active
  public :: local_block_tendency_import_allocation_count
  public :: reset_local_block_tendency_accumulator
  public :: accumulate_local_block_tendency
  public :: begin_local_block_accumulated_tendency_trial
  public :: local_block_tendency_accumulator_state_ready
  public :: local_block_tendency_accumulator_allocation_count
  public :: local_block_tendency_accumulator_statistics
  public :: begin_local_block_tendency_trial
  public :: commit_local_block_tendency_trial
  public :: finalize_local_block_tendency_commit
  public :: restore_local_block_tendency_commit
  public :: local_block_tendency_commit_checkpoint_is_ready
  public :: local_block_tendency_commit_checkpoint_statistics
  public :: rollback_local_block_tendency_trial
  public :: local_block_tendency_trial_is_active
  public :: Local_Block_Hydrostatic_Consumer
  public :: apply_local_block_hydrostatic_consumer
  public :: local_block_hydrostatic_statistics
  public :: install_local_blocks

contains


subroutine get_block_field_layout ( &
     scalar_variable,n_scalar_variable,vector_variable,field_level, &
     n_field_level,scalar_mult,vector_mult)
  ! Return the field layout represented by this Block_Data format.
  ! Centralizing the descriptor prevents builders and consumers from
  ! independently hard-coding the selected variables and sol levels.

  implicit none

  integer, intent(out) :: scalar_variable
  integer, intent(out) :: n_scalar_variable
  integer, intent(out) :: vector_variable
  integer, intent(out) :: field_level
  integer, intent(out) :: n_field_level
  integer, intent(out) :: scalar_mult
  integer, intent(out) :: vector_mult

  scalar_variable = scalars(1)
  n_scalar_variable = scalars(2) - scalars(1) + 1
  vector_variable = S_VELO
  field_level     = zmin
  n_field_level   = zmax - zmin + 1
  scalar_mult     = MULT(scalar_variable)
  vector_mult     = MULT(vector_variable)

end subroutine get_block_field_layout


subroutine get_block_turbulence_layout (tke_level,n_tke_level)
  ! Return the optional node-based TKE layout represented by Block_Data.

  implicit none

  integer, intent(out) :: tke_level
  integer, intent(out) :: n_tke_level

  tke_level = 1
  n_tke_level = 0
  if (vert_diffuse) n_tke_level = zlevels

end subroutine get_block_turbulence_layout


subroutine check_block_storage (block,check_serialization)
  ! Check the allocation and extents of one block without referring to
  ! ownership, catalogues or MPI state. Optionally verify an exact
  ! pack/unpack/pack round trip.

  implicit none

  type(Block_Data), intent(in) :: block
  logical, optional, intent(in) :: check_serialization

  integer :: n_bdry_node
  integer :: n_ghost_node
  integer :: n_node
  integer :: i
  integer :: scalar_variable
  integer :: n_scalar_variable
  integer :: vector_variable
  integer :: field_level
  integer :: n_field_level
  integer :: scalar_mult
  integer :: vector_mult
  integer :: tke_level
  integer :: n_tke_level

  integer(int8), allocatable :: buffer_copy(:)
  integer(int8), allocatable :: buffer_source(:)

  type(Block_Data) :: block_copy

  logical :: serialize

  serialize = .false.
  if (present(check_serialization)) serialize = check_serialization

  call get_block_field_layout( &
       scalar_variable,n_scalar_variable,vector_variable,field_level, &
       n_field_level,scalar_mult,vector_mult)

  call get_block_turbulence_layout(tke_level,n_tke_level)

  if (block%scalar_variable /= scalar_variable .or. &
       block%n_scalar_variable /= n_scalar_variable .or. &
       block%vector_variable /= vector_variable .or. &
       block%field_level /= field_level .or. &
       block%n_field_level /= n_field_level .or. &
       block%scalar_mult /= scalar_mult .or. &
       block%vector_mult /= vector_mult) then

     error stop "check_block_storage: field layout mismatch"

  end if

  if (block%tke_level /= tke_level .or. &
       block%n_tke_level /= n_tke_level) then
     error stop "check_block_storage: turbulence layout mismatch"
  end if

  if (.not. allocated(block%patch) .or. &
       .not. allocated(block%node) .or. &
       .not. allocated(block%scalar) .or. &
       .not. allocated(block%vector) .or. &
       .not. allocated(block%wavelet_scalar) .or. &
       .not. allocated(block%wavelet_vector) .or. &
       .not. allocated(block%scalar_mean) .or. &
       .not. allocated(block%vector_mean) .or. &
       .not. allocated(block%tke) .or. &
       .not. allocated(block%wavelet_tke) .or. &
       .not. allocated(block%topography) .or. &
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%bdry_wavelet_scalar) .or. &
       .not. allocated(block%bdry_wavelet_vector) .or. &
       .not. allocated(block%bdry_scalar_mean) .or. &
       .not. allocated(block%bdry_vector_mean) .or. &
       .not. allocated(block%bdry_tke) .or. &
       .not. allocated(block%bdry_wavelet_tke) .or. &
       .not. allocated(block%bdry_topography) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector) .or. &
       .not. allocated(block%ghost_wavelet_scalar) .or. &
       .not. allocated(block%ghost_wavelet_vector) .or. &
       .not. allocated(block%ghost_scalar_mean) .or. &
       .not. allocated(block%ghost_vector_mean) .or. &
       .not. allocated(block%ghost_tke) .or. &
       .not. allocated(block%ghost_wavelet_tke) .or. &
       .not. allocated(block%ghost_topography)) then

     error stop "check_block_storage: unallocated component"

  end if

  n_node = size(block%patch) * PATCH_SIZE**2

  if (size(block%node) /= n_node .or. &
       size(block%scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_node .or. &
       size(block%vector) /= &
       block%n_field_level*block%vector_mult*n_node .or. &
       size(block%wavelet_scalar) /= size(block%scalar) .or. &
       size(block%wavelet_vector) /= size(block%vector) .or. &
       size(block%scalar_mean) /= size(block%scalar) .or. &
       size(block%vector_mean) /= size(block%vector) .or. &
       size(block%tke) /= block%n_tke_level*n_node .or. &
       size(block%wavelet_tke) /= size(block%tke) .or. &
       size(block%topography) /= n_node) then

     error stop "check_block_storage: interior extent mismatch"

  end if

  if (size(block%neigh_class,1) /= N_BDRY .or. &
       size(block%neigh_class,2) /= size(block%patch) .or. &
       size(block%stencil,1) /= N_BDRY .or. &
       size(block%stencil,2) /= size(block%patch)) then

     error stop "check_block_storage: topology extent mismatch"

  end if

  n_bdry_node = sum(block%bdry_storage%n_node)

  if (size(block%bdry_node) /= n_bdry_node .or. &
       size(block%bdry_scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_bdry_node .or. &
       size(block%bdry_vector) /= &
       block%n_field_level*block%vector_mult*n_bdry_node .or. &
       size(block%bdry_wavelet_scalar) /= size(block%bdry_scalar) .or. &
       size(block%bdry_wavelet_vector) /= size(block%bdry_vector) .or. &
       size(block%bdry_scalar_mean) /= size(block%bdry_scalar) .or. &
       size(block%bdry_vector_mean) /= size(block%bdry_vector) .or. &
       size(block%bdry_tke) /= block%n_tke_level*n_bdry_node .or. &
       size(block%bdry_wavelet_tke) /= size(block%bdry_tke) .or. &
       size(block%bdry_topography) /= n_bdry_node) then

     error stop "check_block_storage: boundary extent mismatch"

  end if

  n_ghost_node = sum(block%ghost_storage%n_node)

  if (size(block%ghost_node) /= n_ghost_node .or. &
       size(block%ghost_scalar) /= &
       block%n_scalar_variable*block%n_field_level* &
       block%scalar_mult*n_ghost_node .or. &
       size(block%ghost_vector) /= &
       block%n_field_level*block%vector_mult*n_ghost_node .or. &
       size(block%ghost_wavelet_scalar) /= size(block%ghost_scalar) .or. &
       size(block%ghost_wavelet_vector) /= size(block%ghost_vector) .or. &
       size(block%ghost_scalar_mean) /= size(block%ghost_scalar) .or. &
       size(block%ghost_vector_mean) /= size(block%ghost_vector) .or. &
       size(block%ghost_tke) /= block%n_tke_level*n_ghost_node .or. &
       size(block%ghost_wavelet_tke) /= size(block%ghost_tke) .or. &
       size(block%ghost_topography) /= n_ghost_node) then

     error stop "check_block_storage: ghost extent mismatch"

  end if

  do i = 1, size(block%ghost_storage)
     if (block%ghost_storage(i)%source_local_patch < 0) then
        error stop "check_block_storage: missing ghost source patch"
     end if
  end do

  if (serialize) then

     call pack_block(block,buffer_source)
     call unpack_block(buffer_source,block_copy)
     call pack_block(block_copy,buffer_copy)

     if (size(buffer_copy) /= size(buffer_source)) then
        error stop &
             "check_block_storage: packed size mismatch"
     end if

     if (any(buffer_copy /= buffer_source)) then
        error stop &
             "check_block_storage: serialization mismatch"
     end if

  end if

end subroutine check_block_storage


subroutine install_local_blocks (n_catalog,local_seen)
  ! Build the final-owner local store from deep copies of the retained
  ! source blocks and validated received blocks. Catalogue ownership is
  ! deliberately checked by the caller, outside this non-MPI module.

  implicit none

  integer, intent(in) :: n_catalog
  integer, allocatable, intent(out) :: local_seen(:)

  integer :: b
  integer :: i
  integer :: ib
  integer :: ilocal
  integer :: n_local

  integer(int8), allocatable :: buffer_local(:)
  integer(int8), allocatable :: buffer_reference(:)

  if (n_catalog < 1) then
     error stop "install_local_blocks: invalid catalogue size"
  end if

  if (.not. allocated(block_source) .or. &
       .not. allocated(block_retained_source_index) .or. &
       .not. allocated(block_source_catalog_index)) then

     error stop "install_local_blocks: retained store missing"

  end if

  if (size(block_source_catalog_index) /= size(block_source)) then
     error stop "install_local_blocks: source map size mismatch"
  end if

  if (.not. allocated(block_received) .or. &
       .not. allocated(block_received_catalog_index)) then

     error stop "install_local_blocks: received store missing"

  end if

  if (size(block_received) /= &
       size(block_received_catalog_index)) then
     error stop "install_local_blocks: received map size mismatch"
  end if

  n_local = size(block_retained_source_index) + &
       size(block_received)

  call clear_local_blocks

  allocate(block_local(n_local))
  allocate(block_local_catalog_index(n_local))
  allocate(block_catalog_local_index(n_catalog))
  allocate(local_seen(n_catalog))

  block_local_catalog_index = -1
  block_catalog_local_index = 0
  local_seen = 0
  ilocal = 0

  do i = 1, size(block_retained_source_index)

     ib = block_retained_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop "install_local_blocks: invalid retained index"
     end if

     b = block_source_catalog_index(ib)

     if (b < 1 .or. b > n_catalog) then
        error stop &
             "install_local_blocks: invalid retained catalogue index"
     end if

     if (local_seen(b) /= 0) then
        error stop "install_local_blocks: duplicate retained block"
     end if

     ilocal = ilocal + 1
     block_local(ilocal) = block_source(ib)
     block_local_catalog_index(ilocal) = b
     block_catalog_local_index(b) = ilocal
     local_seen(b) = 1

     call pack_block(block_local(ilocal),buffer_local)
     call pack_block(block_source(ib),buffer_reference)

     if (size(buffer_local) /= size(buffer_reference)) then
        error stop &
             "install_local_blocks: retained copy size mismatch"
     end if

     if (any(buffer_local /= buffer_reference)) then
        error stop &
             "install_local_blocks: retained deep-copy mismatch"
     end if

  end do

  do i = 1, size(block_received)

     b = block_received_catalog_index(i)

     if (b < 1 .or. b > n_catalog) then
        error stop &
             "install_local_blocks: invalid received catalogue index"
     end if

     if (local_seen(b) /= 0) then
        error stop "install_local_blocks: duplicate received block"
     end if

     ilocal = ilocal + 1
     block_local(ilocal) = block_received(i)
     block_local_catalog_index(ilocal) = b
     block_catalog_local_index(b) = ilocal
     local_seen(b) = 1

     call pack_block(block_local(ilocal),buffer_local)
     call pack_block(block_received(i),buffer_reference)

     if (size(buffer_local) /= size(buffer_reference)) then
        error stop &
             "install_local_blocks: received copy size mismatch"
     end if

     if (any(buffer_local /= buffer_reference)) then
        error stop &
             "install_local_blocks: received deep-copy mismatch"
     end if

  end do

  if (ilocal /= n_local) then
     error stop "install_local_blocks: local fill count mismatch"
  end if

  if (count(local_seen /= 0) /= n_local) then
     error stop "install_local_blocks: local inventory mismatch"
  end if

  if (count(block_catalog_local_index > 0) /= n_local) then
     error stop "install_local_blocks: inverse map count mismatch"
  end if

  do ilocal = 1, n_local

     b = block_local_catalog_index(ilocal)

     if (block_catalog_local_index(b) /= ilocal) then
        error stop "install_local_blocks: inverse map mismatch"
     end if

  end do

  block_store_ready = .true.

end subroutine install_local_blocks


logical function local_block_store_ready () result(ready)
  ! Report whether a complete final-owner local store has been
  ! installed. Allocation checks make a partially modified store
  ! unavailable even if its readiness flag was set previously.

  implicit none

  ready = .false.

  if (.not. block_store_ready) return

  if (.not. allocated(block_local)) return
  if (.not. allocated(block_local_catalog_index)) return
  if (.not. allocated(block_catalog_local_index)) return

  if (size(block_local) /= &
       size(block_local_catalog_index)) return

  if (size(block_catalog_local_index) < 1) return

  ready = .true.

end function local_block_store_ready


integer function n_local_blocks () result(n_block)
  ! Return the number of blocks in the ready final-owner store.

  implicit none

  if (.not. local_block_store_ready()) then
     error stop "n_local_blocks: local store is not ready"
  end if

  n_block = size(block_local)

end function n_local_blocks


integer function local_block_catalog (local_index) result(catalog_index)
  ! Map a valid local block index to its replicated catalogue index.

  implicit none

  integer, intent(in) :: local_index

  if (.not. local_block_store_ready()) then
     error stop "local_block_catalog: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "local_block_catalog: invalid local index"
  end if

  catalog_index = block_local_catalog_index(local_index)

end function local_block_catalog


integer function catalog_local_block (catalog_index) result(local_index)
  ! Map a replicated catalogue index to its local block index. Return
  ! zero when the catalogue block belongs to another rank.

  implicit none

  integer, intent(in) :: catalog_index

  if (.not. local_block_store_ready()) then
     error stop "catalog_local_block: local store is not ready"
  end if

  if (catalog_index < 1 .or. &
       catalog_index > size(block_catalog_local_index)) then
     error stop "catalog_local_block: invalid catalogue index"
  end if

  local_index = block_catalog_local_index(catalog_index)

end function catalog_local_block


subroutine get_local_block_identity ( &
     local_index,id,root_domain,root_patch,level)
  ! Return immutable identity metadata for one local block.

  implicit none

  integer, intent(in)  :: local_index
  integer, intent(out) :: id
  integer, intent(out) :: root_domain
  integer, intent(out) :: root_patch
  integer, intent(out) :: level

  if (.not. local_block_store_ready()) then
     error stop "get_local_block_identity: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "get_local_block_identity: invalid local index"
  end if

  id          = block_local(local_index)%id
  root_domain = block_local(local_index)%root_domain
  root_patch  = block_local(local_index)%root_patch
  level       = block_local(local_index)%level

end subroutine get_local_block_identity


subroutine get_local_block_field_layout ( &
     local_index,scalar_variable,n_scalar_variable, &
     vector_variable,field_level, &
     n_field_level,scalar_mult,vector_mult)
  ! Return the serialized field descriptor for one installed block.

  implicit none

  integer, intent(in)  :: local_index
  integer, intent(out) :: scalar_variable
  integer, intent(out) :: n_scalar_variable
  integer, intent(out) :: vector_variable
  integer, intent(out) :: field_level
  integer, intent(out) :: n_field_level
  integer, intent(out) :: scalar_mult
  integer, intent(out) :: vector_mult

  if (.not. local_block_store_ready()) then
     error stop &
          "get_local_block_field_layout: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop &
          "get_local_block_field_layout: invalid local index"
  end if

  scalar_variable = block_local(local_index)%scalar_variable
  n_scalar_variable = block_local(local_index)%n_scalar_variable
  vector_variable = block_local(local_index)%vector_variable
  field_level     = block_local(local_index)%field_level
  n_field_level   = block_local(local_index)%n_field_level
  scalar_mult     = block_local(local_index)%scalar_mult
  vector_mult     = block_local(local_index)%vector_mult

end subroutine get_local_block_field_layout


subroutine get_local_block_turbulence_layout ( &
     local_index,tke_level,n_tke_level)
  ! Return the serialized optional TKE descriptor for one installed block.

  implicit none

  integer, intent(in)  :: local_index
  integer, intent(out) :: tke_level
  integer, intent(out) :: n_tke_level

  if (.not. local_block_store_ready()) then
     error stop &
          "get_local_block_turbulence_layout: store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop &
          "get_local_block_turbulence_layout: invalid local index"
  end if

  tke_level   = block_local(local_index)%tke_level
  n_tke_level = block_local(local_index)%n_tke_level

end subroutine get_local_block_turbulence_layout


subroutine check_local_block_storage (local_index,check_serialization)
  ! Validate one local block without exposing the private store.

  implicit none

  integer, intent(in) :: local_index
  logical, optional, intent(in) :: check_serialization

  logical :: serialize

  if (.not. local_block_store_ready()) then
     error stop "check_local_block_storage: local store is not ready"
  end if

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop "check_local_block_storage: invalid local index"
  end if

  serialize = .false.
  if (present(check_serialization)) serialize = check_serialization

  call check_block_storage(block_local(local_index),serialize)

end subroutine check_local_block_storage


subroutine apply_local_block_field_consumer (consumer,context)
  ! Apply one caller-supplied read-only production kernel to every local
  ! block. The callback receives the compact interior, boundary and ghost
  ! fields together with their topology without constructing copies.

  implicit none

  procedure(Local_Block_Field_Consumer) :: consumer
  class(*), intent(inout) :: context

  integer :: catalog_index
  integer :: local_index
  integer :: n_block_before

  if (.not. local_block_store_ready()) then
     error stop &
          "apply_local_block_field_consumer: store is not ready"
  end if

  n_block_before = size(block_local)

  do local_index = 1,n_block_before
     catalog_index = block_local_catalog_index(local_index)
     if (catalog_index < 1) then
        error stop &
             "apply_local_block_field_consumer: invalid catalogue index"
     end if

     call consumer(catalog_index,block_local(local_index),context)

     if (.not. local_block_store_ready()) then
        error stop &
             "apply_local_block_field_consumer: store changed"
     end if
     if (size(block_local) /= n_block_before .or. &
          block_local_catalog_index(local_index) /= catalog_index) then
        error stop &
             "apply_local_block_field_consumer: traversal changed"
     end if
  end do

end subroutine apply_local_block_field_consumer


subroutine prepare_local_block_tendency_state
  ! Allocate reusable scalar/vector kernel outputs for the current local
  ! block store. Existing correctly sized allocations are retained.

  implicit none

  integer :: local_index

  if (.not. local_block_store_ready()) then
     error stop &
          "prepare_local_block_tendency_state: store is not ready"
  end if

  if (allocated(block_tendency)) then
     if (size(block_tendency) /= size(block_local)) then
        deallocate(block_tendency)
     end if
  end if
  if (.not. allocated(block_tendency)) then
     allocate(block_tendency(size(block_local)))
     block_tendency_ready = .false.
  end if

  do local_index = 1,size(block_local)
     if (block_tendency(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        if (allocated(block_tendency(local_index)%scalar)) then
           deallocate(block_tendency(local_index)%scalar)
        end if
        if (allocated(block_tendency(local_index)%vector)) then
           deallocate(block_tendency(local_index)%vector)
        end if
        block_tendency(local_index)%catalog_index = &
             block_local_catalog_index(local_index)
        block_tendency(local_index)%ready = .false.
        block_tendency(local_index)%generation = 0_int64
     end if

     if (allocated(block_tendency(local_index)%scalar)) then
        if (size(block_tendency(local_index)%scalar) /= &
             size(block_local(local_index)%scalar)) then
           deallocate(block_tendency(local_index)%scalar)
           block_tendency(local_index)%ready = .false.
        end if
     end if
     if (.not. allocated(block_tendency(local_index)%scalar)) then
        allocate(block_tendency(local_index)%scalar( &
             size(block_local(local_index)%scalar)))
        block_tendency_allocations = block_tendency_allocations + 1_int64
     end if

     if (allocated(block_tendency(local_index)%vector)) then
        if (size(block_tendency(local_index)%vector) /= &
             size(block_local(local_index)%vector)) then
           deallocate(block_tendency(local_index)%vector)
           block_tendency(local_index)%ready = .false.
        end if
     end if
     if (.not. allocated(block_tendency(local_index)%vector)) then
        allocate(block_tendency(local_index)%vector( &
             size(block_local(local_index)%vector)))
        block_tendency_allocations = block_tendency_allocations + 1_int64
     end if
  end do

  block_tendency_ready = all(block_tendency%ready)

end subroutine prepare_local_block_tendency_state


subroutine prepare_local_block_tendency_import_coverage
  ! Allocate reusable per-block import coverage for the current store.

  implicit none

  if (.not. local_block_store_ready()) then
     error stop &
          "prepare_local_block_tendency_import_coverage: store not ready"
  end if

  if (allocated(block_tendency_import_patch_count)) then
     if (size(block_tendency_import_patch_count) /= &
          size(block_local)) then
        deallocate(block_tendency_import_patch_count)
     end if
  end if
  if (.not. allocated(block_tendency_import_patch_count)) then
     allocate(block_tendency_import_patch_count(size(block_local)))
     block_tendency_import_allocations = &
          block_tendency_import_allocations + 1_int64
  end if

  block_tendency_import_patch_count = 0

end subroutine prepare_local_block_tendency_import_coverage


subroutine begin_local_block_tendency_import
  ! Begin a guarded whole-store import into persistent tendency arrays.

  implicit none

  integer :: local_index

  if (block_tendency_import_active) then
     error stop "begin_local_block_tendency_import: import is active"
  end if
  if (block_tendency_trial_active .or. &
       block_tendency_commit_checkpoint_ready) then
     error stop &
          "begin_local_block_tendency_import: transaction pending"
  end if

  call prepare_local_block_tendency_state
  call prepare_local_block_tendency_import_coverage
  block_tendency_ready = .false.
  do local_index = 1,size(block_tendency)
     block_tendency(local_index)%ready = .false.
     block_tendency(local_index)%scalar = 0.0_dp
     block_tendency(local_index)%vector = 0.0_dp
  end do
  block_tendency_import_active = .true.

end subroutine begin_local_block_tendency_import


subroutine prepare_local_block_tendency_workspace
  ! Allocate the reusable production tendency, import-coverage and trial
  ! storage without importing a tendency or activating a transaction.

  implicit none

  if (block_tendency_import_active .or. &
       block_tendency_trial_active .or. &
       block_tendency_commit_checkpoint_ready) then
     error stop &
          "prepare_local_block_tendency_workspace: transaction pending"
  end if

  call prepare_local_block_tendency_state
  call prepare_local_block_tendency_import_coverage
  call prepare_local_block_tendency_trial
  call invalidate_local_block_tendency_products

end subroutine prepare_local_block_tendency_workspace


subroutine set_local_block_tendency_patch_values ( &
     catalog_index,local_patch,scalar_value,vector_value)
  ! Install one compact patch during an active preorder import.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(in) :: scalar_value(:)
  real(dp), intent(in) :: vector_value(:)

  integer :: field_base
  integer :: input_base
  integer :: expected_scalar_nvalue
  integer :: expected_vector_nvalue
  integer :: level_slot
  integer :: local_index
  integer :: n_node
  integer :: n_patch_value
  integer :: patch_start
  integer :: scalar_slot

  if (.not. block_tendency_import_active) then
     error stop &
          "set_local_block_tendency_patch_values: no active import"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_tendency_patch_values: block is not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "set_local_block_tendency_patch_values: invalid patch"
  end if
  if (local_patch /= &
       block_tendency_import_patch_count(local_index)) then
     error stop &
          "set_local_block_tendency_patch_values: patch order mismatch"
  end if
  expected_scalar_nvalue = &
       local_block_scalar_family_patch_nvalue(catalog_index)
  expected_vector_nvalue = &
       local_block_vector_family_patch_nvalue(catalog_index)
  if (size(scalar_value) /= expected_scalar_nvalue .or. &
       size(vector_value) /= expected_vector_nvalue) then
     error stop &
          "set_local_block_tendency_patch_values: payload extent"
  end if

  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start
  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "set_local_block_tendency_patch_values: patch storage"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult*PATCH_SIZE**2
  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)* &
             block_local(local_index)%scalar_mult*n_node
        input_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)*n_patch_value
        block_tendency(local_index)%scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + &
             n_patch_value) = &
             scalar_value(input_base+1:input_base+n_patch_value)
     end do
  end do

  n_patch_value = &
       block_local(local_index)%vector_mult*PATCH_SIZE**2
  do level_slot = 1,block_local(local_index)%n_field_level
     field_base = (level_slot-1)* &
          block_local(local_index)%vector_mult*n_node
     input_base = (level_slot-1)*n_patch_value
     block_tendency(local_index)%vector( &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + &
          n_patch_value) = &
          vector_value(input_base+1:input_base+n_patch_value)
  end do

  block_tendency_import_patch_count(local_index) = &
       block_tendency_import_patch_count(local_index) + 1

end subroutine set_local_block_tendency_patch_values


subroutine finish_local_block_tendency_import
  ! Commit readiness only after every local patch was imported exactly once.

  implicit none

  integer :: local_index

  if (.not. block_tendency_import_active) then
     error stop "finish_local_block_tendency_import: no active import"
  end if

  do local_index = 1,size(block_local)
     if (block_tendency_import_patch_count(local_index) /= &
          size(block_local(local_index)%patch)) then
        error stop &
             "finish_local_block_tendency_import: incomplete block"
     end if
     block_tendency(local_index)%ready = .true.
     block_tendency(local_index)%generation = &
          block_tendency(local_index)%generation + 1_int64
  end do

  block_tendency_import_active = .false.
  block_tendency_ready = all(block_tendency%ready)
  if (.not. local_block_tendency_state_ready()) then
     error stop &
          "finish_local_block_tendency_import: state not ready"
  end if

end subroutine finish_local_block_tendency_import


subroutine get_local_block_tendency_patch_values ( &
     catalog_index,local_patch,scalar_value,vector_value)
  ! Read one imported/kernel-produced tendency patch in compact field order.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(out) :: scalar_value(:)
  real(dp), intent(out) :: vector_value(:)

  integer :: field_base
  integer :: expected_scalar_nvalue
  integer :: expected_vector_nvalue
  integer :: level_slot
  integer :: local_index
  integer :: n_node
  integer :: n_patch_value
  integer :: output_base
  integer :: patch_start
  integer :: scalar_slot

  if (.not. local_block_tendency_state_ready()) then
     error stop &
          "get_local_block_tendency_patch_values: state not ready"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_tendency_patch_values: block is not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_tendency_patch_values: invalid patch"
  end if
  expected_scalar_nvalue = &
       local_block_scalar_family_patch_nvalue(catalog_index)
  expected_vector_nvalue = &
       local_block_vector_family_patch_nvalue(catalog_index)
  if (size(scalar_value) /= expected_scalar_nvalue .or. &
       size(vector_value) /= expected_vector_nvalue) then
     error stop &
          "get_local_block_tendency_patch_values: output extent"
  end if

  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start
  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "get_local_block_tendency_patch_values: patch storage"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult*PATCH_SIZE**2
  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)* &
             block_local(local_index)%scalar_mult*n_node
        output_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)*n_patch_value
        scalar_value(output_base+1:output_base+n_patch_value) = &
             block_tendency(local_index)%scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + &
             n_patch_value)
     end do
  end do

  n_patch_value = &
       block_local(local_index)%vector_mult*PATCH_SIZE**2
  do level_slot = 1,block_local(local_index)%n_field_level
     field_base = (level_slot-1)* &
          block_local(local_index)%vector_mult*n_node
     output_base = (level_slot-1)*n_patch_value
     vector_value(output_base+1:output_base+n_patch_value) = &
          block_tendency(local_index)%vector( &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + &
          n_patch_value)
  end do

end subroutine get_local_block_tendency_patch_values


subroutine assert_local_block_tendency_patch_values ( &
     catalog_index,local_patch,scalar_value,vector_value)
  ! Compare one ready tendency patch directly with an external compact
  ! payload. This avoids allocating diagnostic copies in production paths.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(in) :: scalar_value(:)
  real(dp), intent(in) :: vector_value(:)

  integer :: field_base
  integer :: input_base
  integer :: expected_scalar_nvalue
  integer :: expected_vector_nvalue
  integer :: level_slot
  integer :: local_index
  integer :: n_node
  integer :: n_patch_value
  integer :: patch_start
  integer :: scalar_slot

  if (.not. block_tendency_ready .or. &
       .not. allocated(block_tendency)) then
     error stop &
          "assert_local_block_tendency_patch_values: state not ready"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "assert_local_block_tendency_patch_values: block is not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "assert_local_block_tendency_patch_values: invalid patch"
  end if
  expected_scalar_nvalue = &
       local_block_scalar_family_patch_nvalue(catalog_index)
  expected_vector_nvalue = &
       local_block_vector_family_patch_nvalue(catalog_index)
  if (size(scalar_value) /= expected_scalar_nvalue .or. &
       size(vector_value) /= expected_vector_nvalue) then
     error stop &
          "assert_local_block_tendency_patch_values: payload extent"
  end if

  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start
  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "assert_local_block_tendency_patch_values: patch storage"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult*PATCH_SIZE**2
  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)* &
             block_local(local_index)%scalar_mult*n_node
        input_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1)*n_patch_value
        if (any(abs(block_tendency(local_index)%scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + &
             n_patch_value) - &
             scalar_value(input_base+1:input_base+n_patch_value)) > &
             0.0_dp)) then
           error stop &
                "assert_local_block_tendency_patch_values: scalar mismatch"
        end if
     end do
  end do

  n_patch_value = &
       block_local(local_index)%vector_mult*PATCH_SIZE**2
  do level_slot = 1,block_local(local_index)%n_field_level
     field_base = (level_slot-1)* &
          block_local(local_index)%vector_mult*n_node
     input_base = (level_slot-1)*n_patch_value
     if (any(abs(block_tendency(local_index)%vector( &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + &
          n_patch_value) - &
          vector_value(input_base+1:input_base+n_patch_value)) > &
          0.0_dp)) then
        error stop &
             "assert_local_block_tendency_patch_values: vector mismatch"
     end if
  end do

end subroutine assert_local_block_tendency_patch_values


logical function local_block_tendency_import_is_active () result(active)

  implicit none

  active = block_tendency_import_active

end function local_block_tendency_import_is_active


integer(int64) function local_block_tendency_import_allocation_count () &
     result(n_allocation)

  implicit none

  n_allocation = block_tendency_import_allocations

end function local_block_tendency_import_allocation_count


subroutine apply_local_block_tendency_kernel (kernel,context)
  ! Execute one writable production kernel across every local block.
  ! Output arrays are persistent, zeroed before use and kept separate from
  ! authoritative prognostic and derived-cache storage.

  implicit none

  procedure(Local_Block_Tendency_Kernel) :: kernel
  class(*), intent(inout) :: context

  integer :: catalog_index
  integer :: local_index

  if (block_tendency_import_active) then
     error stop &
          "apply_local_block_tendency_kernel: import is active"
  end if
  if (block_tendency_trial_active) then
     error stop &
          "apply_local_block_tendency_kernel: trial is active"
  end if

  call prepare_local_block_tendency_state
  block_tendency_ready = .false.

  do local_index = 1,size(block_local)
     catalog_index = block_local_catalog_index(local_index)
     block_tendency(local_index)%ready = .false.
     block_tendency(local_index)%scalar = 0.0_dp
     block_tendency(local_index)%vector = 0.0_dp

     call kernel( &
          catalog_index,block_local(local_index), &
          block_tendency(local_index)%scalar, &
          block_tendency(local_index)%vector,context)

     if (size(block_tendency(local_index)%scalar) /= &
          size(block_local(local_index)%scalar) .or. &
          size(block_tendency(local_index)%vector) /= &
          size(block_local(local_index)%vector)) then
        error stop &
             "apply_local_block_tendency_kernel: output extent changed"
     end if

     block_tendency(local_index)%ready = .true.
     block_tendency(local_index)%generation = &
          block_tendency(local_index)%generation + 1_int64
  end do

  block_tendency_ready = all(block_tendency%ready)
  block_tendency_executions = block_tendency_executions + 1_int64

  if (.not. local_block_tendency_state_ready()) then
     error stop &
          "apply_local_block_tendency_kernel: output state not ready"
  end if

end subroutine apply_local_block_tendency_kernel


subroutine apply_local_block_tendency_consumer (consumer,context)
  ! Traverse the persistent tendency outputs without copying or exposing
  ! writable storage. The consumer cannot alter either output family.

  implicit none

  procedure(Local_Block_Tendency_Consumer) :: consumer
  class(*), intent(inout) :: context

  integer :: catalog_index
  integer :: local_index
  integer(int64) :: generation_before

  if (.not. local_block_tendency_state_ready()) then
     error stop &
          "apply_local_block_tendency_consumer: output state not ready"
  end if
  if (block_tendency_trial_active) then
     error stop &
          "apply_local_block_tendency_consumer: trial is active"
  end if

  do local_index = 1,size(block_tendency)
     catalog_index = block_tendency(local_index)%catalog_index
     generation_before = block_tendency(local_index)%generation

     call consumer( &
          catalog_index,block_tendency(local_index)%scalar, &
          block_tendency(local_index)%vector,context)

     if (.not. local_block_tendency_state_ready()) then
        error stop &
             "apply_local_block_tendency_consumer: state changed"
     end if
     if (block_tendency(local_index)%catalog_index /= catalog_index .or. &
          block_tendency(local_index)%generation /= generation_before) then
        error stop &
             "apply_local_block_tendency_consumer: traversal changed"
     end if
  end do

end subroutine apply_local_block_tendency_consumer


subroutine prepare_local_block_tendency_accumulator
  ! Allocate reusable scalar/vector registers for weighted tendency stages.

  implicit none

  integer :: local_index

  if (.not. local_block_store_ready()) then
     error stop &
          "prepare_local_block_tendency_accumulator: store not ready"
  end if
  if (block_tendency_import_active) then
     error stop &
          "prepare_local_block_tendency_accumulator: import active"
  end if
  if (block_tendency_trial_active .or. &
       block_tendency_commit_checkpoint_ready) then
     error stop &
          "prepare_local_block_tendency_accumulator: transaction pending"
  end if

  if (allocated(block_tendency_accumulator)) then
     if (size(block_tendency_accumulator) /= size(block_local)) then
        deallocate(block_tendency_accumulator)
        block_tendency_accumulator_ready = .false.
     end if
  end if
  if (.not. allocated(block_tendency_accumulator)) then
     allocate(block_tendency_accumulator(size(block_local)))
     block_tendency_accumulator_ready = .false.
  end if

  do local_index = 1,size(block_local)
     if (block_tendency_accumulator(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        if (allocated( &
             block_tendency_accumulator(local_index)%scalar)) then
           deallocate(block_tendency_accumulator(local_index)%scalar)
        end if
        if (allocated( &
             block_tendency_accumulator(local_index)%vector)) then
           deallocate(block_tendency_accumulator(local_index)%vector)
        end if
        block_tendency_accumulator(local_index)%catalog_index = &
             block_local_catalog_index(local_index)
        block_tendency_accumulator_ready = .false.
     end if

     if (allocated(block_tendency_accumulator(local_index)%scalar)) then
        if (size(block_tendency_accumulator(local_index)%scalar) /= &
             size(block_local(local_index)%scalar)) then
           deallocate(block_tendency_accumulator(local_index)%scalar)
           block_tendency_accumulator_ready = .false.
        end if
     end if
     if (.not. allocated( &
          block_tendency_accumulator(local_index)%scalar)) then
        allocate(block_tendency_accumulator(local_index)%scalar( &
             size(block_local(local_index)%scalar)))
        block_tendency_accumulator_allocations = &
             block_tendency_accumulator_allocations + 1_int64
     end if

     if (allocated(block_tendency_accumulator(local_index)%vector)) then
        if (size(block_tendency_accumulator(local_index)%vector) /= &
             size(block_local(local_index)%vector)) then
           deallocate(block_tendency_accumulator(local_index)%vector)
           block_tendency_accumulator_ready = .false.
        end if
     end if
     if (.not. allocated( &
          block_tendency_accumulator(local_index)%vector)) then
        allocate(block_tendency_accumulator(local_index)%vector( &
             size(block_local(local_index)%vector)))
        block_tendency_accumulator_allocations = &
             block_tendency_accumulator_allocations + 1_int64
     end if
  end do

end subroutine prepare_local_block_tendency_accumulator


subroutine reset_local_block_tendency_accumulator
  ! Begin a new weighted multi-stage tendency combination.

  implicit none

  integer :: local_index

  call prepare_local_block_tendency_accumulator

  do local_index = 1,size(block_tendency_accumulator)
     block_tendency_accumulator(local_index)%scalar = 0.0_dp
     block_tendency_accumulator(local_index)%vector = 0.0_dp
  end do

  block_tendency_accumulator_stages = 0_int64
  block_tendency_accumulator_ready = .true.

end subroutine reset_local_block_tendency_accumulator


subroutine accumulate_local_block_tendency (weight)
  ! Add one weighted persistent tendency output to the multi-stage register.

  implicit none

  real(dp), intent(in) :: weight

  integer :: local_index

  if (.not. local_block_tendency_state_ready()) then
     error stop "accumulate_local_block_tendency: tendency not ready"
  end if
  if (.not. local_block_tendency_accumulator_state_ready()) then
     error stop "accumulate_local_block_tendency: accumulator not ready"
  end if
  if (block_tendency_trial_active) then
     error stop "accumulate_local_block_tendency: trial is active"
  end if

  do local_index = 1,size(block_local)
     block_tendency_accumulator(local_index)%scalar = &
          block_tendency_accumulator(local_index)%scalar + &
          weight*block_tendency(local_index)%scalar
     block_tendency_accumulator(local_index)%vector = &
          block_tendency_accumulator(local_index)%vector + &
          weight*block_tendency(local_index)%vector
  end do

  block_tendency_accumulator_stages = &
       block_tendency_accumulator_stages + 1_int64

end subroutine accumulate_local_block_tendency


logical function local_block_tendency_accumulator_state_ready () &
     result(ready)
  ! Report whether the current catalogue has a complete stage register.

  implicit none

  integer :: local_index

  ready = .false.
  if (block_tendency_import_active) return
  if (.not. local_block_store_ready()) return
  if (.not. block_tendency_accumulator_ready) return
  if (.not. allocated(block_tendency_accumulator)) return
  if (size(block_tendency_accumulator) /= size(block_local)) return

  do local_index = 1,size(block_local)
     if (block_tendency_accumulator(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) return
     if (.not. allocated( &
          block_tendency_accumulator(local_index)%scalar)) return
     if (.not. allocated( &
          block_tendency_accumulator(local_index)%vector)) return
     if (size(block_tendency_accumulator(local_index)%scalar) /= &
          size(block_local(local_index)%scalar)) return
     if (size(block_tendency_accumulator(local_index)%vector) /= &
          size(block_local(local_index)%vector)) return
  end do

  ready = .true.

end function local_block_tendency_accumulator_state_ready


integer(int64) function local_block_tendency_accumulator_allocation_count () &
     result(n_allocation)
  ! Number of scalar/vector stage-register allocations since installation.

  implicit none

  n_allocation = block_tendency_accumulator_allocations

end function local_block_tendency_accumulator_allocation_count


subroutine local_block_tendency_accumulator_statistics ( &
     scalar_count,vector_count,scalar_changed_block_count, &
     vector_changed_block_count,stage_count,scalar_moment,vector_moment)
  ! Accumulate diagnostics directly from the persistent stage register.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count
  integer(int64), intent(out) :: scalar_changed_block_count
  integer(int64), intent(out) :: vector_changed_block_count
  integer(int64), intent(out) :: stage_count
  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: local_index

  if (.not. local_block_tendency_accumulator_state_ready()) then
     error stop &
          "local_block_tendency_accumulator_statistics: not ready"
  end if

  scalar_count = 0_int64
  vector_count = 0_int64
  scalar_changed_block_count = 0_int64
  vector_changed_block_count = 0_int64
  stage_count = block_tendency_accumulator_stages
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do local_index = 1,size(block_tendency_accumulator)
     scalar_count = scalar_count + int(size( &
          block_tendency_accumulator(local_index)%scalar),int64)
     vector_count = vector_count + int(size( &
          block_tendency_accumulator(local_index)%vector),int64)
     if (maxval(abs( &
          block_tendency_accumulator(local_index)%scalar)) > 0.0_dp) then
        scalar_changed_block_count = &
             scalar_changed_block_count + 1_int64
     end if
     if (maxval(abs( &
          block_tendency_accumulator(local_index)%vector)) > 0.0_dp) then
        vector_changed_block_count = &
             vector_changed_block_count + 1_int64
     end if

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_tendency_accumulator(local_index)%scalar)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_tendency_accumulator(local_index)%scalar))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_tendency_accumulator(local_index)%scalar**2)
     vector_moment(1) = vector_moment(1) + &
          sum(block_tendency_accumulator(local_index)%vector)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_tendency_accumulator(local_index)%vector))
     vector_moment(3) = vector_moment(3) + &
          sum(block_tendency_accumulator(local_index)%vector**2)
  end do

end subroutine local_block_tendency_accumulator_statistics


subroutine prepare_local_block_tendency_trial
  ! Allocate reusable exact snapshots for reversible shadow updates.

  implicit none

  integer :: local_index

  if (block_tendency_trial_active) then
     error stop "prepare_local_block_tendency_trial: trial is active"
  end if
  if (block_tendency_commit_checkpoint_ready) then
     error stop &
          "prepare_local_block_tendency_trial: commit checkpoint pending"
  end if

  if (allocated(block_tendency_trial)) then
     if (size(block_tendency_trial) /= size(block_local)) then
        deallocate(block_tendency_trial)
     end if
  end if
  if (.not. allocated(block_tendency_trial)) then
     allocate(block_tendency_trial(size(block_local)))
  end if

  do local_index = 1,size(block_local)
     if (block_tendency_trial(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        if (allocated(block_tendency_trial(local_index)%scalar)) then
           deallocate(block_tendency_trial(local_index)%scalar)
        end if
        if (allocated(block_tendency_trial(local_index)%vector)) then
           deallocate(block_tendency_trial(local_index)%vector)
        end if
        block_tendency_trial(local_index)%catalog_index = &
             block_local_catalog_index(local_index)
     end if

     if (allocated(block_tendency_trial(local_index)%scalar)) then
        if (size(block_tendency_trial(local_index)%scalar) /= &
             size(block_local(local_index)%scalar)) then
           deallocate(block_tendency_trial(local_index)%scalar)
        end if
     end if
     if (.not. allocated(block_tendency_trial(local_index)%scalar)) then
        allocate(block_tendency_trial(local_index)%scalar( &
             size(block_local(local_index)%scalar)))
     end if

     if (allocated(block_tendency_trial(local_index)%vector)) then
        if (size(block_tendency_trial(local_index)%vector) /= &
             size(block_local(local_index)%vector)) then
           deallocate(block_tendency_trial(local_index)%vector)
        end if
     end if
     if (.not. allocated(block_tendency_trial(local_index)%vector)) then
        allocate(block_tendency_trial(local_index)%vector( &
             size(block_local(local_index)%vector)))
     end if
     block_tendency_trial(local_index)%active = .false.
  end do

end subroutine prepare_local_block_tendency_trial


subroutine begin_local_block_tendency_trial (scale)
  ! Snapshot authoritative interior fields, apply one scaled tendency and
  ! invalidate only hydrostatic caches affected by scalar changes.

  implicit none

  real(dp), intent(in) :: scale

  integer :: local_index

  if (.not. local_block_tendency_state_ready()) then
     error stop "begin_local_block_tendency_trial: tendency not ready"
  end if
  if (block_tendency_trial_active) then
     error stop "begin_local_block_tendency_trial: trial already active"
  end if

  call prepare_local_block_tendency_trial

  do local_index = 1,size(block_local)
     block_tendency_trial(local_index)%scalar = &
          block_local(local_index)%scalar
     block_tendency_trial(local_index)%vector = &
          block_local(local_index)%vector

     block_local(local_index)%scalar = &
          block_tendency_trial(local_index)%scalar + &
          scale*block_tendency(local_index)%scalar
     block_local(local_index)%vector = &
          block_tendency_trial(local_index)%vector + &
          scale*block_tendency(local_index)%vector

     block_tendency_trial(local_index)%active = .true.
     if (abs(scale) > 0.0_dp .and. &
          maxval(abs(block_tendency(local_index)%scalar)) > 0.0_dp) then
        if (maxval(abs(block_local(local_index)%scalar - &
             block_tendency_trial(local_index)%scalar)) <= 0.0_dp) then
           error stop &
                "begin_local_block_tendency_trial: scalar update vanished"
        end if
     end if
     if (maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp) then
        call invalidate_local_block_hydrostatic_block(local_index)
     end if
     if (abs(scale) > 0.0_dp .and. &
          maxval(abs(block_tendency(local_index)%vector)) > 0.0_dp) then
        if (maxval(abs(block_local(local_index)%vector - &
             block_tendency_trial(local_index)%vector)) <= 0.0_dp) then
           error stop &
                "begin_local_block_tendency_trial: vector update vanished"
        end if
     end if
  end do

  block_tendency_trial_active = &
       all(block_tendency_trial%active)

  if (.not. block_tendency_trial_active) then
     error stop "begin_local_block_tendency_trial: incomplete trial"
  end if

end subroutine begin_local_block_tendency_trial


subroutine begin_local_block_accumulated_tendency_trial (scale)
  ! Snapshot authoritative fields and apply the current weighted multi-stage
  ! tendency register as one reversible trial update.

  implicit none

  real(dp), intent(in) :: scale

  integer :: local_index

  if (.not. local_block_tendency_accumulator_state_ready()) then
     error stop &
          "begin_local_block_accumulated_tendency_trial: not ready"
  end if
  if (block_tendency_accumulator_stages < 1_int64) then
     error stop &
          "begin_local_block_accumulated_tendency_trial: no stages"
  end if

  call prepare_local_block_tendency_trial

  do local_index = 1,size(block_local)
     block_tendency_trial(local_index)%scalar = &
          block_local(local_index)%scalar
     block_tendency_trial(local_index)%vector = &
          block_local(local_index)%vector

     block_local(local_index)%scalar = &
          block_tendency_trial(local_index)%scalar + scale* &
          block_tendency_accumulator(local_index)%scalar
     block_local(local_index)%vector = &
          block_tendency_trial(local_index)%vector + scale* &
          block_tendency_accumulator(local_index)%vector

     block_tendency_trial(local_index)%active = .true.
     if (abs(scale) > 0.0_dp .and. maxval(abs( &
          block_tendency_accumulator(local_index)%scalar)) > 0.0_dp) then
        if (maxval(abs(block_local(local_index)%scalar - &
             block_tendency_trial(local_index)%scalar)) <= 0.0_dp) then
           error stop &
                "begin_local_block_accumulated_tendency_trial: scalar vanished"
        end if
     end if
     if (maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp) then
        call invalidate_local_block_hydrostatic_block(local_index)
     end if
     if (abs(scale) > 0.0_dp .and. maxval(abs( &
          block_tendency_accumulator(local_index)%vector)) > 0.0_dp) then
        if (maxval(abs(block_local(local_index)%vector - &
             block_tendency_trial(local_index)%vector)) <= 0.0_dp) then
           error stop &
                "begin_local_block_accumulated_tendency_trial: vector vanished"
        end if
     end if
  end do

  block_tendency_trial_active = all(block_tendency_trial%active)
  if (.not. block_tendency_trial_active) then
     error stop &
          "begin_local_block_accumulated_tendency_trial: incomplete"
  end if

end subroutine begin_local_block_accumulated_tendency_trial


subroutine commit_local_block_tendency_trial
  ! Keep the trial fields as the new authoritative block state. Tendencies
  ! were computed from the preceding state and are therefore marked stale;
  ! their allocations remain available for the next kernel execution.

  implicit none

  integer :: local_index

  if (.not. block_tendency_trial_active) then
     error stop "commit_local_block_tendency_trial: no active trial"
  end if
  if (.not. allocated(block_tendency_trial) .or. &
       .not. allocated(block_tendency)) then
     error stop "commit_local_block_tendency_trial: state missing"
  end if

  do local_index = 1,size(block_local)
     if (.not. block_tendency_trial(local_index)%active .or. &
          block_tendency_trial(local_index)%catalog_index /= &
          block_local_catalog_index(local_index) .or. &
          block_tendency(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        error stop "commit_local_block_tendency_trial: invalid state"
     end if

     block_tendency_trial(local_index)%active = .false.
     block_tendency(local_index)%ready = .false.
  end do

  block_tendency_trial_active = .false.
  block_tendency_commit_checkpoint_ready = .true.
  block_tendency_ready = .false.

end subroutine commit_local_block_tendency_trial


subroutine finalize_local_block_tendency_commit
  ! Accept the committed fields permanently and release the one-level
  ! recovery checkpoint for reuse by the next trial update.

  implicit none

  if (block_tendency_trial_active) then
     error stop "finalize_local_block_tendency_commit: trial is active"
  end if
  if (.not. block_tendency_commit_checkpoint_ready) then
     error stop &
          "finalize_local_block_tendency_commit: checkpoint not ready"
  end if

  block_tendency_commit_checkpoint_ready = .false.

end subroutine finalize_local_block_tendency_commit


subroutine restore_local_block_tendency_commit
  ! Restore the exact pre-commit fields retained in the one-level recovery
  ! checkpoint. Current tendencies and affected hydrostatic caches are stale.

  implicit none

  integer :: local_index
  logical :: scalar_changed

  if (block_tendency_trial_active) then
     error stop "restore_local_block_tendency_commit: trial is active"
  end if
  if (.not. block_tendency_commit_checkpoint_ready .or. &
       .not. allocated(block_tendency_trial)) then
     error stop "restore_local_block_tendency_commit: checkpoint not ready"
  end if

  do local_index = 1,size(block_local)
     if (block_tendency_trial(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        error stop "restore_local_block_tendency_commit: invalid checkpoint"
     end if

     scalar_changed = maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp

     block_local(local_index)%scalar = &
          block_tendency_trial(local_index)%scalar
     block_local(local_index)%vector = &
          block_tendency_trial(local_index)%vector

     if (maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp .or. &
          maxval(abs(block_local(local_index)%vector - &
          block_tendency_trial(local_index)%vector)) > 0.0_dp) then
        error stop "restore_local_block_tendency_commit: restore failed"
     end if

     if (scalar_changed) then
        call invalidate_local_block_hydrostatic_block(local_index)
     end if
     if (allocated(block_tendency)) then
        block_tendency(local_index)%ready = .false.
     end if
  end do

  block_tendency_ready = .false.
  block_tendency_commit_checkpoint_ready = .false.

end subroutine restore_local_block_tendency_commit


logical function local_block_tendency_commit_checkpoint_is_ready () &
     result(ready)
  ! Report whether a committed update can still be restored exactly.

  implicit none

  ready = block_tendency_commit_checkpoint_ready

end function local_block_tendency_commit_checkpoint_is_ready


subroutine local_block_tendency_commit_checkpoint_statistics ( &
     scalar_changed_block_count,vector_changed_block_count, &
     scalar_max_update,vector_max_update)
  ! Measure the accepted update directly against its retained checkpoint.

  implicit none

  integer(int64), intent(out) :: scalar_changed_block_count
  integer(int64), intent(out) :: vector_changed_block_count
  real(dp), intent(out) :: scalar_max_update
  real(dp), intent(out) :: vector_max_update

  integer :: local_index
  real(dp) :: update

  if (.not. block_tendency_commit_checkpoint_ready .or. &
       .not. allocated(block_tendency_trial)) then
     error stop &
          "local_block_tendency_commit_checkpoint_statistics: not ready"
  end if

  scalar_changed_block_count = 0_int64
  vector_changed_block_count = 0_int64
  scalar_max_update = 0.0_dp
  vector_max_update = 0.0_dp

  do local_index = 1,size(block_local)
     if (block_tendency_trial(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        error stop &
             "local_block_tendency_commit_checkpoint_statistics: invalid"
     end if

     update = maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar))
     scalar_max_update = max(scalar_max_update,update)
     if (update > 0.0_dp) then
        scalar_changed_block_count = &
             scalar_changed_block_count + 1_int64
     end if

     update = maxval(abs(block_local(local_index)%vector - &
          block_tendency_trial(local_index)%vector))
     vector_max_update = max(vector_max_update,update)
     if (update > 0.0_dp) then
        vector_changed_block_count = &
             vector_changed_block_count + 1_int64
     end if
  end do

end subroutine local_block_tendency_commit_checkpoint_statistics


subroutine rollback_local_block_tendency_trial
  ! Restore the exact saved fields. Derived hydrostatic values remain stale
  ! until the next explicit ensure, because they may describe trial fields.

  implicit none

  integer :: local_index
  logical :: scalar_changed

  if (.not. block_tendency_trial_active) then
     error stop "rollback_local_block_tendency_trial: no active trial"
  end if
  if (.not. allocated(block_tendency_trial)) then
     error stop "rollback_local_block_tendency_trial: snapshot missing"
  end if

  do local_index = 1,size(block_local)
     if (.not. block_tendency_trial(local_index)%active .or. &
          block_tendency_trial(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) then
        error stop "rollback_local_block_tendency_trial: invalid snapshot"
     end if

     scalar_changed = maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp

     block_local(local_index)%scalar = &
          block_tendency_trial(local_index)%scalar
     block_local(local_index)%vector = &
          block_tendency_trial(local_index)%vector

     if (maxval(abs(block_local(local_index)%scalar - &
          block_tendency_trial(local_index)%scalar)) > 0.0_dp .or. &
          maxval(abs(block_local(local_index)%vector - &
          block_tendency_trial(local_index)%vector)) > 0.0_dp) then
        error stop "rollback_local_block_tendency_trial: restore failed"
     end if

     if (scalar_changed) then
        call invalidate_local_block_hydrostatic_block(local_index)
     end if
     block_tendency_trial(local_index)%active = .false.
  end do

  block_tendency_trial_active = .false.

end subroutine rollback_local_block_tendency_trial


logical function local_block_tendency_trial_is_active () result(active)
  ! Report whether authoritative fields currently contain a trial update.

  implicit none

  active = block_tendency_trial_active

end function local_block_tendency_trial_is_active


subroutine invalidate_local_block_tendency_products
  ! Mark tendency and accumulated-stage products stale after authoritative
  ! block sol is replaced from an external Domain representation.

  implicit none

  integer :: local_index

  if (block_tendency_import_active .or. &
       block_tendency_trial_active .or. &
       block_tendency_commit_checkpoint_ready) then
     error stop &
          "invalidate_local_block_tendency_products: transaction pending"
  end if

  block_tendency_ready = .false.
  if (allocated(block_tendency)) then
     do local_index = 1,size(block_tendency)
        block_tendency(local_index)%ready = .false.
     end do
  end if

  block_tendency_accumulator_ready = .false.
  block_tendency_accumulator_stages = 0_int64

end subroutine invalidate_local_block_tendency_products


subroutine discard_local_block_tendency_output
  ! Mark one completed diagnostic tendency evaluation stale without touching
  ! an enclosing committed-state checkpoint. This supports provisional RK
  ! stage evaluation that must never become an accepted update.

  implicit none

  integer :: local_index

  logical :: tendency_ready

  if (block_tendency_import_active .or. block_tendency_trial_active) then
     error stop &
          "discard_local_block_tendency_output: active operation"
  end if
  if (block_tendency_accumulator_ready) then
     error stop &
          "discard_local_block_tendency_output: accumulator is ready"
  end if

  tendency_ready = local_block_tendency_state_ready()
  if (.not. tendency_ready) then
     error stop &
          "discard_local_block_tendency_output: output is not ready"
  end if

  block_tendency_ready = .false.
  do local_index = 1,size(block_tendency)
     block_tendency(local_index)%ready = .false.
  end do

end subroutine discard_local_block_tendency_output


logical function local_block_tendency_state_ready () result(ready)
  ! Report whether every current local block has a complete tendency output.

  implicit none

  integer :: local_index

  ready = .false.
  if (block_tendency_import_active) return
  if (.not. local_block_store_ready()) return
  if (.not. block_tendency_ready) return
  if (.not. allocated(block_tendency)) return
  if (size(block_tendency) /= size(block_local)) return

  do local_index = 1,size(block_local)
     if (block_tendency(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) return
     if (.not. block_tendency(local_index)%ready) return
     if (.not. allocated(block_tendency(local_index)%scalar)) return
     if (.not. allocated(block_tendency(local_index)%vector)) return
     if (size(block_tendency(local_index)%scalar) /= &
          size(block_local(local_index)%scalar)) return
     if (size(block_tendency(local_index)%vector) /= &
          size(block_local(local_index)%vector)) return
  end do

  ready = .true.

end function local_block_tendency_state_ready


integer(int64) function local_block_tendency_execution_count () &
     result(n_execution)
  ! Number of whole-store tendency-kernel executions since installation.

  implicit none

  n_execution = block_tendency_executions

end function local_block_tendency_execution_count


integer(int64) function local_block_tendency_allocation_count () &
     result(n_allocation)
  ! Number of scalar/vector workspace allocations since installation.

  implicit none

  n_allocation = block_tendency_allocations

end function local_block_tendency_allocation_count


subroutine local_block_tendency_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Accumulate order-independent diagnostics directly from persistent
  ! writable kernel outputs.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count
  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: local_index

  if (.not. local_block_tendency_state_ready()) then
     error stop &
          "local_block_tendency_statistics: output state not ready"
  end if

  scalar_count = 0_int64
  vector_count = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do local_index = 1,size(block_tendency)
     scalar_count = scalar_count + &
          int(size(block_tendency(local_index)%scalar),int64)
     vector_count = vector_count + &
          int(size(block_tendency(local_index)%vector),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_tendency(local_index)%scalar)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_tendency(local_index)%scalar))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_tendency(local_index)%scalar**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_tendency(local_index)%vector)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_tendency(local_index)%vector))
     vector_moment(3) = vector_moment(3) + &
          sum(block_tendency(local_index)%vector**2)
  end do

end subroutine local_block_tendency_statistics


subroutine local_block_field_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Compute order-independent field inventory statistics over the
  ! ready final-owner block store. This is the first read-only field
  ! consumer of the private store; no block data are exposed or
  ! modified.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count

  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_field_statistics: local store is not ready"
  end if

  scalar_count  = 0_int64
  vector_count  = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%scalar) .or. &
          .not. allocated(block_local(i)%vector)) then
        error stop &
             "local_block_field_statistics: field storage missing"
     end if

     scalar_count = scalar_count + &
          int(size(block_local(i)%scalar),int64)
     vector_count = vector_count + &
          int(size(block_local(i)%vector),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_local(i)%scalar)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_local(i)%scalar))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_local(i)%scalar**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_local(i)%vector)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_local(i)%vector))
     vector_moment(3) = vector_moment(3) + &
          sum(block_local(i)%vector**2)

  end do

end subroutine local_block_field_statistics


subroutine local_block_wavelet_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Compute order-independent wav_coeff inventory statistics over the
  ! ready final-owner block store.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count

  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_wavelet_statistics: store is not ready"
  end if

  scalar_count  = 0_int64
  vector_count  = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%wavelet_scalar) .or. &
          .not. allocated(block_local(i)%wavelet_vector)) then
        error stop &
             "local_block_wavelet_statistics: storage missing"
     end if

     scalar_count = scalar_count + &
          int(size(block_local(i)%wavelet_scalar),int64)
     vector_count = vector_count + &
          int(size(block_local(i)%wavelet_vector),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_local(i)%wavelet_scalar)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_local(i)%wavelet_scalar))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_local(i)%wavelet_scalar**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_local(i)%wavelet_vector)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_local(i)%wavelet_vector))
     vector_moment(3) = vector_moment(3) + &
          sum(block_local(i)%wavelet_vector**2)

  end do

end subroutine local_block_wavelet_statistics


subroutine local_block_mean_field_statistics ( &
     scalar_count,vector_count,scalar_moment,vector_moment)
  ! Compute order-independent sol_mean inventory statistics over the
  ! ready final-owner block store.

  implicit none

  integer(int64), intent(out) :: scalar_count
  integer(int64), intent(out) :: vector_count

  real(dp), intent(out) :: scalar_moment(3)
  real(dp), intent(out) :: vector_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_mean_field_statistics: store is not ready"
  end if

  scalar_count  = 0_int64
  vector_count  = 0_int64
  scalar_moment = 0.0_dp
  vector_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%scalar_mean) .or. &
          .not. allocated(block_local(i)%vector_mean)) then
        error stop &
             "local_block_mean_field_statistics: storage missing"
     end if

     scalar_count = scalar_count + &
          int(size(block_local(i)%scalar_mean),int64)
     vector_count = vector_count + &
          int(size(block_local(i)%vector_mean),int64)

     scalar_moment(1) = scalar_moment(1) + &
          sum(block_local(i)%scalar_mean)
     scalar_moment(2) = scalar_moment(2) + &
          sum(abs(block_local(i)%scalar_mean))
     scalar_moment(3) = scalar_moment(3) + &
          sum(block_local(i)%scalar_mean**2)

     vector_moment(1) = vector_moment(1) + &
          sum(block_local(i)%vector_mean)
     vector_moment(2) = vector_moment(2) + &
          sum(abs(block_local(i)%vector_mean))
     vector_moment(3) = vector_moment(3) + &
          sum(block_local(i)%vector_mean**2)

  end do

end subroutine local_block_mean_field_statistics


subroutine local_block_turbulence_statistics ( &
     tke_count,wavelet_count,tke_moment,wavelet_moment)
  ! Compute order-independent tke and wav_tke inventories over the
  ! ready final-owner block store. Counts are zero when vert_diffuse is off.

  implicit none

  integer(int64), intent(out) :: tke_count
  integer(int64), intent(out) :: wavelet_count

  real(dp), intent(out) :: tke_moment(3)
  real(dp), intent(out) :: wavelet_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_turbulence_statistics: store is not ready"
  end if

  tke_count      = 0_int64
  wavelet_count  = 0_int64
  tke_moment     = 0.0_dp
  wavelet_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%tke) .or. &
          .not. allocated(block_local(i)%wavelet_tke)) then
        error stop &
             "local_block_turbulence_statistics: storage missing"
     end if

     tke_count = tke_count + int(size(block_local(i)%tke),int64)
     wavelet_count = wavelet_count + &
          int(size(block_local(i)%wavelet_tke),int64)

     tke_moment(1) = tke_moment(1) + sum(block_local(i)%tke)
     tke_moment(2) = tke_moment(2) + sum(abs(block_local(i)%tke))
     tke_moment(3) = tke_moment(3) + sum(block_local(i)%tke**2)

     wavelet_moment(1) = wavelet_moment(1) + &
          sum(block_local(i)%wavelet_tke)
     wavelet_moment(2) = wavelet_moment(2) + &
          sum(abs(block_local(i)%wavelet_tke))
     wavelet_moment(3) = wavelet_moment(3) + &
          sum(block_local(i)%wavelet_tke**2)

  end do

end subroutine local_block_turbulence_statistics


subroutine local_block_topography_statistics (value_count,value_moment)
  ! Compute an order-independent topography inventory over the ready
  ! final-owner block store.

  implicit none

  integer(int64), intent(out) :: value_count
  real(dp), intent(out) :: value_moment(3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_topography_statistics: store is not ready"
  end if

  value_count  = 0_int64
  value_moment = 0.0_dp

  do i = 1, size(block_local)

     if (.not. allocated(block_local(i)%topography)) then
        error stop &
             "local_block_topography_statistics: storage missing"
     end if

     value_count = value_count + &
          int(size(block_local(i)%topography),int64)

     value_moment(1) = value_moment(1) + &
          sum(block_local(i)%topography)
     value_moment(2) = value_moment(2) + &
          sum(abs(block_local(i)%topography))
     value_moment(3) = value_moment(3) + &
          sum(block_local(i)%topography**2)

  end do

end subroutine local_block_topography_statistics


subroutine source_block_boundary_route_statistics ( &
     link_count,storage_count,node_count)
  ! Inventory the field-independent boundary routes before migration.
  ! These routes are shared by the scalar, rank-one and rank-two
  ! Float_Field update_bdry overloads; field rank only repeats the route.

  implicit none

  integer(int64), intent(out) :: link_count(3)
  integer(int64), intent(out) :: storage_count(2)
  integer(int64), intent(out) :: node_count(2)

  integer :: i

  if (.not. allocated(block_source)) then
     error stop &
          "source_block_boundary_route_statistics: source unavailable"
  end if

  link_count    = 0_int64
  storage_count = 0_int64
  node_count    = 0_int64

  do i = 1, size(block_source)
     call accumulate_block_boundary_route_statistics( &
          block_source(i),link_count,storage_count,node_count)
  end do

end subroutine source_block_boundary_route_statistics


subroutine local_block_boundary_route_statistics ( &
     link_count,storage_count,node_count)
  ! Inventory the identical routes in the installed final-owner store.

  implicit none

  integer(int64), intent(out) :: link_count(3)
  integer(int64), intent(out) :: storage_count(2)
  integer(int64), intent(out) :: node_count(2)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_boundary_route_statistics: store is not ready"
  end if

  link_count    = 0_int64
  storage_count = 0_int64
  node_count    = 0_int64

  do i = 1, size(block_local)
     call accumulate_block_boundary_route_statistics( &
          block_local(i),link_count,storage_count,node_count)
  end do

end subroutine local_block_boundary_route_statistics


subroutine accumulate_block_boundary_route_statistics ( &
     block,link_count,storage_count,node_count)
  ! Validate one field-independent boundary catalogue. Link counters are
  ! inter-block, existing-domain and adaptive. Storage/node counters are
  ! ghost and compact-boundary. AT_EDGE payloads contain EDGE values per
  ! stored node and require the sign convention implemented by update_bdry.

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int64), intent(inout) :: link_count(3)
  integer(int64), intent(inout) :: storage_count(2)
  integer(int64), intent(inout) :: node_count(2)

  integer :: expected_start
  integer :: ghost_id
  integer :: i
  integer :: storage_id

  if (block%scalar_mult /= 1 .or. block%vector_mult /= EDGE) then
     error stop &
          "accumulate_block_boundary_route_statistics: multiplier"
  end if

  expected_start = 0
  do i = 1, size(block%ghost_storage)
     if (block%ghost_storage(i)%local_start /= expected_start .or. &
          block%ghost_storage(i)%n_node <= 0) then
        error stop &
             "accumulate_block_boundary_route_statistics: ghost layout"
     end if
     if (block%ghost_storage(i)%source_domain < 0 .or. &
          block%ghost_storage(i)%source_patch < 0 .or. &
          block%ghost_storage(i)%source_block < 1 .or. &
          block%ghost_storage(i)%source_block_id < 0 .or. &
          block%ghost_storage(i)%source_owner < 0) then
        error stop &
             "accumulate_block_boundary_route_statistics: ghost source"
     end if
     expected_start = expected_start + block%ghost_storage(i)%n_node
  end do

  if (expected_start /= size(block%ghost_node)) then
     error stop &
          "accumulate_block_boundary_route_statistics: ghost extent"
  end if

  expected_start = 0
  do i = 1, size(block%bdry_storage)
     if (block%bdry_storage(i)%local_start /= expected_start .or. &
          block%bdry_storage(i)%n_node <= 0 .or. &
          block%bdry_storage(i)%source_bdry < 1) then
        error stop &
             "accumulate_block_boundary_route_statistics: boundary layout"
     end if
     expected_start = expected_start + block%bdry_storage(i)%n_node
  end do

  if (expected_start /= size(block%bdry_node)) then
     error stop &
          "accumulate_block_boundary_route_statistics: boundary extent"
  end if

  do i = 1, size(block%block_bdry)
     select case (block%block_bdry(i)%class)

     case (NGB_BLOCK)
        ghost_id = block%block_bdry(i)%ghost_id
        if (ghost_id < 1 .or. ghost_id > size(block%ghost_storage)) then
           error stop &
                "accumulate_block_boundary_route_statistics: ghost ID"
        end if
        if (block%ghost_storage(ghost_id)%source_patch /= &
             block%block_bdry(i)%neigh_patch .or. &
             block%ghost_storage(ghost_id)%source_block /= &
             block%block_bdry(i)%source_block .or. &
             block%ghost_storage(ghost_id)%source_block_id /= &
             block%block_bdry(i)%source_block_id .or. &
             block%ghost_storage(ghost_id)%source_owner /= &
             block%block_bdry(i)%source_owner) then
           error stop &
                "accumulate_block_boundary_route_statistics: ghost route"
        end if
        link_count(1) = link_count(1) + 1_int64

     case (NGB_DOMAIN, NGB_ADAPT)
        storage_id = block%block_bdry(i)%storage_id
        if (storage_id < 1 .or. &
             storage_id > size(block%bdry_storage)) then
           error stop &
                "accumulate_block_boundary_route_statistics: storage ID"
        end if
        if (block%bdry_storage(storage_id)%source_bdry /= &
             block%block_bdry(i)%source_bdry) then
           error stop &
                "accumulate_block_boundary_route_statistics: boundary route"
        end if
        if (block%block_bdry(i)%class == NGB_DOMAIN) then
           link_count(2) = link_count(2) + 1_int64
        else
           link_count(3) = link_count(3) + 1_int64
        end if

     case default
        error stop &
             "accumulate_block_boundary_route_statistics: link class"

     end select
  end do

  storage_count(1) = storage_count(1) + &
       int(size(block%ghost_storage),int64)
  storage_count(2) = storage_count(2) + &
       int(size(block%bdry_storage),int64)
  node_count(1) = node_count(1) + int(size(block%ghost_node),int64)
  node_count(2) = node_count(2) + int(size(block%bdry_node),int64)

end subroutine accumulate_block_boundary_route_statistics


subroutine source_block_ghost_source_statistics ( &
     ghost_count,value_count,source_sum)
  ! Inventory source-block addresses before migration. The address is
  ! independent of Float_Field rank and position; position determines
  ! only whether one or EDGE values are packed for each source node.

  implicit none

  integer(int64), intent(out) :: ghost_count
  integer(int64), intent(out) :: value_count(2)
  integer(int64), intent(out) :: source_sum(5)

  integer :: i

  if (.not. allocated(block_source)) then
     error stop &
          "source_block_ghost_source_statistics: source unavailable"
  end if

  ghost_count = 0_int64
  value_count = 0_int64
  source_sum  = 0_int64

  do i = 1, size(block_source)
     call accumulate_block_ghost_source_statistics( &
          block_source(i),ghost_count,value_count,source_sum)
  end do

end subroutine source_block_ghost_source_statistics


subroutine local_block_ghost_source_statistics ( &
     ghost_count,value_count,source_sum)
  ! Inventory source-block addresses after final-owner installation.

  implicit none

  integer(int64), intent(out) :: ghost_count
  integer(int64), intent(out) :: value_count(2)
  integer(int64), intent(out) :: source_sum(5)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_ghost_source_statistics: store is not ready"
  end if

  ghost_count = 0_int64
  value_count = 0_int64
  source_sum  = 0_int64

  do i = 1, size(block_local)
     call accumulate_block_ghost_source_statistics( &
          block_local(i),ghost_count,value_count,source_sum)
  end do

end subroutine local_block_ghost_source_statistics


subroutine accumulate_block_ghost_source_statistics ( &
     block,ghost_count,value_count,source_sum)

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int64), intent(inout) :: ghost_count
  integer(int64), intent(inout) :: value_count(2)
  integer(int64), intent(inout) :: source_sum(5)

  integer :: i

  do i = 1, size(block%ghost_storage)
     if (block%ghost_storage(i)%source_local_patch < 0 .or. &
          block%ghost_storage(i)%n_node /= PATCH_SIZE**2) then
        error stop &
             "accumulate_block_ghost_source_statistics: source address"
     end if

     ghost_count = ghost_count + 1_int64
     value_count(1) = value_count(1) + &
          int(block%ghost_storage(i)%n_node,int64)
     value_count(2) = value_count(2) + &
          int(EDGE*block%ghost_storage(i)%n_node,int64)

     source_sum = source_sum + int([ &
          block%ghost_storage(i)%source_domain, &
          block%ghost_storage(i)%source_block, &
          block%ghost_storage(i)%source_block_id, &
          block%ghost_storage(i)%source_owner, &
          block%ghost_storage(i)%source_local_patch ],int64)
  end do

end subroutine accumulate_block_ghost_source_statistics


subroutine validate_local_block_ghost_sources ( &
     catalog_patch_count,catalog_owner,catalog_id,catalog_domain)
  ! Cross-check each final-owner ghost address against the replicated
  ! block catalogue and the globally assembled source-block patch counts.

  implicit none

  integer, intent(in) :: catalog_patch_count(:)
  integer, intent(in) :: catalog_owner(:)
  integer, intent(in) :: catalog_id(:)
  integer, intent(in) :: catalog_domain(:)

  integer :: b
  integer :: i
  integer :: source_block

  if (.not. local_block_store_ready()) then
     error stop "validate_local_block_ghost_sources: store is not ready"
  end if

  if (size(catalog_patch_count) /= size(catalog_owner) .or. &
       size(catalog_patch_count) /= size(catalog_id) .or. &
       size(catalog_patch_count) /= size(catalog_domain)) then
     error stop "validate_local_block_ghost_sources: catalogue extent"
  end if

  do b = 1, size(block_local)
     do i = 1, size(block_local(b)%ghost_storage)
        source_block = &
             block_local(b)%ghost_storage(i)%source_block

        if (source_block < 1 .or. &
             source_block > size(catalog_patch_count)) then
           error stop &
                "validate_local_block_ghost_sources: source block"
        end if

        if (block_local(b)%ghost_storage(i)%source_local_patch < 0 .or. &
             block_local(b)%ghost_storage(i)%source_local_patch >= &
             catalog_patch_count(source_block)) then
           error stop &
                "validate_local_block_ghost_sources: local patch"
        end if

        if (block_local(b)%ghost_storage(i)%source_owner /= &
             catalog_owner(source_block) .or. &
             block_local(b)%ghost_storage(i)%source_block_id /= &
             catalog_id(source_block) .or. &
             block_local(b)%ghost_storage(i)%source_domain /= &
             catalog_domain(source_block)) then
           error stop &
                "validate_local_block_ghost_sources: source identity"
        end if
     end do
  end do

end subroutine validate_local_block_ghost_sources


subroutine get_local_block_ghost_requests ( &
     source_block,source_local_patch,source_owner, &
     destination_block,destination_ghost)
  ! Export one field-independent request for every ghost record owned by
  ! this rank. Float_Field rank and position affect later payload packing,
  ! not this source/destination manifest.

  implicit none

  integer, allocatable, intent(out) :: source_block(:)
  integer, allocatable, intent(out) :: source_local_patch(:)
  integer, allocatable, intent(out) :: source_owner(:)
  integer, allocatable, intent(out) :: destination_block(:)
  integer, allocatable, intent(out) :: destination_ghost(:)

  integer :: b
  integer :: g
  integer :: i
  integer :: n_request

  if (.not. local_block_store_ready()) then
     error stop "get_local_block_ghost_requests: store is not ready"
  end if

  n_request = 0
  do b = 1, size(block_local)
     n_request = n_request + size(block_local(b)%ghost_storage)
  end do

  allocate(source_block(n_request))
  allocate(source_local_patch(n_request))
  allocate(source_owner(n_request))
  allocate(destination_block(n_request))
  allocate(destination_ghost(n_request))

  i = 0
  do b = 1, size(block_local)
     do g = 1, size(block_local(b)%ghost_storage)
        i = i + 1
        source_block(i) = &
             block_local(b)%ghost_storage(g)%source_block
        source_local_patch(i) = &
             block_local(b)%ghost_storage(g)%source_local_patch
        source_owner(i) = &
             block_local(b)%ghost_storage(g)%source_owner
        destination_block(i) = block_local_catalog_index(b)
        destination_ghost(i) = g
     end do
  end do

  if (i /= n_request) then
     error stop "get_local_block_ghost_requests: request count"
  end if

end subroutine get_local_block_ghost_requests


integer function local_block_patch_count (catalog_index) result(n_patch)
  ! Return the compact interior-patch count for a locally owned block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop "local_block_patch_count: block is not local"
  end if

  n_patch = size(block_local(local_index)%patch)

end function local_block_patch_count


integer function local_block_boundary_count (catalog_index) &
     result(n_boundary)
  ! Return the compact boundary-record count for a locally owned block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop "local_block_boundary_count: block is not local"
  end if

  n_boundary = size(block_local(local_index)%bdry_storage)

end function local_block_boundary_count


subroutine get_local_block_boundary_source ( &
     catalog_index,boundary_index,source_bdry,elts_start,n_node)
  ! Return the authoritative Domain-boundary address represented by one
  ! compact boundary record in a locally owned final block.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  integer, intent(out) :: source_bdry
  integer, intent(out) :: elts_start
  integer, intent(out) :: n_node

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop "get_local_block_boundary_source: block is not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop "get_local_block_boundary_source: invalid boundary"
  end if

  source_bdry = block_local(local_index)% &
       bdry_storage(boundary_index)%source_bdry
  elts_start = block_local(local_index)% &
       bdry_storage(boundary_index)%elts_start
  n_node = block_local(local_index)% &
       bdry_storage(boundary_index)%n_node

  if (source_bdry < 0 .or. elts_start < 0 .or. n_node < 1) then
     error stop "get_local_block_boundary_source: invalid source"
  end if

end subroutine get_local_block_boundary_source


integer function local_block_ghost_count (catalog_index) result(n_ghost)
  ! Return the compact ghost-record count for a locally owned block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop "local_block_ghost_count: block is not local"
  end if

  n_ghost = size(block_local(local_index)%ghost_storage)

end function local_block_ghost_count


integer function local_block_scalar_patch_nvalue (catalog_index) &
     result(n_value)
  ! Number of scalar sol and wav_coeff values carried by one compact
  ! patch. The complete scalar-variable and level ranges are included.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_scalar_patch_nvalue: block is not local"
  end if

  n_value = 2 * block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2

end function local_block_scalar_patch_nvalue


integer function local_block_scalar_family_patch_nvalue ( &
     catalog_index) result(n_value)
  ! Number of values for one scalar payload family on one compact patch.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_scalar_family_patch_nvalue: block is not local"
  end if

  n_value = block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2

end function local_block_scalar_family_patch_nvalue


subroutine get_local_block_scalar_patch_family_values ( &
     catalog_index,local_patch,payload_family,value)
  ! Pack one scalar sol or wav_coeff family from one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: output_base
  integer :: patch_start
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_patch_family_values: block is not local"
  end if

  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_scalar_patch_family_values: invalid patch"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_scalar_patch_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_scalar_family_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_scalar_patch_family_values: output extent"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "get_local_block_scalar_patch_family_values: patch storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_node
        output_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           value(output_base+1:output_base+n_patch_value) = &
                block_local(local_index)%scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + &
                n_patch_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           value(output_base+1:output_base+n_patch_value) = &
                block_local(local_index)%wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + &
                n_patch_value)
        end select
     end do
  end do

end subroutine get_local_block_scalar_patch_family_values


subroutine get_local_block_scalar_patch_values ( &
     catalog_index,local_patch,value)
  ! Pack scalar sol followed by scalar wav_coeff for one compact patch.
  ! local_patch is the zero-based compact source-patch address.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: output_base
  integer :: patch_start
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_patch_values: block is not local"
  end if

  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_scalar_patch_values: invalid local patch"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_scalar_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_scalar_patch_values: output extent"
  end if

  n_node = size(block_local(local_index)%node)
  family_base = &
       block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * n_patch_value
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "get_local_block_scalar_patch_values: patch storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_node
        output_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + &
             n_patch_value)

        value( &
             family_base+output_base+1: &
             family_base+output_base+n_patch_value) = &
             block_local(local_index)%wavelet_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*patch_start + &
             n_patch_value)
     end do
  end do

end subroutine get_local_block_scalar_patch_values


subroutine set_local_block_scalar_patch_family_values ( &
     catalog_index,local_patch,payload_family,value)
  ! Install one scalar sol or wav_coeff family in one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: patch_start
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_scalar_patch_family_values: block not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "set_local_block_scalar_patch_family_values: invalid patch"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_scalar_patch_family_values: invalid family"
  end if
  if (size(value) /= &
       local_block_scalar_family_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_scalar_patch_family_values: input extent"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "set_local_block_scalar_patch_family_values: patch storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_node
        input_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           block_local(local_index)%scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + &
                n_patch_value) = &
                value(input_base+1:input_base+n_patch_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           block_local(local_index)%wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*patch_start + &
                n_patch_value) = &
                value(input_base+1:input_base+n_patch_value)
        end select
     end do
  end do

  if (payload_family == BLOCK_PAYLOAD_SOL) then
     call invalidate_local_block_hydrostatic_block(local_index)
  end if

end subroutine set_local_block_scalar_patch_family_values


subroutine set_local_block_scalar_patch_values ( &
     catalog_index,local_patch,value)
  ! Install scalar sol followed by wav_coeff in one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(in) :: value(:)

  integer :: n_family

  if (size(value) /= local_block_scalar_patch_nvalue( &
       catalog_index)) then
     error stop "set_local_block_scalar_patch_values: input extent"
  end if

  n_family = local_block_scalar_family_patch_nvalue(catalog_index)
  call set_local_block_scalar_patch_family_values( &
       catalog_index,local_patch,BLOCK_PAYLOAD_SOL, &
       value(1:n_family))
  call set_local_block_scalar_patch_family_values( &
       catalog_index,local_patch,BLOCK_PAYLOAD_WAV_COEFF, &
       value(n_family+1:2*n_family))

end subroutine set_local_block_scalar_patch_values


subroutine fill_local_block_scalar_patch_values (value)
  ! Fill scalar sol and wav_coeff patch values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_patch_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%scalar = value
     block_local(b)%wavelet_scalar = value
  end do

  call invalidate_local_block_hydrostatic_state

end subroutine fill_local_block_scalar_patch_values


subroutine fill_local_block_scalar_patch_family_values ( &
     payload_family,value)
  ! Fill one scalar patch payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_patch_family_values: store not ready"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_scalar_patch_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%scalar = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%wavelet_scalar = value
     end select
  end do

  if (payload_family == BLOCK_PAYLOAD_SOL) then
     call invalidate_local_block_hydrostatic_state
  end if

end subroutine fill_local_block_scalar_patch_family_values


integer function local_block_scalar_boundary_nvalue ( &
     catalog_index,boundary_index) result(n_value)
  ! Number of scalar sol and wav_coeff values in one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_scalar_boundary_nvalue: block is not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "local_block_scalar_boundary_nvalue: invalid boundary"
  end if

  n_value = 2 * block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * &
       block_local(local_index)%scalar_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node

end function local_block_scalar_boundary_nvalue


integer function local_block_scalar_family_boundary_nvalue ( &
     catalog_index,boundary_index) result(n_value)
  ! Number of values for one scalar family in one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_scalar_family_boundary_nvalue: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "local_block_scalar_family_boundary_nvalue: invalid boundary"
  end if

  n_value = block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * &
       block_local(local_index)%scalar_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node

end function local_block_scalar_family_boundary_nvalue


subroutine get_local_block_scalar_boundary_values ( &
     catalog_index,boundary_index,value)
  ! Read scalar sol followed by wav_coeff from one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  real(dp), intent(out) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: family_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value
  integer :: output_base
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_boundary_values: block is not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "get_local_block_scalar_boundary_values: invalid boundary"
  end if
  if (size(value) /= local_block_scalar_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "get_local_block_scalar_boundary_values: output extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%scalar_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  family_base = block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level*n_boundary_value
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "get_local_block_scalar_boundary_values: boundary storage"
  end if

  do scalar_slot = 1, block_local(local_index)%n_scalar_variable
     do level_slot = 1, block_local(local_index)%n_field_level
        field_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             block_local(local_index)%scalar_mult*n_boundary_node
        output_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             n_boundary_value

        value(output_base+1:output_base+n_boundary_value) = &
             block_local(local_index)%bdry_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*boundary_start + &
             n_boundary_value)
        value(family_base+output_base+1: &
             family_base+output_base+n_boundary_value) = &
             block_local(local_index)%bdry_wavelet_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*boundary_start + &
             n_boundary_value)
     end do
  end do

end subroutine get_local_block_scalar_boundary_values


subroutine get_local_block_scalar_boundary_family_values ( &
     catalog_index,boundary_index,payload_family,value)
  ! Read one scalar sol or wav_coeff family from a compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value
  integer :: output_base
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_boundary_family_values: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "get_local_block_scalar_boundary_family_values: invalid boundary"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_scalar_boundary_family_values: invalid family"
  end if
  if (size(value) /= local_block_scalar_family_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "get_local_block_scalar_boundary_family_values: output extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%scalar_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "get_local_block_scalar_boundary_family_values: storage"
  end if

  do scalar_slot = 1, block_local(local_index)%n_scalar_variable
     do level_slot = 1, block_local(local_index)%n_field_level
        field_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             block_local(local_index)%scalar_mult*n_boundary_node
        output_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             n_boundary_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           value(output_base+1:output_base+n_boundary_value) = &
                block_local(local_index)%bdry_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + &
                n_boundary_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           value(output_base+1:output_base+n_boundary_value) = &
                block_local(local_index)%bdry_wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + &
                n_boundary_value)
        end select
     end do
  end do

end subroutine get_local_block_scalar_boundary_family_values


subroutine set_local_block_scalar_boundary_family_values ( &
     catalog_index,boundary_index,payload_family,value)
  ! Install one scalar sol or wav_coeff family in a compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_scalar_boundary_family_values: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "set_local_block_scalar_boundary_family_values: invalid boundary"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_scalar_boundary_family_values: invalid family"
  end if
  if (size(value) /= local_block_scalar_family_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "set_local_block_scalar_boundary_family_values: input extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%scalar_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "set_local_block_scalar_boundary_family_values: storage"
  end if

  do scalar_slot = 1, block_local(local_index)%n_scalar_variable
     do level_slot = 1, block_local(local_index)%n_field_level
        field_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             block_local(local_index)%scalar_mult*n_boundary_node
        input_base = ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + level_slot-1) * &
             n_boundary_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           block_local(local_index)%bdry_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + &
                n_boundary_value) = &
                value(input_base+1:input_base+n_boundary_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           block_local(local_index)%bdry_wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*boundary_start + &
                n_boundary_value) = &
                value(input_base+1:input_base+n_boundary_value)
        end select
     end do
  end do

end subroutine set_local_block_scalar_boundary_family_values


subroutine fill_local_block_scalar_boundary_values (value)
  ! Fill scalar sol and wav_coeff boundary values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_boundary_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%bdry_scalar = value
     block_local(b)%bdry_wavelet_scalar = value
  end do

end subroutine fill_local_block_scalar_boundary_values


subroutine fill_local_block_scalar_boundary_family_values ( &
     payload_family,value)
  ! Fill one scalar boundary payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_boundary_family_values: store not ready"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_scalar_boundary_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%bdry_scalar = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%bdry_wavelet_scalar = value
     end select
  end do

end subroutine fill_local_block_scalar_boundary_family_values


subroutine get_local_block_scalar_ghost_values ( &
     catalog_index,ghost_index,value)
  ! Read scalar sol followed by scalar wav_coeff from one compact ghost.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: ghost_start
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: output_base
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_ghost_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "get_local_block_scalar_ghost_values: invalid ghost"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_scalar_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_scalar_ghost_values: output extent"
  end if

  n_ghost_node = size(block_local(local_index)%ghost_node)
  family_base = &
       block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * n_patch_value
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "get_local_block_scalar_ghost_values: ghost storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_ghost_node
        output_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%ghost_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + &
             n_patch_value)

        value( &
             family_base+output_base+1: &
             family_base+output_base+n_patch_value) = &
             block_local(local_index)%ghost_wavelet_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + &
             n_patch_value)
     end do
  end do

end subroutine get_local_block_scalar_ghost_values


subroutine get_local_block_scalar_ghost_family_values ( &
     catalog_index,ghost_index,payload_family,value)
  ! Read one scalar sol or wav_coeff family from one compact ghost patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: ghost_start
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: output_base
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_scalar_ghost_family_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "get_local_block_scalar_ghost_family_values: invalid ghost"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_scalar_ghost_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_scalar_family_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_scalar_ghost_family_values: output extent"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  n_ghost_node = size(block_local(local_index)%ghost_node)
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "get_local_block_scalar_ghost_family_values: ghost storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_ghost_node
        output_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           value(output_base+1:output_base+n_patch_value) = &
                block_local(local_index)%ghost_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + &
                n_patch_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           value(output_base+1:output_base+n_patch_value) = &
                block_local(local_index)%ghost_wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + &
                n_patch_value)
        end select
     end do
  end do

end subroutine get_local_block_scalar_ghost_family_values


subroutine set_local_block_scalar_ghost_family_values ( &
     catalog_index,ghost_index,payload_family,value)
  ! Install one scalar sol or wav_coeff family in one compact ghost patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: ghost_start
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_scalar_ghost_family_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "set_local_block_scalar_ghost_family_values: invalid ghost"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_scalar_ghost_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_scalar_family_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_scalar_ghost_family_values: input extent"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  n_ghost_node = size(block_local(local_index)%ghost_node)
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "set_local_block_scalar_ghost_family_values: ghost storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_ghost_node
        input_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        select case (payload_family)
        case (BLOCK_PAYLOAD_SOL)
           block_local(local_index)%ghost_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + &
                n_patch_value) = &
                value(input_base+1:input_base+n_patch_value)
        case (BLOCK_PAYLOAD_WAV_COEFF)
           block_local(local_index)%ghost_wavelet_scalar( &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + 1: &
                field_base + &
                block_local(local_index)%scalar_mult*ghost_start + &
                n_patch_value) = &
                value(input_base+1:input_base+n_patch_value)
        end select
     end do
  end do

end subroutine set_local_block_scalar_ghost_family_values


subroutine set_local_block_scalar_ghost_values ( &
     catalog_index,ghost_index,value)
  ! Install scalar sol and scalar wav_coeff into compact ghost storage.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: ghost_start
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: scalar_slot

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_scalar_ghost_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "set_local_block_scalar_ghost_values: invalid ghost"
  end if

  n_patch_value = &
       block_local(local_index)%scalar_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_scalar_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_scalar_ghost_values: input extent"
  end if

  n_ghost_node = size(block_local(local_index)%ghost_node)
  family_base = &
       block_local(local_index)%n_scalar_variable * &
       block_local(local_index)%n_field_level * n_patch_value
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "set_local_block_scalar_ghost_values: ghost storage"
  end if

  do scalar_slot = 1, &
       block_local(local_index)%n_scalar_variable
     do level_slot = 1, &
          block_local(local_index)%n_field_level
        field_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * &
             block_local(local_index)%scalar_mult*n_ghost_node
        input_base = &
             ((scalar_slot-1)* &
             block_local(local_index)%n_field_level + &
             level_slot-1) * n_patch_value

        block_local(local_index)%ghost_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + &
             n_patch_value) = &
             value(input_base+1:input_base+n_patch_value)

        block_local(local_index)%ghost_wavelet_scalar( &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%scalar_mult*ghost_start + &
             n_patch_value) = value( &
             family_base+input_base+1: &
             family_base+input_base+n_patch_value)
     end do
  end do

end subroutine set_local_block_scalar_ghost_values


subroutine fill_local_block_scalar_ghost_values (value)
  ! Fill scalar sol and wav_coeff ghost values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_ghost_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%ghost_scalar = value
     block_local(b)%ghost_wavelet_scalar = value
  end do

end subroutine fill_local_block_scalar_ghost_values


subroutine fill_local_block_scalar_ghost_family_values ( &
     payload_family,value)
  ! Fill one scalar ghost payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_scalar_ghost_family_values: store not ready"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_scalar_ghost_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%ghost_scalar = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%ghost_wavelet_scalar = value
     end select
  end do

end subroutine fill_local_block_scalar_ghost_family_values


integer function local_block_vector_patch_nvalue (catalog_index) &
     result(n_value)
  ! Number of vector sol and wav_coeff values carried by one compact
  ! patch. Every field level and vector component is included.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_vector_patch_nvalue: block is not local"
  end if

  n_value = 2 * block_local(local_index)%n_field_level * &
       block_local(local_index)%vector_mult * PATCH_SIZE**2

end function local_block_vector_patch_nvalue


integer function local_block_vector_family_patch_nvalue ( &
     catalog_index) result(n_value)
  ! Number of values for one vector payload family on one compact patch.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_vector_family_patch_nvalue: block is not local"
  end if

  n_value = block_local(local_index)%n_field_level * &
       block_local(local_index)%vector_mult * PATCH_SIZE**2

end function local_block_vector_family_patch_nvalue


subroutine get_local_block_vector_patch_family_values ( &
     catalog_index,local_patch,payload_family,value)
  ! Pack one vector sol or wav_coeff family from one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: output_base
  integer :: patch_start

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_patch_family_values: block is not local"
  end if

  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_vector_patch_family_values: invalid patch"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_vector_patch_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_vector_family_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_vector_patch_family_values: output extent"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "get_local_block_vector_patch_family_values: patch storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_node
     output_base = (level_slot-1)*n_patch_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%vector( &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + &
             n_patch_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + &
             n_patch_value)
     end select
  end do

end subroutine get_local_block_vector_patch_family_values


subroutine get_local_block_vector_patch_values ( &
     catalog_index,local_patch,value)
  ! Pack vector sol followed by vector wav_coeff for one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: output_base
  integer :: patch_start

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_patch_values: block is not local"
  end if

  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_vector_patch_values: invalid local patch"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_vector_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_vector_patch_values: output extent"
  end if

  n_node = size(block_local(local_index)%node)
  family_base = &
       block_local(local_index)%n_field_level*n_patch_value
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "get_local_block_vector_patch_values: patch storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_node
     output_base = (level_slot-1)*n_patch_value

     value(output_base+1:output_base+n_patch_value) = &
          block_local(local_index)%vector( &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + &
          n_patch_value)

     value( &
          family_base+output_base+1: &
          family_base+output_base+n_patch_value) = &
          block_local(local_index)%wavelet_vector( &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*patch_start + &
          n_patch_value)
  end do

end subroutine get_local_block_vector_patch_values


subroutine set_local_block_vector_patch_family_values ( &
     catalog_index,local_patch,payload_family,value)
  ! Install one vector sol or wav_coeff family in one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_node
  integer :: n_patch_value
  integer :: patch_start

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_vector_patch_family_values: block not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "set_local_block_vector_patch_family_values: invalid patch"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_vector_patch_family_values: invalid family"
  end if
  if (size(value) /= &
       local_block_vector_family_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_vector_patch_family_values: input extent"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start

  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "set_local_block_vector_patch_family_values: patch storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_node
     input_base = (level_slot-1)*n_patch_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(local_index)%vector( &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + &
             n_patch_value) = &
             value(input_base+1:input_base+n_patch_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(local_index)%wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*patch_start + &
             n_patch_value) = &
             value(input_base+1:input_base+n_patch_value)
     end select
  end do

end subroutine set_local_block_vector_patch_family_values


subroutine set_local_block_vector_patch_values ( &
     catalog_index,local_patch,value)
  ! Install vector sol followed by wav_coeff in one compact patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch
  real(dp), intent(in) :: value(:)

  integer :: n_family

  if (size(value) /= local_block_vector_patch_nvalue( &
       catalog_index)) then
     error stop "set_local_block_vector_patch_values: input extent"
  end if

  n_family = local_block_vector_family_patch_nvalue(catalog_index)
  call set_local_block_vector_patch_family_values( &
       catalog_index,local_patch,BLOCK_PAYLOAD_SOL, &
       value(1:n_family))
  call set_local_block_vector_patch_family_values( &
       catalog_index,local_patch,BLOCK_PAYLOAD_WAV_COEFF, &
       value(n_family+1:2*n_family))

end subroutine set_local_block_vector_patch_values


subroutine fill_local_block_vector_patch_values (value)
  ! Fill vector sol and wav_coeff patch values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_patch_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%vector = value
     block_local(b)%wavelet_vector = value
  end do

end subroutine fill_local_block_vector_patch_values


subroutine fill_local_block_vector_patch_family_values ( &
     payload_family,value)
  ! Fill one vector patch payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_patch_family_values: store not ready"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_vector_patch_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%vector = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%wavelet_vector = value
     end select
  end do

end subroutine fill_local_block_vector_patch_family_values


integer function local_block_vector_boundary_nvalue ( &
     catalog_index,boundary_index) result(n_value)
  ! Number of vector sol and wav_coeff values in one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_vector_boundary_nvalue: block is not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "local_block_vector_boundary_nvalue: invalid boundary"
  end if

  n_value = 2 * block_local(local_index)%n_field_level * &
       block_local(local_index)%vector_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node

end function local_block_vector_boundary_nvalue


integer function local_block_vector_family_boundary_nvalue ( &
     catalog_index,boundary_index) result(n_value)
  ! Number of values for one vector family in one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index

  integer :: local_index

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_vector_family_boundary_nvalue: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "local_block_vector_family_boundary_nvalue: invalid boundary"
  end if

  n_value = block_local(local_index)%n_field_level * &
       block_local(local_index)%vector_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node

end function local_block_vector_family_boundary_nvalue


subroutine get_local_block_vector_boundary_values ( &
     catalog_index,boundary_index,value)
  ! Read vector sol followed by wav_coeff from one compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  real(dp), intent(out) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: family_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value
  integer :: output_base

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_boundary_values: block is not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "get_local_block_vector_boundary_values: invalid boundary"
  end if
  if (size(value) /= local_block_vector_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "get_local_block_vector_boundary_values: output extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%vector_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  family_base = block_local(local_index)%n_field_level*n_boundary_value
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "get_local_block_vector_boundary_values: boundary storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_boundary_node
     output_base = (level_slot-1)*n_boundary_value

     value(output_base+1:output_base+n_boundary_value) = &
          block_local(local_index)%bdry_vector( &
          field_base + &
          block_local(local_index)%vector_mult*boundary_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*boundary_start + &
          n_boundary_value)
     value(family_base+output_base+1: &
          family_base+output_base+n_boundary_value) = &
          block_local(local_index)%bdry_wavelet_vector( &
          field_base + &
          block_local(local_index)%vector_mult*boundary_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*boundary_start + &
          n_boundary_value)
  end do

end subroutine get_local_block_vector_boundary_values


subroutine get_local_block_vector_boundary_family_values ( &
     catalog_index,boundary_index,payload_family,value)
  ! Read one vector sol or wav_coeff family from a compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value
  integer :: output_base

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_boundary_family_values: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "get_local_block_vector_boundary_family_values: invalid boundary"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_vector_boundary_family_values: invalid family"
  end if
  if (size(value) /= local_block_vector_family_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "get_local_block_vector_boundary_family_values: output extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%vector_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "get_local_block_vector_boundary_family_values: storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_boundary_node
     output_base = (level_slot-1)*n_boundary_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        value(output_base+1:output_base+n_boundary_value) = &
             block_local(local_index)%bdry_vector( &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + &
             n_boundary_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        value(output_base+1:output_base+n_boundary_value) = &
             block_local(local_index)%bdry_wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + &
             n_boundary_value)
     end select
  end do

end subroutine get_local_block_vector_boundary_family_values


subroutine set_local_block_vector_boundary_family_values ( &
     catalog_index,boundary_index,payload_family,value)
  ! Install one vector sol or wav_coeff family in a compact boundary.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: boundary_index
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: boundary_start
  integer :: field_base
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_boundary_node
  integer :: n_boundary_value

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_vector_boundary_family_values: block not local"
  end if
  if (boundary_index < 1 .or. &
       boundary_index > size(block_local(local_index)%bdry_storage)) then
     error stop &
          "set_local_block_vector_boundary_family_values: invalid boundary"
  end if
  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_vector_boundary_family_values: invalid family"
  end if
  if (size(value) /= local_block_vector_family_boundary_nvalue( &
       catalog_index,boundary_index)) then
     error stop &
          "set_local_block_vector_boundary_family_values: input extent"
  end if

  n_boundary_node = size(block_local(local_index)%bdry_node)
  n_boundary_value = block_local(local_index)%vector_mult * &
       block_local(local_index)%bdry_storage(boundary_index)%n_node
  boundary_start = block_local(local_index)% &
       bdry_storage(boundary_index)%local_start

  if (boundary_start < 0 .or. &
       boundary_start + &
       block_local(local_index)%bdry_storage(boundary_index)%n_node > &
       n_boundary_node) then
     error stop &
          "set_local_block_vector_boundary_family_values: storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_boundary_node
     input_base = (level_slot-1)*n_boundary_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(local_index)%bdry_vector( &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + &
             n_boundary_value) = &
             value(input_base+1:input_base+n_boundary_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(local_index)%bdry_wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*boundary_start + &
             n_boundary_value) = &
             value(input_base+1:input_base+n_boundary_value)
     end select
  end do

end subroutine set_local_block_vector_boundary_family_values


subroutine fill_local_block_vector_boundary_values (value)
  ! Fill vector sol and wav_coeff boundary values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_boundary_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%bdry_vector = value
     block_local(b)%bdry_wavelet_vector = value
  end do

end subroutine fill_local_block_vector_boundary_values


subroutine fill_local_block_vector_boundary_family_values ( &
     payload_family,value)
  ! Fill one vector boundary payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_boundary_family_values: store not ready"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_vector_boundary_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%bdry_vector = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%bdry_wavelet_vector = value
     end select
  end do

end subroutine fill_local_block_vector_boundary_family_values


subroutine get_local_block_vector_ghost_values ( &
     catalog_index,ghost_index,value)
  ! Read vector sol followed by vector wav_coeff from one compact ghost.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: ghost_start
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: output_base

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_ghost_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "get_local_block_vector_ghost_values: invalid ghost"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_vector_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_vector_ghost_values: output extent"
  end if

  n_ghost_node = size(block_local(local_index)%ghost_node)
  family_base = &
       block_local(local_index)%n_field_level*n_patch_value
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "get_local_block_vector_ghost_values: ghost storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_ghost_node
     output_base = (level_slot-1)*n_patch_value

     value(output_base+1:output_base+n_patch_value) = &
          block_local(local_index)%ghost_vector( &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + &
          n_patch_value)

     value( &
          family_base+output_base+1: &
          family_base+output_base+n_patch_value) = &
          block_local(local_index)%ghost_wavelet_vector( &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + &
          n_patch_value)
  end do

end subroutine get_local_block_vector_ghost_values


subroutine get_local_block_vector_ghost_family_values ( &
     catalog_index,ghost_index,payload_family,value)
  ! Read one vector sol or wav_coeff family from one compact ghost patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  integer, intent(in) :: payload_family
  real(dp), intent(out) :: value(:)

  integer :: field_base
  integer :: ghost_start
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value
  integer :: output_base

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_vector_ghost_family_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "get_local_block_vector_ghost_family_values: invalid ghost"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "get_local_block_vector_ghost_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_vector_family_patch_nvalue(catalog_index)) then
     error stop &
          "get_local_block_vector_ghost_family_values: output extent"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  n_ghost_node = size(block_local(local_index)%ghost_node)
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "get_local_block_vector_ghost_family_values: ghost storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_ghost_node
     output_base = (level_slot-1)*n_patch_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%ghost_vector( &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + &
             n_patch_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        value(output_base+1:output_base+n_patch_value) = &
             block_local(local_index)%ghost_wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + &
             n_patch_value)
     end select
  end do

end subroutine get_local_block_vector_ghost_family_values


subroutine set_local_block_vector_ghost_family_values ( &
     catalog_index,ghost_index,payload_family,value)
  ! Install one vector sol or wav_coeff family in one compact ghost patch.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: ghost_start
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_vector_ghost_family_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "set_local_block_vector_ghost_family_values: invalid ghost"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "set_local_block_vector_ghost_family_values: invalid family"
  end if

  if (size(value) /= &
       local_block_vector_family_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_vector_ghost_family_values: input extent"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  n_ghost_node = size(block_local(local_index)%ghost_node)
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "set_local_block_vector_ghost_family_values: ghost storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_ghost_node
     input_base = (level_slot-1)*n_patch_value

     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(local_index)%ghost_vector( &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + &
             n_patch_value) = &
             value(input_base+1:input_base+n_patch_value)
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(local_index)%ghost_wavelet_vector( &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + 1: &
             field_base + &
             block_local(local_index)%vector_mult*ghost_start + &
             n_patch_value) = &
             value(input_base+1:input_base+n_patch_value)
     end select
  end do

end subroutine set_local_block_vector_ghost_family_values


subroutine set_local_block_vector_ghost_values ( &
     catalog_index,ghost_index,value)
  ! Install vector sol and vector wav_coeff into compact ghost storage.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: ghost_index
  real(dp), intent(in) :: value(:)

  integer :: field_base
  integer :: family_base
  integer :: ghost_start
  integer :: input_base
  integer :: local_index
  integer :: level_slot
  integer :: n_ghost_node
  integer :: n_patch_value

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "set_local_block_vector_ghost_values: block is not local"
  end if

  if (ghost_index < 1 .or. &
       ghost_index > size(block_local(local_index)%ghost_storage)) then
     error stop &
          "set_local_block_vector_ghost_values: invalid ghost"
  end if

  n_patch_value = &
       block_local(local_index)%vector_mult * PATCH_SIZE**2
  if (size(value) /= &
       local_block_vector_patch_nvalue(catalog_index)) then
     error stop &
          "set_local_block_vector_ghost_values: input extent"
  end if

  n_ghost_node = size(block_local(local_index)%ghost_node)
  family_base = &
       block_local(local_index)%n_field_level*n_patch_value
  ghost_start = block_local(local_index)% &
       ghost_storage(ghost_index)%local_start

  if (ghost_start < 0 .or. &
       ghost_start+PATCH_SIZE**2 > n_ghost_node) then
     error stop &
          "set_local_block_vector_ghost_values: ghost storage"
  end if

  do level_slot = 1, block_local(local_index)%n_field_level
     field_base = (level_slot-1) * &
          block_local(local_index)%vector_mult*n_ghost_node
     input_base = (level_slot-1)*n_patch_value

     block_local(local_index)%ghost_vector( &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + &
          n_patch_value) = value(input_base+1:input_base+n_patch_value)

     block_local(local_index)%ghost_wavelet_vector( &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + 1: &
          field_base + &
          block_local(local_index)%vector_mult*ghost_start + &
          n_patch_value) = value( &
          family_base+input_base+1: &
          family_base+input_base+n_patch_value)
  end do

end subroutine set_local_block_vector_ghost_values


subroutine fill_local_block_vector_ghost_values (value)
  ! Fill vector sol and wav_coeff ghost values in the local block store.

  implicit none

  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_ghost_values: store is not ready"
  end if

  do b = 1, size(block_local)
     block_local(b)%ghost_vector = value
     block_local(b)%ghost_wavelet_vector = value
  end do

end subroutine fill_local_block_vector_ghost_values


subroutine fill_local_block_vector_ghost_family_values ( &
     payload_family,value)
  ! Fill one vector ghost payload family in the local block store.

  implicit none

  integer, intent(in) :: payload_family
  real(dp), intent(in) :: value

  integer :: b

  if (.not. local_block_store_ready()) then
     error stop &
          "fill_local_block_vector_ghost_family_values: store not ready"
  end if

  if (payload_family /= BLOCK_PAYLOAD_SOL .and. &
       payload_family /= BLOCK_PAYLOAD_WAV_COEFF) then
     error stop &
          "fill_local_block_vector_ghost_family_values: invalid family"
  end if

  do b = 1, size(block_local)
     select case (payload_family)
     case (BLOCK_PAYLOAD_SOL)
        block_local(b)%ghost_vector = value
     case (BLOCK_PAYLOAD_WAV_COEFF)
        block_local(b)%ghost_wavelet_vector = value
     end select
  end do

end subroutine fill_local_block_vector_ghost_family_values


subroutine source_block_scalar_stencil_statistics ( &
     address_count,value_count,value_moment)
  ! Exercise the compact scalar stencil reader over every source block
  ! before migration. The resulting inventory is the migration-independent
  ! reference for the final-owner block store.

  implicit none

  integer(int64), intent(out) :: address_count(3)
  integer(int64), intent(out) :: value_count
  real(dp), intent(out) :: value_moment(3,3)

  integer :: i

  if (.not. allocated(block_source)) then
     error stop &
          "source_block_scalar_stencil_statistics: source unavailable"
  end if

  address_count = 0_int64
  value_count   = 0_int64
  value_moment  = 0.0_dp

  do i = 1, size(block_source)
     call accumulate_block_scalar_stencil_statistics( &
          block_source(i),address_count,value_count,value_moment)
  end do

end subroutine source_block_scalar_stencil_statistics


subroutine local_block_scalar_stencil_statistics ( &
     address_count,value_count,value_moment)
  ! Exercise the same compact scalar stencil reader over the installed
  ! final-owner block store.

  implicit none

  integer(int64), intent(out) :: address_count(3)
  integer(int64), intent(out) :: value_count
  real(dp), intent(out) :: value_moment(3,3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_scalar_stencil_statistics: store is not ready"
  end if

  address_count = 0_int64
  value_count   = 0_int64
  value_moment  = 0.0_dp

  do i = 1, size(block_local)
     call accumulate_block_scalar_stencil_statistics( &
          block_local(i),address_count,value_count,value_moment)
  end do

end subroutine local_block_scalar_stencil_statistics


subroutine accumulate_block_scalar_stencil_statistics ( &
     block,address_count,value_count,value_moment)
  ! Read all valid patch-sized windows represented by the explicit compact
  ! stencil catalogue. Field columns are ordered as sol, sol_mean and
  ! wav_coeff. This is a read-only consumer of migrated block storage.

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int64), intent(inout) :: address_count(3)
  integer(int64), intent(inout) :: value_count
  real(dp), intent(inout) :: value_moment(3,3)

  integer :: address_offset
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

  real(dp) :: value(3)

  if (block%scalar_mult /= 1) then
     error stop &
          "accumulate_block_scalar_stencil_statistics: scalar multiplier"
  end if

  if (size(block%stencil,1) /= N_BDRY .or. &
       size(block%stencil,2) /= size(block%patch)) then
     error stop &
          "accumulate_block_scalar_stencil_statistics: stencil extent"
  end if

  do p = 1, size(block%patch)
     do side = 1, N_BDRY

        storage_class = block%stencil(side,p)%storage
        record = block%stencil(side,p)%id
        address_offset = block%stencil(side,p)%offset

        select case (storage_class)

        case (STORE_PATCH)
           if (record < 0 .or. record >= size(block%patch)) then
              error stop &
                   "accumulate_block_scalar_stencil_statistics: patch ID"
           end if
           storage_start = block%patch(record+1)%elts_start
           n_storage_node = PATCH_SIZE**2

        case (STORE_BDRY)
           if (record < 1 .or. record > size(block%bdry_storage)) then
              error stop &
                   "accumulate_block_scalar_stencil_statistics: boundary ID"
           end if
           storage_start = block%bdry_storage(record)%local_start
           n_storage_node = block%bdry_storage(record)%n_node

        case (STORE_GHOST)
           if (record < 1 .or. record > size(block%ghost_storage)) then
              error stop &
                   "accumulate_block_scalar_stencil_statistics: ghost ID"
           end if
           storage_start = block%ghost_storage(record)%local_start
           n_storage_node = block%ghost_storage(record)%n_node

        case default
           error stop &
                "accumulate_block_scalar_stencil_statistics: storage class"

        end select

        do q = 0, PATCH_SIZE**2-1

           if (address_offset+q < 0 .or. &
                address_offset+q >= n_storage_node) cycle

           node_index = storage_start + address_offset + q

           select case (storage_class)
           case (STORE_PATCH)
              if (node_index < 0 .or. &
                   node_index >= size(block%node)) then
                 error stop &
                      "accumulate_block_scalar_stencil_statistics: patch node"
              end if
           case (STORE_BDRY)
              if (node_index < 0 .or. &
                   node_index >= size(block%bdry_node)) then
                 error stop &
                      "accumulate_block_scalar_stencil_statistics: boundary node"
              end if
           case (STORE_GHOST)
              if (node_index < 0 .or. &
                   node_index >= size(block%ghost_node)) then
                 error stop &
                      "accumulate_block_scalar_stencil_statistics: ghost node"
              end if
           end select

           address_count(storage_class) = &
                address_count(storage_class) + 1_int64

           do scalar_slot = 1, block%n_scalar_variable
              do level_slot = 1, block%n_field_level

                 select case (storage_class)
                 case (STORE_PATCH)
                    field_base = &
                         ((scalar_slot-1)*block%n_field_level + &
                         level_slot-1)*size(block%node)
                    field_index = field_base + node_index + 1
                    if (field_index < 1 .or. &
                         field_index > size(block%scalar)) then
                       error stop &
                            "accumulate_block_scalar_stencil_statistics: patch field"
                    end if
                    value = [ &
                         block%scalar(field_index), &
                         block%scalar_mean(field_index), &
                         block%wavelet_scalar(field_index) ]

                 case (STORE_BDRY)
                    field_base = &
                         ((scalar_slot-1)*block%n_field_level + &
                         level_slot-1)*size(block%bdry_node)
                    field_index = field_base + node_index + 1
                    if (field_index < 1 .or. &
                         field_index > size(block%bdry_scalar)) then
                       error stop &
                            "accumulate_block_scalar_stencil_statistics: boundary field"
                    end if
                    value = [ &
                         block%bdry_scalar(field_index), &
                         block%bdry_scalar_mean(field_index), &
                         block%bdry_wavelet_scalar(field_index) ]

                 case (STORE_GHOST)
                    field_base = &
                         ((scalar_slot-1)*block%n_field_level + &
                         level_slot-1)*size(block%ghost_node)
                    field_index = field_base + node_index + 1
                    if (field_index < 1 .or. &
                         field_index > size(block%ghost_scalar)) then
                       error stop &
                            "accumulate_block_scalar_stencil_statistics: ghost field"
                    end if
                    value = [ &
                         block%ghost_scalar(field_index), &
                         block%ghost_scalar_mean(field_index), &
                         block%ghost_wavelet_scalar(field_index) ]
                 end select

                 value_count = value_count + 1_int64
                 value_moment(1,:) = value_moment(1,:) + value
                 value_moment(2,:) = value_moment(2,:) + abs(value)
                 value_moment(3,:) = value_moment(3,:) + value**2

              end do
           end do

        end do

     end do
  end do

end subroutine accumulate_block_scalar_stencil_statistics


subroutine source_block_vector_stencil_statistics ( &
     address_count,value_count,value_moment)
  ! Exercise the compact vector stencil reader over every source block
  ! before migration. The resulting inventory is the migration-independent
  ! reference for the final-owner block store.

  implicit none

  integer(int64), intent(out) :: address_count(3)
  integer(int64), intent(out) :: value_count
  real(dp), intent(out) :: value_moment(3,3)

  integer :: i

  if (.not. allocated(block_source)) then
     error stop &
          "source_block_vector_stencil_statistics: source unavailable"
  end if

  address_count = 0_int64
  value_count   = 0_int64
  value_moment  = 0.0_dp

  do i = 1, size(block_source)
     call accumulate_block_vector_stencil_statistics( &
          block_source(i),address_count,value_count,value_moment)
  end do

end subroutine source_block_vector_stencil_statistics


subroutine local_block_vector_stencil_statistics ( &
     address_count,value_count,value_moment)
  ! Exercise the same compact vector stencil reader over the installed
  ! final-owner block store.

  implicit none

  integer(int64), intent(out) :: address_count(3)
  integer(int64), intent(out) :: value_count
  real(dp), intent(out) :: value_moment(3,3)

  integer :: i

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_vector_stencil_statistics: store is not ready"
  end if

  address_count = 0_int64
  value_count   = 0_int64
  value_moment  = 0.0_dp

  do i = 1, size(block_local)
     call accumulate_block_vector_stencil_statistics( &
          block_local(i),address_count,value_count,value_moment)
  end do

end subroutine local_block_vector_stencil_statistics


subroutine accumulate_block_vector_stencil_statistics ( &
     block,address_count,value_count,value_moment)
  ! Read all valid patch-sized windows represented by the explicit compact
  ! stencil catalogue. Field columns are ordered as sol, sol_mean and
  ! wav_coeff. Every stored edge component and represented level is read.

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int64), intent(inout) :: address_count(3)
  integer(int64), intent(inout) :: value_count
  real(dp), intent(inout) :: value_moment(3,3)

  integer :: address_offset
  integer :: component_slot
  integer :: field_base
  integer :: field_index
  integer :: level_slot
  integer :: n_storage_node
  integer :: node_index
  integer :: p
  integer :: q
  integer :: record
  integer :: side
  integer :: storage_class
  integer :: storage_start

  real(dp) :: value(3)

  if (block%vector_mult < 1) then
     error stop &
          "accumulate_block_vector_stencil_statistics: vector multiplier"
  end if

  if (size(block%stencil,1) /= N_BDRY .or. &
       size(block%stencil,2) /= size(block%patch)) then
     error stop &
          "accumulate_block_vector_stencil_statistics: stencil extent"
  end if

  do p = 1, size(block%patch)
     do side = 1, N_BDRY

        storage_class = block%stencil(side,p)%storage
        record = block%stencil(side,p)%id
        address_offset = block%stencil(side,p)%offset

        select case (storage_class)

        case (STORE_PATCH)
           if (record < 0 .or. record >= size(block%patch)) then
              error stop &
                   "accumulate_block_vector_stencil_statistics: patch ID"
           end if
           storage_start = block%patch(record+1)%elts_start
           n_storage_node = PATCH_SIZE**2

        case (STORE_BDRY)
           if (record < 1 .or. record > size(block%bdry_storage)) then
              error stop &
                   "accumulate_block_vector_stencil_statistics: boundary ID"
           end if
           storage_start = block%bdry_storage(record)%local_start
           n_storage_node = block%bdry_storage(record)%n_node

        case (STORE_GHOST)
           if (record < 1 .or. record > size(block%ghost_storage)) then
              error stop &
                   "accumulate_block_vector_stencil_statistics: ghost ID"
           end if
           storage_start = block%ghost_storage(record)%local_start
           n_storage_node = block%ghost_storage(record)%n_node

        case default
           error stop &
                "accumulate_block_vector_stencil_statistics: storage class"

        end select

        do q = 0, PATCH_SIZE**2-1

           if (address_offset+q < 0 .or. &
                address_offset+q >= n_storage_node) cycle

           node_index = storage_start + address_offset + q

           select case (storage_class)
           case (STORE_PATCH)
              if (node_index < 0 .or. &
                   node_index >= size(block%node)) then
                 error stop &
                      "accumulate_block_vector_stencil_statistics: patch node"
              end if
           case (STORE_BDRY)
              if (node_index < 0 .or. &
                   node_index >= size(block%bdry_node)) then
                 error stop &
                      "accumulate_block_vector_stencil_statistics: boundary node"
              end if
           case (STORE_GHOST)
              if (node_index < 0 .or. &
                   node_index >= size(block%ghost_node)) then
                 error stop &
                      "accumulate_block_vector_stencil_statistics: ghost node"
              end if
           end select

           address_count(storage_class) = &
                address_count(storage_class) + 1_int64

           do level_slot = 1, block%n_field_level
              do component_slot = 1, block%vector_mult

                 select case (storage_class)
                 case (STORE_PATCH)
                    field_base = (level_slot-1)* &
                         block%vector_mult*size(block%node)
                    field_index = field_base + &
                         block%vector_mult*node_index + component_slot
                    if (field_index < 1 .or. &
                         field_index > size(block%vector)) then
                       error stop &
                            "accumulate_block_vector_stencil_statistics: patch field"
                    end if
                    value = [ &
                         block%vector(field_index), &
                         block%vector_mean(field_index), &
                         block%wavelet_vector(field_index) ]

                 case (STORE_BDRY)
                    field_base = (level_slot-1)* &
                         block%vector_mult*size(block%bdry_node)
                    field_index = field_base + &
                         block%vector_mult*node_index + component_slot
                    if (field_index < 1 .or. &
                         field_index > size(block%bdry_vector)) then
                       error stop &
                            "accumulate_block_vector_stencil_statistics: boundary field"
                    end if
                    value = [ &
                         block%bdry_vector(field_index), &
                         block%bdry_vector_mean(field_index), &
                         block%bdry_wavelet_vector(field_index) ]

                 case (STORE_GHOST)
                    field_base = (level_slot-1)* &
                         block%vector_mult*size(block%ghost_node)
                    field_index = field_base + &
                         block%vector_mult*node_index + component_slot
                    if (field_index < 1 .or. &
                         field_index > size(block%ghost_vector)) then
                       error stop &
                            "accumulate_block_vector_stencil_statistics: ghost field"
                    end if
                    value = [ &
                         block%ghost_vector(field_index), &
                         block%ghost_vector_mean(field_index), &
                         block%ghost_wavelet_vector(field_index) ]
                 end select

                 value_count = value_count + 1_int64
                 value_moment(1,:) = value_moment(1,:) + value
                 value_moment(2,:) = value_moment(2,:) + abs(value)
                 value_moment(3,:) = value_moment(3,:) + value**2

              end do
           end do

        end do

     end do
  end do

end subroutine accumulate_block_vector_stencil_statistics


subroutine compute_local_block_hydrostatic_patch ( &
     catalog_index,local_patch,surface_pressure,dynamic_exner, &
     air_temperature)
  ! Production hydrostatic column kernel for one compact block patch.
  ! Results are returned in level-major order for direct consumption by
  ! later block-based dynamics without modifying the installed fields.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch

  real(dp), intent(out) :: surface_pressure(:)
  real(dp), intent(out) :: dynamic_exner(:)
  real(dp), intent(out) :: air_temperature(:)

  integer :: i
  integer :: k
  integer :: local_index
  integer :: mass_index
  integer :: mass_variable_slot
  integer :: n_node
  integer :: node_index
  integer :: output_index
  integer :: patch_start
  integer :: scalar_variable_size
  integer :: temperature_index
  integer :: temperature_variable_slot

  real(dp) :: layer_pressure
  real(dp) :: pressure_factor
  real(dp) :: pressure_lower
  real(dp) :: pressure_upper
  real(dp) :: rho_dz
  real(dp) :: rho_dz_theta

  if (.not. compressible) then
     error stop &
          "compute_local_block_hydrostatic_patch: incompressible case"
  end if
  if (.not. local_block_store_ready()) then
     error stop &
          "compute_local_block_hydrostatic_patch: store is not ready"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "compute_local_block_hydrostatic_patch: block is not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "compute_local_block_hydrostatic_patch: invalid patch"
  end if
  if (size(surface_pressure) /= PATCH_SIZE**2 .or. &
       size(dynamic_exner) /= zlevels*PATCH_SIZE**2 .or. &
       size(air_temperature) /= zlevels*PATCH_SIZE**2) then
     error stop &
          "compute_local_block_hydrostatic_patch: output extent"
  end if

  if (block_local(local_index)%scalar_mult /= 1) then
     error stop &
          "compute_local_block_hydrostatic_patch: scalar multiplier"
  end if
  if (block_local(local_index)%field_level > 1 .or. &
       block_local(local_index)%field_level + &
       block_local(local_index)%n_field_level - 1 < zlevels) then
     error stop &
          "compute_local_block_hydrostatic_patch: level coverage"
  end if

  mass_variable_slot = &
       S_MASS - block_local(local_index)%scalar_variable
  temperature_variable_slot = &
       S_TEMP - block_local(local_index)%scalar_variable

  if (mass_variable_slot < 0 .or. &
       mass_variable_slot >= &
       block_local(local_index)%n_scalar_variable .or. &
       temperature_variable_slot < 0 .or. &
       temperature_variable_slot >= &
       block_local(local_index)%n_scalar_variable) then
     error stop &
          "compute_local_block_hydrostatic_patch: variable coverage"
  end if

  n_node = size(block_local(local_index)%node)
  patch_start = &
       block_local(local_index)%patch(local_patch+1)%elts_start
  if (patch_start < 0 .or. &
       patch_start+PATCH_SIZE**2 > n_node) then
     error stop &
          "compute_local_block_hydrostatic_patch: patch storage"
  end if

  scalar_variable_size = &
       block_local(local_index)%n_field_level*n_node

  do i = 1, PATCH_SIZE**2
     node_index = patch_start + i
     pressure_lower = p_top

     do k = 1, zlevels
        mass_index = &
             mass_variable_slot*scalar_variable_size + &
             (k-block_local(local_index)%field_level)*n_node + &
             node_index
        rho_dz = block_local(local_index)%scalar(mass_index) + &
             block_local(local_index)%scalar_mean(mass_index)
        if (rho_dz <= 0.0_dp) then
           error stop &
                "compute_local_block_hydrostatic_patch: nonpositive mass"
        end if
        pressure_lower = pressure_lower + grav_accel*rho_dz
     end do

     surface_pressure(i) = pressure_lower

     do k = 1, zlevels
        mass_index = &
             mass_variable_slot*scalar_variable_size + &
             (k-block_local(local_index)%field_level)*n_node + &
             node_index
        temperature_index = &
             temperature_variable_slot*scalar_variable_size + &
             (k-block_local(local_index)%field_level)*n_node + &
             node_index

        rho_dz = block_local(local_index)%scalar(mass_index) + &
             block_local(local_index)%scalar_mean(mass_index)
        rho_dz_theta = &
             block_local(local_index)%scalar(temperature_index) + &
             block_local(local_index)%scalar_mean(temperature_index)

        pressure_upper = pressure_lower - grav_accel*rho_dz
        layer_pressure = 0.5_dp*(pressure_lower+pressure_upper)
        if (layer_pressure <= 0.0_dp) then
           error stop &
                "compute_local_block_hydrostatic_patch: nonpositive pressure"
        end if

        pressure_factor = (layer_pressure/p_0)**kappa
        output_index = (k-1)*PATCH_SIZE**2 + i
        dynamic_exner(output_index) = c_p*pressure_factor
        air_temperature(output_index) = &
             rho_dz_theta/rho_dz*pressure_factor
        pressure_lower = pressure_upper
     end do
  end do

end subroutine compute_local_block_hydrostatic_patch


subroutine clear_local_block_hydrostatic_state
  ! Release derived hydrostatic fields. They are regenerated from the
  ! installed compact scalar and scalar-mean storage after migration.

  implicit none

  block_hydrostatic_ready = .false.
  block_hydrostatic_refreshes = 0_int64

  if (allocated(block_hydrostatic)) then
     deallocate(block_hydrostatic)
  end if

  if (allocated(block_hydrostatic)) then
     error stop &
          "clear_local_block_hydrostatic_state: cleanup failed"
  end if

end subroutine clear_local_block_hydrostatic_state


subroutine invalidate_local_block_hydrostatic_block (local_index)
  ! Mark only one block's derived thermodynamic fields stale.

  implicit none

  integer, intent(in) :: local_index

  block_hydrostatic_ready = .false.

  if (.not. allocated(block_hydrostatic)) return
  if (size(block_hydrostatic) /= size(block_local)) return
  if (local_index < 1 .or. local_index > size(block_hydrostatic)) return

  block_hydrostatic(local_index)%ready = .false.

end subroutine invalidate_local_block_hydrostatic_block


subroutine invalidate_local_block_hydrostatic_state
  ! Mark all local blocks' derived thermodynamic fields stale.

  implicit none

  block_hydrostatic_ready = .false.
  if (allocated(block_hydrostatic)) then
     block_hydrostatic%ready = .false.
  end if

end subroutine invalidate_local_block_hydrostatic_state


subroutine prepare_local_block_hydrostatic_state
  ! Ensure the per-block cache catalogue exists without computing fields.

  implicit none

  integer :: local_index

  if (allocated(block_hydrostatic)) then
     if (size(block_hydrostatic) /= size(block_local)) then
        deallocate(block_hydrostatic)
     end if
  end if
  if (.not. allocated(block_hydrostatic)) then
     allocate(block_hydrostatic(size(block_local)))
     do local_index = 1,size(block_local)
        block_hydrostatic(local_index)%catalog_index = &
             block_local_catalog_index(local_index)
        block_hydrostatic(local_index)%ready = .false.
        block_hydrostatic(local_index)%refreshes = 0_int64
     end do
  end if

end subroutine prepare_local_block_hydrostatic_state


subroutine refresh_local_block_hydrostatic_block (local_index)
  ! Refresh one reusable, patch-major thermodynamic block cache.

  implicit none

  integer, intent(in) :: local_index

  integer :: catalog_index
  integer :: column_base
  integer :: local_patch
  integer :: n_column
  integer :: n_patch
  integer :: n_surface
  integer :: surface_base

  if (local_index < 1 .or. local_index > size(block_local)) then
     error stop &
          "refresh_local_block_hydrostatic_block: invalid local index"
  end if

  call prepare_local_block_hydrostatic_state

  catalog_index = block_local_catalog_index(local_index)
  n_patch = size(block_local(local_index)%patch)
  n_surface = n_patch*PATCH_SIZE**2
  n_column = n_patch*zlevels*PATCH_SIZE**2

  block_hydrostatic(local_index)%catalog_index = catalog_index
  block_hydrostatic(local_index)%ready = .false.

  if (allocated( &
       block_hydrostatic(local_index)%surface_pressure)) then
     if (size(block_hydrostatic(local_index)%surface_pressure) /= &
          n_surface) then
        deallocate(block_hydrostatic(local_index)%surface_pressure)
     end if
  end if
  if (.not. allocated( &
       block_hydrostatic(local_index)%surface_pressure)) then
     allocate( &
          block_hydrostatic(local_index)%surface_pressure(n_surface))
  end if

  if (allocated(block_hydrostatic(local_index)%dynamic_exner)) then
     if (size(block_hydrostatic(local_index)%dynamic_exner) /= &
          n_column) then
        deallocate(block_hydrostatic(local_index)%dynamic_exner)
     end if
  end if
  if (.not. allocated( &
       block_hydrostatic(local_index)%dynamic_exner)) then
     allocate(block_hydrostatic(local_index)%dynamic_exner(n_column))
  end if

  if (allocated(block_hydrostatic(local_index)%air_temperature)) then
     if (size(block_hydrostatic(local_index)%air_temperature) /= &
          n_column) then
        deallocate(block_hydrostatic(local_index)%air_temperature)
     end if
  end if
  if (.not. allocated( &
       block_hydrostatic(local_index)%air_temperature)) then
     allocate(block_hydrostatic(local_index)%air_temperature(n_column))
  end if

  do local_patch = 0,n_patch-1
     surface_base = local_patch*PATCH_SIZE**2
     column_base = local_patch*zlevels*PATCH_SIZE**2

     call compute_local_block_hydrostatic_patch( &
          catalog_index,local_patch, &
          block_hydrostatic(local_index)%surface_pressure( &
          surface_base+1:surface_base+PATCH_SIZE**2), &
          block_hydrostatic(local_index)%dynamic_exner( &
          column_base+1:column_base+zlevels*PATCH_SIZE**2), &
          block_hydrostatic(local_index)%air_temperature( &
          column_base+1:column_base+zlevels*PATCH_SIZE**2))
  end do

  block_hydrostatic(local_index)%ready = .true.
  block_hydrostatic(local_index)%refreshes = &
       block_hydrostatic(local_index)%refreshes + 1_int64
  block_hydrostatic_refreshes = &
       block_hydrostatic_refreshes + 1_int64

end subroutine refresh_local_block_hydrostatic_block


subroutine ensure_local_block_hydrostatic_block (catalog_index)
  ! Lazily refresh only the requested local block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  if (.not. compressible) then
     error stop &
          "ensure_local_block_hydrostatic_block: incompressible case"
  end if
  if (.not. local_block_store_ready()) then
     error stop &
          "ensure_local_block_hydrostatic_block: store is not ready"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "ensure_local_block_hydrostatic_block: block not local"
  end if

  call prepare_local_block_hydrostatic_state

  if (.not. block_hydrostatic(local_index)%ready) then
     call refresh_local_block_hydrostatic_block(local_index)
  end if

  block_hydrostatic_ready = all(block_hydrostatic%ready)

end subroutine ensure_local_block_hydrostatic_block


subroutine refresh_local_block_hydrostatic_state
  ! Explicitly recompute every local block while retaining allocations.

  implicit none

  integer :: local_index

  if (.not. compressible) then
     error stop &
          "refresh_local_block_hydrostatic_state: incompressible case"
  end if
  if (.not. local_block_store_ready()) then
     error stop &
          "refresh_local_block_hydrostatic_state: store is not ready"
  end if

  call prepare_local_block_hydrostatic_state
  block_hydrostatic_ready = .false.
  block_hydrostatic%ready = .false.

  do local_index = 1,size(block_local)
     call refresh_local_block_hydrostatic_block(local_index)
  end do

  block_hydrostatic_ready = all(block_hydrostatic%ready)

end subroutine refresh_local_block_hydrostatic_state


subroutine ensure_local_block_hydrostatic_state
  ! Refresh only stale local block caches for whole-store consumers.

  implicit none

  integer :: local_index

  if (.not. compressible) then
     error stop &
          "ensure_local_block_hydrostatic_state: incompressible case"
  end if
  if (.not. local_block_store_ready()) then
     error stop &
          "ensure_local_block_hydrostatic_state: store is not ready"
  end if

  call prepare_local_block_hydrostatic_state

  do local_index = 1,size(block_local)
     if (.not. block_hydrostatic(local_index)%ready) then
        call refresh_local_block_hydrostatic_block(local_index)
     end if
  end do

  block_hydrostatic_ready = all(block_hydrostatic%ready)

end subroutine ensure_local_block_hydrostatic_state


logical function local_block_hydrostatic_state_ready () result(ready)
  ! Report whether complete derived hydrostatic storage corresponds to
  ! the currently installed local block catalogue.

  implicit none

  integer :: local_index
  integer :: n_patch

  ready = .false.

  if (.not. local_block_store_ready()) return
  if (.not. block_hydrostatic_ready) return
  if (.not. allocated(block_hydrostatic)) return
  if (size(block_hydrostatic) /= size(block_local)) return

  do local_index = 1,size(block_local)
     if (block_hydrostatic(local_index)%catalog_index /= &
          block_local_catalog_index(local_index)) return
     if (.not. block_hydrostatic(local_index)%ready) return

     n_patch = size(block_local(local_index)%patch)
     if (.not. allocated( &
          block_hydrostatic(local_index)%surface_pressure)) return
     if (.not. allocated( &
          block_hydrostatic(local_index)%dynamic_exner)) return
     if (.not. allocated( &
          block_hydrostatic(local_index)%air_temperature)) return
     if (size(block_hydrostatic(local_index)%surface_pressure) /= &
          n_patch*PATCH_SIZE**2) return
     if (size(block_hydrostatic(local_index)%dynamic_exner) /= &
          n_patch*zlevels*PATCH_SIZE**2) return
     if (size(block_hydrostatic(local_index)%air_temperature) /= &
          n_patch*zlevels*PATCH_SIZE**2) return
  end do

  ready = .true.

end function local_block_hydrostatic_state_ready


integer(int64) function local_block_hydrostatic_refresh_count () &
     result(n_refresh)
  ! Number of completed per-block refreshes since the current local
  ! block store was installed. This supports cache-lifecycle diagnostics.

  implicit none

  n_refresh = block_hydrostatic_refreshes

end function local_block_hydrostatic_refresh_count


integer(int64) function local_block_hydrostatic_block_refresh_count ( &
     catalog_index) result(n_refresh)
  ! Number of completed refreshes for one local block cache.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  if (.not. local_block_store_ready()) then
     error stop &
          "local_block_hydrostatic_block_refresh_count: store not ready"
  end if
  if (.not. allocated(block_hydrostatic)) then
     error stop &
          "local_block_hydrostatic_block_refresh_count: cache missing"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_hydrostatic_block_refresh_count: block not local"
  end if
  if (block_hydrostatic(local_index)%catalog_index /= catalog_index) then
     error stop &
          "local_block_hydrostatic_block_refresh_count: catalogue mismatch"
  end if

  n_refresh = block_hydrostatic(local_index)%refreshes

end function local_block_hydrostatic_block_refresh_count


integer function local_block_hydrostatic_surface_nvalue ( &
     catalog_index) result(n_value)
  ! Number of persistent surface-pressure values in one local block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  call ensure_local_block_hydrostatic_block(catalog_index)

  if (.not. local_block_store_ready() .or. &
       .not. allocated(block_hydrostatic)) then
     error stop &
          "local_block_hydrostatic_surface_nvalue: state not ready"
  end if
  if (size(block_hydrostatic) /= size(block_local)) then
     error stop &
          "local_block_hydrostatic_surface_nvalue: store extent"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_hydrostatic_surface_nvalue: block not local"
  end if
  if (block_hydrostatic(local_index)%catalog_index /= &
       catalog_index .or. &
       .not. block_hydrostatic(local_index)%ready .or. &
       .not. allocated( &
       block_hydrostatic(local_index)%surface_pressure)) then
     error stop &
          "local_block_hydrostatic_surface_nvalue: local storage"
  end if

  n_value = size(block_hydrostatic(local_index)%surface_pressure)

end function local_block_hydrostatic_surface_nvalue


integer function local_block_hydrostatic_column_nvalue ( &
     catalog_index) result(n_value)
  ! Number of persistent Exner or temperature values in one local block.

  implicit none

  integer, intent(in) :: catalog_index

  integer :: local_index

  call ensure_local_block_hydrostatic_block(catalog_index)

  if (.not. local_block_store_ready() .or. &
       .not. allocated(block_hydrostatic)) then
     error stop &
          "local_block_hydrostatic_column_nvalue: state not ready"
  end if
  if (size(block_hydrostatic) /= size(block_local)) then
     error stop &
          "local_block_hydrostatic_column_nvalue: store extent"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "local_block_hydrostatic_column_nvalue: block not local"
  end if
  if (block_hydrostatic(local_index)%catalog_index /= &
       catalog_index .or. &
       .not. block_hydrostatic(local_index)%ready .or. &
       .not. allocated( &
       block_hydrostatic(local_index)%dynamic_exner)) then
     error stop &
          "local_block_hydrostatic_column_nvalue: local storage"
  end if

  n_value = size(block_hydrostatic(local_index)%dynamic_exner)

end function local_block_hydrostatic_column_nvalue


subroutine get_local_block_hydrostatic_patch_values ( &
     catalog_index,local_patch,surface_pressure,dynamic_exner, &
     air_temperature)
  ! Return one patch from persistent, patch-major hydrostatic storage.

  implicit none

  integer, intent(in) :: catalog_index
  integer, intent(in) :: local_patch

  real(dp), intent(out) :: surface_pressure(:)
  real(dp), intent(out) :: dynamic_exner(:)
  real(dp), intent(out) :: air_temperature(:)

  integer :: column_base
  integer :: local_index
  integer :: surface_base

  call ensure_local_block_hydrostatic_block(catalog_index)

  if (.not. local_block_store_ready() .or. &
       .not. allocated(block_hydrostatic)) then
     error stop &
          "get_local_block_hydrostatic_patch_values: state not ready"
  end if
  if (size(block_hydrostatic) /= size(block_local)) then
     error stop &
          "get_local_block_hydrostatic_patch_values: store extent"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_hydrostatic_patch_values: block not local"
  end if
  if (local_patch < 0 .or. &
       local_patch >= size(block_local(local_index)%patch)) then
     error stop &
          "get_local_block_hydrostatic_patch_values: invalid patch"
  end if
  if (size(surface_pressure) /= PATCH_SIZE**2 .or. &
       size(dynamic_exner) /= zlevels*PATCH_SIZE**2 .or. &
       size(air_temperature) /= zlevels*PATCH_SIZE**2) then
     error stop &
          "get_local_block_hydrostatic_patch_values: output extent"
  end if
  if (block_hydrostatic(local_index)%catalog_index /= &
       catalog_index .or. &
       .not. block_hydrostatic(local_index)%ready) then
     error stop &
          "get_local_block_hydrostatic_patch_values: catalogue mismatch"
  end if

  surface_base = local_patch*PATCH_SIZE**2
  column_base = local_patch*zlevels*PATCH_SIZE**2

  surface_pressure = &
       block_hydrostatic(local_index)%surface_pressure( &
       surface_base+1:surface_base+PATCH_SIZE**2)
  dynamic_exner = block_hydrostatic(local_index)%dynamic_exner( &
       column_base+1:column_base+zlevels*PATCH_SIZE**2)
  air_temperature = &
       block_hydrostatic(local_index)%air_temperature( &
       column_base+1:column_base+zlevels*PATCH_SIZE**2)

end subroutine get_local_block_hydrostatic_patch_values


subroutine get_local_block_hydrostatic_values ( &
     catalog_index,surface_pressure,dynamic_exner,air_temperature)
  ! Return all persistent hydrostatic values for one local block. Arrays
  ! retain patch-major order and level-major order within each patch.

  implicit none

  integer, intent(in) :: catalog_index

  real(dp), intent(out) :: surface_pressure(:)
  real(dp), intent(out) :: dynamic_exner(:)
  real(dp), intent(out) :: air_temperature(:)

  integer :: local_index

  call ensure_local_block_hydrostatic_block(catalog_index)

  if (.not. local_block_store_ready() .or. &
       .not. allocated(block_hydrostatic)) then
     error stop &
          "get_local_block_hydrostatic_values: state not ready"
  end if
  if (size(block_hydrostatic) /= size(block_local)) then
     error stop &
          "get_local_block_hydrostatic_values: store extent"
  end if

  local_index = catalog_local_block(catalog_index)
  if (local_index < 1) then
     error stop &
          "get_local_block_hydrostatic_values: block not local"
  end if
  if (block_hydrostatic(local_index)%catalog_index /= &
       catalog_index .or. &
       .not. block_hydrostatic(local_index)%ready .or. &
       .not. allocated( &
       block_hydrostatic(local_index)%surface_pressure) .or. &
       .not. allocated( &
       block_hydrostatic(local_index)%dynamic_exner) .or. &
       .not. allocated( &
       block_hydrostatic(local_index)%air_temperature)) then
     error stop &
          "get_local_block_hydrostatic_values: local storage"
  end if
  if (size(surface_pressure) /= &
       size(block_hydrostatic(local_index)%surface_pressure) .or. &
       size(dynamic_exner) /= &
       size(block_hydrostatic(local_index)%dynamic_exner) .or. &
       size(air_temperature) /= &
       size(block_hydrostatic(local_index)%air_temperature)) then
     error stop &
          "get_local_block_hydrostatic_values: output extent"
  end if

  surface_pressure = &
       block_hydrostatic(local_index)%surface_pressure
  dynamic_exner = block_hydrostatic(local_index)%dynamic_exner
  air_temperature = &
       block_hydrostatic(local_index)%air_temperature

end subroutine get_local_block_hydrostatic_values


subroutine apply_local_block_hydrostatic_consumer (consumer,context)
  ! Apply one caller-supplied read-only production kernel to every local
  ! thermodynamic block without constructing patch or whole-block copies.

  implicit none

  procedure(Local_Block_Hydrostatic_Consumer) :: consumer
  class(*), intent(inout) :: context

  integer :: local_index
  integer(int64) :: refresh_count_before

  call ensure_local_block_hydrostatic_state

  if (.not. local_block_hydrostatic_state_ready()) then
     error stop &
          "apply_local_block_hydrostatic_consumer: state not ready"
  end if

  refresh_count_before = local_block_hydrostatic_refresh_count()

  do local_index = 1,size(block_local)
     call consumer( &
          block_local_catalog_index(local_index), &
          size(block_local(local_index)%patch), &
          block_hydrostatic(local_index)%surface_pressure, &
          block_hydrostatic(local_index)%dynamic_exner, &
          block_hydrostatic(local_index)%air_temperature,context)

     if (.not. block_hydrostatic(local_index)%ready) then
        error stop &
             "apply_local_block_hydrostatic_consumer: cache invalidated"
     end if
  end do

  if (.not. local_block_hydrostatic_state_ready()) then
     error stop &
          "apply_local_block_hydrostatic_consumer: state changed"
  end if
  if (local_block_hydrostatic_refresh_count() /= &
       refresh_count_before) then
     error stop &
          "apply_local_block_hydrostatic_consumer: unexpected refresh"
  end if

end subroutine apply_local_block_hydrostatic_consumer


subroutine local_block_hydrostatic_statistics ( &
     surface_count,column_count,surface_moment,exner_moment, &
     temperature_moment)
  ! Accumulate diagnostics from persistent production hydrostatic state.
  ! The state is regenerated after migration rather than serialized.

  implicit none

  integer(int64), intent(out) :: surface_count
  integer(int64), intent(out) :: column_count

  real(dp), intent(out) :: surface_moment(3)
  real(dp), intent(out) :: exner_moment(3)
  real(dp), intent(out) :: temperature_moment(3)

  integer :: catalog_index
  integer :: column_nvalue
  integer :: local_index
  integer :: surface_nvalue

  real(dp), allocatable :: surface_pressure(:)
  real(dp), allocatable :: dynamic_exner(:)
  real(dp), allocatable :: air_temperature(:)

  if (.not. compressible) then
     error stop &
          "local_block_hydrostatic_statistics: incompressible case"
  end if

  call ensure_local_block_hydrostatic_state

  if (.not. local_block_hydrostatic_state_ready()) then
     error stop &
          "local_block_hydrostatic_statistics: state is not ready"
  end if

  surface_count      = 0_int64
  column_count       = 0_int64
  surface_moment     = 0.0_dp
  exner_moment       = 0.0_dp
  temperature_moment = 0.0_dp

  do local_index = 1, size(block_local)
     catalog_index = block_local_catalog_index(local_index)
     surface_nvalue = &
          local_block_hydrostatic_surface_nvalue(catalog_index)
     column_nvalue = &
          local_block_hydrostatic_column_nvalue(catalog_index)

     if (allocated(surface_pressure)) then
        if (size(surface_pressure) /= surface_nvalue) then
           deallocate(surface_pressure)
        end if
     end if
     if (.not. allocated(surface_pressure)) then
        allocate(surface_pressure(surface_nvalue))
     end if

     if (allocated(dynamic_exner)) then
        if (size(dynamic_exner) /= column_nvalue) then
           deallocate(dynamic_exner)
        end if
     end if
     if (.not. allocated(dynamic_exner)) then
        allocate(dynamic_exner(column_nvalue))
     end if

     if (allocated(air_temperature)) then
        if (size(air_temperature) /= column_nvalue) then
           deallocate(air_temperature)
        end if
     end if
     if (.not. allocated(air_temperature)) then
        allocate(air_temperature(column_nvalue))
     end if

     call get_local_block_hydrostatic_values( &
          catalog_index,surface_pressure,dynamic_exner,air_temperature)

     surface_count = surface_count + &
          int(size(surface_pressure),int64)
     surface_moment(1) = surface_moment(1) + sum(surface_pressure)
     surface_moment(2) = surface_moment(2) + &
          sum(abs(surface_pressure))
     surface_moment(3) = surface_moment(3) + &
          sum(surface_pressure**2)

     column_count = column_count + int(size(dynamic_exner),int64)
     exner_moment(1) = exner_moment(1) + sum(dynamic_exner)
     exner_moment(2) = exner_moment(2) + sum(abs(dynamic_exner))
     exner_moment(3) = exner_moment(3) + sum(dynamic_exner**2)

     temperature_moment(1) = temperature_moment(1) + &
          sum(air_temperature)
     temperature_moment(2) = temperature_moment(2) + &
          sum(abs(air_temperature))
     temperature_moment(3) = temperature_moment(3) + &
          sum(air_temperature**2)
  end do

end subroutine local_block_hydrostatic_statistics


subroutine clear_local_block_tendency_state
  ! Release derived tendency outputs before the local block store changes.

  implicit none

  block_tendency_ready = .false.
  block_tendency_import_active = .false.
  block_tendency_trial_active = .false.
  block_tendency_commit_checkpoint_ready = .false.
  block_tendency_accumulator_ready = .false.
  block_tendency_executions = 0_int64
  block_tendency_allocations = 0_int64
  block_tendency_import_allocations = 0_int64
  block_tendency_accumulator_allocations = 0_int64
  block_tendency_accumulator_stages = 0_int64

  if (allocated(block_tendency)) deallocate(block_tendency)
  if (allocated(block_tendency_import_patch_count)) then
     deallocate(block_tendency_import_patch_count)
  end if
  if (allocated(block_tendency_trial)) deallocate(block_tendency_trial)
  if (allocated(block_tendency_accumulator)) then
     deallocate(block_tendency_accumulator)
  end if

  if (allocated(block_tendency) .or. &
       allocated(block_tendency_import_patch_count) .or. &
       allocated(block_tendency_trial) .or. &
       allocated(block_tendency_accumulator)) then
     error stop "clear_local_block_tendency_state: cleanup failed"
  end if

end subroutine clear_local_block_tendency_state


subroutine clear_local_blocks
  ! Invalidate and release the persistent final-owner local store.
  ! This routine is deliberately idempotent so it is safe before the
  ! first restart and before every subsequent checkpoint restart.

  implicit none

  call clear_local_block_tendency_state
  call clear_local_block_hydrostatic_state
  block_store_ready = .false.

  if (allocated(block_local)) then
     deallocate(block_local)
  end if

  if (allocated(block_local_catalog_index)) then
     deallocate(block_local_catalog_index)
  end if

  if (allocated(block_catalog_local_index)) then
     deallocate(block_catalog_local_index)
  end if

  if (allocated(block_local) .or. &
       allocated(block_local_catalog_index) .or. &
       allocated(block_catalog_local_index)) then
     error stop "clear_local_blocks: cleanup failed"
  end if

end subroutine clear_local_blocks


subroutine clear_block_staging
  ! Release only temporary source/receive migration stores. The
  ! final-owner block_local store and its catalogue are retained.

  implicit none

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
  if (allocated(block_received)) then
     deallocate(block_received)
  end if
  if (allocated(block_received_catalog_index)) then
     deallocate(block_received_catalog_index)
  end if

  if (allocated(block_source) .or. &
       allocated(block_source_catalog_index) .or. &
       allocated(block_retained_source_index) .or. &
       allocated(block_migrating_source_index) .or. &
       allocated(block_received) .or. &
       allocated(block_received_catalog_index)) then

     error stop "clear_block_staging: cleanup failed"

  end if

end subroutine clear_block_staging

integer function packed_block_nbyte (block) result(nbyte)
  ! Return the exact byte count used by pack_block.

  implicit none

  type(Block_Data), intent(in) :: block

  integer :: nbyte_integer

  if (.not. allocated(block%patch) .or. &
       .not. allocated(block%node) .or. &
       .not. allocated(block%scalar) .or. &
       .not. allocated(block%vector) .or. &
       .not. allocated(block%wavelet_scalar) .or. &
       .not. allocated(block%wavelet_vector) .or. &
       .not. allocated(block%scalar_mean) .or. &
       .not. allocated(block%vector_mean) .or. &
       .not. allocated(block%tke) .or. &
       .not. allocated(block%wavelet_tke) .or. &
       .not. allocated(block%topography) .or. &
       .not. allocated(block%neigh_class) .or. &
       .not. allocated(block%block_bdry) .or. &
       .not. allocated(block%bdry_storage) .or. &
       .not. allocated(block%stencil) .or. &
       .not. allocated(block%bdry_node) .or. &
       .not. allocated(block%bdry_scalar) .or. &
       .not. allocated(block%bdry_vector) .or. &
       .not. allocated(block%bdry_wavelet_scalar) .or. &
       .not. allocated(block%bdry_wavelet_vector) .or. &
       .not. allocated(block%bdry_scalar_mean) .or. &
       .not. allocated(block%bdry_vector_mean) .or. &
       .not. allocated(block%bdry_tke) .or. &
       .not. allocated(block%bdry_wavelet_tke) .or. &
       .not. allocated(block%bdry_topography) .or. &
       .not. allocated(block%ghost_storage) .or. &
       .not. allocated(block%ghost_node) .or. &
       .not. allocated(block%ghost_scalar) .or. &
       .not. allocated(block%ghost_vector) .or. &
       .not. allocated(block%ghost_wavelet_scalar) .or. &
       .not. allocated(block%ghost_wavelet_vector) .or. &
       .not. allocated(block%ghost_scalar_mean) .or. &
       .not. allocated(block%ghost_vector_mean) .or. &
       .not. allocated(block%ghost_tke) .or. &
       .not. allocated(block%ghost_wavelet_tke) .or. &
       .not. allocated(block%ghost_topography)) then

     error stop "packed_block_nbyte: unallocated component"

  end if

  nbyte_integer = storage_size(0) / 8

  nbyte = BLOCK_PACK_HEADER_SIZE * nbyte_integer

  nbyte = nbyte + &
       size(block%patch) * storage_size(block%patch) / 8

  nbyte = nbyte + &
       size(block%node) * storage_size(block%node) / 8

  nbyte = nbyte + &
       size(block%scalar) * storage_size(block%scalar) / 8

  nbyte = nbyte + &
       size(block%vector) * storage_size(block%vector) / 8

  nbyte = nbyte + &
       size(block%wavelet_scalar) * &
       storage_size(block%wavelet_scalar) / 8

  nbyte = nbyte + &
       size(block%wavelet_vector) * &
       storage_size(block%wavelet_vector) / 8

  nbyte = nbyte + &
       size(block%scalar_mean) * storage_size(block%scalar_mean) / 8

  nbyte = nbyte + &
       size(block%vector_mean) * storage_size(block%vector_mean) / 8

  nbyte = nbyte + &
       size(block%tke) * storage_size(block%tke) / 8

  nbyte = nbyte + &
       size(block%wavelet_tke) * storage_size(block%wavelet_tke) / 8

  nbyte = nbyte + &
       size(block%topography) * storage_size(block%topography) / 8

  nbyte = nbyte + size(block%neigh_class) * nbyte_integer

  nbyte = nbyte + &
       size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8

  nbyte = nbyte + &
       size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8

  nbyte = nbyte + &
       size(block%stencil) * storage_size(block%stencil) / 8

  nbyte = nbyte + &
       size(block%bdry_node) * storage_size(block%bdry_node) / 8

  nbyte = nbyte + &
       size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8

  nbyte = nbyte + &
       size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8

  nbyte = nbyte + &
       size(block%bdry_wavelet_scalar) * &
       storage_size(block%bdry_wavelet_scalar) / 8

  nbyte = nbyte + &
       size(block%bdry_wavelet_vector) * &
       storage_size(block%bdry_wavelet_vector) / 8

  nbyte = nbyte + &
       size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8

  nbyte = nbyte + &
       size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8

  nbyte = nbyte + &
       size(block%bdry_tke) * storage_size(block%bdry_tke) / 8

  nbyte = nbyte + &
       size(block%bdry_wavelet_tke) * &
       storage_size(block%bdry_wavelet_tke) / 8

  nbyte = nbyte + &
       size(block%bdry_topography) * &
       storage_size(block%bdry_topography) / 8

  nbyte = nbyte + &
       size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8

  nbyte = nbyte + &
       size(block%ghost_node) * storage_size(block%ghost_node) / 8

  nbyte = nbyte + &
       size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8

  nbyte = nbyte + &
       size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8

  nbyte = nbyte + &
       size(block%ghost_wavelet_scalar) * &
       storage_size(block%ghost_wavelet_scalar) / 8

  nbyte = nbyte + &
       size(block%ghost_wavelet_vector) * &
       storage_size(block%ghost_wavelet_vector) / 8

  nbyte = nbyte + &
       size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8

  nbyte = nbyte + &
       size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8

  nbyte = nbyte + &
       size(block%ghost_tke) * storage_size(block%ghost_tke) / 8

  nbyte = nbyte + &
       size(block%ghost_wavelet_tke) * &
       storage_size(block%ghost_wavelet_tke) / 8

  nbyte = nbyte + &
       size(block%ghost_topography) * &
       storage_size(block%ghost_topography) / 8

end function packed_block_nbyte


subroutine pack_block (block,buffer)
  ! Serialize a Block_Data into a versioned contiguous byte buffer.
  ! This raw representation is intended for homogeneous MPI jobs.

  implicit none

  type(Block_Data), intent(in) :: block
  integer(int8), allocatable, intent(out) :: buffer(:)

  integer :: header(BLOCK_PACK_HEADER_SIZE)
  integer :: n
  integer :: nbyte
  integer :: pos

  nbyte = packed_block_nbyte(block)

  allocate(buffer(nbyte))

  header = [ &
       BLOCK_PACK_MAGIC, &
       BLOCK_PACK_VERSION, &
       block%id, &
       block%root_domain, &
       block%root_patch, &
       block%level, &
       block%scalar_variable, &
       block%n_scalar_variable, &
       block%vector_variable, &
       block%field_level, &
       block%n_field_level, &
       block%scalar_mult, &
       block%vector_mult, &
       size(block%patch), &
       size(block%node), &
       size(block%scalar), &
       size(block%vector), &
       size(block%scalar_mean), &
       size(block%vector_mean), &
       size(block%neigh_class,1), &
       size(block%neigh_class,2), &
       size(block%block_bdry), &
       size(block%bdry_storage), &
       size(block%stencil,1), &
       size(block%stencil,2), &
       size(block%bdry_node), &
       size(block%bdry_scalar), &
       size(block%bdry_vector), &
       size(block%bdry_scalar_mean), &
       size(block%bdry_vector_mean), &
       size(block%ghost_storage), &
       size(block%ghost_node), &
       size(block%ghost_scalar), &
       size(block%ghost_vector), &
       size(block%ghost_scalar_mean), &
       size(block%ghost_vector_mean), &
       block%tke_level, &
       block%n_tke_level, &
       size(block%tke), &
       size(block%bdry_tke), &
       size(block%ghost_tke), &
       size(block%topography), &
       size(block%bdry_topography), &
       size(block%ghost_topography) ]

  pos = 0

  n = size(header) * storage_size(header) / 8
  buffer(pos+1:pos+n) = transfer(header,0_int8,n)
  pos = pos + n

  n = size(block%patch) * storage_size(block%patch) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%patch,0_int8,n)
     pos = pos + n
  end if

  n = size(block%node) * storage_size(block%node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%scalar) * storage_size(block%scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%vector) * storage_size(block%vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%wavelet_scalar) * &
       storage_size(block%wavelet_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%wavelet_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%wavelet_vector) * &
       storage_size(block%wavelet_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%wavelet_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%scalar_mean) * &
       storage_size(block%scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%vector_mean) * &
       storage_size(block%vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%vector_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%tke) * storage_size(block%tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%wavelet_tke) * storage_size(block%wavelet_tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%wavelet_tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%topography) * storage_size(block%topography) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%topography,0_int8,n)
     pos = pos + n
  end if

  n = size(block%neigh_class) * storage_size(0) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%neigh_class,0_int8,n)
     pos = pos + n
  end if

  n = size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%block_bdry,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_storage,0_int8,n)
     pos = pos + n
  end if

  n = size(block%stencil) * storage_size(block%stencil) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%stencil,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_node) * storage_size(block%bdry_node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_scalar) * &
       storage_size(block%bdry_wavelet_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_wavelet_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_vector) * &
       storage_size(block%bdry_wavelet_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_wavelet_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%bdry_vector_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_tke) * storage_size(block%bdry_tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_tke) * &
       storage_size(block%bdry_wavelet_tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_wavelet_tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%bdry_topography) * &
       storage_size(block%bdry_topography) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%bdry_topography,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_storage,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_node) * storage_size(block%ghost_node) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_node,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_scalar) * &
       storage_size(block%ghost_wavelet_scalar) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_wavelet_scalar,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_vector) * &
       storage_size(block%ghost_wavelet_vector) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_wavelet_vector,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_scalar_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = &
          transfer(block%ghost_vector_mean,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_tke) * storage_size(block%ghost_tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_tke) * &
       storage_size(block%ghost_wavelet_tke) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_wavelet_tke,0_int8,n)
     pos = pos + n
  end if

  n = size(block%ghost_topography) * &
       storage_size(block%ghost_topography) / 8
  if (n > 0) then
     buffer(pos+1:pos+n) = transfer(block%ghost_topography,0_int8,n)
     pos = pos + n
  end if

  if (pos /= size(buffer)) then
     error stop "pack_block: final byte count mismatch"
  end if

end subroutine pack_block


subroutine unpack_block (buffer,block)
  ! Reconstruct a Block_Data from pack_block's byte stream.

  implicit none

  integer(int8), intent(in) :: buffer(:)
  type(Block_Data), intent(out) :: block

  integer :: header(BLOCK_PACK_HEADER_SIZE)
  integer :: n
  integer :: pos
  integer :: scalar_variable
  integer :: n_scalar_variable
  integer :: vector_variable
  integer :: field_level
  integer :: n_field_level
  integer :: scalar_mult
  integer :: vector_mult
  integer :: tke_level
  integer :: n_tke_level

  n = size(header) * storage_size(header) / 8

  if (size(buffer) < n) then
     error stop "unpack_block: truncated header"
  end if

  header = transfer(buffer(1:n),header,size(header))
  pos = n

  if (header(1) /= BLOCK_PACK_MAGIC) then
     error stop "unpack_block: invalid pack magic"
  end if

  if (header(2) /= BLOCK_PACK_VERSION) then
     error stop "unpack_block: unsupported pack version"
  end if

  if (any(header(14:BLOCK_PACK_HEADER_SIZE) < 0)) then
     error stop "unpack_block: negative component extent"
  end if

  call get_block_field_layout( &
       scalar_variable,n_scalar_variable,vector_variable,field_level, &
       n_field_level,scalar_mult,vector_mult)

  call get_block_turbulence_layout(tke_level,n_tke_level)

  if (header(7) /= scalar_variable .or. &
       header(8) /= n_scalar_variable .or. &
       header(9) /= vector_variable .or. &
       header(10) /= field_level .or. &
       header(11) /= n_field_level .or. &
       header(12) /= scalar_mult .or. &
       header(13) /= vector_mult) then

     error stop "unpack_block: unsupported field layout"

  end if

  if (header(37) /= tke_level .or. &
       header(38) /= n_tke_level) then
     error stop "unpack_block: unsupported turbulence layout"
  end if

  if (header(20) /= N_BDRY .or. &
       header(21) /= header(14) .or. &
       header(24) /= N_BDRY .or. &
       header(25) /= header(14)) then

     error stop "unpack_block: invalid topology extents"

  end if

  if (header(15) /= header(14)*PATCH_SIZE**2) then
     error stop "unpack_block: invalid interior node extent"
  end if

  if (header(39) /= header(38)*header(15) .or. &
       header(40) /= header(38)*header(26) .or. &
       header(41) /= header(38)*header(32)) then
     error stop "unpack_block: invalid turbulence extents"
  end if

  if (header(42) /= header(15) .or. &
       header(43) /= header(26) .or. &
       header(44) /= header(32)) then
     error stop "unpack_block: invalid topography extents"
  end if

  block%id          = header(3)
  block%root_domain = header(4)
  block%root_patch  = header(5)
  block%level       = header(6)
  block%scalar_variable = header(7)
  block%n_scalar_variable = header(8)
  block%vector_variable = header(9)
  block%field_level     = header(10)
  block%n_field_level   = header(11)
  block%scalar_mult     = header(12)
  block%vector_mult     = header(13)
  block%tke_level       = header(37)
  block%n_tke_level     = header(38)

  allocate(block%patch(header(14)))
  allocate(block%node(header(15)))
  allocate(block%scalar(header(16)))
  allocate(block%vector(header(17)))
  allocate(block%wavelet_scalar(header(16)))
  allocate(block%wavelet_vector(header(17)))
  allocate(block%scalar_mean(header(18)))
  allocate(block%vector_mean(header(19)))
  allocate(block%tke(header(39)))
  allocate(block%wavelet_tke(header(39)))
  allocate(block%topography(header(42)))
  allocate(block%neigh_class(header(20),header(21)))
  allocate(block%block_bdry(header(22)))
  allocate(block%bdry_storage(header(23)))
  allocate(block%stencil(header(24),header(25)))
  allocate(block%bdry_node(header(26)))
  allocate(block%bdry_scalar(header(27)))
  allocate(block%bdry_vector(header(28)))
  allocate(block%bdry_wavelet_scalar(header(27)))
  allocate(block%bdry_wavelet_vector(header(28)))
  allocate(block%bdry_scalar_mean(header(29)))
  allocate(block%bdry_vector_mean(header(30)))
  allocate(block%bdry_tke(header(40)))
  allocate(block%bdry_wavelet_tke(header(40)))
  allocate(block%bdry_topography(header(43)))
  allocate(block%ghost_storage(header(31)))
  allocate(block%ghost_node(header(32)))
  allocate(block%ghost_scalar(header(33)))
  allocate(block%ghost_vector(header(34)))
  allocate(block%ghost_wavelet_scalar(header(33)))
  allocate(block%ghost_wavelet_vector(header(34)))
  allocate(block%ghost_scalar_mean(header(35)))
  allocate(block%ghost_vector_mean(header(36)))
  allocate(block%ghost_tke(header(41)))
  allocate(block%ghost_wavelet_tke(header(41)))
  allocate(block%ghost_topography(header(44)))

  n = size(block%patch) * storage_size(block%patch) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated patch data"
  end if
  if (n > 0) then
     block%patch = transfer( &
          buffer(pos+1:pos+n),block%patch,size(block%patch))
     pos = pos + n
  end if

  n = size(block%node) * storage_size(block%node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated node data"
  end if
  if (n > 0) then
     block%node = transfer( &
          buffer(pos+1:pos+n),block%node,size(block%node))
     pos = pos + n
  end if

  n = size(block%scalar) * storage_size(block%scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated scalar data"
  end if
  if (n > 0) then
     block%scalar = transfer( &
          buffer(pos+1:pos+n),block%scalar,size(block%scalar))
     pos = pos + n
  end if

  n = size(block%vector) * storage_size(block%vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated vector data"
  end if
  if (n > 0) then
     block%vector = transfer( &
          buffer(pos+1:pos+n),block%vector,size(block%vector))
     pos = pos + n
  end if

  n = size(block%wavelet_scalar) * &
       storage_size(block%wavelet_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated wavelet scalar data"
  end if
  if (n > 0) then
     block%wavelet_scalar = transfer( &
          buffer(pos+1:pos+n),block%wavelet_scalar, &
          size(block%wavelet_scalar))
     pos = pos + n
  end if

  n = size(block%wavelet_vector) * &
       storage_size(block%wavelet_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated wavelet vector data"
  end if
  if (n > 0) then
     block%wavelet_vector = transfer( &
          buffer(pos+1:pos+n),block%wavelet_vector, &
          size(block%wavelet_vector))
     pos = pos + n
  end if

  n = size(block%scalar_mean) * &
       storage_size(block%scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean scalar data"
  end if
  if (n > 0) then
     block%scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%scalar_mean, &
          size(block%scalar_mean))
     pos = pos + n
  end if

  n = size(block%vector_mean) * &
       storage_size(block%vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean vector data"
  end if
  if (n > 0) then
     block%vector_mean = transfer( &
          buffer(pos+1:pos+n),block%vector_mean, &
          size(block%vector_mean))
     pos = pos + n
  end if

  n = size(block%tke) * storage_size(block%tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated tke data"
  end if
  if (n > 0) then
     block%tke = transfer( &
          buffer(pos+1:pos+n),block%tke,size(block%tke))
     pos = pos + n
  end if

  n = size(block%wavelet_tke) * storage_size(block%wavelet_tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated wavelet tke data"
  end if
  if (n > 0) then
     block%wavelet_tke = transfer( &
          buffer(pos+1:pos+n),block%wavelet_tke, &
          size(block%wavelet_tke))
     pos = pos + n
  end if

  n = size(block%topography) * storage_size(block%topography) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated topography data"
  end if
  if (n > 0) then
     block%topography = transfer( &
          buffer(pos+1:pos+n),block%topography, &
          size(block%topography))
     pos = pos + n
  end if

  n = size(block%neigh_class) * storage_size(0) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated neighbour data"
  end if
  if (n > 0) then
     block%neigh_class = reshape( &
          transfer(buffer(pos+1:pos+n),0, &
          size(block%neigh_class)),shape(block%neigh_class))
     pos = pos + n
  end if

  n = size(block%block_bdry) * &
       storage_size(block%block_bdry) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated block boundary data"
  end if
  if (n > 0) then
     block%block_bdry = transfer( &
          buffer(pos+1:pos+n),block%block_bdry, &
          size(block%block_bdry))
     pos = pos + n
  end if

  n = size(block%bdry_storage) * &
       storage_size(block%bdry_storage) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary catalogue"
  end if
  if (n > 0) then
     block%bdry_storage = transfer( &
          buffer(pos+1:pos+n),block%bdry_storage, &
          size(block%bdry_storage))
     pos = pos + n
  end if

  n = size(block%stencil) * storage_size(block%stencil) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated stencil data"
  end if
  if (n > 0) then
     block%stencil = reshape( &
          transfer(buffer(pos+1:pos+n),block%stencil(1,1), &
          size(block%stencil)),shape(block%stencil))
     pos = pos + n
  end if

  n = size(block%bdry_node) * storage_size(block%bdry_node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary nodes"
  end if
  if (n > 0) then
     block%bdry_node = transfer( &
          buffer(pos+1:pos+n),block%bdry_node, &
          size(block%bdry_node))
     pos = pos + n
  end if

  n = size(block%bdry_scalar) * &
       storage_size(block%bdry_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary scalar"
  end if
  if (n > 0) then
     block%bdry_scalar = transfer( &
          buffer(pos+1:pos+n),block%bdry_scalar, &
          size(block%bdry_scalar))
     pos = pos + n
  end if

  n = size(block%bdry_vector) * &
       storage_size(block%bdry_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary vector"
  end if
  if (n > 0) then
     block%bdry_vector = transfer( &
          buffer(pos+1:pos+n),block%bdry_vector, &
          size(block%bdry_vector))
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_scalar) * &
       storage_size(block%bdry_wavelet_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary wavelet scalar"
  end if
  if (n > 0) then
     block%bdry_wavelet_scalar = transfer( &
          buffer(pos+1:pos+n),block%bdry_wavelet_scalar, &
          size(block%bdry_wavelet_scalar))
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_vector) * &
       storage_size(block%bdry_wavelet_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary wavelet vector"
  end if
  if (n > 0) then
     block%bdry_wavelet_vector = transfer( &
          buffer(pos+1:pos+n),block%bdry_wavelet_vector, &
          size(block%bdry_wavelet_vector))
     pos = pos + n
  end if

  n = size(block%bdry_scalar_mean) * &
       storage_size(block%bdry_scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean boundary scalar"
  end if
  if (n > 0) then
     block%bdry_scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%bdry_scalar_mean, &
          size(block%bdry_scalar_mean))
     pos = pos + n
  end if

  n = size(block%bdry_vector_mean) * &
       storage_size(block%bdry_vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean boundary vector"
  end if
  if (n > 0) then
     block%bdry_vector_mean = transfer( &
          buffer(pos+1:pos+n),block%bdry_vector_mean, &
          size(block%bdry_vector_mean))
     pos = pos + n
  end if

  n = size(block%bdry_tke) * storage_size(block%bdry_tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary tke"
  end if
  if (n > 0) then
     block%bdry_tke = transfer( &
          buffer(pos+1:pos+n),block%bdry_tke, &
          size(block%bdry_tke))
     pos = pos + n
  end if

  n = size(block%bdry_wavelet_tke) * &
       storage_size(block%bdry_wavelet_tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary wavelet tke"
  end if
  if (n > 0) then
     block%bdry_wavelet_tke = transfer( &
          buffer(pos+1:pos+n),block%bdry_wavelet_tke, &
          size(block%bdry_wavelet_tke))
     pos = pos + n
  end if

  n = size(block%bdry_topography) * &
       storage_size(block%bdry_topography) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated boundary topography"
  end if
  if (n > 0) then
     block%bdry_topography = transfer( &
          buffer(pos+1:pos+n),block%bdry_topography, &
          size(block%bdry_topography))
     pos = pos + n
  end if

  n = size(block%ghost_storage) * &
       storage_size(block%ghost_storage) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost catalogue"
  end if
  if (n > 0) then
     block%ghost_storage = transfer( &
          buffer(pos+1:pos+n),block%ghost_storage, &
          size(block%ghost_storage))
     pos = pos + n
  end if

  n = size(block%ghost_node) * storage_size(block%ghost_node) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost nodes"
  end if
  if (n > 0) then
     block%ghost_node = transfer( &
          buffer(pos+1:pos+n),block%ghost_node, &
          size(block%ghost_node))
     pos = pos + n
  end if

  n = size(block%ghost_scalar) * &
       storage_size(block%ghost_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost scalar"
  end if
  if (n > 0) then
     block%ghost_scalar = transfer( &
          buffer(pos+1:pos+n),block%ghost_scalar, &
          size(block%ghost_scalar))
     pos = pos + n
  end if

  n = size(block%ghost_vector) * &
       storage_size(block%ghost_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost vector"
  end if
  if (n > 0) then
     block%ghost_vector = transfer( &
          buffer(pos+1:pos+n),block%ghost_vector, &
          size(block%ghost_vector))
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_scalar) * &
       storage_size(block%ghost_wavelet_scalar) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost wavelet scalar"
  end if
  if (n > 0) then
     block%ghost_wavelet_scalar = transfer( &
          buffer(pos+1:pos+n),block%ghost_wavelet_scalar, &
          size(block%ghost_wavelet_scalar))
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_vector) * &
       storage_size(block%ghost_wavelet_vector) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost wavelet vector"
  end if
  if (n > 0) then
     block%ghost_wavelet_vector = transfer( &
          buffer(pos+1:pos+n),block%ghost_wavelet_vector, &
          size(block%ghost_wavelet_vector))
     pos = pos + n
  end if

  n = size(block%ghost_scalar_mean) * &
       storage_size(block%ghost_scalar_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean ghost scalar"
  end if
  if (n > 0) then
     block%ghost_scalar_mean = transfer( &
          buffer(pos+1:pos+n),block%ghost_scalar_mean, &
          size(block%ghost_scalar_mean))
     pos = pos + n
  end if

  n = size(block%ghost_vector_mean) * &
       storage_size(block%ghost_vector_mean) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated mean ghost vector"
  end if
  if (n > 0) then
     block%ghost_vector_mean = transfer( &
          buffer(pos+1:pos+n),block%ghost_vector_mean, &
          size(block%ghost_vector_mean))
     pos = pos + n
  end if

  n = size(block%ghost_tke) * storage_size(block%ghost_tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost tke"
  end if
  if (n > 0) then
     block%ghost_tke = transfer( &
          buffer(pos+1:pos+n),block%ghost_tke, &
          size(block%ghost_tke))
     pos = pos + n
  end if

  n = size(block%ghost_wavelet_tke) * &
       storage_size(block%ghost_wavelet_tke) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost wavelet tke"
  end if
  if (n > 0) then
     block%ghost_wavelet_tke = transfer( &
          buffer(pos+1:pos+n),block%ghost_wavelet_tke, &
          size(block%ghost_wavelet_tke))
     pos = pos + n
  end if

  n = size(block%ghost_topography) * &
       storage_size(block%ghost_topography) / 8
  if (pos+n > size(buffer)) then
     error stop "unpack_block: truncated ghost topography"
  end if
  if (n > 0) then
     block%ghost_topography = transfer( &
          buffer(pos+1:pos+n),block%ghost_topography, &
          size(block%ghost_topography))
     pos = pos + n
  end if

  if (pos /= size(buffer)) then
     error stop "unpack_block: final byte count mismatch"
  end if

end subroutine unpack_block

end module parallel_block_mod
