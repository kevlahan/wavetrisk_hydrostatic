module parallel_block_build_mod

  use, intrinsic :: iso_fortran_env, only : int8

  use kind_mod, only : dp

  use shared_mod, only : BDRY_THICKNESS, Coord, EDGE, N_BDRY, &
       N_CHDRN

  use ops_mod,   only : comp_offs3
  use patch_mod, only : Patch, PATCH_SIZE

  use parallel_block_mod, only : Block_Data, &
       STORE_NONE, STORE_PATCH, STORE_BDRY, STORE_GHOST, &
       NGB_INTERNAL, NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT, NGB_OTHER, &
       block_source, block_source_catalog_index, &
       block_retained_source_index, block_migrating_source_index, &
       pack_block, unpack_block, check_block_storage, &
       get_block_field_layout, get_block_turbulence_layout

  use arch_mod, only : block_catalog, loc_id, n_process, owner, rank

  use domain_mod, only : grid, sol, sol_mean, tke, topography, &
       wav_coeff, wav_tke, &
       count_subtree_patches_Domain, &
       extract_subtree_patches_Domain, subtree_depth_Domain, &
       compact_subtree_storage_Domain, copy_subtree_nodes_Domain, &
       copy_subtree_field_Domain, renumber_subtree_neigh_Domain, &
       get_bdry_dims_Domain

  implicit none

  private
  public :: build_source_blocks

contains


subroutine build_source_blocks (verbose)
  ! Build and verify every candidate block whose source Domain is
  ! currently local to this rank. Optional diagnostics describe one
  ! representative block and the rank totals; validation always runs.

  implicit none

  logical, optional, intent(in) :: verbose

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

  logical :: print_summary

  print_summary = .false.
  if (present(verbose)) print_summary = verbose

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
     error stop "build_source_blocks: no source-local blocks"
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
  if (print_summary .and. rank == 0) then

  do b = 1, size(block_catalog)

     if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

     d = loc_id(block_catalog(b)%root_domain+1) + 1

     if (d < 1 .or. d > size(grid)) then
        error stop &
             "build_source_blocks: invalid local domain mapping"
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
  if (print_summary .and. rank == 0 .and. b_verbose < 1) then

     do b = 1, size(block_catalog)

        if (owner(block_catalog(b)%root_domain+1) /= rank) cycle

        d = loc_id(block_catalog(b)%root_domain+1) + 1

        if (d < 1 .or. d > size(grid)) then
           error stop &
                "build_source_blocks: invalid local domain mapping"
        end if

        b_verbose = b
        exit


     end do

  end if

  if (print_summary .and. rank == 0 .and. b_verbose < 1) then
     error stop &
          "build_source_blocks: no local source-domain block"
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
             "build_source_blocks: invalid local domain mapping"
     end if

     i_block = i_block + 1

     block_source_catalog_index(i_block) = b

     call build_one_source_block( &
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
     error stop "build_source_blocks: no blocks constructed"
  end if

  if (n_block_owned + n_block_migrating /= n_block_built) then
     error stop &
          "build_source_blocks: block migration count mismatch"
  end if

  if (i_block /= size(block_source)) then
     error stop &
          "build_source_blocks: persistent block count mismatch"
  end if

  call check_source_blocks


  call check_migrating_block_serialization( &
       n_pack_block, n_pack_byte_total, n_pack_byte_max)

  if (print_summary) then

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

  end if

end subroutine build_source_blocks


subroutine check_source_blocks
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
          "check_source_blocks: block store is not allocated"
  end if

  if (.not. allocated(block_source_catalog_index)) then
     error stop &
          "check_source_blocks: catalog map is not allocated"
  end if

  if (.not. allocated(block_retained_source_index) .or. &
       .not. allocated(block_migrating_source_index)) then
     error stop &
          "check_source_blocks: ownership sets not allocated"
  end if

  if (size(block_source_catalog_index) /= size(block_source)) then
     error stop &
          "check_source_blocks: catalog map size mismatch"
  end if

  allocate(seen(size(block_source)))
  seen = 0

  do i = 1, size(block_source)

     b = block_source_catalog_index(i)

     if (b < 1 .or. b > size(block_catalog)) then
        error stop &
             "check_source_blocks: invalid catalog index"
     end if

     if (owner(block_catalog(b)%root_domain+1) /= rank) then
        error stop &
             "check_source_blocks: source owner mismatch"
     end if

     if (block_source(i)%id /= block_catalog(b)%id .or. &
          block_source(i)%root_domain /= &
          block_catalog(b)%root_domain .or. &
          block_source(i)%root_patch /= &
          block_catalog(b)%root_patch .or. &
          block_source(i)%level /= block_catalog(b)%level) then

        error stop &
             "check_source_blocks: block identity mismatch"

     end if

     call check_block_storage(block_source(i))

  end do

  do i = 1, size(block_retained_source_index)

     ib = block_retained_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "check_source_blocks: invalid retained index"
     end if

     if (seen(ib) /= 0) then
        error stop &
             "check_source_blocks: duplicate retained index"
     end if

     b = block_source_catalog_index(ib)

     if (block_catalog(b)%owner /= rank) then
        error stop &
             "check_source_blocks: retained owner mismatch"
     end if

     seen(ib) = 1

  end do

  do i = 1, size(block_migrating_source_index)

     ib = block_migrating_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "check_source_blocks: invalid migrating index"
     end if

     if (seen(ib) /= 0) then
        error stop &
             "check_source_blocks: duplicate migrating index"
     end if

     b = block_source_catalog_index(ib)

     if (block_catalog(b)%owner == rank) then
        error stop &
             "check_source_blocks: migrating owner mismatch"
     end if

     seen(ib) = 1

  end do

  if (any(seen /= 1)) then
     error stop &
          "check_source_blocks: ownership partition mismatch"
  end if

  deallocate(seen)

end subroutine check_source_blocks


subroutine check_migrating_block_serialization ( &
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
          "check_migrating_block_serialization: block store not allocated"

  end if

  n_block_out    = 0
  n_byte_out     = 0
  n_byte_max_out = 0

  do i = 1, size(block_migrating_source_index)

     ib = block_migrating_source_index(i)

     if (ib < 1 .or. ib > size(block_source)) then
        error stop &
             "check_migrating_block_serialization: invalid block index"
     end if

     call pack_block(block_source(ib),buffer_source)
     call unpack_block(buffer_source,block_copy)
     call pack_block(block_copy,buffer_copy)

     if (size(buffer_copy) /= size(buffer_source)) then
        error stop &
             "check_migrating_block_serialization: packed size mismatch"
     end if

     if (any(buffer_copy /= buffer_source)) then
        error stop &
             "check_migrating_block_serialization: round-trip mismatch"
     end if

     nbyte = size(buffer_source)

     n_block_out = n_block_out + 1
     n_byte_out = n_byte_out + nbyte
     n_byte_max_out = max(n_byte_max_out,nbyte)

  end do

  if (n_block_out /= size(block_migrating_source_index)) then
     error stop &
          "check_migrating_block_serialization: tested count mismatch"
  end if

end subroutine check_migrating_block_serialization


subroutine build_one_source_block ( &
     b_catalog, block_out, verbose, &
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

  integer, intent(in) :: b_catalog
  type(Block_Data), intent(out) :: block_out
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
  integer :: next_source_patch
  integer :: source_local_patch
  integer :: source_patch_count

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
  integer :: first_field_level
  integer :: mult_scalar
  integer :: mult_vector
  integer :: n_field_level
  integer :: n_scalar_variable
  integer :: tke_level
  integer :: n_tke_level
  integer :: tke_storage_size
  integer :: field_base
  integer :: field_level
  integer :: level_slot
  integer :: scalar_base
  integer :: scalar_id
  integer :: scalar_slot
  integer :: scalar_storage_size
  integer :: scalar_variable_size
  integer :: vector_storage_size

  integer :: offs_src(0:N_BDRY)
  integer :: dims_src(2,N_BDRY)

  integer, allocatable :: old_to_new(:)
  integer, allocatable :: old_elts_start(:)

  integer, allocatable :: bdry_required(:)
  integer, allocatable :: bdry_closure(:)
  integer, allocatable :: bdry_level(:)
  integer, allocatable :: ghost_patch(:)

  real(dp) :: val_src
  real(dp) :: val_blk

  real(dp), allocatable :: scalar_copy(:)
  real(dp), allocatable :: wavelet_scalar_copy(:)
  real(dp), allocatable :: wavelet_scalar_one(:)
  real(dp), allocatable :: scalar_mean_copy(:)
  real(dp), allocatable :: scalar_mean_one(:)
  real(dp), allocatable :: scalar_one(:)
  real(dp), allocatable :: vector_copy(:)
  real(dp), allocatable :: wavelet_vector_copy(:)
  real(dp), allocatable :: wavelet_vector_one(:)
  real(dp), allocatable :: vector_mean_copy(:)
  real(dp), allocatable :: vector_mean_one(:)
  real(dp), allocatable :: vector_one(:)
  real(dp), allocatable :: tke_copy(:)
  real(dp), allocatable :: tke_one(:)
  real(dp), allocatable :: wavelet_tke_copy(:)
  real(dp), allocatable :: wavelet_tke_one(:)
  real(dp), allocatable :: topography_copy(:)

  type(Patch), allocatable :: patch_copy(:)
  type(Coord), allocatable :: node_copy(:)

  logical :: already_present

  if (owner(block_catalog(b_catalog)%root_domain+1) /= rank) then
     error stop &
          "build_source_blocks: source domain is not local"
  end if

  d = loc_id(block_catalog(b_catalog)%root_domain+1) + 1

  if (d < 1 .or. d > size(grid)) then
     error stop "build_source_blocks: invalid local domain"
  end if

  p_root   = block_catalog(b_catalog)%root_patch
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
     error stop "build_source_blocks: patch count mismatch"
  end if

  if (old_to_new(p_root) /= 0) then
     error stop "build_source_blocks: extracted root is not patch zero"
  end if

  allocate(old_elts_start(size(patch_copy)))

  do i = 1, size(patch_copy)
     old_elts_start(i) = patch_copy(i)%elts_start
  end do

  do i = 1, size(patch_copy)

     if (old_elts_start(i) < 0) then
        error stop "build_source_blocks: invalid original elts_start"
     end if

     if (old_elts_start(i) + PATCH_SIZE**2 > grid(d)%node%length) then
        error stop &
             "build_source_blocks: original node storage out of bounds"
     end if

  end do

  call compact_subtree_storage_Domain(patch_copy)

  do i = 1, size(patch_copy)

     if (patch_copy(i)%elts_start /= (i-1)*PATCH_SIZE**2) then
        error stop &
             "build_source_blocks: incorrect compact elts_start"
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
          "build_source_blocks: incorrect copied node storage size"
  end if

  do i = 1, size(patch_copy)

     old_start = old_elts_start(i)
     new_start = patch_copy(i)%elts_start

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%x - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%x)) > 0.0_dp) then

        error stop &
             "build_source_blocks: node x-coordinate mismatch"

     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%y - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%y)) > 0.0_dp) then

        error stop &
             "build_source_blocks: node y-coordinate mismatch"

     end if

     if (maxval(abs( &
          node_copy(new_start+1:new_start+PATCH_SIZE**2)%z - &
          grid(d)%node%elts( &
          old_start+1:old_start+PATCH_SIZE**2)%z)) > 0.0_dp) then

        error stop &
             "build_source_blocks: node z-coordinate mismatch"

     end if

  end do

  !
  ! ===============================================================
  ! Copy and verify every scalar field over the prognostic column.
  ! ===============================================================
  !
  call get_block_field_layout( &
       v_scalar,n_scalar_variable,v_vector,first_field_level, &
       n_field_level,mult_scalar,mult_vector)

  call get_block_turbulence_layout(tke_level,n_tke_level)

  if (mult_scalar /= 1) then
     error stop &
          "build_source_blocks: unexpected scalar multiplier"
  end if

  if (n_scalar_variable < 1) then
     error stop "build_source_blocks: no scalar variables"
  end if

  if (n_field_level < 1) then
     error stop "build_source_blocks: no prognostic field levels"
  end if

  scalar_storage_size = mult_scalar * n_node_storage
  scalar_variable_size = n_field_level * scalar_storage_size
  allocate(scalar_copy(n_scalar_variable*scalar_variable_size))
  allocate(wavelet_scalar_copy( &
       n_scalar_variable*scalar_variable_size))
  allocate(scalar_mean_copy( &
       n_scalar_variable*scalar_variable_size))

  n_patch_field = mult_scalar * PATCH_SIZE**2

  do scalar_slot = 1, n_scalar_variable

     scalar_id = v_scalar + scalar_slot - 1

     do level_slot = 1, n_field_level

        field_level = first_field_level + level_slot - 1
        scalar_base = (scalar_slot-1) * scalar_variable_size + &
             (level_slot-1) * scalar_storage_size

        call copy_subtree_field_Domain( &
             patch_copy, old_elts_start, mult_scalar, &
             sol(scalar_id,field_level)%data(d)%elts, scalar_one)

        if (size(scalar_one) /= scalar_storage_size) then
           error stop &
                "build_source_blocks: incorrect scalar storage size"
        end if

        scalar_copy( &
             scalar_base+1:scalar_base+scalar_storage_size) = &
             scalar_one

        call copy_subtree_field_Domain( &
             patch_copy, old_elts_start, mult_scalar, &
             wav_coeff(scalar_id,field_level)%data(d)%elts, &
             wavelet_scalar_one)

        if (size(wavelet_scalar_one) /= scalar_storage_size) then
           error stop &
                "build_source_blocks: incorrect scalar wavelet storage"
        end if

        wavelet_scalar_copy( &
             scalar_base+1:scalar_base+scalar_storage_size) = &
             wavelet_scalar_one

        call copy_subtree_field_Domain( &
             patch_copy, old_elts_start, mult_scalar, &
             sol_mean(scalar_id,field_level)%data(d)%elts, &
             scalar_mean_one)

        if (size(scalar_mean_one) /= scalar_storage_size) then
           error stop &
                "build_source_blocks: incorrect scalar mean storage size"
        end if

        scalar_mean_copy( &
             scalar_base+1:scalar_base+scalar_storage_size) = &
             scalar_mean_one

        do i = 1, size(patch_copy)

           old_start = mult_scalar * old_elts_start(i)
           new_start = scalar_base + &
                mult_scalar * patch_copy(i)%elts_start

           if (maxval(abs( &
                scalar_copy(new_start+1:new_start+n_patch_field) - &
                sol(scalar_id,field_level)%data(d)%elts( &
                old_start+1:old_start+n_patch_field))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: scalar field copy mismatch"

           end if

           if (maxval(abs( &
                wavelet_scalar_copy( &
                new_start+1:new_start+n_patch_field) - &
                wav_coeff(scalar_id,field_level)%data(d)%elts( &
                old_start+1:old_start+n_patch_field))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: scalar wavelet copy mismatch"

           end if

           if (maxval(abs( &
                scalar_mean_copy( &
                new_start+1:new_start+n_patch_field) - &
                sol_mean(scalar_id,field_level)%data(d)%elts( &
                old_start+1:old_start+n_patch_field))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: scalar mean copy mismatch"

           end if

        end do

     end do

  end do

  !
  ! ===============================================================
  ! Copy and verify velocity over the complete sol level range.
  ! ===============================================================
  !
  if (mult_vector /= EDGE) then
     error stop &
          "build_source_blocks: unexpected vector multiplier"
  end if

  vector_storage_size = mult_vector * n_node_storage
  allocate(vector_copy(n_field_level*vector_storage_size))
  allocate(wavelet_vector_copy(n_field_level*vector_storage_size))
  allocate(vector_mean_copy(n_field_level*vector_storage_size))
  n_patch_field = mult_vector * PATCH_SIZE**2

  do level_slot = 1, n_field_level

     field_level = first_field_level + level_slot - 1
     field_base = (level_slot-1) * vector_storage_size

     call copy_subtree_field_Domain( &
          patch_copy, old_elts_start, mult_vector, &
          sol(v_vector,field_level)%data(d)%elts, vector_one)

     if (size(vector_one) /= vector_storage_size) then
        error stop &
             "build_source_blocks: incorrect vector storage size"
     end if

     vector_copy(field_base+1:field_base+vector_storage_size) = &
          vector_one

     call copy_subtree_field_Domain( &
          patch_copy, old_elts_start, mult_vector, &
          wav_coeff(v_vector,field_level)%data(d)%elts, &
          wavelet_vector_one)

     if (size(wavelet_vector_one) /= vector_storage_size) then
        error stop &
             "build_source_blocks: incorrect vector wavelet storage"
     end if

     wavelet_vector_copy( &
          field_base+1:field_base+vector_storage_size) = &
          wavelet_vector_one

     call copy_subtree_field_Domain( &
          patch_copy, old_elts_start, mult_vector, &
          sol_mean(v_vector,field_level)%data(d)%elts, &
          vector_mean_one)

     if (size(vector_mean_one) /= vector_storage_size) then
        error stop &
             "build_source_blocks: incorrect vector mean storage size"
     end if

     vector_mean_copy(field_base+1:field_base+vector_storage_size) = &
          vector_mean_one

     do i = 1, size(patch_copy)

        old_start = mult_vector * old_elts_start(i)
        new_start = field_base + &
             mult_vector * patch_copy(i)%elts_start

        if (maxval(abs( &
             vector_copy(new_start+1:new_start+n_patch_field) - &
             sol(v_vector,field_level)%data(d)%elts( &
             old_start+1:old_start+n_patch_field))) > 0.0_dp) then

           error stop &
                "build_source_blocks: vector field copy mismatch"

        end if

        if (maxval(abs( &
             wavelet_vector_copy( &
             new_start+1:new_start+n_patch_field) - &
             wav_coeff(v_vector,field_level)%data(d)%elts( &
             old_start+1:old_start+n_patch_field))) > 0.0_dp) then

           error stop &
                "build_source_blocks: vector wavelet copy mismatch"

        end if

        if (maxval(abs( &
             vector_mean_copy( &
             new_start+1:new_start+n_patch_field) - &
             sol_mean(v_vector,field_level)%data(d)%elts( &
             old_start+1:old_start+n_patch_field))) > 0.0_dp) then

           error stop &
                "build_source_blocks: vector mean copy mismatch"

        end if

     end do

  end do

  !
  ! ===============================================================
  ! Copy and verify optional turbulent kinetic energy fields.
  ! ===============================================================
  !
  tke_storage_size = n_node_storage
  allocate(tke_copy(n_tke_level*tke_storage_size))
  allocate(wavelet_tke_copy(n_tke_level*tke_storage_size))
  n_patch_field = PATCH_SIZE**2

  do level_slot = 1, n_tke_level

     field_level = tke_level + level_slot - 1
     field_base = (level_slot-1) * tke_storage_size

     call copy_subtree_field_Domain( &
          patch_copy,old_elts_start,1, &
          tke(field_level)%data(d)%elts,tke_one)

     call copy_subtree_field_Domain( &
          patch_copy,old_elts_start,1, &
          wav_tke(field_level)%data(d)%elts,wavelet_tke_one)

     if (size(tke_one) /= tke_storage_size .or. &
          size(wavelet_tke_one) /= tke_storage_size) then
        error stop "build_source_blocks: incorrect tke storage size"
     end if

     tke_copy(field_base+1:field_base+tke_storage_size) = tke_one
     wavelet_tke_copy( &
          field_base+1:field_base+tke_storage_size) = wavelet_tke_one

     do i = 1, size(patch_copy)

        old_start = old_elts_start(i)
        new_start = field_base + patch_copy(i)%elts_start

        if (maxval(abs( &
             tke_copy(new_start+1:new_start+n_patch_field) - &
             tke(field_level)%data(d)%elts( &
             old_start+1:old_start+n_patch_field))) > 0.0_dp) then
           error stop "build_source_blocks: tke copy mismatch"
        end if

        if (maxval(abs( &
             wavelet_tke_copy(new_start+1:new_start+n_patch_field) - &
             wav_tke(field_level)%data(d)%elts( &
             old_start+1:old_start+n_patch_field))) > 0.0_dp) then
           error stop "build_source_blocks: wavelet tke copy mismatch"
        end if

     end do

  end do

  !
  ! ===============================================================
  ! Copy and verify topography over the complete extracted subtree.
  ! ===============================================================
  !
  call copy_subtree_field_Domain( &
       patch_copy,old_elts_start,1, &
       topography%data(d)%elts,topography_copy)

  if (size(topography_copy) /= n_node_storage) then
     error stop "build_source_blocks: incorrect topography storage size"
  end if

  n_patch_field = PATCH_SIZE**2

  do i = 1, size(patch_copy)

     old_start = old_elts_start(i)
     new_start = patch_copy(i)%elts_start

     if (maxval(abs( &
          topography_copy(new_start+1:new_start+n_patch_field) - &
          topography%data(d)%elts( &
          old_start+1:old_start+n_patch_field))) > 0.0_dp) then
        error stop "build_source_blocks: topography copy mismatch"
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
                   "build_source_blocks: unexpected child link"
           end if

        else if (grid(d)%patch%elts(p_chd_old+1)%deleted) then

           if (p_chd_new /= 0) then
              error stop &
                   "build_source_blocks: deleted child copied"
           end if

        else

           if (old_to_new(p_chd_old) < 0) then
              error stop &
                   "build_source_blocks: child missing from map"
           end if

           if (p_chd_new /= old_to_new(p_chd_old)) then
              error stop &
                   "build_source_blocks: incorrect child renumbering"
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
     error stop "build_source_blocks: leaf count mismatch"
  end if

  depth_new = copied_depth(0)

  if (depth_new /= depth_old) then
     error stop "build_source_blocks: subtree depth mismatch"
  end if

  !
  ! ===============================================================
  ! Classify source neighbour links.
  ! ===============================================================
  !
  allocate(block_out%neigh_class(N_BDRY,size(patch_copy)))

  block_out%neigh_class = NGB_OTHER

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

        if (p_ngb_old > 0) then

           if (p_ngb_old >= grid(d)%patch%length) then
              error stop &
                   "build_source_blocks: invalid positive neighbour"
           end if

           if (old_to_new(p_ngb_old) >= 0) then

              block_out%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_INTERNAL

           else

              block_out%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_BLOCK

           end if

        else if (p_ngb_old < 0) then

           b_src = -p_ngb_old

           if (b_src >= grid(d)%bdry_patch%length) then
              error stop &
                   "build_source_blocks: invalid source boundary"
           end if

           if (grid(d)%bdry_patch%elts(b_src+1)%side > 0) then

              block_out%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_DOMAIN

           else

              block_out%neigh_class( &
                   c,old_to_new(p_old)+1) = NGB_ADAPT

           end if

        else

           block_out%neigh_class( &
                c,old_to_new(p_old)+1) = NGB_OTHER

        end if

     end do

  end do

  n_ngb_internal = count(block_out%neigh_class == NGB_INTERNAL)
  n_ngb_block    = count(block_out%neigh_class == NGB_BLOCK)
  n_ngb_domain   = count(block_out%neigh_class == NGB_DOMAIN)
  n_ngb_adapt    = count(block_out%neigh_class == NGB_ADAPT)
  n_ngb_other    = count(block_out%neigh_class == NGB_OTHER)

  if (n_ngb_internal + n_ngb_block + n_ngb_domain + &
       n_ngb_adapt + n_ngb_other /= &
       N_BDRY*size(patch_copy)) then

     error stop &
          "build_source_blocks: neighbour count mismatch"

  end if

  !
  ! ===============================================================
  ! Build explicit local boundary-link catalogue.
  ! ===============================================================
  !
  n_bdry_local = n_ngb_block + n_ngb_domain + n_ngb_adapt

  allocate(block_out%block_bdry(n_bdry_local))

  ib = 0

  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        select case (block_out%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL, NGB_OTHER)

           cycle

        case (NGB_BLOCK)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (p_ngb_old <= 0 .or. &
                p_ngb_old >= grid(d)%patch%length) then

              error stop &
                   "build_source_blocks: invalid inter-block neighbour"

           end if

           ib = ib + 1

           block_out%block_bdry(ib)%patch = old_to_new(p_old)
           block_out%block_bdry(ib)%side  = c
           block_out%block_bdry(ib)%class = NGB_BLOCK

           block_out%block_bdry(ib)%root_domain = &
                block_catalog(b_catalog)%root_domain

           block_out%block_bdry(ib)%neigh_patch = p_ngb_old

           block_out%block_bdry(ib)%source_block = -1
           block_out%block_bdry(ib)%source_block_id = -1
           block_out%block_bdry(ib)%source_owner = -1
           block_out%block_bdry(ib)%ghost_id     = -1
           block_out%block_bdry(ib)%storage_id = -1

        case (NGB_DOMAIN, NGB_ADAPT)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (p_ngb_old >= 0) then
              error stop &
                   "build_source_blocks: existing boundary is not negative"
           end if

           b_src = -p_ngb_old

           if (b_src >= grid(d)%bdry_patch%length) then
              error stop &
                   "build_source_blocks: invalid source bdry_patch"
           end if

           ib = ib + 1

           block_out%block_bdry(ib)%patch = old_to_new(p_old)
           block_out%block_bdry(ib)%side  = c

           block_out%block_bdry(ib)%class = &
                block_out%neigh_class(c,old_to_new(p_old)+1)

           block_out%block_bdry(ib)%root_domain = &
                block_catalog(b_catalog)%root_domain

           block_out%block_bdry(ib)%source_bdry = b_src

           block_out%block_bdry(ib)%elts_start = &
                grid(d)%bdry_patch%elts(b_src+1)%elts_start

           block_out%block_bdry(ib)%bdry_side = &
                grid(d)%bdry_patch%elts(b_src+1)%side

           block_out%block_bdry(ib)%bdry_neigh = &
                grid(d)%bdry_patch%elts(b_src+1)%neigh

           call get_bdry_dims_Domain( &
                grid(d), b_src, block_out%block_bdry(ib)%dims)

           block_out%block_bdry(ib)%n_node = &
                BDRY_THICKNESS * PATCH_SIZE

           block_out%block_bdry(ib)%storage_id = -1

        case default

           error stop &
                "build_source_blocks: unexpected neighbour class"

        end select

     end do

  end do

  if (ib /= n_bdry_local) then
     error stop &
          "build_source_blocks: local boundary count mismatch"
  end if

  !
  ! ===============================================================
  ! Map every inter-block link to its unique source patch, candidate
  ! source block, and prospective source owner.
  ! ===============================================================
  !
  n_block_source_local  = 0
  n_block_source_remote = 0

  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_block  = -1
     n_source_match = 0

     do i = 1, size(block_catalog)

        if (block_catalog(i)%root_domain /= &
             block_catalog(b_catalog)%root_domain) cycle

        if (.not. patch_in_subtree( &
             block_catalog(i)%root_patch, &
             block_out%block_bdry(ib)%neigh_patch)) cycle

        source_block   = i
        n_source_match = n_source_match + 1

     end do

     if (n_source_match /= 1) then
        error stop &
             "build_source_blocks: nonunique inter-block source"
     end if

     if (source_block == b_catalog) then
        error stop &
             "build_source_blocks: inter-block source is current block"
     end if

     block_out%block_bdry(ib)%source_block = source_block
     block_out%block_bdry(ib)%source_block_id = &
          block_catalog(source_block)%id
     block_out%block_bdry(ib)%source_owner = &
          block_catalog(source_block)%owner

     if (block_out%block_bdry(ib)%source_owner < 0 .or. &
          block_out%block_bdry(ib)%source_owner >= n_process) then

        error stop &
             "build_source_blocks: invalid inter-block source owner"

     end if

     if (block_out%block_bdry(ib)%source_owner == &
          block_catalog(b_catalog)%owner) then
        n_block_source_local = n_block_source_local + 1
     else
        n_block_source_remote = n_block_source_remote + 1
     end if

  end do

  if (n_block_source_local + n_block_source_remote /= &
       n_ngb_block) then

     error stop &
          "build_source_blocks: inter-block mapping count mismatch"

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
  allocate(bdry_level(max(1,grid(d)%bdry_patch%length)))

  bdry_required = -1
  bdry_closure  = -1
  bdry_level = -1

  n_bdry_direct   = 0
  n_bdry_closure  = 0
  n_bdry_required = 0

  !
  ! Distinct directly referenced boundaries.
  !
  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_DOMAIN .and. &
          block_out%block_bdry(ib)%class /= NGB_ADAPT) cycle

     b_src = block_out%block_bdry(ib)%source_bdry

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
        bdry_level(n_bdry_required) = &
             patch_copy(block_out%block_bdry(ib)%patch+1)%level

     else if (bdry_level(j) /= &
          patch_copy(block_out%block_bdry(ib)%patch+1)%level) then

        error stop &
             "build_source_blocks: direct boundary level differs"

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

        select case (block_out%neigh_class( &
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

        if (already_present) then

           if (bdry_level(j) /= grid(d)%patch%elts(p_old+1)%level) then
              error stop &
                   "build_source_blocks: closure boundary level differs"
           end if

           cycle

        end if

        n_bdry_required = n_bdry_required + 1
        n_bdry_closure  = n_bdry_closure + 1

        bdry_required(n_bdry_required) = b_missing
        bdry_closure(n_bdry_closure)   = b_missing
        bdry_level(n_bdry_required) = &
             grid(d)%patch%elts(p_old+1)%level

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
             "build_source_blocks: invalid closure boundary"

     end if

     do j = 1, i-1

        if (bdry_closure(j) == b_src) then
           error stop &
                "build_source_blocks: duplicate closure boundary"
        end if

     end do

     do j = 1, n_bdry_direct

        if (bdry_required(j) == b_src) then
           error stop &
                "build_source_blocks: closure boundary already direct"
        end if

     end do

  end do

  if (n_bdry_required /= n_bdry_direct + n_bdry_closure) then
     error stop &
          "build_source_blocks: required boundary count mismatch"
  end if

  !
  ! ===============================================================
  ! Construct complete compact boundary-storage catalogue.
  ! ===============================================================
  !
  allocate(block_out%bdry_storage(n_bdry_required))

  n_bdry_node_unique = 0

  do is = 1, n_bdry_required

     b_src = bdry_required(is)

     block_out%bdry_storage(is)%source_bdry = b_src
     block_out%bdry_storage(is)%level = bdry_level(is)

     block_out%bdry_storage(is)%elts_start = &
          grid(d)%bdry_patch%elts(b_src+1)%elts_start

     call get_bdry_dims_Domain( &
          grid(d), b_src, block_out%bdry_storage(is)%dims)

     block_out%bdry_storage(is)%n_node = &
          BDRY_THICKNESS * PATCH_SIZE

     if (block_out%bdry_storage(is)%n_node <= 0) then
        error stop &
             "build_source_blocks: invalid boundary storage size"
     end if

     block_out%bdry_storage(is)%local_start = &
          n_bdry_node_unique

     n_bdry_node_unique = n_bdry_node_unique + &
          block_out%bdry_storage(is)%n_node

  end do

  !
  ! Assign storage IDs to directly referenced boundary links.
  !
  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class == NGB_BLOCK) cycle

     do is = 1, size(block_out%bdry_storage)

        if (block_out%bdry_storage(is)%source_bdry == &
             block_out%block_bdry(ib)%source_bdry) then

           block_out%block_bdry(ib)%storage_id = is
           exit

        end if

     end do

     if (block_out%block_bdry(ib)%storage_id < 1) then
        error stop &
             "build_source_blocks: missing boundary storage ID"
     end if

  end do

  !
  ! Verify compact storage layout and uniqueness.
  !
  do is = 1, size(block_out%bdry_storage)

     if (is == 1) then

        if (block_out%bdry_storage(is)%local_start /= 0) then
           error stop &
                "build_source_blocks: invalid first boundary offset"
        end if

     else

        if (block_out%bdry_storage(is)%local_start /= &
             block_out%bdry_storage(is-1)%local_start + &
             block_out%bdry_storage(is-1)%n_node) then

           error stop &
                "build_source_blocks: noncompact boundary storage"

        end if

     end if

     if (block_out%bdry_storage(is)%elts_start < 0 .or. &
          block_out%bdry_storage(is)%elts_start + &
          block_out%bdry_storage(is)%n_node > &
          grid(d)%node%length) then

        error stop &
             "build_source_blocks: boundary source range invalid"

     end if

     do j = is+1, size(block_out%bdry_storage)

        if (block_out%bdry_storage(is)%source_bdry == &
             block_out%bdry_storage(j)%source_bdry) then

           error stop &
                "build_source_blocks: duplicate boundary storage"

        end if

     end do

  end do

  !
  ! ===============================================================
  ! Copy complete compact boundary data.
  ! ===============================================================
  !
  allocate(block_out%bdry_node(n_bdry_node_unique))

  scalar_storage_size = mult_scalar * n_bdry_node_unique
  scalar_variable_size = n_field_level * scalar_storage_size
  vector_storage_size = mult_vector * n_bdry_node_unique

  allocate(block_out%bdry_scalar( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%bdry_wavelet_scalar( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%bdry_scalar_mean( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%bdry_vector( &
       n_field_level*vector_storage_size))

  allocate(block_out%bdry_wavelet_vector( &
       n_field_level*vector_storage_size))

  allocate(block_out%bdry_vector_mean( &
       n_field_level*vector_storage_size))

  tke_storage_size = n_bdry_node_unique
  allocate(block_out%bdry_tke(n_tke_level*tke_storage_size))
  allocate(block_out%bdry_wavelet_tke( &
       n_tke_level*tke_storage_size))

  allocate(block_out%bdry_topography(n_bdry_node_unique))

  do is = 1, size(block_out%bdry_storage)

     old_start = block_out%bdry_storage(is)%elts_start
     new_start = block_out%bdry_storage(is)%local_start

     block_out%bdry_node( &
          new_start+1 : &
          new_start+block_out%bdry_storage(is)%n_node) = &
          grid(d)%node%elts( &
          old_start+1 : &
          old_start+block_out%bdry_storage(is)%n_node)

     block_out%bdry_topography( &
          new_start+1 : &
          new_start+block_out%bdry_storage(is)%n_node) = &
          topography%data(d)%elts( &
          old_start+1 : &
          old_start+block_out%bdry_storage(is)%n_node)

     do scalar_slot = 1, n_scalar_variable

        scalar_id = v_scalar + scalar_slot - 1

        do level_slot = 1, n_field_level

           field_level = first_field_level + level_slot - 1
           scalar_base = (scalar_slot-1) * scalar_variable_size + &
                (level_slot-1) * scalar_storage_size

           block_out%bdry_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) = &
                sol(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node))

           block_out%bdry_wavelet_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) = &
                wav_coeff(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node))

           block_out%bdry_scalar_mean( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) = &
                sol_mean(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node))

        end do

     end do

     do level_slot = 1, n_field_level

        field_level = first_field_level + level_slot - 1
        field_base = (level_slot-1) * vector_storage_size

        block_out%bdry_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) = &
             sol(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node))

        block_out%bdry_wavelet_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) = &
             wav_coeff(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node))

        block_out%bdry_vector_mean( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) = &
             sol_mean(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node))

     end do

     do level_slot = 1, n_tke_level

        field_level = tke_level + level_slot - 1
        field_base = (level_slot-1) * tke_storage_size

        block_out%bdry_tke( &
             field_base+new_start+1 : &
             field_base+new_start+ &
             block_out%bdry_storage(is)%n_node) = &
             tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+ &
             block_out%bdry_storage(is)%n_node)

        block_out%bdry_wavelet_tke( &
             field_base+new_start+1 : &
             field_base+new_start+ &
             block_out%bdry_storage(is)%n_node) = &
             wav_tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+ &
             block_out%bdry_storage(is)%n_node)

     end do

  end do

  !
  ! Verify copied boundary data.
  !
  do is = 1, size(block_out%bdry_storage)

     old_start = block_out%bdry_storage(is)%elts_start
     new_start = block_out%bdry_storage(is)%local_start

     if (maxval(abs( &
          block_out%bdry_topography( &
          new_start+1 : &
          new_start+block_out%bdry_storage(is)%n_node) - &
          topography%data(d)%elts( &
          old_start+1 : &
          old_start+block_out%bdry_storage(is)%n_node))) > 0.0_dp) then
        error stop "build_source_blocks: boundary topography mismatch"
     end if

     if (maxval(abs( &
          block_out%bdry_node( &
          new_start+1 : &
          new_start+block_out%bdry_storage(is)%n_node)%x - &
          grid(d)%node%elts( &
          old_start+1 : &
          old_start+block_out%bdry_storage(is)%n_node)%x)) > &
          0.0_dp) then

        error stop &
             "build_source_blocks: boundary coordinate mismatch"

     end if

     do scalar_slot = 1, n_scalar_variable

        scalar_id = v_scalar + scalar_slot - 1

        do level_slot = 1, n_field_level

           field_level = first_field_level + level_slot - 1
           scalar_base = (scalar_slot-1) * scalar_variable_size + &
                (level_slot-1) * scalar_storage_size

           if (maxval(abs( &
                block_out%bdry_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) - &
                sol(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: boundary scalar mismatch"

           end if

           if (maxval(abs( &
                block_out%bdry_wavelet_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) - &
                wav_coeff(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: boundary scalar wavelet mismatch"

           end if

           if (maxval(abs( &
                block_out%bdry_scalar_mean( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start + &
                block_out%bdry_storage(is)%n_node)) - &
                sol_mean(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start + &
                block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: boundary scalar mean mismatch"

           end if

        end do

     end do

     do level_slot = 1, n_field_level

        field_level = first_field_level + level_slot - 1
        field_base = (level_slot-1) * vector_storage_size

        if (maxval(abs( &
             block_out%bdry_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) - &
             sol(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: boundary vector mismatch"

        end if

        if (maxval(abs( &
             block_out%bdry_wavelet_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) - &
             wav_coeff(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: boundary vector wavelet mismatch"

        end if

        if (maxval(abs( &
             block_out%bdry_vector_mean( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start + &
             block_out%bdry_storage(is)%n_node)) - &
             sol_mean(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start + &
             block_out%bdry_storage(is)%n_node)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: boundary vector mean mismatch"

        end if

     end do

     do level_slot = 1, n_tke_level

        field_level = tke_level + level_slot - 1
        field_base = (level_slot-1) * tke_storage_size

        if (maxval(abs( &
             block_out%bdry_tke( &
             field_base+new_start+1 : &
             field_base+new_start+ &
             block_out%bdry_storage(is)%n_node) - &
             tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+ &
             block_out%bdry_storage(is)%n_node))) > 0.0_dp) then
           error stop "build_source_blocks: boundary tke mismatch"
        end if

        if (maxval(abs( &
             block_out%bdry_wavelet_tke( &
             field_base+new_start+1 : &
             field_base+new_start+ &
             block_out%bdry_storage(is)%n_node) - &
             wav_tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+ &
             block_out%bdry_storage(is)%n_node))) > 0.0_dp) then
           error stop "build_source_blocks: boundary wav_tke mismatch"
        end if

     end do

  end do

  !
  ! Boundary-storage statistics for directly referenced links.
  !
  n_bdry_node_total = 0
  n_bdry_node_max   = 0

  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_DOMAIN .and. &
          block_out%block_bdry(ib)%class /= NGB_ADAPT) cycle

     n_bdry_node_total = n_bdry_node_total + &
          block_out%block_bdry(ib)%n_node

     n_bdry_node_max = max( &
          n_bdry_node_max, &
          block_out%block_bdry(ib)%n_node)

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

        select case (block_out%neigh_class( &
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

        do is = 1, size(block_out%bdry_storage)

           old_start = block_out%bdry_storage(is)%elts_start

           if (idx_src >= old_start .and. &
                idx_src < old_start + &
                block_out%bdry_storage(is)%n_node) then

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
                "build_source_blocks: ghost source is inside block"
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
  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_patch = block_out%block_bdry(ib)%neigh_patch
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
  allocate(block_out%ghost_storage(n_ghost))

  n_ghost_node = 0

  do ghost_id = 1, n_ghost

     source_block   = -1
     n_source_match = 0

     do i = 1, size(block_catalog)

        if (block_catalog(i)%root_domain /= &
             block_catalog(b_catalog)%root_domain) cycle

        if (.not. patch_in_subtree( &
             block_catalog(i)%root_patch, &
             ghost_patch(ghost_id))) cycle

        source_block   = i
        n_source_match = n_source_match + 1

     end do

     if (n_source_match /= 1) then
        error stop &
             "build_source_blocks: nonunique ghost source block"
     end if

     if (source_block == b_catalog) then
        error stop &
             "build_source_blocks: ghost source is current block"
     end if

     next_source_patch = 0
     source_local_patch = -1

     call find_subtree_patch_index( &
          block_catalog(source_block)%root_patch, &
          ghost_patch(ghost_id),next_source_patch,source_local_patch)

     if (source_local_patch < 0) then
        error stop &
             "build_source_blocks: invalid compact ghost source patch"
     end if

     source_patch_count = count_subtree_patches_Domain( &
          grid(d),block_catalog(source_block)%root_patch)

     if (source_local_patch >= source_patch_count) then
        error stop &
             "build_source_blocks: compact ghost source patch out of range"
     end if

     block_out%ghost_storage(ghost_id)%source_domain = &
          block_catalog(b_catalog)%root_domain

     block_out%ghost_storage(ghost_id)%source_patch = &
          ghost_patch(ghost_id)

     block_out%ghost_storage(ghost_id)%source_block = &
          source_block

     block_out%ghost_storage(ghost_id)%source_block_id = &
          block_catalog(source_block)%id

     block_out%ghost_storage(ghost_id)%source_owner = &
          block_catalog(source_block)%owner

     block_out%ghost_storage(ghost_id)%source_local_patch = &
          source_local_patch

     if (block_out%ghost_storage(ghost_id)%source_owner < 0 .or. &
          block_out%ghost_storage(ghost_id)%source_owner >= &
          n_process) then

        error stop &
             "build_source_blocks: invalid ghost source owner"

     end if

     block_out%ghost_storage(ghost_id)%elts_start = &
          grid(d)%patch%elts(ghost_patch(ghost_id)+1)%elts_start

     block_out%ghost_storage(ghost_id)%local_start = &
          n_ghost_node

     block_out%ghost_storage(ghost_id)%n_node = &
          PATCH_SIZE**2

     n_ghost_node = n_ghost_node + PATCH_SIZE**2

  end do

  !
  ! Assign every nominal inter-block link its compact ghost ID.
  !
  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     do ghost_id = 1, n_ghost

        if (block_out%ghost_storage(ghost_id)%source_patch == &
             block_out%block_bdry(ib)%neigh_patch) then

           block_out%block_bdry(ib)%ghost_id = ghost_id
           exit

        end if

     end do

     if (block_out%block_bdry(ib)%ghost_id < 1) then
        error stop &
             "build_source_blocks: missing inter-block ghost ID"
     end if

     ghost_id = block_out%block_bdry(ib)%ghost_id

     if (block_out%ghost_storage(ghost_id)%source_block /= &
          block_out%block_bdry(ib)%source_block .or. &
          block_out%ghost_storage(ghost_id)%source_block_id /= &
          block_out%block_bdry(ib)%source_block_id .or. &
          block_out%ghost_storage(ghost_id)%source_owner /= &
          block_out%block_bdry(ib)%source_owner) then

        error stop &
             "build_source_blocks: ghost source metadata mismatch"

     end if

  end do

  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     source_block = block_out%block_bdry(ib)%source_block

     if (source_block < 1 .or. &
          source_block > size(block_catalog)) then

        error stop &
             "build_source_blocks: invalid source catalog index"

     end if

     if (block_out%block_bdry(ib)%source_block_id /= &
          block_catalog(source_block)%id) then

        error stop &
             "build_source_blocks: source block ID mismatch"

     end if

  end do

  do ghost_id = 1, size(block_out%ghost_storage)

     source_block = &
          block_out%ghost_storage(ghost_id)%source_block

     if (source_block < 1 .or. &
          source_block > size(block_catalog)) then

        error stop &
             "build_source_blocks: invalid ghost catalog index"

     end if

     if (block_out%ghost_storage(ghost_id)%source_block_id /= &
          block_catalog(source_block)%id) then

        error stop &
             "build_source_blocks: ghost block ID mismatch"

     end if

  end do

  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     do j = ib+1, size(block_out%block_bdry)

        if (block_out%block_bdry(j)%class /= NGB_BLOCK) cycle

        if ((block_out%block_bdry(ib)%neigh_patch == &
             block_out%block_bdry(j)%neigh_patch) .neqv. &
             (block_out%block_bdry(ib)%ghost_id == &
             block_out%block_bdry(j)%ghost_id)) then

           error stop &
                "build_source_blocks: inconsistent ghost deduplication"

        end if

     end do

  end do

  allocate(block_out%ghost_node(n_ghost_node))

  scalar_storage_size = mult_scalar * n_ghost_node
  scalar_variable_size = n_field_level * scalar_storage_size
  vector_storage_size = mult_vector * n_ghost_node

  allocate(block_out%ghost_scalar( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%ghost_wavelet_scalar( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%ghost_scalar_mean( &
       n_scalar_variable*scalar_variable_size))

  allocate(block_out%ghost_vector( &
       n_field_level*vector_storage_size))

  allocate(block_out%ghost_wavelet_vector( &
       n_field_level*vector_storage_size))

  allocate(block_out%ghost_vector_mean( &
       n_field_level*vector_storage_size))

  tke_storage_size = n_ghost_node
  allocate(block_out%ghost_tke(n_tke_level*tke_storage_size))
  allocate(block_out%ghost_wavelet_tke( &
       n_tke_level*tke_storage_size))

  allocate(block_out%ghost_topography(n_ghost_node))

  !
  ! Copy temporary ghost data from the source domain. Eventually the
  ! field data will be supplied by inter-block communication.
  !
  do ghost_id = 1, n_ghost

     old_start = &
          block_out%ghost_storage(ghost_id)%elts_start

     new_start = &
          block_out%ghost_storage(ghost_id)%local_start

     block_out%ghost_node( &
          new_start+1 : new_start+PATCH_SIZE**2) = &
          grid(d)%node%elts( &
          old_start+1 : old_start+PATCH_SIZE**2)

     block_out%ghost_topography( &
          new_start+1 : new_start+PATCH_SIZE**2) = &
          topography%data(d)%elts( &
          old_start+1 : old_start+PATCH_SIZE**2)

     do scalar_slot = 1, n_scalar_variable

        scalar_id = v_scalar + scalar_slot - 1

        do level_slot = 1, n_field_level

           field_level = first_field_level + level_slot - 1
           scalar_base = (scalar_slot-1) * scalar_variable_size + &
                (level_slot-1) * scalar_storage_size

           block_out%ghost_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) = &
                sol(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2))

           block_out%ghost_wavelet_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) = &
                wav_coeff(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2))

           block_out%ghost_scalar_mean( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) = &
                sol_mean(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2))

        end do

     end do

     do level_slot = 1, n_field_level

        field_level = first_field_level + level_slot - 1
        field_base = (level_slot-1) * vector_storage_size

        block_out%ghost_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) = &
             sol(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2))

        block_out%ghost_wavelet_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) = &
             wav_coeff(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2))

        block_out%ghost_vector_mean( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) = &
             sol_mean(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2))

     end do

     do level_slot = 1, n_tke_level

        field_level = tke_level + level_slot - 1
        field_base = (level_slot-1) * tke_storage_size

        block_out%ghost_tke( &
             field_base+new_start+1 : &
             field_base+new_start+PATCH_SIZE**2) = &
             tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+PATCH_SIZE**2)

        block_out%ghost_wavelet_tke( &
             field_base+new_start+1 : &
             field_base+new_start+PATCH_SIZE**2) = &
             wav_tke(field_level)%data(d)%elts( &
             old_start+1 : old_start+PATCH_SIZE**2)

     end do

  end do

  !
  ! Verify temporary ghost copies.
  !
  do ghost_id = 1, n_ghost

     old_start = &
          block_out%ghost_storage(ghost_id)%elts_start

     new_start = &
          block_out%ghost_storage(ghost_id)%local_start

     if (maxval(abs( &
          block_out%ghost_topography( &
          new_start+1 : new_start+PATCH_SIZE**2) - &
          topography%data(d)%elts( &
          old_start+1 : old_start+PATCH_SIZE**2))) > 0.0_dp) then
        error stop "build_source_blocks: ghost topography mismatch"
     end if

     if (maxval(abs( &
          block_out%ghost_node( &
          new_start+1 : new_start+PATCH_SIZE**2)%x - &
          grid(d)%node%elts( &
          old_start+1 : old_start+PATCH_SIZE**2)%x)) > 0.0_dp) then

        error stop &
             "build_source_blocks: ghost coordinate mismatch"

     end if

     do scalar_slot = 1, n_scalar_variable

        scalar_id = v_scalar + scalar_slot - 1

        do level_slot = 1, n_field_level

           field_level = first_field_level + level_slot - 1
           scalar_base = (scalar_slot-1) * scalar_variable_size + &
                (level_slot-1) * scalar_storage_size

           if (maxval(abs( &
                block_out%ghost_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) - &
                sol(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: ghost scalar mismatch"

           end if

           if (maxval(abs( &
                block_out%ghost_wavelet_scalar( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) - &
                wav_coeff(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: ghost scalar wavelet mismatch"

           end if

           if (maxval(abs( &
                block_out%ghost_scalar_mean( &
                scalar_base+mult_scalar*new_start+1 : &
                scalar_base+mult_scalar*(new_start+PATCH_SIZE**2)) - &
                sol_mean(scalar_id,field_level)%data(d)%elts( &
                mult_scalar*old_start+1 : &
                mult_scalar*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

              error stop &
                   "build_source_blocks: ghost scalar mean mismatch"

           end if

        end do

     end do

     do level_slot = 1, n_field_level

        field_level = first_field_level + level_slot - 1
        field_base = (level_slot-1) * vector_storage_size

        if (maxval(abs( &
             block_out%ghost_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) - &
             sol(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: ghost vector mismatch"

        end if

        if (maxval(abs( &
             block_out%ghost_wavelet_vector( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) - &
             wav_coeff(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: ghost vector wavelet mismatch"

        end if

        if (maxval(abs( &
             block_out%ghost_vector_mean( &
             field_base+mult_vector*new_start+1 : &
             field_base+mult_vector*(new_start+PATCH_SIZE**2)) - &
             sol_mean(v_vector,field_level)%data(d)%elts( &
             mult_vector*old_start+1 : &
             mult_vector*(old_start+PATCH_SIZE**2)))) > 0.0_dp) then

           error stop &
                "build_source_blocks: ghost vector mean mismatch"

        end if

     end do

     do level_slot = 1, n_tke_level

        field_level = tke_level + level_slot - 1
        field_base = (level_slot-1) * tke_storage_size

        if (maxval(abs( &
             block_out%ghost_tke( &
             field_base+new_start+1 : &
             field_base+new_start+PATCH_SIZE**2) - &
             tke(field_level)%data(d)%elts( &
             old_start+1:old_start+PATCH_SIZE**2))) > 0.0_dp) then
           error stop "build_source_blocks: ghost tke mismatch"
        end if

        if (maxval(abs( &
             block_out%ghost_wavelet_tke( &
             field_base+new_start+1 : &
             field_base+new_start+PATCH_SIZE**2) - &
             wav_tke(field_level)%data(d)%elts( &
             old_start+1:old_start+PATCH_SIZE**2))) > 0.0_dp) then
           error stop "build_source_blocks: ghost wav_tke mismatch"
        end if

     end do

  end do

  !
  ! ===============================================================
  ! Renumber neighbours and convert all external links to local
  ! block_bdry references.
  ! ===============================================================
  !
  call renumber_subtree_neigh_Domain( &
       grid(d), patch_copy, old_to_new)

  do ib = 1, size(block_out%block_bdry)

     patch_copy(block_out%block_bdry(ib)%patch+1)%neigh( &
          block_out%block_bdry(ib)%side) = -ib

  end do

  !
  ! Verify local neighbour representation.
  !
  do p_old = 0, grid(d)%patch%length-1

     if (old_to_new(p_old) < 0) cycle

     do c = 1, N_BDRY

        select case (block_out%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_INTERNAL)

           p_ngb_old = grid(d)%patch%elts(p_old+1)%neigh(c)

           if (patch_copy(old_to_new(p_old)+1)%neigh(c) /= &
                old_to_new(p_ngb_old)) then

              error stop &
                   "build_source_blocks: internal neighbour mismatch"

           end if

        case (NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT)

           ib = -patch_copy(old_to_new(p_old)+1)%neigh(c)

           if (ib < 1 .or. &
                ib > size(block_out%block_bdry)) then

              error stop &
                   "build_source_blocks: invalid boundary reference"

           end if

           if (block_out%block_bdry(ib)%patch /= &
                old_to_new(p_old) .or. &
                block_out%block_bdry(ib)%side /= c) then

              error stop &
                   "build_source_blocks: incorrect boundary reference"

           end if

        case (NGB_OTHER)

           if (patch_copy(old_to_new(p_old)+1)%neigh(c) /= 0) then
              error stop &
                   "build_source_blocks: zero neighbour changed"
           end if

        case default

           error stop &
                "build_source_blocks: unexpected neighbour class"

        end select

     end do

  end do

  !
  ! ===============================================================
  ! Package extracted interior arrays.
  ! ===============================================================
  !
  block_out%id          = block_catalog(b_catalog)%id
  block_out%root_domain = block_catalog(b_catalog)%root_domain
  block_out%root_patch  = block_catalog(b_catalog)%root_patch
  block_out%level       = block_catalog(b_catalog)%level
  block_out%scalar_variable = v_scalar
  block_out%n_scalar_variable = n_scalar_variable
  block_out%vector_variable = v_vector
  block_out%field_level     = first_field_level
  block_out%n_field_level   = n_field_level
  block_out%scalar_mult     = mult_scalar
  block_out%vector_mult     = mult_vector
  block_out%tke_level       = tke_level
  block_out%n_tke_level     = n_tke_level

  call move_alloc(patch_copy,  block_out%patch)
  call move_alloc(node_copy,   block_out%node)
  call move_alloc(scalar_copy, block_out%scalar)
  call move_alloc(wavelet_scalar_copy, block_out%wavelet_scalar)
  call move_alloc(scalar_mean_copy, block_out%scalar_mean)
  call move_alloc(vector_copy, block_out%vector)
  call move_alloc(wavelet_vector_copy, block_out%wavelet_vector)
  call move_alloc(vector_mean_copy, block_out%vector_mean)
  call move_alloc(tke_copy, block_out%tke)
  call move_alloc(wavelet_tke_copy, block_out%wavelet_tke)
  call move_alloc(topography_copy, block_out%topography)

  !
  ! ===============================================================
  ! Build explicit compact stencil addresses.
  !
  ! NGB_BLOCK addresses use the unified compact ghost catalogue and
  ! the orientation-adjusted effective address from comp_offs3.
  ! ===============================================================
  !
  allocate(block_out%stencil( &
       N_BDRY,size(block_out%patch)))

  block_out%stencil%storage = STORE_NONE
  block_out%stencil%id      = -1
  block_out%stencil%offset  = 0

  do i = 1, size(block_out%stencil,2)
     do c = 1, size(block_out%stencil,1)
        block_out%stencil(c,i)%dims = 0
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

        select case (block_out%neigh_class( &
             c,old_to_new(p_old)+1))

        case (NGB_OTHER)

           cycle

        case (NGB_INTERNAL, NGB_BLOCK, NGB_DOMAIN, NGB_ADAPT)

           idx_src = grid(d)%patch%elts(p_old+1)%elts_start + &
                offs_src(c)

        case default

           error stop &
                "build_source_blocks: unexpected stencil class"

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

           do is = 1, size(block_out%bdry_storage)

              old_start = block_out%bdry_storage(is)%elts_start

              if (idx_src >= old_start .and. &
                   idx_src < old_start + &
                   block_out%bdry_storage(is)%n_node) then

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

           do i = 1, size(block_out%ghost_storage)

              old_start = &
                   block_out%ghost_storage(i)%elts_start

              if (idx_src >= old_start .and. &
                   idx_src < old_start + &
                   block_out%ghost_storage(i)%n_node) then

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
        if (block_out%neigh_class( &
             c,old_to_new(p_old)+1) == NGB_BLOCK) then

           ib = -block_out%patch( &
                old_to_new(p_old)+1)%neigh(c)

           if (ib < 1 .or. &
                ib > size(block_out%block_bdry)) then

              error stop &
                   "build_source_blocks: invalid NGB_BLOCK record"

           end if

           if (block_out%block_bdry(ib)%class /= NGB_BLOCK) then
              error stop &
                   "build_source_blocks: incorrect NGB_BLOCK record"
           end if

           ghost_id = block_out%block_bdry(ib)%ghost_id

           if (ghost_id < 1 .or. &
                ghost_id > size(block_out%ghost_storage)) then

              error stop &
                   "build_source_blocks: invalid NGB_BLOCK ghost ID"

           end if

           if (block_out%ghost_storage(ghost_id)%source_patch /= &
                block_out%block_bdry(ib)%neigh_patch) then

              error stop &
                   "build_source_blocks: incorrect NGB_BLOCK ghost"

           end if

           target_patch = -1
           target_bdry  = -1

           old_start = &
                block_out%ghost_storage(ghost_id)%elts_start

           target_offset = idx_src - old_start

           if (target_offset > &
                block_out%ghost_storage(ghost_id)%n_node-1 .or. &
                target_offset + PATCH_SIZE**2-1 < 0) then

              error stop &
                   "build_source_blocks: empty NGB_BLOCK mapping"

           end if

        end if

        !
        ! Resolved interior target.
        !
        if (target_patch >= 0) then

           block_out%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_PATCH

           block_out%stencil( &
                c,old_to_new(p_old)+1)%id = target_patch

           block_out%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_out%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_patch = n_stencil_patch + 1
           n_stencil_built = n_stencil_built + 1

        !
        ! Resolved compact boundary target.
        !
        else if (target_bdry >= 0) then

           block_out%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_BDRY

           block_out%stencil( &
                c,old_to_new(p_old)+1)%id = target_bdry

           block_out%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_out%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_bdry  = n_stencil_bdry + 1
           n_stencil_built = n_stencil_built + 1

        !
        ! Resolved compact ghost-patch target.
        !
        else if (ghost_id >= 0) then

           block_out%stencil( &
                c,old_to_new(p_old)+1)%storage = STORE_GHOST

           block_out%stencil( &
                c,old_to_new(p_old)+1)%id = ghost_id

           block_out%stencil( &
                c,old_to_new(p_old)+1)%offset = target_offset

           block_out%stencil( &
                c,old_to_new(p_old)+1)%dims = dims_src(:,c)

           n_stencil_ghost = n_stencil_ghost + 1
           n_stencil_built = n_stencil_built + 1

           if (block_out%neigh_class( &
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
        select case (block_out%stencil( &
             c,old_to_new(p_old)+1)%storage)

        case (STORE_PATCH)

           target_patch = block_out%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_out%stencil( &
                c,old_to_new(p_old)+1)%offset

           val_src = sol(v_scalar,first_field_level)%data(d)%elts( &
                mult_scalar*idx_src + 1)

           val_blk = block_out%scalar( &
                mult_scalar * &
                (block_out%patch(target_patch+1)%elts_start + &
                target_offset) + 1)

        case (STORE_BDRY)

           target_bdry = block_out%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_out%stencil( &
                c,old_to_new(p_old)+1)%offset

           val_src = sol(v_scalar,first_field_level)%data(d)%elts( &
                mult_scalar*idx_src + 1)

           val_blk = block_out%bdry_scalar( &
                mult_scalar * &
                (block_out%bdry_storage(target_bdry)%local_start + &
                target_offset) + 1)

        case (STORE_GHOST)

           ghost_id = block_out%stencil( &
                c,old_to_new(p_old)+1)%id

           target_offset = block_out%stencil( &
                c,old_to_new(p_old)+1)%offset

           if (block_out%neigh_class( &
                c,old_to_new(p_old)+1) == NGB_BLOCK) then

              n_mapped = 0

              do q = 0, PATCH_SIZE**2-1

                 ghost_offset = target_offset + q

                 if (ghost_offset < 0 .or. &
                      ghost_offset >= &
                      block_out%ghost_storage(ghost_id)%n_node) cycle

                 val_src = sol(v_scalar,first_field_level)%data(d)%elts( &
                      mult_scalar*(idx_src+q) + 1)

                 val_blk = block_out%ghost_scalar( &
                      mult_scalar * &
                      (block_out%ghost_storage( &
                      ghost_id)%local_start + ghost_offset) + 1)

                 if (abs(val_blk-val_src) > 0.0_dp) then
                    error stop &
                         "build_source_blocks: NGB_BLOCK value mismatch"
                 end if

                 n_mapped = n_mapped + 1

              end do

              if (n_mapped < 1) then
                 error stop &
                      "build_source_blocks: empty NGB_BLOCK validation"
              end if

              n_block_value_checked = &
                   n_block_value_checked + n_mapped

              ! Suppress the single-address comparison below; every
              ! physically mapped value has already been checked.
              val_src = 0.0_dp
              val_blk = 0.0_dp

           else

              val_src = sol(v_scalar,first_field_level)%data(d)%elts( &
                   mult_scalar*idx_src + 1)

              val_blk = block_out%ghost_scalar( &
                   mult_scalar * &
                   (block_out%ghost_storage( &
                   ghost_id)%local_start + target_offset) + 1)

           end if

        case default

           error stop &
                "build_source_blocks: invalid stencil storage"

        end select

        if (abs(val_blk-val_src) > 0.0_dp) then
           error stop &
                "build_source_blocks: explicit scalar stencil mismatch"
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
          "build_source_blocks: stencil count mismatch"

  end if

  if (n_stencil_patch + n_stencil_bdry + n_stencil_ghost /= &
       n_stencil_built) then

     error stop &
          "build_source_blocks: stencil class count mismatch"

  end if

  if (n_stencil_block /= n_ngb_block) then
     error stop &
          "build_source_blocks: NGB_BLOCK stencil count mismatch"
  end if

  if (n_block_value_checked < n_stencil_block) then
     error stop &
          "build_source_blocks: incomplete NGB_BLOCK value check"
  end if

  if (n_unresolved_patch + n_unresolved_bdry + &
       n_unresolved_none /= n_stencil_unresolved) then

     error stop &
          "build_source_blocks: unresolved count mismatch"

  end if

  !
  ! The boundary closure must now be complete.
  !
  if (n_unresolved_bdry /= 0) then
     error stop &
          "build_source_blocks: boundary-storage closure incomplete"
  end if

  !
  ! ===============================================================
  ! Diagnostic output.
  ! ===============================================================
  !
  if (verbose) then

  write(6,'(/,a,i0,a,i0)') &
       "Source block extraction: domain ", &
       block_catalog(b_catalog)%root_domain, &
       ", root patch ", p_root

  write(6,'(a,i0)') &
       "  block ID = ", block_catalog(b_catalog)%id

  write(6,'(a,i0)') &
       "  block weight = ", block_catalog(b_catalog)%weight

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
       "  interior coordinate/sol/sol_mean/wav_coeff copy checks passed"

  write(6,'(a)') &
       "  interior topography copy checks passed"

  if (n_tke_level > 0) then
     write(6,'(a)') &
          "  interior tke/wav_tke copy checks passed"
  end if

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

  do ib = 1, size(block_out%block_bdry)

     if (block_out%block_bdry(ib)%class /= NGB_BLOCK) cycle

     write(6,'(a,i0,a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "    link ", ib, &
          ": source patch = ", &
          block_out%block_bdry(ib)%neigh_patch, &
          ", block index = ", &
          block_out%block_bdry(ib)%source_block, &
          ", block ID = ", &
          block_out%block_bdry(ib)%source_block_id, &
          ", owner = ", &
          block_out%block_bdry(ib)%source_owner, &
          ", ghost ID = ", &
          block_out%block_bdry(ib)%ghost_id

  end do

  write(6,'(a)') &
       "  unique inter-block source mapping check passed"

  write(6,'(/,a,i0)') &
       "  total local boundary records = ", &
       size(block_out%block_bdry)

  write(6,'(a,i0)') &
       "  directly referenced boundary patches = ", &
       n_bdry_direct

  write(6,'(a,i0)') &
       "  stencil-closure boundary patches = ", &
       n_bdry_closure

  write(6,'(a,i0)') &
       "  total compact boundary patches = ", &
       size(block_out%bdry_storage)

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
       "  boundary coordinate/sol/sol_mean/wav_coeff copy checks passed"

  write(6,'(a)') &
       "  boundary topography copy checks passed"

  if (n_tke_level > 0) then
     write(6,'(a)') &
          "  boundary tke/wav_tke copy checks passed"
  end if

  write(6,'(/,a,i0)') &
       "  compact ghost source patches = ", &
       size(block_out%ghost_storage)

  write(6,'(a,i0)') &
       "  compact ghost node storage = ", &
       size(block_out%ghost_node)

  do ghost_id = 1, size(block_out%ghost_storage)

     write(6,'(a,i0,a,i0,a,i0,a,i0,a,i0)') &
          "    ghost ", ghost_id, &
          ": source patch = ", &
          block_out%ghost_storage(ghost_id)%source_patch, &
          ", block index = ", &
          block_out%ghost_storage(ghost_id)%source_block, &
          ", block ID = ", &
          block_out%ghost_storage(ghost_id)%source_block_id, &
          ", owner = ", &
          block_out%ghost_storage(ghost_id)%source_owner

  end do

  write(6,'(a)') &
       "  unified ghost sol/sol_mean/wav_coeff copy checks passed"

  write(6,'(a)') &
       "  unified ghost topography copy checks passed"

  if (n_tke_level > 0) then
     write(6,'(a)') &
          "  unified ghost tke/wav_tke copy checks passed"
  end if

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
  deallocate(bdry_level)
  deallocate(old_to_new)
  deallocate(old_elts_start)

contains

  recursive subroutine find_subtree_patch_index ( &
       p,p_target,next_patch,local_patch)
    ! Reproduce extract_subtree_patches_Domain preorder numbering for a
    ! target patch in another candidate block of the same root domain.

    implicit none

    integer, intent(in) :: p
    integer, intent(in) :: p_target
    integer, intent(inout) :: next_patch
    integer, intent(inout) :: local_patch

    integer :: c
    integer :: p_chd

    if (local_patch >= 0) return

    if (p < 0 .or. p >= grid(d)%patch%length) then
       error stop "find_subtree_patch_index: invalid patch"
    end if

    if (grid(d)%patch%elts(p+1)%deleted) return

    if (p == p_target) then
       local_patch = next_patch
       return
    end if

    next_patch = next_patch + 1

    do c = 1, N_CHDRN
       p_chd = grid(d)%patch%elts(p+1)%children(c)
       if (p_chd <= 0) cycle
       call find_subtree_patch_index( &
            p_chd,p_target,next_patch,local_patch)
       if (local_patch >= 0) return
    end do

  end subroutine find_subtree_patch_index

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

end subroutine build_one_source_block


end module parallel_block_build_mod
