! Historical full-record update retained as an independent regression reference.
    subroutine reference_boundary_record ( &
         d,source_bdry,boundary_level,n_node,scalar_slot,level_slot,grid_level, &
         capture_direct,local_index,n_boundary_node,node_start)

      implicit none

      integer, intent(in) :: d
      integer, intent(in) :: source_bdry
      integer, intent(in) :: boundary_level
      integer, intent(in) :: n_node
      integer, intent(in) :: scalar_slot
      integer, intent(in) :: level_slot
      integer, intent(in) :: grid_level
      logical, intent(in) :: capture_direct
      integer, intent(in) :: local_index
      integer, intent(in) :: n_boundary_node
      integer, intent(in) :: node_start

      integer :: data_start
      integer :: id
      integer :: node
      integer :: sample

      if (boundary_level /= grid_level) return
      id = grid(d)%bdry_patch%elts(source_bdry+1)%elts_start
      do node = 0,n_node-1
         sample = ((scalar_slot-1)*n_field_level + level_slot-1)* &
              n_boundary_node + node_start + node
         data_start = BLOCK_SCALAR_DIVERGENCE_INPUT_COUNT*sample + 1
         block
           real(dp) :: record(BLOCK_SCALAR_DIVERGENCE_INPUT_COUNT)
           record=scalar_read_range(block_scalar_tendency(local_index)%bdry, &
                data_start,data_start+BLOCK_SCALAR_DIVERGENCE_INPUT_COUNT-1)
           call fill_boundary_node(d,id+node,v_scalar+scalar_slot-1,capture_direct,record)
           call scalar_write_range(block_scalar_tendency(local_index)%bdry, &
                data_start,data_start+BLOCK_SCALAR_DIVERGENCE_INPUT_COUNT-1,record)
         end block
      end do

    end subroutine reference_boundary_record
