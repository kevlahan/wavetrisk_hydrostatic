module parallel_block_rhs_layout_mod
  ! Compile the actual resident primitive input footprint. No state/geometry copies.
  ! This is address infrastructure; it does not claim derived-field or mask readiness.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use shared_mod, only : N_BDRY, EDGE, S_MASS, S_TEMP, WEST, EAST, SOUTH, NORTH, &
       SOUTHWEST, NORTHWEST, SOUTHEAST, NORTHEAST
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data, STORE_PATCH, STORE_BDRY, STORE_GHOST
  implicit none
  private
  public :: Rhs_Input_Layout, compile_rhs_input_layout, rhs_input_values, rhs_input_layout_bytes
  public :: rhs_layout_slot, rhs_layout_extent
  integer, parameter :: dx(7)=[0,0,1,0,-1,1,-1],dy(7)=[0,1,0,-1,0,1,-1]
  integer, parameter :: vnode(10)=[1,2,5,1,4,7,5,1,3,4],vedge(10)=[1,1,1,2,2,2,2,3,3,3]
  type :: Rhs_Input_Layout
     private
     integer(int64) :: generation=-1_int64
     integer :: identity(8)=0
     integer :: nodes(3)=0
     integer, allocatable :: patch(:,:),boundary(:,:),ghost(:,:)
     ! One shared seven-node footprint for all layers, scalars and velocity components.
     integer, allocatable :: storage(:,:),node(:,:)
     integer(int64) :: rebuilds=0_int64,reuses=0_int64
  end type
contains
  subroutine compile_rhs_input_layout(block,generation,plan,rebuilt)
    type(Block_Data), intent(in) :: block
    integer(int64), intent(in) :: generation
    type(Rhs_Input_Layout), intent(inout) :: plan
    logical, intent(out) :: rebuilt
    integer :: identity(8),nodes(3),p,i,j,n,t
    integer, allocatable :: patches(:,:),boundaries(:,:),ghosts(:,:)
    if (generation<0_int64) error stop 'RHS layout generation invalid'
    if (.not.allocated(block%patch).or..not.allocated(block%node)) error stop 'RHS layout missing owned topology'
    if (block%scalar_mult/=1.or.block%vector_mult/=EDGE) error stop 'RHS layout field multiplicity invalid'
    identity=[block%id,block%root_domain,block%root_patch,block%scalar_variable,block%n_scalar_variable, &
         block%field_level,block%n_field_level,block%vector_variable]
    nodes=[size(block%node),0,0]
    if (allocated(block%bdry_node)) nodes(2)=size(block%bdry_node)
    if (allocated(block%ghost_node)) nodes(3)=size(block%ghost_node)
    call validate_fields(block,nodes)
    ! Exact structural snapshots catch changed mappings even if a caller reuses an epoch.
    allocate(patches(N_BDRY+2,size(block%patch)))
    do p=1,size(block%patch)
       if (block%patch(p)%deleted) error stop 'RHS layout contains deleted owned patch'
       patches(:,p)=[block%patch(p)%elts_start,block%patch(p)%level,block%patch(p)%neigh]
    end do
    n=0
    if (allocated(block%block_bdry)) n=size(block%block_bdry)
    allocate(boundaries(6,n))
    do p=1,n
       boundaries(:,p)=[block%block_bdry(p)%storage_id,block%block_bdry(p)%ghost_id, &
            block%block_bdry(p)%dims,0,block%block_bdry(p)%bdry_side]
       i=block%block_bdry(p)%ghost_id
       if (i>0) then
          if (.not.allocated(block%ghost_storage)) error stop 'RHS layout missing ghost descriptor'
          if (i>size(block%ghost_storage)) error stop 'RHS layout ghost descriptor invalid'
          boundaries(5,p)=block%ghost_storage(i)%local_start
       else
          i=block%block_bdry(p)%storage_id
          if (.not.allocated(block%bdry_storage)) error stop 'RHS layout missing boundary descriptor'
          if (i==0) cycle
          if (i<1.or.i>size(block%bdry_storage)) error stop 'RHS layout boundary descriptor invalid'
          boundaries(5,p)=block%bdry_storage(i)%local_start
       end if
    end do
    n=0
    if (allocated(block%ghost_storage)) n=size(block%ghost_storage)
    allocate(ghosts(8,n))
    do p=1,n
       associate(g=>block%ghost_storage(p))
       ghosts(:,p)=[g%source_domain,g%source_block,g%source_block_id,g%source_owner, &
            g%source_local_patch,g%local_start,g%n_node,g%source_level]
       end associate
    end do
    rebuilt=.true.
    if (plan%generation==generation.and.allocated(plan%storage)) then
       if (any(plan%identity/=identity).or.any(plan%nodes/=nodes)) &
            error stop 'RHS layout identity changed without new generation'
       if (.not.same(plan%patch,patches).or..not.same(plan%boundary,boundaries).or..not.same(plan%ghost,ghosts)) &
            error stop 'RHS layout topology changed without new generation'
       plan%reuses=plan%reuses+1_int64
       rebuilt=.false.
       return
    end if
    plan%identity=identity;plan%nodes=nodes;plan%generation=generation
    call move_alloc(patches,plan%patch)
    call move_alloc(boundaries,plan%boundary)
    call move_alloc(ghosts,plan%ghost)
    if (allocated(plan%storage)) deallocate(plan%storage,plan%node)
    n=size(block%patch)*PATCH_SIZE**2
    allocate(plan%storage(7,n),plan%node(7,n))
    do p=1,size(block%patch)
       do j=0,PATCH_SIZE-1
          do i=0,PATCH_SIZE-1
             n=(p-1)*PATCH_SIZE**2+j*PATCH_SIZE+i+1
             do t=1,7
                call resolve(block,p,i+dx(t),j+dy(t),plan%storage(t,n),plan%node(t,n))
                if (plan%node(t,n)<1.or.plan%node(t,n)>nodes(plan%storage(t,n))) &
                     error stop 'RHS layout node outside resident storage'
             end do
          end do
       end do
    end do
    plan%rebuilds=plan%rebuilds+1_int64
  end subroutine

  logical function same(a,b)
    integer, allocatable, intent(in) :: a(:,:),b(:,:)
    same=.false.
    if (.not.allocated(a).or..not.allocated(b)) return
    if (any(shape(a)/=shape(b))) return
    same=all(a==b)
  end function

  subroutine validate_fields(block,nodes)
    type(Block_Data), intent(in) :: block
    integer, intent(in) :: nodes(3)
    integer :: ns,nv
    ns=block%n_scalar_variable*block%n_field_level;nv=EDGE*block%n_field_level
    if (block%scalar_variable>S_MASS.or.block%scalar_variable+block%n_scalar_variable<=S_TEMP.or. &
         block%n_field_level<1) error stop 'RHS layout scalar family invalid'
    if (.not.allocated(block%scalar).or..not.allocated(block%scalar_mean).or..not.allocated(block%vector)) &
         error stop 'RHS layout missing owned fields'
    if (size(block%scalar)/=ns*nodes(1).or.size(block%scalar_mean)/=ns*nodes(1).or. &
         size(block%vector)/=nv*nodes(1)) error stop 'RHS layout owned field extent invalid'
    if (nodes(2)>0) then
       if (.not.allocated(block%bdry_scalar).or..not.allocated(block%bdry_scalar_mean).or. &
            .not.allocated(block%bdry_vector)) error stop 'RHS layout missing boundary fields'
       if (size(block%bdry_scalar)/=ns*nodes(2).or.size(block%bdry_scalar_mean)/=ns*nodes(2).or. &
            size(block%bdry_vector)/=nv*nodes(2)) error stop 'RHS layout boundary field extent invalid'
    end if
    if (nodes(3)>0) then
       if (.not.allocated(block%ghost_scalar).or..not.allocated(block%ghost_scalar_mean).or. &
            .not.allocated(block%ghost_vector)) error stop 'RHS layout missing ghost fields'
       if (size(block%ghost_scalar)/=ns*nodes(3).or.size(block%ghost_scalar_mean)/=ns*nodes(3).or. &
            size(block%ghost_vector)/=nv*nodes(3)) error stop 'RHS layout ghost field extent invalid'
    end if
  end subroutine

  subroutine resolve(block,p,i,j,storage,node)
    type(Block_Data), intent(in) :: block
    integer, intent(in) :: p,i,j
    integer, intent(out) :: storage,node
    integer :: side,neighbor,link,record,dims(2),first,x,y
    if (i>=0.and.i<PATCH_SIZE.and.j>=0.and.j<PATCH_SIZE) then
       storage=STORE_PATCH;node=block%patch(p)%elts_start+j*PATCH_SIZE+i+1
       return
    end if
    if (i<0) then
       side=WEST
       if (j<0) side=SOUTHWEST
       if (j>=PATCH_SIZE) side=NORTHWEST
    else if (i>=PATCH_SIZE) then
       side=EAST
       if (j<0) side=SOUTHEAST
       if (j>=PATCH_SIZE) side=NORTHEAST
    else if (j<0) then
       side=SOUTH
    else
       side=NORTH
    end if
    neighbor=block%patch(p)%neigh(side)
    if (neighbor>0) then
       if (neighbor+1>size(block%patch)) error stop 'RHS layout owned neighbor invalid'
       storage=STORE_PATCH;first=block%patch(neighbor+1)%elts_start;dims=PATCH_SIZE
    else if (neighbor<0) then
       link=-neighbor
       if (.not.allocated(block%block_bdry)) error stop 'RHS layout missing boundary link'
       if (link>size(block%block_bdry)) error stop 'RHS layout boundary link invalid'
       record=block%block_bdry(link)%ghost_id
       if (record>0) then
          storage=STORE_GHOST;first=block%ghost_storage(record)%local_start;dims=PATCH_SIZE
       else
          record=block%block_bdry(link)%storage_id
          if (record<1.or..not.allocated(block%bdry_storage)) error stop 'RHS layout missing boundary storage'
          if (record>size(block%bdry_storage)) error stop 'RHS layout boundary storage invalid'
          storage=STORE_BDRY;first=block%bdry_storage(record)%local_start;dims=block%block_bdry(link)%dims
       end if
    else
       error stop 'RHS layout missing neighbor'
    end if
    x=i;y=j
    if (i<0) x=i+dims(1)
    if (i>=PATCH_SIZE) x=i-PATCH_SIZE
    if (j<0) y=j+dims(2)
    if (j>=PATCH_SIZE) y=j-PATCH_SIZE
    if (x<0.or.x>=dims(1).or.y<0.or.y>=dims(2)) error stop 'RHS layout neighbor coordinate invalid'
    node=first+y*dims(1)+x+1
  end subroutine

  subroutine rhs_input_values(block,plan,generation,patch,node,layer,values,classes)
    ! Development consumer for exact comparison against accepted operands.
    ! Final production loops bind patch arrays directly; this gather is not a performance claim.
    type(Block_Data), intent(in) :: block
    type(Rhs_Input_Layout), intent(in) :: plan
    integer(int64), intent(in) :: generation
    integer, intent(in) :: patch,node,layer
    real(dp), intent(out) :: values(16)
    integer, intent(out) :: classes(16)
    integer :: n,t,s,a,k,offset,extent
    if (plan%generation/=generation.or..not.allocated(plan%storage)) error stop 'RHS input stale layout'
    if (patch<1.or.patch>size(block%patch).or.node<1.or.node>PATCH_SIZE**2) error stop 'RHS input target invalid'
    k=layer-block%field_level
    if (k<0.or.k>=block%n_field_level) error stop 'RHS input layer invalid'
    n=(patch-1)*PATCH_SIZE**2+node
    do t=1,6
       s=plan%storage(t,n);a=plan%node(t,n);extent=plan%nodes(s)
       offset=((S_MASS-block%scalar_variable)*block%n_field_level+k)*extent+a
       classes(t)=s
       select case(s)
       case(STORE_PATCH)
          values(t)=block%scalar(offset)+block%scalar_mean(offset)
       case(STORE_BDRY)
          values(t)=block%bdry_scalar(offset)+block%bdry_scalar_mean(offset)
       case(STORE_GHOST)
          values(t)=block%ghost_scalar(offset)+block%ghost_scalar_mean(offset)
       end select
    end do
    do t=1,10
       s=plan%storage(vnode(t),n);a=plan%node(vnode(t),n);extent=plan%nodes(s)
       offset=k*EDGE*extent+EDGE*(a-1)+vedge(t)
       classes(6+t)=s
       select case(s)
       case(STORE_PATCH)
          values(6+t)=block%vector(offset)
       case(STORE_BDRY)
          values(6+t)=block%bdry_vector(offset)
       case(STORE_GHOST)
          values(6+t)=block%ghost_vector(offset)
       end select
    end do
  end subroutine

  integer function rhs_layout_extent(plan,generation) result(n)
    type(Rhs_Input_Layout), intent(in) :: plan
    integer(int64), intent(in) :: generation
    if (plan%generation/=generation.or..not.allocated(plan%storage)) error stop 'RHS slot stale layout'
    n=sum(plan%nodes)
  end function

  integer function rhs_layout_slot(block,plan,generation,p,i,j) result(slot)
    ! Topology-time mapping to one flat derived-field workspace. Kernels use
    ! these integer offsets directly; there is no per-value storage dispatch.
    type(Block_Data), intent(in) :: block
    type(Rhs_Input_Layout), intent(in) :: plan
    integer(int64), intent(in) :: generation
    integer, intent(in) :: p,i,j
    integer :: storage,node
    if (plan%generation/=generation.or..not.allocated(plan%storage)) error stop 'RHS slot stale layout'
    if (p<1.or.p>size(block%patch)) error stop 'RHS slot patch invalid'
    call resolve(block,p,i,j,storage,node)
    if (node<1.or.node>plan%nodes(storage)) error stop 'RHS slot outside resident storage'
    slot=sum(plan%nodes(:storage-1))+node
  end function

  integer(int64) function rhs_input_layout_bytes(plan) result(bytes)
    type(Rhs_Input_Layout), intent(in) :: plan
    bytes=0_int64
    if (allocated(plan%storage)) bytes=int(size(plan%storage),int64)+size(plan%node,kind=int64)+ &
         size(plan%patch,kind=int64)+size(plan%boundary,kind=int64)+size(plan%ghost,kind=int64)
    bytes=bytes*int(storage_size(0)/8,int64)
  end function
end module
