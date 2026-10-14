module parallel_block_rhs_schedule_mod
  ! Compile the actual step1 traversal, including positive physical boundaries
  ! and immediate southwest writes. The caller supplies accepted patch order
  ! and pentagon identities; block DFS order and geometric corners are NOT an
  ! alternative source of numerical ordering. No Domain or MPI dependency.
  use iso_fortran_env, only : int64
  use shared_mod, only : NORTH,EAST,SOUTH,WEST,SOUTHEAST,IJPLUS,IJMINUS,IPLUSJMINUS,IMINUSJPLUS
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data
  use parallel_block_rhs_layout_mod, only : Rhs_Input_Layout,rhs_layout_slot,rhs_layout_extent
  implicit none
  private
  public :: Rhs_Producer_Schedule,compile_rhs_producer_schedule,rhs_schedule_last_writes,rhs_schedule_bytes
  public :: rhs_schedule_definitions
  integer, parameter :: dx(7)=[0,0,1,0,-1,1,-1],dy(7)=[0,1,0,-1,0,1,-1]
  type :: Rhs_Producer_Schedule
     integer(int64) :: generation=-1_int64
     integer :: block_id=-1,level=-1,nslot=0
     ! C,N,E,S,W,NE,SW; integers address a flat owned/boundary/ghost workspace.
     integer, allocatable :: node(:,:),patch(:),xy(:,:)
     logical, allocatable :: southwest(:)
     integer, allocatable :: pentagon_node(:,:),pentagon_patch(:),corner(:)
  end type
contains
  subroutine compile_rhs_producer_schedule(block,layout,generation,level,patch_order, &
       pentagon_patch,corner,program)
    type(Block_Data), intent(in) :: block
    type(Rhs_Input_Layout), intent(in) :: layout
    integer(int64), intent(in) :: generation
    integer, intent(in) :: level,patch_order(:),pentagon_patch(:),corner(:)
    type(Rhs_Producer_Schedule), intent(out) :: program
    integer :: pass,n,s,p,i,j,t,c
    logical :: boundary(4)
    logical, allocatable :: seen(:)
    program%nslot=rhs_layout_extent(layout,generation)
    program%generation=generation;program%block_id=block%id;program%level=level
    allocate(seen(size(block%patch)));seen=.false.
    do s=1,size(patch_order)
       p=patch_order(s)
       if (p<1.or.p>size(block%patch)) error stop 'RHS schedule patch invalid'
       if (seen(p)) error stop 'RHS schedule duplicate patch'
       if (block%patch(p)%level/=level) error stop 'RHS schedule level mismatch'
       seen(p)=.true.
    end do
    if (any(seen.neqv.(block%patch%level==level))) error stop 'RHS schedule incomplete level'
    if (size(corner)/=size(pentagon_patch)) error stop 'RHS schedule pentagon extent mismatch'
    do s=1,size(corner)
       p=pentagon_patch(s);c=corner(s)
       if (p<1.or.p>size(block%patch)) error stop 'RHS schedule pentagon patch invalid'
       if (.not.seen(p)) error stop 'RHS schedule pentagon outside level'
       if (c<IJPLUS.or.c>IMINUSJPLUS) error stop 'RHS schedule corner invalid'
    end do
    do s=2,size(corner)
       if (corner(s)<=corner(s-1)) error stop 'RHS schedule pentagon order invalid'
    end do
    ! Each direct call retains its immediately following SW call. Physical
    ! Domain boundaries are identified by the migrated boundary descriptor;
    ! cross-block ghost boundaries must not manufacture extra producers.
    do pass=1,2
       n=0
       do s=1,size(patch_order)
          p=patch_order(s)
          do t=1,4
             boundary(t)=physical_boundary(block,p,t)
          end do
          do j=0,PATCH_SIZE-1
             do i=0,PATCH_SIZE-1
                call event(i,j,(i==0.and.boundary(WEST)).or.(j==0.and.boundary(SOUTH)))
             end do
          end do
          if (boundary(NORTH)) then
             do i=0,PATCH_SIZE-1
                call event(i,PATCH_SIZE,i==0.and.boundary(WEST))
             end do
             call event(PATCH_SIZE,PATCH_SIZE,.false.)
          end if
          if (boundary(EAST)) then
             do j=0,PATCH_SIZE-1
                call event(PATCH_SIZE,j,j==0.and.boundary(SOUTH))
             end do
          end if
       end do
       if (pass==1) allocate(program%node(7,n),program%patch(n),program%xy(2,n),program%southwest(n))
    end do
    program%pentagon_patch=pentagon_patch;program%corner=corner
    allocate(program%pentagon_node(7,size(corner)))
    do s=1,size(corner)
       p=pentagon_patch(s);i=0;j=0
       if (corner(s)==IJPLUS.or.corner(s)==IPLUSJMINUS) i=PATCH_SIZE
       if (corner(s)==IJPLUS.or.corner(s)==IMINUSJPLUS) j=PATCH_SIZE
       do t=1,7
          ! Only addresses used by this corner are necessary. In particular
          ! IJPLUS does not read the far NE diagonal; do not demand that halo.
          if ((corner(s)==IJPLUS.and.t>=6).or.(corner(s)==IJMINUS.and.t==6).or. &
               (corner(s)==IPLUSJMINUS.and.t==2).or.(corner(s)==IMINUSJPLUS.and.t==3)) then
             program%pentagon_node(t,s)=rhs_layout_slot(block,layout,generation,p,i,j)
          else
             program%pentagon_node(t,s)=rhs_layout_slot(block,layout,generation,p,i+dx(t),j+dy(t))
          end if
       end do
    end do
  contains
    subroutine event(x,y,sw)
      integer, intent(in) :: x,y
      logical, intent(in) :: sw
      integer :: q,south_dims(2),southeast_dims(2)
      n=n+1
      if (pass==1) return
      program%patch(n)=p;program%xy(:,n)=[x,y];program%southwest(n)=sw
      do q=1,7
         if (x==PATCH_SIZE.and.y==0.and.q==4) then
            ! comp_offs3 intentionally uses SOUTH's height in the SE offset.
            ! With a 4x3 SOUTH strip and 3x4 SE strip this is row 2, not
            ! the last SE row selected by geometric idx(PATCH_SIZE,-1).
            ! Preserve the accepted operator's address, including this alias.
            south_dims=neighbor_dims(block,p,SOUTH)
            southeast_dims=neighbor_dims(block,p,SOUTHEAST)
            program%node(q,n)=rhs_layout_slot(block,layout,generation,p,x, &
                 south_dims(2)-southeast_dims(2)-1)
         else
            program%node(q,n)=rhs_layout_slot(block,layout,generation,p,x+dx(q),y+dy(q))
         end if
      end do
    end subroutine
  end subroutine

  function neighbor_dims(block,p,side) result(dims)
    type(Block_Data), intent(in) :: block
    integer, intent(in) :: p,side
    integer :: dims(2),link
    dims=PATCH_SIZE
    link=-block%patch(p)%neigh(side)
    if (link>0) then
       if (block%block_bdry(link)%ghost_id<=0) dims=block%block_bdry(link)%dims
    end if
  end function

  logical function physical_boundary(block,p,side) result(physical)
    type(Block_Data), intent(in) :: block
    integer, intent(in) :: p,side
    integer :: link
    physical=.false.
    link=-block%patch(p)%neigh(side)
    if (link<=0) return
    if (.not.allocated(block%block_bdry)) error stop 'RHS schedule missing boundary'
    if (link>size(block%block_bdry)) error stop 'RHS schedule boundary invalid'
    physical=block%block_bdry(link)%bdry_side>0
  end function

  subroutine rhs_schedule_last_writes(program,flux_writer,pv_writer,vort_writer)
    ! Local write provenance for route compilation. Positive ordinal = direct
    ! event, negative ordinal = its SW phase; n_event+s denotes pentagon s.
    ! These are local aliases only: remote final-writer arbitration is required
    ! before a derived halo can be declared ready. The caller owns scratch.
    type(Rhs_Producer_Schedule), intent(in) :: program
    integer, intent(out) :: flux_writer(:,:),pv_writer(:,:),vort_writer(:,:)
    integer :: s,n(7),w
    if (any(shape(flux_writer)/=[3,program%nslot]).or.any(shape(pv_writer)/=[3,program%nslot]).or. &
         any(shape(vort_writer)/=[2,program%nslot])) &
         error stop 'RHS schedule writer extent invalid'
    flux_writer=0;vort_writer=0
    do s=1,size(program%patch)
       n=program%node(:,s)
       flux_writer(:,n(1))=s
       vort_writer(:,n(1))=s
       vort_writer(1,n(5))=s;vort_writer(2,n(4))=s
       if (program%southwest(s)) then
          flux_writer(1,n(5))=-s;flux_writer(2,n(7))=-s;flux_writer(3,n(4))=-s
          vort_writer(:,n(7))=-s
       end if
    end do
    pv_writer=flux_writer
    ! Flux is unaffected by post_step1; PV has different final writers.
    do s=1,size(program%corner)
       n=program%pentagon_node(:,s);w=size(program%patch)+s
       select case(program%corner(s))
       case(IJMINUS)
          vort_writer(:,n(7))=w
          pv_writer(1,n(5))=w;pv_writer(3,n(4))=w
       case(IPLUSJMINUS)
          vort_writer(1,n(7))=w;vort_writer(2,n(4))=w
          pv_writer(1,n(1))=w;pv_writer(2,n(7))=w
       case(IMINUSJPLUS)
          vort_writer(2,n(7))=w;vort_writer(1,n(5))=w
          pv_writer(3,n(1))=w;pv_writer(2,n(7))=w
       case(IJPLUS)
          vort_writer(:,n(1))=w
          pv_writer(1,n(1))=w;pv_writer(3,n(1))=w
       end select
    end do
  end subroutine

  subroutine rhs_schedule_definitions(program,patch_ordinal,last_patch,physical_key,field,workspace_base, &
       definition_key,definition_order,definition_slot)
    ! Bridge the ordered producer schedule to compile_rhs_exchange. This
    ! function cannot infer physical alias identity or accepted global order:
    ! both must come from the topology compiler, without coordinate merging.
    ! field: 1=flux, 2=PV, 3=vorticity. Each uses a separate phase workspace.
    type(Rhs_Producer_Schedule), intent(in) :: program
    integer, intent(in) :: patch_ordinal(:),last_patch,physical_key(:,:),field,workspace_base
    integer, allocatable, intent(out) :: definition_key(:,:),definition_order(:,:),definition_slot(:)
    integer, allocatable :: flux(:,:),pv(:,:),vort(:,:),writer(:,:),event_order(:,:)
    integer :: i,n,c,w,p,previous,within,q,m
    if (field<1.or.field>3.or.workspace_base<0.or.last_patch<1.or.last_patch==huge(0)) &
         error stop 'RHS writer field or extent invalid'
    if (any(shape(physical_key)/=[3,program%nslot]).or.any(physical_key<=0)) &
         error stop 'RHS writer physical identity invalid'
    allocate(event_order(2,size(program%patch)))
    previous=0;within=0
    do i=1,size(program%patch)
       p=program%patch(i)
       if (p<1.or.p>size(patch_ordinal)) error stop 'RHS writer patch identity invalid'
       q=patch_ordinal(p)
       if (q<1.or.q>last_patch.or.q<previous) error stop 'RHS writer global patch order invalid'
       if (q/=previous) within=0
       within=within+1;event_order(:,i)=[q,2*within-1];previous=q
    end do
    allocate(flux(3,program%nslot),pv(3,program%nslot),vort(2,program%nslot))
    call rhs_schedule_last_writes(program,flux,pv,vort)
    select case(field)
    case(1)
       call move_alloc(flux,writer)
    case(2)
       call move_alloc(pv,writer)
    case(3)
       call move_alloc(vort,writer)
    end select
    m=size(writer,1);n=count(writer/=0)
    allocate(definition_key(4,n),definition_order(2,n),definition_slot(n));i=0
    do p=1,program%nslot
       do c=1,m
          w=writer(c,p)
          if (w==0) cycle
          i=i+1;definition_key(:,i)=[physical_key(:,p),c]
          definition_slot(i)=workspace_base+m*(p-1)+c
          if (abs(w)<=size(program%patch)) then
             definition_order(:,i)=event_order(:,abs(w))
             if (w<0) definition_order(2,i)=definition_order(2,i)+1
          else
             ! post_step1 runs after ALL patches, including other owners'
             ! patches. A local maximum would incorrectly precede a peer.
             q=w-size(program%patch)
             definition_order(:,i)=[last_patch+1,program%corner(q)]
          end if
       end do
    end do
  end subroutine

  integer(int64) function rhs_schedule_bytes(program) result(bytes)
    type(Rhs_Producer_Schedule), intent(in) :: program
    bytes=0_int64
    if (.not.allocated(program%node)) return
    bytes=int(storage_size(0)/8,int64)*(size(program%node,kind=int64)+size(program%patch,kind=int64)+ &
         size(program%xy,kind=int64)+size(program%pentagon_node,kind=int64)+ &
         size(program%pentagon_patch,kind=int64)+size(program%corner,kind=int64)) + &
         int(storage_size(.false.)/8,int64)*size(program%southwest,kind=int64)
  end function
end module
