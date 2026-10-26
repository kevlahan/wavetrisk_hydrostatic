module parallel_block_transform_placement_mod
  ! Preparation-only placement of fixed kernels. Physical storage identities
  ! are preserved across the cut; runtime execution uses compact integer slots.
  use iso_fortran_env,only : int8,int64
  use mpi_f08
  use kind_mod,only : dp
  use arch_mod,only : comm,rank,n_process,glo_id,loc_id
  use domain_mod,only : grid
  use shared_mod,only : EDGE,level_start,level_end
  use parallel_block_transform_kernels_mod
  use parallel_block_rhs_placement_mod,only : compile_rhs_owner_map
  implicit none
  private
  public :: Transform_Placement,distribute_transform_geometry
  type :: Transform_Placement
     type(Transform_Geometry) :: program
     integer,allocatable :: key(:,:),owner(:),state_key(:,:),level(:),mask(:,:)
     integer,allocatable :: scalar_order(:),inner_order(:),outer_order(:)
     logical,allocatable :: owned(:)
     integer(int64) :: sent_bytes=0,received_bytes=0
  end type
  type :: Node_Record
     integer :: key(4)=0,owner=-1,state_key(4)=0,level=0,mask(4)=0
  end type
  type :: Domain_Identity
     integer,allocatable :: owner(:),key(:,:)
  end type
contains
  subroutine prefix(count,displacement)
    integer,intent(in) :: count(:)
    integer,intent(out) :: displacement(:)
    integer :: r
    displacement=0
    do r=2,size(count)
       displacement(r)=displacement(r-1)+count(r-1)
    end do
  end subroutine

  subroutine distribute_transform_geometry(input,result,domain_extent)
    integer,intent(in),optional :: domain_extent(:)
    type(Transform_Geometry),intent(in) :: input
    type(Transform_Placement),intent(out) :: result
    type(Domain_Identity),allocatable :: identity(:)
    integer,allocatable :: scalar_owner(:),inner_owner(:),outer_owner(:),sc(:),sd(:),rc(:),rd(:),cursor(:)
    integer(int8),allocatable :: sent(:),received(:)
    type(Node_Record) :: node
    type(Transform_Scalar) :: sop
    type(Transform_Inner) :: iop
    type(Transform_Outer) :: oop
    integer :: d,i,j,pass,p,tag,width,ierr,capacity,used,b
    integer,allocatable :: table(:),keys(:,:),domain_base(:),touched(:)
    integer :: original_capacity,key_domain
    character(8) :: option
    integer :: counts(4),sizes(4)
    sizes=[storage_size(node),storage_size(sop),storage_size(iop),storage_size(oop)]/8
    allocate(identity(size(grid)))
    do d=1,size(grid)
       call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),identity(d)%owner,identity(d)%key)
    end do
    allocate(scalar_owner(size(input%scalar)),inner_owner(size(input%inner)),outer_owner(size(input%outer)))
    do i=1,size(input%scalar)
       d=input%scalar(i)%d
       scalar_owner(i)=identity(d)%owner(input%scalar(i)%target)
    end do
    do i=1,size(input%inner)
       d=input%inner(i)%d
       inner_owner(i)=identity(d)%owner((input%inner(i)%target(1)-1)/EDGE+1)
    end do
    ! The outer stencil is in-place. Its true intra-phase dependency components
    ! stay together; a runtime per-edge/dependency-depth interpreter is forbidden.
    call place_outer_components(input%outer,identity,outer_owner)
    width=storage_size(i)/8
    allocate(sc(n_process),sd(n_process),rc(n_process),rd(n_process),cursor(n_process))
    sc=0
    do pass=1,2
       if(pass==2) then
          call prefix(sc,sd);cursor=sd
          allocate(sent(sum(sc)))
       end if
       do d=1,size(grid)
          do i=1,size(identity(d)%owner)
             node%key=[glo_id(rank+1,d)+1,i,1,1]
             node%owner=identity(d)%owner(i);node%state_key=identity(d)%key(:,i)
             node%level=grid(d)%level%elts(i)
             node%mask=[grid(d)%mask_n%elts(i),grid(d)%mask_e%elts(EDGE*(i-1)+1:EDGE*i)]
             call emit(1,node%owner,transfer(node,[0_int8]))
          end do
       end do
       do i=1,size(input%scalar)
          sop=input%scalar(i);sop%d=i
          call emit(2,scalar_owner(i),transfer(sop,[0_int8]))
       end do
       do i=1,size(input%inner)
          iop=input%inner(i);iop%d=i
          call emit(3,inner_owner(i),transfer(iop,[0_int8]))
       end do
       do i=1,size(input%outer)
          oop=input%outer(i)
          call emit(4,outer_owner(i),transfer(oop,[0_int8]),i)
       end do
    end do
    if(any(cursor/=sd+sc)) error stop 'transform metadata send extent'
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform metadata counts'
    call prefix(rc,rd);allocate(received(sum(rc)))
    call MPI_Alltoallv(sent,sc,sd,MPI_BYTE,received,rc,rd,MPI_BYTE,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform metadata exchange'
    result%sent_bytes=sum(sc)-sc(rank+1);result%received_bytes=sum(rc)-rc(rank+1)
    deallocate(sent,identity)
    counts=0;p=1
    do while(p<=size(received))
       tag=transfer(received(p:p+width-1),tag);p=p+2*width
       if(tag<1.or.tag>4) error stop 'transform metadata tag'
       counts(tag)=counts(tag)+1;p=p+sizes(tag)
    end do
    if(p/=size(received)+1) error stop 'transform metadata frame extent'
    allocate(result%program%scalar(counts(2)),result%program%inner(counts(3)),result%program%outer(counts(4)), &
         result%scalar_order(counts(2)),result%inner_order(counts(3)),result%outer_order(counts(4)))
    ! All local canonical nodes plus the operation footprint. Open addressing
    ! compacts physical keys once; no owner/key lookup occurs during arithmetic.
    original_capacity=counts(1)+13*counts(2)+38*counts(3)+16*counts(4)
    capacity=original_capacity
    if(present(domain_extent)) then
       allocate(domain_base(size(domain_extent)),touched(size(domain_extent)));touched=0;domain_base=-1
       p=1
       do while(p<=size(received))
          tag=transfer(received(p:p+width-1),tag);p=p+2*width
          select case(tag)
          case(1)
             node=transfer(received(p:p+sizes(tag)-1),node);key_domain=node%key(1)
          case(2)
             sop=transfer(received(p:p+sizes(tag)-1),sop);key_domain=sop%global_domain
          case(3)
             iop=transfer(received(p:p+sizes(tag)-1),iop);key_domain=iop%global_domain
          case(4)
             oop=transfer(received(p:p+sizes(tag)-1),oop);key_domain=oop%global_domain
          case default
             error stop 'transform footprint record tag'
          end select
          if(key_domain<1.or.key_domain>size(touched)) error stop 'transform footprint namespace'
          touched(key_domain)=1;p=p+sizes(tag)
       end do
       capacity=0
       do d=1,size(touched)
          if(touched(d)==0) cycle
          domain_base(d)=capacity;capacity=capacity+domain_extent(d)
       end do
       allocate(table(max(1,capacity)));table=0;used=0
       ! First count the unique physical slots; then allocate keys exactly.
       p=1
       do while(p<=size(received))
          tag=transfer(received(p:p+width-1),tag);p=p+2*width
          select case(tag)
          case(1)
             node=transfer(received(p:p+sizes(tag)-1),node)
             call mark(node%key(1),node%key(2))
          case(2)
             sop=transfer(received(p:p+sizes(tag)-1),sop)
             call mark(sop%global_domain,sop%target);call mark(sop%global_domain,sop%base)
             do j=1,sop%n
                call mark(sop%global_domain,sop%source(j))
             end do
          case(3)
             iop=transfer(received(p:p+sizes(tag)-1),iop)
             do j=1,32
                if(iop%source(j)>0) call mark(iop%global_domain,(iop%source(j)-1)/EDGE+1)
             end do
             do j=1,6
                call mark(iop%global_domain,(iop%target(j)-1)/EDGE+1)
             end do
          case(4)
             oop=transfer(received(p:p+sizes(tag)-1),oop)
             do j=1,13
                if(oop%source(j)>0) call mark(oop%global_domain,(oop%source(j)-1)/EDGE+1)
             end do
             do j=1,3
                if(oop%target(j)>0) call mark(oop%global_domain,(oop%target(j)-1)/EDGE+1)
             end do
          end select
          p=p+sizes(tag)
       end do
       allocate(keys(4,used));table=0
    else
       allocate(keys(4,max(1,capacity)),table(max(1,2*capacity+1)))
    end if
    table=0;used=0;counts=0;p=1
    do while(p<=size(received))
       tag=transfer(received(p:p+width-1),tag);p=p+width
       b=transfer(received(p:p+width-1),b);p=p+width
       counts(tag)=counts(tag)+1;i=counts(tag)
       select case(tag)
       case(1)
          node=transfer(received(p:p+sizes(tag)-1),node)
          j=slot(node%key(1),node%key(2))
       case(2)
          sop=transfer(received(p:p+sizes(tag)-1),sop)
          result%scalar_order(i)=sop%d;sop%d=0
          sop%target=slot(sop%global_domain,sop%target);sop%base=slot(sop%global_domain,sop%base)
          do j=1,sop%n
             sop%source(j)=slot(sop%global_domain,sop%source(j))
          end do
          result%program%scalar(i)=sop
       case(3)
          iop=transfer(received(p:p+sizes(tag)-1),iop)
          result%inner_order(i)=iop%d;iop%d=0
          do j=1,32
             if(iop%source(j)>0) iop%source(j)=edge_slot(iop%global_domain,iop%source(j))
          end do
          do j=1,6
             iop%target(j)=edge_slot(iop%global_domain,iop%target(j))
          end do
          result%program%inner(i)=iop
       case(4)
          oop=transfer(received(p:p+sizes(tag)-1),oop);result%outer_order(i)=b
          do j=1,13
             if(oop%source(j)>0) oop%source(j)=edge_slot(oop%global_domain,oop%source(j))
          end do
          do j=1,3
             if(oop%target(j)>0) oop%target(j)=edge_slot(oop%global_domain,oop%target(j))
          end do
          result%program%outer(i)=oop
       end select
       p=p+sizes(tag)
    end do
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
    if(trim(option)=='1') write(*,*) 'Structure284 footprint:',rank,original_capacity,used, &
         4_int64*(4_int64*original_capacity+2_int64*original_capacity+1), &
         4_int64*(size(keys,kind=int64)+size(table,kind=int64))
    result%key=keys(:,:used)
    allocate(result%owned(used),result%owner(used),result%state_key(4,used),result%level(used),result%mask(4,used))
    result%owned=.false.;result%owner=-1;result%state_key=0;result%level=-1;result%mask=0
    p=1
    do while(p<=size(received))
       tag=transfer(received(p:p+width-1),tag);p=p+2*width
       if(tag==1) then
          node=transfer(received(p:p+sizes(tag)-1),node);j=slot(node%key(1),node%key(2))
          if(result%owned(j)) error stop 'transform duplicate owned node'
          result%owned(j)=.true.;result%owner(j)=node%owner;result%state_key(:,j)=node%state_key
          result%level(j)=node%level;result%mask(:,j)=node%mask
       end if
       p=p+sizes(tag)
    end do
  contains
    subroutine mark(domain,node)
      integer,intent(in) :: domain,node
      integer :: h
      if(node<1.or.node>domain_extent(domain)) error stop 'transform footprint node'
      h=domain_base(domain)+node
      if(table(h)==0) then
         used=used+1;table(h)=used
      end if
    end subroutine
    subroutine emit(kind,destination,data,ordinal)
      integer,intent(in) :: kind,destination
      integer(int8),intent(in) :: data(:)
      integer,optional,intent(in) :: ordinal
      integer :: a,t,z
      a=destination+1;t=0
      if(present(ordinal)) t=ordinal
      if(pass==1) then
         sc(a)=sc(a)+2*width+size(data)
      else
         z=cursor(a)+1
         sent(z:z+width-1)=transfer(kind,[0_int8]);z=z+width
         sent(z:z+width-1)=transfer(t,[0_int8]);z=z+width
         sent(z:z+size(data)-1)=data;cursor(a)=z+size(data)-1
      end if
    end subroutine
    integer function slot(domain,node) result(index)
      integer,intent(in) :: domain,node
      integer :: h
      integer(int64) :: bits
      if(domain<1.or.node<1) error stop 'transform physical key'
      if(present(domain_extent)) then
         if(domain>size(domain_extent)) error stop 'transform physical namespace'
         if(node>domain_extent(domain).or.domain_base(domain)<0) error stop 'transform physical slot'
         h=domain_base(domain)+node;index=table(h)
         if(index==0) then
            used=used+1
            if(used>size(keys,2)) error stop 'transform exact footprint capacity'
            index=used;table(h)=index;keys(:,index)=[domain,node,1,1]
         end if
         return
      end if
      bits=modulo(int(node,int64)*104729_int64+int(domain,int64)*65599_int64,2147483647_int64)
      bits=ieor(bits,shiftr(bits,16));h=int(modulo(bits,int(size(table),int64)))+1
      do
         index=table(h)
         if(index==0) then
            used=used+1
            if(used>size(keys,2)) error stop 'transform footprint capacity'
            index=used;table(h)=index;keys(:,index)=[domain,node,1,1];return
         end if
         if(keys(1,index)==domain.and.keys(2,index)==node) return
         h=modulo(h,size(table))+1
      end do
    end function
    integer function edge_slot(domain,edge_code) result(index)
      integer,intent(in) :: domain,edge_code
      index=EDGE*(slot(domain,(edge_code-1)/EDGE+1)-1)+mod(edge_code-1,EDGE)+1
    end function
  end subroutine

  subroutine place_outer_components(operations,identity,destination)
    type(Transform_Outer),intent(in) :: operations(:)
    type(Domain_Identity),intent(in) :: identity(:)
    integer,intent(out) :: destination(:)
    integer,allocatable :: parent(:),last_writer(:),first_reader(:),next_reader(:),reader_operation(:),choice(:)
    integer :: d,l,i,j,a,n,nread,root,home
    allocate(parent(size(operations)),choice(size(operations)))
    parent=[(i,i=1,size(operations))];choice=-1;destination=-1
    do d=1,size(grid)
       n=EDGE*grid(d)%node%length
       allocate(last_writer(n),first_reader(n))
       allocate(next_reader(16*size(operations)),reader_operation(16*size(operations)))
       do l=level_start-1,level_end-1
          last_writer=0;first_reader=0;nread=0
          do i=1,size(operations)
             associate(op=>operations(i))
             if(op%global_domain/=glo_id(rank+1,d)+1.or.op%level/=l) cycle
             do j=1,13
                a=op%source(j)
                if(a==0) cycle
                call read_address(a,i)
             end do
             if(op%kind/=0) then
                call read_address(op%target(2),i);call read_address(op%target(3),i)
             end if
             ! Include both writes, ordered as in the fixed kernel.
             call write_address(op%target(3),i);call write_address(op%target(2),i)
             end associate
          end do
       end do
       deallocate(last_writer,first_reader,next_reader,reader_operation)
    end do
    ! Choose a real canonical target owner where possible. Components with
    ! only alias targets retain the deterministic numerical scaffold owner.
    do i=1,size(operations)
       d=loc_id(operations(i)%global_domain)+1;root=find(i)
       do j=2,3
          a=(operations(i)%target(j)-1)/EDGE+1
          home=identity(d)%owner(a)
          if(choice(root)<0.or.identity(d)%key(1,a)>0) choice(root)=home
          if(identity(d)%key(1,a)>0) exit
       end do
    end do
    do i=1,size(operations)
       destination(i)=choice(find(i))
       if(destination(i)<0.or.destination(i)>=n_process) error stop 'transform component placement'
    end do
  contains
    integer function find(a) result(r)
      integer,intent(in) :: a
      r=a
      do while(parent(r)/=r)
         parent(r)=parent(parent(r));r=parent(r)
      end do
    end function
    subroutine join(a,b)
      integer,intent(in) :: a,b
      integer :: x,y
      if(a==0.or.b==0) return
      x=find(a);y=find(b)
      if(x/=y) parent(y)=x
    end subroutine
    subroutine read_address(a,op)
      integer,intent(in) :: a,op
      call join(op,last_writer(a))
      nread=nread+1
      if(nread>size(reader_operation)) error stop 'transform dependency capacity'
      reader_operation(nread)=op;next_reader(nread)=first_reader(a);first_reader(a)=nread
    end subroutine
    subroutine write_address(a,op)
      integer,intent(in) :: a,op
      integer :: r
      call join(op,last_writer(a));r=first_reader(a)
      do while(r>0)
         call join(op,reader_operation(r));r=next_reader(r)
      end do
      first_reader(a)=0;last_writer(a)=op
    end subroutine
  end subroutine
end module
