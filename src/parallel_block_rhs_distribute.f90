module parallel_block_rhs_distribute_mod
  ! Topology placement and metadata delivery. Numerical stages use precompiled
  ! peer routes; no global metadata collective is executed by this module then.
  use iso_fortran_env, only : int8,int64
  use kind_mod, only : dp
  use mpi_f08
  use arch_mod, only : comm,rank,n_process,glo_id
  use domain_mod, only : grid
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators,merge_rhs_domain_operators
  use parallel_block_rhs_placement_mod, only : compile_rhs_owner_map,select_rhs_placement
  use parallel_block_rhs_wire_mod, only : pack_rhs_placement,unpack_rhs_placement
  use parallel_block_rhs_exchange_mod, only : Rhs_Exchange,compile_rhs_exchange,exchange_rhs_values,release_rhs_exchange
  implicit none
  private
  public :: Rhs_Placement,distribute_rhs_geometry,deliver_rhs_metadata
  type :: Rhs_Placement
     type(Rhs_Domain_Operators) :: program
     integer,allocatable :: key(:,:),base(:)
     logical,allocatable :: owned(:)
     integer(int64) :: sent_metadata_bytes=0,received_metadata_bytes=0
  end type
  type :: Bytes
     integer(int8),allocatable :: value(:)
  end type
contains
  subroutine distribute_rhs_geometry(input,result)
    type(Rhs_Domain_Operators),intent(in) :: input(:)
    type(Rhs_Placement),intent(out) :: result
    type(Rhs_Domain_Operators),allocatable :: received(:)
    type :: Node_Identity
       integer,allocatable :: original(:)
       logical,allocatable :: owned(:)
    end type
    type(Node_Identity),allocatable :: identity(:)
    type(Bytes),allocatable :: packet_list(:)
    integer,allocatable :: destination(:),sc(:),sd(:),rc(:),rd(:),cursor(:)
    integer(int8),allocatable :: sent(:),incoming(:)
    integer :: r,d,n,p,width,ierr,i,j,used,total,first,last
    if(size(input)/=size(grid)) error stop 'RHS geometry source extent'
    width=storage_size(n)/8
    allocate(sc(n_process),sd(n_process),rc(n_process),rd(n_process),cursor(n_process))
    ! At most one packet per source Domain/destination. Retain each once,
    ! then pack directly into its final peer segment without growing copies.
    allocate(packet_list(size(input)*n_process),destination(size(input)*n_process));sc=0;used=0
    block
      type(Rhs_Domain_Operators) :: part
      integer,allocatable :: owners(:),state_key(:,:),original(:)
      logical,allocatable :: owned(:),active_rank(:)
      integer(int8),allocatable :: packet(:)
      allocate(active_rank(0:n_process-1))
      do d=1,size(input)
         call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),owners,state_key)
         active_rank=.false.
         do i=1,size(owners)
            if(owners(i)<0.or.owners(i)>=n_process) error stop 'RHS geometry owner range'
            active_rank(owners(i))=.true.
         end do
         do r=0,n_process-1
            if(.not.active_rank(r)) cycle
            call select_rhs_placement(input(d),owners,r,part,original,owned)
            part%domain_id=glo_id(rank+1,d)
            call pack_rhs_placement(part,original,owned,packet)
            used=used+1;destination(used)=r+1;sc(r+1)=sc(r+1)+width+size(packet)
            call move_alloc(packet,packet_list(used)%value)
         end do
      end do
    end block
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'RHS geometry metadata counts'
    sd=0;rd=0
    do r=2,n_process
       sd(r)=sd(r-1)+sc(r-1);rd(r)=rd(r-1)+rc(r-1)
    end do
    allocate(sent(sum(sc)),incoming(sum(rc)));cursor=sd
    do i=1,used
       r=destination(i);p=cursor(r)+1;n=size(packet_list(i)%value)
       sent(p:p+width-1)=transfer(n,[0_int8]);p=p+width
       sent(p:p+n-1)=packet_list(i)%value;cursor(r)=p+n-1
    end do
    if(any(cursor/=sd+sc)) error stop 'RHS geometry packed extent'
    deallocate(packet_list,destination)
    call MPI_Alltoallv(sent,sc,sd,MPI_BYTE,incoming,rc,rd,MPI_BYTE,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'RHS geometry metadata delivery'
    deallocate(sent)
    ! Count framed packets before unpacking. An append of a derived type with
    ! allocatable members deep-copied every previously received program.
    total=0
    do r=1,n_process
       p=rd(r)+1
       do while(p<=rd(r)+rc(r))
          if(rd(r)+rc(r)-p+1<width) error stop 'RHS geometry truncated header'
          n=transfer(incoming(p:p+width-1),n);p=p+width
          if(n<0.or.n>rd(r)+rc(r)-p+1) error stop 'RHS geometry truncated plan'
          p=p+n;total=total+1
       end do
    end do
    allocate(received(total),identity(total));used=0
    do r=1,n_process
       p=rd(r)+1
       do while(p<=rd(r)+rc(r))
          n=transfer(incoming(p:p+width-1),n);p=p+width;used=used+1
          call unpack_rhs_placement(incoming(p:p+n-1),received(used), &
               identity(used)%original,identity(used)%owned)
          p=p+n
          do i=1,used-1
             if(received(i)%domain_id==received(used)%domain_id) error stop 'RHS duplicate physical namespace'
          end do
          if(size(identity(used)%original)/=received(used)%nnode.or. &
               size(identity(used)%owned)/=received(used)%nnode) error stop 'RHS geometry identity extent'
       end do
    end do
    deallocate(incoming)
    total=sum(received%nnode);allocate(result%key(4,total),result%owned(total));last=0
    do i=1,size(received)
       first=last+1;last=last+received(i)%nnode
       do j=1,received(i)%nnode
          result%key(:,first+j-1)=[received(i)%domain_id+1,identity(i)%original(j),1,1]
       end do
       result%owned(first:last)=identity(i)%owned
    end do
    deallocate(identity)
    call merge_rhs_domain_operators(received,result%base,result%program,with_dynamics=.true.)
    if(size(result%owned)/=result%program%nnode) error stop 'RHS geometry merged extent'
    result%sent_metadata_bytes=int(sum(sc)-sc(rank+1),int64)
    result%received_metadata_bytes=int(sum(rc)-rc(rank+1),int64)
  end subroutine

  subroutine deliver_rhs_metadata(source_key,source_values,request_key,values,topology,phase)
    ! Topology-only keyed transfer. Integer metadata are exactly representable
    ! in binary64; callers validate round trips before converting back.
    integer,intent(in) :: source_key(:,:),request_key(:,:),phase
    real(dp),intent(in) :: source_values(:,:)
    real(dp),allocatable,intent(out) :: values(:,:)
    integer(int64),intent(in) :: topology
    type(Rhs_Exchange) :: route
    integer,allocatable :: order(:,:),source(:),target(:),sign(:)
    real(dp),allocatable :: work(:,:)
    integer :: ns,nq,nc,i
    ns=size(source_key,2);nq=size(request_key,2);nc=size(source_values,2)
    if(size(source_values,1)/=ns) error stop 'RHS metadata source extent'
    allocate(order(2,ns),source(ns),target(nq),sign(nq),work(ns+nq,nc))
    order=1;sign=1;source=[(i,i=1,ns)];target=[(ns+i,i=1,nq)]
    work=0.0_dp;work(:ns,:)=source_values
    call compile_rhs_exchange(route,comm,topology,0_int64,phase,ns+nq,nc,source_key,order,source, &
         request_key,target,sign)
    call exchange_rhs_values(route,topology,0_int64,phase,1_int64,work)
    allocate(values(nq,nc));values=work(ns+1:,:)
    call release_rhs_exchange(route)
  end subroutine
end module
