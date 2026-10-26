module parallel_block_transform_adapter_mod
  use audit316_mod, only: work325_enabled,work325_record
  ! Bounded raw-RK boundary version adapter. Register only actual readers;
  ! no full Domain numeric image, runtime lookup or global numeric collective.
  use iso_fortran_env,only : int64
  use kind_mod,only : dp
  use mpi_f08
  use arch_mod,only : comm,rank,n_process,owner,loc_id
  use domain_mod,only : Float_Field,grid
  use shared_mod,only : N_VARIABLE,zlevels,scalars,S_VELO
  implicit none
  private
  public :: Transform_Adapter,compile_transform_adapter,read_transform_adapter
  type :: Transform_Adapter
     integer :: component=0,tag=0,columns=0
     integer,allocatable :: sc(:),sd(:),rc(:),rd(:),source(:,:),index(:),target(:),sign(:)
     integer,allocatable :: send_peer(:),receive_peer(:)
     real(dp),allocatable :: sent(:,:),received(:,:)
     type(MPI_Request),allocatable :: request(:)
  end type
contains
  subroutine compile_transform_adapter(plan,component,tag,key,target,sign)
    type(Transform_Adapter),intent(out) :: plan
    integer,intent(in) :: component,tag,key(:,:),target(:),sign(:)
    integer,allocatable :: unique(:,:),table(:),map(:),position(:),cursor(:),outgoing(:,:),incoming(:,:)
    integer :: i,j,h,n,r,a,ierr
    integer(int64) :: bits
    if(component/=1.and.component/=3) error stop 'transform adapter component'
    if(size(key,1)/=4.or.size(key,2)/=size(target).or.size(sign)/=size(target)) error stop 'transform adapter extent'
    plan%component=component;plan%tag=tag;plan%columns=merge(2*zlevels,zlevels,component==1)
    allocate(unique(4,size(target)),table(max(1,2*size(target)+1)),map(size(target)))
    table=0;n=0
    do i=1,size(target)
       bits=0
       do j=1,4
          bits=modulo(bits*65599_int64+key(j,i),2147483647_int64)
       end do
       h=int(modulo(bits,int(size(table),int64)))+1
       do
          a=table(h)
          if(a==0) then
             n=n+1;a=n;table(h)=a;unique(:,a)=key(:,i);exit
          end if
          if(all(unique(:,a)==key(:,i))) exit
          h=modulo(h,size(table))+1
       end do
       map(i)=a
    end do
    allocate(plan%sc(n_process),plan%sd(n_process),plan%rc(n_process),plan%rd(n_process), &
         cursor(n_process),position(n),outgoing(4,n))
    plan%sc=0
    do i=1,n
       r=owner(unique(1,i))+1;plan%sc(r)=plan%sc(r)+1
    end do
    call MPI_Alltoall(plan%sc,1,MPI_INTEGER,plan%rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform adapter request counts'
    plan%sd=0;plan%rd=0
    do r=2,n_process
       plan%sd(r)=plan%sd(r-1)+plan%sc(r-1)
       plan%rd(r)=plan%rd(r-1)+plan%rc(r-1)
    end do
    cursor=plan%sd
    do i=1,n
       r=owner(unique(1,i))+1;cursor(r)=cursor(r)+1;a=cursor(r)
       outgoing(:,a)=unique(:,i);position(i)=a
    end do
    allocate(incoming(4,sum(plan%rc)),plan%source(3,sum(plan%rc)))
    call MPI_Alltoallv(outgoing,4*plan%sc,4*plan%sd,MPI_INTEGER,incoming,4*plan%rc,4*plan%rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform adapter request delivery'
    do i=1,size(incoming,2)
       if(owner(incoming(1,i))/=rank) error stop 'transform adapter source owner'
       r=loc_id(incoming(1,i))+1;a=incoming(2,i)
       if(r<1.or.r>size(grid)) error stop 'transform adapter source Domain'
       if(a<1.or.a>grid(r)%node%length) error stop 'transform adapter source node'
       plan%source(:,i)=[r,a,incoming(3,i)]
    end do
    plan%index=position(map);plan%target=target;plan%sign=sign
    plan%send_peer=pack([(r,r=1,n_process)],plan%rc>0)
    plan%receive_peer=pack([(r,r=1,n_process)],plan%sc>0)
    allocate(plan%sent(plan%columns,sum(plan%rc)),plan%received(plan%columns,n), &
         plan%request(size(plan%send_peer)+size(plan%receive_peer)))
  end subroutine

  subroutine read_transform_adapter(plan,context,fields,values)
    type(Transform_Adapter),intent(inout) :: plan
    type(MPI_Comm),intent(in) :: context
    type(Float_Field),intent(in) :: fields(N_VARIABLE,zlevels)
    real(dp),intent(inout) :: values(:,:)
    integer :: i,j,k,d,p,e,r,q,nreq,ierr
    if(work325_enabled()) call work325_record('transform.adapter',plan%tag, &
         [int(size(plan%source,2),int64),int(size(plan%target),int64), &
         int(plan%columns,int64)*int(sum(plan%rc),int64),int(plan%columns,int64)*int(plan%sc(rank+1),int64), &
         int(plan%columns,int64)*int(size(plan%target),int64),int(size(plan%sent),int64)+int(size(plan%received),int64), &
         int(size(plan%request),int64)], &
         8_int64*int(plan%columns,int64)*int(plan%rc,int64), &
         8_int64*int(plan%columns,int64)*int(plan%sc,int64))
    do i=1,size(plan%source,2)
       d=plan%source(1,i);p=plan%source(2,i);e=plan%source(3,i)
       do k=1,zlevels
          if(plan%component==1) then
             plan%sent(2*k-1,i)=fields(scalars(1),k)%data(d)%elts(p)
             plan%sent(2*k,i)=fields(scalars(2),k)%data(d)%elts(p)
          else
             plan%sent(k,i)=fields(S_VELO,k)%data(d)%elts(3*(p-1)+e)
          end if
       end do
    end do
    nreq=0
    do j=1,size(plan%receive_peer)
       r=plan%receive_peer(j)
       if(r==rank+1) cycle
       nreq=nreq+1
       call MPI_Irecv(plan%received(:,plan%sd(r)+1:plan%sd(r)+plan%sc(r)),plan%columns*plan%sc(r), &
            MPI_DOUBLE_PRECISION,r-1,plan%tag,context,plan%request(nreq),ierr)
       if(ierr/=MPI_SUCCESS) error stop 'transform adapter receive'
    end do
    do j=1,size(plan%send_peer)
       r=plan%send_peer(j)
       if(r==rank+1) cycle
       nreq=nreq+1
       call MPI_Isend(plan%sent(:,plan%rd(r)+1:plan%rd(r)+plan%rc(r)),plan%columns*plan%rc(r), &
            MPI_DOUBLE_PRECISION,r-1,plan%tag,context,plan%request(nreq),ierr)
       if(ierr/=MPI_SUCCESS) error stop 'transform adapter send'
    end do
    r=rank+1
    if(plan%sc(r)/=plan%rc(r)) error stop 'transform adapter self extent'
    plan%received(:,plan%sd(r)+1:plan%sd(r)+plan%sc(r))=plan%sent(:,plan%rd(r)+1:plan%rd(r)+plan%rc(r))
    call MPI_Waitall(nreq,plan%request(:nreq),MPI_STATUSES_IGNORE,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform adapter completion'
    do i=1,size(plan%target)
       p=plan%target(i);q=plan%index(i)
       if(plan%sign(i)==1) then
          values(p,:)=plan%received(:,q)
       else
          values(p,:)=-plan%received(:,q)
       end if
    end do
  end subroutine
end module
