module parallel_block_rhs_exchange_mod
  use audit316_mod, only : audit316_i,audit316_i64,audit316_l,audit316_begin,audit316_end
  use parallel_block_profile_mod, only : detail_enter,detail_leave,rhs_detail_enter,rhs_detail_leave, &
       DP_RHS_EXCHANGE,DP_RHS_PACK,DP_RHS_POST,DP_RHS_WAIT,DP_RHS_INSTALL
  ! Topology-time final-writer selection and phase-bounded derived-field halos.
  ! Identities describe physical storage aliases, NOT geometric coordinates.
  ! The caller supplies the accepted patch/event ordering and orientation.
  ! Runtime contains no dependency interpreter or global collective.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  implicit none
  private
  public :: Rhs_Exchange, compile_rhs_exchange, exchange_rhs_values, release_rhs_exchange
  public :: rhs_exchange_cost,rhs_exchange_source_slots,deliver_rhs_indexed_metadata,rhs_exchange_buffer_values
  public :: Rhs_Owner_Directory,compile_rhs_owner_directory,lookup_rhs_owner_directory
  public :: rhs_exchange_backward_demand,rhs_exchange_restrict,rhs_exchange_target_count
  public :: audit_rhs_exchange316
  integer, parameter :: KEY_SIZE=4, DEFINITION_SIZE=8
  type :: Rhs_Owner_Directory
     private
     ! Borrowed only during preparation; allocatables expire with its scope.
     type(MPI_Comm) :: comm=MPI_COMM_NULL
     integer :: me=-1,nrank=0
     integer,allocatable :: record(:,:),table(:)
  end type
  type :: Rhs_Exchange
     private
     type(MPI_Comm) :: comm=MPI_COMM_NULL
     integer(int64) :: topology=-1_int64, masks=-1_int64, last_state=-1_int64
     integer :: phase=-1,nslot=0,nlayer=0,me=-1,nrank=0,tag=0,profile_region=DP_RHS_EXCHANGE
     logical :: owns_comm=.true.,variable_columns=.false.
     integer :: last_columns=0
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),source(:),target(:),address(:),orientation(:)
     integer, allocatable :: send_peer(:),receive_peer(:),target_peer(:)
     real(dp), allocatable :: send(:),received(:)
     type(MPI_Request), allocatable :: requests(:)
  end type
contains
  subroutine audit_rhs_exchange316(plan,label)
    ! Final committed route only; requests/communicator handles are process-local
    ! MPI implementation identities, not canonical peers or numerical addresses.
    type(Rhs_Exchange),intent(in) :: plan
    character(*),optional,intent(in) :: label
    if(present(label)) call audit316_begin(label,plan%masks)
    call audit316_l('route.allocated',[allocated(plan%sc)])
    if(.not.allocated(plan%sc)) then
       if(present(label)) call audit316_end()
       return
    end if
    call audit316_i64('route.epoch',[plan%topology,plan%masks,plan%last_state])
    call audit316_i('route.header',[plan%phase,plan%nslot,plan%nlayer,plan%me,plan%nrank,plan%tag,plan%last_columns])
    call audit316_l('route.flags',[plan%owns_comm,plan%variable_columns])
    call audit316_i('route.sc',plan%sc)
    call audit316_i('route.sd',plan%sd)
    call audit316_i('route.rc',plan%rc)
    call audit316_i('route.rd',plan%rd)
    call audit316_i('route.source',plan%source)
    call audit316_i('route.target',plan%target)
    call audit316_i('route.address',plan%address)
    call audit316_i('route.orientation',plan%orientation)
    call audit316_i('route.send_peer',plan%send_peer)
    call audit316_i('route.receive_peer',plan%receive_peer)
    call audit316_i('route.target_peer',plan%target_peer)
    call audit316_i('route.allocations',[size(plan%send),size(plan%received),size(plan%requests)])
    if(present(label)) call audit316_end()
  end subroutine

  subroutine die(comm,message)
    type(MPI_Comm), intent(in) :: comm
    character(*), intent(in) :: message
    integer :: ierr
    write(*,'(2a)') 'Block RHS exchange: ',message
    call MPI_Abort(comm,1,ierr)
    error stop 'RHS exchange abort returned'
  end subroutine

  subroutine mpi_check(comm,ierr,message)
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: ierr
    character(*), intent(in) :: message
    if (ierr/=MPI_SUCCESS) call die(comm,message)
  end subroutine

  subroutine prefix(count,disp)
    integer, intent(in) :: count(:)
    integer, intent(out) :: disp(:)
    integer :: i
    disp=0
    do i=2,size(count)
       disp(i)=disp(i-1)+count(i-1)
    end do
  end subroutine

  integer function bucket(key,n) result(h)
    integer, intent(in) :: key(:),n
    integer(int64) :: bits
    integer :: i
    bits=0_int64
    do i=1,size(key)
       ! Bound every intermediate to avoid compiler-dependent signed overflow.
       bits=modulo(bits*65599_int64+int(key(i),int64),2147483647_int64)
    end do
    h=int(modulo(bits,int(n,int64)))+1
  end function

  subroutine deliver(comm,record,destination,received)
    ! Metadata only, used during compilation. No replicated global directory.
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in) :: record(:,:),destination(:)
    integer, allocatable, intent(out) :: received(:,:)
    integer, allocatable :: sc(:),sd(:),rc(:),rd(:),cursor(:),packed(:,:)
    integer :: n,width,i,r,ierr
    call MPI_Comm_size(comm,n,ierr);call mpi_check(comm,ierr,'metadata communicator')
    width=size(record,1)
    if (size(record,2)/=size(destination)) call die(comm,'metadata extent')
    if (any(destination<0).or.any(destination>=n)) call die(comm,'metadata rank')
    allocate(sc(n),sd(n),rc(n),rd(n),cursor(n),packed(width,size(destination)))
    sc=0
    do i=1,size(destination)
       r=destination(i)+1;sc(r)=sc(r)+1
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    call mpi_check(comm,ierr,'metadata counts')
    call prefix(sc,sd);call prefix(rc,rd);cursor=sd
    do i=1,size(destination)
       r=destination(i)+1;cursor(r)=cursor(r)+1
       packed(:,cursor(r))=record(:,i)
    end do
    allocate(received(width,sum(rc)))
    call MPI_Alltoallv(packed,width*sc,width*sd,MPI_INTEGER,received,width*rc,width*rd,MPI_INTEGER,comm,ierr)
    call mpi_check(comm,ierr,'metadata delivery')
  end subroutine

  integer function rank_slot_bucket(key,n) result(h)
    integer,intent(in) :: key(2),n
    integer(int64) :: bits
    ! A rank prefix plus an unmodified dense slot creates overlapping linear
    ! probe clusters. Mix the slot before reduction to the local table size.
    ! All intermediates stay below 2**49 for nonnegative default integers.
    bits=modulo(int(key(2),int64)*104729_int64+int(key(1),int64)*65599_int64,2147483647_int64)
    bits=ieor(bits,shiftr(bits,16))
    bits=modulo(bits*104729_int64,2147483647_int64)
    h=int(modulo(bits,int(n,int64)))+1
  end function

  subroutine deliver_rhs_indexed_metadata(comm,source,home,values,topology,phase,remote_values)
    ! Preparation-only pull from an already proved rank/slot directory. There
    ! is one producer per slot: no writer election or persistent halo is needed.
    ! Buffers and requests expire on return; all four collectives are blocking.
    type(MPI_Comm),intent(in) :: comm
    real(dp),intent(in) :: source(:,:)
    integer,intent(in) :: home(:,:)
    real(dp),allocatable,intent(out) :: values(:,:)
    integer(int64),intent(in) :: topology
    integer,intent(in) :: phase
    integer(int64),optional,intent(out) :: remote_values(2)
    integer(int64) :: signature(6),checked(6)
    integer,allocatable :: sc(:),sd(:),rc(:),rd(:),cursor(:),request(:),incoming(:),index(:)
    real(dp),allocatable :: reply(:,:),received(:,:)
    integer,allocatable :: table(:),unique(:,:),mapping(:),packed_index(:)
    integer :: nr,nq,nc,i,r,a,ierr,h,j,u,me
    call MPI_Comm_size(comm,nr,ierr);call mpi_check(comm,ierr,'indexed metadata communicator')
    if(size(home,1)/=2) call die(comm,'indexed metadata home extent')
    if(any(home(1,:)<0).or.any(home(1,:)>=nr).or.any(home(2,:)<1)) &
         call die(comm,'indexed metadata home address')
    nq=size(home,2);nc=size(source,2)
    if(nc<1.or.topology<0.or.phase<1) call die(comm,'indexed metadata epoch or extent')
    signature(:3)=[topology,int(phase,int64),int(nc,int64)];signature(4:)=-signature(:3)
    call MPI_Allreduce(signature,checked,6,MPI_INTEGER8,MPI_MIN,comm,ierr)
    call mpi_check(comm,ierr,'indexed metadata signature')
    if(any(checked(:3)/=-checked(4:))) call die(comm,'inconsistent indexed metadata epoch')
    ! Match the original halo compiler's wire-volume bound: repeated aliases
    ! request one record per producer rank/slot, then expand locally in order.
    allocate(table(max(1,2*nq+1)),unique(2,nq),mapping(nq));table=0;u=0
    do i=1,nq
       h=rank_slot_bucket(home(:,i),size(table))
       do
          j=table(h)
          if(j==0) then
             u=u+1;j=u;table(h)=j;unique(:,j)=home(:,i);exit
          end if
          if(all(unique(:,j)==home(:,i))) exit
          h=modulo(h,size(table))+1
       end do
       mapping(i)=j
    end do
    deallocate(table)
    allocate(sc(nr),sd(nr),rc(nr),rd(nr),cursor(nr),request(u),index(u),packed_index(nq));sc=0
    do i=1,u
       r=unique(1,i)+1;sc(r)=sc(r)+1
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    call mpi_check(comm,ierr,'indexed metadata request counts')
    call prefix(sc,sd);call prefix(rc,rd);cursor=sd
    do i=1,u
       r=unique(1,i)+1;cursor(r)=cursor(r)+1;a=cursor(r)
       request(a)=unique(2,i);index(i)=a
    end do
    packed_index=index(mapping)
    deallocate(unique,mapping,index)
    allocate(incoming(sum(rc)))
    call MPI_Alltoallv(request,sc,sd,MPI_INTEGER,incoming,rc,rd,MPI_INTEGER,comm,ierr)
    call mpi_check(comm,ierr,'indexed metadata request slots')
    if(any(incoming<1).or.any(incoming>size(source,1))) call die(comm,'indexed metadata source address')
    allocate(reply(nc,size(incoming)),received(nc,u),values(nq,nc))
    do i=1,size(incoming)
       reply(:,i)=source(incoming(i),:)
    end do
    call MPI_Alltoallv(reply,nc*rc,nc*rd,MPI_DOUBLE_PRECISION,received,nc*sc,nc*sd, &
         MPI_DOUBLE_PRECISION,comm,ierr)
    call mpi_check(comm,ierr,'indexed metadata reply')
    if(present(remote_values)) then
       call MPI_Comm_rank(comm,me,ierr);call mpi_check(comm,ierr,'indexed metadata rank')
       remote_values=int(nc,int64)*[int(sum(rc)-rc(me+1),int64),int(sum(sc)-sc(me+1),int64)]
    end if
    do i=1,nq
       values(i,:)=received(:,packed_index(i))
    end do
  end subroutine

  subroutine compile_rhs_owner_directory(directory,comm,key,slot)
    ! Canonical physical ownership is common to every field phase. Publish it
    ! once per executor preparation, instead of rediscovering all interiors.
    type(Rhs_Owner_Directory),intent(out) :: directory
    type(MPI_Comm),intent(in) :: comm
    integer,intent(in) :: key(:,:),slot(:)
    integer,allocatable :: record(:,:),destination(:)
    integer :: i,j,h,n,ierr,bad,total
    directory%comm=comm
    call MPI_Comm_rank(comm,directory%me,ierr);call mpi_check(comm,ierr,'owner directory rank')
    call MPI_Comm_size(comm,directory%nrank,ierr);call mpi_check(comm,ierr,'owner directory size')
    if(any(shape(key)/=[KEY_SIZE,size(slot)]).or.any(key<=0).or.any(slot<1)) &
         call die(comm,'invalid canonical owner declaration')
    allocate(record(KEY_SIZE+2,size(slot)),destination(size(slot)))
    do i=1,size(slot)
       record(:,i)=[key(:,i),directory%me,slot(i)]
       destination(i)=bucket(key(:,i),directory%nrank)-1
    end do
    call deliver(comm,record,destination,directory%record)
    n=size(directory%record,2);allocate(directory%table(max(1,2*n+1)));directory%table=0;bad=0
    do i=1,n
       h=bucket(directory%record(:KEY_SIZE,i),size(directory%table))
       do
          j=directory%table(h)
          if(j==0) then
             directory%table(h)=i;exit
          end if
          if(all(directory%record(:KEY_SIZE,j)==directory%record(:KEY_SIZE,i))) then
             bad=bad+1;exit
          end if
          h=modulo(h,size(directory%table))+1
       end do
    end do
    call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call mpi_check(comm,ierr,'canonical owner uniqueness')
    if(total/=0) call die(comm,'duplicate canonical owner')
  end subroutine

  subroutine lookup_rhs_owner_directory(directory,key,home)
    type(Rhs_Owner_Directory),intent(in) :: directory
    integer,intent(in) :: key(:,:)
    integer,allocatable,intent(out) :: home(:,:)
    integer,allocatable :: query(:,:),destination(:),incoming(:,:),reply(:,:),answer(:,:),seen(:)
    integer :: i,j,h,n,ierr,bad,total
    n=size(key,2)
    if(size(key,1)/=KEY_SIZE.or.any(key<=0)) call die(directory%comm,'invalid canonical owner request')
    allocate(query(KEY_SIZE+2,n),destination(n))
    do i=1,n
       query(:,i)=[key(:,i),directory%me,i];destination(i)=bucket(key(:,i),directory%nrank)-1
    end do
    call deliver(directory%comm,query,destination,incoming);deallocate(query,destination)
    allocate(reply(3,size(incoming,2)),destination(size(incoming,2)));bad=0
    do i=1,size(incoming,2)
       h=bucket(incoming(:KEY_SIZE,i),size(directory%table))
       do
          j=directory%table(h)
          if(j==0) exit
          if(all(directory%record(:KEY_SIZE,j)==incoming(:KEY_SIZE,i))) exit
          h=modulo(h,size(directory%table))+1
       end do
       reply(:,i)=[incoming(KEY_SIZE+2,i),-1,0];destination(i)=incoming(KEY_SIZE+1,i)
       if(j==0) then
          bad=bad+1
       else
          reply(2:3,i)=directory%record(KEY_SIZE+1:KEY_SIZE+2,j)
       end if
    end do
    call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_SUM,directory%comm,ierr)
    call mpi_check(directory%comm,ierr,'canonical owner lookup')
    if(total/=0) call die(directory%comm,'missing canonical owner')
    call deliver(directory%comm,reply,destination,answer)
    if(size(answer,2)/=n) call die(directory%comm,'canonical owner reply extent')
    allocate(home(2,n),seen(n));seen=0
    do i=1,n
       j=answer(1,i)
       if(j<1.or.j>n) call die(directory%comm,'canonical owner reply index')
       if(seen(j)/=0) call die(directory%comm,'duplicate canonical owner reply')
       seen(j)=1;home(:,j)=answer(2:3,i)
    end do
  end subroutine

  subroutine resolve_rhs_home_writers(comm,me,nrank,nslot,definition_home,definition_order,definition_slot, &
       consumer_home,canonical_slot,answer)
    ! Only nonlocal writes and consumer requests cross ranks. The canonical
    ! owner chooses the same final writer by the accepted event ordering.
    ! All private interior definitions/queries are resolved by direct indexing.
    type(MPI_Comm),intent(in) :: comm
    integer,intent(in) :: me,nrank,nslot,definition_home(:,:),definition_order(:,:),definition_slot(:), &
         consumer_home(:,:),canonical_slot(:)
    integer,allocatable,intent(out) :: answer(:,:)
    integer,allocatable :: directory(:,:),record(:,:),destination(:),incoming(:,:),reply(:,:),returned(:,:),seen(:)
    logical,allocatable :: canonical(:)
    integer :: nd,nc,i,j,n,s,bad,total,ierr
    nd=size(definition_slot);nc=size(consumer_home,2)
    if(any(shape(definition_home)/=[2,nd]).or.size(consumer_home,1)/=2) call die(comm,'home binding extent')
    if(any(definition_home(1,:)<0).or.any(definition_home(1,:)>=nrank).or.any(definition_home(2,:)<1).or. &
         any(consumer_home(1,:)<0).or.any(consumer_home(1,:)>=nrank).or.any(consumer_home(2,:)<1)) &
         call die(comm,'home binding address')
    if(any(canonical_slot<1).or.any(canonical_slot>nslot)) call die(comm,'canonical slot extent')
    allocate(canonical(nslot),directory(4,nslot));canonical=.false.;directory=0
    do i=1,size(canonical_slot)
       s=canonical_slot(i)
       if(canonical(s)) call die(comm,'duplicate canonical slot')
       canonical(s)=.true.
    end do
    n=count(definition_home(1,:)/=me);allocate(record(5,n),destination(n));n=0;bad=0
    do i=1,nd
       if(definition_home(1,i)==me) then
          call define(definition_home(2,i),definition_order(:,i),me,definition_slot(i))
       else
          n=n+1;record(:,n)=[definition_home(2,i),definition_order(:,i),me,definition_slot(i)]
          destination(n)=definition_home(1,i)
       end if
    end do
    call deliver(comm,record,destination,incoming);deallocate(record,destination)
    do i=1,size(incoming,2)
       call define(incoming(1,i),incoming(2:3,i),incoming(4,i),incoming(5,i))
    end do
    call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call mpi_check(comm,ierr,'home writer ambiguity check')
    if(total/=0) call die(comm,'ambiguous final writer')
    deallocate(incoming)
    n=count(consumer_home(1,:)/=me);allocate(record(3,n),destination(n),answer(2,nc),seen(nc));n=0;bad=0;seen=0
    do i=1,nc
       if(consumer_home(1,i)==me) then
          answer(:,i)=lookup(consumer_home(2,i));seen(i)=1
       else
          n=n+1;record(:,n)=[consumer_home(2,i),me,i];destination(n)=consumer_home(1,i)
       end if
    end do
    call deliver(comm,record,destination,incoming);deallocate(record,destination)
    allocate(reply(3,size(incoming,2)),destination(size(incoming,2)))
    do i=1,size(incoming,2)
       reply(1,i)=incoming(3,i);reply(2:3,i)=lookup(incoming(1,i));destination(i)=incoming(2,i)
    end do
    call MPI_Allreduce(bad,total,1,MPI_INTEGER,MPI_SUM,comm,ierr)
    call mpi_check(comm,ierr,'home missing producer check')
    if(total/=0) call die(comm,'missing phase producer')
    call deliver(comm,reply,destination,returned)
    do i=1,size(returned,2)
       j=returned(1,i)
       if(j<1.or.j>nc) call die(comm,'home consumer reply index')
       if(seen(j)/=0) call die(comm,'duplicate home consumer reply')
       seen(j)=1;answer(:,j)=returned(2:3,i)
    end do
    if(any(seen/=1)) call die(comm,'missing home consumer reply')
  contains
    subroutine define(slot,order,writer,address)
      integer,intent(in) :: slot,order(2),writer,address
      integer :: t,choice
      if(slot<1.or.slot>nslot) call die(comm,'home definition outside workspace')
      if(.not.canonical(slot)) call die(comm,'home definition is not canonical')
      choice=0
      do t=1,2
         if(order(t)==directory(t,slot)) cycle
         if(order(t)>directory(t,slot)) choice=1
         if(order(t)<directory(t,slot)) choice=-1
         exit
      end do
      if(choice==1) directory(:,slot)=[order,writer,address]
      if(choice==0.and.any(directory(3:4,slot)/=[writer,address])) bad=bad+1
    end subroutine
    function lookup(slot) result(provider)
      integer,intent(in) :: slot
      integer :: provider(2)
      if(slot<1.or.slot>nslot) call die(comm,'home consumer outside workspace')
      if(.not.canonical(slot)) call die(comm,'home consumer is not canonical')
      provider=directory(3:4,slot)
      if(directory(1,slot)==0) bad=bad+1
    end function
  end subroutine

  subroutine compile_rhs_exchange(plan,comm,topology,masks,phase,nslot,nlayer, &
       definition_key,definition_order,definition_slot,consumer_key,consumer_slot,orientation,shared_tag,variable_columns, &
       definition_home,consumer_home,canonical_slot,profile_region,shared_buffers,consumer_provider)
    ! Each definition is the FINAL local write for its physical storage key.
    ! order=(accepted patch ordinal, accepted assignment ordinal) compares
    ! lexicographically. All pentagon writes must follow all direct/SW writes.
    ! Different keys remain distinct even when they have equal coordinates.
    ! Consumers explicitly name the required producer alias; orientation only
    ! changes sign. Vector-component permutations belong in consumer_key.
    type(Rhs_Exchange), intent(inout) :: plan
    type(MPI_Comm), intent(in) :: comm
    integer(int64), intent(in) :: topology,masks
    integer, intent(in) :: phase,nslot,nlayer
    integer, intent(in) :: definition_key(:,:),definition_order(:,:),definition_slot(:)
    integer, intent(in) :: consumer_key(:,:),consumer_slot(:),orientation(:)
    logical, optional, intent(in) :: shared_buffers
    integer, optional, intent(in) :: profile_region
    integer, optional, intent(in) :: shared_tag
    logical, optional, intent(in) :: variable_columns
    integer, optional, intent(in) :: definition_home(:,:),consumer_home(:,:),canonical_slot(:),consumer_provider(:,:)
    integer, allocatable :: outgoing(:,:),destination(:),definitions(:,:),table(:),query(:,:),answer(:,:),active(:)
    integer, allocatable :: incoming(:,:),reply(:,:),replies(:,:),slots(:),unique(:,:),mapping(:),order(:),cursor(:)
    integer :: ierr,nd,nc,i,j,h,t,n,u,r,a,failures,total_failures,choice
    integer(int64) :: signature(7),lower(7),upper(7)
    if (plan%comm/=MPI_COMM_NULL) call die(comm,'release previous plan before recompiling')
    if (topology<0.or.masks<0.or.phase<1.or.nslot<0.or.nlayer<1) call die(comm,'invalid plan epoch or extent')
    if(present(profile_region)) plan%profile_region=profile_region
    nd=size(definition_slot);nc=size(consumer_slot)
    if (any(shape(definition_key)/=[KEY_SIZE,nd]).or.any(shape(definition_order)/=[2,nd]).or. &
         any(shape(consumer_key)/=[KEY_SIZE,nc]).or.size(orientation)/=nc) call die(comm,'plan argument extent')
    if (any(definition_slot<1).or.any(definition_slot>nslot).or.any(consumer_slot<1).or. &
         any(consumer_slot>nslot).or.any(abs(orientation)/=1)) call die(comm,'plan slot or orientation')
    if (any(definition_key<=0).or.any(consumer_key<=0).or.any(definition_order<=0)) &
         call die(comm,'invalid producer identity or order')
    allocate(slots(nslot));slots=0
    do i=1,nc
       if (slots(consumer_slot(i))/=0) call die(comm,'duplicate consumer destination')
       slots(consumer_slot(i))=i
    end do
    deallocate(slots)
    ! An isolated communicator prevents any overlap with accepted halo tags.
    ! It lives for one topology/mask plan and is explicitly released on rebuild.
    plan%owns_comm=.not.present(shared_tag);plan%tag=0
    if(present(shared_tag)) then
       ! Caller owns an isolated context until every borrowing route is released.
       ! Distinct phase/channel tags stay within MPI's guaranteed portable range.
       if(shared_tag<0.or.shared_tag>32767) call die(comm,'shared route tag outside portable range')
       plan%comm=comm;plan%tag=shared_tag
    else
       call MPI_Comm_dup(comm,plan%comm,ierr);call mpi_check(comm,ierr,'plan communicator duplication')
    end if
    call MPI_Comm_rank(plan%comm,plan%me,ierr);call mpi_check(comm,ierr,'plan rank')
    call MPI_Comm_size(plan%comm,plan%nrank,ierr);call mpi_check(comm,ierr,'plan size')
    plan%variable_columns=.false.
    if(present(variable_columns)) plan%variable_columns=variable_columns
    plan%last_columns=nlayer
    signature=[topology,masks,int(phase,int64),int(nlayer,int64),int(plan%tag,int64), &
         int(merge(1,0,plan%variable_columns),int64), &
         int(merge(2,merge(1,0,present(definition_home)),present(consumer_provider)),int64)]
    call MPI_Allreduce(signature,lower,7,MPI_INTEGER8,MPI_MIN,plan%comm,ierr)
    call mpi_check(comm,ierr,'plan signature minimum')
    call MPI_Allreduce(signature,upper,7,MPI_INTEGER8,MPI_MAX,plan%comm,ierr)
    call mpi_check(comm,ierr,'plan signature maximum')
    if (any(lower/=upper)) call die(plan%comm,'inconsistent distributed plan epoch')
    plan%topology=topology;plan%masks=masks;plan%phase=phase;plan%nslot=nslot;plan%nlayer=nlayer
    if(present(consumer_provider)) then
       ! The caller supplies providers from a UNIQUE canonical ownership
       ! directory for a phase with no competing writes. This is not a writer
       ! election shortcut for alias/in-place phases: they retain that proof.
       if(nd/=0.or.present(definition_home).or.present(consumer_home).or.present(canonical_slot)) &
            call die(comm,'prebound provider with competing definitions')
       if(any(shape(consumer_provider)/=[2,nc])) call die(comm,'prebound provider extent')
       if(any(consumer_provider(1,:)<0).or.any(consumer_provider(1,:)>=plan%nrank).or. &
            any(consumer_provider(2,:)<1)) call die(comm,'prebound provider address')
       answer=consumer_provider
    else if(present(definition_home)) then
       if(.not.present(consumer_home).or..not.present(canonical_slot)) call die(comm,'incomplete home binding')
       if(size(consumer_home,2)/=nc) call die(comm,'home consumer extent')
       call resolve_rhs_home_writers(plan%comm,plan%me,plan%nrank,nslot,definition_home, &
            definition_order,definition_slot,consumer_home,canonical_slot,answer)
    else
       if(present(consumer_home).or.present(canonical_slot)) call die(comm,'incomplete home binding')
       allocate(outgoing(DEFINITION_SIZE,nd),destination(nd))
       do i=1,nd
          outgoing(:,i)=[definition_key(:,i),definition_order(:,i),plan%me,definition_slot(i)]
          destination(i)=bucket(definition_key(:,i),plan%nrank)-1
       end do
       call deliver(plan%comm,outgoing,destination,definitions)
       deallocate(outgoing,destination)
       n=size(definitions,2);allocate(table(max(1,2*n+1)));table=0;failures=0
       do i=1,n
          h=bucket(definitions(:KEY_SIZE,i),size(table))
          do
             j=table(h)
             if (j==0) then
                table(h)=i;exit
             end if
             if (all(definitions(:KEY_SIZE,j)==definitions(:KEY_SIZE,i))) then
                choice=0
                do t=KEY_SIZE+1,KEY_SIZE+2
                   if (definitions(t,i)==definitions(t,j)) cycle
                   if (definitions(t,i)>definitions(t,j)) choice=1
                   if (definitions(t,i)<definitions(t,j)) choice=-1
                   exit
                end do
                if (choice==1) table(h)=i
                if (choice==0.and.any(definitions(7:8,i)/=definitions(7:8,j))) failures=failures+1
                exit
             end if
             h=modulo(h,size(table))+1
          end do
       end do
       call MPI_Allreduce(failures,total_failures,1,MPI_INTEGER,MPI_SUM,plan%comm,ierr)
       call mpi_check(comm,ierr,'writer ambiguity check')
       if (total_failures/=0) call die(plan%comm,'ambiguous final writer')
       allocate(query(KEY_SIZE+2,nc),destination(nc))
       do i=1,nc
          query(:,i)=[consumer_key(:,i),plan%me,i]
          destination(i)=bucket(consumer_key(:,i),plan%nrank)-1
       end do
       call deliver(plan%comm,query,destination,incoming)
       deallocate(query,destination)
       allocate(reply(3,size(incoming,2)),destination(size(incoming,2)))
       failures=0
       do i=1,size(incoming,2)
          h=bucket(incoming(:KEY_SIZE,i),size(table))
          do
             j=table(h)
             if (j==0) exit
             if (all(definitions(:KEY_SIZE,j)==incoming(:KEY_SIZE,i))) exit
             h=modulo(h,size(table))+1
          end do
          reply(:,i)=[incoming(KEY_SIZE+2,i),-1,0]
          if (j==0) then
             failures=failures+1
          else
             reply(2:3,i)=definitions(7:8,j)
          end if
          destination(i)=incoming(KEY_SIZE+1,i)
       end do
       call MPI_Allreduce(failures,total_failures,1,MPI_INTEGER,MPI_SUM,plan%comm,ierr)
       call mpi_check(comm,ierr,'missing producer check')
       if (total_failures/=0) call die(plan%comm,'missing phase producer')
       call deliver(plan%comm,reply,destination,replies)
       deallocate(reply,destination,incoming,definitions,table)
       if (size(replies,2)/=nc) call die(plan%comm,'consumer reply count')
       allocate(answer(2,nc),slots(nc));slots=0
       do i=1,nc
          j=replies(1,i)
          if (j<1.or.j>nc) call die(plan%comm,'consumer reply index')
          if (slots(j)/=0) call die(plan%comm,'duplicate consumer reply')
          slots(j)=1;answer(:,j)=replies(2:3,i)
       end do
       deallocate(slots,replies)
    end if
    ! A canonical slot already written on this rank needs no staging copy.
    ! Prune only true identity deliveries; signed aliases and cycles still use
    ! the complete pre-install snapshot, including local sources needed remotely.
    active=pack([(i,i=1,nc)],answer(1,:)/=plan%me.or.answer(2,:)/=consumer_slot.or.orientation/=1)
    answer=answer(:,active);nc=size(active)
    ! Deduplicate requests by physical producer slot on each peer. Several
    ! transformed consumer aliases share one value on the wire.
    allocate(table(max(1,2*nc+1)),unique(2,nc),mapping(nc),order(nc))
    table=0;u=0
    do i=1,nc
       h=rank_slot_bucket(answer(:,i),size(table))
       do
          j=table(h)
          if (j==0) then
             u=u+1;j=u;table(h)=j;unique(:,j)=answer(:,i);exit
          end if
          if (all(unique(:,j)==answer(:,i))) exit
          h=modulo(h,size(table))+1
       end do
       mapping(i)=j
    end do
    allocate(plan%sc(plan%nrank),plan%sd(plan%nrank),plan%rc(plan%nrank),plan%rd(plan%nrank),cursor(plan%nrank))
    plan%rc=0
    do i=1,u
       r=unique(1,i)+1;plan%rc(r)=plan%rc(r)+1
    end do
    call prefix(plan%rc,plan%rd);cursor=plan%rd
    allocate(slots(u))
    do i=1,u
       r=unique(1,i)+1;cursor(r)=cursor(r)+1;a=cursor(r)
       slots(a)=unique(2,i);order(i)=a
    end do
    call MPI_Alltoall(plan%rc,1,MPI_INTEGER,plan%sc,1,MPI_INTEGER,plan%comm,ierr)
    call mpi_check(comm,ierr,'phase request counts')
    call prefix(plan%sc,plan%sd)
    plan%send_peer=pack([(i,i=1,plan%nrank)],plan%sc>0)
    plan%receive_peer=pack([(i,i=1,plan%nrank)],plan%rc>0)
    allocate(plan%source(sum(plan%sc)))
    call MPI_Alltoallv(slots,plan%rc,plan%rd,MPI_INTEGER,plan%source,plan%sc,plan%sd,MPI_INTEGER,plan%comm,ierr)
    call mpi_check(comm,ierr,'phase request addresses')
    if (any(plan%source<1).or.any(plan%source>nslot)) call die(plan%comm,'remote source outside workspace')
    plan%target=consumer_slot(active);plan%orientation=orientation(active)
    allocate(plan%address(nc),plan%target_peer(nc))
    do i=1,nc
       r=answer(1,i)+1
       plan%address(i)=order(mapping(i))-1;plan%target_peer(i)=r
    end do
    ! One int64 state stamp is embedded in each peer message. Even otherwise
    ! valid local transactions must reject a peer sending a different stage.
    allocate(plan%requests(size(plan%send_peer)+size(plan%receive_peer)))
    if(present(shared_buffers)) then
       if(shared_buffers) then
          allocate(plan%send(0),plan%received(0))
          return
       end if
    end if
    allocate(plan%send(nlayer*sum(plan%sc)+plan%nrank),plan%received(nlayer*sum(plan%rc)+plan%nrank))
  end subroutine

  integer function rhs_exchange_target_count(plan) result(n)
    type(Rhs_Exchange),intent(in) :: plan
    if(plan%comm==MPI_COMM_NULL) error stop 'demand requires compiled route'
    n=size(plan%target)
  end function

  subroutine select_demand_requests(plan,wanted,requested,exported)
    ! Reverse the already-proved provider binding. There is no writer election,
    ! numerical value, new communicator or runtime work in this preparation pass.
    type(Rhs_Exchange),intent(in) :: plan
    logical,intent(in) :: wanted(:,:)
    logical,allocatable,intent(out) :: requested(:,:),exported(:,:)
    integer :: i,c,nc,ierr
    if(plan%comm==MPI_COMM_NULL.or.plan%last_state>=0) error stop 'demand route lifetime'
    if(size(wanted,1)/=size(plan%target).or.size(wanted,2)<1) error stop 'demand consumer extent'
    nc=size(wanted,2)
    allocate(requested(nc,sum(plan%rc)),exported(nc,sum(plan%sc)));requested=.false.
    do i=1,size(plan%target)
       do c=1,nc
          if(wanted(i,c)) requested(c,plan%address(i)+1)=.true.
       end do
    end do
    call MPI_Alltoallv(requested,nc*plan%rc,nc*plan%rd,MPI_LOGICAL, &
         exported,nc*plan%sc,nc*plan%sd,MPI_LOGICAL,plan%comm,ierr)
    call mpi_check(plan%comm,ierr,'backward demand requests')
  end subroutine

  subroutine rhs_exchange_backward_demand(plan,live,keep)
    ! The halo is one simultaneous pre-install snapshot. Kill all destinations
    ! before marking requested producer versions; this preserves signed cycles.
    type(Rhs_Exchange),intent(in) :: plan
    logical,intent(inout) :: live(:,:),keep(:)
    logical,allocatable :: wanted(:,:),requested(:,:),exported(:,:)
    integer :: i,c
    if(size(live,1)/=plan%nslot.or.size(keep)/=size(plan%target)) error stop 'backward demand extent'
    allocate(wanted(size(plan%target),size(live,2)))
    wanted(:,:)=live(plan%target,:)
    keep=keep.or.any(wanted,dim=2)
    call select_demand_requests(plan,wanted,requested,exported)
    live(plan%target,:)=.false.
    do i=1,size(plan%source)
       do c=1,size(live,2)
          if(exported(c,i)) live(plan%source(i),c)=.true.
       end do
    end do
  end subroutine

  subroutine rhs_exchange_restrict(plan,keep)
    ! Retain a subset of established consumer bindings. One reverse flag exchange
    ! selects matching exports; existing ordering supplies reciprocal new counts.
    ! No second compiler, consensus, directory or new route context is needed.
    type(Rhs_Exchange),intent(inout) :: plan
    logical,intent(in) :: keep(:)
    logical,allocatable :: wanted(:,:),requested(:,:),exported(:,:)
    integer,allocatable :: position(:),sc(:),rc(:),selected(:)
    integer :: i,r,n
    logical :: shared
    if(size(keep)/=size(plan%target)) error stop 'restricted demand extent'
    allocate(wanted(size(keep),1));wanted(:,1)=keep
    call select_demand_requests(plan,wanted,requested,exported)
    allocate(sc(plan%nrank),rc(plan%nrank),position(size(requested,2)));position=0;n=0
    do i=1,size(position)
       if(.not.requested(1,i)) cycle
       n=n+1;position(i)=n
    end do
    do r=1,plan%nrank
       sc(r)=count(exported(1,plan%sd(r)+1:plan%sd(r)+plan%sc(r)))
       rc(r)=count(requested(1,plan%rd(r)+1:plan%rd(r)+plan%rc(r)))
    end do
    selected=pack([(i,i=1,size(keep))],keep)
    plan%address=position(plan%address(selected)+1)-1
    if(any(plan%address<0)) error stop 'restricted request missing provider'
    plan%target=plan%target(selected);plan%orientation=plan%orientation(selected)
    plan%target_peer=plan%target_peer(selected)
    plan%source=pack(plan%source,exported(1,:))
    plan%sc=sc;plan%rc=rc;call prefix(plan%sc,plan%sd);call prefix(plan%rc,plan%rd)
    plan%send_peer=pack([(r,r=1,plan%nrank)],plan%sc>0)
    plan%receive_peer=pack([(r,r=1,plan%nrank)],plan%rc>0)
    deallocate(plan%requests);allocate(plan%requests(size(plan%send_peer)+size(plan%receive_peer)))
    shared=size(plan%send)==0.and.size(plan%received)==0
    if(.not.shared) then
       deallocate(plan%send,plan%received)
       allocate(plan%send(plan%nlayer*sum(plan%sc)+plan%nrank), &
            plan%received(plan%nlayer*sum(plan%rc)+plan%nrank))
    end if
  end subroutine

  subroutine exchange_rhs_values(plan,topology,masks,phase,state,values,scratch)
    type(Rhs_Exchange), intent(inout), target, asynchronous :: plan
    integer(int64), intent(in) :: topology,masks,state
    integer, intent(in) :: phase
    real(dp), intent(inout) :: values(:,:)
    real(dp),optional,intent(inout),target,asynchronous :: scratch(:)
    real(dp),pointer,asynchronous :: sent(:),received(:)
    integer :: i,k,r,n,ierr,a,b,tag,peer,columns,ns,nr
    if (plan%comm==MPI_COMM_NULL) error stop 'RHS exchange missing compiled plan'
    if (topology/=plan%topology.or.masks/=plan%masks.or.phase/=plan%phase) &
         call die(plan%comm,'stale topology, mask or phase')
    if (state<=plan%last_state) call die(plan%comm,'stale state or duplicate phase')
    columns=size(values,2)
    if(size(values,1)/=plan%nslot.or.columns<1.or.columns>plan%nlayer) &
         call die(plan%comm,'runtime workspace extent')
    if(.not.plan%variable_columns.and.columns/=plan%nlayer) call die(plan%comm,'runtime workspace extent')
    ns=columns*sum(plan%sc)+plan%nrank;nr=columns*sum(plan%rc)+plan%nrank
    if(present(scratch)) then
       if(size(plan%send)/=0.or.size(scratch)<ns+nr) call die(plan%comm,'shared exchange scratch extent')
       sent=>scratch(:ns);received=>scratch(ns+1:ns+nr)
    else
       if(size(plan%send)<ns.or.size(plan%received)<nr) call die(plan%comm,'missing exchange scratch')
       sent=>plan%send;received=>plan%received
    end if
    call detail_enter(plan%profile_region)
    call rhs_detail_enter(DP_RHS_PACK)
    ! Snapshot ALL sends before modifying any destination, including local
    ! cycles and aliases which are both consumers and remote producers.
    do peer=1,size(plan%send_peer)
       r=plan%send_peer(peer)
       a=columns*plan%sd(r)+r
       sent(a)=transfer(state,0.0_dp)
       do i=plan%sd(r)+1,plan%sd(r)+plan%sc(r)
          do k=1,columns
             sent((i-1)*columns+r+k)=values(plan%source(i),k)
          end do
       end do
    end do
    call rhs_detail_leave(DP_RHS_PACK)
    call rhs_detail_enter(DP_RHS_POST)
    n=0;tag=plan%tag
    do peer=1,size(plan%receive_peer)
       r=plan%receive_peer(peer)
       if (r==plan%me+1) cycle
       n=n+1;a=columns*plan%rd(r)+r;b=a+columns*plan%rc(r)
       call MPI_Irecv(received(a:b),b-a+1,MPI_DOUBLE_PRECISION, &
            r-1,tag,plan%comm,plan%requests(n),ierr)
       call mpi_check(plan%comm,ierr,'phase receive')
    end do
    do peer=1,size(plan%send_peer)
       r=plan%send_peer(peer)
       if (r==plan%me+1) cycle
       n=n+1;a=columns*plan%sd(r)+r;b=a+columns*plan%sc(r)
       call MPI_Isend(sent(a:b),b-a+1,MPI_DOUBLE_PRECISION, &
            r-1,tag,plan%comm,plan%requests(n),ierr)
       call mpi_check(plan%comm,ierr,'phase send')
    end do
    r=plan%me+1
    a=columns*plan%rd(r)+r;b=columns*plan%sd(r)+r
    if (plan%rc(r)>0) received(a:a+columns*plan%rc(r))=sent(b:b+columns*plan%sc(r))
    call rhs_detail_leave(DP_RHS_POST)
    call rhs_detail_enter(DP_RHS_WAIT)
    call MPI_Waitall(n,plan%requests(:n),MPI_STATUSES_IGNORE,ierr)
    call mpi_check(plan%comm,ierr,'phase completion')
    call rhs_detail_leave(DP_RHS_WAIT)
    call rhs_detail_enter(DP_RHS_INSTALL)
    do peer=1,size(plan%receive_peer)
       r=plan%receive_peer(peer)
       a=columns*plan%rd(r)+r
       if (transfer(received(a),0_int64)/=state) call die(plan%comm,'peer state epoch mismatch')
    end do
    do i=1,size(plan%target)
       do k=1,columns
          values(plan%target(i),k)=plan%orientation(i)*received(plan%address(i)*columns+plan%target_peer(i)+k)
       end do
    end do
    plan%last_state=state;plan%last_columns=columns
    call rhs_detail_leave(DP_RHS_INSTALL)
    call detail_leave(plan%profile_region)
  end subroutine

  subroutine rhs_exchange_source_slots(plan,source)
    ! Topology-time provider binding. Repeated peer requests can name the same
    ! source; the state adapter deduplicates these before reading resident data.
    type(Rhs_Exchange),intent(in) :: plan
    integer,allocatable,intent(out) :: source(:)
    if(plan%comm==MPI_COMM_NULL) error stop 'RHS source binding needs a compiled route'
    source=plan%source
  end subroutine

  integer function rhs_exchange_buffer_values(plan) result(n)
    type(Rhs_Exchange),intent(in) :: plan
    if(plan%comm==MPI_COMM_NULL) error stop 'buffer capacity requires compiled route'
    n=plan%nlayer*(sum(plan%sc)+sum(plan%rc))+2*plan%nrank
  end function

  subroutine rhs_exchange_cost(plan,send_values,receive_values,send_peers,receive_peers,storage_bytes,header_bytes)
    type(Rhs_Exchange), intent(in) :: plan
    integer(int64), intent(out) :: send_values,receive_values,storage_bytes
    integer, intent(out) :: send_peers,receive_peers
    integer(int64), optional, intent(out) :: header_bytes
    integer :: r
    if (plan%comm==MPI_COMM_NULL) error stop 'RHS exchange missing cost plan'
    r=plan%me+1
    send_values=int(sum(plan%sc)-plan%sc(r),int64)*plan%last_columns
    receive_values=int(sum(plan%rc)-plan%rc(r),int64)*plan%last_columns
    send_peers=count(plan%sc>0);receive_peers=count(plan%rc>0)
    if (plan%sc(r)>0) send_peers=send_peers-1
    if (plan%rc(r)>0) receive_peers=receive_peers-1
    if (present(header_bytes)) header_bytes=int(storage_size(0_int64)/8,int64)*send_peers
    ! Allocated payload only: excludes MPI communicator/request implementation,
    ! descriptors and transient compilation directory. Not a peak RSS claim.
    storage_bytes=int(storage_size(0)/8,int64)*(4_int64*plan%nrank+size(plan%source,kind=int64)+ &
         size(plan%target,kind=int64)+size(plan%address,kind=int64)+size(plan%orientation,kind=int64)+ &
         size(plan%target_peer,kind=int64)+ &
         size(plan%send_peer,kind=int64)+size(plan%receive_peer,kind=int64)) + &
         int(storage_size(0.0_dp)/8,int64)*(size(plan%send,kind=int64)+size(plan%received,kind=int64))
  end subroutine

  subroutine release_rhs_exchange(plan)
    type(Rhs_Exchange), intent(inout) :: plan
    integer :: ierr
    if (plan%comm/=MPI_COMM_NULL.and.plan%owns_comm) then
       call MPI_Comm_free(plan%comm,ierr)
       if (ierr/=MPI_SUCCESS) error stop 'RHS exchange communicator release'
    end if
    if (allocated(plan%sc)) deallocate(plan%sc,plan%sd,plan%rc,plan%rd,plan%source, &
         plan%target,plan%address,plan%orientation,plan%send,plan%received,plan%requests, &
         plan%send_peer,plan%receive_peer,plan%target_peer)
    plan%comm=MPI_COMM_NULL;plan%owns_comm=.true.;plan%tag=0;plan%variable_columns=.false.;plan%last_columns=0
    plan%topology=-1_int64;plan%masks=-1_int64;plan%last_state=-1_int64
    plan%phase=-1;plan%nslot=0;plan%nlayer=0;plan%me=-1;plan%nrank=0
  end subroutine
end module
