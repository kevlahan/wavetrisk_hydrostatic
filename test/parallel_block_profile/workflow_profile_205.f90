module workflow_profile_195_mod
  ! Stage201: PMU snapshots at common context boundaries; identical coarse contexts and selected MPI scopes; integer work counters.
  ! No batch tape, per-kernel/access clocks, barriers or timed-step I/O.
  use iso_fortran_env, only: int64,real64
  use iso_c_binding, only: c_double,c_int
  implicit none
  private
  integer,parameter :: ns=256,np=6,nt=2,ne=8192
  integer :: mode=0,n=0,phase=0,events=0,depth=0,stack(32)=0,owners(32)=0
  logical :: active=.false.,ready=.false.
  real(real64) :: epoch,last,begin_event(32),wall(ns)=0,parts(0:np-1,ns)=0,transactions(nt,0:np-1,ns)=0
  integer :: calls(nt,0:np-1,ns)=0,flags(2,ns)=0,event_info(4,ne)=0
  real(real64) :: event_time(2,ne)=0
  integer,parameter :: nr205=29
  integer :: mode205=0,depth205=0,stack205(128)=0
  real(real64) :: last205=0,seconds205(0:nr205,0:8,ns)=0
  real(real64) :: inclusive205(0:nr205,0:8,ns)=0
  integer(int64) :: calls205(0:nr205,0:8,ns)=0
  real(real64) :: mpi205(0:nr205,26:29,0:8,ns)=0
  integer,parameter :: max_life205=100000
  integer :: life_n205=0,rk205=0,event_index205(128)=0
  real(real64) :: scope_begin205(128)=0
  real(real64),allocatable :: life_time205(:,:)
  integer(int64),allocatable :: life205(:,:)
  integer(int64) :: producer205(2)=0
  public :: lifecycle_event205,lifecycle_enabled205,lifecycle_selected205
  public :: lifecycle_producer_begin205,lifecycle_producer_end205
  public :: measure_enter205,measure_leave205
  public :: wf_init,wf_begin,wf_end,wf_phase,wf_enter,wf_leave,wf_flag,wf_finish
  integer,parameter :: nr196=135,nc196=42,nv196=56
  integer :: omode=0,od=0,os(128)=0,context196=0
  logical :: oa=.false.
  real(real64) :: ot=0,osec(0:nr196,0:8,ns)=0
  integer(int64) :: ocall(0:nr196,0:8,ns)=0,ocount(nc196,0:8,ns)=0
  integer(int64) :: ovolume(2,nv196,0:8,ns)=0
  logical,public :: count199=.false.,census199=.false.
  integer(int64),public :: leaf199(42)=0_int64
  public :: sample_begin199,sample_end199,batch_begin199,batch_end199
  public :: op_enter,op_leave,op_count,op_volume,op_enabled,op_dependency,op_post,op_wait_begin,op_wait_end
  integer,parameter :: mpi199(11)=[13,18,19,34,37,38,39,40,60,72,75]
  interface
    function monotonic199() bind(C,name='workflow_clock199') result(t)
      import c_double
      real(c_double) :: t
    end function
    function affinity200() bind(C,name='workflow_affinity200') result(cpu)
      use iso_c_binding, only: c_int
      integer(c_int) :: cpu
    end function
    subroutine hw_init201() bind(C,name='workflow_hw_init201')
    end subroutine
    subroutine hw_begin201() bind(C,name='workflow_hw_begin201')
    end subroutine
    subroutine hw_tick201(step,context) bind(C,name='workflow_hw_tick201')
      import c_int
      integer(c_int),value :: step,context
    end subroutine
    subroutine hw_finish201(rank) bind(C,name='workflow_hw_finish201')
      import c_int
      integer(c_int),value :: rank
    end subroutine
  end interface
contains
  real(real64) function clock_now() result(t)
    t=monotonic199()
    if(t<0) error stop 'operations199: monotonic clock unavailable'
  end function
  subroutine wf_init
    character(16) :: value
    integer :: status,length
    integer :: expected_cpu
    if(ready) error stop 'workflow195: repeated init'
    call get_environment_variable('WAVETRISK_WORKFLOW195',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'workflow195: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'012')/=0) error stop 'workflow195: invalid selector'
      read(value,*) mode
    end if
    call get_environment_variable('WAVETRISK_SCALING200',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations199: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'01')/=0) error stop 'operations199: invalid selector'
      read(value,*) omode
    end if
    if(mode/=1) error stop 'operations199: common coarse observer required'
    call get_environment_variable('WAVETRISK_PROFILE_BLOCK_DETAIL',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations199: old detail profiler must be off'
    call get_environment_variable('WAVETRISK_PROFILE_PARALLEL_BLOCKS',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations199: old block profiler must be off'
    call get_environment_variable('WAVETRISK_EXPECT_CPU200',value,length=length,status=status)
    if(status==0) then
      read(value,*) expected_cpu
      if(affinity200()/=expected_cpu) error stop 'Stage200 MPI initialization changed rank CPU affinity'
    end if
    if(omode/=0) error stop 'Stage201 uses coarse contexts only'
    call get_environment_variable('WAVETRISK_LIFECYCLE205',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'Stage205 malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'01')/=0) error stop 'Stage205 invalid selector'
      read(value,*) mode205
    end if
    if(mode205==1) allocate(life205(10,max_life205),life_time205(2,max_life205))
    call hw_init201()
    ready=.true.
  end subroutine
  subroutine wf_begin
    if(.not.ready.or.active.or.n>=ns) error stop 'workflow195: invalid step begin/overflow'
    call hw_begin201()
    n=n+1;phase=0;active=.true.;epoch=clock_now();last=epoch
    last205=epoch;rk205=0
    context196=5;ot=epoch;oa=.true.;count199=omode==1;census199=.false.
  end subroutine
  subroutine wf_phase(p)
    integer,intent(in) :: p
    real(real64) :: t
    if(.not.active.or.mode==0) return
    if(p/=phase+1.or.p>=np.or.depth/=0) error stop 'workflow195: phase order or nested transaction'
    t=clock_now();parts(phase,n)=parts(phase,n)+t-last
    call event(0,phase,last,t)
    if(od/=0) error stop 'operations199: scope crosses phase boundary'
    call flush199
    call op_tick(t)
    phase=p;last=t
    select case(p)
    case(1)
      context196=0
    case(2)
      context196=6
    case(3)
      context196=4
    case(4)
      context196=7
    case(5)
      context196=8
    end select
  end subroutine
  subroutine wf_enter(tag)
    integer,intent(in) :: tag
    if(.not.active.or.mode==0) return
    if(tag<1.or.tag>nt.or.depth>=size(stack)) error stop 'workflow195: bad transaction'
    depth=depth+1;stack(depth)=tag;owners(depth)=phase;begin_event(depth)=clock_now()
    if(oa) then
      if(depth/=1) error stop 'operations199: nested workflow transaction'
      call flush199
      call op_tick(begin_event(depth))
      if(phase==1) then
        context196=tag
      else if(phase==3.and.tag==1) then
        context196=3
      end if
    end if
  end subroutine
  subroutine wf_leave(tag)
    integer,intent(in) :: tag
    real(real64) :: t
    if(.not.active.or.mode==0) return
    if(depth<1) error stop 'workflow195: empty transaction'
    if(stack(depth)/=tag.or.owners(depth)/=phase) error stop 'workflow195: mismatched transaction'
    t=clock_now()
    if(oa) then
      call flush199
      call op_tick(t)
      if(phase==1) context196=0
      if(phase==3) context196=4
    end if
    transactions(tag,phase,n)=transactions(tag,phase,n)+t-begin_event(depth)
    calls(tag,phase,n)=calls(tag,phase,n)+1
    call event(tag,phase,begin_event(depth),t)
    depth=depth-1
  end subroutine
  subroutine wf_flag(tag)
    integer,intent(in) :: tag
    if(.not.active) return
    if(tag<1.or.tag>2) error stop 'workflow195: invalid flag'
    flags(tag,n)=flags(tag,n)+1
  end subroutine
  subroutine event(tag,p,a,b)
    integer,intent(in) :: tag,p
    real(real64),intent(in) :: a,b
    if(mode<2) return
    if(events>=ne) error stop 'workflow195: event overflow'
    events=events+1;event_info(:,events)=[n,tag,p,events]
    event_time(:,events)=[a-epoch,b-epoch]
  end subroutine
  subroutine wf_end
    real(real64) :: t
    if(.not.active.or.depth/=0) error stop 'workflow195: invalid end'
    t=clock_now();wall(n)=t-epoch
    if(mode>0) then
      if(phase/=np-1) error stop 'workflow195: missing phase'
      parts(phase,n)=parts(phase,n)+t-last;call event(0,phase,last,t)
    else
      parts(0,n)=wall(n)
    end if
    if(depth205/=0) error stop 'Stage205 unfinished scope'
    if(od/=0) error stop 'operations199: unfinished operation at step end'
    call flush199
    call op_tick(t)
    oa=.false.;count199=.false.;census199=.false.;active=.false.
  end subroutine
  subroutine wf_finish(rank)
    integer,intent(in) :: rank
    integer :: u,s,p,t,e
    character(80) :: filename
    if(active.or.depth/=0.or..not.ready) error stop 'workflow195: unfinished'
    write(filename,'(a,i0,a)') 'workflow195-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',1,mode,rank
    do s=1,n
      write(u,'(a,3(1x,i0),1x,es24.16)') 'step',s,flags(:,s),wall(s)
      do p=0,np-1
        write(u,'(a,2(1x,i0),1x,es24.16)') 'phase',s,p,parts(p,s)
        do t=1,nt
          write(u,'(a,4(1x,i0),1x,es24.16)') 'transaction',s,p,t,calls(t,p,s),transactions(t,p,s)
        end do
      end do
    end do
    do e=1,events
      write(u,'(a,4(1x,i0),2(1x,es24.16))') 'event',event_info(:,e),event_time(:,e)
    end do
    close(u)
    call hw_finish201(int(rank,c_int))
    call op_finish(rank)
    call measure_finish205(rank)
    call lifecycle_finish205(rank)
  end subroutine

  logical function op_selected(p) result(yes)
    integer,intent(in) :: p
    yes=omode==1.and.any(mpi199==p)
  end function

  subroutine op_tick(t)
    real(real64),intent(in) :: t
    integer :: id
    if(.not.oa) return
    call hw_tick201(int(n,c_int),int(context196,c_int))
    call measure_tick205(t)
    if(t<ot) then
      write(*,*) 'operations199 clock reversal: step/context/depth/current/previous/delta',n,context196,od,t,ot,t-ot
      error stop 'operations199: nonmonotonic clock'
    end if
    id=0
    if(od>0) id=os(od)
    osec(id,context196,n)=osec(id,context196,n)+t-ot
    ot=t
  end subroutine

  subroutine op_enter(p)
    integer,intent(in) :: p
    if(.not.oa) return
    if(.not.op_selected(p)) return
    if(od>=size(os)) error stop 'operations199: stack overflow'
    call op_tick(clock_now())
    od=od+1;os(od)=p
    ocall(p,context196,n)=ocall(p,context196,n)+1_int64
  end subroutine

  subroutine op_leave(p)
    integer,intent(in) :: p
    if(.not.oa) return
    if(.not.op_selected(p)) return
    if(od<1) error stop 'operations199: empty stack'
    if(os(od)/=p) then
      write(*,*) 'operations199: mismatched scopes',os(od),p
      error stop 'operations199: scope mismatch'
    end if
    call op_tick(clock_now());od=od-1
  end subroutine

  subroutine op_count(p,value)
    integer,intent(in) :: p
    integer(int64),intent(in) :: value
    if(.not.count199) return
    if(p<1.or.p>nc196.or.value<0) error stop 'operations199: invalid counter'
    ocount(p,context196,n)=ocount(p,context196,n)+value
  end subroutine

  subroutine op_volume(p,messages,bytes)
    integer,intent(in) :: p
    integer(int64),intent(in) :: messages,bytes
    if(.not.count199) return
    if(p<1.or.p>nv196.or.messages<0.or.bytes<0) error stop 'operations199: invalid volume'
    ovolume(:,p,context196,n)=ovolume(:,p,context196,n)+[messages,bytes]
  end subroutine

  subroutine op_finish(rank)
    integer,intent(in) :: rank
    integer :: u,s,c,p
    character(80) :: filename
    if(oa.or.od/=0) error stop 'operations199: incomplete scopes'
    write(filename,'(a,i0,a)') 'scaling200-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',1,omode,rank
    do s=1,n
      write(u,'(a,3(1x,i0),1x,es24.16)') 'step',s,flags(:,s),wall(s)
      do c=0,8
        do p=0,nr196
          if(p/=0.and..not.any(mpi199==p)) cycle
          write(u,'(a,3(1x,i0),1x,i0,1x,es24.16)') 'self',s,c,p,ocall(p,c,s),osec(p,c,s)
        end do
        do p=1,nc196
          write(u,'(a,3(1x,i0),1x,i0)') 'count',s,c,p,ocount(p,c,s)
        end do
        do p=1,nv196
          if(any(ovolume(:,p,c,s)/=0_int64)) &
            write(u,'(a,3(1x,i0),2(1x,i0))') 'volume',s,c,p,ovolume(:,p,c,s)
        end do
      end do
    end do
    close(u)
  end subroutine
  logical function op_enabled() result(enabled)
    enabled=active.and.omode==1
  end function
  subroutine flush199
    if(count199) ocount(:,context196,n)=ocount(:,context196,n)+leaf199
    leaf199=0_int64
  end subroutine
  subroutine op_dependency(kind,component,family,first,last)
    integer,intent(in)::kind,component,family
    integer,optional,intent(in)::first,last
  end subroutine
  subroutine op_post(requests,sent,received,tag)
    integer,intent(in)::requests,sent,received,tag
  end subroutine
  subroutine op_wait_begin
  end subroutine
  subroutine op_wait_end
  end subroutine
  subroutine batch_begin199(id,k,l,pos,work)
    integer,intent(in)::id
    integer,optional,intent(in)::k,l,pos,work
  end subroutine
  subroutine batch_end199(id)
    integer,intent(in)::id
  end subroutine
  logical function sample_begin199(id) result(take)
    integer,intent(in)::id
    take=.false.
  end function
  subroutine sample_end199(id,take)
    integer,intent(in)::id
    logical,intent(in)::take
    if(take) error stop 'Stage200 forbids per-access clocks'
  end subroutine
  subroutine measure_tick205(t)
    real(real64),intent(in) :: t
    integer :: id,parent,j
    if(.not.oa) return
    if(t<last205) error stop 'Stage205 clock reversal'
    id=0
    if(depth205>0) id=stack205(depth205)
    seconds205(id,context196,n)=seconds205(id,context196,n)+t-last205
    do j=1,depth205
      if(any(stack205(:j-1)==stack205(j))) cycle
      inclusive205(stack205(j),context196,n)=inclusive205(stack205(j),context196,n)+t-last205
    end do
    if(id>=26.and.id<=29) then
      parent=0
      if(depth205>1) parent=stack205(depth205-1)
      mpi205(parent,id,context196,n)=mpi205(parent,id,context196,n)+t-last205
    end if
    last205=t
  end subroutine
  subroutine measure_enter205(id)
    integer,intent(in) :: id
    if(.not.oa.or.mode205==0) return
    if(id>=26.and.depth205==0) return
    if(id<1.or.id>nr205.or.depth205>=size(stack205)) error stop 'Stage205 invalid scope'
    call measure_tick205(clock_now())
    depth205=depth205+1;stack205(depth205)=id
    event_index205(depth205)=0;scope_begin205(depth205)=last205
    calls205(id,context196,n)=calls205(id,context196,n)+1_int64
  end subroutine
  subroutine measure_leave205(id)
    integer,intent(in) :: id
    real(real64) :: now
    if(.not.oa.or.mode205==0) return
    if(id>=26.and.depth205==0) return
    if(depth205<1) error stop 'Stage205 empty stack'
    if(stack205(depth205)/=id) error stop 'Stage205 mismatched scope'
    now=clock_now()
    call measure_tick205(now)
    if(event_index205(depth205)>0) &
         life_time205(2,event_index205(depth205))=now-scope_begin205(depth205)
    depth205=depth205-1
  end subroutine
  subroutine measure_finish205(rank)
    integer,intent(in) :: rank
    integer :: u,s,c,p,parent
    character(80) :: filename
    write(filename,'(a,i0,a)') 'operations205-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',1,mode205,rank
    do s=1,n
      write(u,'(a,3(1x,i0),1x,es24.16)') 'step',s,flags(:,s),wall(s)
      do c=0,8
        do p=0,nr205
          write(u,'(a,4(1x,i0),2(1x,es24.16))') 'self',s,c,p,calls205(p,c,s), &
               seconds205(p,c,s),inclusive205(p,c,s)
        end do
        do p=26,29
          do parent=0,nr205
            if(mpi205(parent,p,c,s)==0) cycle
            write(u,'(a,4(1x,i0),1x,es24.16)') 'mpi',s,c,p,parent,mpi205(parent,p,c,s)
          end do
        end do
      end do
    end do
    close(u)
  end subroutine
  logical function lifecycle_enabled205() result(yes)
    yes=oa.and.mode205==1
  end function
  logical function lifecycle_selected205() result(yes)
    yes=oa.and.mode205==1.and.depth205>0
  end function
  subroutine lifecycle_event205(kind,object,generation,a,b,c)
    integer,intent(in) :: kind,object
    integer(int64),intent(in) :: generation,a,b,c
    integer :: owner
    if(.not.lifecycle_enabled205()) return
    if(life_n205>=max_life205) error stop 'Stage205 lifecycle event overflow'
    owner=0
    if(depth205>0) owner=stack205(depth205)
    if(kind==40.and.depth205>1) owner=stack205(depth205-1)
    if(kind==8.or.kind==9.or.kind==15) rk205=int(a)
    life_n205=life_n205+1
    life_time205(:,life_n205)=[-1.0_real64,0.0_real64]
    if(kind==40) then
      if(depth205==0) error stop 'Stage205 collective outside selected scope'
      event_index205(depth205)=life_n205
      life_time205(1,life_n205)=scope_begin205(depth205)-epoch
    end if
    life205(:,life_n205)=[int(n,int64),int(context196,int64),int(owner,int64),int(rk205,int64), &
         int(kind,int64),int(object,int64),generation,a,b,c]
  end subroutine
  subroutine lifecycle_producer_begin205(values,geometry)
    integer(int64),intent(in) :: values,geometry
    if(.not.lifecycle_enabled205()) return
    producer205=[values,geometry]
  end subroutine
  subroutine lifecycle_producer_end205(generation,values,geometry)
    integer(int64),intent(in) :: generation,values,geometry
    if(.not.lifecycle_enabled205()) return
    if(values<producer205(1).or.geometry<producer205(2)) error stop 'Stage205 producer counter reversed'
    call lifecycle_event205(13,1,generation,values-producer205(1),geometry-producer205(2),0_int64)
  end subroutine
  subroutine lifecycle_finish205(rank)
    integer,intent(in) :: rank
    integer :: u,i
    character(80) :: filename
    write(filename,'(a,i0,a)') 'lifecycle205-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,4(1x,i0))') 'schema',1,mode205,rank,life_n205
    do i=1,life_n205
      write(u,'(a,11(1x,i0),2(1x,es24.16))') 'event',i,life205(:,i),life_time205(:,i)
    end do
    close(u)
  end subroutine
end module
