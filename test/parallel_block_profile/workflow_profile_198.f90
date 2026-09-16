module workflow_profile_195_mod
  ! Identical observer in both implementations. No MPI or I/O in timed steps.
  use iso_fortran_env, only: int64,real64
  use iso_c_binding, only: c_double
  implicit none
  private
  integer,parameter :: ns=256,np=6,nt=2,ne=8192
  integer :: mode=0,n=0,phase=0,events=0,depth=0,stack(32)=0,owners(32)=0
  logical :: active=.false.,ready=.false.
  real(real64) :: epoch,last,begin_event(32),wall(ns)=0,parts(0:np-1,ns)=0,transactions(nt,0:np-1,ns)=0
  integer :: calls(nt,0:np-1,ns)=0,flags(2,ns)=0,event_info(4,ne)=0
  real(real64) :: event_time(2,ne)=0
  public :: wf_init,wf_begin,wf_end,wf_phase,wf_enter,wf_leave,wf_flag,wf_finish
  integer,parameter :: nr196=100,nc196=42,nv196=56
  integer :: omode=0,od=0,os(128)=0,context196=0
  logical :: oa=.false.
  real(real64) :: ot=0,osec(0:nr196,0:8,ns)=0
  integer(int64) :: ocall(0:nr196,0:8,ns)=0,ocount(nc196,0:8,ns)=0
  integer(int64) :: ovolume(2,nv196,0:8,ns)=0
  integer,parameter :: nd197=16384
  integer :: depn=0,depkey(5)=0,depinfo(12,nd197)=0
  logical :: post_pending=.false.,wait_pending=.false.
  real(real64) :: post_at=0,wait_at=0,deptimes(2,nd197)=0
  logical,public :: count198=.false.
  integer(int64),public :: leaf198(42)=0_int64
  integer,parameter :: nprobe=4,stride198=257
  integer :: offset198=0
  integer(int64) :: pcalls(nprobe,0:8,ns)=0,pns(nprobe,0:8,ns)=0
  real(real64) :: ptime(nprobe,0:8,ns)=0,psquare(nprobe,0:8,ns)=0,pstart(nprobe)=0
  public :: sample_begin198,sample_end198
  public :: op_enter,op_leave,op_count,op_volume,op_enabled,op_dependency,op_post,op_wait_begin,op_wait_end
  interface
    function monotonic198() bind(C,name='workflow_clock198') result(t)
      import c_double
      real(c_double) :: t
    end function
  end interface
contains
  real(real64) function clock_now() result(t)
    t=monotonic198()
    if(t<0) error stop 'operations198: monotonic clock unavailable'
  end function
  subroutine wf_init
    character(16) :: value
    integer :: status,length
    if(ready) error stop 'workflow195: repeated init'
    call get_environment_variable('WAVETRISK_WORKFLOW195',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'workflow195: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'012')/=0) error stop 'workflow195: invalid selector'
      read(value,*) mode
    end if
    call get_environment_variable('WAVETRISK_OPERATIONS198',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations198: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'0123')/=0) error stop 'operations198: invalid selector'
      read(value,*) omode
    end if
    call get_environment_variable('WAVETRISK_SAMPLE_OFFSET198',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations198: sample offset unreadable'
    if(status==0) then
      if(length<1.or.verify(trim(value),'0123456789')/=0) error stop 'operations198: bad offset'
      read(value,*) offset198
      if(offset198<0.or.offset198>=stride198) error stop 'operations198: offset out of range'
    end if
    if(mode/=1) error stop 'operations198: common coarse observer required'
    call get_environment_variable('WAVETRISK_PROFILE_BLOCK_DETAIL',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations198: old detail profiler must be off'
    call get_environment_variable('WAVETRISK_PROFILE_PARALLEL_BLOCKS',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations198: old block profiler must be off'
    ready=.true.
  end subroutine
  subroutine wf_begin
    if(.not.ready.or.active.or.n>=ns) error stop 'workflow195: invalid step begin/overflow'
    n=n+1;phase=0;active=.true.;epoch=clock_now();last=epoch
    context196=5;ot=epoch;oa=omode>0;count198=omode>=2
  end subroutine
  subroutine wf_phase(p)
    integer,intent(in) :: p
    real(real64) :: t
    if(.not.active.or.mode==0) return
    if(p/=phase+1.or.p>=np.or.depth/=0) error stop 'workflow195: phase order or nested transaction'
    t=clock_now();parts(phase,n)=parts(phase,n)+t-last
    call event(0,phase,last,t)
    if(od/=0) error stop 'operations198: scope crosses phase boundary'
    call flush198
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
      if(depth/=1) error stop 'operations198: nested workflow transaction'
      call flush198
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
      call flush198
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
    if(od/=0) error stop 'operations198: unfinished operation at step end'
    call flush198
    call op_tick(t)
    oa=.false.;count198=.false.;active=.false.
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
    call op_finish(rank)
  end subroutine

  logical function op_selected(p) result(yes)
    integer,intent(in) :: p
    yes=omode>0.and.p>=2.and.p<=nr196
  end function

  subroutine op_tick(t)
    real(real64),intent(in) :: t
    integer :: id
    if(.not.oa) return
    if(t<ot) then
      write(*,*) 'operations198 clock reversal: step/context/depth/current/previous/delta',n,context196,od,t,ot,t-ot
      error stop 'operations198: nonmonotonic clock'
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
    if(od>=size(os)) error stop 'operations198: stack overflow'
    call op_tick(clock_now())
    od=od+1;os(od)=p
    ocall(p,context196,n)=ocall(p,context196,n)+1_int64
  end subroutine

  subroutine op_leave(p)
    integer,intent(in) :: p
    if(.not.oa) return
    if(.not.op_selected(p)) return
    if(od<1) error stop 'operations198: empty stack'
    if(os(od)/=p) then
      write(*,*) 'operations198: mismatched scopes',os(od),p
      error stop 'operations198: scope mismatch'
    end if
    call op_tick(clock_now());od=od-1
  end subroutine

  subroutine op_count(p,value)
    integer,intent(in) :: p
    integer(int64),intent(in) :: value
    if(.not.count198) return
    if(p<1.or.p>nc196.or.value<0) error stop 'operations198: invalid counter'
    ocount(p,context196,n)=ocount(p,context196,n)+value
  end subroutine

  subroutine op_volume(p,messages,bytes)
    integer,intent(in) :: p
    integer(int64),intent(in) :: messages,bytes
    if(.not.oa) return
    if(p<1.or.p>nv196.or.messages<0.or.bytes<0) error stop 'operations198: invalid volume'
    ovolume(:,p,context196,n)=ovolume(:,p,context196,n)+[messages,bytes]
  end subroutine

  subroutine op_finish(rank)
    integer,intent(in) :: rank
    integer :: u,s,c,p
    character(80) :: filename
    if(post_pending.or.wait_pending) error stop 'operations198: incomplete dependency'
    if(oa.or.od/=0) error stop 'operations198: incomplete scopes'
    write(filename,'(a,i0,a)') 'operations198-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',2,omode,rank
    do s=1,n
      write(u,'(a,3(1x,i0),1x,es24.16)') 'step',s,flags(:,s),wall(s)
      do c=0,8
        do p=0,nr196
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
    write(filename,'(a,i0,a)') 'storage198-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,4(1x,i0))') 'schema',1,rank,stride198,offset198
    do s=1,n
      do c=0,8
        do p=1,nprobe
          write(u,'(a,3(1x,i0),2(1x,i0),2(1x,es24.16))') &
            'sample',s,c,p,pcalls(p,c,s),pns(p,c,s),ptime(p,c,s),psquare(p,c,s)
        end do
      end do
    end do
    close(u)
    write(filename,'(a,i0,a)') 'dependencies198-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,2(1x,i0))') 'schema',1,rank
    do s=1,depn
      write(u,'(12(i0,1x),2(es24.16,1x))') depinfo(:,s),deptimes(:,s)
    end do
    close(u)
  end subroutine
  logical function op_enabled() result(enabled)
    enabled=oa
  end function
  subroutine op_dependency(kind,component,family,first,last)
    integer,intent(in)::kind,component,family
    integer,optional,intent(in)::first,last
    if(.not.oa.or.omode<2) return
    if(post_pending) error stop 'operations198: dependency changed while pending'
    depkey=[kind,component,family,-999,-999]
    if(present(first)) depkey(4)=first
    if(present(last)) depkey(5)=last
  end subroutine
  subroutine op_post(requests,sent,received,tag)
    integer,intent(in)::requests,sent,received,tag
    if(.not.oa.or.omode<2) return
    if(post_pending.or.wait_pending.or.depn>=nd197) error stop 'operations198: post state/overflow'
    depn=depn+1;post_pending=.true.;post_at=clock_now()
    depinfo(:,depn)=[n,context196,depkey,requests,sent,received,tag,depn]
  end subroutine
  subroutine op_wait_begin
    if(.not.oa.or.omode<2) return
    if(.not.post_pending.or.wait_pending) error stop 'operations198: unmatched wait begin'
    wait_pending=.true.;wait_at=clock_now()
    deptimes(1,depn)=wait_at-post_at
  end subroutine
  subroutine op_wait_end
    if(.not.oa.or.omode<2) return
    if(.not.wait_pending) error stop 'operations198: unmatched wait end'
    deptimes(2,depn)=clock_now()-wait_at
    post_pending=.false.;wait_pending=.false.
  end subroutine
  subroutine flush198
    if(count198) ocount(:,context196,n)=ocount(:,context196,n)+leaf198
    leaf198=0_int64
  end subroutine
  logical function sample_begin198(id) result(take)
    integer,intent(in)::id
    take=.false.
    if(.not.oa.or.omode/=3) return
    if(id<1.or.id>nprobe) error stop 'operations198: invalid sample id'
    pcalls(id,context196,n)=pcalls(id,context196,n)+1_int64
    take=mod(pcalls(id,context196,n)+offset198,int(stride198,int64))==0
    if(take) pstart(id)=clock_now()
  end function
  subroutine sample_end198(id,take)
    integer,intent(in)::id
    logical,intent(in)::take
    real(real64)::dt
    if(.not.take) return
    dt=clock_now()-pstart(id)
    if(dt<0) error stop 'operations198: invalid sample interval'
    pns(id,context196,n)=pns(id,context196,n)+1_int64
    ptime(id,context196,n)=ptime(id,context196,n)+dt
    psquare(id,context196,n)=psquare(id,context196,n)+dt*dt
  end subroutine
end module
