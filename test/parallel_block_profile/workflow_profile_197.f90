module workflow_profile_195_mod
  ! Identical observer in both implementations. No MPI or I/O in timed steps.
  use iso_fortran_env, only: int64,real64
  implicit none
  private
  integer,parameter :: ns=256,np=6,nt=2,ne=8192
  integer :: mode=0,n=0,phase=0,events=0,depth=0,stack(32)=0,owners(32)=0
  logical :: active=.false.,ready=.false.
  real(real64) :: epoch,last,begin_event(32),wall(ns)=0,parts(0:np-1,ns)=0,transactions(nt,0:np-1,ns)=0
  integer :: calls(nt,0:np-1,ns)=0,flags(2,ns)=0,event_info(4,ne)=0
  real(real64) :: event_time(2,ne)=0
  public :: wf_init,wf_begin,wf_end,wf_phase,wf_enter,wf_leave,wf_flag,wf_finish
  integer,parameter :: nr196=82,nc196=35,nv196=56
  integer :: omode=0,od=0,os(128)=0,context196=0
  logical :: oa=.false.
  real(real64) :: ot=0,osec(0:nr196,0:2,ns)=0
  integer(int64) :: ocall(0:nr196,0:2,ns)=0,ocount(nc196,0:2,ns)=0
  integer(int64) :: ovolume(2,nv196,0:2,ns)=0
  integer,parameter :: nd197=16384
  integer :: depn=0,depkey(5)=0,depinfo(12,nd197)=0
  logical :: post_pending=.false.,wait_pending=.false.
  real(real64) :: post_at=0,wait_at=0,deptimes(2,nd197)=0
  public :: op_enter,op_leave,op_count,op_volume,op_enabled,op_dependency,op_post,op_wait_begin,op_wait_end
contains
  real(real64) function clock_now() result(t)
    integer(int64) :: c,r
    call system_clock(c,r)
    if(r<=0) error stop 'workflow195: no clock'
    t=real(c,real64)/real(r,real64)
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
    call get_environment_variable('WAVETRISK_OPERATIONS197',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations197: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'012')/=0) error stop 'operations197: invalid selector'
      read(value,*) omode
    end if
    if(mode/=1) error stop 'operations197: common coarse observer required'
    call get_environment_variable('WAVETRISK_PROFILE_BLOCK_DETAIL',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations197: old detail profiler must be off'
    call get_environment_variable('WAVETRISK_PROFILE_PARALLEL_BLOCKS',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations197: old block profiler must be off'
    ready=.true.
  end subroutine
  subroutine wf_begin
    if(.not.ready.or.active.or.n>=ns) error stop 'workflow195: invalid step begin/overflow'
    n=n+1;phase=0;active=.true.;epoch=clock_now();last=epoch
  end subroutine
  subroutine wf_phase(p)
    integer,intent(in) :: p
    real(real64) :: t
    if(.not.active.or.mode==0) return
    if(p/=phase+1.or.p>=np.or.depth/=0) error stop 'workflow195: phase order or nested transaction'
    t=clock_now();parts(phase,n)=parts(phase,n)+t-last
    call event(0,phase,last,t)
    if(phase==1) then
      if(od/=0) error stop 'operations197: open scope at dynamics end'
      call op_tick(t)
      oa=.false.
    end if
    phase=p;last=t
    if(p==1) then
      oa=omode>0;ot=t;context196=0
    end if
  end subroutine
  subroutine wf_enter(tag)
    integer,intent(in) :: tag
    if(.not.active.or.mode==0) return
    if(tag<1.or.tag>nt.or.depth>=size(stack)) error stop 'workflow195: bad transaction'
    depth=depth+1;stack(depth)=tag;owners(depth)=phase;begin_event(depth)=clock_now()
    if(oa) then
      if(depth/=1) error stop 'operations197: nested workflow transaction'
      call op_tick(begin_event(depth));context196=tag
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
      call op_tick(t);context196=0
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
    active=.false.
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
    yes=.false.
    if(omode==2) then
      yes=p>=2.and.p<=nr196
    else if(omode==1) then
      select case(p)
      case(15:20,24,28,32:34,37:40,58:63,68:77)
        yes=.true.
      end select
    end if
  end function

  subroutine op_tick(t)
    real(real64),intent(in) :: t
    integer :: id
    if(.not.oa) return
    if(t<ot) error stop 'operations197: nonmonotonic clock'
    id=0
    if(od>0) id=os(od)
    osec(id,context196,n)=osec(id,context196,n)+t-ot
    ot=t
  end subroutine

  subroutine op_enter(p)
    integer,intent(in) :: p
    if(.not.oa) return
    if(.not.op_selected(p)) return
    if(od>=size(os)) error stop 'operations197: stack overflow'
    call op_tick(clock_now())
    od=od+1;os(od)=p
    ocall(p,context196,n)=ocall(p,context196,n)+1_int64
  end subroutine

  subroutine op_leave(p)
    integer,intent(in) :: p
    if(.not.oa) return
    if(.not.op_selected(p)) return
    if(od<1) error stop 'operations197: empty stack'
    if(os(od)/=p) then
      write(*,*) 'operations197: mismatched scopes',os(od),p
      error stop 'operations197: scope mismatch'
    end if
    call op_tick(clock_now());od=od-1
  end subroutine

  subroutine op_count(p,value)
    integer,intent(in) :: p
    integer(int64),intent(in) :: value
    if(.not.oa) return
    if(p<1.or.p>nc196.or.value<0) error stop 'operations197: invalid counter'
    ocount(p,context196,n)=ocount(p,context196,n)+value
  end subroutine

  subroutine op_volume(p,messages,bytes)
    integer,intent(in) :: p
    integer(int64),intent(in) :: messages,bytes
    if(.not.oa) return
    if(p<1.or.p>nv196.or.messages<0.or.bytes<0) error stop 'operations197: invalid volume'
    ovolume(:,p,context196,n)=ovolume(:,p,context196,n)+[messages,bytes]
  end subroutine

  subroutine op_finish(rank)
    integer,intent(in) :: rank
    integer :: u,s,c,p
    character(80) :: filename
    if(post_pending.or.wait_pending) error stop 'operations197: incomplete dependency'
    if(oa.or.od/=0) error stop 'operations197: incomplete scopes'
    write(filename,'(a,i0,a)') 'operations197-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',1,omode,rank
    do s=1,n
      write(u,'(a,3(1x,i0),1x,es24.16)') 'step',s,flags(:,s),parts(1,s)
      do c=0,2
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
    write(filename,'(a,i0,a)') 'dependencies197-rank-',rank,'.txt'
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
    if(.not.oa) return
    if(post_pending) error stop 'operations197: dependency changed while pending'
    depkey=[kind,component,family,-999,-999]
    if(present(first)) depkey(4)=first
    if(present(last)) depkey(5)=last
  end subroutine
  subroutine op_post(requests,sent,received,tag)
    integer,intent(in)::requests,sent,received,tag
    if(.not.oa) return
    if(post_pending.or.wait_pending.or.depn>=nd197) error stop 'operations197: post state/overflow'
    depn=depn+1;post_pending=.true.;post_at=clock_now()
    depinfo(:,depn)=[n,context196,depkey,requests,sent,received,tag,depn]
  end subroutine
  subroutine op_wait_begin
    if(.not.oa) return
    if(.not.post_pending.or.wait_pending) error stop 'operations197: unmatched wait begin'
    wait_pending=.true.;wait_at=clock_now()
    deptimes(1,depn)=wait_at-post_at
  end subroutine
  subroutine op_wait_end
    if(.not.oa) return
    if(.not.wait_pending) error stop 'operations197: unmatched wait end'
    deptimes(2,depn)=clock_now()-wait_at
    post_pending=.false.;wait_pending=.false.
  end subroutine
end module
