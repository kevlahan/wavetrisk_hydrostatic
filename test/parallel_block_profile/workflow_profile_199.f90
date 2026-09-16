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
  integer,parameter :: nr196=135,nc196=42,nv196=56
  integer :: omode=0,od=0,os(128)=0,context196=0
  logical :: oa=.false.
  real(real64) :: ot=0,osec(0:nr196,0:8,ns)=0
  integer(int64) :: ocall(0:nr196,0:8,ns)=0,ocount(nc196,0:8,ns)=0
  integer(int64) :: ovolume(2,nv196,0:8,ns)=0
  integer,parameter :: nd197=16384
  integer :: depn=0,depkey(5)=0,depinfo(12,nd197)=0
  logical :: post_pending=.false.,wait_pending=.false.
  real(real64) :: post_at=0,wait_at=0,deptimes(2,nd197)=0
  logical,public :: count199=.false.,census199=.false.
  integer(int64),public :: leaf199(42)=0_int64
  integer,parameter :: nprobe=4,stride199=257
  integer :: offset199=0
  integer(int64) :: pcalls(nprobe,0:8,ns)=0,pns(nprobe,0:8,ns)=0
  real(real64) :: ptime(nprobe,0:8,ns)=0,psquare(nprobe,0:8,ns)=0
  public :: sample_begin199,sample_end199, batch_begin199,batch_end199
  public :: op_enter,op_leave,op_count,op_volume,op_enabled,op_dependency,op_post,op_wait_begin,op_wait_end
  integer,parameter :: maxbatch199=262144
  integer :: nb199=0,bd199=0,bs199(64)=0
  integer :: bi199(9,maxbatch199)=0
  integer(int64) :: bw199(maxbatch199)=0,bc199(7,maxbatch199)=0,bl199(7,64)=0
  real(real64) :: bt199(3,maxbatch199)=0,bm199(64)=0
  integer,parameter :: mpi199(11)=[13,18,19,34,37,38,39,40,60,72,75]
  interface
    function monotonic199() bind(C,name='workflow_clock199') result(t)
      import c_double
      real(c_double) :: t
    end function
  end interface
contains
  real(real64) function clock_now() result(t)
    t=monotonic199()
    if(t<0) error stop 'operations199: monotonic clock unavailable'
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
    call get_environment_variable('WAVETRISK_OPERATIONS199',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations199: malformed selector'
    if(status==0) then
      if(length/=1.or.verify(trim(value),'0123')/=0) error stop 'operations199: invalid selector'
      read(value,*) omode
    end if
    call get_environment_variable('WAVETRISK_SAMPLE_OFFSET199',value,length=length,status=status)
    if(status/=0.and.status/=1) error stop 'operations199: sample offset unreadable'
    if(status==0) then
      if(length<1.or.verify(trim(value),'0123456789')/=0) error stop 'operations199: bad offset'
      read(value,*) offset199
      if(offset199<0.or.offset199>=stride199) error stop 'operations199: offset out of range'
    end if
    if(mode/=1) error stop 'operations199: common coarse observer required'
    call get_environment_variable('WAVETRISK_PROFILE_BLOCK_DETAIL',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations199: old detail profiler must be off'
    call get_environment_variable('WAVETRISK_PROFILE_PARALLEL_BLOCKS',value,status=status)
    if(status==0.and.trim(value)/='0') error stop 'operations199: old block profiler must be off'
    ready=.true.
  end subroutine
  subroutine wf_begin
    if(.not.ready.or.active.or.n>=ns) error stop 'workflow195: invalid step begin/overflow'
    n=n+1;phase=0;active=.true.;epoch=clock_now();last=epoch
    context196=5;ot=epoch;oa=omode>0;count199=omode>=2;census199=omode==3
  end subroutine
  subroutine wf_phase(p)
    integer,intent(in) :: p
    real(real64) :: t
    if(.not.active.or.mode==0) return
    if(p/=phase+1.or.p>=np.or.depth/=0) error stop 'workflow195: phase order or nested transaction'
    t=clock_now();parts(phase,n)=parts(phase,n)+t-last
    call event(0,phase,last,t)
    if(bd199/=0) error stop 'stage199: unfinished batch'
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
    call batch_finish199(rank)
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
    if(.not.oa) return
    if(p<1.or.p>nv196.or.messages<0.or.bytes<0) error stop 'operations199: invalid volume'
    ovolume(:,p,context196,n)=ovolume(:,p,context196,n)+[messages,bytes]
  end subroutine

  subroutine op_finish(rank)
    integer,intent(in) :: rank
    integer :: u,s,c,p
    character(80) :: filename
    if(post_pending.or.wait_pending) error stop 'operations199: incomplete dependency'
    if(oa.or.od/=0) error stop 'operations199: incomplete scopes'
    write(filename,'(a,i0,a)') 'operations199-rank-',rank,'.txt'
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
    write(filename,'(a,i0,a)') 'storage199-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,4(1x,i0))') 'schema',1,rank,stride199,offset199
    do s=1,n
      do c=0,8
        do p=1,nprobe
          write(u,'(a,3(1x,i0),2(1x,i0),2(1x,es24.16))') &
            'sample',s,c,p,pcalls(p,c,s),pns(p,c,s),ptime(p,c,s),psquare(p,c,s)
        end do
      end do
    end do
    close(u)
    write(filename,'(a,i0,a)') 'dependencies199-rank-',rank,'.txt'
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
    if(post_pending) error stop 'operations199: dependency changed while pending'
    depkey=[kind,component,family,-999,-999]
    if(present(first)) depkey(4)=first
    if(present(last)) depkey(5)=last
  end subroutine
  subroutine op_post(requests,sent,received,tag)
    integer,intent(in)::requests,sent,received,tag
    if(.not.oa.or.omode<2) return
    if(post_pending.or.wait_pending.or.depn>=nd197) error stop 'operations199: post state/overflow'
    depn=depn+1;post_pending=.true.;post_at=clock_now()
    depinfo(:,depn)=[n,context196,depkey,requests,sent,received,tag,depn]
  end subroutine
  subroutine op_wait_begin
    if(.not.oa.or.omode<2) return
    if(.not.post_pending.or.wait_pending) error stop 'operations199: unmatched wait begin'
    wait_pending=.true.;wait_at=clock_now()
    deptimes(1,depn)=wait_at-post_at
  end subroutine
  subroutine op_wait_end
    if(.not.oa.or.omode<2) return
    if(.not.wait_pending) error stop 'operations199: unmatched wait end'
    deptimes(2,depn)=clock_now()-wait_at
    post_pending=.false.;wait_pending=.false.
  end subroutine
  subroutine flush199
    if(count199) ocount(:,context196,n)=ocount(:,context196,n)+leaf199
    leaf199=0_int64
  end subroutine
  ! Retained ABI for old hooks; Stage199 injects no sampled-access calls.
  logical function sample_begin199(id) result(take)
    integer,intent(in)::id
    take=.false.
    if(id<1) error stop 'invalid unused sample id'
  end function
  subroutine sample_end199(id,take)
    integer,intent(in)::id
    logical,intent(in)::take
    if(take.or.id<1) error stop 'Stage199 has no access clocks'
  end subroutine

  subroutine batch_begin199(id,k,l,pos,work)
    integer,intent(in)::id
    integer,optional,intent(in)::k,l,pos,work
    integer :: parent,kk,ll,pp
    call op_enter(id)
    if(.not.oa.or.omode<2.or.id==102.or.id==103.or.id==128) return
    if(nb199>=maxbatch199.or.bd199>=size(bs199)) error stop 'stage199: batch overflow'
    parent=0;kk=0;ll=0;pp=0
    if(bd199>0) then
      parent=bs199(bd199);kk=bi199(6,parent);ll=bi199(7,parent);pp=bi199(8,parent)
    end if
    if(present(k)) kk=k
    if(present(l)) ll=l
    if(present(pos)) pp=pos
    nb199=nb199+1;bd199=bd199+1;bs199(bd199)=nb199
    bi199(:,nb199)=[n,context196,nb199,parent,id,kk,ll,pp,0]
    if(present(work)) bw199(nb199)=int(work,int64)
    bl199(:,bd199)=leaf199(36:42)+ocount(36:42,context196,n)
    bm199(bd199)=sum(osec(mpi199,context196,n))
    bt199(1,nb199)=clock_now()-epoch
  end subroutine
  subroutine batch_end199(id)
    integer,intent(in)::id
    integer :: q
    if(oa.and.omode>=2.and.id/=102.and.id/=103.and.id/=128) then
      if(bd199<1) error stop 'stage199: empty batch stack'
      q=bs199(bd199)
      if(bi199(5,q)/=id.or.bi199(1,q)/=n.or.bi199(2,q)/=context196) &
        error stop 'stage199: mismatched batch'
      bt199(2,q)=clock_now()-epoch
      bt199(3,q)=sum(osec(mpi199,context196,n))-bm199(bd199)
      bc199(:,q)=leaf199(36:42)+ocount(36:42,context196,n)-bl199(:,bd199)
      if(any(bc199(:,q)<0)) error stop 'stage199: negative batch work'
      bi199(9,q)=1;bd199=bd199-1
    end if
    call op_leave(id)
  end subroutine
  subroutine batch_finish199(rank)
    integer,intent(in)::rank
    integer :: u,q
    character(80)::filename
    if(bd199/=0) error stop 'stage199: unfinished batch at output'
    write(filename,'(a,i0,a)') 'batches199-rank-',rank,'.txt'
    open(newunit=u,file=trim(filename),status='new',action='write')
    write(u,'(a,3(1x,i0))') 'schema',1,omode,rank
    do q=1,nb199
      write(u,'(9(i0,1x),8(i0,1x),3(es24.16,1x))') bi199(:,q),bw199(q),bc199(:,q),bt199(:,q)
    end do
    close(u)
  end subroutine
end module
