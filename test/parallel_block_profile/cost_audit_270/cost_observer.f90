module cost_observer_mod
  use iso_fortran_env, only : int64, real64
  implicit none
  private
  integer, parameter :: max_steps=64,max_depth=32
  integer, parameter, public :: C_STEP=1, C_RK=2, C_TENDENCY=3, C_WT=4, &
       C_ADAPT=5, C_DOMAIN_INVERSE=6, C_NATIVE_INVERSE=7, C_INVERSE_SETUP=8, &
       C_PUBLICATION=9, C_INSTALL=10, C_PROVISIONAL=11, C_HYDRO=12, &
       C_RECOMPOSE=13, C_CATALOG=14, C_SOURCE=15, C_MIGRATE=16, &
       C_POSTGRID=17, C_COMPRESSION=18, C_FORWARD_SCALAR=19, C_FORWARD_VECTOR=20, &
       C_GHOST=21, C_VSOURCE=22, C_VGRADIENT=23, C_PRESSURE=24, &
       C_BOUNDARY=25, C_ACCEPTED=26, C_INVERSE_MPI=27, C_COMPAT=28, &
       C_CAPTURE_BEGIN=29, C_CAPTURE_PUBLISH=30, C_VELOCITY_REMAINDER=31, C_RK_ASSEMBLE=32, &
       C_SCALAR_RESTRICT=33, C_SCALAR_LOCAL=34, C_RESTRICT_GHOST=35, C_RESTRICT_PACK=36, &
       C_RESTRICT_MPI=37, C_RESTRICT_INSTALL=38, C_TEMP_EDGES=39, C_TEMP_WAIT=40, &
       C_EXNER_TRANSACTION=41, C_EXNER_KERNEL=42, C_ACCUM_RETAIN=43, C_EDGE_METRIC=44, &
       C_DOMAIN_PAYLOAD=45, C_MIG_MANIFEST=46, C_MIG_PACK=47, C_MIG_UNPACK=48, &
       C_MIG_PAYLOAD=49, C_MIG_SIZES=50, C_STORE_INSTALL=51, C_GHOST_PLAN=52, &
       C_WRITEBACK_PLAN=53, C_PRODUCTION_PREP=54, C_CAPTURE_MPI=55, C_VELOCITY_PLAN=56, &
       C_MASS_PLAN=57, C_MIG_TRANSFER_MPI=58, C_MIG_CHECKSUM=59, C_MIG_PROOF_MPI=60, C_COUNT=60
  character(48), parameter, public :: cost_name(C_COUNT)=[character(48) :: &
       'step', 'RK', 'tendency context', &
       'wavelet transaction', 'adapt algorithm', 'Domain inverse', &
       'native inverse', 'inverse setup', 'canonical publication', &
       'block installation', 'provisional refresh', 'block hydrostatic cache', &
       'native recomposition', 'catalogue', 'source blocks', &
       'migration', 'post-grid refresh', 'native compression', &
       'native scalar forward', 'native vector forward', 'block ghosts', &
       'retired velocity source', 'retired velocity gradient', 'Domain surface pressure', &
       'retired Domain boundary', 'accepted state refresh', 'inverse MPI wait/traffic', &
       'shared compatibility producer', 'scalar capture preparation', 'scalar capture publication', &
       'velocity remainder preparation', 'RK state assembly', 'scalar restriction transaction', &
       'restriction local kernels', 'restriction ghost transaction', 'restriction packing/local copies', &
       'restriction MPI post/wait', 'restriction received installation', 'temperature edge transaction', &
       'temperature edge MPI wait', 'final physical transaction', 'final physical kernel batch', &
       'tendency accumulator retention', 'edge metric preparation', 'Domain to block payload', &
       'migration manifest/check', 'migration packing batch', 'migration unpacking batch', &
       'migration payload exchange', 'migration size exchange', 'complete store installation', &
       'ghost plan preparation', 'writeback plan preparation', 'production workspace preparation', &
       'scalar capture MPI', 'native velocity plan preparation', 'native mass plan preparation', &
       'migration payload MPI', 'migration byte checksum scans', 'migration proof reductions']
  logical :: enabled=.false.,configured=.false.,active=.false.
  integer :: depth=0,nsteps=0,stack(max_depth)=0
  real(real64) :: started(2,max_depth)=0,last(2)=0
  real(real64) :: values(4,C_COUNT,max_steps)=0
  integer(int64) :: calls(C_COUNT,max_steps)=0,work(3,C_COUNT,max_steps)=0
  public :: cost_begin_step,cost_end_step,cost_enter,cost_leave,cost_work,cost_flush
contains
  subroutine stamp(t)
    real(real64), intent(out) :: t(2)
    integer(int64) :: tick,rate
    call system_clock(tick,rate)
    if(rate<=0) error stop 'Stage270 unavailable monotonic clock'
    t(1)=real(tick,real64)/real(rate,real64)
    call cpu_time(t(2))
  end subroutine
  subroutine account(t)
    real(real64), intent(in) :: t(2)
    if(depth>0) values(1:2,stack(depth),nsteps)=values(1:2,stack(depth),nsteps)+max(0.0_real64,t-last)
    last=t
  end subroutine
  subroutine cost_begin_step
    character(8) :: setting
    integer :: status,length
    if(.not.configured) then
       call get_environment_variable('WAVETRISK_COST_AUDIT',setting,length,status)
       if(status/=1) then
          if(status/=0.or.length/=1) error stop 'Stage270 COST_AUDIT must be 0 or 1'
          if(setting(1:1)/='0'.and.setting(1:1)/='1') error stop 'Stage270 invalid COST_AUDIT'
          enabled=setting(1:1)=='1'
       end if
       configured=.true.
    end if
    if(.not.enabled) return
    if(active.or.depth/=0.or.nsteps>=max_steps) error stop 'Stage270 step lifetime/overflow'
    nsteps=nsteps+1
    active=.true.
    call cost_enter(C_STEP)
  end subroutine
  subroutine cost_end_step
    if(.not.enabled) return
    if(.not.active.or.depth/=1) error stop 'Stage270 incomplete scopes'
    call cost_leave(C_STEP)
    active=.false.
  end subroutine
  subroutine cost_enter(p)
    integer, intent(in) :: p
    real(real64) :: t(2)
    if(.not.active) return
    if(p<1.or.p>C_COUNT.or.depth>=max_depth) error stop 'Stage270 invalid scope'
    call stamp(t)
    call account(t)
    depth=depth+1
    stack(depth)=p
    started(:,depth)=t
    calls(p,nsteps)=calls(p,nsteps)+1_int64
  end subroutine
  subroutine cost_leave(p)
    integer, intent(in) :: p
    real(real64) :: t(2)
    if(.not.active) return
    if(depth<1) error stop 'Stage270 empty stack'
    if(stack(depth)/=p) error stop 'Stage270 scope mismatch'
    call stamp(t)
    call account(t)
    values(3:4,p,nsteps)=values(3:4,p,nsteps)+max(0.0_real64,t-started(:,depth))
    depth=depth-1
  end subroutine
  subroutine cost_work(p,a,b,c)
    integer, intent(in) :: p
    integer(int64), intent(in) :: a,b,c
    if(.not.active) return
    if(p<1.or.p>C_COUNT.or.min(a,b,c)<0_int64) error stop 'Stage270 invalid work'
    work(:,p,nsteps)=work(:,p,nsteps)+[a,b,c]
  end subroutine
  subroutine cost_flush(rank)
    integer, intent(in) :: rank
    integer :: unit,s,p
    character(40) :: filename
    if(.not.enabled) return
    if(active.or.depth/=0.or.nsteps<1) error stop 'Stage270 invalid final observer state'
    write(filename,'(a,i6.6,a)') 'cost-rank-',rank,'.txt'
    open(newunit=unit,file=trim(filename),status='new',action='write')
    write(unit,'(a,3(1x,i0))') 'STAGE270_COST_V1',rank,nsteps,C_COUNT
    do s=1,nsteps
       do p=1,C_COUNT
          write(unit,'(3(i0,1x),4(es25.17e3,1x),3(i0,1x))') s,p,calls(p,s),values(:,p,s),work(:,p,s)
       end do
    end do
    close(unit)
  end subroutine
end module cost_observer_mod
