module parallel_block_profile_mod
  ! Optional coarse-region attribution. No MPI, allocation or clocks when off.
  ! CPU_TIME is process CPU (including MPI spinning), NOT useful-work time.
  use iso_fortran_env, only : int64, real64
  implicit none
  private
  integer, parameter, public :: DP_STEP=1, DP_DYNAMICS=2, DP_ADAPT=3, DP_OUTPUT=4, &
       DP_SOURCE=5, DP_SOURCE_PROOF=6, DP_RECEIVE_PROOF=7, DP_INSTALL=8, DP_COPY_PROOF=9, &
       DP_MIGRATE=10, DP_PACK=11, DP_UNPACK=12, DP_MIGRATE_MPI=13, DP_INVERSE=14, &
       DP_INVERSE_PLAN=15, DP_INVERSE_PACK=16, DP_INVERSE_LOCAL=17, DP_INVERSE_POST=18, &
       DP_INVERSE_WAIT=19, DP_INVERSE_INSTALL=20, DP_OUTER_KERNEL=21, DP_SCALAR_SETUP=22, &
       DP_GEOMETRY_EXPAND=23, DP_SCALAR_TRANSPORT=24, DP_SCALAR_REPLAY=25, DP_SHARED=26, &
       DP_BASIC=27, DP_DOMAIN_BOUNDARY=28, DP_PRIMITIVE=29, DP_NATIVE_MASS=30, &
       DP_NATIVE_VELOCITY=31, DP_WRITEBACK=32, DP_GHOST=33, DP_GHOST_MPI=34, &
       DP_CATALOG=35, DP_INVERSE_KERNEL=36, DP_SCALAR_MPI=37, DP_DOMAIN_WAIT=38, &
       DP_MASS_WAIT=39, DP_WRITEBACK_MPI=40, DP_REMOTE_GEOMETRY=41, DP_INVERSE_SEED=42, &
       DP_WAVELET=43, DP_COMPRESSION=44, DP_RK_ASSEMBLE=45, DP_RETAIN=46, DP_TENDENCY=47, &
       DP_GRID_CHANGE=48, DP_TREND_REFRESH=49, DP_STATE_REFRESH=50, DP_CAPTURE=51, DP_VELOCITY_PREP=52, &
       DP_RK_COMPAT=53, DP_PHYSICS=54, DP_REMAP=55, DP_RESTART=56, DP_ORACLE=57, &
       DP_BDRY_PACK=58, DP_BDRY_ROUTE=59, DP_BDRY_POST=60, DP_BDRY_INSTALL=61, DP_BDRY_LOCAL=62, &
       DP_BDRY_STORAGE=63, DP_OUTSIDE=64, DP_WT_DRIVER=65, DP_COUNT=65
  character(32), parameter, public :: detail_name(DP_COUNT)=[character(32) :: &
       'timestep residual', 'dynamics residual', 'adaptation residual', 'output/checkpoint', &
       'source extraction', 'source serialization proof', 'receive serialization proof', &
       'local install/deep copy', 'local copy byte proof', 'migration residual', &
       'migration pack', 'migration unpack', 'migration MPI', 'inverse residual', &
       'inverse plan', 'inverse pack', 'inverse local aliases', 'inverse MPI post/self copy', &
       'inverse MPI wait', 'inverse install', 'outer inverse arithmetic', 'scalar setup/storage', &
       'scalar geometry expansion', 'scalar transport/install', 'scalar replay', &
       'shared producer residual', 'shared basic operators', 'Domain boundary interface', &
       'primitive/physics preparation', 'native mass', 'native velocity', &
       'Domain writeback', 'block ghost pack/install', 'block ghost MPI', &
       'catalog construction', 'scalar/inner inverse kernels', 'scalar transport MPI', &
       'Domain boundary MPI wait', 'native mass MPI wait', 'Domain writeback MPI', &
       'received geometry expansion', 'inverse scaffold staging', &
       'native wavelet work', 'native compression', 'native RK assembly', 'RK retention/publication', &
       'tendency assembly residual', 'grid-change synchronization', 'trend state refresh', &
       'prognostic state refresh', 'scalar capture lifecycle', 'velocity remainder preparation', &
       'RK compatibility closure', 'split physics consumers', 'vertical remap', 'restart bootstrap', &
       'validation oracle', 'boundary scan and pack', 'boundary receive route scan', &
       'boundary MPI posting', 'boundary scan and install', 'boundary local copies', &
       'boundary buffer allocation/zero', 'outside instrumented context', 'wavelet driver residual']
  logical, public, save :: detail_enabled=.false.
  ! self wall, self process CPU, inclusive wall, inclusive process CPU.
  real(real64), public, save :: detail_time(4,DP_COUNT)=0.0_real64
  integer(int64), public, save :: detail_calls(DP_COUNT)=0_int64
  integer, parameter, public :: DC_SOURCE_BYTES=1, DC_RECEIVE_BYTES=2, DC_COPY_BYTES=3, &
       DC_GEOMETRY_WRITES=4, DC_FIELD_SAMPLES=5, DC_REMOTE_GEOMETRY_WRITES=6, &
       DC_INVERSE_BUILD=7, DC_INVERSE_REUSE=8, DC_BDRY_SEND=9, DC_BDRY_RECV=10, &
       DC_BDRY_REQUESTS=11, DC_BDRY_PACK_VISITS=12, DC_BDRY_ROUTE_VISITS=13, &
       DC_BDRY_INSTALL_VISITS=14, DC_BDRY_ZERO=15, DC_BDRY_FIELDS=16, DC_COUNT=16
  character(32), parameter, public :: detail_counter_name(DC_COUNT)=[character(32) :: &
       'source proof payload bytes', 'receive proof payload bytes', 'copy proof payload bytes', &
       'local geometry doubles written', 'local scalar field samples', 'received geometry doubles copied', &
       'inverse plan builds', 'inverse plan reuses', 'boundary doubles sent', 'boundary doubles received', &
       'boundary posted requests', 'boundary pack candidates', 'boundary route candidates', &
       'boundary install candidates', 'boundary buffer doubles zeroed', 'boundary fields per start']
  integer(int64), public, save :: detail_count(DC_COUNT)=0_int64
  integer, parameter, public :: DM_SCALAR_BYTES=1, DM_GEOMETRY_BYTES=2, DM_SHARED_BYTES=3, DM_COUNT=3
  character(32), parameter, public :: detail_memory_name(DM_COUNT)=[character(32) :: &
       'scalar record storage bytes', 'geometry subset bytes', 'shared-once geometry estimate']
  integer(int64), public, save :: detail_peak(DM_COUNT)=0_int64
  integer, parameter :: STACK_SIZE=128
  integer :: depth=0, phase_stack(STACK_SIZE)=0
  integer :: boundary_owner(STACK_SIZE)=0
  ! Cross-cut attribution, NOT extra self time. Each boundary event is charged
  ! to its nearest instrumented non-boundary caller, including nested waits.
  real(real64), public :: detail_boundary(2,DP_COUNT)=0
  integer(int64), public :: detail_boundary_calls(DP_COUNT)=0
  ! Bounded per-step records, flushed only at the existing window report.
  ! No per-step MPI or I/O. Overflow is explicit, never silently extrapolated.
  integer, parameter, public :: DS_MAX=256
  integer, public :: detail_nsteps=0, detail_dropped_steps=0
  real(real64), public :: detail_steps(2,DP_COUNT+1,DS_MAX)=0
  integer(int64), public :: detail_step_info(3,DS_MAX)=0
  integer(int64) :: step_sequence=0, step_calls_before(DP_COUNT)=0
  real(real64) :: step_before(4,DP_COUNT)=0
  real(real64) :: start_wall(STACK_SIZE), start_cpu(STACK_SIZE), last_wall=0, last_cpu=0
  public :: detail_enter, detail_leave, detail_reset, detail_add, detail_memory, detail_idle
contains
  subroutine stamp(w,c)
    real(real64), intent(out) :: w,c
    integer(int64) :: ticks,rate
    call system_clock(ticks,rate)
    if (rate<=0) error stop 'detail profiler: unavailable monotonic clock'
    w=real(ticks,real64)/real(rate,real64)
    call cpu_time(c)
  end subroutine

  subroutine account(w,c)
    real(real64), intent(in) :: w,c
    integer :: p
    if (depth>0) then
       p=phase_stack(depth)
       detail_time(1,p)=detail_time(1,p)+max(0.0_real64,w-last_wall)
       detail_time(2,p)=detail_time(2,p)+max(0.0_real64,c-last_cpu)
       if (boundary_owner(depth)>0) then
          p=boundary_owner(depth)
          detail_boundary(1,p)=detail_boundary(1,p)+max(0.0_real64,w-last_wall)
          detail_boundary(2,p)=detail_boundary(2,p)+max(0.0_real64,c-last_cpu)
       end if
    end if
    last_wall=w
    last_cpu=c
  end subroutine

  subroutine detail_enter(p)
    integer, intent(in) :: p
    real(real64) :: w,c
    if (.not.detail_enabled) return
    if (p<1.or.p>DP_COUNT.or.depth>=STACK_SIZE) error stop 'detail profiler: invalid enter'
    if (p==DP_STEP) then
       if (depth/=0) error stop 'detail profiler: nested timestep'
       step_before=detail_time
       step_calls_before=detail_calls
    end if
    call stamp(w,c)
    call account(w,c)
    depth=depth+1
    boundary_owner(depth)=0
    if (depth>1) boundary_owner(depth)=boundary_owner(depth-1)
    if (p==DP_DOMAIN_BOUNDARY.and.boundary_owner(depth)==0) then
       boundary_owner(depth)=DP_OUTSIDE
       if (depth>1) boundary_owner(depth)=phase_stack(depth-1)
       detail_boundary_calls(boundary_owner(depth))=detail_boundary_calls(boundary_owner(depth))+1_int64
    end if
    phase_stack(depth)=p
    start_wall(depth)=w
    start_cpu(depth)=c
    detail_calls(p)=detail_calls(p)+1_int64
  end subroutine

  subroutine detail_leave(p)
    integer, intent(in) :: p
    real(real64) :: w,c
    if (.not.detail_enabled) return
    if (depth<1) error stop 'detail profiler: empty stack'
    if (phase_stack(depth)/=p) error stop 'detail profiler: mismatched scope'
    call stamp(w,c)
    call account(w,c)
    detail_time(3,p)=detail_time(3,p)+max(0.0_real64,w-start_wall(depth))
    detail_time(4,p)=detail_time(4,p)+max(0.0_real64,c-start_cpu(depth))
    depth=depth-1
    if (p==DP_STEP) then
       step_sequence=step_sequence+1_int64
       if (detail_nsteps<DS_MAX) then
          detail_nsteps=detail_nsteps+1
          detail_steps(:,1:DP_COUNT,detail_nsteps)=detail_time(1:2,:)-step_before(1:2,:)
          detail_steps(:,DP_COUNT+1,detail_nsteps)=detail_time(3:4,DP_STEP)-step_before(3:4,DP_STEP)
          detail_step_info(:,detail_nsteps)=[step_sequence, &
               detail_calls(DP_RESTART)-step_calls_before(DP_RESTART), &
               detail_calls(DP_REMAP)-step_calls_before(DP_REMAP)]
       else
          detail_dropped_steps=detail_dropped_steps+1
       end if
    end if
  end subroutine

  subroutine detail_add(p,n)
    integer, intent(in) :: p
    integer(int64), intent(in) :: n
    if (.not.detail_enabled) return
    if (p<1.or.p>DC_COUNT.or.n<0) error stop 'detail profiler: invalid counter'
    detail_count(p)=detail_count(p)+n
  end subroutine

  subroutine detail_memory(values)
    integer(int64), intent(in) :: values(DM_COUNT)
    if (.not.detail_enabled) return
    detail_peak=max(detail_peak,values)
  end subroutine

  logical function detail_idle()
    detail_idle=depth==0
  end function

  subroutine detail_reset
    if (depth/=0) error stop 'detail profiler: reset inside scope'
    detail_time=0
    detail_calls=0_int64
    detail_count=0_int64
    detail_peak=0_int64
    detail_boundary=0
    detail_boundary_calls=0_int64
    detail_nsteps=0
    detail_dropped_steps=0
  end subroutine
end module parallel_block_profile_mod
