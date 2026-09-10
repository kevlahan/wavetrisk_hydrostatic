program test_profile
  use iso_fortran_env, only : real64,int64
  use parallel_block_profile_mod
  implicit none
  real(real64) :: total
  call detail_enter(-1) ! Disabled instrumentation does not inspect arguments.
  call detail_leave(-1)
  if (.not.detail_idle()) error stop 'disabled stack changed'
  if (any(detail_calls/=0)) error stop 'disabled counters changed'
  detail_enabled=.true.
  call detail_enter(DP_STEP)
  call work
  call detail_enter(DP_SOURCE)
  call work
  call detail_enter(DP_SOURCE_PROOF)
  call work
  call detail_leave(DP_SOURCE_PROOF)
  call detail_enter(DP_DOMAIN_BOUNDARY)
  call detail_enter(DP_BDRY_PACK)
  call work
  call detail_leave(DP_BDRY_PACK)
  call detail_enter(DP_DOMAIN_WAIT)
  call work
  call detail_leave(DP_DOMAIN_WAIT)
  call detail_leave(DP_DOMAIN_BOUNDARY)
  call detail_enter(DP_SOURCE) ! Recursive scopes also retain exclusive accounting.
  call work
  call detail_leave(DP_SOURCE)
  call detail_leave(DP_SOURCE)
  call detail_leave(DP_STEP)
  if (.not.detail_idle()) error stop 'stack not closed'
  total=sum(detail_time(1,:))
  if (abs(total-detail_time(3,DP_STEP))>1.0e-7_real64) error stop 'self wall conservation'
  if (abs(sum(detail_time(2,:))-detail_time(4,DP_STEP))>1.0e-7_real64) error stop 'self CPU conservation'
  if (detail_calls(DP_SOURCE)/=2) error stop 'recursive calls'
  if (detail_boundary_calls(DP_SOURCE)/=1) error stop 'boundary caller'
  if (abs(detail_boundary(1,DP_SOURCE)-detail_time(3,DP_DOMAIN_BOUNDARY))>1.0e-7_real64) &
       error stop 'boundary attribution conservation'
  if (detail_nsteps/=1) error stop 'missing timestep snapshot'
  if (abs(sum(detail_steps(1,1:DP_COUNT,1))-detail_steps(1,DP_COUNT+1,1))>1.0e-7_real64) &
       error stop 'step snapshot conservation'
  call detail_add(DC_SOURCE_BYTES,100_int64)
  call detail_add(DC_SOURCE_BYTES,50_int64)
  call detail_memory([20_int64,10_int64,5_int64])
  call detail_memory([10_int64,15_int64,2_int64])
  if (detail_count(DC_SOURCE_BYTES)/=150) error stop 'counter sum'
  if (any(detail_peak/=[20_int64,15_int64,5_int64])) error stop 'peak samples'
  call detail_reset
  if (any(detail_calls/=0).or.any(detail_count/=0).or.any(detail_peak/=0)) error stop 'reset'
  if (detail_nsteps/=0.or.any(detail_boundary_calls/=0)) error stop 'snapshot reset'
  call detail_enter(DP_STEP)
  call detail_enter(DP_RESTART)
  call detail_leave(DP_RESTART)
  call detail_leave(DP_STEP)
  if (any(detail_step_info(:,1)/=[2_int64,1_int64,0_int64])) error stop 'step sequence across reset'
  do while (detail_nsteps<DS_MAX)
     call detail_enter(DP_STEP)
     call detail_leave(DP_STEP)
  end do
  call detail_enter(DP_STEP)
  call detail_leave(DP_STEP)
  if (detail_dropped_steps/=1) error stop 'snapshot overflow not reported'
  print *, 'PASS: nesting, boundary callers, timestep snapshots/overflow, disabled mode and reset'
contains
  subroutine work
    integer(int64) :: first,current,rate
    call system_clock(first,rate)
    do
       call system_clock(current)
       if (real(current-first,real64)/rate>0.005_real64) exit
    end do
  end subroutine
end program
