program profile278_unit
  use iso_fortran_env,only:int64
  use parallel_block_profile_mod
  implicit none
  integer :: i
  character(16) :: mode,option
  call get_command_argument(1,mode)
  call rhs_detail_enter(DP_RHS_PACK)
  call rhs_detail_leave(DP_RHS_PACK)
  if(any(detail_calls/=0_int64).or.detail_nsteps/=0) error stop 'disabled profile recorded work'
  if(trim(mode)=='disabled') then
     print *, 'PASS disabled RHS profile'
     stop
  end if
  detail_enabled=.true.
  call detail_enter(DP_STEP)
  call detail_enter(DP_RHS_EXCHANGE)
  call rhs_detail_enter(DP_RHS_PACK)
  call rhs_detail_leave(DP_RHS_PACK)
  call detail_leave(DP_RHS_EXCHANGE)
  call detail_leave(DP_STEP)
  call get_environment_variable('WAVETRISK_PROFILE_RHS_FINE',option)
  if(trim(option)=='1') then
     if(detail_calls(DP_RHS_PACK)/=1) error stop 'missing fine region'
  else
     if(detail_calls(DP_RHS_PACK)/=0) error stop 'fine region recorded in coarse mode'
  end if
  if(detail_calls(DP_RHS_EXCHANGE)/=1.or..not.detail_idle()) error stop 'coarse region or stack invalid'
  if(detail_nsteps/=1) error stop 'missing profile step'
  call detail_reset
  do i=1,DS_MAX+2
     call detail_enter(DP_STEP)
     call detail_leave(DP_STEP)
  end do
  if(detail_nsteps/=DS_MAX.or.detail_dropped_steps/=2) error stop 'unbounded profile buffer'
  print *, 'PASS coarse/fine RHS profile and bounded output'
end program
