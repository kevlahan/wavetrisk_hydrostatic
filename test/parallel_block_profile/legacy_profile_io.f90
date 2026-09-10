module legacy_profile_io_mod
  ! Experimental archive only. No MPI and no output inside a measured timestep.
  use parallel_block_profile_mod
  implicit none
contains
  subroutine legacy_profile_begin
    character(16) :: value
    integer :: status
    call get_environment_variable('WAVETRISK_PROFILE_BLOCK_DETAIL',value,status=status)
    detail_enabled=status==0.and.trim(value)=='1'
    call detail_reset
  end subroutine

  subroutine legacy_profile_end(rank)
    integer, intent(in) :: rank
    integer :: unit, p, s
    character(80) :: filename
    if (.not.detail_enabled) return
    if (.not.detail_idle()) error stop 'legacy profile: unclosed scope'
    if (detail_dropped_steps/=0) error stop 'legacy profile: step overflow'
    write(filename,'(a,i0,a)') 'legacy-detail-rank-',rank,'.txt'
    open(newunit=unit,file=trim(filename),status='new',action='write')
    write(unit,'(a)') '# Experimental legacy main: IDs reuse block labels, not identical implementation boundaries.'
    do p=1,DP_COUNT
       write(unit,'(a,1x,i0,1x,a)') 'name',p,trim(detail_name(p))
       write(unit,'(a,1x,i0,4(1x,es24.16),1x,i0)') 'region',p,detail_time(:,p),detail_calls(p)
       write(unit,'(a,1x,i0,2(1x,es24.16),1x,i0)') &
            'boundary',p,detail_boundary(:,p),detail_boundary_calls(p)
    end do
    do p=1,DC_COUNT
       write(unit,'(a,1x,i0,1x,i0)') 'counter',p,detail_count(p)
    end do
    do s=1,detail_nsteps
       write(unit,'(a,3(1x,i0),2(1x,es24.16))') &
            'step',detail_step_info(:,s),detail_steps(:,DP_COUNT+1,s)
       do p=1,DP_COUNT
          write(unit,'(a,2(1x,i0),2(1x,es24.16))') 'self',s,p,detail_steps(:,p,s)
       end do
    end do
    close(unit)
  end subroutine
end module
