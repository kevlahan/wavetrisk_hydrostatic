module runtime322_mod
  use iso_c_binding, only: c_int,c_int64_t
  implicit none
  private
  public :: init_runtime322,observe_runtime322,start_window322,end_window322,finish_runtime322
  logical :: placement322=.false.,window322=.false.,qualification322=.false.
  integer(c_int64_t) :: begin322=0_c_int64_t,end322=0_c_int64_t
  interface
     function monotonic322() bind(C) result(value)
       import c_int64_t
       integer(c_int64_t) :: value
     end function
     function observe322(stage) bind(C) result(code)
       import c_int
       integer(c_int),value :: stage
       integer(c_int) :: code
     end function
     function flush_runtime322(rank,first,last) bind(C) result(code)
       import c_int,c_int64_t
       integer(c_int),value :: rank
       integer(c_int64_t),value :: first,last
       integer(c_int) :: code
     end function
  end interface
contains
  logical function enabled322(name) result(enabled)
    character(*),intent(in) :: name
    character(16) :: option
    integer :: status
    option='';call get_environment_variable(name,option,status=status)
    if(status/=0.and.status/=1) error stop 'Stage322 runtime selector unreadable'
    if(trim(option)/=''.and.trim(option)/='0'.and.trim(option)/='1') &
         error stop 'Stage322 runtime selector must be 0 or 1'
    enabled=trim(option)=='1'
  end function
  subroutine init_runtime322()
    placement322=enabled322('WAVETRISK_PLACEMENT322')
    window322=enabled322('WAVETRISK_WINDOW322')
    qualification322=enabled322('WAVETRISK_EXECUTION_AUDIT316')
    if(window322.and..not.placement322) error stop 'Stage322 window requires placement'
    if(window322.and.qualification322) error stop 'Stage322 bare timing cannot enable audit'
    call observe_runtime322(0)
  end subroutine
  subroutine observe_runtime322(stage)
    integer,intent(in) :: stage
    integer(c_int) :: code
    if(.not.placement322) return
    if(stage>=2.and.stage<=5.and..not.qualification322) return
    code=observe322(int(stage,c_int))
    if(code/=0) then
       write(*,*) 'Stage322 solver placement error',stage,code
       error stop 'Stage322 solver placement gate failed'
    end if
  end subroutine
  subroutine start_window322()
    if(.not.window322) return
    if(begin322/=0) error stop 'Stage322 window began twice'
    begin322=monotonic322()
    if(begin322<=0) error stop 'Stage322 start clock failed'
  end subroutine
  subroutine end_window322()
    if(.not.window322) return
    if(begin322<=0.or.end322/=0) error stop 'Stage322 window lifecycle'
    end322=monotonic322()
    if(end322<=begin322) error stop 'Stage322 end clock failed'
  end subroutine
  subroutine finish_runtime322(rank,steps)
    integer,intent(in) :: rank,steps
    integer(c_int) :: code
    if(.not.placement322) return
    if(steps/=4) error stop 'Stage322 native runtime requires four steps'
    if(window322.and.(begin322<=0.or.end322<=begin322)) error stop 'Stage322 missing window'
    call observe_runtime322(6)
    code=flush_runtime322(int(rank,c_int),begin322,end322)
    if(code/=0) error stop 'Stage322 runtime flush failed'
  end subroutine
end module runtime322_mod
