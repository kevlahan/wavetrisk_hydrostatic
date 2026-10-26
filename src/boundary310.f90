module boundary310_mod
  use, intrinsic :: iso_fortran_env, only: real64,int64
  implicit none
  private
  integer,parameter :: groups=12,capacity=256,stack_capacity=8,stamp_limit=160
  logical,save :: enabled=.false.,initialized=.false.,flushed=.false.,active=.false.
  integer,save :: rows=0,depth=0,stack(stack_capacity)=0,stamps=0
  integer,save :: calls(groups,capacity)=0,stamp_count(capacity)=0
  real(real64),save :: elapsed(groups,capacity)=0,totals(capacity)=0,last=0,started=0
  public :: b310_init,b310_begin,b310_enter,b310_leave,b310_end,b310_flush
contains
  subroutine b310_init
    character(16) :: option
    integer :: length,status
    if(initialized) error stop '310 duplicate initialization'
    option=''
    call get_environment_variable('WAVETRISK_BOUNDARY310',option,length,status)
    if(status/=0.and.status/=1) error stop '310 invalid recorder option'
    if(len_trim(option)>0.and.trim(option)/='0'.and.trim(option)/='1') error stop '310 invalid recorder option'
    enabled=trim(option)=='1';initialized=.true.
    ! Pre-touch bounded storage outside the numerical trajectory.
    elapsed=0;totals=0;calls=0;stamp_count=0
  end subroutine
  function clock310() result(t)
    real(real64) :: t
    integer(int64) :: tick,rate
    call system_clock(tick,rate)
    if(rate<=0) error stop '310 unavailable monotonic clock'
    stamps=stamps+1
    if(stamps>stamp_limit) error stop '310 timestamp budget exceeded'
    t=real(tick,real64)/real(rate,real64)
  end function
  subroutine charge(t)
    real(real64),intent(in) :: t
    if(t<last) error stop '310 reversed clock'
    elapsed(stack(depth),rows)=elapsed(stack(depth),rows)+t-last
    last=t
  end subroutine
  subroutine b310_begin
    if(.not.initialized) error stop '310 recorder not initialized'
    if(.not.enabled) return
    if(active.or.flushed.or.rows==capacity) error stop '310 invalid begin or capacity exceeded'
    rows=rows+1;depth=1;stack(1)=1;calls(1,rows)=1;stamps=0;active=.true.
    started=clock310();last=started
  end subroutine
  subroutine b310_enter(group)
    integer,intent(in) :: group
    real(real64) :: t
    ! Setup outside measured steps never starts a partial interval.
    if(.not.active) return
    if(group<2.or.group>groups.or.depth==stack_capacity) error stop '310 invalid enter'
    t=clock310();call charge(t)
    depth=depth+1;stack(depth)=group;calls(group,rows)=calls(group,rows)+1
  end subroutine
  subroutine b310_leave(group)
    integer,intent(in) :: group
    real(real64) :: t
    if(.not.active) return
    if(depth<2) error stop '310 empty scope'
    if(stack(depth)/=group) error stop '310 mismatched scope'
    t=clock310();call charge(t);depth=depth-1
  end subroutine
  subroutine b310_end
    real(real64) :: t
    if(.not.enabled) return
    if(.not.active.or.depth/=1) error stop '310 invalid end'
    t=clock310();call charge(t)
    totals(rows)=t-started;stamp_count(rows)=stamps;active=.false.;depth=0
  end subroutine
  subroutine b310_flush(rank)
    integer,intent(in) :: rank
    integer :: unit,j
    character(80) :: filename
    if(.not.enabled) return
    if(active.or.flushed.or.rows<1) error stop '310 invalid flush'
    write(filename,'(a,i0,a)') 'boundary310-rank-',rank,'.txt'
    open(newunit=unit,file=trim(filename),status='new',action='write')
    write(unit,*) 310,rank,rows,groups,stamp_limit
    do j=1,rows
       write(unit,'(3(i0,1x),es25.17,1x,12(i0,1x),12(es25.17,1x))') &
            j,rank,stamp_count(j),totals(j),calls(:,j),elapsed(:,j)
    end do
    close(unit);flushed=.true.
    if(rank==0) write(*,'(a,i0)') 'Stage310 boundary records after trajectory: ',rows
  end subroutine
end module
