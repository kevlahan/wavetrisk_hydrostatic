module campaign312_mod
  use iso_fortran_env, only: int64,real64
  implicit none
  integer :: steps312=-1,step312=0,nrhs312=0
  integer(int64) :: rhs312(12,32)=0
  real(real64) :: timing312(4)=0,physical312(4)=0,dt312(4)=0
  integer(int64) :: active312(4)=0
contains
  subroutine configure312
    character(32) :: option
    integer :: ios,status
    call get_environment_variable('WAVETRISK_CAMPAIGN312_STEPS',option,status=status)
    if(status/=0) error stop 'Stage312 requires explicit step limit'
    read(option,*,iostat=ios) steps312
    if(ios/=0.or.steps312<0.or.steps312>4) error stop 'Stage312 step limit must be 0..4'
  end subroutine
  subroutine record_rhs312(values)
    integer(int64),intent(in) :: values(11)
    if(step312<1.or.step312>steps312) error stop 'Stage312 RHS outside timestep'
    nrhs312=nrhs312+1
    if(nrhs312>size(rhs312,2)) error stop 'Stage312 RHS record overflow'
    rhs312(:,nrhs312)=[int(step312,int64),values]
  end subroutine
  subroutine flush312(rank,startup)
    integer,intent(in) :: rank
    real(real64),intent(in) :: startup
    integer :: u,i
    character(80) :: name
    write(name,'(a,i0,a)') 'campaign312-rank-',rank,'.txt'
    open(newunit=u,file=trim(name),status='new')
    write(u,*) 'H',rank,steps312,nrhs312,startup
    do i=1,steps312
       write(u,*) 'T',i,timing312(i),physical312(i),dt312(i),active312(i)
    end do
    do i=1,nrhs312
       write(u,*) 'R',rhs312(:,i)
    end do
    close(u)
  end subroutine
end module
