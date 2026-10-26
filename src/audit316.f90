module audit316_mod
  ! Qualification-only component records, losslessly compressed in this process.
  ! No payload replay, collective, clock or numerical substitution is used.
  use iso_fortran_env, only: int8,int32,int64,real64
  use iso_c_binding, only: c_char,c_int,c_int64_t,c_ptr,c_loc,c_null_char
  implicit none
  private
  public :: audit316_init,audit316_finish,audit316_enabled,audit316_begin,audit316_end
  public :: audit316_i,audit316_i64,audit316_r,audit316_l
  integer :: depth316=0
  integer(int64) :: records316=0_int64,bytes316=0_int64
  integer(int64),parameter :: limit316=512_int64*1024_int64*1024_int64
  integer,parameter :: max_depth316=32,max_tag316=256
  character(max_tag316) :: scopes316(max_depth316)
  integer(int64) :: epochs316(max_depth316)
  logical :: initialized316=.false.,enabled316=.false.,finished316=.false.
  interface
     integer(c_int) function compress_open322(path) bind(C,name='audit_compress322_open')
       import :: c_int,c_char
       character(kind=c_char),intent(in) :: path(*)
     end function
     integer(c_int) function compress_write322(address,length) bind(C,name='audit_compress322_write')
       import :: c_int,c_ptr,c_int64_t
       type(c_ptr),value :: address
       integer(c_int64_t),value :: length
     end function
     integer(c_int) function compress_close322() bind(C,name='audit_compress322_close')
       import :: c_int
     end function
  end interface
contains
  subroutine write_bytes322(address,length)
    type(c_ptr),intent(in) :: address
    integer(int64),intent(in) :: length
    if(compress_write322(address,int(length,c_int64_t))/=0_c_int) error stop 'Stage322 compressed audit write failed'
  end subroutine

  subroutine audit316_init(rank)
    integer,intent(in) :: rank
    character(16) :: option
    character(80) :: filename
    character(8),target :: magic
    integer(int32),target :: header(5)
    integer :: status
    if(initialized316) error stop 'Stage316 audit initialized twice'
    initialized316=.true.
    option='';call get_environment_variable('WAVETRISK_EXECUTION_AUDIT316',option,status=status)
    if(status/=0.and.status/=1) error stop 'Stage316 audit environment read failed'
    if(trim(option)/=''.and.trim(option)/='0'.and.trim(option)/='1') &
         error stop 'Stage316 audit requires 0 or 1'
    enabled316=trim(option)=='1'
    if(.not.enabled316) return
    if(rank<0.or.storage_size(0)/=32.or.storage_size(0_int64)/=64.or.storage_size(0.0_real64)/=64.or. &
         storage_size('a')/=8.or.kind('a')/=c_char) error stop 'Stage316 unsupported integer or real representation'
    write(filename,'(a,i0,a)') 'execution316-rank-',rank,'.bin.gz'
    if(compress_open322(trim(filename)//c_null_char)/=0_c_int) error stop 'Stage322 compressed audit open failed'
    magic='RK316V01'
    header=[int(1,int32),int(rank,int32),int(4,int32),int(8,int32),int(z'01020304',int32)]
    call write_bytes322(c_loc(magic(1:1)),8_int64)
    call write_bytes322(c_loc(header(1)),20_int64)
    bytes316=28_int64
  end subroutine

  logical function audit316_enabled() result(on)
    on=enabled316.and..not.finished316
  end function

  subroutine header316(name,kind,n,width)
    character(*),intent(in) :: name
    integer,intent(in) :: kind,width
    integer(int64),intent(in) :: n
    integer :: k
    integer(int64) :: added
    integer(int32),target :: descriptor(2)
    integer(int64),target :: count_value(1)
    character(max_tag316),target :: tag
    k=len_trim(name)
    if(.not.initialized316.or..not.enabled316.or.finished316) error stop 'Stage316 inactive write'
    if(k<1.or.k>max_tag316.or.n<0.or.width<0) error stop 'Stage316 invalid record descriptor'
    if(n>limit316/max(1,width)) error stop 'Stage316 record size exceeds bound'
    added=16_int64+k+n*width
    if(bytes316>limit316-added) error stop 'Stage316 audit file exceeds 512MiB bound'
    descriptor=[int(k,int32),int(kind,int32)];count_value=[n];tag=name
    call write_bytes322(c_loc(descriptor(1)),8_int64)
    call write_bytes322(c_loc(count_value(1)),8_int64)
    call write_bytes322(c_loc(tag(1:1)),int(k,int64))
    records316=records316+1_int64;bytes316=bytes316+added
  end subroutine

  subroutine audit316_begin(name,epoch)
    character(*),intent(in) :: name
    integer(int64),intent(in) :: epoch
    integer(int64),target :: epoch_value(1)
    if(.not.audit316_enabled()) return
    if(depth316==max_depth316) error stop 'Stage316 scope nesting exceeds bound'
    call header316(name,10,1_int64,8)
    epoch_value=[epoch]
    call write_bytes322(c_loc(epoch_value(1)),8_int64)
    depth316=depth316+1;scopes316(depth316)=name;epochs316(depth316)=epoch
  end subroutine

  subroutine audit316_end()
    integer(int64),target :: epoch_value(1)
    if(.not.audit316_enabled()) return
    if(depth316<1) error stop 'Stage316 unmatched scope end'
    call header316(trim(scopes316(depth316)),11,1_int64,8)
    epoch_value=[epochs316(depth316)]
    call write_bytes322(c_loc(epoch_value(1)),8_int64)
    depth316=depth316-1
  end subroutine

  subroutine audit316_i(name,value)
    character(*),intent(in) :: name
    integer,intent(in),target :: value(:)
    integer,allocatable,target :: packed(:)
    if(.not.audit316_enabled()) return
    call header316(name,1,int(size(value),int64),4)
    if(size(value)>0) then
       if(is_contiguous(value)) then
          call write_bytes322(c_loc(value(1)),4_int64*size(value))
       else
          allocate(packed(size(value)));packed=value
          call write_bytes322(c_loc(packed(1)),4_int64*size(packed))
       end if
    end if
  end subroutine

  subroutine audit316_i64(name,value)
    character(*),intent(in) :: name
    integer(int64),intent(in),target :: value(:)
    integer(int64),allocatable,target :: packed(:)
    if(.not.audit316_enabled()) return
    call header316(name,2,int(size(value),int64),8)
    if(size(value)>0) then
       if(is_contiguous(value)) then
          call write_bytes322(c_loc(value(1)),8_int64*size(value))
       else
          allocate(packed(size(value)));packed=value
          call write_bytes322(c_loc(packed(1)),8_int64*size(packed))
       end if
    end if
  end subroutine

  subroutine audit316_r(name,value)
    character(*),intent(in) :: name
    real(real64),intent(in),target :: value(:)
    real(real64),allocatable,target :: packed(:)
    if(.not.audit316_enabled()) return
    call header316(name,3,int(size(value),int64),8)
    if(size(value)>0) then
       if(is_contiguous(value)) then
          call write_bytes322(c_loc(value(1)),8_int64*size(value))
       else
          allocate(packed(size(value)));packed=value
          call write_bytes322(c_loc(packed(1)),8_int64*size(packed))
       end if
    end if
  end subroutine

  subroutine audit316_l(name,value)
    character(*),intent(in) :: name
    logical,intent(in) :: value(:)
    integer(int8),allocatable,target :: encoded(:)
    if(.not.audit316_enabled()) return
    call header316(name,4,int(size(value),int64),1)
    if(size(value)>0) then
       allocate(encoded(size(value)))
       encoded=int(merge(1,0,value),int8)
       call write_bytes322(c_loc(encoded(1)),int(size(encoded),int64))
    end if
  end subroutine

  subroutine audit316_finish()
    integer(int64),target :: counts(2)
    if(.not.initialized316.or.finished316) error stop 'Stage316 audit finish lifecycle'
    if(enabled316) then
       if(depth316/=0) error stop 'Stage316 unfinished execution scope'
       counts=[records316,bytes316]
       call header316('complete',12,2_int64,8)
       call write_bytes322(c_loc(counts(1)),16_int64)
       if(compress_close322()/=0_c_int) error stop 'Stage322 compressed audit close or fsync failed'
    end if
    finished316=.true.
  end subroutine
end module audit316_mod
