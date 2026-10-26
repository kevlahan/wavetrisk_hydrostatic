module audit316_mod
  ! Qualification records: exact integer metadata and finite real SHA-256 digests.
  ! No payload replay, collective, clock or numerical substitution is used.
  use iso_fortran_env, only: int8,int32,int64,real64
  use iso_c_binding, only: c_char,c_int,c_int64_t,c_ptr,c_loc,c_null_char,c_null_ptr
  implicit none
  private
  public :: audit316_init,audit316_finish,audit316_enabled,audit316_begin,audit316_end
  public :: audit316_i,audit316_i64,audit316_r,audit316_l
  integer :: depth316=0
  integer(int64) :: records316=0_int64,bytes316=0_int64
  integer(int64),parameter :: limit316=128_int64*1024_int64*1024_int64
  integer,parameter :: max_depth316=32,max_tag316=256
  character(max_tag316) :: scopes316(max_depth316)
  integer(int64) :: epochs316(max_depth316)
  logical :: initialized316=.false.,enabled316=.false.,finished316=.false.
  ! Stage325 qualification ledger: counters only, never a timing instrument.
  public :: work325_enabled,work325_record
  integer,parameter :: work325_limit=4096,work325_width=24
  integer :: work325_n=0,work325_rank=-1
  logical :: work325_on=.false.
  type :: Work325_Entry
     character(48) :: tag=''
     integer :: phase=0,width=0
     integer(int64) :: calls=0,total(work325_width)=0,peak(work325_width)=0
     integer(int64),allocatable :: sent(:),received(:),send_peak(:),receive_peak(:)
  end type
  type(Work325_Entry),allocatable :: work325(:)
  interface
     integer(c_int) function digest_write325(address,length) bind(C,name='audit_digest325_write')
       import :: c_int,c_ptr,c_int64_t
       type(c_ptr),value :: address
       integer(c_int64_t),value :: length
     end function
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
  subroutine work325_init(rank)
    integer,intent(in) :: rank
    character(16) :: option
    integer :: status
    option='';call get_environment_variable('WAVETRISK_WORK_LEDGER325',option,status=status)
    if(status/=0.and.status/=1) error stop 'Stage325 ledger environment read'
    if(trim(option)/=''.and.trim(option)/='0'.and.trim(option)/='1') error stop 'Stage325 ledger requires 0 or 1'
    work325_on=trim(option)=='1';work325_rank=rank
    if(work325_on) allocate(work325(work325_limit))
  end subroutine

  logical function work325_enabled() result(on)
    on=work325_on.and..not.finished316
  end function

  subroutine work325_record(tag,phase,values,sent,received)
    character(*),intent(in) :: tag
    integer,intent(in) :: phase
    integer(int64),intent(in) :: values(:)
    integer(int64),optional,intent(in) :: sent(:),received(:)
    integer :: i,n,k
    if(.not.work325_enabled()) return
    n=size(values)
    if(n<1.or.n>work325_width.or.len_trim(tag)>48.or.len_trim(tag)<1.or.phase<0) &
         error stop 'Stage325 ledger descriptor'
    if(any(values<0)) error stop 'Stage325 negative work count'
    k=0
    do i=1,work325_n
       if(work325(i)%tag==tag.and.work325(i)%phase==phase) then
          k=i;exit
       end if
    end do
    if(k==0) then
       if(work325_n==work325_limit) error stop 'Stage325 ledger entry bound'
       work325_n=work325_n+1;k=work325_n
       work325(k)%tag=tag;work325(k)%phase=phase;work325(k)%width=n
    end if
    associate(w=>work325(k))
    if(w%width/=n) error stop 'Stage325 changed ledger width'
    if(any(values>huge(0_int64)-w%total(:n))) error stop 'Stage325 ledger overflow'
    w%calls=w%calls+1;w%total(:n)=w%total(:n)+values;w%peak(:n)=max(w%peak(:n),values)
    if(present(sent).neqv.present(received)) error stop 'Stage325 incomplete peer ledger'
    if(present(sent).and.present(received)) then
       if(size(sent)/=size(received).or.size(sent)>1024.or.size(sent)<1) error stop 'Stage325 peer extent'
       if(any(sent<0).or.any(received<0)) error stop 'Stage325 negative peer count'
       if(.not.allocated(w%sent)) then
          if(w%calls/=1) error stop 'Stage325 late peer ledger'
          allocate(w%sent(size(sent)),w%received(size(sent)),w%send_peak(size(sent)),w%receive_peak(size(sent)))
          w%sent=0;w%received=0;w%send_peak=0;w%receive_peak=0
       end if
       if(size(w%sent)/=size(sent)) error stop 'Stage325 peer extent changed'
       if(any(sent>huge(0_int64)-w%sent).or.any(received>huge(0_int64)-w%received)) &
            error stop 'Stage325 peer overflow'
       w%sent=w%sent+sent;w%received=w%received+received
       w%send_peak=max(w%send_peak,sent);w%receive_peak=max(w%receive_peak,received)
    else if(allocated(w%sent)) then
       error stop 'Stage325 missing peer ledger'
    end if
    end associate
  end subroutine

  subroutine work325_finish()
    integer :: unit,status,i,p,n
    integer(int64) :: extent,rows
    character(80) :: filename
    character(1200) :: line
    if(.not.work325_on) return
    write(filename,'(a,i0,a)') 'work325-rank-',work325_rank,'.tsv'
    open(newunit=unit,file=trim(filename),status='new',action='write',iostat=status)
    if(status/=0) error stop 'Stage325 ledger open'
    extent=0
    write(line,'(a,1x,i0)') 'RK325WORK1',work325_rank
    call work325_line(unit,line,extent)
    rows=0
    do i=1,work325_n
       associate(w=>work325(i))
       n=w%width
       write(line,'(a,1x,a,1x,*(i0,1x))') 'C',trim(w%tag),w%phase,w%calls,n,w%total(:n),w%peak(:n)
       call work325_line(unit,line,extent)
       rows=rows+1
       if(allocated(w%sent)) then
          do p=1,size(w%sent)
             if(w%sent(p)==0.and.w%received(p)==0) cycle
             write(line,'(a,1x,a,1x,*(i0,1x))') 'P',trim(w%tag),w%phase,p-1, &
                  w%sent(p),w%received(p),w%send_peak(p),w%receive_peak(p)
             call work325_line(unit,line,extent)
             rows=rows+1
          end do
       end if
       end associate
    end do
    write(line,'(a,1x,*(i0,1x))') 'COMPLETE',work325_n,rows
    call work325_line(unit,line,extent)
    flush(unit,iostat=status)
    if(status/=0) error stop 'Stage325 ledger flush'
    inquire(unit=unit,size=extent,iostat=status)
    if(status/=0.or.extent>8_int64*1024*1024) error stop 'Stage325 ledger 8MiB bound'
    close(unit,iostat=status)
    if(status/=0) error stop 'Stage325 ledger close'
    deallocate(work325)
  end subroutine

  subroutine work325_line(unit,line,extent)
    integer,intent(in) :: unit
    character(*),intent(in) :: line
    integer(int64),intent(inout) :: extent
    integer :: status
    extent=extent+int(len_trim(line),int64)+1_int64
    if(extent>8_int64*1024*1024) error stop 'Stage325 ledger 8MiB bound'
    write(unit,'(a)',iostat=status) trim(line)
    if(status/=0) error stop 'Stage325 ledger write'
  end subroutine
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
    call work325_init(rank)
    option='';call get_environment_variable('WAVETRISK_EXECUTION_AUDIT316',option,status=status)
    if(status/=0.and.status/=1) error stop 'Stage316 audit environment read failed'
    if(trim(option)/=''.and.trim(option)/='0'.and.trim(option)/='1') &
         error stop 'Stage316 audit requires 0 or 1'
    enabled316=trim(option)=='1'
    if(.not.enabled316) return
    option='';call get_environment_variable('WAVETRISK_COMPACT_AUDIT325',option,status=status)
    if(status/=0.or.trim(option)/='1') error stop 'Stage325 audit requires compact mode1'
    if(rank<0.or.storage_size(0)/=32.or.storage_size(0_int64)/=64.or.storage_size(0.0_real64)/=64.or. &
         storage_size('a')/=8.or.kind('a')/=c_char) error stop 'Stage316 unsupported integer or real representation'
    write(filename,'(a,i0,a)') 'execution316-rank-',rank,'.bin.gz'
    if(compress_open322(trim(filename)//c_null_char)/=0_c_int) error stop 'Stage322 compressed audit open failed'
    magic='RK325V01'
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
    if(kind==13) added=16_int64+k+32_int64
    if(bytes316>limit316-added) error stop 'Stage325 compact audit file exceeds 128MiB bound'
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
    integer(c_int) :: status
    if(.not.audit316_enabled()) return
    ! Kind13 preserves the original logical element count; payload is32 bytes.
    call header316(name,13,int(size(value),int64),0)
    if(size(value)==0) then
       status=digest_write325(c_null_ptr,0_c_int64_t)
    else if(is_contiguous(value)) then
       status=digest_write325(c_loc(value(1)),8_c_int64_t*size(value))
    else
       allocate(packed(size(value)));packed=value
       status=digest_write325(c_loc(packed(1)),8_c_int64_t*size(packed))
    end if
    if(status/=0_c_int) error stop 'Stage325 nonfinite real or digest write failed'
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
    call work325_finish()
    finished316=.true.
  end subroutine
end module audit316_mod
