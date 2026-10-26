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
