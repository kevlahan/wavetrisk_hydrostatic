! Diagnostic only: mimic checkpoint metadata + collective payload writes.
! Compare ordering, synchronization, collective/independent writes and a single-writer control.
program checkpoint_io_probe
  use mpi_f08
  use iso_fortran_env, only: int8,int64
  implicit none
  type(MPI_File) :: fh
  type(MPI_Status) :: status
  type(MPI_Info) :: hints,active_hints
  integer,parameter :: domains=160,bytes_per_rank=8192,trials=8,data_pos=24+20*domains
  integer :: ierr,rank,nrank,trial,mode,bad_header,bad_directory,bad_payload,transfer_count,fid,r,ios,q,first_bad,host_length
  integer :: loads(domains),read_loads(domains),library_length
  integer(int64) :: header(3),got(3),offsets(domains),lengths(domains),read_offsets(domains),read_lengths(domains)
  integer(int8) :: payload(bytes_per_rank),read_payload(bytes_per_rank)
  integer(int8),allocatable :: gathered(:,:)
  logical :: found
  character(MPI_MAX_INFO_VAL) :: hint_value
  character(MPI_MAX_PROCESSOR_NAME) :: host
  double precision :: start
  character(100) :: filename
  character(len=MPI_MAX_LIBRARY_VERSION_STRING) :: library_version
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,rank,ierr)
  call MPI_Comm_size(MPI_COMM_WORLD,nrank,ierr)
  if(nrank<2)call MPI_Abort(MPI_COMM_WORLD,1,ierr)
  call MPI_Get_library_version(library_version,library_length,ierr)
  if(rank==0)then
    write(*,'(a)')trim(library_version(:library_length))
    write(*,'(a,i0)')'MPI ranks = ',nrank
  end if
  header=[int(z'5741564554524953',int64),1_int64,int(domains,int64)]
  loads=1
  offsets=[(int(data_pos+16*r,int64),r=0,domains-1)]
  lengths=16_int64
  payload=int(mod(rank,100)+1,int8)
  allocate(gathered(bytes_per_rank,nrank))
  call MPI_Get_processor_name(host,host_length,ierr)
  if(rank==0.or.rank==nrank-1)write(*,'(a,i0,2a)')'reader rank=',rank,' host=',trim(host)
  do mode=0,5
    bad_header=0;bad_directory=0;bad_payload=0
    do trial=1,trials
      write(filename,'(a,i0,a,i0,a)')'probe-',mode,'-',trial,'.bin'
      if(mode==5)then
        ! Small diagnostic control only; not a proposal to gather climate checkpoints.
        call MPI_Gather(payload,bytes_per_rank,MPI_BYTE,gathered,bytes_per_rank,MPI_BYTE,0,MPI_COMM_WORLD,ierr)
        call check('gather single-writer control')
        if(rank==0)then
          open(newunit=fid,file=trim(filename),access='stream',form='unformatted',status='replace',iostat=ios)
          if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,7,ierr)
          write(fid)header,loads,offsets,lengths,gathered
          close(fid)
        end if
      else
      hints=MPI_INFO_NULL
      if(mode==3)then
        call MPI_Info_create(hints,ierr)
        call check('create hints')
        call MPI_Info_set(hints,'romio_cb_write','disable',ierr)
        call check('disable collective buffering hint')
        call MPI_Info_set(hints,'romio_ds_write','disable',ierr)
        call check('disable data sieving hint')
      end if
      call MPI_File_open(MPI_COMM_WORLD,trim(filename),MPI_MODE_CREATE+MPI_MODE_WRONLY,hints,fh,ierr)
      call check('open')
      if(mode==3)then
        call MPI_Info_free(hints,ierr)
        if(trial==1)then
          call MPI_File_get_info(fh,active_hints,ierr)
          call check('get active hints')
          hint_value='not reported'
          call MPI_Info_get(active_hints,'romio_cb_write',MPI_MAX_INFO_VAL,hint_value,found,ierr)
          call check('get cb hint')
          if(rank==0)write(*,'(a,l1,2a)')'romio_cb_write reported=',found,' value=',trim(hint_value)
          hint_value='not reported'
          call MPI_Info_get(active_hints,'romio_ds_write',MPI_MAX_INFO_VAL,hint_value,found,ierr)
          call check('get ds hint')
          if(rank==0)write(*,'(a,l1,2a)')'romio_ds_write reported=',found,' value=',trim(hint_value)
          call MPI_Info_free(active_hints,ierr)
        end if
      end if
      ! Let rank zero arrive first; the collective itself may synchronize on
      ! some MPI/filesystem implementations. No failure is forced or assumed.
      if(rank==1+mod(trial-1,nrank-1))then
        start=MPI_Wtime()
        do while(MPI_Wtime()-start<0.10d0)
        end do
      end if
      call MPI_File_set_size(fh,0_MPI_OFFSET_KIND,ierr)
      call check('set_size')
      if(mode>=1)then
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call check('post-truncation barrier')
      end if
      if(rank==0)then
        call MPI_File_write_at(fh,0_MPI_OFFSET_KIND,header,3,MPI_INTEGER8,status,ierr)
        call transferred('write header',MPI_INTEGER8,3)
        call MPI_File_write_at(fh,24_MPI_OFFSET_KIND,loads,domains,MPI_INTEGER,status,ierr)
        call transferred('write loads',MPI_INTEGER,domains)
        call MPI_File_write_at(fh,int(24+4*domains,MPI_OFFSET_KIND),offsets,domains,MPI_INTEGER8,status,ierr)
        call transferred('write offsets',MPI_INTEGER8,domains)
        call MPI_File_write_at(fh,int(24+12*domains,MPI_OFFSET_KIND),lengths,domains,MPI_INTEGER8,status,ierr)
        call transferred('write lengths',MPI_INTEGER8,domains)
      end if
      if(mode>=2)then
        call MPI_File_sync(fh,ierr)
        call check('sync metadata')
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call check('pre-payload barrier')
      if(mode==4)then
        call MPI_File_write_at(fh,int(data_pos+bytes_per_rank*rank,MPI_OFFSET_KIND), &
             payload,bytes_per_rank,MPI_BYTE,status,ierr)
      else
        call MPI_File_write_at_all(fh,int(data_pos+bytes_per_rank*rank,MPI_OFFSET_KIND), &
             payload,bytes_per_rank,MPI_BYTE,status,ierr)
      end if
      call transferred('write payload',MPI_BYTE,bytes_per_rank)
      if(mode>=2)then
        call MPI_File_sync(fh,ierr)
        call check('sync payload')
        call MPI_Barrier(MPI_COMM_WORLD,ierr)
        call check('post-sync barrier')
      end if
      call MPI_File_close(fh,ierr)
      call check('close')
      end if ! MPI writes versus single-writer control
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
      call check('pre-readback barrier')
      if(rank==0.or.rank==nrank-1)then
        open(newunit=fid,file=trim(filename),access='stream',form='unformatted',status='old',iostat=ios)
        if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,2,ierr)
        read(fid,iostat=ios)got,read_loads,read_offsets,read_lengths
        if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,3,ierr)
        if(any(got/=header))bad_header=bad_header+1
        if(any(read_loads/=loads).or.any(read_offsets/=offsets).or.any(read_lengths/=lengths)) &
             bad_directory=bad_directory+1
        do r=0,nrank-1
          read(fid,pos=data_pos+bytes_per_rank*r+1,iostat=ios)read_payload
          if(ios/=0)call MPI_Abort(MPI_COMM_WORLD,4,ierr)
          if(any(read_payload/=int(mod(r,100)+1,int8)))then
            bad_payload=bad_payload+1
            first_bad=0
            do q=1,bytes_per_rank
              if(read_payload(q)==int(mod(r,100)+1,int8))cycle
              if(first_bad==0)first_bad=q
            end do
            write(*,'(a,8(i0,1x))')'BAD payload: mode trial reader block offset count expected observed = ', &
                 mode,trial,rank,r,data_pos+bytes_per_rank*r+first_bad-1, &
                 count(read_payload/=int(mod(r,100)+1,int8)),mod(r,100)+1,int(read_payload(first_bad))
          end if
        end do
        close(fid)
        if(any(got/=header))write(*,'(a,3(i0,1x),a,3(z16.16,1x))') &
             'BAD header: mode trial reader = ',mode,trial,rank,' values=',got
      end if
      call MPI_Barrier(MPI_COMM_WORLD,ierr)
    end do
    if(rank==0.or.rank==nrank-1)then
      write(*,'(a,i0,a,i0,a,3(i0,1x))')'mode=',mode,' reader=',rank, &
           ' bad header/directory/payload = ',bad_header,bad_directory,bad_payload
      flush(6)
    end if
  end do
  if(rank==0)then
    write(*,'(a)')'Modes: 0=current; 1=barrier; 2=barrier+sync; 3=sync+hints; 4=sync+independent; 5=single writer.'
    write(*,'(a)')'Probe complete. Inspect verification counts; completion alone is not a passing result.'
  end if
  call MPI_Finalize(ierr)
contains
  subroutine check(context)
    character(*),intent(in) :: context
    if(ierr/=MPI_SUCCESS)then
      write(*,'(a,i0,2a,i0)')'rank ',rank,': '//context//' MPI error ',': ',ierr
      call MPI_Abort(MPI_COMM_WORLD,5,ierr)
    end if
  end subroutine check
  subroutine transferred(context,datatype,expected)
    character(*),intent(in) :: context
    type(MPI_Datatype),intent(in) :: datatype
    integer,intent(in) :: expected
    call check(context)
    call MPI_Get_count(status,datatype,transfer_count,ierr)
    call check('get_count')
    if(transfer_count/=expected)then
      write(*,'(a,i0,2a,2i12)')'rank ',rank,': '//context,' expected/actual=',expected,transfer_count
      call MPI_Abort(MPI_COMM_WORLD,6,ierr)
    end if
  end subroutine transferred
end program checkpoint_io_probe
