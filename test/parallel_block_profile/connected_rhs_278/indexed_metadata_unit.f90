program indexed_metadata_unit
  use iso_fortran_env,only:int64
  use kind_mod,only:dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  integer :: me,np,ierr,epoch,owners,ns,nq,nc,i,k,phase,width
  integer(int64) :: topology,remote_values(2),sent,received,bytes
  integer :: sp,rp
  integer,allocatable :: home(:,:),dk(:,:),ck(:,:),ord(:,:),slot(:),target(:),sign(:)
  real(dp),allocatable :: source(:,:),actual(:,:),expected(:,:)
  type(Rhs_Exchange) :: reference
  character(32) :: mode
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,mode)
  owners=max(1,np-1)
  do epoch=1,40
     width=3+modulo(epoch,3);ns=merge(width,0,me<owners)
     nq=9
     if(modulo(epoch,2)==0.and.me==0) nq=0
     nc=1+modulo(epoch,5);phase=8000+epoch;topology=int(epoch,int64)
     allocate(home(2,nq),dk(4,ns),ck(4,nq),ord(2,ns),slot(ns),target(nq),sign(nq),source(ns,nc))
     do i=1,ns
        dk(:,i)=[me+1,i,1,1];slot(i)=i
        do k=1,nc
           source(i,k)=real(100000*epoch+1000*me+10*i+k,dp)
        end do
     end do
     do i=1,nq
        ! Duplicate remote reads and local reads share the same exact semantics.
        home(:,i)=[modulo(me+i/2+epoch,owners),1+modulo(i/2+epoch,width)]
        ck(:,i)=[home(1,i)+1,home(2,i),1,1];target(i)=ns+i
     end do
     if(trim(mode)=='bad-rank'.and.me==0.and.nq>0) home(1,1)=np
     if(trim(mode)=='bad-slot'.and.me==0.and.nq>0) home(2,1)=width+1
     if(trim(mode)=='epoch'.and.me==0) topology=topology+1
     if(trim(mode)=='phase'.and.me==0) phase=phase+1
     if(trim(mode)=='columns'.and.me==0) then
        deallocate(source);allocate(source(ns,nc+1));source=1
     end if
     call deliver_rhs_indexed_metadata(MPI_COMM_WORLD,source,home,actual,topology,phase,remote_values)
     allocate(expected(ns+nq,nc));expected=0;expected(:ns,:)=source;ord=1;sign=1
     call compile_rhs_exchange(reference,MPI_COMM_WORLD,topology,0_int64,phase,ns+nq,nc, &
          dk,ord,slot,ck,target,sign)
     call exchange_rhs_values(reference,topology,0_int64,phase,1_int64,expected)
     if(any(transfer(actual,[0_int64],nq*nc)/=transfer(expected(ns+1:,:),[0_int64],nq*nc))) &
          error stop 'indexed metadata differs from generic writer reference'
     call rhs_exchange_cost(reference,sent,received,sp,rp,bytes)
     if(any(remote_values/=[sent,received])) error stop 'indexed metadata duplicated wire values'
     call release_rhs_exchange(reference)
     deallocate(home,dk,ck,ord,slot,target,sign,source,actual,expected)
  end do
  if(me==0) print *, 'PASS indexed metadata exact reference, changing layouts and bounded lifetime'
  call MPI_Finalize(ierr)
end program
