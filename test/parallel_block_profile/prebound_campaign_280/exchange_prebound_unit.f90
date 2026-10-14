program test_exchange
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  type(Rhs_Exchange) :: plan
  integer :: me,np,ierr,i,k,j,nc,nd,pass,ns,nr
  integer, allocatable :: dk(:,:),priority(:,:),ds(:),ck(:,:),cs(:),signs(:),demand(:),provider(:,:)
  integer(int64) :: sv,rv,bytes,headers,mask_epoch,state_epoch,totals(4),global(4)
  real(dp),allocatable :: scratch(:)
  real(dp) :: values(8,3),expected(8,3),a,b
  character(40) :: mode
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,mode)
  if (np/=1.and.np/=4) error stop 'test uses one or four ranks'
  do pass=1,3
     nd=0;nc=0
     if (me<2) nd=2
     if (me<3) nc=6
     allocate(dk(4,nd),priority(2,nd),ds(nd),ck(4,nc),cs(nc),signs(nc))
     if (nd>0) then
        dk(:,1)=[1,1,1,1];dk(:,2)=[1,1,1,2];ds=[1,2]
        ! The same geometric node's two components have DIFFERENT final writers.
        if (me==0) then
           priority(:,1)=[1,2];priority(:,2)=[2,1]
        else
           priority(:,1)=[2,1];priority(:,2)=[1,2]
        end if
     end if
     if (nc>0) then
        ck(:,1)=[1,1,1,1];ck(:,2)=ck(:,1);ck(:,3)=[1,1,1,2]
        ck(:,4)=ck(:,3);ck(:,5)=ck(:,1);ck(:,6)=ck(:,3)
        cs=[3,4,5,1,2,6];signs=[1,-1,1,1,1,-1]
     end if
     if (trim(mode)=='missing'.and.nc>0) ck(3,1)=99
     if (trim(mode)=='ambiguous'.and.me==1) priority(:,1)=[1,2]
     if (trim(mode)=='destination'.and.nc>0) cs(6)=cs(5)
     if (trim(mode)=='orientation'.and.nc>0) signs(1)=0
     if (trim(mode)=='identity'.and.nd>0) dk(1,1)=0
     if (trim(mode)=='extent'.and.nd>0) ds(1)=9
     ! Reverse definition and consumer presentation order in the second plan.
     ! Only the supplied numerical priorities may decide the final writer.
     if (pass==2) then
        dk=dk(:,nd:1:-1);priority=priority(:,nd:1:-1);ds=ds(nd:1:-1)
        ck=ck(:,nc:1:-1);cs=cs(nc:1:-1);signs=signs(nc:1:-1)
     end if
     allocate(provider(2,nc))
     do i=1,nc
        if(ck(4,i)==1) then
           provider(:,i)=[min(1,np-1),1]
        else
           provider(:,i)=[0,2]
        end if
     end do
     if(trim(mode)=='provider-rank'.and.nc>0) provider(1,1)=np
     if(trim(mode)=='provider-slot'.and.nc>0) provider(2,1)=9
     if(trim(mode)=='binding-mode'.and.me==1) then
        call compile_rhs_exchange(plan,MPI_COMM_WORLD,int(pass,int64),int(2*pass,int64),17,8,3, &
             dk,priority,ds,ck,cs,signs,shared_buffers=.true.)
     end if
     if(trim(mode)=='mixed'.and.nd>0) then
        call compile_rhs_exchange(plan,MPI_COMM_WORLD,int(pass,int64),int(2*pass,int64),17,8,3, &
             dk,priority,ds,ck,cs,signs,shared_buffers=.true.,consumer_provider=provider)
     end if
     mask_epoch=int(2*pass,int64)
     if (trim(mode)=='distributed-mask'.and.me==1) mask_epoch=mask_epoch+1
     call compile_rhs_exchange(plan,MPI_COMM_WORLD,int(pass,int64),mask_epoch,17,8,3, &
          dk(:,:0),priority(:,:0),ds(:0),ck,cs,signs,shared_buffers=.true.,consumer_provider=provider)
     if (trim(mode)=='recompile') call compile_rhs_exchange(plan,MPI_COMM_WORLD,9_int64,9_int64,17,8,3, &
          dk(:,:0),priority(:,:0),ds(:0),ck,cs,signs,shared_buffers=.true.,consumer_provider=provider)
     if(.not.allocated(scratch)) allocate(scratch(max(rhs_exchange_buffer_values(plan),6+2*np)))
     if(trim(mode)=='scratch') then
        deallocate(scratch);allocate(scratch(0))
     end if
     call rhs_exchange_source_slots(plan,demand)
     if(any(demand<1).or.any(demand>2)) error stop 'provider demand includes unused slots'
     if((me<2).neqv.(size(demand)>0)) error stop 'provider demand owner differs'
     call rhs_exchange_cost(plan,sv,rv,ns,nr,bytes,headers)
     if (headers/=8_int64*ns) error stop 'state header cost differs'
     if (sv<0.or.rv<0.or.ns<0.or.nr<0.or.bytes<0) error stop 'negative cost'
     totals=[sv,rv,int(ns,int64),int(nr,int64)]
     call MPI_Allreduce(totals,global,4,MPI_INTEGER8,MPI_SUM,MPI_COMM_WORLD,ierr)
     if (global(1)/=global(2).or.global(3)/=global(4)) error stop 'traffic totals differ'
     if (np==4) then
        ! Unique A goes rank1->rank0/rank2 and B rank0->rank1/rank2.
        ! Three layers per edge, four directed peers; six aliases per consumer.
        if (any(global/=[12_int64,12_int64,4_int64,4_int64])) error stop 'request deduplication differs'
     else
        if (any(global/=0)) error stop 'local aliases used MPI payload'
     end if
     do j=1,5
        do k=1,3
           do i=1,8
              values(i,k)=value(me,i,k,j,pass)
           end do
        end do
        expected=values
        if (nc>0) then
           do k=1,3
              a=value(min(1,np-1),1,k,j,pass);b=value(0,2,k,j,pass)
              expected(1,k)=b;expected(2,k)=a;expected(3,k)=a
              expected(4,k)=-a;expected(5,k)=b;expected(6,k)=-b
           end do
        end if
        if (trim(mode)=='topology') call exchange_rhs_values(plan,99_int64,int(2*pass,int64),17,int(j,int64),values,scratch)
        if (trim(mode)=='mask') call exchange_rhs_values(plan,int(pass,int64),99_int64,17,int(j,int64),values,scratch)
        if (trim(mode)=='phase') call exchange_rhs_values(plan,int(pass,int64),int(2*pass,int64),99,int(j,int64),values,scratch)
        if (trim(mode)=='state') call exchange_rhs_values(plan,int(pass,int64),int(2*pass,int64),17,-1_int64,values,scratch)
        state_epoch=int(j,int64)
        if (trim(mode)=='peer-state'.and.me==1) state_epoch=state_epoch+1
        call exchange_rhs_values(plan,int(pass,int64),int(2*pass,int64),17,state_epoch,values,scratch)
        if (trim(mode)=='duplicate') call exchange_rhs_values(plan,int(pass,int64),int(2*pass,int64),17,int(j,int64),values,scratch)
        if (any(transfer(values,[0_int64],24)/=transfer(expected,[0_int64],24))) &
             error stop 'final writer, orientation, local cycle or state differs'
     end do
     call release_rhs_exchange(plan)
     ! Release is idempotent, including an all-empty rank's plan.
     call release_rhs_exchange(plan)
     deallocate(provider,dk,priority,ds,ck,cs,signs)
  end do
  ! Identity copies are removed, but an in-place sign reversal remains a
  ! real consumer and must still declare its provider demand.
  do pass=1,2
     call compile_rhs_exchange(plan,MPI_COMM_WORLD,int(90+pass,int64),90_int64,17,8,3, &
          reshape([integer::],[4,0]),reshape([integer::],[2,0]),[integer::], &
          reshape([me+1,1,1,1],[4,1]),[1],[merge(1,-1,pass==1)],shared_buffers=.true.,consumer_provider=reshape([me,1],[2,1]))
     if(.not.allocated(scratch)) allocate(scratch(max(rhs_exchange_buffer_values(plan),6+2*np)))
     if(trim(mode)=='scratch') then
        deallocate(scratch);allocate(scratch(0))
     end if
     call rhs_exchange_source_slots(plan,demand)
     if(size(demand)/=pass-1) error stop 'identity pruning provider demand differs'
     values=2.0_dp;expected=values
     if(pass==2) expected(1,:)=-2.0_dp
     call exchange_rhs_values(plan,int(90+pass,int64),90_int64,17,1_int64,values,scratch)
     if(any(transfer(values,[0_int64],24)/=transfer(expected,[0_int64],24))) &
          error stop 'identity/sign pruning differs'
     call release_rhs_exchange(plan)
  end do
  ! All ranks can also have empty producer, consumer and numerical storage.
  ! This checks zero-count collectives during compilation and Waitall(0).
  call compile_rhs_exchange(plan,MPI_COMM_WORLD,99_int64,99_int64,17,0,3, &
       reshape([integer::],[4,0]),reshape([integer::],[2,0]),[integer::], &
       reshape([integer::],[4,0]),[integer::],[integer::],shared_buffers=.true., &
       consumer_provider=reshape([integer::],[2,0]))
  call exchange_rhs_values(plan,99_int64,99_int64,17,1_int64,values(:0,:),scratch)
  call rhs_exchange_cost(plan,sv,rv,ns,nr,bytes,headers)
  if (sv/=0.or.rv/=0.or.ns/=0.or.nr/=0.or.headers/=0) error stop 'empty plan traffic'
  call release_rhs_exchange(plan)
  if (me==0) print *, 'PASS RHS exchange: 3 epochs, 15 transactions, exact aliases, bounded deduplicated traffic'
  call MPI_Finalize(ierr)
contains
  real(dp) function value(rank,slot,layer,state,epoch)
    integer, intent(in) :: rank,slot,layer,state,epoch
    value=real(100000*epoch+10000*state+1000*rank+100*slot+layer,dp)/8
  end function
end program
