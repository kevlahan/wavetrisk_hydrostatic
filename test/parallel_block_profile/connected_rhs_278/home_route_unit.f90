program home_route_unit
  use iso_fortran_env,only:int64
  use kind_mod,only:dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  type(MPI_Comm) :: context
  type(Rhs_Owner_Directory) :: directory
  type(Rhs_Exchange) :: home,reference
  integer :: me,np,ierr,owners,owner,previous,next,epoch,state,nc,k,nd,pass,i
  integer,allocatable :: owned_key(:,:),canonical(:),dh(:,:),ch(:,:)
  integer :: dk(4,6),ck(4,6),ord(2,6),ds(6),orientation(6),target(6)
  integer :: sp,rp,sp_ref,rp_ref
  integer(int64) :: sent,received,bytes,sent_ref,received_ref,bytes_ref
  real(dp) :: actual(10,2),expected(10,2)
  character(32) :: mode
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,mode)
  owners=max(1,np-1);owner=modulo(me,owners);previous=modulo(owner-1,owners);next=modulo(owner+1,owners)
  target=[1,2,3,4,5,6];orientation=[1,1,1,-1,1,-1]
  if(me<owners) then
     canonical=[1,2];allocate(owned_key(4,2))
     owned_key(:,1)=key(me,1);owned_key(:,2)=key(me,2)
  else
     allocate(canonical(0),owned_key(4,0))
  end if
  if(trim(mode)=='duplicate-owner'.and.me==0) owned_key(:,2)=owned_key(:,1)
  do epoch=1,6
     call MPI_Comm_dup(MPI_COMM_WORLD,context,ierr)
     call compile_rhs_owner_directory(directory,context,owned_key,canonical)
     ck(:,1)=key(owner,2);ck(:,2)=key(owner,1);ck(:,3)=key(previous,1)
     ck(:,4)=ck(:,3);ck(:,5)=key(next,2);ck(:,6)=ck(:,2)
     if(trim(mode)=='missing-owner') ck(:,6)=key(np+1,1)
     call lookup_rhs_owner_directory(directory,ck,ch)
     if(any(ch(1,:)/=[owner,owner,previous,previous,next,owner]).or.any(ch(2,:)/=[2,1,1,1,2,1])) &
          error stop 'wrong canonical binding'
     do pass=1,2
        nd=0
        if(me<owners) then
           call define(me,1,1,[1,1]);call define(me,2,2,[1,1])
        end if
        if(pass==2) then
           call define(previous,1,7,[2,me+1]);call define(previous,1,8,[3,me+1])
           if(me<owners) call define(me,2,9,[4,1])
        end if
        if(trim(mode)=='ambiguous'.and.me==0) call define(0,1,10,[1,1])
        call lookup_rhs_owner_directory(directory,dk(:,:nd),dh)
        if(trim(mode)=='bad-home'.and.me==0) dh(2,1)=10
        if(trim(mode)=='missing-producer') then
           nd=0;deallocate(dh);allocate(dh(2,0))
        end if
        call compile_rhs_exchange(home,context,int(epoch,int64),2_int64,17,10,2, &
             dk(:,:nd),ord(:,:nd),ds(:nd),ck,target,orientation,shared_tag=53,variable_columns=.true., &
             definition_home=dh,consumer_home=ch,canonical_slot=canonical)
        call compile_rhs_exchange(reference,context,int(epoch,int64),2_int64,17,10,2, &
             dk(:,:nd),ord(:,:nd),ds(:nd),ck,target,orientation,shared_tag=54,variable_columns=.true.)
        do state=1,12
           nc=modulo(state,2)+1
           do k=1,2
              do i=1,10
                 actual(i,k)=real(10000*k+1000*me+100*epoch+10*state+i,dp)
              end do
           end do
           expected=actual
           call exchange_rhs_values(home,int(epoch,int64),2_int64,17,int(state,int64),actual(:,:nc))
           call exchange_rhs_values(reference,int(epoch,int64),2_int64,17,int(state,int64),expected(:,:nc))
           if(any(transfer(actual,[0_int64],20)/=transfer(expected,[0_int64],20))) &
                error stop 'home route differs from hash reference'
           call rhs_exchange_cost(home,sent,received,sp,rp,bytes)
           call rhs_exchange_cost(reference,sent_ref,received_ref,sp_ref,rp_ref,bytes_ref)
           if(sent/=sent_ref.or.received/=received_ref.or.sp/=sp_ref.or.rp/=rp_ref.or.bytes/=bytes_ref) &
                error stop 'home route differs from reference cost'
        end do
        call release_rhs_exchange(home);call release_rhs_exchange(home);call release_rhs_exchange(reference)
     end do
     call MPI_Comm_free(context,ierr)
  end do
  if(me==0) print *, 'PASS home routes, exact hash equivalence and traffic ledger'
  call MPI_Finalize(ierr)
contains
  function key(rank,slot) result(value)
    integer,intent(in) :: rank,slot
    integer :: value(4)
    value=[rank+1,slot,1,1]
  end function
  subroutine define(rank,slot,address,order)
    integer,intent(in) :: rank,slot,address,order(2)
    nd=nd+1;dk(:,nd)=key(rank,slot);ord(:,nd)=order;ds(nd)=address
  end subroutine
end program
