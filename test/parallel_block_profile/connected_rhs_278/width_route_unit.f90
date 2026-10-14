program width_route_unit
  use iso_fortran_env,only:int64
  use kind_mod,only:dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  type(MPI_Comm) :: context
  type(Rhs_Exchange) :: flexible,fixed(2)
  integer :: me,np,ierr,epoch,state,nc,k,r,tag
  integer :: dk(4,1),ck(4,2),ord(2,1),sp,rp
  integer(int64) :: sent,received,bytes
  logical :: varying
  real(dp) :: actual(3,2),expected(3,2),invalid(3,3)
  character(32) :: mode
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,mode)
  dk(:,1)=[me+1,1,1,1];ck(:,1)=[modulo(me-1,np)+1,1,1,1];ck(:,2)=ck(:,1);ord(:,1)=[1,1]
  do epoch=1,8
     call MPI_Comm_dup(MPI_COMM_WORLD,context,ierr)
     varying=.true.
     if(trim(mode)=='variable-mismatch'.and.me==1) varying=.false.
     call compile_rhs_exchange(flexible,context,int(epoch,int64),2_int64,17,3,2, &
          dk,ord,[1],ck,[2,3],[1,-1],shared_tag=53,variable_columns=varying)
     do nc=1,2
        tag=53+nc
        call compile_rhs_exchange(fixed(nc),context,int(epoch,int64),2_int64,17,3,nc, &
             dk,ord,[1],ck,[2,3],[1,-1],shared_tag=tag)
     end do
     if(trim(mode)=='bad-width') then
        invalid=0
        call exchange_rhs_values(flexible,int(epoch,int64),2_int64,17,1_int64,invalid)
        error stop 'invalid width accepted'
     end if
     do state=1,12
        nc=modulo(state,2)+1
        actual=-999
        do k=1,2
           actual(1,k)=real(100*k+me+epoch+state,dp)
        end do
        expected=actual
        call exchange_rhs_values(flexible,int(epoch,int64),2_int64,17,int(state,int64),actual(:,:nc))
        call exchange_rhs_values(fixed(nc),int(epoch,int64),2_int64,17,int(state,int64),expected(:,:nc))
        if(any(transfer(actual,[0_int64],6)/=transfer(expected,[0_int64],6))) &
             error stop 'variable width differs from fixed route'
        r=modulo(me-1,np)
        do k=1,nc
           if(transfer(actual(2,k),0_int64)/=transfer(real(100*k+r+epoch+state,dp),0_int64)) &
                error stop 'variable width source mismatch'
           if(transfer(actual(3,k),0_int64)/=transfer(-actual(2,k),0_int64)) &
                error stop 'variable width signed alias mismatch'
        end do
        call rhs_exchange_cost(flexible,sent,received,sp,rp,bytes)
        if(np>1) then
           if(sent/=nc.or.received/=nc.or.sp/=1.or.rp/=1) error stop 'variable width cost ledger'
        else
           if(sent/=0.or.received/=0.or.sp/=0.or.rp/=0) error stop 'local width cost ledger'
        end if
     end do
     call release_rhs_exchange(flexible);call release_rhs_exchange(flexible)
     do nc=1,2
        call release_rhs_exchange(fixed(nc))
     end do
     call MPI_Comm_free(context,ierr)
  end do
  if(me==0) print *, 'PASS variable widths, exact fixed-route equivalence and traffic ledger'
  call MPI_Finalize(ierr)
end program
