program shared_route_unit
  use iso_fortran_env,only:int64
  use kind_mod,only:dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  type(MPI_Comm) :: context
  type(Rhs_Exchange) :: a,b
  integer :: me,np,ierr,epoch,state,n,tag
  integer :: dk(4,1),ck(4,1),order(2,1)
  real(dp) :: va(2,2),vb(2,2)
  character(32) :: mode
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  call get_command_argument(1,mode)
  dk(:,1)=[me+1,1,1,1];ck(:,1)=[modulo(me-1,np)+1,1,1,1];order(:,1)=[1,1]
  do epoch=1,8
     call MPI_Comm_dup(MPI_COMM_WORLD,context,ierr)
     tag=53
     if(trim(mode)=='tag-mismatch'.and.me==1) tag=54
     if(trim(mode)=='tag-range') tag=32768
     call compile_rhs_exchange(a,context,int(epoch,int64),2_int64,17,2,2, &
          dk,order,[1],ck,[2],[1],shared_tag=tag)
     call compile_rhs_exchange(b,context,int(epoch,int64),2_int64,18,2,2, &
          dk,order,[1],ck,[2],[-1],shared_tag=56)
     do state=1,5
        va(1,:)=[real(me+epoch+state,dp),real(100+me+epoch+state,dp)];va(2,:)=-999
        vb(1,:)=va(1,:)*2;vb(2,:)=-999
        call exchange_rhs_values(a,int(epoch,int64),2_int64,17,int(state,int64),va)
        call exchange_rhs_values(b,int(epoch,int64),2_int64,18,int(state,int64),vb)
        n=modulo(me-1,np)+epoch+state
        if(any(transfer(va(2,:),[0_int64,0_int64])/=transfer([real(n,dp),real(100+n,dp)],[0_int64,0_int64]))) &
             error stop 'shared route payload mismatch'
        if(any(transfer(vb(2,:),[0_int64,0_int64])/=transfer(-2*va(2,:),[0_int64,0_int64]))) &
             error stop 'shared route payload/tag mismatch'
     end do
     call release_rhs_exchange(a);call release_rhs_exchange(a)
     ! Releasing a borrowing route must neither free the context nor disrupt b.
     call MPI_Comm_size(context,n,ierr)
     if(ierr/=MPI_SUCCESS.or.n/=np) error stop 'borrowed context lifetime'
     vb(1,:)=[real(me,dp),real(me+1,dp)]
     call exchange_rhs_values(b,int(epoch,int64),2_int64,18,6_int64,vb)
     if(any(transfer(vb(2,:),[0_int64,0_int64])/= &
          transfer(-[real(modulo(me-1,np),dp),real(modulo(me-1,np)+1,dp)],[0_int64,0_int64]))) &
          error stop 'surviving shared route failed'
     call release_rhs_exchange(b)
     call MPI_Comm_free(context,ierr)
     if(ierr/=MPI_SUCCESS) error stop 'shared context teardown'
  end do
  if(me==0) print *, 'PASS shared routes, distinct tags, repeated epochs and lifetime'
  call MPI_Finalize(ierr)
end program
