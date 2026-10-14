program exchange_stress
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  use parallel_block_rhs_exchange_mod
  implicit none
  integer, parameter :: n=8192,nlayer=30
  type(Rhs_Exchange) :: plan
  integer :: me,np,ierr,i,k,state,previous,ns,nr
  integer :: dk(4,n),orders(2,n),ds(n),ck(4,n),cs(n),signs(n)
  integer(int64) :: sv,rv,bytes,headers
  real(dp), allocatable :: values(:,:)
  real(dp) :: expected
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,me,ierr);call MPI_Comm_size(MPI_COMM_WORLD,np,ierr)
  if (np/=4) error stop 'stress requires four ranks'
  previous=modulo(me-1,np)
  allocate(values(2*n,nlayer))
  do i=1,n
     dk(:,i)=[me+1,1,i,1];ck(:,i)=[previous+1,1,i,1]
     orders(:,i)=[1,i];ds(i)=i;cs(i)=n+i;signs(i)=merge(1,-1,mod(i,2)==0)
  end do
  call compile_rhs_exchange(plan,MPI_COMM_WORLD,91_int64,73_int64,4,2*n,nlayer,dk,orders,ds,ck,cs,signs)
  call rhs_exchange_cost(plan,sv,rv,ns,nr,bytes,headers)
  if (sv/=int(n,int64)*nlayer.or.rv/=sv.or.ns/=1.or.nr/=1.or.headers/=8) error stop 'stress traffic differs'
  do state=1,3
     values=-999.0_dp
     do k=1,nlayer
        do i=1,n
           values(i,k)=real(state*10000000+(me+1)*1000000+k*10000+i,dp)+0.125_dp
        end do
     end do
     call exchange_rhs_values(plan,91_int64,73_int64,4,int(state,int64),values)
     do k=1,nlayer
        do i=1,n
           expected=signs(i)*(real(state*10000000+(previous+1)*1000000+k*10000+i,dp)+0.125_dp)
           if (transfer(values(n+i,k),0_int64)/=transfer(expected,0_int64)) error stop 'rendezvous payload differs'
        end do
     end do
  end do
  call release_rhs_exchange(plan)
  if (me==0) print *, 'PASS RHS exchange rendezvous: 4 ranks, 3 states, 8192 values x 30 layers per peer'
  call MPI_Finalize(ierr)
end program
