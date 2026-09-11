program benchmark
use kind_mod, only: dp
use parallel_block_scalar_storage_mod
use iso_fortran_env, only: int64
implicit none
type(Scalar_Record_Storage)::store
integer,parameter::span=128,nk=41,ns=2,groups=64,repeats=1
integer,parameter::order(8)=[0,1,1,0,1,0,0,1]
integer::trial,rep,sample,q,batches
integer(int64)::start_clock,end_clock,rate
real(dp)::record(50),values(6),start_cpu,end_cpu,expected
logical::rebuilt
call scalar_allocate(store,groups*span*nk*ns,span,nk,ns,-10,30,.true.,rebuilt)
values=[1.0_dp,2.0_dp,3.0_dp,1.0_dp,2.0_dp,3.0_dp]
batches=store%samples/span
expected=real(store%samples,dp)*(sum(values)-11.0_dp)
call system_clock(count_rate=rate)
do trial=1,size(order)
 call scalar_fill(store,-1.0_dp)
 call cpu_time(start_cpu)
 call system_clock(start_clock)
 if(order(trial)==0)then
  do rep=1,repeats
   do sample=0,store%samples-1
    record=scalar_read_range(store,50*sample+1,50*(sample+1))
    record(15:20)=values
    call scalar_write_range(store,50*sample+1,50*(sample+1),record)
   end do
  end do
 else
  do rep=1,repeats
   do sample=0,store%samples-1
    record=1.0e290_dp
    record(15:20)=values
    call scalar_write_range(store,50*sample+15,50*sample+17,record(15:17))
    call scalar_write_range(store,50*sample+18,50*sample+20,record(18:20))
   end do
  end do
 end if
 call system_clock(end_clock)
 call cpu_time(end_cpu)
 if(abs(sum(store%field)-expected)>0.0_dp)error stop 'benchmark result differs'
 write(*,'(2(i0,1x),2(es24.16,1x))')trial,order(trial),end_cpu-start_cpu, &
      real(end_clock-start_clock,dp)/real(rate,dp)
end do
end program
