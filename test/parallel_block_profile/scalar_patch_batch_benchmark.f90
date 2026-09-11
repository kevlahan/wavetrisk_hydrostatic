program benchmark
use kind_mod, only: dp
use parallel_block_scalar_storage_mod
use iso_fortran_env, only: int64
implicit none
type(Scalar_Record_Storage)::store
integer,parameter::span=16,nk=41,ns=2,groups=512,repeats=4
integer,parameter::order(8)=[0,1,1,0,1,0,0,1]
integer::trial,rep,sample,q,batches
integer(int64)::start_clock,end_clock,rate
real(dp)::values(3,span),start_cpu,end_cpu,expected
logical::rebuilt
call scalar_allocate(store,groups*span*nk*ns,span,nk,ns,-10,30,.true.,rebuilt)
values=reshape([(real(q,dp),q=1,3*span)],shape(values))
batches=store%samples/span
expected=real(batches,dp)*(sum(values)-real(14*span,dp))
call system_clock(count_rate=rate)
do trial=1,size(order)
 call scalar_fill(store,-1.0_dp)
 call cpu_time(start_cpu)
 call system_clock(start_clock)
 if(order(trial)==0)then
  do rep=1,repeats
   do sample=0,store%samples-1,span
    do q=1,span
     call scalar_write_range(store,50*(sample+q-1)+9,50*(sample+q-1)+11,values(:,q))
    end do
   end do
  end do
 else
  do rep=1,repeats
   do sample=0,store%samples-1,span
    call scalar_write_field_records(store,sample,9,values)
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
