program check
  use fixture
  implicit none
  type(Block_Data) :: block
  type(Block_Scalar_Restriction_Cursor) :: cursor
  real(dp) :: sm(30,4),values(30),reference,poison
  integer :: scalar,first,count,i,j,k,e,trial,rep
  integer, parameter :: order(8)=[0,1,1,0,1,0,0,1]
  integer(int64) :: t0,t1,rate
  real(dp) :: checksum,expected,started,elapsed
  character(16) :: argument
  call get_command_argument(1,argument)
  call seed(block_scalar_tendency(1)%patch,PATCH_SIZE**2,2,argument/='full',1)
  call seed(block_scalar_tendency(1)%bdry,stride,1,argument/='full',2)
  call seed(block_scalar_tendency(1)%ghost,stride,1,argument/='full',3)
  if(argument=='poison'.or.argument=='nan')then
     poison=0.125_dp*huge(0.0_dp)
     if(argument=='nan')poison=ieee_value(0.0_dp,ieee_quiet_nan)
     call scalar_fill(block_scalar_tendency(1)%patch,poison)
     call scalar_fill(block_scalar_tendency(1)%bdry,poison)
     call scalar_fill(block_scalar_tendency(1)%ghost,poison)
  end if
  if(argument=='benchmark')then
     first=12; count=30
     sm=0.0_dp
     block_profile=.false.
     call system_clock(count_rate=rate)
     do trial=1,size(order)
        checksum=0.0_dp
        call cpu_time(started)
        call system_clock(t0)
        do rep=1,6000
           i=2*modulo(rep,2);j=2*modulo(rep/2,2);scalar=modulo(rep/4,2);e=modulo(rep,3)
           if(order(trial)==0)then
              do k=1,count
                 values(k)=complete_block_coarse_flux(sm(k,:),block,1,2,scalar,first+k-1,0,0,i,j,e,cursor)
              end do
           else
              values=complete_block_coarse_flux_columns(sm,block,1,2,scalar,first,count,0,0,i,j,e,cursor)
           end if
           checksum=checksum+sum(values)
        end do
        call system_clock(t1)
        call cpu_time(elapsed)
        if(trial==1)expected=checksum
        if(transfer(checksum,0_int64)/=transfer(expected,0_int64))error stop 'benchmark checksum differs'
        write(*,'(2(i0,1x),3(es24.16,1x))')trial,order(trial),elapsed-started,real(t1-t0,dp)/real(rate,dp),checksum
     end do
     stop
  end if
  do scalar=0,1
     do first=12,29,17
        count=nk-first+1
        do j=0,2,2
           do i=0,2,2
              sm=0.0_dp
              if(argument=='valid')sm(:count,:)=interpolate_block_scalar_flux_columns( &
                   block,1,2,scalar,first,count,i,j,BLOCK_SCALAR_DIRECT_FLUX_START)
              do e=RT,UP
                 inverse_address_work=0
                 values(:count)=complete_block_coarse_flux_columns( &
                      sm(:count,:),block,1,2,scalar,merge(0,first,argument=='extent'),count,0,0,i,j,e,cursor)
                 do k=1,count
                    reference=complete_block_coarse_flux(sm(k,:),block,1,2,scalar,first+k-1, &
                         0,0,i,j,e,cursor)
                    if(transfer(reference,0_int64)/=transfer(values(k),0_int64))then
                       print *,scalar,first,i,j,e,k,reference,values(k)
                       error stop 'closure bits differ'
                    end if
                 end do
                 if(inverse_address_work(9)/=1)error stop 'missing closure counter'
                 if(inverse_address_work(11)==0)error stop 'no reused input columns'
                 if(inverse_address_work(12)==0)error stop 'missing scalar fallback coverage'
              end do
           end do
        end do
     end do
  end do
  print *,'PASS restriction closure'
end program
