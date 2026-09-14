module fixture
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use ieee_arithmetic
  use parallel_block_scalar_storage_mod
  implicit none
  include 'constants.inc'
  integer, parameter :: PATCH_SIZE=4,nk=41,ns=2,stride=64
  type :: Block_Data
     integer :: n_field_level=nk
     integer :: bdry_node(stride),ghost_node(stride)
  end type
  type :: tendency
     type(Scalar_Record_Storage) :: patch,bdry,ghost
  end type
  type :: plan
     integer :: scalar_patch_nvalue=nk*ns*PATCH_SIZE**2
  end type
  type(tendency) :: block_scalar_tendency(1)
  type(plan) :: block_writeback_plan
  integer(int64) :: inverse_address_work(8)=0,scalar_reads=0
  logical :: block_profile=.true.
contains
  subroutine fail(message)
    character(*),intent(in) :: message
    print *,message
    error stop 1
  end subroutine
  subroutine locate_block_scalar_record(block,local_index,p,x,y,storage_class,record,node)
    type(Block_Data),intent(in) :: block
    integer,intent(in) :: local_index,p,x,y
    integer,intent(out) :: storage_class,record,node
    storage_class=modulo(x+2*y,3)+1
    record=modulo(p,2)
    node=modulo(7*x+11*y,merge(PATCH_SIZE**2,stride,storage_class==STORE_PATCH))
  end subroutine
  function block_scalar_record_value(block,local_index,p,scalar_slot,level_slot,x,y,slot) result(value)
    type(Block_Data),intent(in) :: block
    integer,intent(in) :: local_index,p,scalar_slot,level_slot,x,y,slot
    real(dp) :: value
    integer :: storage_class,record,node,base
    call locate_block_scalar_record(block,local_index,p,x,y,storage_class,record,node)
    base=scalar_slot*nk+level_slot-1
    select case(storage_class)
    case(STORE_PATCH)
       base=record*block_writeback_plan%scalar_patch_nvalue+base*PATCH_SIZE**2+node
       value=scalar_read(block_scalar_tendency(1)%patch,50*base+slot)
    case(STORE_BDRY)
       value=scalar_read(block_scalar_tendency(1)%bdry,50*(base*stride+node)+slot)
    case(STORE_GHOST)
       value=scalar_read(block_scalar_tendency(1)%ghost,50*(base*stride+node)+slot)
    end select
    scalar_reads=scalar_reads+1
  end function
  subroutine seed(store,nodes,groups,compact,tag)
    type(Scalar_Record_Storage),intent(inout) :: store
    integer,intent(in) :: nodes,groups,tag
    logical,intent(in) :: compact
    logical :: rebuilt
    integer :: sample,slot,layer,node,field
    real(dp) :: values(50)
    call scalar_allocate(store,nodes*nk*ns*groups,nodes,nk,ns,-10,30,compact,rebuilt)
    call scalar_fill(store,0.0_dp)
    do sample=0,store%samples-1
       field=modulo(sample/nodes,nk*ns)
       layer=modulo(field,nk)
       node=modulo(sample,nodes)+nodes*(sample/(nodes*nk*ns))
       do slot=1,50
          ! Distinct signs and magnitudes, non-binary fractions. Compact
          ! geometry is shared; full mode deliberately varies it by layer.
          values(slot)=(-1.0_dp)**modulo(slot+node,2)* &
               real(1+slot*7+node*11+tag*13,dp)/real(31+slot+tag,dp)
          if (any(slot==[1,2,3,4,5,6,9,10,11,15,16,17,18,19,20,34,35]).or..not.compact) &
               values(slot)=values(slot)*real(1+field,dp)/7.0_dp
       end do
       call scalar_write_range(store,50*sample+1,50*(sample+1),values)
    end do
  end subroutine
  include 'indices.inc'
  include 'interpolate_block_scalar_flux.inc'
  include 'block_restriction_edge_flux.inc'
  include 'interpolate_block_scalar_flux_columns.inc'
end module
program check
  use fixture
  implicit none
  type(Block_Data) :: block
  real(dp) :: columns(30,4),reference(4),poison
  integer :: mode,scalar,k,i,j,first,count,trial,rep
  integer, parameter :: order(8)=[0,1,1,0,1,0,0,1]
  real(dp) :: started,elapsed,checksum,expected
  integer(int64) :: t0,t1,rate
  character(16) :: argument
  call get_command_argument(1,argument)
  do mode=0,1
     call seed(block_scalar_tendency(1)%patch,PATCH_SIZE**2,2,mode==0,1)
     call seed(block_scalar_tendency(1)%bdry,stride,1,mode==0,2)
     call seed(block_scalar_tendency(1)%ghost,stride,1,mode==0,3)
     if(argument=='benchmark')then
        first=12
        count=30
        expected=0.0_dp
        call system_clock(count_rate=rate)
        do trial=1,size(order)
           checksum=0.0_dp
           call cpu_time(started)
           call system_clock(t0)
           do rep=1,30000
              i=2*modulo(rep,2)
              j=2*modulo(rep/2,2)
              scalar=modulo(rep/4,2)
              if(order(trial)==0)then
                 do k=1,count
                    columns(k,:)=interpolate_block_scalar_flux( &
                         block,1,2,scalar,first+k-1,i,j,BLOCK_SCALAR_DIRECT_FLUX_START)
                 end do
              else
                 columns=interpolate_block_scalar_flux_columns( &
                      block,1,2,scalar,first,count,i,j,BLOCK_SCALAR_DIRECT_FLUX_START)
              end if
              checksum=checksum+sum(columns)
           end do
           call system_clock(t1)
           call cpu_time(elapsed)
           if(trial==1)expected=checksum
           if(transfer(checksum,0_int64)/=transfer(expected,0_int64))error stop 'benchmark checksum differs'
           write(*,'(2(i0,1x),3(es24.16,1x))')trial,order(trial),elapsed-started,real(t1-t0,dp)/real(rate,dp),checksum
        end do
        exit
     end if
     if(argument=='poison'.or.argument=='nan')then
        poison=0.125_dp*huge(0.0_dp)
        if(argument=='nan')poison=ieee_value(0.0_dp,ieee_quiet_nan)
        call scalar_fill(block_scalar_tendency(1)%patch,poison)
        call scalar_fill(block_scalar_tendency(1)%bdry,poison)
        call scalar_fill(block_scalar_tendency(1)%ghost,poison)
     end if
     do scalar=0,1
        do first=12,19,7
           count=nk-first+1
           if(argument=='extent')count=count+1
           do j=0,2,2
              do i=0,2,2
                 inverse_address_work=0
                 scalar_reads=0
                 columns(:count,:)=interpolate_block_scalar_flux_columns( &
                      block,1,2,scalar,first,count,i,j,BLOCK_SCALAR_DIRECT_FLUX_START)
                 do k=1,count
                    reference=interpolate_block_scalar_flux( &
                         block,1,2,scalar,first+k-1,i,j,BLOCK_SCALAR_DIRECT_FLUX_START)
                    if(any(transfer(reference,[0_int64],4)/=transfer(columns(k,:),[0_int64],4))) &
                         error stop 'interpolation bits differ'
                 end do
                 if(scalar_reads/=inverse_address_work(8))error stop 'work count differs'
                 if(scalar_reads/=count*inverse_address_work(7))error stop 'column scaling differs'
              end do
           end do
        end do
     end do
  end do
  print *,'PASS restriction columns'
end program
