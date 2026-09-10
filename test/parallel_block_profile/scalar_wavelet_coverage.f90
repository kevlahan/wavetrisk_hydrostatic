! Traversal regression: producer.inc is extracted from the production source.
! Stubs isolate coverage, masks and transform range, not stencil accuracy.
module coverage_test
  use iso_fortran_env, only: int64,real64
  use ieee_arithmetic, only: ieee_is_finite
  implicit none
  integer, parameter :: dp=real64,PATCH_SIZE=4,ADJZONE=8,zlevels=2
  integer, parameter :: BLOCK_SCALAR_WAVELET_MASK_INDEX=1
  type :: patch_type
     integer :: level=5,children(4)=0
  end type
  type :: Block_Data
     type(patch_type) :: patch(2)
     integer :: n_scalar_variable=2,n_field_level=2,field_level=1
     real(dp) :: values(0:3,0:3,2,0:1,2)=-1
  end type
  type :: readiness
     logical :: ready=.true.
     integer :: catalog_index=1
  end type
  type(readiness) :: block_scalar_tendency(1)
  include 'context.inc'
contains
  integer function catalog_local_block(i)
    integer,intent(in)::i
    catalog_local_block=i
  end function
  subroutine fail(message)
    character(*),intent(in)::message
    print *,message
    error stop 1
  end subroutine
  real(dp) function block_scalar_record_value(b,l,p,s,k,i,j,field)
    type(Block_Data),intent(in)::b
    integer,intent(in)::l,p,s,k,i,j,field
    block_scalar_record_value=real(ADJZONE,dp)
    if(i==1.and.j==0) block_scalar_record_value=0
  end function
  real(dp) function compute_block_scalar_wavelet(b,l,p,s,k,i,j,si,sj,scale) result(value)
    type(Block_Data),intent(in)::b
    integer,intent(in)::l,p,s,k,i,j,si(4),sj(4)
    real(dp),intent(out)::scale
    value=real(100*p+10*k+2*s+i+j,dp)
    scale=value
  end function
  subroutine set_block_scalar_wavelet_value(b,p,s,k,i,j,value)
    type(Block_Data),intent(inout)::b
    integer,intent(in)::p,s,k,i,j
    real(dp),intent(in)::value
    b%values(i,j,k,s,p)=value
  end subroutine
  include 'producer.inc'
end module

program test
  use coverage_test
  implicit none
  type(Block_Data)::b
  type(Block_Scalar_Wavelet_Context)::stats
  integer::i,j,k,s,p
  b%patch(2)%level=6
  b%patch(1)%children(1)=1
  stats%first_level=4
  stats%last_level=5
  call produce_block_scalar_wavelets(1,b,stats)
  if(stats%target_patch_count/=2.or.stats%candidate_count/=96) error stop 'missing root coverage'
  if(stats%produced_count/=96.or.stats%active_count/=88) error stop 'site/mask coverage'
  do p=1,2
     do s=0,1
        do k=1,2
           do j=0,3
              do i=0,3
                 if(mod(i,2)==0.and.mod(j,2)==0) then
                    if(b%values(i,j,k,s,p)/=-1) error stop 'scaling site overwritten'
                 else if(i==1.and.j==0) then
                    if(b%values(i,j,k,s,p)/=0) error stop 'inactive site not zero'
                 else
                    if(b%values(i,j,k,s,p)/=real(100*p+10*k+2*s+i+j,dp)) &
                         error stop 'missing/wrong target'
                 end if
              end do
           end do
        end do
     end do
  end do
  stats=Block_Scalar_Wavelet_Context()
  stats%first_level=5
  stats%last_level=5
  b%values=-1
  call produce_block_scalar_wavelets(1,b,stats)
  if(stats%target_patch_count/=1.or.stats%candidate_count/=48) error stop 'provisional coverage'
  if(any(b%values(:,:,:,:,1)/=-1)) error stop 'provisional root modified'
  ! An unrefined root must not depend on a child or a retained parent.
  stats=Block_Scalar_Wavelet_Context()
  stats%first_level=4
  stats%last_level=4
  b%patch(1)%children=0
  b%values=-1
  call produce_block_scalar_wavelets(1,b,stats)
  if(stats%target_patch_count/=1.or.stats%candidate_count/=48) error stop 'leaf root omitted'
  if(any(b%values(:,:,:,:,2)/=-1)) error stop 'out-of-range patch modified'
  print *, 'scalar wavelet root/provisional/mask coverage passed'
end program
