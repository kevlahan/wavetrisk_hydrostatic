module parallel_block_scalar_storage_mod
  ! Scalar execution records: field values vary with variable/vertical level;
  ! horizontal geometry does not. Keep the oracle's full records independent.
  use kind_mod, only : dp
  use iso_fortran_env, only : int64
  implicit none
  private
  public :: Scalar_Record_Storage, scalar_allocate, scalar_release, scalar_extent, scalar_capacity
  public :: scalar_geometry_capacity, scalar_read, scalar_read_range, scalar_write, scalar_write_range, scalar_fill
  public :: scalar_is_allocated
  public :: scalar_seed_patch
  public :: scalar_install_geometry, scalar_share_inactive, scalar_fill_fields
  public :: SCALAR_RECORD_WIDTH, SCALAR_FIELD_SLOTS, SCALAR_SHARED_SLOTS

  integer, parameter :: SCALAR_RECORD_WIDTH=50, SCALAR_FIELD_SLOTS=17, SCALAR_SHARED_SLOTS=33
  integer, parameter :: slot_map(50)=[1,2,3,4,5,6,-1,-2,7,8,9,-3,-4,-5,10,11,12,13,14,15, &
       -6,-7,-8,-9,-10,-11,-12,-13,-14,-15,-16,-17,-18,16,17, &
       -19,-20,-21,-22,-23,-24,-25,-26,-27,-28,-29,-30,-31,-32,-33]

  type :: Scalar_Record_Storage
     integer :: samples=0, stride=0, nlevel=0, nfield=0, first_level=0, last_physical=0
     integer :: physical_field=0, inactive_field=-1
     logical :: compact=.false.
     real(dp), allocatable :: full(:)
     real(dp), allocatable :: field(:,:)
     ! A separate inactive class preserves local soil/surface zero records
     ! versus transported records carrying geometry on every field level.
     real(dp), allocatable :: geometry(:,:,:)
  end type Scalar_Record_Storage

  interface scalar_read
     module procedure read_value,read_indices
  end interface
  interface scalar_write
     module procedure write_value,write_indices
  end interface
  interface scalar_write_range
     module procedure write_range_values,fill_range
  end interface

contains

  subroutine scalar_release(store)
    type(Scalar_Record_Storage), intent(inout) :: store
    if (allocated(store%full)) deallocate(store%full)
    if (allocated(store%field)) deallocate(store%field)
    if (allocated(store%geometry)) deallocate(store%geometry)
    store%samples=0
    store%nfield=0
  end subroutine scalar_release

  subroutine scalar_allocate(store,samples,stride,nlevel,nscalar,first_level,last_physical,compact,rebuilt)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: samples,stride,nlevel,nscalar,first_level,last_physical
    logical, intent(in) :: compact
    logical, intent(out) :: rebuilt
    integer :: f
    if (samples<0.or.stride<0.or.nlevel<1.or.nscalar<1) error stop 'scalar storage dimensions invalid'
    if (samples>0) then
       if (stride==0) error stop 'scalar storage stride invalid'
       if (mod(samples,stride*nlevel*nscalar)/=0) error stop 'scalar storage shape invalid'
    end if
    if (1-first_level<0.or.1-first_level>=nlevel) error stop 'scalar storage physical field absent'
    rebuilt=.false.
    if (allocated(store%full).or.allocated(store%field)) then
       if (store%samples==samples.and.store%stride==stride.and.store%nlevel==nlevel.and. &
            store%nfield==nlevel*nscalar.and.store%first_level==first_level.and. &
            store%last_physical==last_physical.and.(store%compact.eqv.compact)) return
    end if
    call scalar_release(store)
    store%samples=samples
    store%stride=stride
    store%nlevel=nlevel
    store%nfield=nlevel*nscalar
    store%first_level=first_level
    store%last_physical=last_physical
    store%physical_field=1-first_level
    store%inactive_field=-1
    do f=0,nlevel-1
       if (f+first_level>=1.and.f+first_level<=last_physical) cycle
       store%inactive_field=f
       exit
    end do
    store%compact=compact
    if (compact) then
       allocate(store%field(SCALAR_FIELD_SLOTS,samples))
       allocate(store%geometry(SCALAR_SHARED_SLOTS,2,samples/store%nfield))
    else
       allocate(store%full(SCALAR_RECORD_WIDTH*samples))
    end if
    rebuilt=.true.
  end subroutine scalar_allocate

  integer function scalar_extent(store) result(n)
    type(Scalar_Record_Storage), intent(in) :: store
    n=SCALAR_RECORD_WIDTH*store%samples
  end function scalar_extent

  logical function scalar_is_allocated(store) result(ready)
    type(Scalar_Record_Storage), intent(in) :: store
    ready=allocated(store%full).or.allocated(store%field)
  end function scalar_is_allocated

  integer(int64) function scalar_capacity(store) result(n)
    type(Scalar_Record_Storage), intent(in) :: store
    n=0_int64
    if (allocated(store%full)) n=n+size(store%full,kind=int64)
    if (allocated(store%field)) n=n+size(store%field,kind=int64)
    if (allocated(store%geometry)) n=n+size(store%geometry,kind=int64)
    n=n*int(storage_size(0.0_dp)/8,int64)
  end function scalar_capacity

  integer(int64) function scalar_geometry_capacity(store) result(n)
    type(Scalar_Record_Storage), intent(in) :: store
    n=0_int64
    if (allocated(store%full)) n=int(store%samples,int64)*SCALAR_SHARED_SLOTS
    if (allocated(store%geometry)) n=size(store%geometry,kind=int64)
    n=n*int(storage_size(0.0_dp)/8,int64)
  end function scalar_geometry_capacity

  subroutine shared_address(store,sample,node,zone,field)
    type(Scalar_Record_Storage), intent(in) :: store
    integer, intent(in) :: sample
    integer, intent(out) :: node,zone,field
    integer :: level
    field=mod(sample/store%stride,store%nfield)
    node=(sample/(store%stride*store%nfield))*store%stride+mod(sample,store%stride)+1
    level=mod(field,store%nlevel)+store%first_level
    zone=2
    if (level>=1.and.level<=store%last_physical) zone=1
  end subroutine shared_address

  real(dp) function read_value(store,address) result(value)
    type(Scalar_Record_Storage), intent(in) :: store
    integer, intent(in) :: address
    integer :: sample,slot,node,zone,field
    if (address<1.or.address>SCALAR_RECORD_WIDTH*store%samples) error stop 'scalar storage read index invalid'
    if (.not.store%compact) then
       value=store%full(address)
       return
    end if
    sample=(address-1)/SCALAR_RECORD_WIDTH
    slot=slot_map(mod(address-1,SCALAR_RECORD_WIDTH)+1)
    if (slot>0) then
       value=store%field(slot,sample+1)
    else
       call shared_address(store,sample,node,zone,field)
       value=store%geometry(-slot,zone,node)
    end if
  end function read_value

  function read_indices(store,address) result(value)
    type(Scalar_Record_Storage), intent(in) :: store
    integer, intent(in) :: address(:)
    real(dp) :: value(size(address))
    integer :: i,sample,node,zone,field,mapped(size(address))
    if (store%compact.and.size(address)>0) then
       if (minval(address)<1.or.maxval(address)>SCALAR_RECORD_WIDTH*store%samples) &
            error stop 'scalar storage indexed read invalid'
       sample=(address(1)-1)/SCALAR_RECORD_WIDTH
       if (all((address-1)/SCALAR_RECORD_WIDTH==sample)) then
          mapped=slot_map(mod(address-1,SCALAR_RECORD_WIDTH)+1)
          if (all(mapped>0)) then
             value=store%field(mapped,sample+1)
             return
          else if (all(mapped<0)) then
             call shared_address(store,sample,node,zone,field)
             value=store%geometry(-mapped,zone,node)
             return
          end if
       end if
    end if
    do i=1,size(address)
       value(i)=read_value(store,address(i))
    end do
  end function read_indices

  function scalar_read_range(store,first,last) result(value)
    type(Scalar_Record_Storage), intent(in) :: store
    integer, intent(in) :: first,last
    real(dp) :: value(max(0,last-first+1))
    integer :: i,sample,a,b
    if (first<1.or.last>SCALAR_RECORD_WIDTH*store%samples) error stop 'scalar storage read range invalid'
    if (.not.store%compact) then
       value=store%full(first:last)
    else
       sample=(first-1)/SCALAR_RECORD_WIDTH
       if (size(value)>0.and.(last-1)/SCALAR_RECORD_WIDTH==sample) then
          a=slot_map(mod(first-1,SCALAR_RECORD_WIDTH)+1)
          b=slot_map(mod(last-1,SCALAR_RECORD_WIDTH)+1)
          if (a>0.and.b-a==last-first) then
             value=store%field(a:b,sample+1)
             return
          end if
       end if
       do i=1,size(value)
          value(i)=read_value(store,first+i-1)
       end do
    end if
  end function scalar_read_range

  subroutine write_value(store,address,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: address
    real(dp), intent(in) :: value
    integer :: sample,slot,node,zone,field
    if (address<1.or.address>SCALAR_RECORD_WIDTH*store%samples) error stop 'scalar storage write index invalid'
    if (.not.store%compact) then
       store%full(address)=value
       return
    end if
    sample=(address-1)/SCALAR_RECORD_WIDTH
    slot=slot_map(mod(address-1,SCALAR_RECORD_WIDTH)+1)
    if (slot>0) then
       store%field(slot,sample+1)=value
    else
       call shared_address(store,sample,node,zone,field)
       ! Geometry installation is performed by the canonical field once.
       ! Later zeroing of another record must not erase installed geometry.
       if (zone==1.and.field/=store%physical_field) return
       if (zone==2.and.field/=store%inactive_field) return
       store%geometry(-slot,zone,node)=value
    end if
  end subroutine write_value

  subroutine write_indices(store,address,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: address(:)
    real(dp), intent(in) :: value(:)
    integer :: i,sample,mapped(size(address))
    if (size(address)/=size(value)) error stop 'scalar storage indexed write shape invalid'
    if (store%compact.and.size(address)>0) then
       if (minval(address)<1.or.maxval(address)>SCALAR_RECORD_WIDTH*store%samples) &
            error stop 'scalar storage indexed write invalid'
       sample=(address(1)-1)/SCALAR_RECORD_WIDTH
       if (all((address-1)/SCALAR_RECORD_WIDTH==sample)) then
          mapped=slot_map(mod(address-1,SCALAR_RECORD_WIDTH)+1)
          if (all(mapped>0)) then
             store%field(mapped,sample+1)=value
             return
          end if
       end if
    end if
    do i=1,size(address)
       call write_value(store,address(i),value(i))
    end do
  end subroutine write_indices

  subroutine write_range_values(store,first,last,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: first,last
    real(dp), intent(in) :: value(:)
    integer :: i,sample,a,b
    if (first<1.or.last>SCALAR_RECORD_WIDTH*store%samples.or.size(value)/=last-first+1) &
         error stop 'scalar storage write range invalid'
    if (.not.store%compact) then
       store%full(first:last)=value
    else
       sample=(first-1)/SCALAR_RECORD_WIDTH
       if (size(value)>0.and.(last-1)/SCALAR_RECORD_WIDTH==sample) then
          a=slot_map(mod(first-1,SCALAR_RECORD_WIDTH)+1)
          b=slot_map(mod(last-1,SCALAR_RECORD_WIDTH)+1)
          if (a>0.and.b-a==last-first) then
             store%field(a:b,sample+1)=value
             return
          end if
       end if
       do i=1,size(value)
          call write_value(store,first+i-1,value(i))
       end do
    end if
  end subroutine write_range_values

  subroutine fill_range(store,first,last,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: first,last
    real(dp), intent(in) :: value
    integer :: i
    if (first<1.or.last>SCALAR_RECORD_WIDTH*store%samples) error stop 'scalar storage fill range invalid'
    if (.not.store%compact) then
       store%full(first:last)=value
    else
       do i=first,last
          call write_value(store,i,value)
       end do
    end if
  end subroutine fill_range

  subroutine scalar_fill(store,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    real(dp), intent(in) :: value
    if (allocated(store%full)) store%full=value
    if (allocated(store%field)) store%field=value
    if (allocated(store%geometry)) store%geometry=value
  end subroutine scalar_fill

  subroutine scalar_seed_patch(store,first,geometry)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: first
    real(dp), intent(in) :: geometry(:,:)
    integer :: node,last
    if (.not.store%compact.or.size(geometry,1)/=SCALAR_SHARED_SLOTS.or.size(geometry,2)/=store%stride) &
         error stop 'scalar geometry seed shape invalid'
    if (mod(first-1,store%stride*store%nfield)/=0) error stop 'scalar geometry seed alignment invalid'
    last=first+store%stride*store%nfield-1
    if (first<1.or.last>store%samples) error stop 'scalar geometry seed extent invalid'
    node=(first-1)/store%nfield+1
    store%field(:,first:last)=0.0_dp
    store%geometry(:,1,node:node+store%stride-1)=geometry
    store%geometry(:,2,node:node+store%stride-1)=0.0_dp
  end subroutine scalar_seed_patch

  subroutine scalar_install_geometry(store,first,geometry)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: first
    real(dp), intent(in) :: geometry(:,:)
    integer :: last
    if (.not.store%compact.or.size(geometry,1)/=SCALAR_SHARED_SLOTS) error stop 'scalar shared install shape invalid'
    last=first+size(geometry,2)-1
    if (first<1.or.last>size(store%geometry,3)) error stop 'scalar shared install extent invalid'
    store%geometry(:,1,first:last)=geometry
    store%geometry(:,2,first:last)=geometry
  end subroutine scalar_install_geometry

  subroutine scalar_share_inactive(store)
    type(Scalar_Record_Storage), intent(inout) :: store
    if (.not.store%compact) error stop 'scalar inactive sharing requires compact storage'
    store%geometry(:,2,:)=store%geometry(:,1,:)
  end subroutine scalar_share_inactive

  subroutine scalar_fill_fields(store,first,last,value)
    type(Scalar_Record_Storage), intent(inout) :: store
    integer, intent(in) :: first,last
    real(dp), intent(in) :: value
    if (.not.store%compact.or.first<1.or.last>store%samples) error stop 'scalar field fill extent invalid'
    store%field(:,first:last)=value
  end subroutine scalar_fill_fields

end module parallel_block_scalar_storage_mod
