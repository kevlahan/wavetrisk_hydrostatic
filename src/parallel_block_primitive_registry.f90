module parallel_block_primitive_registry_mod
  ! Geometry-side identities only. No numerical primitive is stored here.
  use iso_fortran_env, only: int64
  use domain_mod, only: grid
  use arch_mod, only: rank,abort_run
  implicit none
  private
  public :: declare_primitive_result,registry_needs_declarations,finish_primitive_declarations
  public :: reset_primitive_registry,register_primitive_result,primitive_result_source
  type :: Primitive_Definitions
     integer,allocatable :: pv(:,:,:),vort(:,:,:)
     logical,allocatable :: pv_written(:,:),vort_written(:,:)
  end type
  type(Primitive_Definitions),allocatable :: definition(:)
  integer(int64) :: layout=-1_int64
  integer :: layer=0
  logical :: registry_needs_declarations=.true.
contains
  subroutine reset_primitive_registry(k,generation)
    integer,intent(in)::k
    integer(int64),intent(in)::generation
    integer::d,n,i,e
    if(layout/=generation) then
       if(allocated(definition)) deallocate(definition)
       allocate(definition(size(grid)))
       do d=1,size(grid)
          n=grid(d)%node%length
          allocate(definition(d)%pv(7,3,n),definition(d)%vort(7,2,n), &
               definition(d)%pv_written(3,n),definition(d)%vort_written(2,n))
          definition(d)%pv_written=.false.;definition(d)%vort_written=.false.
       end do
       layout=generation;registry_needs_declarations=.true.
    end if
    do d=1,size(definition)
       definition(d)%pv=0;definition(d)%vort=0
       do i=1,size(definition(d)%pv,3)
          do e=1,3
             if(.not.definition(d)%pv_written(e,i)) definition(d)%pv(:,e,i)=[4,0,0,0,1,0,0]
          end do
          do e=1,2
             if(.not.definition(d)%vort_written(e,i)) definition(d)%vort(:,e,i)=[4,0,0,0,1,0,0]
          end do
       end do
    end do
    layer=k
  end subroutine
  subroutine declare_primitive_result(d,id,e,vorticity)
    integer,intent(in)::d,id,e
    logical,intent(in)::vorticity
    if(vorticity) then
       definition(d)%vort_written(e,id+1)=.true.;definition(d)%vort(:,e,id+1)=0
    else
       definition(d)%pv_written(e,id+1)=.true.;definition(d)%pv(:,e,id+1)=0
    end if
  end subroutine
  subroutine finish_primitive_declarations
    registry_needs_declarations=.false.
  end subroutine
  subroutine register_primitive_result(d,id,component,vorticity,key,k)
    integer,intent(in)::d,id,component,key(7),k
    logical,intent(in)::vorticity
    if(layer/=k) call die('registration layer')
    if(vorticity) then
       definition(d)%vort(:,component,id+1)=key
    else
       definition(d)%pv(:,component,id+1)=key
    end if
  end subroutine
  function primitive_result_source(d,id,component,vorticity,k,required) result(key)
    integer,intent(in)::d,id,component,k
    logical,intent(in)::vorticity
    integer::key(7)
    logical,intent(in),optional::required
    logical::must_exist
    if(layer/=k) call die('lookup layer')
    if(vorticity) then
       key=definition(d)%vort(:,component,id+1)
    else
       key=definition(d)%pv(:,component,id+1)
    end if
    must_exist=.true.
    if(present(required)) must_exist=required
    if(key(1)==0.and.must_exist) call die('missing current primitive producer')
  end function
  subroutine die(message)
    character(*),intent(in)::message
    write(*,'(a,i0,2a)') 'Rank ',rank,': primitive registry: ',message
    call abort_run
    error stop 'primitive registry abort returned'
  end subroutine
end module
