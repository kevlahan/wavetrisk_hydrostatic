module parallel_block_rhs_wire_mod
  ! Homogeneous MPI topology metadata only. Fixed numerical records are copied
  ! as bytes; allocatable descriptors are never transmitted. Runtime values do
  ! not use this format. The header rejects incompatible version/kind sizes.
  use iso_fortran_env, only : int8
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  implicit none
  private
  public :: pack_rhs_placement,unpack_rhs_placement
contains
  subroutine pack_rhs_placement(p,original,owned,bytes)
    type(Rhs_Domain_Operators),intent(in) :: p
    integer,intent(in) :: original(:)
    logical,intent(in) :: owned(:)
    integer(int8),allocatable,intent(out) :: bytes(:)
    integer :: pass,pos,n,nb,l,first,last,header(8)
    if(.not.allocated(p%dynamics)) error stop 'RHS wire missing dynamics'
    if(size(original)/=p%nnode.or.size(owned)/=p%nnode) error stop 'RHS wire slot extent'
    first=lbound(p%diffusion,1);last=ubound(p%diffusion,1)
    header=[278001,storage_size(0),storage_size(.false.),storage_size(0.0d0),p%domain_id,p%nnode,first,last]
    do pass=1,2
       pos=1
       nb=storage_size(header(1))/8*size(header)
       if(pass==2) bytes(pos:pos+nb-1)=transfer(header,[0_int8],nb)
       pos=pos+nb
       n=size(p%surface_visits)
       nb=storage_size(n)/8
       if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
       pos=pos+nb
       nb=storage_size(p%surface_visits)/8*size(p%surface_visits)
       if(pass==2) bytes(pos:pos+nb-1)=transfer(p%surface_visits,[0_int8],nb)
       pos=pos+nb
       n=size(original)
       nb=storage_size(n)/8
       if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
       pos=pos+nb
       nb=storage_size(original)/8*size(original)
       if(pass==2) bytes(pos:pos+nb-1)=transfer(original,[0_int8],nb)
       pos=pos+nb
       n=size(owned)
       nb=storage_size(n)/8
       if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
       pos=pos+nb
       nb=storage_size(owned)/8*size(owned)
       if(pass==2) bytes(pos:pos+nb-1)=transfer(owned,[0_int8],nb)
       pos=pos+nb
       n=size(p%gradient)
       nb=storage_size(n)/8
       if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
       pos=pos+nb
       nb=storage_size(p%gradient)/8*size(p%gradient)
       if(pass==2) bytes(pos:pos+nb-1)=transfer(p%gradient,[0_int8],nb)
       pos=pos+nb
       do l=first,last
          n=size(p%pressure(l)%patch)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%pressure(l)%patch)/8*size(p%pressure(l)%patch)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%pressure(l)%patch,[0_int8],nb)
          pos=pos+nb
          n=size(p%pressure(l)%node)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%pressure(l)%node)/8*size(p%pressure(l)%node)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%pressure(l)%node,[0_int8],nb)
          pos=pos+nb
          n=size(p%pressure(l)%producer_ordinal)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%pressure(l)%producer_ordinal)/8*size(p%pressure(l)%producer_ordinal)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%pressure(l)%producer_ordinal,[0_int8],nb)
          pos=pos+nb
          n=size(p%diffusion(l)%producer)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%diffusion(l)%producer)/8*size(p%diffusion(l)%producer)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%diffusion(l)%producer,[0_int8],nb)
          pos=pos+nb
          n=size(p%diffusion(l)%restriction)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%diffusion(l)%restriction)/8*size(p%diffusion(l)%restriction)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%diffusion(l)%restriction,[0_int8],nb)
          pos=pos+nb
          n=size(p%diffusion(l)%divergence)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%diffusion(l)%divergence)/8*size(p%diffusion(l)%divergence)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%diffusion(l)%divergence,[0_int8],nb)
          pos=pos+nb
          n=size(p%rotation(l)%curl)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%rotation(l)%curl)/8*size(p%rotation(l)%curl)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%rotation(l)%curl,[0_int8],nb)
          pos=pos+nb
          n=size(p%rotation(l)%pentagon)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%rotation(l)%pentagon)/8*size(p%rotation(l)%pentagon)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%rotation(l)%pentagon,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%divergence)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%divergence)/8*size(p%dynamics(l)%divergence)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%divergence,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%physics)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%physics)/8*size(p%dynamics(l)%physics)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%physics,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%injection,2)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%injection)/8*size(p%dynamics(l)%injection)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%injection,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%source%action)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%source%action)/8*size(p%dynamics(l)%source%action)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%source%action,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%source%operand)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%source%operand)/8*size(p%dynamics(l)%source%operand)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%source%operand,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%source%direct)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%source%direct)/8*size(p%dynamics(l)%source%direct)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%source%direct,[0_int8],nb)
          pos=pos+nb
          n=size(p%dynamics(l)%source%restriction)
          nb=storage_size(n)/8
          if(pass==2) bytes(pos:pos+nb-1)=transfer(n,[0_int8],nb)
          pos=pos+nb
          nb=storage_size(p%dynamics(l)%source%restriction)/8*size(p%dynamics(l)%source%restriction)
          if(pass==2) bytes(pos:pos+nb-1)=transfer(p%dynamics(l)%source%restriction,[0_int8],nb)
          pos=pos+nb
       end do
       if(pass==1) allocate(bytes(pos-1))
    end do
  end subroutine

  subroutine unpack_rhs_placement(bytes,p,original,owned)
    integer(int8),intent(in) :: bytes(:)
    type(Rhs_Domain_Operators),intent(out) :: p
    integer,allocatable,intent(out) :: original(:)
    logical,allocatable,intent(out) :: owned(:)
    integer :: pos,n,nb,l,first,last,header(8)
    pos=1;nb=storage_size(header(1))/8*size(header)
    call require_bytes(nb)
    header=transfer(bytes(pos:pos+nb-1),header);pos=pos+nb
    if(any(header(:4)/=[278001,storage_size(0),storage_size(.false.),storage_size(0.0d0)])) &
         error stop 'RHS wire incompatible header'
    p%domain_id=header(5);p%nnode=header(6);first=header(7);last=header(8)
    if(p%nnode<0.or.first<0.or.last<first.or.last>100) error stop 'RHS wire invalid extent'
    allocate(p%diffusion(first:last),p%pressure(first:last),p%rotation(first:last),p%dynamics(first:last))
    call read_count(n)
    allocate(p%surface_visits(n))
    nb=storage_size(p%surface_visits)/8*size(p%surface_visits)
    call require_bytes(nb)
    p%surface_visits=transfer(bytes(pos:pos+nb-1),p%surface_visits)
    pos=pos+nb
    call read_count(n)
    allocate(original(n))
    nb=storage_size(original)/8*size(original)
    call require_bytes(nb)
    original=transfer(bytes(pos:pos+nb-1),original)
    pos=pos+nb
    call read_count(n)
    allocate(owned(n))
    nb=storage_size(owned)/8*size(owned)
    call require_bytes(nb)
    owned=transfer(bytes(pos:pos+nb-1),owned)
    pos=pos+nb
    call read_count(n)
    allocate(p%gradient(n))
    nb=storage_size(p%gradient)/8*size(p%gradient)
    call require_bytes(nb)
    p%gradient=transfer(bytes(pos:pos+nb-1),p%gradient)
    pos=pos+nb
    do l=first,last
       p%diffusion(l)%level=l
       call read_count(n)
       allocate(p%pressure(l)%patch(n))
       nb=storage_size(p%pressure(l)%patch)/8*size(p%pressure(l)%patch)
       call require_bytes(nb)
       p%pressure(l)%patch=transfer(bytes(pos:pos+nb-1),p%pressure(l)%patch)
       pos=pos+nb
       call read_count(n)
       allocate(p%pressure(l)%node(n))
       nb=storage_size(p%pressure(l)%node)/8*size(p%pressure(l)%node)
       call require_bytes(nb)
       p%pressure(l)%node=transfer(bytes(pos:pos+nb-1),p%pressure(l)%node)
       pos=pos+nb
       call read_count(n)
       allocate(p%pressure(l)%producer_ordinal(n))
       nb=storage_size(p%pressure(l)%producer_ordinal)/8*size(p%pressure(l)%producer_ordinal)
       call require_bytes(nb)
       p%pressure(l)%producer_ordinal=transfer(bytes(pos:pos+nb-1),p%pressure(l)%producer_ordinal)
       pos=pos+nb
       call read_count(n)
       allocate(p%diffusion(l)%producer(n))
       nb=storage_size(p%diffusion(l)%producer)/8*size(p%diffusion(l)%producer)
       call require_bytes(nb)
       p%diffusion(l)%producer=transfer(bytes(pos:pos+nb-1),p%diffusion(l)%producer)
       pos=pos+nb
       call read_count(n)
       allocate(p%diffusion(l)%restriction(n))
       nb=storage_size(p%diffusion(l)%restriction)/8*size(p%diffusion(l)%restriction)
       call require_bytes(nb)
       p%diffusion(l)%restriction=transfer(bytes(pos:pos+nb-1),p%diffusion(l)%restriction)
       pos=pos+nb
       call read_count(n)
       allocate(p%diffusion(l)%divergence(n))
       nb=storage_size(p%diffusion(l)%divergence)/8*size(p%diffusion(l)%divergence)
       call require_bytes(nb)
       p%diffusion(l)%divergence=transfer(bytes(pos:pos+nb-1),p%diffusion(l)%divergence)
       pos=pos+nb
       call read_count(n)
       allocate(p%rotation(l)%curl(n))
       nb=storage_size(p%rotation(l)%curl)/8*size(p%rotation(l)%curl)
       call require_bytes(nb)
       p%rotation(l)%curl=transfer(bytes(pos:pos+nb-1),p%rotation(l)%curl)
       pos=pos+nb
       call read_count(n)
       allocate(p%rotation(l)%pentagon(n))
       nb=storage_size(p%rotation(l)%pentagon)/8*size(p%rotation(l)%pentagon)
       call require_bytes(nb)
       p%rotation(l)%pentagon=transfer(bytes(pos:pos+nb-1),p%rotation(l)%pentagon)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%divergence(n))
       nb=storage_size(p%dynamics(l)%divergence)/8*size(p%dynamics(l)%divergence)
       call require_bytes(nb)
       p%dynamics(l)%divergence=transfer(bytes(pos:pos+nb-1),p%dynamics(l)%divergence)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%physics(n))
       nb=storage_size(p%dynamics(l)%physics)/8*size(p%dynamics(l)%physics)
       call require_bytes(nb)
       p%dynamics(l)%physics=transfer(bytes(pos:pos+nb-1),p%dynamics(l)%physics)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%injection(2,n))
       nb=storage_size(p%dynamics(l)%injection)/8*size(p%dynamics(l)%injection)
       call require_bytes(nb)
       p%dynamics(l)%injection=reshape(transfer(bytes(pos:pos+nb-1), &
            reshape(p%dynamics(l)%injection,[size(p%dynamics(l)%injection)])),[2,n])
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%source%action(n))
       nb=storage_size(p%dynamics(l)%source%action)/8*size(p%dynamics(l)%source%action)
       call require_bytes(nb)
       p%dynamics(l)%source%action=transfer(bytes(pos:pos+nb-1),p%dynamics(l)%source%action)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%source%operand(n))
       nb=storage_size(p%dynamics(l)%source%operand)/8*size(p%dynamics(l)%source%operand)
       call require_bytes(nb)
       p%dynamics(l)%source%operand=transfer(bytes(pos:pos+nb-1),p%dynamics(l)%source%operand)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%source%direct(n))
       nb=storage_size(p%dynamics(l)%source%direct)/8*size(p%dynamics(l)%source%direct)
       call require_bytes(nb)
       p%dynamics(l)%source%direct=transfer(bytes(pos:pos+nb-1),p%dynamics(l)%source%direct)
       pos=pos+nb
       call read_count(n)
       allocate(p%dynamics(l)%source%restriction(n))
       nb=storage_size(p%dynamics(l)%source%restriction)/8*size(p%dynamics(l)%source%restriction)
       call require_bytes(nb)
       p%dynamics(l)%source%restriction=transfer(bytes(pos:pos+nb-1), &
            p%dynamics(l)%source%restriction)
       pos=pos+nb
    end do
    if(pos/=size(bytes)+1) error stop 'RHS wire trailing bytes'
    if(size(original)/=p%nnode.or.size(owned)/=p%nnode.or.size(p%surface_visits)/=p%nnode) &
         error stop 'RHS wire decoded slot extent'
  contains
    subroutine require_bytes(count)
      integer,intent(in) :: count
      if(count<0.or.count>size(bytes)-pos+1) error stop 'RHS wire truncated packet'
    end subroutine
    subroutine read_count(count)
      integer,intent(out) :: count
      integer :: width
      width=storage_size(count)/8
      call require_bytes(width)
      count=transfer(bytes(pos:pos+width-1),count);pos=pos+width
      if(count<0.or.count>size(bytes)-pos+1) error stop 'RHS wire invalid record count'
    end subroutine
  end subroutine
end module
