module parallel_block_rhs_placement_mod
  ! Topology-time selection and address compaction. Operation selection follows
  ! the target's numerical owner, while shared writes retain original ordinals.
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  use parallel_block_velocity_mod, only : VELOCITY_DIRECT,VELOCITY_RESTRICT
  implicit none
  private
  public :: select_rhs_placement,compile_rhs_owner_map
contains
  subroutine compile_rhs_owner_map(dom,global_domain,numerical_owner,state_key)
    use domain_mod, only : Domain
    use arch_mod, only : block_catalog
    use shared_mod, only : N_CHDRN
    use patch_mod, only : PATCH_SIZE
    use parallel_block_mod, only : STORE_PATCH
    type(Domain),intent(in) :: dom
    integer,intent(in) :: global_domain
    integer,allocatable,intent(out) :: numerical_owner(:),state_key(:,:)
    integer :: b,first,next_patch
    allocate(numerical_owner(dom%node%length),state_key(4,dom%node%length))
    numerical_owner=-1;state_key=0;first=0
    do b=1,size(block_catalog)
       if(block_catalog(b)%root_domain/=global_domain) cycle
       if(first==0) first=b
       next_patch=0
       call visit(block_catalog(b)%root_patch)
    end do
    if(first==0) error stop 'RHS owner map missing Domain catalog'
    ! Boundary/scaffold numerical work has a deterministic final block owner.
    ! Its state source is bound separately; it is never assumed to be a Domain
    ! numerical producer merely because it lacks an interior catalogue key.
    where(numerical_owner<0) numerical_owner=block_catalog(first)%owner
  contains
    recursive subroutine visit(p)
      integer,intent(in) :: p
      integer :: a,c,q
      if(p<0.or.p>=dom%patch%length) error stop 'RHS owner map patch outside Domain'
      if(dom%patch%elts(p+1)%deleted) return
      next_patch=next_patch+1;a=dom%patch%elts(p+1)%elts_start
      do q=1,PATCH_SIZE**2
         if(state_key(1,a+q)/=0) error stop 'RHS owner map overlapping block interiors'
         state_key(:,a+q)=[b,STORE_PATCH,next_patch,q-1]
         numerical_owner(a+q)=block_catalog(b)%owner
      end do
      do c=1,N_CHDRN
         q=dom%patch%elts(p+1)%children(c)
         if(q>0) call visit(q)
      end do
    end subroutine
  end subroutine

  subroutine select_rhs_placement(input,owner,rank,output,original,owned)
    type(Rhs_Domain_Operators), intent(in) :: input
    integer, intent(in) :: owner(:),rank
    type(Rhs_Domain_Operators), intent(out) :: output
    integer, allocatable, intent(out) :: original(:)
    logical, allocatable, intent(out) :: owned(:)
    integer, allocatable :: index(:),map(:),dmap(:),rmap(:),actions(:),operands(:)
    logical, allocatable :: seen(:),keep(:)
    integer :: l,s,j,a,b,n,c,d,i,pass,nold,np,nn,first,last
    if(size(owner)/=input%nnode.or.any(owner<0).or.rank<0) error stop 'RHS placement owner extent'
    if(.not.allocated(input%dynamics)) error stop 'RHS placement incomplete dynamics'
    nold=input%nnode;output%domain_id=input%domain_id;output%nnode=nold
    first=lbound(input%diffusion,1);last=ubound(input%diffusion,1)
    allocate(output%diffusion(first:last),output%rotation(first:last), &
         output%pressure(first:last),output%dynamics(first:last))
    output%surface_visits=input%surface_visits
    where(owner/=rank) output%surface_visits=0
    allocate(keep(size(input%gradient)))
    do s=1,size(keep)
       keep(s)=owner(input%gradient(s)%node(1))==rank
    end do
    output%gradient=pack(input%gradient,keep)
    do l=lbound(input%diffusion,1),ubound(input%diffusion,1)
       associate(src=>input%diffusion(l),dst=>output%diffusion(l), &
            sp=>input%pressure(l),dp=>output%pressure(l))
       dst%level=src%level;allocate(dp%patch(size(sp%patch)))
       do pass=1,2
          nn=0;np=0
          do j=1,size(sp%patch)
             a=sp%patch(j)%first_pressure;b=sp%patch(j)%last_pressure
             c=sp%patch(j)%first_producer;d=sp%patch(j)%last_producer
             dp%patch(j)%first_pressure=nn+1;dp%patch(j)%first_producer=np+1
             n=count(owner(sp%node(a:b))==rank)
             if(pass==2) dp%node(nn+1:nn+n)=pack(sp%node(a:b),owner(sp%node(a:b))==rank)
             nn=nn+n
             keep=owner(src%producer(c:d)%target)==rank;n=count(keep)
             if(pass==2) then
                dst%producer(np+1:np+n)=pack(src%producer(c:d),keep)
                dp%producer_ordinal(np+1:np+n)=pack(sp%producer_ordinal(c:d),keep)
             end if
             np=np+n;dp%patch(j)%last_pressure=nn;dp%patch(j)%last_producer=np
          end do
          if(pass==1) allocate(dp%node(nn),dst%producer(np),dp%producer_ordinal(np))
       end do
       dst%restriction=pack(src%restriction,owner(src%restriction%target)==rank)
       dst%divergence=pack(src%divergence,owner(src%divergence%target)==rank)
       end associate
       associate(src=>input%rotation(l),dst=>output%rotation(l))
       dst%curl=pack(src%curl,owner(src%curl%target)==rank)
       if(allocated(keep)) deallocate(keep)
       allocate(keep(size(src%pentagon)))
       do s=1,size(keep)
          keep(s)=owner(src%pentagon(s)%mass(1))==rank
       end do
       dst%pentagon=pack(src%pentagon,keep)
       end associate
       associate(src=>input%dynamics(l),dst=>output%dynamics(l))
       dst%divergence=pack(src%divergence,owner(src%divergence%target)==rank)
       if(allocated(keep)) deallocate(keep)
       allocate(keep(size(src%physics)))
       do s=1,size(keep)
          keep(s)=owner(src%physics(s)%node(1))==rank
       end do
       dst%physics=pack(src%physics,keep)
       if(allocated(index)) deallocate(index)
       allocate(index(size(src%injection,2)))
       do s=1,size(index)
          index(s)=s
       end do
       index=pack(index,owner(src%injection(1,:))==rank)
       dst%injection=src%injection(:,index)
       if(allocated(dmap)) deallocate(dmap,rmap)
       allocate(dmap(size(src%source%direct)),rmap(size(src%source%restriction)))
       dmap=0;rmap=0;n=0
       do s=1,size(dmap)
          if(owner(src%source%direct(s)%node(1))/=rank) cycle
          n=n+1;dmap(s)=n
       end do
       dst%source%direct=pack(src%source%direct,dmap>0)
       n=0
       do s=1,size(rmap)
          if(owner(src%source%restriction(s)%target)/=rank) cycle
          n=n+1;rmap(s)=n
       end do
       dst%source%restriction=pack(src%source%restriction,rmap>0)
       if(allocated(actions)) deallocate(actions,operands)
       allocate(actions(size(src%source%action)),operands(size(src%source%operand)))
       n=0
       do s=1,size(src%source%action)
          i=src%source%operand(s)
          select case(src%source%action(s))
          case(VELOCITY_DIRECT)
             i=dmap(i)
          case(VELOCITY_RESTRICT)
             i=rmap(i)
          case default
             error stop 'RHS placement unknown source action'
          end select
          if(i==0) cycle
          n=n+1;actions(n)=src%source%action(s);operands(n)=i
       end do
       dst%source%action=actions(:n);dst%source%operand=operands(:n)
       end associate
    end do
    allocate(seen(nold),map(nold));seen=owner==rank;map=0
    do pass=1,2
       do s=1,size(output%gradient)
          output%gradient(s)%node=remap(output%gradient(s)%node,1)
       end do
       do l=lbound(output%diffusion,1),ubound(output%diffusion,1)
          output%pressure(l)%node=remap(output%pressure(l)%node,1)
          do s=1,size(output%diffusion(l)%producer)
             associate(p=>output%diffusion(l)%producer(s))
             p%target=remap(p%target,1);p%mass=remap(p%mass,1)
             p%velocity=remap(p%velocity,3);p%vorticity=remap(p%vorticity,2)
             if(p%southwest) then
                p%sw_mass=remap(p%sw_mass,1);p%sw_velocity=remap(p%sw_velocity,3)
                p%sw_edge(2,:)=remap(p%sw_edge(2,:),1);p%sw_vorticity=remap(p%sw_vorticity,2)
             end if
             end associate
          end do
          do s=1,size(output%diffusion(l)%restriction)
             associate(p=>output%diffusion(l)%restriction(s))
             p%target=remap(p%target,1);p%small=remap(p%small,3)
             do i=1,3
                if(.not.p%edge(i)) cycle
                p%partial_flux(:,i)=remap(p%partial_flux(:,i),3)
                p%partial_node(:,i)=remap(p%partial_node(:,i),1)
                p%coarse_node(:,i)=remap(p%coarse_node(:,i),1)
             end do
             end associate
          end do
          do s=1,size(output%diffusion(l)%divergence)
             associate(p=>output%diffusion(l)%divergence(s))
             p%target=remap(p%target,1);p%edge=remap(p%edge,3)
             end associate
          end do
          do s=1,size(output%rotation(l)%curl)
             associate(p=>output%rotation(l)%curl(s))
             p%target=remap(p%target,1);p%vorticity=remap(p%vorticity,2)
             end associate
          end do
          do s=1,size(output%rotation(l)%pentagon)
             associate(p=>output%rotation(l)%pentagon(s))
             p%mass=remap(p%mass,1);p%velocity=remap(p%velocity,3)
             p%vorticity=remap(p%vorticity,2);p%pv_node=remap(p%pv_node,1)
             end associate
          end do
          associate(p=>output%dynamics(l))
          p%injection=remap(p%injection,1)
          do s=1,size(p%divergence)
             p%divergence(s)%target=remap(p%divergence(s)%target,1)
             if(p%divergence(s)%active) p%divergence(s)%edge=remap(p%divergence(s)%edge,3)
          end do
          do s=1,size(p%physics)
             p%physics(s)%node=remap(p%physics(s)%node,1)
             p%physics(s)%vorticity=remap(p%physics(s)%vorticity,2)
          end do
          do s=1,size(p%source%direct)
             p%source%direct(s)%node=remap(p%source%direct(s)%node,1)
          end do
          do s=1,size(p%source%restriction)
             p%source%restriction(s)%target=remap(p%source%restriction(s)%target,1)
             p%source%restriction(s)%child=remap(p%source%restriction(s)%child,1)
             p%source%restriction(s)%neighbor=remap(p%source%restriction(s)%neighbor,1)
          end do
          end associate
       end do
       if(pass==1) then
          n=count(seen);allocate(original(n),owned(n));j=0
          do i=1,nold
             if(.not.seen(i)) cycle
             j=j+1;map(i)=j;original(j)=i;owned(j)=owner(i)==rank
          end do
       end if
    end do
    output%nnode=size(original)
    output%surface_visits=output%surface_visits(original)
  contains
    impure elemental integer function remap(address,mult) result(value)
      integer,intent(in) :: address,mult
      integer :: node
      value=address
      if(address==0) return
      node=(address-1)/mult+1
      if(address<0.or.node>nold) error stop 'RHS placement invalid physical address'
      if(pass==1) then
         seen(node)=.true.
      else
         if(map(node)==0) error stop 'RHS placement missing compact address'
         value=mult*(map(node)-1)+mod(address-1,mult)+1
      end if
    end function
  end subroutine
end module
