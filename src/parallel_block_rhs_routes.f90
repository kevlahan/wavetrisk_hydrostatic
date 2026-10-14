module parallel_block_rhs_routes_mod
  use iso_fortran_env, only : int64
  use mpi_f08, only : MPI_Comm
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  use parallel_block_rhs_exchange_mod, only : Rhs_Exchange,compile_rhs_exchange
  use parallel_block_rhs_dynamics_mod
  use parallel_block_rhs_diffusion_mod
  implicit none
  private
  public :: compile_rhs_phase_route
contains
  subroutine compile_rhs_phase_route(plan,key,owned,phase,l,ncolumn,comm,topology,masks, &
       consumer_key,orientation,route,shared_context,producer_home,consumer_home)
    type(Rhs_Domain_Operators),intent(in) :: plan
    integer,intent(in) :: key(:,:),phase,l,ncolumn,consumer_key(:,:),orientation(:)
    logical,intent(in) :: owned(:)
    type(MPI_Comm),intent(in) :: comm
    integer(int64),intent(in) :: topology,masks
    type(Rhs_Exchange),intent(inout) :: route
    logical,optional,intent(in) :: shared_context
    integer,optional,intent(in) :: producer_home(:,:),consumer_home(:,:)
    logical,allocatable :: needed(:)
    integer,allocatable :: definition(:,:),order(:,:),slot(:),target(:),canonical(:)
    integer :: mult,n,s,c,i,pass,ordinal,pkey(4),write_order(2),nslot
    if(size(key,1)/=4.or.size(key,2)/=plan%nnode.or.size(owned)/=plan%nnode) &
         error stop 'RHS phase directory node extent'
    mult=1
    select case(phase)
    case(RHS_DIFFUSION_EDGES,RHS_ROTATION_EDGES,RHS_FINAL_PV,RHS_SCALAR_FLUX,RHS_VELOCITY_SOURCE, &
         RHS_DIFFUSION_DIRECT_EDGES,RHS_DIRECT_SCALAR_FLUX)
       mult=3
    case(RHS_PRIMITIVE_VORT,RHS_ROTATION_TRIANGLES,RHS_ROTATION_DIRECT_TRIANGLES)
       mult=2
    end select
    nslot=mult*plan%nnode
    if(any(shape(consumer_key)/=[4,nslot]).or.size(orientation)/=nslot) error stop 'RHS phase consumer extent'
    call phase_demand(plan,owned,phase,mult,needed)
    do pass=1,2
       n=0
       ! A single numerical owner preserves slots not written in this phase.
       ! Direct/alias definitions supersede these baseline values by ordinal.
       do s=1,plan%nnode
          if(.not.owned(s)) cycle
          do c=1,mult
             call define(s,c,[0,1])
          end do
       end do
       select case(phase)
       case(RHS_DIFFUSION_EDGES,RHS_SCALAR_FLUX,RHS_FINAL_PV,RHS_PRIMITIVE_VORT,RHS_ROTATION_TRIANGLES, &
            RHS_DIFFUSION_DIRECT_EDGES,RHS_DIRECT_SCALAR_FLUX,RHS_ROTATION_DIRECT_TRIANGLES)
          do s=1,size(plan%diffusion(l)%producer)
             associate(p=>plan%diffusion(l)%producer(s))
             ordinal=plan%pressure(l)%producer_ordinal(s)
             write_order=[1,2*ordinal]
             if(mult==3) then
                do c=1,3
                   call define(p%target,c,write_order)
                end do
             else
                do c=1,4
                   call define_component(p%vorticity(c),write_order)
                end do
             end if
             if(p%southwest) then
                write_order=[1,2*ordinal+1]
                if(mult==3) then
                   do c=1,3
                      call define(p%sw_edge(2,c),p%sw_edge(1,c),write_order)
                   end do
                else
                   do c=1,2
                      call define_component(p%sw_vorticity(c),write_order)
                   end do
                end if
             end if
             end associate
          end do
          if(phase==RHS_DIFFUSION_EDGES.or.phase==RHS_SCALAR_FLUX) then
             do s=1,size(plan%diffusion(l)%restriction)
                do c=1,3
                   if(.not.plan%diffusion(l)%restriction(s)%edge(c)) cycle
                   call define(plan%diffusion(l)%restriction(s)%target,c,[2,s])
                end do
             end do
          end if
          if(phase==RHS_FINAL_PV.or.phase==RHS_ROTATION_TRIANGLES) then
             do s=1,size(plan%rotation(l)%pentagon)
                associate(p=>plan%rotation(l)%pentagon(s))
                write_order=[3,p%corner]
                select case(p%corner)
                case(7)
                   if(mult==3) then
                      call define(p%pv_node(5),1,write_order);call define(p%pv_node(4),3,write_order)
                   else
                      call define_component(p%vorticity(1,7),write_order)
                      call define_component(p%vorticity(2,7),write_order)
                   end if
                case(6)
                   if(mult==3) then
                      call define(p%pv_node(1),1,write_order);call define(p%pv_node(7),2,write_order)
                   else
                      call define_component(p%vorticity(1,7),write_order)
                      call define_component(p%vorticity(2,4),write_order)
                   end if
                case(8)
                   if(mult==3) then
                      call define(p%pv_node(1),3,write_order);call define(p%pv_node(7),2,write_order)
                   else
                      call define_component(p%vorticity(2,7),write_order)
                      call define_component(p%vorticity(1,5),write_order)
                   end if
                case(5)
                   if(mult==3) then
                      call define(p%pv_node(1),1,write_order);call define(p%pv_node(1),3,write_order)
                   else
                      call define_component(p%vorticity(1,1),write_order)
                      call define_component(p%vorticity(2,1),write_order)
                   end if
                case default
                   error stop 'RHS phase invalid pentagon'
                end select
                end associate
             end do
          end if
       case(RHS_ROTATION_EDGES)
          do s=1,size(plan%rotation(l)%curl)
             do c=1,3
                call define(plan%rotation(l)%curl(s)%target,c,[1,s])
             end do
          end do
       case(RHS_DIFFUSION_NODES,RHS_DIFFUSION_DIRECT_NODES)
          do s=1,size(plan%diffusion(l)%divergence)
             call define(plan%diffusion(l)%divergence(s)%target,1,[1,s])
          end do
       case(RHS_SCALAR_TENDENCY)
          do s=1,size(plan%dynamics(l)%divergence)
             call define(plan%dynamics(l)%divergence(s)%target,1,[1,s])
          end do
       case(RHS_VELOCITY_SOURCE,RHS_BERNOULLI,RHS_EXNER,RHS_INJECTED_BERNOULLI,RHS_INJECTED_EXNER)
          ! These operators write only their target owner's canonical slot.
       case default
          error stop 'RHS phase directory unknown phase'
       end select
       if(pass==1) allocate(definition(4,n),order(2,n),slot(n))
    end do
    needed(slot)=.true.
    allocate(target(count(needed)));target=pack([(i,i=1,nslot)],needed)
    if(present(producer_home)) then
       if(.not.present(consumer_home).or..not.present(shared_context)) error stop 'RHS incomplete route home binding'
       if(.not.shared_context) error stop 'RHS route home context must be owned by executor'
       if(any(shape(producer_home)/=[2,nslot]).or.any(shape(consumer_home)/=[2,nslot])) &
            error stop 'RHS route home extent'
       canonical=pack([(i,i=1,nslot)],[(owned((i-1)/mult+1),i=1,nslot)])
       call compile_rhs_exchange(route,comm,topology,masks,100*l+phase,nslot,ncolumn,definition,order,slot, &
            consumer_key(:,target),target,orientation(target),shared_tag=3*(100*l+phase)+ncolumn, &
            variable_columns=.true.,definition_home=producer_home(:,slot),consumer_home=consumer_home(:,target), &
            canonical_slot=canonical)
    else if(present(shared_context)) then
       if(.not.shared_context) error stop 'RHS shared context must be owned by executor'
       call compile_rhs_exchange(route,comm,topology,masks,100*l+phase,nslot,ncolumn,definition,order,slot, &
            consumer_key(:,target),target,orientation(target),shared_tag=3*(100*l+phase)+ncolumn,variable_columns=.true.)
    else
       call compile_rhs_exchange(route,comm,topology,masks,100*l+phase,nslot,ncolumn,definition,order,slot, &
            consumer_key(:,target),target,orientation(target))
    end if
  contains
    subroutine define(node,component,sequence)
      integer,intent(in) :: node,component,sequence(2)
      if(node<1.or.node>plan%nnode.or.component<1.or.component>mult) error stop 'RHS phase definition address'
      n=n+1
      if(pass==1) return
      pkey=key(:,node);pkey(3)=component
      definition(:,n)=pkey;order(:,n)=sequence+[1,0];slot(n)=mult*(node-1)+component
    end subroutine
    subroutine define_component(address,sequence)
      integer,intent(in) :: address,sequence(2)
      call define((address-1)/mult+1,mod(address-1,mult)+1,sequence)
    end subroutine
  end subroutine
  subroutine phase_demand(plan,owned,phase,mult,needed)
    ! A field's reader/writer union across levels is conservative about later
    ! coarse reads and aliases. Unlike the old full-program refresh, unrelated
    ! fields' ghost slots are not dependencies. Keep canonical owners coherent
    ! because they are fallback providers in later phases, even without readers.
    type(Rhs_Domain_Operators),intent(in) :: plan
    logical,intent(in) :: owned(:)
    integer,intent(in) :: phase,mult
    logical,allocatable,intent(out) :: needed(:)
    integer :: l,s,c,n
    allocate(needed(mult*plan%nnode));needed=.false.
    do n=1,plan%nnode
       if(owned(n)) call nodes([n])
    end do
    do l=lbound(plan%diffusion,1),ubound(plan%diffusion,1)
       select case(phase)
       case(RHS_BERNOULLI,RHS_EXNER,RHS_INJECTED_BERNOULLI,RHS_INJECTED_EXNER)
          call nodes(reshape(plan%dynamics(l)%injection,[size(plan%dynamics(l)%injection)]))
       case(RHS_FINAL_PV)
          do s=1,size(plan%dynamics(l)%source%direct)
             if(plan%dynamics(l)%source%direct(s)%active) call nodes(plan%dynamics(l)%source%direct(s)%node)
          end do
          do s=1,size(plan%diffusion(l)%producer)
             associate(p=>plan%diffusion(l)%producer(s))
             call nodes([p%target])
             if(p%southwest) call slots(3*(p%sw_edge(2,:)-1)+p%sw_edge(1,:))
             end associate
          end do
          do s=1,size(plan%rotation(l)%pentagon)
             call nodes(plan%rotation(l)%pentagon(s)%pv_node)
          end do
       case(RHS_VELOCITY_SOURCE)
          do s=1,size(plan%dynamics(l)%source%restriction)
             associate(p=>plan%dynamics(l)%source%restriction(s))
             do c=1,3
                if(p%edge(c)) call slots([3*(p%target-1)+c,3*(p%child-1)+c,3*(p%neighbor(c)-1)+c])
             end do
             end associate
          end do
       case(RHS_PRIMITIVE_VORT,RHS_ROTATION_TRIANGLES,RHS_ROTATION_DIRECT_TRIANGLES)
          do s=1,size(plan%rotation(l)%curl)
             call slots(plan%rotation(l)%curl(s)%vorticity)
          end do
          do s=1,size(plan%dynamics(l)%physics)
             call slots(plan%dynamics(l)%physics(s)%vorticity)
          end do
          do s=1,size(plan%diffusion(l)%producer)
             call slots(plan%diffusion(l)%producer(s)%vorticity)
             if(plan%diffusion(l)%producer(s)%southwest) call slots(plan%diffusion(l)%producer(s)%sw_vorticity)
          end do
          do s=1,size(plan%rotation(l)%pentagon)
             call slots(reshape(plan%rotation(l)%pentagon(s)%vorticity,[14]))
          end do
       case(RHS_ROTATION_EDGES)
          do s=1,size(plan%diffusion(l)%producer)
             call slots(plan%diffusion(l)%producer(s)%velocity)
             if(plan%diffusion(l)%producer(s)%southwest) call slots(plan%diffusion(l)%producer(s)%sw_velocity)
          end do
          do s=1,size(plan%rotation(l)%pentagon)
             call slots(reshape(plan%rotation(l)%pentagon(s)%velocity,[21]))
          end do
          do s=1,size(plan%rotation(l)%curl)
             call nodes([plan%rotation(l)%curl(s)%target])
          end do
       case(RHS_DIFFUSION_EDGES,RHS_DIFFUSION_DIRECT_EDGES,RHS_SCALAR_FLUX,RHS_DIRECT_SCALAR_FLUX)
          do s=1,size(plan%diffusion(l)%restriction)
             associate(p=>plan%diffusion(l)%restriction(s))
             call slots(reshape(p%small,[size(p%small)]))
             do c=1,3
                if(.not.p%edge(c)) cycle
                call slots(p%partial_flux(:,c));call slots([3*(p%target-1)+c])
             end do
             end associate
          end do
          do s=1,size(plan%diffusion(l)%divergence)
             call slots(plan%diffusion(l)%divergence(s)%edge)
          end do
          if(phase==RHS_SCALAR_FLUX.or.phase==RHS_DIRECT_SCALAR_FLUX) then
             do s=1,size(plan%dynamics(l)%divergence)
                if(plan%dynamics(l)%divergence(s)%active) call slots(plan%dynamics(l)%divergence(s)%edge)
             end do
             do s=1,size(plan%dynamics(l)%source%direct)
                if(plan%dynamics(l)%source%direct(s)%active) call nodes(plan%dynamics(l)%source%direct(s)%node)
             end do
          end if
          do s=1,size(plan%diffusion(l)%producer)
             associate(p=>plan%diffusion(l)%producer(s))
             call nodes([p%target])
             if(p%southwest) call slots(3*(p%sw_edge(2,:)-1)+p%sw_edge(1,:))
             end associate
          end do
       case(RHS_DIFFUSION_NODES,RHS_DIFFUSION_DIRECT_NODES,RHS_SCALAR_TENDENCY)
          do s=1,size(plan%diffusion(l)%restriction)
             associate(p=>plan%diffusion(l)%restriction(s))
             do c=1,3
                if(.not.p%edge(c)) cycle
                call nodes(p%partial_node(:,c));call nodes(p%coarse_node(:,c))
             end do
             end associate
          end do
          if(phase/=RHS_SCALAR_TENDENCY) then
             do s=1,size(plan%diffusion(l)%producer)
                associate(p=>plan%diffusion(l)%producer(s))
                call nodes(p%mass([1,3,6,2]))
                if(p%southwest) call nodes([p%mass(1),p%mass(5),p%sw_mass,p%mass(4)])
                end associate
             end do
             do s=1,size(plan%dynamics(l)%physics)
                call nodes(plan%dynamics(l)%physics(s)%node)
             end do
          end if
       case default
          error stop 'RHS demand compiler unknown phase'
       end select
    end do
    if(phase==RHS_BERNOULLI.or.phase==RHS_EXNER.or.phase==RHS_INJECTED_BERNOULLI.or. &
         phase==RHS_INJECTED_EXNER.or.phase==RHS_VELOCITY_SOURCE) then
       do s=1,size(plan%gradient)
          if(.not.plan%gradient(s)%active) cycle
          if(phase==RHS_VELOCITY_SOURCE) then
             call nodes([plan%gradient(s)%node(1)])
          else
             call nodes(plan%gradient(s)%node)
          end if
       end do
    end if
  contains
    subroutine nodes(address)
      integer,intent(in) :: address(:)
      integer :: i,j
      do i=1,size(address)
         j=address(i)
         if(j<1.or.j>plan%nnode) error stop 'RHS demand node extent'
         needed(mult*(j-1)+1:mult*j)=.true.
      end do
    end subroutine
    subroutine slots(address)
      integer,intent(in) :: address(:)
      if(any(address<1).or.any(address>size(needed))) error stop 'RHS demand component extent'
      needed(address)=.true.
    end subroutine
  end subroutine
end module
