module parallel_block_rhs_distributed_mod
  use audit316_mod, only : audit316_enabled,audit316_i,audit316_i64,audit316_r,audit316_l, &
       audit316_begin,audit316_end
  use parallel_block_profile_mod, only : detail_enter,detail_leave,rhs_detail_enter,rhs_detail_leave, &
       DP_RHS_GEOMETRY,DP_RHS_ROUTES,DP_RHS_COST
  ! Connected direct execution on target block owners. The state adapter must
  ! supply this compact input layout at the declared epoch before run entry.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  use arch_mod, only : comm,rank
  use shared_mod, only : AT_NODE,AT_EDGE,level_start,level_end,Laplace_sclr,Laplace_divu
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators,audit_rhs_program316
  use parallel_block_rhs_distribute_mod
  use parallel_block_rhs_boundary_mod
  use parallel_block_rhs_routes_mod
  use parallel_block_rhs_exchange_mod
  use parallel_block_rhs_diffusion_mod
  use parallel_block_rhs_dynamics_mod
  implicit none
  private
  public :: Rhs_Distributed_Executor,prepare_rhs_distributed,run_rhs_distributed,rhs_distributed_halo,release_rhs_distributed
  public :: audit_rhs_distributed316
  type :: Rhs_Route_Context
     type(MPI_Comm) :: comm=MPI_COMM_NULL
     type(Rhs_Exchange),allocatable :: route(:,:)
     logical,allocatable :: bound(:,:)
     integer,allocatable :: route_phase(:,:),route_level(:,:)
     real(dp),allocatable :: scratch(:)
     integer,allocatable :: audit_event316(:,:)
     integer :: audit_used316=0,audit_layer316=0
     integer(int64) :: topology=-1,masks=-1,state=0,exchanges=0,remote_values=0,remote_messages=0
  end type
  type :: Rhs_Distributed_Executor
     type(Rhs_Placement) :: placement
     type(Rhs_Dynamics_Workspace) :: work
     type(Rhs_Route_Context) :: transport
     real(dp),allocatable :: scalar_rhs(:,:,:),velocity_rhs(:,:)
     integer(int64) :: columns=0,primitive_events=0,scalar_targets=0,velocity_targets=0
  end type
  logical,save :: demand_option_initialized=.false.,demand_enabled=.false.
contains
  subroutine audit_rhs_distributed316(executor)
    type(Rhs_Distributed_Executor),intent(in) :: executor
    integer :: l,p
    call audit316_i('placement.key.shape',shape(executor%placement%key))
    call audit316_i('placement.key',reshape(executor%placement%key,[size(executor%placement%key)]))
    call audit316_i('placement.base',executor%placement%base)
    call audit316_l('placement.owned',executor%placement%owned)
    call audit_rhs_program316(executor%placement%program)
    associate(t=>executor%transport)
    call audit316_i64('transport.epoch',[t%topology,t%masks,t%state])
    call audit316_i('transport.route.bounds',[lbound(t%route,2),ubound(t%route,2)])
    call audit316_l('transport.bound',reshape(t%bound,[size(t%bound)]))
    call audit316_i('transport.route_phase',reshape(t%route_phase,[size(t%route_phase)]))
    call audit316_i('transport.route_level',reshape(t%route_level,[size(t%route_level)]))
    call audit316_i('transport.scratch.shape',shape(t%scratch))
    do l=lbound(t%route,2),ubound(t%route,2)
       do p=1,size(t%route,1)
          if(.not.t%bound(p,l)) cycle
          call audit316_i('transport.route.coordinate',[p,l])
          call audit_rhs_exchange316(t%route(p,l))
       end do
    end do
    end associate
    call audit316_i('executor.scalar_rhs.shape',shape(executor%scalar_rhs))
    call audit316_i('executor.velocity_rhs.shape',shape(executor%velocity_rhs))
    call audit316_i('workspace.lower.shape',shape(executor%work%lower))
    call audit316_i('workspace.phi.shape',shape(executor%work%phi))
    call audit316_i('workspace.pressure.shape',shape(executor%work%pressure))
    call audit316_i('workspace.midpoint.shape',shape(executor%work%midpoint))
    call audit316_i('workspace.exner.shape',shape(executor%work%exner))
    call audit316_i('workspace.thermodynamic.shape',shape(executor%work%thermodynamic))
    call audit316_i('workspace.flux.shape',shape(executor%work%flux))
    call audit316_i('workspace.pv.shape',shape(executor%work%pv))
    call audit316_i('workspace.ke.shape',shape(executor%work%ke))
    call audit316_i('workspace.bernoulli.shape',shape(executor%work%bernoulli))
    call audit316_i('workspace.vorticity.shape',shape(executor%work%vorticity))
    call audit316_i('workspace.diffusion_flux.shape',shape(executor%work%diffusion_flux))
    call audit316_i('workspace.laplacian.shape',shape(executor%work%laplacian))
    call audit316_i('workspace.divu.shape',shape(executor%work%divu))
    call audit316_i('workspace.lap_divu.shape',shape(executor%work%lap_divu))
    call audit316_i('workspace.rotation.shape',shape(executor%work%rotation))
    call audit316_i('workspace.source.shape',shape(executor%work%source))
    call audit316_i('workspace.physics.shape',shape(executor%work%physics))
    call audit316_i('workspace.scalar_rhs.shape',shape(executor%work%scalar_rhs))
    call audit316_i('workspace.velocity_rhs.shape',shape(executor%work%velocity_rhs))
    call audit316_i('workspace.scalar_physics.shape',shape(executor%work%scalar_physics))
    call audit316_i('workspace.sw_physics.shape',shape(executor%work%sw_physics))
  end subroutine

  logical function prefix_demand_enabled(context) result(enabled)
    type(MPI_Comm),intent(in) :: context
    integer :: local(2),global(2),ierr
    character(8) :: option
    if(.not.demand_option_initialized) then
       call get_environment_variable('WAVETRISK_RHS_PREFIX_DEMAND',option)
       if(len_trim(option)>0.and.trim(option)/='0'.and.trim(option)/='1') error stop 'prefix demand selector'
       demand_enabled=trim(option)=='1'
       local=[merge(1,0,demand_enabled),merge(0,1,demand_enabled)]
       call MPI_Allreduce(local,global,2,MPI_INTEGER,MPI_MAX,context,ierr)
       if(ierr/=MPI_SUCCESS.or.all(global==1)) error stop 'inconsistent prefix demand selector'
       demand_option_initialized=.true.
    end if
    enabled=demand_enabled
  end function

  subroutine compile_prefix_demand(executor)
    ! Backward live-variable analysis of the COMPLETE ordered prefix, including
    ! all actual numerical reads, aliases, fallback providers and retained values
    ! across layers/RHS calls. The resulting sparse routes execute unchanged.
    type(Rhs_Distributed_Executor),intent(inout) :: executor
    type :: Demand
       logical,allocatable :: target(:)
    end type
    type(Demand) :: keep(20,level_start:level_end)
    logical,allocatable :: seed(:,:,:),carry_node(:,:,:),live_node(:,:,:),carry_flux(:,:),live_flux(:,:)
    logical :: changed,global_changed,converged
    integer :: n,p,l,s,c,field,nc,iteration,ierr,a
    integer(int64) :: before_target,after_target,before_remote,after_remote,sv,rv,bytes,temp_bound
    integer :: sp,pr
    character(8) :: option
    n=executor%placement%program%nnode
    allocate(seed(n,2,3),carry_node(n,2,3),live_node(n,2,3),carry_flux(3*n,2),live_flux(3*n,2))
    seed=.false.;carry_node=.false.;carry_flux=.false.;before_target=0;before_remote=0
    temp_bound=4_int64*(30_int64*n)
    associate(t=>executor%transport,plan=>executor%placement%program)
    do l=level_start,level_end
       do p=1,20
          if(.not.t%bound(p,l)) cycle
          if(p/=RHS_DIFFUSION_EDGES.and.p/=RHS_DIFFUSION_NODES.and. &
               p/=RHS_DIFFUSION_DIRECT_EDGES.and.p/=RHS_DIFFUSION_DIRECT_NODES) cycle
          allocate(keep(p,l)%target(rhs_exchange_target_count(t%route(p,l))));keep(p,l)%target=.false.
          before_target=before_target+size(keep(p,l)%target)
          temp_bound=temp_bound+4_int64*size(keep(p,l)%target)
          call rhs_exchange_cost(t%route(p,l),sv,rv,sp,pr,bytes);before_remote=before_remote+sv
       end do
       if(Laplace_sclr==2) then
          do s=1,size(plan%diffusion(l)%producer)
             associate(q=>plan%diffusion(l)%producer(s))
             seed(q%mass([1,3,6,2]),:,2)=.true.
             if(q%southwest) seed([q%mass(1),q%mass(5),q%sw_mass,q%mass(4)],:,2)=.true.
             end associate
          end do
       end if
       if(Laplace_divu/=0) then
          field=merge(1,3,Laplace_divu==1)
          do s=1,size(plan%dynamics(l)%physics)
             seed(plan%dynamics(l)%physics(s)%node,1,field)=.true.
          end do
       end if
    end do
    converged=.false.
    do iteration=1,12
       live_node=seed.or.carry_node;live_flux=carry_flux
       ! Reverse third Laplacian(divu), second scalar Laplacian, first divu.
       do field=3,1,-1
          if(field==3.and.Laplace_divu/=2) cycle
          if(field==2.and.Laplace_sclr/=2) cycle
          if(field==1.and.Laplace_divu==0) cycle
          nc=merge(2,1,field==2)
          do l=level_start,level_end
             p=RHS_DIFFUSION_NODES
             if(field==1) p=RHS_DIFFUSION_DIRECT_NODES
             call reverse_halo(p,l,live_node(:,:nc,field))
             ! Divergence overwrites its targets, reading all six incident edges.
             do s=1,size(plan%diffusion(l)%divergence)
                live_node(plan%diffusion(l)%divergence(s)%target,:nc,field)=.false.
                live_flux(plan%diffusion(l)%divergence(s)%edge,:nc)=.true.
             end do
             call reverse_halo(RHS_DIFFUSION_EDGES,l,live_flux(:,:nc))
             ! Restriction is in-place: reverse edge writes, then its earlier
             ! small-flux reads. Unused arithmetic still receives exact inputs.
             do s=size(plan%diffusion(l)%restriction),1,-1
                associate(q=>plan%diffusion(l)%restriction(s))
                do c=3,1,-1
                   if(.not.q%edge(c)) cycle
                   live_flux(3*(q%target-1)+c,:nc)=.false.
                   live_flux(q%partial_flux(:,c),:nc)=.true.
                   live_node(q%partial_node(:,c),:nc,field)=.true.
                   live_node(q%coarse_node(:,c),:nc,field)=.true.
                end do
                live_flux(reshape(q%small,[size(q%small)]),:nc)=.true.
                end associate
             end do
             call reverse_halo(RHS_DIFFUSION_DIRECT_EDGES,l,live_flux(:,:nc))
             ! Direct/SW edge writes use resident inputs, or the completed divu
             ! field for the third sweep. They do not read the reused flux array.
             do s=1,size(plan%diffusion(l)%producer)
                associate(q=>plan%diffusion(l)%producer(s))
                a=3*(q%target-1);live_flux(a+1:a+3,:nc)=.false.
                if(q%southwest) then
                   do c=1,3
                      a=3*(q%sw_edge(2,c)-1)+q%sw_edge(1,c);live_flux(a,:nc)=.false.
                   end do
                end if
                if(field==3) then
                   live_node(q%mass([1,3,6,2]),1,1)=.true.
                   if(q%southwest) live_node([q%mass(1),q%mass(5),q%sw_mass,q%mass(4)],1,1)=.true.
                end if
                end associate
             end do
          end do
       end do
       changed=any(live_node.and..not.carry_node).or.any(live_flux.and..not.carry_flux)
       call MPI_Allreduce(changed,global_changed,1,MPI_LOGICAL,MPI_LOR,t%comm,ierr)
       if(ierr/=MPI_SUCCESS) error stop 'prefix retained-version closure'
       if(.not.global_changed) then
          converged=.true.;exit
       end if
       carry_node=carry_node.or.live_node;carry_flux=carry_flux.or.live_flux
    end do
    after_target=before_target;after_remote=before_remote
    if(converged) then
       after_target=0;after_remote=0
       do l=level_start,level_end
          do p=1,20
             if(.not.allocated(keep(p,l)%target)) cycle
             call rhs_exchange_restrict(t%route(p,l),keep(p,l)%target)
             after_target=after_target+rhs_exchange_target_count(t%route(p,l))
             call rhs_exchange_cost(t%route(p,l),sv,rv,sp,pr,bytes);after_remote=after_remote+sv
          end do
       end do
    end if
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
    if(trim(option)=='1') write(*,'(a,8(i0,1x))') 'Demand287 preparation: ',rank,min(iteration,12), &
         merge(1,0,converged),before_target,after_target,before_remote,after_remote,temp_bound
    end associate
  contains
    subroutine reverse_halo(phase,level,live)
      integer,intent(in) :: phase,level
      logical,intent(inout) :: live(:,:)
      integer :: physical_phase,physical_level
      physical_phase=executor%transport%route_phase(phase,level)
      physical_level=executor%transport%route_level(phase,level)
      call rhs_exchange_backward_demand(executor%transport%route(physical_phase,physical_level),live, &
           keep(physical_phase,physical_level)%target)
    end subroutine
  end subroutine

  subroutine release_rhs_distributed(executor)
    type(Rhs_Distributed_Executor),intent(inout) :: executor
    integer :: p,l,ierr
    if(allocated(executor%transport%route)) then
       do l=lbound(executor%transport%route,2),ubound(executor%transport%route,2)
          do p=1,size(executor%transport%route,1)
             if(executor%transport%bound(p,l)) call release_rhs_exchange(executor%transport%route(p,l))
          end do
       end do
    end if
    if(executor%transport%comm/=MPI_COMM_NULL) then
       call MPI_Comm_free(executor%transport%comm,ierr)
       if(ierr/=MPI_SUCCESS) error stop 'RHS context release'
    end if
    executor=Rhs_Distributed_Executor()
  end subroutine

  subroutine prepare_rhs_distributed(input,base,topology,masks,nlayer,executor,directory)
    type(Rhs_Domain_Operators),intent(in) :: input(:)
    integer,intent(in) :: base(:),nlayer
    integer(int64),intent(in) :: topology,masks
    type(Rhs_Distributed_Executor),intent(inout) :: executor
    type(Rhs_Geometry_Directory),optional,intent(in) :: directory
    type(Rhs_Boundary_Map) :: mapping
    type :: Route_Keys
       integer,allocatable :: key(:,:),orientation(:),home(:,:)
    end type
    type(Route_Keys) :: direct(3),boundary(2)
    type(Rhs_Owner_Directory) :: owner_directory
    integer,allocatable :: owner_slot(:),ghost(:),home(:,:),resolved(:,:),request(:,:),changed(:)
    real(dp),allocatable :: metadata(:,:),received(:,:)
    integer,parameter :: phases(15)=[1,2,3,4,5,6,7,11,12,13,14,15,16,18,20]
    integer :: p,l,c,m,i,j,n,location,nc,ierr,a,rp,rl,capacity,logical_buffers,direct_routes
    character(len=8) :: audit_option
    call release_rhs_distributed(executor)
    call detail_enter(DP_RHS_GEOMETRY)
    call distribute_rhs_geometry(input,executor%placement)
    call detail_leave(DP_RHS_GEOMETRY)
    call detail_enter(DP_RHS_ROUTES)
    n=executor%placement%program%nnode
    call initialize_rhs_dynamics(executor%placement%program,nlayer,executor%work)
    allocate(executor%scalar_rhs(n,2,nlayer),executor%velocity_rhs(3*n,nlayer))
    executor%scalar_rhs=0;executor%velocity_rhs=0
    associate(t=>executor%transport)
    call MPI_Comm_dup(comm,t%comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'RHS context creation'
    allocate(t%route(20,level_start:level_end),t%bound(20,level_start:level_end), &
         t%route_phase(20,level_start:level_end),t%route_level(20,level_start:level_end))
    t%route_phase=0;t%route_level=0;capacity=0;logical_buffers=0;direct_routes=0
    t%bound=.false.;t%topology=topology;t%masks=masks
    owner_slot=pack([(j,j=1,n)],executor%placement%owned)
    call compile_rhs_owner_directory(owner_directory,t%comm,executor%placement%key(:,owner_slot),owner_slot)
    allocate(home(2,n));home(1,:)=rank;home(2,:)=[(j,j=1,n)]
    ghost=pack([(j,j=1,n)],.not.executor%placement%owned)
    call lookup_rhs_owner_directory(owner_directory,executor%placement%key(:,ghost),resolved)
    home(:,ghost)=resolved
    do m=1,3
       allocate(direct(m)%key(4,m*n),direct(m)%orientation(m*n),direct(m)%home(2,m*n))
       direct(m)%orientation=1
       do j=1,n
          do c=1,m
             direct(m)%key(:,m*(j-1)+c)=executor%placement%key(:,j)
             direct(m)%key(3,m*(j-1)+c)=c
             direct(m)%home(:,m*(j-1)+c)=[home(1,j),m*(home(2,j)-1)+c]
          end do
       end do
    end do
    do l=level_start,level_end
       ! The accepted node/edge alias transaction depends on level and topology,
       ! not on field phase or channel width. Resolve each exactly once here.
       do i=1,2
          m=1;location=AT_NODE
          if(i==2) then
             m=3;location=AT_EDGE
          end if
          call compile_rhs_boundary_map(location,l,base,mapping)
          if(allocated(metadata)) deallocate(metadata)
          allocate(metadata(size(mapping%key,2),5))
          metadata(:,:4)=real(transpose(mapping%key),dp);metadata(:,5)=1.0_dp
          do j=1,size(mapping%target)
             metadata(mapping%target(j),:4)=real(mapping%source(:,j),dp)
             metadata(mapping%target(j),5)=real(mapping%orientation(j),dp)
          end do
          if(present(directory)) then
             call deliver_rhs_metadata(mapping%key,metadata,direct(m)%key,received,topology,9900+10*l+i,directory,m)
          else
             call deliver_rhs_metadata(mapping%key,metadata,direct(m)%key,received,topology,9900+10*l+i)
          end if
          if(any(abs(real(nint(received),dp)-received)>0.0_dp)) error stop 'RHS alias metadata lost integer identity'
          boundary(i)%key=transpose(nint(received(:,:4)));boundary(i)%orientation=nint(received(:,5))
          boundary(i)%home=direct(m)%home
          do j=1,m*n
             boundary(i)%home(2,j)=m*((direct(m)%home(2,j)-1)/m)+boundary(i)%key(3,j)
          end do
          ! Unchanged interior identities already have a home. Only physical
          ! aliases outside that identity need an additional directory query.
          changed=pack([(j,j=1,m*n)],any(boundary(i)%key([1,2,4],:)/=direct(m)%key([1,2,4],:),dim=1))
          request=boundary(i)%key(:,changed);request(3,:)=1
          call lookup_rhs_owner_directory(owner_directory,request,resolved)
          do j=1,size(changed)
             a=changed(j)
             boundary(i)%home(:,a)=[resolved(1,j),m*(resolved(2,j)-1)+boundary(i)%key(3,a)]
          end do
       end do
       do i=1,size(phases)
          p=phases(i)
          ! Demand for these direct canonical fields is the union over all
          ! levels. The provider layout and widths are identical before/after
          ! injection and across levels; bind once within this preparation.
          rp=p;rl=l
          select case(p)
          case(RHS_BERNOULLI,RHS_INJECTED_BERNOULLI)
             rp=RHS_BERNOULLI;rl=level_start
          case(RHS_VELOCITY_SOURCE,RHS_DIFFUSION_DIRECT_NODES)
             rl=level_start
          end select
          t%route_phase(p,l)=rp;t%route_level(p,l)=rl
          if(rp/=p.or.rl/=l) then
             if(.not.t%bound(rp,rl)) error stop 'RHS shared route constructed out of order'
             logical_buffers=logical_buffers+rhs_exchange_buffer_values(t%route(rp,rl))
             cycle
          end if
          m=1;location=0
          select case(p)
          case(RHS_DIFFUSION_EDGES,RHS_ROTATION_EDGES,RHS_SCALAR_FLUX)
             m=3;location=2
          case(RHS_DIFFUSION_NODES,RHS_SCALAR_TENDENCY)
             location=1
          case(RHS_FINAL_PV,RHS_VELOCITY_SOURCE,RHS_DIFFUSION_DIRECT_EDGES,RHS_DIRECT_SCALAR_FLUX)
             m=3
          case(RHS_PRIMITIVE_VORT,RHS_ROTATION_TRIANGLES,RHS_ROTATION_DIRECT_TRIANGLES)
             m=2
          end select
          nc=1
          if(p==RHS_DIFFUSION_EDGES.or.p==RHS_DIFFUSION_NODES.or.p==RHS_SCALAR_FLUX.or. &
               p==RHS_SCALAR_TENDENCY.or.p==RHS_DIFFUSION_DIRECT_EDGES.or.p==RHS_DIRECT_SCALAR_FLUX.or. &
               p==RHS_BERNOULLI.or.p==RHS_INJECTED_BERNOULLI) nc=2
          ! Width does not change any producer or consumer identity. One route
          ! serves the one/two-field calls, with epoch stamps on every message.
          if(location==0) then
             call compile_rhs_phase_route(executor%placement%program,executor%placement%key,executor%placement%owned, &
                  p,l,nc,t%comm,topology,masks,direct(m)%key,direct(m)%orientation,t%route(p,l), &
                  shared_context=.true.,producer_home=direct(m)%home, &
                  consumer_home=direct(m)%home,shared_buffers=.true.)
          else
             call compile_rhs_phase_route(executor%placement%program,executor%placement%key,executor%placement%owned, &
                  p,l,nc,t%comm,topology,masks,boundary(location)%key,boundary(location)%orientation, &
                  t%route(p,l),shared_context=.true.,producer_home=direct(m)%home, &
                  consumer_home=boundary(location)%home,shared_buffers=.true.)
          end if
          t%bound(p,l)=.true.
          a=rhs_exchange_buffer_values(t%route(p,l));capacity=max(capacity,a)
          logical_buffers=logical_buffers+a
          if(rhs_phase_canonical_writers(p)) direct_routes=direct_routes+1
       end do
    end do
    if(prefix_demand_enabled(t%comm)) then
       call compile_prefix_demand(executor)
       capacity=0;logical_buffers=0
       do l=level_start,level_end
          do i=1,size(phases)
             p=phases(i);rp=t%route_phase(p,l);rl=t%route_level(p,l)
             a=rhs_exchange_buffer_values(t%route(rp,rl));capacity=max(capacity,a)
             logical_buffers=logical_buffers+a
          end do
       end do
    end if
    ! Every route completes all requests before another callback can run.
    ! One maximum-sized buffer replaces one payload allocation per route.
    allocate(t%scratch(capacity))
    if(audit316_enabled()) allocate(t%audit_event316(7,128*nlayer*(level_end-level_start+1)))
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',audit_option)
    if(trim(audit_option)=='1') write(*,'(a,6(i0,1x))') 'RHS281 plans: ',rank, &
         size(phases)*(level_end-level_start+1),count(t%bound),direct_routes,capacity,logical_buffers
    end associate
    call detail_leave(DP_RHS_ROUTES)
  end subroutine

  subroutine run_rhs_distributed(executor,scalar,mean,velocity,surface,order,nu,nu_pentagon,gravity,cp,p0,kappa,ptop)
    type(Rhs_Distributed_Executor),intent(inout) :: executor
    real(dp),intent(in) :: scalar(:,:,:),mean(:,:,:),velocity(:,:),surface(:),nu(:,:),nu_pentagon(:)
    real(dp),intent(in) :: gravity,cp,p0,kappa,ptop
    integer,intent(in) :: order(3)
    integer :: k,l,n,nlayer
    n=executor%placement%program%nnode;nlayer=size(scalar,3)
    if(any(shape(scalar)/=[n,2,nlayer]).or.any(shape(mean)/=shape(scalar)).or. &
         any(shape(velocity)/=[3*n,nlayer]).or.any(shape(nu)/=[4,nlayer]).or.size(nu_pentagon)/=nlayer.or. &
         size(surface)/=n.or.any(shape(executor%velocity_rhs)/=[3*n,nlayer])) error stop 'RHS distributed input extent'
    if(audit316_enabled()) then
       executor%transport%audit_used316=0
       call audit316_i('execution.shape',[n,nlayer])
       call audit316_i('execution.order',order)
       call audit316_r('execution.constants',[gravity,cp,p0,kappa,ptop])
       call audit316_i('execution.nu.shape',shape(nu))
       call audit316_r('execution.nu',reshape(nu,[size(nu)]))
       call audit316_r('execution.nu_pentagon',nu_pentagon)
       call audit316_i64('execution.state.begin',[executor%transport%state])
    end if
    call rhs_dynamics_start(executor%placement%program,scalar(:,1,:),mean(:,1,:),surface,gravity,ptop,executor%work)
    do k=1,nlayer
       if(audit316_enabled()) executor%transport%audit_layer316=k
       call rhs_dynamics_layer(executor%placement%program,k,scalar(:,:,k),mean(:,:,k),velocity(:,k:k),order,nu(:,k), &
            nu_pentagon(k),gravity,cp,p0,kappa,executor%work,rhs_distributed_halo,executor%transport)
       executor%scalar_rhs(:,:,k)=executor%work%scalar_rhs
       executor%velocity_rhs(:,k)=reshape(executor%work%velocity_rhs,[3*n])
    end do
    if(audit316_enabled()) then
       associate(t=>executor%transport)
       call audit316_i('execution.halo.shape',[7,t%audit_used316])
       call audit316_i('execution.halo.order',reshape(t%audit_event316(:,:t%audit_used316),[7*t%audit_used316]))
       call audit316_i64('execution.state.end',[t%state])
       end associate
    end if
    executor%columns=executor%columns+count(executor%placement%program%surface_visits>0)
    do l=level_start,level_end
       executor%primitive_events=executor%primitive_events+ &
            int(size(executor%placement%program%diffusion(l)%producer),int64)*nlayer
       executor%scalar_targets=executor%scalar_targets+ &
            int(size(executor%placement%program%dynamics(l)%divergence),int64)*nlayer
    end do
    executor%velocity_targets=executor%velocity_targets+int(size(executor%placement%program%gradient),int64)*nlayer
  end subroutine

  subroutine rhs_distributed_halo(phase,l,values,context)
    integer,intent(in) :: phase,l
    real(dp),intent(inout) :: values(:,:)
    class(*),intent(inout) :: context
    integer(int64) :: send,receive,bytes,halo_state
    integer :: c,sp,rp,p,level
    select type(t=>context)
    type is(Rhs_Route_Context)
       c=size(values,2)
       if(c<1.or.c>2.or.phase<1.or.phase>20.or.l<level_start.or.l>level_end) &
            error stop 'RHS distributed phase coordinates'
       p=t%route_phase(phase,l);level=t%route_level(phase,l)
       if(p<1.or.level<level_start.or.level>level_end) error stop 'RHS phase has no route binding'
       if(.not.t%bound(p,level)) error stop 'RHS distributed phase is not bound'
       t%state=t%state+1
       ! Shared physical routes must still reject mismatched logical phases.
       ! The transaction stride exceeds every level/phase identity, so stamps
       ! remain strictly increasing on every reused route.
       halo_state=t%state*(20_int64*(level_end-level_start+1)+1)+20_int64*(l-level_start)+phase
       if(audit316_enabled()) call audit_event(0,c)
       call exchange_rhs_values(t%route(p,level),t%topology,t%masks,100*level+p,halo_state,values, &
            scratch=t%scratch)
       call rhs_detail_enter(DP_RHS_COST)
       call rhs_exchange_cost(t%route(p,level),send,receive,sp,rp,bytes)
       call rhs_detail_leave(DP_RHS_COST)
       t%exchanges=t%exchanges+1;t%remote_values=t%remote_values+send;t%remote_messages=t%remote_messages+sp
       ! Preserve the existing native mass's second boundary application.
       ! It is distinct from temperature and must be retired by an exact test.
       if(phase==RHS_SCALAR_FLUX.and.c==2) then
          t%state=t%state+1
          halo_state=t%state*(20_int64*(level_end-level_start+1)+1)+20_int64*(l-level_start)+phase
          if(audit316_enabled()) call audit_event(1,1)
          call exchange_rhs_values(t%route(p,level),t%topology,t%masks,100*level+p,halo_state,values(:,:1), &
               scratch=t%scratch)
          call rhs_detail_enter(DP_RHS_COST)
          call rhs_exchange_cost(t%route(p,level),send,receive,sp,rp,bytes)
          call rhs_detail_leave(DP_RHS_COST)
          t%exchanges=t%exchanges+1;t%remote_values=t%remote_values+send;t%remote_messages=t%remote_messages+sp
       end if
    class default
       error stop 'RHS distributed transport context'
    end select
  contains
    subroutine audit_event(repeated,columns)
      integer,intent(in) :: repeated,columns
      select type(t=>context)
      type is(Rhs_Route_Context)
         if(.not.allocated(t%audit_event316)) error stop 'Stage316 missing event storage'
         t%audit_used316=t%audit_used316+1
         if(t%audit_used316>size(t%audit_event316,2)) error stop 'Stage316 event bound exceeded'
         t%audit_event316(:,t%audit_used316)=[phase,l,columns,p,level,repeated,t%audit_layer316]
      end select
    end subroutine
  end subroutine
end module
