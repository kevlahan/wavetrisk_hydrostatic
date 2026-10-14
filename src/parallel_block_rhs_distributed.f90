module parallel_block_rhs_distributed_mod
  use parallel_block_profile_mod, only : detail_enter,detail_leave,rhs_detail_enter,rhs_detail_leave, &
       DP_RHS_GEOMETRY,DP_RHS_ROUTES,DP_RHS_COST
  ! Connected direct execution on target block owners. The state adapter must
  ! supply this compact input layout at the declared epoch before run entry.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  use arch_mod, only : comm,rank
  use shared_mod, only : AT_NODE,AT_EDGE,level_start,level_end
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  use parallel_block_rhs_distribute_mod
  use parallel_block_rhs_boundary_mod
  use parallel_block_rhs_routes_mod
  use parallel_block_rhs_exchange_mod
  use parallel_block_rhs_diffusion_mod
  use parallel_block_rhs_dynamics_mod
  implicit none
  private
  public :: Rhs_Distributed_Executor,prepare_rhs_distributed,run_rhs_distributed,rhs_distributed_halo,release_rhs_distributed
  type :: Rhs_Route_Context
     type(MPI_Comm) :: comm=MPI_COMM_NULL
     type(Rhs_Exchange),allocatable :: route(:,:)
     logical,allocatable :: bound(:,:)
     integer(int64) :: topology=-1,masks=-1,state=0,exchanges=0,remote_values=0,remote_messages=0
  end type
  type :: Rhs_Distributed_Executor
     type(Rhs_Placement) :: placement
     type(Rhs_Dynamics_Workspace) :: work
     type(Rhs_Route_Context) :: transport
     real(dp),allocatable :: scalar_rhs(:,:,:),velocity_rhs(:,:)
     integer(int64) :: columns=0,primitive_events=0,scalar_targets=0,velocity_targets=0
  end type
contains
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
    integer :: p,l,c,m,i,j,n,location,nc,ierr,a
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
    allocate(t%route(20,level_start:level_end),t%bound(20,level_start:level_end))
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
          p=phases(i);m=1;location=0
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
                  shared_context=.true.,producer_home=direct(m)%home,consumer_home=direct(m)%home)
          else
             call compile_rhs_phase_route(executor%placement%program,executor%placement%key,executor%placement%owned, &
                  p,l,nc,t%comm,topology,masks,boundary(location)%key,boundary(location)%orientation, &
                  t%route(p,l),shared_context=.true.,producer_home=direct(m)%home,consumer_home=boundary(location)%home)
          end if
          t%bound(p,l)=.true.
       end do
    end do
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
    call rhs_dynamics_start(executor%placement%program,scalar(:,1,:),mean(:,1,:),surface,gravity,ptop,executor%work)
    do k=1,nlayer
       call rhs_dynamics_layer(executor%placement%program,k,scalar(:,:,k),mean(:,:,k),velocity(:,k:k),order,nu(:,k), &
            nu_pentagon(k),gravity,cp,p0,kappa,executor%work,rhs_distributed_halo,executor%transport)
       executor%scalar_rhs(:,:,k)=executor%work%scalar_rhs
       executor%velocity_rhs(:,k)=reshape(executor%work%velocity_rhs,[3*n])
    end do
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
    integer(int64) :: send,receive,bytes
    integer :: c,sp,rp
    select type(t=>context)
    type is(Rhs_Route_Context)
       c=size(values,2)
       if(c<1.or.c>2.or.phase<1.or.phase>20.or.l<level_start.or.l>level_end) &
            error stop 'RHS distributed phase coordinates'
       if(.not.t%bound(phase,l)) error stop 'RHS distributed phase is not bound'
       t%state=t%state+1
       call exchange_rhs_values(t%route(phase,l),t%topology,t%masks,100*l+phase,t%state,values)
       call rhs_detail_enter(DP_RHS_COST)
       call rhs_exchange_cost(t%route(phase,l),send,receive,sp,rp,bytes)
       call rhs_detail_leave(DP_RHS_COST)
       t%exchanges=t%exchanges+1;t%remote_values=t%remote_values+send;t%remote_messages=t%remote_messages+sp
       ! Preserve the existing native mass's second boundary application.
       ! It is distinct from temperature and must be retired by an exact test.
       if(phase==RHS_SCALAR_FLUX.and.c==2) then
          t%state=t%state+1
          call exchange_rhs_values(t%route(phase,l),t%topology,t%masks,100*l+phase,t%state,values(:,:1))
          call rhs_detail_enter(DP_RHS_COST)
          call rhs_exchange_cost(t%route(phase,l),send,receive,sp,rp,bytes)
          call rhs_detail_leave(DP_RHS_COST)
          t%exchanges=t%exchanges+1;t%remote_values=t%remote_values+send;t%remote_messages=t%remote_messages+sp
       end if
    class default
       error stop 'RHS distributed transport context'
    end select
  end subroutine
end module
