module parallel_block_rhs_backend_mod
  use parallel_block_profile_mod, only : detail_enter,detail_leave, &
       DP_RHS_INPUT_PLAN,DP_RHS_RETAIN_PLAN,DP_RHS_INPUT,DP_RHS_KERNEL,DP_RHS_RETAIN,DP_RHS_RK
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use mpi_f08
  use arch_mod, only : rank,comm,glo_id,block_catalog
  use shared_mod, only : N_VARIABLE,zlevels,scalars,S_VELO,AT_NODE,AT_EDGE,level_start,level_end, &
       Laplace_sclr,Laplace_divu,Laplace_rotu,grav_accel,c_p,p_0,kappa,p_top
  use domain_mod, only : grid,Float_Field
  use parallel_block_rhs_compile_mod, only : Rhs_Domain_Operators
  use parallel_block_rhs_placement_mod, only : compile_rhs_owner_map
  use parallel_block_rhs_distribute_mod, only : deliver_rhs_metadata
  use parallel_block_rhs_distributed_mod
  use parallel_block_rhs_exchange_mod
  use parallel_block_rhs_state_mod
  use parallel_block_rhs_retention_mod
  implicit none
  private
  public :: connected_rhs_enabled,prepare_connected_rhs,evaluate_connected_rhs,retain_connected_rhs
  public :: prepare_connected_rhs_rk,install_connected_rhs_rk,connected_rhs_needs_prepare
  real(dp),save :: evaluate_seconds=0.0_dp,rk_seconds=0.0_dp
  type(Rhs_Distributed_Executor),save,target :: executor
  type(Rhs_Resident_Input),save :: input,rk_input
  type(Rhs_Retention),save :: retention
  type(Rhs_Exchange),save :: adapter_route
  real(dp),allocatable,save :: scalar(:,:,:),mean(:,:,:),velocity(:,:),surface(:),adapter_value(:,:), &
       rk_scalar(:,:,:),rk_mean(:,:,:),rk_velocity(:,:),rk_value(:,:)
  integer,allocatable,save :: adapter_domain(:),adapter_node(:),rk_domain(:),rk_node(:),owned(:)
  type :: Mask_State
     integer,allocatable :: node(:),edge(:)
  end type
  type(Mask_State),allocatable,save :: masks(:)
  integer(int64),save :: generation=-1,mask_epoch=0,state_epoch=0,rk_epoch=0
  integer(int64),save :: rhs_calls=0,rk_calls=0,adapter_rk_values=0,resident_rk_values=0
  logical,save :: option_read=.false.,enabled=.false.,completed=.false.,rk_prepared=.false.
contains
  logical function connected_rhs_enabled() result(value)
    character(len=8) :: option
    if(.not.option_read) then
       call get_environment_variable('WAVETRISK_CONNECTED_RHS',option)
       if(len_trim(option)>0.and.trim(option)/='0'.and.trim(option)/='1') error stop 'invalid WAVETRISK_CONNECTED_RHS'
       enabled=trim(option)=='1';option_read=.true.
    end if
    value=enabled
  end function

  logical function connected_rhs_needs_prepare(topology) result(rebuild)
    integer(int64),intent(in) :: topology
    logical :: local
    integer :: d,ierr
    local=generation/=topology
    if(.not.allocated(masks)) then
       local=.true.
    else if(size(masks)/=size(grid)) then
       local=.true.
    else if(.not.local) then
       do d=1,size(grid)
          if(size(masks(d)%node)/=grid(d)%node%length.or.size(masks(d)%edge)/=3*grid(d)%node%length) then
             local=.true.;exit
          end if
          if(any(masks(d)%node/=grid(d)%mask_n%elts(:grid(d)%node%length)).or. &
               any(masks(d)%edge/=grid(d)%mask_e%elts(:3*grid(d)%node%length))) local=.true.
       end do
    end if
    ! A single consensus also covers mask-only changes on geometry ranks with
    ! no local block. It is explicit preparation cost, not a numeric halo.
    call MPI_Allreduce(local,rebuild,1,MPI_LOGICAL,MPI_LOR,comm,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'RHS mask epoch consensus'
  end function

  subroutine prepare_connected_rhs(plan,base,topology,geometry_surface)
    type(Rhs_Domain_Operators),intent(in) :: plan(:)
    integer,intent(in) :: base(:)
    integer(int64),intent(in) :: topology
    real(dp),intent(in) :: geometry_surface(:)
    integer,allocatable :: physical(:,:),key(:,:),owner(:),definition(:,:),request(:,:),source(:), &
         target(:),order(:,:),orientation(:),rk_key(:,:),domains(:),nodes(:),adapter_index(:),rk_index(:)
    real(dp),allocatable :: delivered(:,:)
    real(dp) :: started
    logical,allocatable :: reader(:),is_adapter(:),is_reader(:),footprint(:)
    integer :: d,b,n,ng,nd,i,g,p,np,na,nr,a
    started=MPI_Wtime()
    call release_rhs_exchange(adapter_route)
    mask_epoch=mask_epoch+1;generation=topology
    call prepare_rhs_distributed(plan,base,topology,mask_epoch,zlevels,executor)
    call state_footprint(executor%placement%program,footprint)
    call detail_enter(DP_RHS_INPUT_PLAN)
    call prepare_rhs_resident_input(executor%placement%key,base,topology,mask_epoch,input,required=footprint)
    call detail_leave(DP_RHS_INPUT_PLAN)
    call detail_enter(DP_RHS_RETAIN_PLAN)
    call prepare_rhs_retention(executor,base,topology,retention)
    call detail_leave(DP_RHS_RETAIN_PLAN)
    if(allocated(masks)) deallocate(masks)
    allocate(masks(size(grid)))
    ng=base(size(base));nd=executor%placement%program%nnode
    allocate(physical(4,ng))
    do d=1,size(grid)
       b=base(d);n=grid(d)%node%length
       masks(d)%node=grid(d)%mask_n%elts(:n);masks(d)%edge=grid(d)%mask_e%elts(:3*n)
       do i=1,n
          physical(:,b+i)=[glo_id(rank+1,d)+1,i,1,1]
       end do
    end do
    call deliver_rhs_metadata(physical,reshape(geometry_surface,[ng,1]),executor%placement%key,delivered,topology,8501)
    if(allocated(scalar)) deallocate(scalar,mean,velocity,surface)
    allocate(scalar(nd,2,zlevels),mean(nd,2,zlevels),velocity(3*nd,zlevels))
    surface=delivered(:,1)
    ! Catalogue interiors are read only by the accepted halo pack lists and
    ! the remaining fixed-coarse velocity restriction. Uncatalogued slots are
    ! the explicit inverse/boundary adapter; they never feed interior RHS work
    ! when a resident catalogue identity exists.
    allocate(domains(ng),nodes(ng),is_adapter(ng),is_reader(ng))
    is_adapter=.false.;is_reader=.false.
    if(allocated(adapter_domain)) deallocate(adapter_domain,adapter_node,rk_domain,rk_node)
    do d=1,size(grid)
       b=base(d);n=grid(d)%node%length
       call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),owner,key)
       allocate(reader(n));reader=.false.
       do g=1,size(grid(d)%pack,2)
          do p=1,grid(d)%pack(AT_NODE,g)%length
             i=grid(d)%pack(AT_NODE,g)%elts(p)+1;reader(i)=.true.
          end do
          do p=1,grid(d)%pack(AT_EDGE,g)%length
             i=grid(d)%pack(AT_EDGE,g)%elts(p)/3+1;reader(i)=.true.
          end do
       end do
       ! Only level_start children are read by the fixed-coarse restriction;
       ! later tightening can enumerate its exact masked edge footprint.
       where(grid(d)%level%elts(:n)==level_start) reader=.true.
       a=grid(d)%patch%elts(3)%elts_start+1
       do i=1,n
          domains(b+i)=d;nodes(b+i)=i
          if(i<a) cycle
          is_adapter(b+i)=key(1,i)==0
          is_reader(b+i)=key(1,i)>0.and.reader(i)
       end do
       deallocate(reader)
    end do
    adapter_index=pack([(i,i=1,ng)],is_adapter.or.is_reader)
    rk_index=pack([(i,i=1,ng)],is_reader)
    request=physical(:,adapter_index);rk_key=physical(:,rk_index)
    adapter_domain=domains(adapter_index);adapter_node=nodes(adapter_index)
    rk_domain=domains(rk_index);rk_node=nodes(rk_index)
    owned=pack([(i,i=1,nd)],executor%placement%owned)
    np=size(owned);na=size(adapter_node);nr=size(rk_node)
    definition=executor%placement%key(:,owned)
    source=owned;target=[(nd+i,i=1,na)]
    allocate(order(2,np),orientation(na));order=1;orientation=1
    call compile_rhs_exchange(adapter_route,comm,topology,mask_epoch,8502,nd+na,5*zlevels, &
         definition,order,source,request,target,orientation)
    if(allocated(adapter_value)) deallocate(adapter_value,rk_value,rk_scalar,rk_mean,rk_velocity)
    allocate(adapter_value(nd+na,5*zlevels),rk_value(na,5*zlevels))
    allocate(rk_scalar(nr,2,zlevels),rk_mean(nr,0,zlevels),rk_velocity(3*nr,zlevels))
    adapter_value=0;rk_value=0
    call detail_enter(DP_RHS_INPUT_PLAN)
    call prepare_rhs_resident_input(rk_key,base,topology,mask_epoch,rk_input,apply_boundary=.false.,with_means=.false.)
    call detail_leave(DP_RHS_INPUT_PLAN)
    completed=.false.;rk_prepared=.false.
    call audit_plan(MPI_Wtime()-started)
  end subroutine

  subroutine state_footprint(plan,required)
    type(Rhs_Domain_Operators),intent(in) :: plan
    logical,allocatable,intent(out) :: required(:)
    integer :: l,j,n
    allocate(required(plan%nnode));required=plan%surface_visits>0
    do l=level_start,level_end
       call nodes(plan%pressure(l)%node)
       do j=1,size(plan%diffusion(l)%producer)
          associate(p=>plan%diffusion(l)%producer(j))
          call nodes(p%mass);call edges(p%velocity)
          if(p%southwest) then
             call nodes([p%sw_mass]);call edges(p%sw_velocity)
          end if
          end associate
       end do
       do j=1,size(plan%rotation(l)%pentagon)
          call nodes(plan%rotation(l)%pentagon(j)%mass)
          call edges(reshape(plan%rotation(l)%pentagon(j)%velocity,[21]))
       end do
    end do
    do j=1,size(plan%gradient)
       if(plan%gradient(j)%active) call nodes(plan%gradient(j)%node)
    end do
  contains
    subroutine nodes(address)
      integer,intent(in) :: address(:)
      integer :: i
      do i=1,size(address)
         n=address(i)
         if(n==0) cycle
         if(n<1.or.n>size(required)) error stop 'RHS state node footprint'
         required(n)=.true.
      end do
    end subroutine
    subroutine edges(address)
      integer,intent(in) :: address(:)
      integer :: i
      do i=1,size(address)
         if(address(i)==0) cycle
         n=(address(i)-1)/3+1
         if(n<1.or.n>size(required)) error stop 'RHS state vector footprint'
         required(n)=.true.
      end do
    end subroutine
  end subroutine

  subroutine evaluate_connected_rhs(q,nu,pentagon)
    type(Float_Field),intent(in) :: q(N_VARIABLE,zlevels)
    real(dp),intent(in) :: nu(:,:),pentagon(:)
    real(dp) :: started
    started=MPI_Wtime()
    if(rhs_calls==0.and.rank==0) write(*,'(a)') 'Stage278 connected resident RHS = on'
    state_epoch=state_epoch+1
    call detail_enter(DP_RHS_INPUT)
    call read_rhs_resident_input(input,q,state_epoch,scalar,mean,velocity)
    call detail_leave(DP_RHS_INPUT)
    call detail_enter(DP_RHS_KERNEL)
    call run_rhs_distributed(executor,scalar,mean,velocity,surface,[Laplace_sclr,Laplace_divu,Laplace_rotu], &
         nu,pentagon,grav_accel,c_p,p_0,kappa,p_top)
    call detail_leave(DP_RHS_KERNEL)
    evaluate_seconds=evaluate_seconds+MPI_Wtime()-started
    completed=.true.;rhs_calls=rhs_calls+1
  end subroutine

  subroutine retain_connected_rhs
    if(.not.completed) error stop 'RHS retention before evaluation'
    call detail_enter(DP_RHS_RETAIN)
    call retain_rhs_tendency(executor,retention)
    call detail_leave(DP_RHS_RETAIN)
  end subroutine

  subroutine prepare_connected_rhs_rk(base,h)
    type(Float_Field),intent(in) :: base(N_VARIABLE,zlevels)
    real(dp),intent(in) :: h
    integer :: n,i,d,j,k,v,c
    real(dp) :: started
    started=MPI_Wtime()
    if(.not.completed.or.rk_prepared) error stop 'RHS RK preparation sequence'
    call detail_enter(DP_RHS_RK)
    n=executor%placement%program%nnode;c=0
    do k=1,zlevels
       adapter_value(:n,c+1:c+2)=executor%scalar_rhs(:,:,k)
       adapter_value(:n,c+3:c+5)=transpose(reshape(executor%velocity_rhs(:,k),[3,n]));c=c+5
    end do
    call exchange_rhs_values(adapter_route,generation,mask_epoch,8502,state_epoch,adapter_value)
    do i=1,size(adapter_node)
       d=adapter_domain(i);j=adapter_node(i);c=0
       do k=1,zlevels
          do v=1,2
             rk_value(i,c+v)=base(scalars(1)+v-1,k)%data(d)%elts(j)+h*adapter_value(n+i,c+v)
          end do
          do v=1,3
             rk_value(i,c+2+v)=base(S_VELO,k)%data(d)%elts(3*(j-1)+v)+h*adapter_value(n+i,c+2+v)
          end do
          c=c+5
       end do
    end do
    ! RK has already committed the resident stage at this call site. Only the
    ! named catalogue readers are delivered. Only their temperature is installed;
    ! accepted mass/velocity aliases use raw RK above. All five delivered fields
    ! are included in the resident adapter payload counter (means excluded).
    rk_epoch=rk_epoch+1
    call read_rhs_resident_input(rk_input,base,rk_epoch,rk_scalar,rk_mean,rk_velocity)
    rk_seconds=rk_seconds+MPI_Wtime()-started
    rk_prepared=.true.
    call detail_leave(DP_RHS_RK)
  end subroutine

  subroutine audit_plan(seconds)
    ! Quantify the actual replacement layout/fanout only in validation runs.
    ! Payload memory excludes descriptors, compilation peaks, Domain geometry
    ! adapters, retained accepted workspaces, and MPI implementation storage.
    real(dp),intent(in) :: seconds
    character(len=8) :: option
    character(len=64) :: path
    integer :: unit,p,l,i,j,sp,rp
    integer(int64) :: sends,receives,bytes,route_bytes,values,rounding_bytes
    integer(int64) :: input_send,input_receive,input_messages,rk_send,rk_receive,rk_messages
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
    if(trim(option)/='1') return
    route_bytes=0;input_send=0;input_receive=0;input_messages=0
    call input_cost(input%scalar_route);call input_cost(input%vector_route);call input_cost(input%mean_route)
    rk_send=0;rk_receive=0;rk_messages=0
    call rk_cost(adapter_route);call rk_cost(rk_input%scalar_route);call rk_cost(rk_input%vector_route)
    do l=lbound(executor%transport%route,2),ubound(executor%transport%route,2)
       do p=1,size(executor%transport%route,1)
          if(.not.executor%transport%bound(p,l)) cycle
          call rhs_exchange_cost(executor%transport%route(p,l),sends,receives,sp,rp,bytes)
          route_bytes=route_bytes+bytes
       end do
    end do
    values=0
    values=values+size(scalar,kind=int64)
    values=values+size(mean,kind=int64)
    values=values+size(velocity,kind=int64)
    values=values+size(surface,kind=int64)
    values=values+size(adapter_value,kind=int64)
    values=values+size(rk_value,kind=int64)
    values=values+size(rk_scalar,kind=int64)
    values=values+size(rk_mean,kind=int64)
    values=values+size(rk_velocity,kind=int64)
    values=values+size(executor%scalar_rhs,kind=int64)
    values=values+size(executor%velocity_rhs,kind=int64)
    values=values+size(executor%work%lower,kind=int64)
    values=values+size(executor%work%phi,kind=int64)
    values=values+size(executor%work%pressure,kind=int64)
    values=values+size(executor%work%midpoint,kind=int64)
    values=values+size(executor%work%exner,kind=int64)
    values=values+size(executor%work%thermodynamic,kind=int64)
    values=values+size(executor%work%flux,kind=int64)
    values=values+size(executor%work%pv,kind=int64)
    values=values+size(executor%work%ke,kind=int64)
    values=values+size(executor%work%bernoulli,kind=int64)
    values=values+size(executor%work%vorticity,kind=int64)
    values=values+size(executor%work%diffusion_flux,kind=int64)
    values=values+size(executor%work%laplacian,kind=int64)
    values=values+size(executor%work%divu,kind=int64)
    values=values+size(executor%work%lap_divu,kind=int64)
    values=values+size(executor%work%rotation,kind=int64)
    values=values+size(executor%work%source,kind=int64)
    values=values+size(executor%work%physics,kind=int64)
    values=values+size(executor%work%scalar_rhs,kind=int64)
    values=values+size(executor%work%velocity_rhs,kind=int64)
    values=values+size(executor%work%scalar_physics,kind=int64)
    values=values+size(executor%work%sw_physics,kind=int64)
    values=values+size(input%scalar,kind=int64)
    values=values+size(input%mean,kind=int64)
    values=values+size(input%vector,kind=int64)
    values=values+size(input%read_scalar,kind=int64)
    values=values+size(input%read_mean,kind=int64)
    values=values+size(input%read_vector,kind=int64)
    values=values+size(rk_input%scalar,kind=int64)
    values=values+size(rk_input%vector,kind=int64)
    values=values+size(rk_input%read_scalar,kind=int64)
    values=values+size(rk_input%read_vector,kind=int64)
    rounding_bytes=0
    do i=1,size(retention%block)
       do j=1,3
          rounding_bytes=rounding_bytes+8_int64*(size(retention%block(i)%column(j)%exner,kind=int64)+ &
               size(retention%block(i)%column(j)%theta,kind=int64))+ &
               int(storage_size(.false.)/8,int64)*size(retention%block(i)%column(j)%ready,kind=int64)
       end do
    end do
    write(path,'(a,i0,a)') 'rhs278-cost-rank-',rank,'.txt'
    open(newunit=unit,file=trim(path),status='unknown',position='append',action='write')
    write(unit,*) generation,mask_epoch,size(owned),executor%placement%program%nnode,input%nprovider,rk_input%nprovider, &
         input_send,input_receive,input_messages,rk_send,rk_receive,rk_messages,route_bytes,8_int64*values, &
         rounding_bytes,executor%placement%sent_metadata_bytes,executor%placement%received_metadata_bytes,seconds
    close(unit)
  contains
    subroutine input_cost(route)
      type(Rhs_Exchange),intent(in) :: route
      call rhs_exchange_cost(route,sends,receives,sp,rp,bytes)
      input_send=input_send+sends;input_receive=input_receive+receives;input_messages=input_messages+sp
      route_bytes=route_bytes+bytes
    end subroutine
    subroutine rk_cost(route)
      type(Rhs_Exchange),intent(in) :: route
      call rhs_exchange_cost(route,sends,receives,sp,rp,bytes)
      rk_send=rk_send+sends;rk_receive=rk_receive+receives;rk_messages=rk_messages+sp
      route_bytes=route_bytes+bytes
    end subroutine
  end subroutine

  subroutine audit
    character(len=8) :: option
    character(len=64) :: path
    integer :: unit
    call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
    if(trim(option)/='1') return
    write(path,'(a,i0,a)') 'rhs278-production-rank-',rank,'.txt'
    open(newunit=unit,file=trim(path),status='unknown',position='append',action='write')
    write(unit,*) rhs_calls,rk_calls,generation,mask_epoch,size(grid),size(retention%block), &
         executor%primitive_events,executor%scalar_targets,executor%velocity_targets,retention%retained_stages, &
         input%adapter_values,input%resident_values,executor%transport%remote_values,executor%transport%remote_messages, &
         adapter_rk_values,resident_rk_values,evaluate_seconds,rk_seconds
    close(unit)
  end subroutine

  subroutine install_connected_rhs_rk(dest,validate)
    type(Float_Field),intent(inout) :: dest(N_VARIABLE,zlevels)
    logical,intent(in) :: validate
    integer :: i,d,j,k,v,c
    if(.not.rk_prepared) error stop 'RHS RK install before preparation'
    do i=1,size(adapter_node)
       d=adapter_domain(i);j=adapter_node(i);c=0
       do k=1,zlevels
          do v=1,2
             call install(dest(scalars(1)+v-1,k)%data(d)%elts(j),rk_value(i,c+v))
          end do
          do v=1,3
             call install(dest(S_VELO,k)%data(d)%elts(3*(j-1)+v),rk_value(i,c+2+v))
          end do
          c=c+5
       end do
    end do
    do i=1,size(rk_node)
       d=rk_domain(i);j=rk_node(i)
       do k=1,zlevels
          v=2
          call install(dest(scalars(1)+v-1,k)%data(d)%elts(j),rk_scalar(i,v,k))
       end do
    end do
    dest%bdry_uptodate=.false.
    adapter_rk_values=adapter_rk_values+int(size(adapter_node),int64)*5*zlevels
    resident_rk_values=resident_rk_values+int(size(rk_node),int64)*5*zlevels
    rk_prepared=.false.;completed=.false.;rk_calls=rk_calls+1
    call audit
  contains
    subroutine install(target,value)
      real(dp),intent(inout) :: target
      real(dp),intent(in) :: value
      if(validate) then
         if(transfer(target,0_int64)/=transfer(value,0_int64)) then
            write(*,*) 'RHS RK adapter mismatch',rank,d,j,k,v,target,value
            error stop 'RHS RK adapter differs bit-for-bit'
         end if
      end if
      target=value
    end subroutine
  end subroutine
end module
