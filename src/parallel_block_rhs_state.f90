module parallel_block_rhs_state_mod
  use audit316_mod, only : audit316_i,audit316_i64,audit316_l,audit316_enabled
  ! Resident interior inputs plus the explicitly bounded compatibility alias
  ! adapter. Geometric state is read only for uncatalogued boundary/scaffold
  ! sources requested by the compiled post-halo footprint.
  use iso_fortran_env, only : int64
  use kind_mod, only : dp
  use arch_mod, only : comm,rank,glo_id,block_catalog
  use shared_mod, only : N_VARIABLE,zlevels,scalars,S_VELO,AT_NODE,AT_EDGE,NONE
  use domain_mod, only : grid,Float_Field,sol_mean
  use patch_mod, only : PATCH_SIZE
  use parallel_block_mod, only : Block_Data,STORE_PATCH,apply_local_block_field_consumer, &
       compile_local_block_inverse_routes,gather_local_block_rhs_state
  use parallel_block_rhs_placement_mod, only : compile_rhs_owner_map
  use parallel_block_rhs_boundary_mod, only : Rhs_Boundary_Map,compile_rhs_boundary_map
  use parallel_block_rhs_distribute_mod, only : deliver_rhs_metadata,Rhs_Geometry_Directory
  use parallel_block_rhs_exchange_mod
  implicit none
  private
  public :: Rhs_Resident_Input,prepare_rhs_resident_input,read_rhs_resident_input,release_rhs_resident_input
  public :: audit_rhs_input316
  type :: Provider_Inventory
     integer :: pass=0,used=0
     integer,allocatable :: key(:,:)
  end type
  type :: Rhs_Resident_Input
     type(Rhs_Exchange) :: scalar_route,mean_route,vector_route
     integer(int64) :: topology=-1,masks=-1,adapter_values=0,resident_values=0
     integer :: nprovider=0,nslot=0
     integer,allocatable :: audit_provider_key325(:,:),audit_physical_key325(:,:)
     logical :: means=.true.
     integer,allocatable :: local_provider(:),address(:,:),adapter_provider(:),adapter_domain(:),adapter_node(:)
     real(dp),allocatable :: scalar(:,:),mean(:,:),vector(:,:),read_scalar(:,:),read_mean(:,:),read_vector(:,:)
  end type
contains
  subroutine audit_rhs_input316(input)
    type(Rhs_Resident_Input),intent(in) :: input
    call audit316_i64('input.epoch',[input%topology,input%masks])
    call audit316_i('input.extent',[input%nprovider,input%nslot])
    call audit316_i('input.provider_key325.shape',shape(input%audit_provider_key325))
    call audit316_i('input.provider_key325',reshape(input%audit_provider_key325,[size(input%audit_provider_key325)]))
    call audit316_i('input.physical_key325.shape',shape(input%audit_physical_key325))
    call audit316_i('input.physical_key325',reshape(input%audit_physical_key325,[size(input%audit_physical_key325)]))
    call audit316_l('input.means',[input%means])
    call audit316_i('input.local_provider',input%local_provider)
    call audit316_i('input.address.shape',shape(input%address))
    call audit316_i('input.address',reshape(input%address,[size(input%address)]))
    call audit316_i('input.adapter_provider',input%adapter_provider)
    call audit316_i('input.adapter_domain',input%adapter_domain)
    call audit316_i('input.adapter_node',input%adapter_node)
    call audit_rhs_exchange316(input%scalar_route,'input.scalar_route')
    call audit_rhs_exchange316(input%vector_route,'input.vector_route')
    call audit_rhs_exchange316(input%mean_route,'input.mean_route')
    call audit316_i('input.scalar.shape',shape(input%scalar))
    call audit316_i('input.mean.shape',shape(input%mean))
    call audit316_i('input.vector.shape',shape(input%vector))
    call audit316_i('input.read_scalar.shape',shape(input%read_scalar))
    call audit316_i('input.read_mean.shape',shape(input%read_mean))
    call audit316_i('input.read_vector.shape',shape(input%read_vector))
  end subroutine

  subroutine collect_providers(catalog,block,context)
    integer,intent(in) :: catalog
    type(Block_Data),intent(in) :: block
    class(*),intent(inout) :: context
    integer :: p,i
    select type(c=>context)
    type is(Provider_Inventory)
       do p=1,size(block%patch)
          do i=1,PATCH_SIZE**2
             c%used=c%used+1
             if(c%pass==2) c%key(:,c%used)=[catalog,STORE_PATCH,p,i]
          end do
       end do
    class default
       error stop 'RHS resident inventory context'
    end select
  end subroutine

  subroutine release_rhs_resident_input(input)
    type(Rhs_Resident_Input),intent(inout) :: input
    call release_rhs_exchange(input%scalar_route)
    call release_rhs_exchange(input%mean_route)
    call release_rhs_exchange(input%vector_route)
    input=Rhs_Resident_Input()
  end subroutine

  subroutine prepare_rhs_resident_input(physical_key,base,topology,masks,input,apply_boundary,required,with_means,directory)
    integer,intent(in) :: physical_key(:,:),base(:)
    integer(int64),intent(in) :: topology,masks
    type(Rhs_Resident_Input),intent(inout) :: input
    logical,intent(in),optional :: apply_boundary,with_means
    logical,intent(in),optional :: required(:)
    type(Rhs_Geometry_Directory),optional,intent(in) :: directory
    logical :: with_halo
    type(Provider_Inventory) :: inventory
    type(Rhs_Boundary_Map) :: mapping
    integer,allocatable :: geo_key(:,:),geo_state(:,:),owners(:),keys(:,:),provider_key(:,:), &
         provider_domain(:),provider_node(:),request(:,:),query(:,:),sign(:),definition(:,:), &
         order(:,:),source(:),target(:),demand(:),index(:),address_key(:,:),needed(:)
    real(dp),allocatable :: metadata(:,:),resolved(:,:),aliases(:,:)
    logical,allocatable :: active(:)
    integer :: d,n,i,j,b,p,m,mode,ng,nd,np,nlocal,na,c,nc,pass,nq,t
    if(scalars(2)-scalars(1)+1/=2) error stop 'RHS resident path requires mass and temperature'
    call release_rhs_resident_input(input)
    with_halo=.true.
    if(present(apply_boundary)) with_halo=apply_boundary
    if(present(with_means)) input%means=with_means
    ng=base(size(base));nd=size(physical_key,2)
    needed=[(i,i=1,nd)]
    if(present(required)) then
       if(size(required)/=nd) error stop 'RHS state footprint extent'
       needed=pack(needed,required)
    end if
    nq=size(needed)
    allocate(geo_key(4,ng),geo_state(4,ng))
    do d=1,size(grid)
       b=base(d);n=grid(d)%node%length
       call compile_rhs_owner_map(grid(d),glo_id(rank+1,d),owners,keys)
       do i=1,n
          geo_key(:,b+i)=[glo_id(rank+1,d)+1,i,1,1]
          if(keys(1,i)>0) then
             geo_state(:,b+i)=keys(:,i);geo_state(4,b+i)=keys(4,i)+1
          else
             geo_state(:,b+i)=[size(block_catalog)+glo_id(rank+1,d)+1,4,i,1]
          end if
       end do
    end do
    do pass=1,2
       inventory%pass=pass;inventory%used=0
       call apply_local_block_field_consumer(collect_providers,inventory)
       if(inventory%pass==1) allocate(inventory%key(4,inventory%used))
    end do
    nlocal=inventory%used;na=count(geo_state(1,:)>size(block_catalog));np=nlocal+na
    allocate(provider_key(4,np),provider_domain(np),provider_node(np),active(np))
    provider_key(:,:nlocal)=inventory%key;provider_domain=0;provider_node=0;active=.false.
    p=nlocal
    do d=1,size(grid)
       b=base(d);n=grid(d)%node%length
       do i=1,n
          if(geo_state(1,b+i)<=size(block_catalog)) cycle
          p=p+1;provider_key(:,p)=geo_state(:,b+i);provider_domain(p)=d;provider_node(p)=i
       end do
    end do
    input%nprovider=np;input%nslot=nd;input%topology=topology;input%masks=masks
    if(audit316_enabled()) then
       input%audit_provider_key325=provider_key;input%audit_physical_key325=physical_key
    end if
    do mode=1,3
       if(mode==3.and..not.input%means) cycle
       ! Scalar perturbation and velocity use the exact accepted NONE halo.
       ! Means retain their own stored epoch and must not silently share it.
       m=1;nc=2*zlevels
       if(mode==2) then
          m=3;nc=zlevels
       end if
       if(allocated(query)) deallocate(query,request,sign,definition,source,target,order)
       allocate(query(4,m*nq),request(4,m*nq),sign(m*nq),definition(4,m*np), &
            source(m*np),target(m*nq),order(2,m*np))
       do t=1,nq
          i=needed(t)
          do c=1,m
             query(:,m*(t-1)+c)=physical_key(:,i);query(3,m*(t-1)+c)=c
             target(m*(t-1)+c)=m*np+m*(i-1)+c
          end do
       end do
       sign=1
       if(mode/=3.and.with_halo) then
          call compile_rhs_boundary_map(merge(AT_EDGE,AT_NODE,m==3),NONE,base,mapping)
          if(allocated(metadata)) deallocate(metadata)
          allocate(metadata(size(mapping%key,2),5))
          metadata(:,:4)=real(transpose(mapping%key),dp);metadata(:,5)=1.0_dp
          do i=1,size(mapping%target)
             metadata(mapping%target(i),:4)=real(mapping%source(:,i),dp)
             metadata(mapping%target(i),5)=real(mapping%orientation(i),dp)
          end do
          if(present(directory)) then
             call deliver_rhs_metadata(mapping%key,metadata,query,aliases,topology,8100+mode,directory,m)
          else
             call deliver_rhs_metadata(mapping%key,metadata,query,aliases,topology,8100+mode)
          end if
          if(any(abs(real(nint(aliases),dp)-aliases)>0.0_dp)) error stop 'RHS state alias identity'
          query=transpose(nint(aliases(:,:4)));sign=nint(aliases(:,5))
       end if
       ! Resolve every physical source to either a resident catalogue node or
       ! an explicitly uncatalogued compatibility slot on its geometry host.
       index=query(3,:);query(3,:)=1
       if(present(directory)) then
          call deliver_rhs_metadata(geo_key,real(transpose(geo_state),dp),query,resolved,topology,8200+mode,directory,1)
       else
          call deliver_rhs_metadata(geo_key,real(transpose(geo_state),dp),query,resolved,topology,8200+mode)
       end if
       if(any(abs(real(nint(resolved),dp)-resolved)>0.0_dp)) error stop 'RHS state provider identity'
       request=transpose(nint(resolved))
       request(4,:)=m*(request(4,:)-1)+index
       do i=1,np
          do c=1,m
             j=m*(i-1)+c;definition(:,j)=provider_key(:,i)
             definition(4,j)=m*(provider_key(4,i)-1)+c
          end do
       end do
       source=[(i,i=1,m*np)];order=1
       select case(mode)
       case(1)
          call compile_rhs_exchange(input%scalar_route,comm,topology,masks,8301,m*(np+nd),nc, &
               definition,order,source,request,target,sign)
          call rhs_exchange_source_slots(input%scalar_route,demand)
          allocate(input%scalar(np+nd,nc));input%scalar=0
       case(2)
          call compile_rhs_exchange(input%vector_route,comm,topology,masks,8302,m*(np+nd),nc, &
               definition,order,source,request,target,sign)
          call rhs_exchange_source_slots(input%vector_route,demand)
          allocate(input%vector(3*(np+nd),nc));input%vector=0
       case(3)
          call compile_rhs_exchange(input%mean_route,comm,topology,masks,8303,m*(np+nd),nc, &
               definition,order,source,request,target,sign)
          call rhs_exchange_source_slots(input%mean_route,demand)
          allocate(input%mean(np+nd,nc));input%mean=0
       end select
       if(any(demand<1).or.any(demand>m*np)) error stop 'RHS state provider demand outside registered sources'
       do i=1,size(demand)
          active((demand(i)-1)/m+1)=.true.
       end do
    end do
    input%local_provider=pack([(i,i=1,nlocal)],active(:nlocal))
    address_key=provider_key(:,input%local_provider);address_key(4,:)=address_key(4,:)-1
    call compile_local_block_inverse_routes(address_key,input%address)
    input%adapter_provider=pack([(i,i=nlocal+1,np)],active(nlocal+1:))
    input%adapter_domain=provider_domain(input%adapter_provider)
    input%adapter_node=provider_node(input%adapter_provider)
    n=size(input%local_provider)
    allocate(input%read_scalar(n,2*zlevels),input%read_mean(n,merge(2*zlevels,0,input%means)), &
         input%read_vector(n,3*zlevels))
    if(.not.input%means) allocate(input%mean(np+nd,0))
  end subroutine

  subroutine read_rhs_resident_input(input,boundary_stage,state,scalar,mean,velocity)
    type(Rhs_Resident_Input),intent(inout) :: input
    type(Float_Field),intent(in) :: boundary_stage(N_VARIABLE,zlevels)
    integer(int64),intent(in) :: state
    real(dp),intent(out) :: scalar(:,:,:),mean(:,:,:),velocity(:,:)
    integer :: i,p,k,e,d,n,np,nd
    np=input%nprovider;nd=input%nslot
    if(any(shape(scalar)/=[nd,2,zlevels]).or.any(shape(mean)/=[nd,merge(2,0,input%means),zlevels]).or. &
         any(shape(velocity)/=[3*nd,zlevels])) error stop 'RHS resident output extent'
    call gather_local_block_rhs_state(input%address,1,input%read_scalar)
    call gather_local_block_rhs_state(input%address,2,input%read_vector)
    if(input%means) call gather_local_block_rhs_state(input%address,3,input%read_mean)
    do i=1,size(input%local_provider)
       p=input%local_provider(i)
       input%scalar(p,:)=input%read_scalar(i,:)
       if(input%means) input%mean(p,:)=input%read_mean(i,:)
       do k=1,zlevels
          do e=1,3
             input%vector(3*(p-1)+e,k)=input%read_vector(i,3*(k-1)+e)
          end do
       end do
    end do
    do i=1,size(input%adapter_provider)
       p=input%adapter_provider(i);d=input%adapter_domain(i);n=input%adapter_node(i)
       do k=1,zlevels
          do e=1,2
             input%scalar(p,2*(k-1)+e)=boundary_stage(scalars(1)+e-1,k)%data(d)%elts(n)
             if(input%means) input%mean(p,2*(k-1)+e)=sol_mean(scalars(1)+e-1,k)%data(d)%elts(n)
          end do
          do e=1,3
             input%vector(3*(p-1)+e,k)=boundary_stage(S_VELO,k)%data(d)%elts(3*(n-1)+e)
          end do
       end do
    end do
    input%adapter_values=input%adapter_values+int(size(input%adapter_provider),int64)*merge(7,5,input%means)*zlevels
    input%resident_values=input%resident_values+int(size(input%local_provider),int64)*merge(7,5,input%means)*zlevels
    call exchange_rhs_values(input%scalar_route,input%topology,input%masks,8301,state,input%scalar)
    call exchange_rhs_values(input%vector_route,input%topology,input%masks,8302,state,input%vector)
    if(input%means) call exchange_rhs_values(input%mean_route,input%topology,input%masks,8303,state,input%mean)
    do k=1,zlevels
       scalar(:,:,k)=input%scalar(np+1:,2*k-1:2*k)
       if(input%means) mean(:,:,k)=input%mean(np+1:,2*k-1:2*k)
       velocity(:,k)=input%vector(3*np+1:,k)
    end do
  end subroutine
end module
