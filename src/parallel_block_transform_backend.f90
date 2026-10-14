module parallel_block_transform_backend_mod
  use iso_fortran_env,only : int64
  use kind_mod,only : dp
  use parallel_block_profile_mod,only : detail_enter,detail_leave,DP_TRANSFORM_PLAN,DP_TRANSFORM_INPUT, &
       DP_TRANSFORM_FORWARD,DP_TRANSFORM_INVERSE,DP_TRANSFORM_OUTPUT,DP_TRANSFORM_EXCHANGE
  use mpi_f08
  use arch_mod,only : comm,rank,n_process,glo_id,loc_id,block_catalog,geometry_owner=>owner
  use shared_mod,only : EDGE,N_VARIABLE,zlevels,scalars,S_VELO,AT_NODE,AT_EDGE,NONE,ADJZONE,level_start,level_end
  use domain_mod,only : grid,sol,Float_Field
  use patch_mod,only : PATCH_SIZE
  use parallel_block_transform_kernels_mod
  use parallel_block_transform_adapter_mod
  use parallel_block_transform_placement_mod
  use parallel_block_rhs_placement_mod,only : compile_rhs_owner_map
  use parallel_block_rhs_distribute_mod,only : Rhs_Geometry_Directory,prepare_rhs_geometry_directory,deliver_rhs_metadata
  use parallel_block_rhs_exchange_mod
  use parallel_block_rhs_boundary_mod
  use parallel_block_rhs_backend_mod,only : connected_rhs_mask_generation
  use parallel_block_mod,only : Block_Data,STORE_PATCH,STORE_BDRY,BLOCK_PAYLOAD_SOL,BLOCK_PAYLOAD_WAV_COEFF, &
       compile_local_block_inverse_routes,transfer_local_block_inverse_routes, &
       apply_local_block_field_consumer
  implicit none
  private
  public :: check_transform_vector_stage
  public :: check_resident_transform, resident_transform_check_enabled, install_checked_transform
  public :: resident_transform_enabled,prepare_resident_transform,execute_resident_transform,release_backend
  type :: Boundary_Identity
     integer,allocatable :: key(:,:),sign(:)
  end type
  type :: Transform_Route
     type(Rhs_Exchange) :: exchange
     integer :: phase=0
  end type
  type :: Transform_Homes
     integer,allocatable :: key(:,:),home(:,:),table(:)
     integer :: used=0
  end type
  type :: Boundary_Inventory
     integer :: pass=0,used=0
     integer,allocatable :: physical(:,:),state_key(:,:)
  end type
  type(Transform_Placement),save :: placement
  type(MPI_Comm),save :: context=MPI_COMM_NULL
  type(Transform_Route),allocatable,save :: scalar_halo(:),vector_halo(:),outer_halo(:)
  type(Transform_Adapter),save :: scalar_adapter,vector_adapter
  type(Transform_Route),save :: scalar_initial,vector_initial
  type(Transform_Route),save :: scalar_wave_halo,vector_wave_halo,scalar_publication,vector_publication
  type(Rhs_Exchange),save :: fixed_scalar_route,fixed_vector_route
  type(Boundary_Identity),allocatable,save :: scalar_boundary(:),vector_boundary(:)
  real(dp),allocatable,save :: scalar(:,:),vector(:,:),scalar_wavelet(:,:),vector_wavelet(:,:), &
       fixed_scalar(:,:),fixed_vector(:,:),wire(:),exchange_scratch(:)
  integer,allocatable,save :: info(:,:),save_node(:),save_index(:),fixed_node(:),fixed_index(:), &
       interior_node(:),interior_address(:,:),boundary_address(:,:),publication_domain(:),publication_node(:),publication_row(:), &
       forward_inner(:),base(:),fixed_domain(:),fixed_physical_node(:)
  integer(int64),save :: generation=-1,masks=0,state=0,completed=0
  real(dp),allocatable,save :: checked_forward_vector(:,:),checked_wavelet_vector(:,:)
  integer,save :: nnode=0,nrow=0,nfixed=0,npublication=0,nboundary=0
  logical,save :: option_read=.false.,enabled=.false.,audit_on=.false.
  integer(int64),save :: route_bytes=0,buffer_bytes=0,sent_values=0,received_values=0,messages=0,rounds=0
  integer(int64),save :: builds=0
  real(dp),save :: preparation_seconds=0

contains
  logical function resident_transform_enabled() result(value)
    character(8) :: option
    if(.not.option_read) then
       call get_environment_variable('WAVETRISK_RESIDENT_TRANSFORM',option)
       if(len_trim(option)>0.and.trim(option)/='0'.and.trim(option)/='1') error stop 'invalid resident transform selector'
       enabled=trim(option)=='1';option_read=.true.
       call get_environment_variable('WAVETRISK_CONNECTED_RHS_AUDIT',option)
       audit_on=trim(option)=='1'
       if(enabled.and.rank==0) write(*,'(a)') 'Stage279 resident provisional transforms = on'
    end if
    value=enabled
  end function

  logical function resident_transform_check_enabled() result(value)
    character(8) :: option
    call get_environment_variable('WAVETRISK_RESIDENT_TRANSFORM_CHECK',option)
    if(len_trim(option)>0.and.trim(option)/='0'.and.trim(option)/='1') error stop 'invalid transform check selector'
    value=trim(option)=='1'
  end function

  subroutine check_transform_vector_stage(values,displacement,ncolumns,first_column,coefficients)
    real(dp),intent(in) :: values(:)
    integer,intent(in) :: displacement(:),ncolumns,first_column
    logical,intent(in) :: coefficients
    integer,allocatable :: physical(:,:),request(:,:)
    real(dp),allocatable :: source(:,:),reference(:,:)
    logical,allocatable :: needed(:)
    type(Rhs_Geometry_Directory) :: directory
    integer :: d,i,j,k,e,a,b,n,stride,inputs(18)
    real(dp) :: expected,observed
    if(.not.resident_transform_check_enabled()) return
    if(.not.allocated(checked_forward_vector)) return
    n=base(size(base));allocate(physical(4,3*n),source(3*n,zlevels),request(4,3*nnode),needed(3*nnode))
    do d=1,size(grid)
       stride=(displacement(d+1)-displacement(d))/ncolumns
       do i=1,grid(d)%node%length
          do e=1,3
             a=3*(base(d)+i-1)+e;physical(:,a)=[glo_id(rank+1,d)+1,i,e,1]
             do k=1,zlevels
                b=displacement(d)+(k-first_column)*stride+3*(i-1)+e
                source(a,k)=values(b)
             end do
          end do
       end do
    end do
    do i=1,nnode
       do e=1,3
          request(:,3*(i-1)+e)=[placement%key(1:2,i),e,1]
       end do
    end do
    call prepare_rhs_geometry_directory(base,generation,directory)
    call deliver_rhs_metadata(physical,source,request,reference,generation,merge(8952,8951,coefficients),directory,3)
    needed=.false.
    if(coefficients) then
       do i=1,nnode
          if(placement%owned(i).and.info(1,i)==1.and.info(2,i)>level_start) needed(3*i-2:3*i)=.true.
       end do
    else
       inputs=[1,2,3,5,6,7,10,11,14,15,18,19,22,23,26,27,30,31]
       do i=1,size(placement%program%outer)
          associate(op=>placement%program%outer(i))
          if(op%level<level_start.or.op%forward_active==0) cycle
          do j=1,13
             if(op%source(j)>0) needed(op%source(j))=.true.
          end do
          needed(op%target(2:3))=.true.
          end associate
       end do
       do i=1,size(placement%program%inner)
          associate(op=>placement%program%inner(i))
          needed(op%source(inputs))=.true.;needed(op%target)=.true.
          end associate
       end do
    end if
    do i=1,3*nnode
       if(.not.needed(i)) cycle
       do k=1,zlevels
          observed=checked_forward_vector(i,k)
          if(coefficients) observed=checked_wavelet_vector(i,k)
          expected=reference(i,k)
          if(transfer(expected,0_int64)/=transfer(observed,0_int64)) then
             write(*,*) 'Resident forward stage difference ',coefficients,rank,request(:,i),k,expected,observed
             error stop 'resident forward input or coefficient differs'
          end if
       end do
    end do
  end subroutine

  subroutine check_resident_transform(scaling)
    type(Float_Field),intent(in) :: scaling(N_VARIABLE,zlevels)
    integer :: i,p,q,k,v,d,np,component,j
    real(dp) :: expected,observed
    do component=1,2
       call transfer_local_block_inverse_routes(interior_address,BLOCK_PAYLOAD_SOL,component,.false.,wire)
       do i=1,size(interior_node)
          p=interior_node(i);q=p
          if(save_index(p)>0) q=save_index(p)
          do k=1,zlevels
             do v=1,merge(2,3,component==1)
                if(component==1) then
                   expected=scalar(q,2*(k-1)+v);observed=wire((i-1)*2*zlevels+2*(k-1)+v)
                else
                   expected=vector(3*(q-1)+v,k);observed=wire((i-1)*3*zlevels+3*(k-1)+v)
                end if
                if(transfer(expected,0_int64)/=transfer(observed,0_int64)) then
                   write(*,*) 'Resident transform canonical difference ',rank,component,placement%key(:,p),k,v, &
                        expected,observed
                   error stop 'resident transform differs from independent inverse'
                end if
             end do
          end do
       end do
    end do
    np=size(publication_node)
    do i=1,np
       d=publication_domain(i);p=publication_node(i);q=publication_row(i)
       do k=1,zlevels
          do v=1,5
             if(v<=2) then
                expected=scalar(q,2*(k-1)+v);j=scalars(1)+v-1;observed=scaling(j,k)%data(d)%elts(p)
             else
                expected=vector(3*(q-1)+v-2,k);observed=scaling(S_VELO,k)%data(d)%elts(3*(p-1)+v-2)
             end if
             if(transfer(expected,0_int64)/=transfer(observed,0_int64)) then
                write(*,*) 'Resident transform adapter difference ',rank,d,p,k,v,expected,observed
                error stop 'resident transform adapter differs'
             end if
          end do
       end do
    end do
    do component=1,2
       call transfer_local_block_inverse_routes(boundary_address,BLOCK_PAYLOAD_SOL,component,.false.,wire)
       do i=1,nboundary
          q=publication_row(np+i)
          do k=1,zlevels
             do v=1,merge(2,3,component==1)
                if(component==1) then
                   expected=scalar(q,2*k-2+v);observed=wire((i-1)*2*zlevels+2*k-2+v)
                else
                   expected=vector(3*(q-1)+v,k);observed=wire((i-1)*3*zlevels+3*k-3+v)
                end if
                if(transfer(expected,0_int64)/=transfer(observed,0_int64)) &
                     error stop 'resident transform block boundary differs'
             end do
          end do
       end do
    end do
    if(rank==0) write(*,'(a,i0)') 'Resident transform shadow passed: ',completed
  end subroutine

  subroutine release_route(route)
    type(Transform_Route),intent(inout) :: route
    call release_rhs_exchange(route%exchange)
    route%phase=0
  end subroutine

  subroutine release_backend
    integer :: i,ierr
    call release_rhs_exchange(fixed_scalar_route);call release_rhs_exchange(fixed_vector_route)
    call release_route(scalar_initial);call release_route(vector_initial)
    scalar_adapter=Transform_Adapter();vector_adapter=Transform_Adapter()
    call release_route(scalar_wave_halo);call release_route(vector_wave_halo)
    call release_route(scalar_publication);call release_route(vector_publication)
    if(allocated(scalar_halo)) then
       do i=lbound(scalar_halo,1),ubound(scalar_halo,1)
          call release_route(scalar_halo(i));call release_route(vector_halo(i))
          call release_route(outer_halo(i))
       end do
       deallocate(scalar_halo,vector_halo,outer_halo,scalar_boundary,vector_boundary)
    end if
    if(context/=MPI_COMM_NULL) then
       call MPI_Comm_free(context,ierr)
       if(ierr/=MPI_SUCCESS) error stop 'transform context release'
    end if
    if(allocated(scalar)) deallocate(scalar,vector,scalar_wavelet,vector_wavelet, &
         fixed_scalar,fixed_vector,wire,exchange_scratch,info,save_node,save_index,fixed_node,fixed_index, &
         interior_node,interior_address, &
         boundary_address,publication_domain,publication_node,publication_row,forward_inner,base,fixed_domain,fixed_physical_node)
    if(allocated(checked_forward_vector)) deallocate(checked_forward_vector,checked_wavelet_vector)
    placement=Transform_Placement();generation=-1
  end subroutine

  subroutine prepare_resident_transform(topology)
    integer(int64),intent(in) :: topology
    type(Transform_Geometry) :: geometry
    type(Rhs_Geometry_Directory) :: directory
    type(Rhs_Boundary_Map) :: mapping
    type(Boundary_Inventory) :: inventory
    type(Transform_Homes) :: homes
    integer,allocatable :: physical(:,:),metadata(:,:),owners(:),keys(:,:),request(:,:),selected(:),pkeys(:,:)
    real(dp),allocatable :: delivered(:,:),boundary_metadata(:,:)
    integer :: d,i,j,n,g,b,l,m,p,a,e,ng,ns,nf,np,pass,ierr,phase
    logical :: rebuild
    real(dp) :: started,stamp,segments(8)
    rebuild=generation/=topology.or.masks/=connected_rhs_mask_generation()
    if(.not.rebuild) return
    started=MPI_Wtime()
    call detail_enter(DP_TRANSFORM_PLAN)
    call release_backend
    generation=topology;masks=connected_rhs_mask_generation()
    call MPI_Comm_dup(comm,context,ierr)
    if(ierr/=MPI_SUCCESS) error stop 'transform context create'
    allocate(base(size(grid)+1));base(1)=0
    do d=1,size(grid)
       base(d+1)=base(d)+grid(d)%node%length
    end do
    ng=base(size(base))
    call prepare_rhs_geometry_directory(base,topology,directory)
    stamp=MPI_Wtime();segments(1)=stamp-started
    call compile_transform_geometry(geometry)
    segments(2)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    call distribute_transform_geometry(geometry,placement)
    segments(3)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    nnode=size(placement%key,2)
    allocate(physical(4,ng),metadata(10,ng));metadata=0
    do d=1,size(grid)
       b=base(d);n=grid(d)%node%length;g=glo_id(rank+1,d)+1
       call compile_rhs_owner_map(grid(d),g-1,owners,keys)
       do i=1,n
          physical(:,b+i)=[g,i,1,1]
          metadata(1,b+i)=merge(1,3,keys(1,i)>0)
          metadata(2,b+i)=grid(d)%level%elts(i)
          metadata(3:6,b+i)=[grid(d)%mask_n%elts(i),grid(d)%mask_e%elts(3*i-2:3*i)]
       end do
       do p=1,grid(d)%patch%length
          if(grid(d)%patch%elts(p)%deleted) cycle
          a=grid(d)%patch%elts(p)%elts_start
          if(keys(1,a+1)>0) cycle
          metadata(1,b+a+1:b+a+PATCH_SIZE**2)=2
       end do
    end do
    allocate(scalar_boundary(level_start-2:level_end),vector_boundary(level_start-2:level_end), &
         scalar_halo(level_start-2:level_end),vector_halo(level_start-2:level_end),outer_halo(level_start-2:level_end))
    ! Two all-level ranges: NONE, and level_start:level_end for the initial
    ! inverse. Local ordered copies still span all levels in both cases.
    do l=level_start-2,level_end
       do m=1,3,2
          if(l==level_start-2) then
             call compile_rhs_boundary_map(merge(AT_NODE,AT_EDGE,m==1),NONE,base,mapping)
          else if(l==level_start-1) then
             call compile_rhs_boundary_map(merge(AT_NODE,AT_EDGE,m==1),NONE,base,mapping,first_level=level_start)
          else
             call compile_rhs_boundary_map(merge(AT_NODE,AT_EDGE,m==1),l,base,mapping)
          end if
          if(allocated(boundary_metadata)) deallocate(boundary_metadata,request)
          allocate(boundary_metadata(m*ng,5),request(4,m*nnode))
          boundary_metadata(:,:4)=real(transpose(mapping%key),dp);boundary_metadata(:,5)=1
          do i=1,size(mapping%target)
             a=mapping%target(i);boundary_metadata(a,:4)=real(mapping%source(:,i),dp)
             boundary_metadata(a,5)=mapping%orientation(i)
             metadata(7,(a-1)/m+1)=1
          end do
          ! Local alias destinations are all-level even for a single-level
          ! refresh. They are identified separately for outer canonical saves.
          if(m==3) then
             do d=1,size(grid)
                do j=1,size(grid)
                   g=glo_id(rank+1,j)+1
                   do i=1,grid(d)%pack(AT_EDGE,g)%length
                      a=abs(grid(j)%unpk(AT_EDGE,glo_id(rank+1,d)+1)%elts(i))/3+1
                      metadata(8,base(j)+a)=1
                   end do
                end do
             end do
          end if
          do i=1,nnode
             do e=1,m
                request(:,m*(i-1)+e)=placement%key(:,i);request(3,m*(i-1)+e)=e
             end do
          end do
          phase=8700+10*(l-level_start+2)+m
          call deliver_rhs_metadata(mapping%key,boundary_metadata,request,delivered,topology,phase,directory,m)
          if(m==1) then
             scalar_boundary(l)%key=transpose(nint(delivered(:,:4)));scalar_boundary(l)%sign=nint(delivered(:,5))
          else
             vector_boundary(l)%key=transpose(nint(delivered(:,:4)));vector_boundary(l)%sign=nint(delivered(:,5))
          end if
       end do
    end do
    do i=1,size(geometry%outer)
       associate(op=>geometry%outer(i))
       d=loc_id(op%global_domain)+1
       do j=2,3
          a=(op%target(j)-1)/EDGE+1
          metadata(9,base(d)+a)=1
          if(grid(d)%level%elts(a)/=op%level+1) then
             metadata(7,base(d)+a)=1;metadata(8,base(d)+a)=1
          end if
       end do
       end associate
    end do
    call deliver_rhs_metadata(physical,real(transpose(metadata),dp),placement%key,delivered,topology,8790,directory,1)
    info=transpose(nint(delivered))
    selected=[(i,i=1,nnode)]
    save_node=pack(selected,placement%owned.and.info(1,:)==1.and.info(7,:)==1)
    fixed_node=pack(selected,placement%owned.and.info(1,:)==2)
    interior_node=pack(selected,placement%owned.and.info(1,:)==1)
    ns=size(save_node);nf=size(fixed_node)
    allocate(save_index(nnode),fixed_index(nnode));save_index=0;fixed_index=0
    do i=1,ns
       save_index(save_node(i))=nnode+i
    end do
    do i=1,nf
       fixed_index(fixed_node(i))=nnode+ns+i
    end do
    call compile_local_block_inverse_routes(placement%state_key(:,interior_node),interior_address)
    segments(4)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    ! Canonical input ownership is already local and bound by interior_address.
    ! Do not rediscover providers or compile a general exchange for local reads.
    segments(5)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    ! Publication requests name only surviving compatibility aliases/scaffold
    ! and compact block boundary readers. No complete canonical Domain mirror.
    do pass=1,2
       inventory%pass=pass;inventory%used=0
       call apply_local_block_field_consumer(collect_boundary,inventory)
       if(pass==1) allocate(inventory%physical(4,inventory%used),inventory%state_key(4,inventory%used))
    end do
    nboundary=inventory%used
    selected=[(i,i=1,ng)]
    selected=pack(selected,metadata(1,:)/=1)
    np=size(selected);npublication=np+nboundary
    allocate(publication_domain(np),publication_node(np),pkeys(4,npublication))
    pkeys(:,:np)=physical(:,selected);pkeys(:,np+1:)=inventory%physical
    do i=1,np
       publication_domain(i)=loc_id(pkeys(1,i))+1;publication_node(i)=pkeys(2,i)
    end do
    call compile_local_block_inverse_routes(inventory%state_key,boundary_address)
    call compile_fixed_inputs(physical,metadata,directory)
    segments(6)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    call compile_transform_homes(pkeys,homes)
    call compile_publication_rows(pkeys,homes)
    segments(7)=MPI_Wtime()-stamp;stamp=MPI_Wtime()
    allocate(scalar(nrow,2*zlevels),vector(3*nrow,zlevels), &
         scalar_wavelet(nnode,2*zlevels),vector_wavelet(3*nnode,zlevels))
    scalar=0;vector=0;scalar_wavelet=0;vector_wavelet=0
    if(resident_transform_check_enabled()) allocate(checked_forward_vector(3*nnode,zlevels), &
         checked_wavelet_vector(3*nnode,zlevels))
    allocate(wire(max(1,3*zlevels*max(size(interior_node),nboundary))))
    do l=level_start-2,level_end
       call bind_halo(scalar_halo(l),1,l,.false.,homes)
       call bind_halo(vector_halo(l),3,l,.false.,homes)
       if(l>level_start) call bind_halo(outer_halo(l),3,l,.true.,homes)
    end do
    call bind_halo(scalar_wave_halo,1,level_start-1,.false.,homes,wave=.true.)
    call bind_halo(vector_wave_halo,3,level_start-1,.false.,homes,wave=.true.)
    call bind_publication(scalar_publication,1,pkeys,homes)
    call bind_publication(vector_publication,3,pkeys,homes)
    selected=[(i,i=1,nnode)];selected=pack(selected,placement%owned.and.info(1,:)/=1)
    call compile_transform_adapter(scalar_adapter,1,9501,scalar_boundary(level_start-2)%key(:,selected), &
         selected,scalar_boundary(level_start-2)%sign(selected))
    selected=reshape(spread(3*selected-3,1,3)+spread([1,2,3],2,size(selected)),[3*size(selected)])
    call compile_transform_adapter(vector_adapter,3,9503,vector_boundary(level_start-2)%key(:,selected), &
         selected,vector_boundary(level_start-2)%sign(selected))
    call bind_halo(scalar_initial,1,level_start-2,.false.,homes,identity=.true.)
    call bind_halo(vector_initial,3,level_start-2,.false.,homes,identity=.true.)
    call prepare_exchange_scratch
    call compile_forward_order
    segments(8)=MPI_Wtime()-stamp
    if(audit_on) write(*,'(a,3(i0,1x),8(es13.5,1x))') 'Transform279 preparation: ',rank,generation,masks,segments
    builds=builds+1;preparation_seconds=MPI_Wtime()-started
    call audit_transform_plan
    call detail_leave(DP_TRANSFORM_PLAN)
    if(resident_transform_check_enabled().and.rank==0) write(*,'(a)') 'Resident provisional transform plan prepared'
  end subroutine

  subroutine collect_boundary(catalog,block,context)
    integer,intent(in) :: catalog
    type(Block_Data),intent(in) :: block
    class(*),intent(inout) :: context
    integer :: p,q,i
    select type(c=>context)
    type is(Boundary_Inventory)
       do p=1,size(block%bdry_storage)
          do q=0,block%bdry_storage(p)%n_node-1
             c%used=c%used+1;i=c%used
             if(c%pass/=2) cycle
             c%physical(:,i)=[block%root_domain+1,block%bdry_storage(p)%elts_start+q+1,1,1]
             c%state_key(:,i)=[catalog,STORE_BDRY,p,q]
          end do
       end do
    class default
       error stop 'transform boundary inventory context'
    end select
  end subroutine

  integer function transform_key_bucket(key,n) result(h)
    integer,intent(in) :: key(:),n
    integer(int64) :: bits
    bits=modulo(int(key(2),int64)*104729_int64+int(key(1),int64)*65599_int64,2147483647_int64)
    bits=ieor(bits,shiftr(bits,16))
    bits=modulo(bits*104729_int64,2147483647_int64)
    h=int(modulo(bits,int(n,int64)))+1
  end function

  subroutine compile_transform_homes(publication,homes)
    integer,intent(in) :: publication(:,:)
    type(Transform_Homes),intent(out) :: homes
    type(Rhs_Owner_Directory) :: directory
    integer,allocatable :: selected(:)
    integer :: l,n
    ! Resolve unique physical nodes once per epoch. Component homes are derived
    ! locally; every phase reuses v7's direct final-writer binding path.
    n=nnode+size(publication,2)
    do l=lbound(scalar_boundary,1),ubound(scalar_boundary,1)
       n=n+size(scalar_boundary(l)%key,2)+size(vector_boundary(l)%key,2)
    end do
    allocate(homes%key(4,n),homes%table(max(1,2*n+1)));homes%table=0
    call insert(placement%key);call insert(publication)
    do l=lbound(scalar_boundary,1),ubound(scalar_boundary,1)
       call insert(scalar_boundary(l)%key);call insert(vector_boundary(l)%key)
    end do
    selected=pack([(l,l=1,nnode)],placement%owned)
    call compile_rhs_owner_directory(directory,context,placement%key(:,selected),selected)
    call lookup_rhs_owner_directory(directory,homes%key(:,:homes%used),homes%home)
  contains
    subroutine insert(key)
      integer,intent(in) :: key(:,:)
      integer :: i,h,a
      do i=1,size(key,2)
         h=transform_key_bucket(key(:,i),size(homes%table))
         do
            a=homes%table(h)
            if(a==0) then
               homes%used=homes%used+1;a=homes%used;homes%table(h)=a
               homes%key(:,a)=[key(1:2,i),1,1];exit
            end if
            if(all(homes%key(1:2,a)==key(1:2,i))) exit
            h=modulo(h,size(homes%table))+1
         end do
      end do
    end subroutine
  end subroutine

  integer function transform_home_index(homes,key) result(a)
    type(Transform_Homes),intent(in) :: homes
    integer,intent(in) :: key(:)
    integer :: h
    h=transform_key_bucket(key,size(homes%table))
    do
       a=homes%table(h)
       if(a==0) error stop 'transform home was not registered'
       if(all(homes%key(1:2,a)==key(1:2))) return
       h=modulo(h,size(homes%table))+1
    end do
  end function

  subroutine compile_publication_rows(physical,homes)
    integer,intent(in) :: physical(:,:)
    type(Transform_Homes),intent(in) :: homes
    integer,allocatable :: imported(:)
    integer :: i,a
    allocate(publication_row(size(physical,2)),imported(homes%used));imported=0
    nrow=nnode+size(save_node)+size(fixed_node)
    do i=1,size(physical,2)
       a=transform_home_index(homes,physical(:,i))
       if(homes%home(1,a)==rank) then
          publication_row(i)=homes%home(2,a)
          if(publication_row(i)<1.or.publication_row(i)>nnode) error stop 'local publication address'
       else
          if(imported(a)==0) then
             nrow=nrow+1;imported(a)=nrow
          end if
          publication_row(i)=imported(a)
       end if
    end do
    if(audit_on) write(*,'(a,4(i0,1x))') 'Transform280 publication: ',rank,size(publication_row), &
         count(imported>0),count(publication_row<=nnode)
  end subroutine

  subroutine transform_homes_for(homes,key,m,result)
    type(Transform_Homes),intent(in) :: homes
    integer,intent(in) :: key(:,:),m
    integer,allocatable,intent(out) :: result(:,:)
    integer :: i,a
    allocate(result(2,size(key,2)))
    do i=1,size(key,2)
       a=transform_home_index(homes,key(:,i))
       result(:,i)=[homes%home(1,a),m*(homes%home(2,a)-1)+key(3,i)]
    end do
  end subroutine

  subroutine bind_halo(route,m,level,outer,homes,wave,identity)
    type(Transform_Route),intent(inout) :: route
    integer,intent(in) :: m,level
    logical,intent(in) :: outer
    type(Transform_Homes),intent(in) :: homes
    integer,allocatable :: definition_home(:,:),consumer_home(:,:),canonical(:)
    logical,optional,intent(in) :: wave,identity
    integer,allocatable :: priority(:),selected(:),definition(:,:),order(:,:),source(:),target(:),request(:,:),sign(:)
    integer :: i,e,a,p,j,n
    logical,allocatable :: needed(:)
    logical :: is_wave
    is_wave=.false.
    if(present(wave)) is_wave=wave
    if(outer) then
       allocate(priority(m*nnode));priority=0
       do i=1,nnode
          if(placement%owned(i)) priority(m*(i-1)+1:m*i)=1
       end do
       do i=1,size(placement%program%outer)
          associate(op=>placement%program%outer(i))
          if(op%level/=level-1) cycle
          do j=2,3
             a=op%target(j)
             if(info(1,(a-1)/3+1)==1.and..not.placement%owned((a-1)/3+1)) &
                  error stop 'transform component crosses canonical block ownership'
             priority(a)=max(priority(a),placement%outer_order(i)+1)
          end do
          end associate
       end do
       selected=[(i,i=1,m*nnode)];selected=pack(selected,priority>0);n=size(selected)
       allocate(definition(4,n),order(2,n),source(n))
       do i=1,n
          a=selected(i);p=(a-1)/m+1;e=mod(a-1,m)+1
          definition(:,i)=placement%key(:,p);definition(3,i)=e
          order(:,i)=[1,priority(a)];source(i)=a
       end do
    else
       allocate(definition(4,0),order(2,0),source(0))
    end if
    allocate(target(m*nnode))
    target=[(i,i=1,m*nnode)]
    if(m==1) then
       request=scalar_boundary(level)%key;sign=scalar_boundary(level)%sign
    else
       request=vector_boundary(level)%key;sign=vector_boundary(level)%sign
    end if
    if(present(identity)) then
       if(identity) then
          do i=1,m*nnode
             p=(i-1)/m+1;e=mod(i-1,m)+1
             request(:,i)=[placement%key(1:2,p),e,1]
          end do
          sign=1
       end if
    end if
    allocate(needed(m*nnode));needed=.false.
    do i=1,m*nnode
       p=(i-1)/m+1;e=mod(i-1,m)+1
       needed(i)=.not.placement%owned(p).or.any(request(:,i)/=[placement%key(1:2,p),e,1]).or.sign(i)/=1
       if(outer.and.info(9,p)==1.and.info(1,p)/=1) needed(i)=.true.
    end do
    selected=pack(target,needed);request=request(:,selected);sign=sign(selected);target=selected
    route%phase=9000+10*(level-level_start+2)+m
    if(outer) route%phase=route%phase+100
    if(is_wave) route%phase=route%phase+200
    if(present(identity)) then
       if(identity) route%phase=route%phase+500
    end if
    call transform_homes_for(homes,request,m,consumer_home)
    if(outer) then
       call transform_homes_for(homes,definition,m,definition_home)
       canonical=pack([(i,i=1,nnode)],placement%owned)
       canonical=reshape(spread(m*canonical-m,1,m)+spread([(i,i=1,m)],2,size(canonical)),[m*size(canonical)])
       call compile_rhs_exchange(route%exchange,context,generation,masks,route%phase,m*nrow, &
            merge(2*zlevels,zlevels,m==1),definition,order,source,request,target,sign, &
            shared_tag=route%phase,profile_region=DP_TRANSFORM_EXCHANGE,shared_buffers=.true., &
            definition_home=definition_home,consumer_home=consumer_home,canonical_slot=canonical)
    else
       call compile_rhs_exchange(route%exchange,context,generation,masks,route%phase,m*merge(nnode,nrow,is_wave), &
            merge(2*zlevels,zlevels,m==1),definition,order,source,request,target,sign, &
            shared_tag=route%phase,profile_region=DP_TRANSFORM_EXCHANGE,shared_buffers=.true., &
            consumer_provider=consumer_home)
    end if
  end subroutine

  subroutine bind_publication(route,m,physical,homes)
    type(Transform_Route),intent(inout) :: route
    integer,intent(in) :: m,physical(:,:)
    type(Transform_Homes),intent(in) :: homes
    integer,allocatable :: consumer_home(:,:),definition(:,:),source(:),order(:,:),request(:,:),target(:),sign(:)
    logical,allocatable :: seen(:)
    integer :: i,e,n,a,p,pass
    allocate(seen(nrow));n=0
    do pass=1,2
       seen=.false.;a=0
       do i=1,size(publication_row)
          p=publication_row(i)
          if(p<=nnode.or.seen(p)) cycle
          seen(p)=.true.
          do e=1,m
             a=a+1
             if(pass==1) cycle
             request(:,a)=physical(:,i);request(3,a)=e;target(a)=m*(p-1)+e
          end do
       end do
       if(pass==1) then
          n=a;allocate(request(4,n),target(n),sign(n),definition(4,0),order(2,0),source(0))
       end if
    end do
    sign=1;route%phase=9300+m
    call transform_homes_for(homes,request,m,consumer_home)
    call compile_rhs_exchange(route%exchange,context,generation,masks,route%phase,m*nrow, &
         merge(2*zlevels,zlevels,m==1),definition,order,source,request,target,sign, &
         shared_tag=route%phase,profile_region=DP_TRANSFORM_EXCHANGE,shared_buffers=.true., &
         consumer_provider=consumer_home)
  end subroutine

  subroutine compile_fixed_inputs(physical,metadata,directory)
    integer,intent(in) :: physical(:,:),metadata(:,:)
    type(Rhs_Geometry_Directory),intent(in) :: directory
    integer,allocatable :: selected(:),definition(:,:),request(:,:),source(:),target(:),order(:,:),sign(:)
    integer :: i,m,e,a,nf,nslot,phase
    ! Explicit allocation also handles an empty selection without descriptor
    ! self-reallocation diagnostics in the GCC 13 checked build.
    nfixed=count(metadata(1,:)==2)
    allocate(selected(nfixed));a=0
    do i=1,size(physical,2)
       if(metadata(1,i)/=2) cycle
       a=a+1;selected(a)=i
    end do
    nf=size(fixed_node);nslot=nfixed+nf
    allocate(fixed_domain(nfixed),fixed_physical_node(nfixed))
    do i=1,nfixed
       fixed_domain(i)=loc_id(physical(1,selected(i)))+1
       fixed_physical_node(i)=physical(2,selected(i))
    end do
    allocate(fixed_scalar(nslot,2*zlevels),fixed_vector(3*nslot,zlevels))
    do m=1,3,2
       if(allocated(definition)) deallocate(definition,request,source,target,order,sign)
       allocate(definition(4,m*nfixed),request(4,m*nf),source(m*nfixed),target(m*nf),order(2,m*nfixed),sign(m*nf))
       do i=1,nfixed
          do e=1,m
             a=m*(i-1)+e;definition(:,a)=physical(:,selected(i));definition(3,a)=e;source(a)=a
          end do
       end do
       do i=1,nf
          do e=1,m
             a=m*(i-1)+e;request(:,a)=placement%key(:,fixed_node(i));request(3,a)=e
             target(a)=m*(nfixed+i-1)+e
          end do
       end do
       order=1;sign=1;phase=9400+m
       if(m==1) then
          call compile_rhs_exchange(fixed_scalar_route,context,generation,masks,phase,nslot,2*zlevels, &
               definition,order,source,request,target,sign,shared_tag=phase, &
               profile_region=DP_TRANSFORM_EXCHANGE,shared_buffers=.true.)
       else
          call compile_rhs_exchange(fixed_vector_route,context,generation,masks,phase,3*nslot,zlevels, &
               definition,order,source,request,target,sign,shared_tag=phase, &
               profile_region=DP_TRANSFORM_EXCHANGE,shared_buffers=.true.)
       end if
    end do
  end subroutine

  subroutine compile_forward_order
    integer,allocatable :: table(:)
    integer :: i,h,j,n
    integer(int64) :: key
    n=size(placement%program%inner)
    allocate(table(max(1,2*n+1)),forward_inner(size(placement%program%outer)));table=0;forward_inner=0
    do i=1,n
       associate(op=>placement%program%inner(i))
       key=int(op%global_domain,int64)*65599_int64+int(op%cell,int64)*104729_int64
       h=int(modulo(key,int(size(table),int64)))+1
       do while(table(h)/=0)
          h=modulo(h,size(table))+1
       end do
       table(h)=i
       end associate
    end do
    do i=1,size(forward_inner)
       associate(op=>placement%program%outer(i))
       if(op%level<level_start.or.op%kind/=0.or.op%cell==0.or.op%component/=EDGE) cycle
       key=int(op%global_domain,int64)*65599_int64+int(op%cell,int64)*104729_int64
       h=int(modulo(key,int(size(table),int64)))+1
       do
          j=table(h)
          if(j==0) error stop 'transform forward cell has a different owner'
          if(placement%program%inner(j)%global_domain==op%global_domain.and. &
               placement%program%inner(j)%cell==op%cell) exit
          h=modulo(h,size(table))+1
       end do
       forward_inner(i)=j
       end associate
    end do
    if(count(forward_inner>0)/=n) error stop 'transform forward inner coverage'
  end subroutine

  subroutine prepare_exchange_scratch
    integer :: n,l
    ! Every exchange drains its requests before returning. Numerical send and
    ! receive storage is bounded by the largest route, rather than their sum.
    n=max(rhs_exchange_buffer_values(fixed_scalar_route),rhs_exchange_buffer_values(fixed_vector_route), &
         rhs_exchange_buffer_values(scalar_initial%exchange),rhs_exchange_buffer_values(vector_initial%exchange), &
         rhs_exchange_buffer_values(scalar_wave_halo%exchange),rhs_exchange_buffer_values(vector_wave_halo%exchange), &
         rhs_exchange_buffer_values(scalar_publication%exchange),rhs_exchange_buffer_values(vector_publication%exchange))
    do l=lbound(scalar_halo,1),ubound(scalar_halo,1)
       n=max(n,rhs_exchange_buffer_values(scalar_halo(l)%exchange),rhs_exchange_buffer_values(vector_halo(l)%exchange))
       if(l>level_start) n=max(n,rhs_exchange_buffer_values(outer_halo(l)%exchange))
    end do
    allocate(exchange_scratch(n))
  end subroutine

  subroutine communicate(route,values)
    type(Transform_Route),intent(inout) :: route
    real(dp),intent(inout) :: values(:,:)
    state=state+1
    call exchange_rhs_values(route%exchange,generation,masks,route%phase,state,values,exchange_scratch)
    call count_route(route%exchange)
  end subroutine

  subroutine reset_aliases(component,include_fixed)
    integer,intent(in) :: component
    logical,intent(in) :: include_fixed
    integer :: i,p,q
    do i=1,size(save_node)
       p=save_node(i);q=save_index(p)
       if(component==1) then
          scalar(p,:)=scalar(q,:)
       else
          vector(3*p-2:3*p,:)=vector(3*q-2:3*q,:)
       end if
    end do
    if(.not.include_fixed) return
    do i=1,size(fixed_node)
       p=fixed_node(i);q=fixed_index(p)
       if(component==1) then
          scalar(p,:)=scalar(q,:)
       else
          vector(3*p-2:3*p,:)=vector(3*q-2:3*q,:)
       end if
    end do
  end subroutine

  subroutine save_scalar_target(p)
    integer,intent(in) :: p
    if(save_index(p)>0) scalar(save_index(p),:)=scalar(p,:)
  end subroutine

  subroutine save_vector_target(edge_code)
    integer,intent(in) :: edge_code
    integer :: p,e,q
    p=(edge_code-1)/3+1;e=mod(edge_code-1,3)+1;q=save_index(p)
    if(q>0) vector(3*(q-1)+e,:)=vector(edge_code,:)
  end subroutine

  subroutine execute_resident_transform(scaling,wavelet,commit)
    type(Float_Field),intent(inout) :: scaling(N_VARIABLE,zlevels),wavelet(N_VARIABLE,zlevels)
    logical,intent(in) :: commit
    integer :: i,j,k,p,q,d,e,l
    if(generation<0.or.context==MPI_COMM_NULL) error stop 'unprepared resident transform'
    sent_values=0;received_values=0;messages=0;rounds=0
    call detail_enter(DP_TRANSFORM_INPUT)
    state=state+1
    scalar(:nnode,:)=0.0_dp;vector(:3*nnode,:)=0.0_dp
    call transfer_local_block_inverse_routes(interior_address,BLOCK_PAYLOAD_SOL,1,.false.,wire)
    do i=1,size(interior_node)
       p=interior_node(i)
       do k=1,zlevels
          scalar(p,2*k-1:2*k)=wire((i-1)*2*zlevels+2*k-1:(i-1)*2*zlevels+2*k)
       end do
    end do
    call transfer_local_block_inverse_routes(interior_address,BLOCK_PAYLOAD_SOL,2,.false.,wire)
    do i=1,size(interior_node)
       p=interior_node(i)
       do k=1,zlevels
          vector(3*p-2:3*p,k)=wire((i-1)*3*zlevels+3*k-2:(i-1)*3*zlevels+3*k)
       end do
    end do
    do i=1,size(save_node)
       p=save_node(i);q=save_index(p)
       scalar(q,:)=scalar(p,:);vector(3*q-2:3*q,:)=vector(3*p-2:3*p,:)
    end do
    do i=1,nfixed
       d=fixed_domain(i);p=fixed_physical_node(i)
       do k=1,zlevels
          fixed_scalar(i,2*k-1)=sol(scalars(1),k)%data(d)%elts(p)
          fixed_scalar(i,2*k)=sol(scalars(2),k)%data(d)%elts(p)
          fixed_vector(3*i-2:3*i,k)=sol(S_VELO,k)%data(d)%elts(3*p-2:3*p)
       end do
    end do
    state=state+1
    call exchange_rhs_values(fixed_scalar_route,generation,masks,9401,state,fixed_scalar,exchange_scratch)
    call exchange_rhs_values(fixed_vector_route,generation,masks,9403,state,fixed_vector,exchange_scratch)
    call count_route(fixed_scalar_route);call count_route(fixed_vector_route)
    do i=1,size(fixed_node)
       p=fixed_index(fixed_node(i));q=nfixed+i
       scalar(p,:)=fixed_scalar(q,:);vector(3*p-2:3*p,:)=fixed_vector(3*q-2:3*q,:)
    end do
    call detail_enter(DP_TRANSFORM_EXCHANGE)
    call read_transform_adapter(scalar_adapter,context,scaling,scalar)
    call read_transform_adapter(vector_adapter,context,scaling,vector)
    call count_adapter(scalar_adapter);call count_adapter(vector_adapter)
    call detail_leave(DP_TRANSFORM_EXCHANGE)
    call communicate(scalar_initial,scalar)
    call communicate(vector_initial,vector)
    if(allocated(checked_forward_vector)) checked_forward_vector=vector(:3*nnode,:)
    call detail_leave(DP_TRANSFORM_INPUT)
    call detail_enter(DP_TRANSFORM_FORWARD)
    scalar_wavelet=0.0_dp;vector_wavelet=0.0_dp
    do i=1,size(placement%program%scalar)
       call execute_transform_scalar(placement%program%scalar(i),scalar,scalar_wavelet,.false.)
    end do
    do l=level_end-1,level_start,-1
       do i=1,size(placement%program%outer)
          if(placement%program%outer(i)%level/=l) cycle
          call execute_transform_outer(placement%program%outer(i),vector,vector_wavelet,.false.)
          j=forward_inner(i)
          if(j>0) call execute_transform_inner(placement%program%inner(j),vector,vector_wavelet,.false.)
       end do
    end do
    if(allocated(checked_wavelet_vector)) checked_wavelet_vector=vector_wavelet(:3*nnode,:)
    ! Only canonical patch coefficients survive forward activation. Alias
    ! values are established by the exact initial inverse boundary routes.
    do p=1,nnode
       if(info(1,p)/=1) then
          scalar_wavelet(p,:)=0.0_dp;vector_wavelet(3*p-2:3*p,:)=0.0_dp
       else if(info(2,p)>level_start) then
          if(info(3,p)<ADJZONE) scalar_wavelet(p,:)=0.0_dp
          do e=1,3
             if(info(3+e,p)<ADJZONE) vector_wavelet(3*(p-1)+e,:)=0.0_dp
          end do
       end if
    end do
    if(commit) call commit_interiors(BLOCK_PAYLOAD_WAV_COEFF)
    call detail_leave(DP_TRANSFORM_FORWARD)
    call detail_enter(DP_TRANSFORM_INVERSE)
    call reset_aliases(1,.false.);call reset_aliases(2,.false.)
    call communicate(scalar_halo(level_start-1),scalar)
    call communicate(vector_halo(level_start-1),vector)
    call communicate(scalar_wave_halo,scalar_wavelet)
    call communicate(vector_wave_halo,vector_wavelet)
    call reset_aliases(1,.true.)
    call communicate(scalar_halo(level_start),scalar)
    do l=level_start,level_end-1
       do i=1,size(placement%program%scalar)
          associate(op=>placement%program%scalar(i))
          if(op%level/=l.or.op%phase/=1) cycle
          call execute_transform_scalar(op,scalar,scalar_wavelet,.true.)
          if(op%inverse_active/=0) call save_scalar_target(op%target)
          end associate
       end do
       call reset_aliases(1,.true.);call communicate(scalar_halo(l+1),scalar)
       call reset_aliases(2,.true.);call communicate(vector_halo(l),vector)
       do i=1,size(placement%program%outer)
          if(placement%program%outer(i)%level/=l) cycle
          call execute_transform_outer(placement%program%outer(i),vector,vector_wavelet,.true.)
       end do
       call communicate(outer_halo(l+1),vector)
       do i=1,size(save_node)
          p=save_node(i);q=save_index(p)
          if(info(2,p)==l.or.info(2,p)==l+1.or.info(8,p)==1) vector(3*q-2:3*q,:)=vector(3*p-2:3*p,:)
       end do
       do i=1,size(placement%program%scalar)
          associate(op=>placement%program%scalar(i))
          if(op%level/=l.or.op%phase/=2) cycle
          call execute_transform_scalar(op,scalar,scalar_wavelet,.true.)
          call save_scalar_target(op%target)
          end associate
       end do
       call reset_aliases(1,.true.);call communicate(scalar_halo(l+1),scalar)
       do i=1,size(placement%program%inner)
          associate(op=>placement%program%inner(i))
          if(op%level/=l) cycle
          call execute_transform_inner(op,vector,vector_wavelet,.true.)
          do j=1,6
             call save_vector_target(op%target(j))
          end do
          end associate
       end do
       call reset_aliases(2,.true.);call communicate(vector_halo(l+1),vector)
    end do
    call reset_aliases(1,.true.);call reset_aliases(2,.true.)
    call communicate(scalar_halo(level_start-2),scalar)
    call communicate(vector_halo(level_start-2),vector)
    call detail_leave(DP_TRANSFORM_INVERSE)
    call detail_enter(DP_TRANSFORM_OUTPUT)
    if(commit) call commit_interiors(BLOCK_PAYLOAD_SOL)
    call communicate(scalar_publication,scalar);call communicate(vector_publication,vector)
    completed=completed+1
    if(.not.commit) then
       call detail_leave(DP_TRANSFORM_OUTPUT)
       return
    end if
    call install_transform_outputs(scaling)
    call audit_transform_execution(.false.)
    call detail_leave(DP_TRANSFORM_OUTPUT)
  end subroutine

  subroutine install_checked_transform(scaling)
    ! The independent path has checked every canonical/adapter result. Restore
    ! the candidate's outputs so subsequent RHS readers consume this backend.
    type(Float_Field),intent(inout) :: scaling(N_VARIABLE,zlevels)
    call commit_interiors(BLOCK_PAYLOAD_SOL)
    call install_transform_outputs(scaling)
    call audit_transform_execution(.true.)
  end subroutine

  subroutine install_transform_outputs(scaling)
    type(Float_Field),intent(inout) :: scaling(N_VARIABLE,zlevels)
    integer :: np,i,d,p,q,k
    np=size(publication_node)
    do i=1,np
       d=publication_domain(i);p=publication_node(i);q=publication_row(i)
       do k=1,zlevels
          scaling(scalars(1),k)%data(d)%elts(p)=scalar(q,2*k-1)
          scaling(scalars(2),k)%data(d)%elts(p)=scalar(q,2*k)
          scaling(S_VELO,k)%data(d)%elts(3*p-2:3*p)=vector(3*q-2:3*q,k)
       end do
    end do
    do i=1,nboundary
       q=publication_row(np+i)
       do k=1,zlevels
          wire((i-1)*2*zlevels+2*k-1:(i-1)*2*zlevels+2*k)=scalar(q,2*k-1:2*k)
       end do
    end do
    call transfer_local_block_inverse_routes(boundary_address,BLOCK_PAYLOAD_SOL,1,.true.,wire)
    do i=1,nboundary
       q=publication_row(np+i)
       do k=1,zlevels
          wire((i-1)*3*zlevels+3*k-2:(i-1)*3*zlevels+3*k)=vector(3*q-2:3*q,k)
       end do
    end do
    call transfer_local_block_inverse_routes(boundary_address,BLOCK_PAYLOAD_SOL,2,.true.,wire)
    scaling%bdry_uptodate=.true.
  end subroutine

  subroutine count_route(route)
    type(Rhs_Exchange),intent(in) :: route
    integer(int64) :: sends,receives,bytes
    integer :: sp,rp
    if(.not.audit_on) return
    call rhs_exchange_cost(route,sends,receives,sp,rp,bytes)
    sent_values=sent_values+sends;received_values=received_values+receives;messages=messages+sp;rounds=rounds+1
  end subroutine

  subroutine count_adapter(adapter)
    type(Transform_Adapter),intent(in) :: adapter
    if(.not.audit_on) return
    sent_values=sent_values+int(sum(adapter%rc)-adapter%rc(rank+1),int64)*adapter%columns
    received_values=received_values+int(sum(adapter%sc)-adapter%sc(rank+1),int64)*adapter%columns
    messages=messages+count(adapter%rc>0)-merge(1,0,adapter%rc(rank+1)>0);rounds=rounds+1
  end subroutine

  subroutine audit_transform_execution(reference)
    logical,intent(in) :: reference
    character(64) :: path
    integer :: unit
    if(.not.audit_on) return
    write(path,'(a,i0,a)') 'transform279-production-rank-',rank,'.txt'
    open(newunit=unit,file=trim(path),status='unknown',position='append',action='write')
    write(unit,*) completed,generation,masks,size(grid),size(interior_node),size(placement%program%scalar), &
         size(placement%program%inner),size(placement%program%outer),merge(1,0,reference), &
         sent_values,received_values,messages,rounds,route_bytes,buffer_bytes,1,builds, &
         5_int64*zlevels*size(publication_node),5_int64*zlevels*nboundary
    close(unit)
  end subroutine

  subroutine audit_transform_plan
    character(64) :: path
    integer :: unit,l,i,p,b,sp,rp
    integer,allocatable :: work(:,:)
    integer(int64) :: values,bytes,sends,receives
    if(.not.audit_on) return
    route_bytes=0
    call add(fixed_scalar_route);call add(fixed_vector_route)
    call add(scalar_initial%exchange);call add(vector_initial%exchange)
    call add(scalar_wave_halo%exchange);call add(vector_wave_halo%exchange)
    call add(scalar_publication%exchange);call add(vector_publication%exchange)
    do l=lbound(scalar_halo,1),ubound(scalar_halo,1)
       call add(scalar_halo(l)%exchange);call add(vector_halo(l)%exchange)
       if(l>level_start) call add(outer_halo(l)%exchange)
    end do
    values=size(scalar,kind=int64)+size(vector,kind=int64)+size(scalar_wavelet,kind=int64)+ &
         size(vector_wavelet,kind=int64)+ &
         size(fixed_scalar,kind=int64)+size(fixed_vector,kind=int64)+size(wire,kind=int64)+ &
         size(exchange_scratch,kind=int64)+ &
         size(scalar_adapter%sent,kind=int64)+ &
         size(scalar_adapter%received,kind=int64)+size(vector_adapter%sent,kind=int64)+ &
         size(vector_adapter%received,kind=int64)
    buffer_bytes=8_int64*values
    write(path,'(a,i0,a)') 'transform279-cost-rank-',rank,'.txt'
    open(newunit=unit,file=trim(path),status='unknown',position='append',action='write')
    write(unit,*) generation,masks,nnode,count(placement%owned),size(save_node),size(fixed_node), &
         npublication,route_bytes,buffer_bytes,placement%sent_bytes,placement%received_bytes, &
         1,builds,preparation_seconds
    close(unit)
    ! Per-block operation vectors, with unowned scaffold/alias work in bucket0.
    ! This is a work census, not a cycle-weight or a runtime speedup estimate.
    allocate(work(4,0:size(block_catalog)));work=0
    do i=1,nnode
       if(.not.placement%owned(i)) cycle
       b=max(0,placement%state_key(1,i));work(1,b)=work(1,b)+1
    end do
    do i=1,size(placement%program%scalar)
       p=placement%program%scalar(i)%target;b=max(0,placement%state_key(1,p));work(2,b)=work(2,b)+1
    end do
    do i=1,size(placement%program%inner)
       p=(placement%program%inner(i)%target(1)-1)/3+1
       b=max(0,placement%state_key(1,p));work(3,b)=work(3,b)+1
    end do
    do i=1,size(placement%program%outer)
       b=0
       do l=2,3
          p=(placement%program%outer(i)%target(l)-1)/3+1
          if(placement%state_key(1,p)>0) b=placement%state_key(1,p)
       end do
       work(4,b)=work(4,b)+1
    end do
    write(path,'(a,i0,a)') 'transform279-work-rank-',rank,'.txt'
    open(newunit=unit,file=trim(path),status='unknown',position='append',action='write')
    do b=0,size(block_catalog)
       if(any(work(:,b)>0)) write(unit,*) generation,masks,b,work(:,b)
    end do
    close(unit)
  contains
    subroutine add(route)
      type(Rhs_Exchange),intent(in) :: route
      call rhs_exchange_cost(route,sends,receives,sp,rp,bytes)
      route_bytes=route_bytes+bytes
    end subroutine
  end subroutine

  subroutine commit_interiors(family)
    integer,intent(in) :: family
    integer :: i,p,q,k
    do i=1,size(interior_node)
       p=interior_node(i);q=p
       if(family==BLOCK_PAYLOAD_SOL.and.save_index(p)>0) q=save_index(p)
       do k=1,zlevels
          if(family==BLOCK_PAYLOAD_SOL) then
             wire((i-1)*2*zlevels+2*k-1:(i-1)*2*zlevels+2*k)=scalar(q,2*k-1:2*k)
          else
             wire((i-1)*2*zlevels+2*k-1:(i-1)*2*zlevels+2*k)=scalar_wavelet(q,2*k-1:2*k)
          end if
       end do
    end do
    call transfer_local_block_inverse_routes(interior_address,family,1,.true.,wire)
    do i=1,size(interior_node)
       p=interior_node(i);q=p
       if(family==BLOCK_PAYLOAD_SOL.and.save_index(p)>0) q=save_index(p)
       do k=1,zlevels
          if(family==BLOCK_PAYLOAD_SOL) then
             wire((i-1)*3*zlevels+3*k-2:(i-1)*3*zlevels+3*k)=vector(3*q-2:3*q,k)
          else
             wire((i-1)*3*zlevels+3*k-2:(i-1)*3*zlevels+3*k)=vector_wavelet(3*q-2:3*q,k)
          end if
       end do
    end do
    call transfer_local_block_inverse_routes(interior_address,family,2,.true.,wire)
  end subroutine
end module
