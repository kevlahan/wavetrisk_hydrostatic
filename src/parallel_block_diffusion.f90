module parallel_block_diffusion_mod
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use ieee_arithmetic, only: ieee_value,ieee_quiet_nan
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run,glo_id,block_catalog
  use shared_mod, only: n_domain,level_start,level_end,N_BDRY,AT_NODE,AT_EDGE, &
       S_MASS,S_TEMP,S_DIVU,zlevels,Laplace_divu,Laplace_sclr,z_null,TRSK
  use domain_mod, only: Domain,grid,idx,Laplacian_scalar,Laplacian_vector
  use utils_mod, only: nu_scale
  use domain_ops_mod, only: apply_onescale_to_patch
  use parallel_block_primitive_walk_mod, only: walk_primitive_patch
  use parallel_block_columns_mod, only: owner_tendency_without_reference
  use parallel_block_columns_mod, only: columns_active,columns_oracle
  use parallel_block_inverse_mod, only: inverse_operand_generation,block_rk_enabled, &
       rk_alias_nodes,set_rk_mass_alias_layer
  use parallel_block_mod, only: local_block_patch_count
  use patch_mod, only: PATCH_SIZE
  use parallel_block_resident_operands_mod, only: resident_operand_key
  use parallel_block_mass_mod, only: native_mass
  use parallel_block_program_mod
  implicit none
  private
  public :: owned_mass_node
  public :: owner_mass_source,owner_diffusion_source,owner_diffusion_value
  public :: mass_owner_active,mass_owner_oracle,install_owner_mass
  public :: owner_scalar_flux
  public :: begin_block_diffusion,install_block_diffusion,diffusion_active,diffusion_oracle
  type :: Definition_Field
     integer, allocatable :: ref(:,:)
  end type
  type :: Diffusion_Map
     type(Definition_Field) :: flux,divu,mass,temp,lapdiv,tendency,initial_rhs
     type(Definition_Field),allocatable :: mass_flux(:),mass_rhs(:)
     integer, allocatable :: primary(:,:,:),scalar_slot(:,:,:,:),scalar_ref(:,:,:,:,:)
  end type
  type :: Mass_Delivery
     type(Owner_Read) :: input
     integer,allocatable :: target(:,:)
  end type
  type(Mass_Delivery),allocatable :: mass_output(:)
  type :: Owned_Mass_Map
     integer,allocatable::slot(:,:)
  end type
  type(Owned_Mass_Map),allocatable::owned_mass(:)
  type(Owner_Read)::owned_mass_read,alias_mass_read
  real(dp),allocatable::owned_mass_value(:,:)
  logical,allocatable::owned_mass_ready(:)
  logical :: mass_owner_active=.false.,mass_owner_oracle=.false.
  integer(int64)::mass_checked=0_int64
  integer::zero_definition(5)=0
  type(Diffusion_Map),allocatable :: map(:)
  type(Owner_Program) :: program
  type(Owner_Read) :: delivery,scalar_delivery
  integer,allocatable :: destination(:,:)
  integer(int64) :: generation=-1_int64,checked=0_int64
  logical :: diffusion_active=.false.,diffusion_oracle=.false.
  integer :: compiling_domain=0,compiling_field=0,compiling_level=0
  logical :: coverage_pass=.false.
  type(Definition_Field),allocatable :: initial_flux(:)
contains
  subroutine die(message)
    character(*),intent(in)::message
    write(*,'(a,i0,2a)') 'Rank ',rank,': block diffusion: ',message
    call abort_run
    error stop 'block diffusion abort returned'
  end subroutine
  integer function target_owner(d,id) result(owner)
    integer,intent(in)::d,id
    integer::key(5),b
    key=resident_operand_key(d,id)
    if(key(2)>0) then
       owner=key(1);return
    end if
    do b=1,size(block_catalog)
       if(block_catalog(b)%root_domain==glo_id(rank+1,d)) exit
    end do
    if(b>size(block_catalog)) call die('alias without parent block')
    owner=block_catalog(b)%owner
  end function
  function primary(d,id,component) result(ref)
    integer,intent(in)::d,id,component
    integer::ref(5),empty(5,0),key(5)
    if(map(d)%primary(3,1,id+1)==0) then
       key=resident_operand_key(d,id)
       map(d)%primary(:,:,id+1)=program_add(program,target_owner(d,id),OP_RESIDENT,empty,[real(dp)::],7,key)
    end if
    ref=map(d)%primary(:,component,id+1)
  end function
  function operand(d,id) result(ref)
    integer,intent(in)::d,id
    integer::ref(5)
    select case(compiling_field)
    case(2,3)
       ref=primary(d,id,compiling_field+4)
    case(4)
       ref=map(d)%divu%ref(:,id+1)
    case default
       call die('scalar input selector')
    end select
  end function
  subroutine compile_flux(ids,southwest)
    integer,intent(in)::ids(0:7)
    logical,intent(in)::southwest
    integer::d,t,edges(3),nodes(4),refs(5,4),result(5,3),code
    real(dp)::g(6)
    d=compiling_domain
    edges=3*ids(0)+[1,2,3];nodes=ids([0,2,5,1]);code=OP_GRADIENT
    if(southwest) then
       edges=3*ids([4,7,3])+[1,2,3];nodes=ids([0,4,7,3]);code=OP_GRADIENT_SW
    end if
    if(coverage_pass) then
       map(d)%flux%ref(3,edges)=1;return
    end if
    g(1:3)=grid(d)%len%elts(edges);g(4:6)=grid(d)%pedlen%elts(edges)
    if(compiling_field==1) then
       do t=1,3
          refs(:,t)=primary(d,(edges(t)-1)/3,t+2)
       end do
       result=program_add(program,target_owner(d,ids(0)),OP_VELOCITY_FLUX,refs(:,1:3),g(4:6),3)
    else
       do t=1,4
          refs(:,t)=operand(d,nodes(t))
       end do
       result=program_add(program,target_owner(d,ids(0)),code,refs,g,3)
    end if
    map(d)%flux%ref(:,edges)=result
  end subroutine
  subroutine compile_divergence(dom,i,j,k,offs,dims)
    type(Domain),intent(inout)::dom
    integer,intent(in)::i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer::d,id,edges(6),r(5,1)
    d=compiling_domain;id=idx(i,j,offs,dims)
    if(coverage_pass) then
       map(d)%divu%ref(3,id+1)=1;map(d)%mass%ref(3,id+1)=1
       return
    end if
    edges=[3*id+1,3*idx(i-1,j,offs,dims)+1,3*idx(i-1,j-1,offs,dims)+2,3*id+2, &
         3*id+3,3*idx(i,j-1,offs,dims)+3]
    if(compiling_field==5) then
       r(:,1)=zero_definition
       if(dom%mask_n%elts(id+1)>=TRSK) r=program_add(program,target_owner(d,id), &
            OP_MASS_DIVERGENCE,map(d)%flux%ref(:,edges),[dom%areas%elts(id+1)%hex_inv],1)
    else
       r=program_add(program,target_owner(d,id),OP_DIVERGENCE,map(d)%flux%ref(:,edges), &
            [dom%areas%elts(id+1)%hex_inv],1)
    end if
    select case(compiling_field)
    case(1)
       map(d)%divu%ref(:,id+1)=r(:,1)
    case(2)
       map(d)%mass%ref(:,id+1)=r(:,1)
    case(3)
       map(d)%temp%ref(:,id+1)=r(:,1)
    case(4)
       map(d)%lapdiv%ref(:,id+1)=r(:,1)
    case(5)
       map(d)%tendency%ref(:,id+1)=r(:,1)
    end select
  end subroutine
  function rhs(d,ids) result(refs)
    integer,intent(in)::d,ids(:)
    integer::refs(5,size(ids))
    select case(compiling_field)
    case(1)
       refs=map(d)%divu%ref(:,ids)
    case(2)
       refs=map(d)%mass%ref(:,ids)
    case(3)
       refs=map(d)%temp%ref(:,ids)
    case(4)
       refs=map(d)%lapdiv%ref(:,ids)
    case(5)
       refs=map(d)%tendency%ref(:,ids)
    end select
  end function
  subroutine compile_restrictions(d,l)
    integer,intent(in)::d,l
    integer::s,e,t,n,refs(5,84),result(5,3),small(36),seed(5)
    real(dp)::g(60)
    do s=1,size(native_mass(d)%level(l)%restriction)
       associate(p=>native_mass(d)%level(l)%restriction(s))
       small=reshape(p%small,[36])
       do e=1,3
          if(.not.p%edge(e)) cycle
          seed=map(d)%flux%ref(:,small(9*(e-1)+1))
          do t=1,84
             refs(:,t)=seed
          end do
          do n=e,e+1
             refs(:,9*(n-1)+1:9*n)=map(d)%flux%ref(:,small(9*(n-1)+1:9*n))
          end do
          refs(:,36+6*(e-1)+1:36+6*e)=map(d)%flux%ref(:,p%partial_flux(:,e))
          refs(:,54+4*(e-1)+1:54+4*e)=rhs(d,p%partial_node(:,e))
          refs(:,66+6*(e-1)+1:66+6*e)=rhs(d,p%coarse_node(:,e))
          g=[reshape(p%weight,[24]),reshape(p%area,[6]),reshape(p%overlap,[12]),reshape(p%coarse,[15]), &
               0.0_dp,0.0_dp,0.0_dp]
          g(57+e)=1.0_dp
          result=program_add(program,target_owner(d,p%target-1),OP_RESTRICTION,refs,g,3)
          map(d)%flux%ref(:,3*(p%target-1)+e)=result(:,e)
       end do
       end associate
    end do
  end subroutine
  function boundary_ref(d,i,pos) result(ref)
    integer,intent(in)::d,i,pos
    integer::ref(5)
    if(pos==AT_EDGE) then
       ref=map(d)%flux%ref(:,i)
    else
       ref=reshape(rhs(d,[i]),[5])
    end if
  end function
  subroutine set_boundary_ref(d,i,pos,ref,orientation)
    integer,intent(in)::d,i,pos,ref(5),orientation
    integer::value(5)
    value=ref;value(4)=orientation*value(4)
    if(coverage_pass) value(3)=1
    if(pos==AT_EDGE) then
       map(d)%flux%ref(:,i)=value
       return
    end if
    select case(compiling_field)
    case(1)
       map(d)%divu%ref(:,i)=value
    case(2)
       map(d)%mass%ref(:,i)=value
    case(3)
       map(d)%temp%ref(:,i)=value
    case(4)
       map(d)%lapdiv%ref(:,i)=value
    case(5)
       map(d)%tendency%ref(:,i)=value
    end select
  end subroutine
  subroutine boundary(pos,l)
    ! Exchange definitions, not values. Snapshot remote sends before ordered
    ! same-rank copies; publish remote definitions after those copies, exactly
    ! matching update_bdry (including its all-level same-rank copy contract).
    integer,intent(in)::pos,l
    integer::sc(n_process),sd(n_process),rc(n_process),rd(n_process), &
         r,ds,dd,g,i,id,n,mult,a,b,ns,nr,pass,ierr,orientation
    integer,allocatable::send(:,:),recv(:,:),target(:,:)
    mult=1
    if(pos==AT_EDGE) mult=3
    do pass=1,2
       ns=0;nr=0;sc=0;rc=0
       do r=1,n_process
          sd(r)=ns;rd(r)=nr
          if(r==rank+1) cycle
          do ds=1,size(grid)
             do dd=1,n_domain(r)
                g=glo_id(r,dd)+1
                do i=1,grid(ds)%pack(pos,g)%length
                   id=grid(ds)%pack(pos,g)%elts(i)
                   if(grid(ds)%level%elts(id/mult+1)/=l) cycle
                   ns=ns+1
                   if(pass==2) send(:,ns)=boundary_ref(ds,id+1,pos)
                end do
             end do
          end do
          do ds=1,n_domain(r)
             g=glo_id(r,ds)+1
             do dd=1,size(grid)
                do i=1,grid(dd)%unpk(pos,g)%length
                   id=grid(dd)%unpk(pos,g)%elts(i)
                   if(grid(dd)%level%elts(abs(id)/mult+1)/=l) cycle
                   nr=nr+1
                   if(pass==2) target(:,nr)=[dd,abs(id)+1,merge(-1,1,id<0.and.pos==AT_EDGE)]
                end do
             end do
          end do
          sc(r)=ns-sd(r);rc(r)=nr-rd(r)
       end do
       if(pass==1) allocate(send(5,ns),recv(5,nr),target(3,nr))
    end do
    call MPI_Alltoallv(send,5*sc,5*sd,MPI_INTEGER,recv,5*rc,5*rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('boundary definitions')
    do ds=1,size(grid)
       do dd=1,size(grid)
          g=glo_id(rank+1,dd)+1
          do i=1,grid(ds)%pack(pos,g)%length
             a=grid(ds)%pack(pos,g)%elts(i)+1
             id=grid(dd)%unpk(pos,glo_id(rank+1,ds)+1)%elts(i);b=abs(id)+1
             orientation=merge(-1,1,id<0.and.pos==AT_EDGE)
             call set_boundary_ref(dd,b,pos,boundary_ref(ds,a,pos),orientation)
          end do
       end do
    end do
    do n=1,nr
       dd=target(1,n);b=target(2,n)
       call set_boundary_ref(dd,b,pos,recv(:,n),target(3,n))
    end do
  end subroutine
  subroutine seed_untouched_zeros
    ! init_Field and patch extension explicitly initialize scratch to zero.
    ! Preserve that definition only for slots outside EVERY write destination
    ! in this topology. Written slots with missing inputs remain undefined.
    integer::d,l,j,i,s,e,empty(5,0),zero(5,1)
    coverage_pass=.true.;compiling_field=2
    do l=level_end,level_start,-1
       do d=1,size(grid)
          compiling_domain=d
          do j=1,grid(d)%lev(l)%length
             call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(j),compile_flux)
             call apply_onescale_to_patch(compile_divergence,grid(d),grid(d)%lev(l)%elts(j),z_null,0,1)
          end do
          do s=1,size(native_mass(d)%level(l)%restriction)
             do e=1,3
                if(native_mass(d)%level(l)%restriction(s)%edge(e)) &
                     map(d)%flux%ref(3,3*(native_mass(d)%level(l)%restriction(s)%target-1)+e)=1
             end do
          end do
       end do
       call boundary(AT_EDGE,l)
       call boundary(AT_NODE,l)
    end do
    coverage_pass=.false.
    zero=program_add(program,block_catalog(1)%owner,OP_ZERO,empty,[real(dp)::],1)
    zero_definition=zero(:,1)
    allocate(initial_flux(size(grid)))
    do d=1,size(grid)
       do i=1,size(map(d)%flux%ref,2)
          if(map(d)%flux%ref(3,i)==0) then
             map(d)%flux%ref(:,i)=zero(:,1)
          else
             map(d)%flux%ref(:,i)=0
          end if
       end do
       initial_flux(d)%ref=map(d)%flux%ref
       do i=1,size(map(d)%mass%ref,2)
          if(map(d)%mass%ref(3,i)==0) then
             map(d)%mass%ref(:,i)=zero(:,1)
          else
             map(d)%mass%ref(:,i)=0
          end if
          if(map(d)%divu%ref(3,i)==0) then
             map(d)%divu%ref(:,i)=zero(:,1)
          else
             map(d)%divu%ref(:,i)=0
          end if
       end do
       map(d)%temp%ref=map(d)%mass%ref;map(d)%lapdiv%ref=map(d)%mass%ref
       map(d)%initial_rhs%ref=map(d)%divu%ref
    end do
  end subroutine
  subroutine compile_plan
    integer::d,l,j,n,f,i,pass,cursor,ref(5),ff,t
    integer,parameter::field_order(4)=[1,2,4,3]
    integer,allocatable::refs(:,:)
    if(allocated(map)) deallocate(map,destination,initial_flux)
    program=Owner_Program()
    allocate(map(size(grid)))
    do d=1,size(grid)
       n=grid(d)%node%length
       allocate(map(d)%primary(5,7,n),map(d)%flux%ref(5,3*n),map(d)%divu%ref(5,n), &
            map(d)%mass%ref(5,n),map(d)%temp%ref(5,n),map(d)%lapdiv%ref(5,n))
       allocate(map(d)%scalar_slot(6,2,2,n),map(d)%scalar_ref(5,6,2,2,n))
       map(d)%scalar_slot=0;map(d)%scalar_ref=0
       map(d)%primary=0;map(d)%flux%ref=0;map(d)%divu%ref=0;map(d)%mass%ref=0;map(d)%temp%ref=0;map(d)%lapdiv%ref=0
    end do
    call seed_untouched_zeros
    do ff=1,4
       f=field_order(ff)
       ! Temperature has separate flux scratch; do not inherit mass definitions.
       if(f==3) then
          do d=1,size(map)
             map(d)%flux%ref=initial_flux(d)%ref
          end do
       end if
       if(f==1.and.Laplace_divu==0) cycle
       if((f==2.or.f==3).and.Laplace_sclr/=2) cycle
       if(f==4.and.Laplace_divu/=2) cycle
       compiling_field=f
       do l=level_end,level_start,-1
          compiling_level=l
          do d=1,size(grid)
             compiling_domain=d
             do j=1,grid(d)%lev(l)%length
                call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(j),compile_flux)
             end do
             if(l<level_end) call compile_restrictions(d,l)
          end do
          call boundary(AT_EDGE,l)
          do d=1,size(grid)
             compiling_domain=d
             do j=1,grid(d)%lev(l)%length
                call apply_onescale_to_patch(compile_divergence,grid(d),grid(d)%lev(l)%elts(j),z_null,0,1)
             end do
          end do
          select case(f)
          case(2)
             call boundary(AT_NODE,l)
          case(3)
             call boundary(AT_NODE,l)
          case(4)
             call boundary(AT_NODE,l)
          end select
       end do
    end do
    do l=level_end,level_start,-1
       do d=1,size(grid)
          compiling_domain=d
          do j=1,grid(d)%lev(l)%length
             call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(j),compile_scalar_flux)
          end do
       end do
    end do
    if(mass_owner_active) call compile_mass_program
    call program_compile(program)
    do pass=1,2
       cursor=0
       do f=1,4
          do d=1,size(grid)
             do i=1,grid(d)%node%length
                select case(f)
                case(1)
                   ref=map(d)%divu%ref(:,i)
                case(2)
                   ref=map(d)%mass%ref(:,i)
                case(3)
                   ref=map(d)%temp%ref(:,i)
                case(4)
                   ref=map(d)%lapdiv%ref(:,i)
                end select
                if(ref(3)==0) cycle
                cursor=cursor+1
                if(pass==2) then
                   refs(:,cursor)=ref;destination(:,cursor)=[f,d,i]
                end if
             end do
          end do
       end do
       if(pass==1) allocate(refs(5,cursor),destination(3,cursor))
    end do
    call program_read_compile(program,delivery,refs)
    deallocate(refs)
    do pass=1,2
       cursor=0
       do d=1,size(grid)
          do i=1,grid(d)%node%length
             do l=1,2
                do f=1,2
                   do t=1,6
                      if(map(d)%scalar_ref(3,t,f,l,i)==0) cycle
                      cursor=cursor+1
                      if(pass==2) then
                         map(d)%scalar_slot(t,f,l,i)=cursor
                         refs(:,cursor)=map(d)%scalar_ref(:,t,f,l,i)
                      end if
                   end do
                end do
             end do
          end do
       end do
       if(pass==1) allocate(refs(5,cursor))
    end do
    call program_read_compile(program,scalar_delivery,refs)
    if(mass_owner_active.and.(.not.block_rk_enabled.or.mass_owner_oracle)) call compile_mass_delivery
    if(block_rk_enabled) call compile_owned_mass
    do d=1,size(map)
       deallocate(map(d)%primary,map(d)%scalar_ref)
    end do
    generation=inverse_operand_generation()
  end subroutine
  subroutine compile_scalar_flux(ids,southwest)
    integer,intent(in)::ids(0:7)
    logical,intent(in)::southwest
    integer::d,t,v,sw,nodes(4),edges(3),refs(5,11),part(5,2)
    real(dp)::g(9)
    d=compiling_domain;sw=merge(2,1,southwest)
    if(map(d)%scalar_ref(3,1,1,sw,ids(0)+1)/=0) return
    nodes=ids([0,2,5,1]);edges=3*ids(0)+[1,2,3]
    if(southwest) then
       nodes=ids([0,4,7,3]);edges=3*ids([4,7,3])+[1,2,3]
    end if
    g(1:3)=grid(d)%len%elts(edges);g(4:6)=grid(d)%pedlen%elts(edges)
    if(southwest) g(1:3)=-g(1:3)
    g(7)=Laplace_sclr
    do v=1,2
       g(8)=v
       do t=1,4
          refs(:,t)=primary(d,nodes(t),v)
          select case(Laplace_sclr)
          case(0,1)
             refs(:,t+7)=primary(d,nodes(t),v+5)
          case(2)
             if(v==1) refs(:,t+7)=map(d)%mass%ref(:,nodes(t)+1)
             if(v==2) refs(:,t+7)=map(d)%temp%ref(:,nodes(t)+1)
          case default
             call die('unsupported scalar diffusion order')
          end select
       end do
       do t=1,3
          refs(:,t+4)=primary(d,(edges(t)-1)/3,t+2)
       end do
       do t=1,3
          g(9)=t
          part=program_add(program,target_owner(d,ids(0)),OP_SCALAR_FLUX,refs(:,[1,t+1,t+4,8,t+8]),g,2)
          map(d)%scalar_ref(:,t,v,sw,ids(0)+1)=part(:,1)
          map(d)%scalar_ref(:,t+3,v,sw,ids(0)+1)=part(:,2)
       end do
    end do
  end subroutine
  subroutine owner_scalar_flux(d,id,v,k,southwest,physics,flux,defined)
    integer,intent(in)::d,id,v,k
    logical,intent(in)::southwest
    real(dp),intent(out)::physics(3),flux(3)
    integer::a,t,sw,vv
    logical,intent(out)::defined(3)
    if(.not.diffusion_active.or.program%layer/=k) call die('scalar flux epoch')
    sw=merge(2,1,southwest);vv=v-S_MASS+1
    physics=ieee_value(0.0_dp,ieee_quiet_nan);flux=physics;defined=.false.
    do t=1,3
       a=map(d)%scalar_slot(t,vv,sw,id+1)
       if(a<=0) cycle
       physics(t)=scalar_delivery%value(a)
       a=map(d)%scalar_slot(t+3,vv,sw,id+1)
       if(a<=0) call die('scalar physics without matching advective producer')
       flux(t)=scalar_delivery%value(a);defined(t)=.true.
    end do
  end subroutine
  function owner_mass_source(d,edge,l) result(key)
    integer,intent(in)::d,edge,l
    integer::key(5)
    if(.not.mass_owner_active) call die('mass source requested without owner execution')
    key=map(d)%mass_flux(l)%ref(:,edge)
    if(key(3)==0) call die('missing live mass flux producer')
  end function
  function owner_diffusion_source(d,id,field) result(key)
    integer,intent(in)::d,id,field
    integer::key(5)
    select case(field)
    case(1)
       key=map(d)%divu%ref(:,id)
    case(2)
       key=map(d)%mass%ref(:,id)
    case(3)
       key=map(d)%temp%ref(:,id)
    case(4)
       key=map(d)%lapdiv%ref(:,id)
    case default
       call die('external diffusion field')
    end select
    if(key(3)==0) call die('missing live diffusion producer')
  end function
  function owner_diffusion_value(key,k) result(value)
    integer,intent(in)::key(5),k
    real(dp)::value
    value=program_value(program,key,k)
  end function
  subroutine compile_mass_flux(ids,southwest)
    integer,intent(in)::ids(0:7)
    logical,intent(in)::southwest
    integer::edges(3),d
    d=compiling_domain
    edges=3*ids(0)+[1,2,3]
    if(southwest) edges=3*ids([4,7,3])+[1,2,3]
    map(d)%flux%ref(:,edges)=map(d)%scalar_ref(:,4:6,1,merge(2,1,southwest),ids(0)+1)
  end subroutine
  subroutine compile_mass_program
    integer::d,l,j
    compiling_field=5
    do d=1,size(map)
       map(d)%flux%ref=initial_flux(d)%ref
       map(d)%tendency%ref=map(d)%initial_rhs%ref
       allocate(map(d)%mass_flux(level_start:level_end),map(d)%mass_rhs(level_start:level_end))
    end do
    do l=level_end,level_start,-1
       compiling_level=l
       if(l<level_end) call boundary(AT_NODE,l+1)
       do d=1,size(grid)
          compiling_domain=d
          do j=1,grid(d)%lev(l)%length
             call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(j),compile_mass_flux)
          end do
          call compile_restrictions(d,l)
       end do
       if(level_start/=level_end) call boundary(AT_EDGE,l)
       call boundary(AT_EDGE,l)
       do d=1,size(grid)
          compiling_domain=d
          do j=1,grid(d)%lev(l)%length
             call apply_onescale_to_patch(compile_divergence,grid(d),grid(d)%lev(l)%elts(j),z_null,0,1)
          end do
          map(d)%mass_flux(l)%ref=map(d)%flux%ref
          map(d)%mass_rhs(l)%ref=map(d)%tendency%ref
       end do
    end do
  end subroutine
  subroutine compile_mass_delivery
    integer::l,d,i,n,pass,f,ref(5)
    integer,allocatable::refs(:,:)
    if(allocated(mass_output)) deallocate(mass_output)
    allocate(mass_output(level_start:level_end))
    do l=level_start,level_end
       do pass=1,2
          n=0
          do d=1,size(map)
             do f=0,1
                do i=1,merge(size(map(d)%mass_flux(l)%ref,2),size(map(d)%mass_rhs(l)%ref,2),f==0)
                   if(f==0) then
                      ref=map(d)%mass_flux(l)%ref(:,i)
                   else
                      ref=map(d)%mass_rhs(l)%ref(:,i)
                   end if
                   if(ref(3)==0) cycle
                   n=n+1
                   if(pass==2) then
                      refs(:,n)=ref;mass_output(l)%target(:,n)=[d,i,f]
                   end if
                end do
             end do
          end do
          if(pass==1) allocate(refs(5,n),mass_output(l)%target(3,n))
       end do
       call program_read_compile(program,mass_output(l)%input,refs)
       deallocate(refs)
    end do
  end subroutine
  subroutine compile_owned_mass
    integer::d,id,key(5),ref(5),b,p,n,peer,ierr,t,pass,i
    integer,allocatable::sc(:),rc(:),sd(:),rd(:),cursor(:),send(:,:),recv(:,:),aliases(:,:),refs(:,:)
    if(allocated(owned_mass)) deallocate(owned_mass,owned_mass_value,owned_mass_ready)
    allocate(sc(n_process),rc(n_process),sd(n_process),rd(n_process),cursor(n_process))
    sc=0
    do pass=1,2
       cursor=0
       do d=1,size(map)
          do id=0,grid(d)%node%length-1
             key=resident_operand_key(d,id)
             if(key(2)==0) cycle
             peer=key(1)+1
             if(pass==1) then
                sc(peer)=sc(peer)+1
             else
                p=sd(peer)+cursor(peer)+1;cursor(peer)=cursor(peer)+1
                ref=map(d)%mass_rhs(level_start)%ref(:,id+1)
                if(ref(3)==0) call die('missing retained mass definition')
                send(:,p)=[key(2),key(4),key(5),ref]
             end if
          end do
       end do
       if(pass==1) then
          call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
          if(ierr/=MPI_SUCCESS) call die('retained mass counts')
          sd(1)=0;rd(1)=0
          do peer=2,n_process
             sd(peer)=sd(peer-1)+sc(peer-1);rd(peer)=rd(peer-1)+rc(peer-1)
          end do
          allocate(send(8,sum(sc)),recv(8,sum(rc)))
       end if
    end do
    call MPI_Alltoallv(send,8*sc,8*sd,MPI_INTEGER,recv,8*rc,8*rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('retained mass keys')
    n=size(recv,2)
    allocate(owned_mass(size(block_catalog)),owned_mass_value(n,zlevels),owned_mass_ready(zlevels))
    owned_mass_ready=.false.
    do b=1,size(block_catalog)
       if(block_catalog(b)%owner/=rank) cycle
       allocate(owned_mass(b)%slot(PATCH_SIZE**2,local_block_patch_count(b)))
       owned_mass(b)%slot=0
    end do
    do t=1,n
       b=recv(1,t);p=recv(2,t);id=recv(3,t)+1
       if(.not.allocated(owned_mass(b)%slot)) call die('retained mass wrong owner')
       if(owned_mass(b)%slot(id,p)/=0) call die('retained mass duplicate')
       owned_mass(b)%slot(id,p)=t
    end do
    do b=1,size(owned_mass)
       if(.not.allocated(owned_mass(b)%slot)) cycle
       if(any(owned_mass(b)%slot==0)) call die('retained mass missing node')
    end do
    call program_read_compile(program,owned_mass_read,recv(4:8,:))
    call rk_alias_nodes(aliases)
    allocate(refs(5,size(aliases,2)))
    do i=1,size(aliases,2)
       d=aliases(1,i);id=aliases(2,i)
       ref=map(d)%mass_rhs(level_start)%ref(:,id+1)
       ! Only untouched slots may have no live producer. The compiler's
       ! initial_rhs mask distinguishes them from live missing definitions.
       if(ref(3)==0) then
          if(map(d)%initial_rhs%ref(3,id+1)==0) call die('missing RK alias mass producer')
          ref=map(d)%initial_rhs%ref(:,id+1)
       end if
       refs(:,i)=ref
    end do
    call program_read_compile(program,alias_mass_read,refs)
  end subroutine
  subroutine retain_owned_mass(k)
    integer,intent(in)::k
    if(k==1) owned_mass_ready=.false.
    call program_read(program,owned_mass_read)
    owned_mass_value(:,k)=owned_mass_read%value;owned_mass_ready(k)=.true.
    call program_read(program,alias_mass_read)
    call set_rk_mass_alias_layer(k,alias_mass_read%value)
  end subroutine
  real(dp) function owned_mass_node(b,p,node,k) result(x)
    integer,intent(in)::b,p,node,k
    integer::t
    if(.not.block_rk_enabled.or..not.all(owned_mass_ready)) call die('retained mass not ready')
    if(generation/=inverse_operand_generation().or.k<1.or.k>zlevels) call die('retained mass epoch/layer')
    if(.not.allocated(owned_mass(b)%slot)) call die('retained mass non-owner read')
    t=owned_mass(b)%slot(node+1,p)
    if(t==0) call die('retained mass undefined node')
    x=owned_mass_value(t,k)
  end function

  subroutine install_owner_mass(k,l)
    integer,intent(in)::k,l
    integer::s,d,i,f
    real(dp)::actual,expected
    if(.not.mass_owner_active.or.program%layer/=k) call die('mass producer epoch')
    if(block_rk_enabled.and.l==level_start) call retain_owned_mass(k)
    if(block_rk_enabled.and..not.mass_owner_oracle) return
    associate(q=>mass_output(l))
    call program_read(program,q%input)
    if(mass_owner_oracle) then
       do s=1,size(q%target,2)
          d=q%target(1,s);i=q%target(2,s);f=q%target(3,s);actual=q%input%value(s)
          if(f==0) then
             expected=native_mass(d)%flux(i)
          else
             expected=native_mass(d)%tendency(i,k)
          end if
          ! A level has not yet defined coarse RHS; those entries are absent.
          if(transfer(actual,0_int64)/=transfer(expected,0_int64)) then
             write(*,*) 'owner mass mismatch domain/id/k/l/field/actual/expected ',d,i,k,l,f,actual,expected
             call die('independent mass program mismatch')
          end if
          mass_checked=mass_checked+1
       end do
    end if
    do d=1,size(native_mass)
       native_mass(d)%flux=ieee_value(0.0_dp,ieee_quiet_nan)
       ! The RK compatibility return spans untouched scaffold slots too.
       ! They retain zero_float's initial +0, unlike undefined *operands* in
       ! the numerical graph, which remain absent and are rejected if read.
       native_mass(d)%tendency(:,k)=0.0_dp
    end do
    do s=1,size(q%target,2)
       d=q%target(1,s);i=q%target(2,s);f=q%target(3,s)
       if(f==0) then
          native_mass(d)%flux(i)=q%input%value(s)
       else
          native_mass(d)%tendency(i,k)=q%input%value(s)
       end if
    end do
    end associate
    if(k==1.and.l==level_start) write(*,'(a,i0,a,i0)') 'Block mass owner rank=',rank,' checked=',mass_checked
  end subroutine
  subroutine begin_block_diffusion(k)
    integer,intent(in)::k
    diffusion_active=.false.;diffusion_oracle=.false.
#ifdef BLOCK_DIFFUSION
    diffusion_active=columns_active;diffusion_oracle=columns_oracle
#ifdef BLOCK_DIFFUSION_ORACLE
    diffusion_oracle=.true.
#endif
#endif
    mass_owner_active=.false.;mass_owner_oracle=diffusion_oracle
#ifdef BLOCK_MASS_OWNER
    mass_owner_active=diffusion_active
#ifdef BLOCK_MASS_OWNER_ORACLE
    mass_owner_oracle=.true.
#endif
#endif
    if(.not.diffusion_active) return
    if(generation/=inverse_operand_generation()) call compile_plan
    call program_execute(program,k,[nu_scale(S_MASS,k),nu_scale(S_TEMP,k)])
    call program_read(program,scalar_delivery)
    if(.not.owner_tendency_without_reference()) call program_read(program,delivery)
  end subroutine
  subroutine install_block_diffusion(k)
    integer,intent(in)::k
    integer::s,f,d,i
    real(dp)::old,value
    if(.not.diffusion_active) return
    if(program%layer/=k) call die('wrong layer install')
    if(owner_tendency_without_reference()) return
    do s=1,size(destination,2)
       f=destination(1,s);d=destination(2,s);i=destination(3,s);value=delivery%value(s)
       if(diffusion_oracle) then
          select case(f)
          case(1)
             old=grid(d)%divu%elts(i)
          case(2,3)
             old=Laplacian_scalar(S_MASS+f-2)%data(d)%elts(i)
          case(4)
             old=Laplacian_vector(S_DIVU)%data(d)%elts(i)
          end select
          if(transfer(old,0_int64)/=transfer(value,0_int64)) then
             write(*,*) 'diffusion oracle field/domain/node/layer/actual/expected ',f,d,i,k,value,old
             call die('independent diffusion mismatch')
          end if
          checked=checked+1
       end if
       select case(f)
       case(1)
          grid(d)%divu%elts(i)=value
       case(2,3)
          Laplacian_scalar(S_MASS+f-2)%data(d)%elts(i)=value
       case(4)
          Laplacian_vector(S_DIVU)%data(d)%elts(i)=value
       end select
    end do
    if(.not.diffusion_oracle) then
    do d=1,size(map)
       i=grid(d)%node%length
       if(Laplace_divu/=0) then
          where(map(d)%divu%ref(3,:)==0) grid(d)%divu%elts(1:i)=ieee_value(0.0_dp,ieee_quiet_nan)
       end if
       if(Laplace_sclr==2) then
          where(map(d)%mass%ref(3,:)==0) &
               Laplacian_scalar(S_MASS)%data(d)%elts(1:i)=ieee_value(0.0_dp,ieee_quiet_nan)
          where(map(d)%temp%ref(3,:)==0) &
               Laplacian_scalar(S_TEMP)%data(d)%elts(1:i)=ieee_value(0.0_dp,ieee_quiet_nan)
       end if
       if(Laplace_divu==2) then
          where(map(d)%lapdiv%ref(3,:)==0) &
               Laplacian_vector(S_DIVU)%data(d)%elts(1:i)=ieee_value(0.0_dp,ieee_quiet_nan)
       end if
    end do
    end if
    if(k==1) write(*,'(a,i0,a,2(i0,1x))') 'Block diffusion rank=',rank,' delivered/checked=',size(destination,2),checked
  end subroutine
end module
