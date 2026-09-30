module parallel_block_velocity_owner_mod
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run,block_catalog,glo_id
  use shared_mod, only: level_start,level_end,N_BDRY,z_null,Laplace_divu,Laplace_rotu,S_DIVU,S_ROTU, &
       eps,radius,C_visc,min_level,dx_avg,penta_node,IJMINUS,IPLUSJMINUS,IMINUSJPLUS,IJPLUS
  use domain_mod, only: Domain,grid,idx
  use domain_ops_mod, only: apply_onescale_to_patch,apply_to_penta_d
  use patch_mod, only: PATCH_SIZE
  use geom_mod, only: geodesic
  use utils_mod, only: nu_scale
  use parallel_block_resident_operands_mod, only: resident_operand_key
  use parallel_block_inverse_mod, only: inverse_operand_generation
  use parallel_block_primitive_walk_mod, only: walk_primitive_patch
  use parallel_block_primitive_registry_mod
  use parallel_block_primitives_mod, only: primitive_owner_value
  use parallel_block_pentagons_mod, only: pentagon_owner_value
  use parallel_block_diffusion_mod, only: mass_owner_active,mass_owner_oracle,owner_mass_source, &
       owner_diffusion_source,owner_diffusion_value
  use parallel_block_velocity_mod, only: native_velocity,VELOCITY_DIRECT,VELOCITY_RESTRICT
  use parallel_block_program_mod
  use parallel_block_definition_exchange_mod
  implicit none
  private
  public :: velocity_owner_enabled,velocity_owner_check,execute_owner_velocity,install_owner_velocity
  public :: owner_velocity_source_key,owner_velocity_source_value
  type :: Source_Map
     integer,allocatable :: source(:,:,:),rot(:,:,:)
     logical,allocatable :: written(:,:),rot_written(:,:)
  end type
  type :: Source_Plan
     type(Owner_Program) :: code
     type(Owner_Read) :: returned,rotation
     integer,allocatable :: target(:,:),rot_target(:,:),rot_key(:,:)
  end type
  type(Source_Map),allocatable :: map(:)
  type(Owner_Definitions),allocatable :: curl(:)
  type(Source_Plan),allocatable :: plan(:)
  integer(int64) :: generation=-1_int64,checked=0_int64
  integer :: compiling_level=0,compiling_domain=0,compiling_layer=0,zero(5)
  logical :: coverage=.false.
contains
  subroutine die(message)
    character(*),intent(in)::message
    write(*,'(a,i0,2a)') 'Rank ',rank,': owner velocity: ',message
    call abort_run
    error stop 'owner velocity abort returned'
  end subroutine
  logical function velocity_owner_enabled()
    velocity_owner_enabled=.false.
#ifdef BLOCK_VELOCITY_OWNER
    velocity_owner_enabled=mass_owner_active
#endif
  end function
  logical function velocity_owner_check()
    velocity_owner_check=mass_owner_oracle
#ifdef BLOCK_VELOCITY_OWNER_ORACLE
    velocity_owner_check=.true.
#endif
  end function
  integer function owner(d,id) result(r)
    integer,intent(in)::d,id
    integer::key(5),b
    key=resident_operand_key(d,id)
    r=key(1)
    if(key(2)>0) return
    do b=1,size(block_catalog)
       if(block_catalog(b)%root_domain==glo_id(rank+1,d)) exit
    end do
    if(b>size(block_catalog)) call die('alias parent block')
    r=block_catalog(b)%owner
  end function
  function result_key(ref) result(key)
    integer,intent(in)::ref(5)
    integer::key(7)
    key=0
    if(ref(3)>0) key=[3,ref(1:4),compiling_level,ref(5)]
  end function
  function reference(key) result(ref)
    integer,intent(in)::key(7)
    integer::ref(5),empty(5,0),r(5,1),code
    ref=0
    if(key(1)==0) return
    if(key(1)==4) then
       ref=zero;ref(4)=sign(1,key(5))*ref(4);return
    end if
    if(key(1)==3.and.key(6)==compiling_level) then
       ref=[key(2:5),key(7)];return
    end if
    code=100+key(1)
    r=program_add(plan(compiling_level)%code,key(2),code,empty,[real(dp)::],1,key(2:6))
    ref=r(:,1)
  end function
  function external(code,key,k) result(value)
    integer,intent(in)::code,key(5),k
    real(dp)::value
    select case(code)
    case(101)
       value=primitive_owner_value(key,k)
    case(102)
       value=pentagon_owner_value(key,k)
    case(103)
       value=program_value(plan(key(5))%code,key,k)
    case(105)
       value=owner_diffusion_value(key,k)
    case default
       call die('external producer kind')
    end select
  end function
  function parent_reference(key) result(ref)
    integer,intent(in)::key(5)
    integer::ref(5),r(5,1),empty(5,0)
    r=program_add(plan(compiling_level)%code,key(1),105,empty,[real(dp)::],1,key)
    ref=r(:,1)
  end function
  subroutine compile_curl(dom,i,j,k,offs,dims)
    type(Domain),intent(inout)::dom
    integer,intent(in)::i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer::d,id,s,w,e,refs(5,2),r(5,1),nodes(2,3),components(2,3)
    real(dp)::ped
    d=compiling_domain;id=idx(i,j,offs,dims)
    if(coverage) then
       curl(d)%key(1,3*id+1:3*id+3)=99;return
    end if
    s=idx(i,j-1,offs,dims);w=idx(i-1,j,offs,dims)
    nodes(:,1)=[id,s];nodes(:,2)=[id,id];nodes(:,3)=[w,id]
    components(:,1)=[1,2];components(:,2)=[1,2];components(:,3)=[1,2]
    do e=1,3
       ped=dom%pedlen%elts(3*id+e)
       r(:,1)=zero
       if(e/=2.or.ped>eps(radius)) then
          refs(:,1)=reference(primitive_result_source(d,nodes(1,e),components(1,e),.true.,compiling_layer,.false.))
          refs(:,2)=reference(primitive_result_source(d,nodes(2,e),components(2,e),.true.,compiling_layer,.false.))
          r=program_add(plan(compiling_level)%code,owner(d,id),OP_CURL,refs,[ped],1)
       end if
       curl(d)%key(:,3*id+e)=result_key(r(:,1))
    end do
  end subroutine
  subroutine circulate(id,component,edges,mode)
    integer,intent(in)::id,component,edges(3),mode
    integer::d,t,refs(5,3),r(5,1)
    real(dp)::g(5)
    d=compiling_domain
    do t=1,3
       refs(:,t)=reference(curl(d)%key(:,edges(t)))
    end do
    g=[grid(d)%len%elts(edges),grid(d)%triarea%elts(2*id+component),real(mode,dp)]
    r=program_add(plan(compiling_level)%code,owner(d,id),OP_CIRCULATION,refs,g,1)
    map(d)%rot(:,component,id+1)=result_key(r(:,1))
    map(d)%rot_written(component,id+1)=.true.
  end subroutine
  subroutine compile_circulation(ids,southwest)
    integer,intent(in)::ids(0:7)
    logical,intent(in)::southwest
    if(southwest) then
       call circulate(ids(7),1,[3*ids(7)+1,3*ids(3)+3,3*ids(7)+2],1)
       call circulate(ids(7),2,[3*ids(4)+1,3*ids(7)+2,3*ids(7)+3],2)
    else
       call circulate(ids(0),1,[3*ids(0)+1,3*ids(2)+3,3*ids(0)+2],1)
       call circulate(ids(0),2,[3*ids(0)+2,3*ids(0)+3,3*ids(1)+1],2)
       call circulate(ids(4),1,[3*ids(4)+1,3*ids(0)+3,3*ids(4)+2],1)
       call circulate(ids(3),2,[3*ids(0)+1,3*ids(3)+2,3*ids(3)+3],2)
    end if
  end subroutine
  subroutine compile_pentagon(dom,p,c,offs,dims,k)
    type(Domain),intent(inout)::dom
    integer,intent(in)::p,c,offs(N_BDRY+1),dims(2,N_BDRY+1),k
    integer::d,id,s,w,sw,n
    d=compiling_domain
    select case(c)
    case(IJMINUS)
       s=idx(0,-1,offs,dims);sw=idx(-1,-1,offs,dims);w=idx(-1,0,offs,dims)
       call circulate(sw,1,[3*s+3,3*w+1,3*sw+1],3)
       map(d)%rot(:,2,sw+1)=map(d)%rot(:,1,sw+1)
       map(d)%rot_written(2,sw+1)=.true.
    case(IPLUSJMINUS)
       id=idx(PATCH_SIZE,0,offs,dims);s=idx(PATCH_SIZE,-1,offs,dims);sw=idx(PATCH_SIZE-1,-1,offs,dims)
       call circulate(sw,1,[3*id+1,3*sw+1,3*sw+2],4)
       map(d)%rot(:,2,s+1)=map(d)%rot(:,1,sw+1)
       map(d)%rot_written(2,s+1)=.true.
    case(IMINUSJPLUS)
       id=idx(0,PATCH_SIZE,offs,dims);sw=idx(-1,PATCH_SIZE-1,offs,dims);w=idx(-1,PATCH_SIZE,offs,dims)
       call circulate(sw,2,[3*id+3,3*sw+2,3*sw+3],5)
       map(d)%rot(:,1,w+1)=map(d)%rot(:,2,sw+1)
       map(d)%rot_written(1,w+1)=.true.
    case(IJPLUS)
       id=idx(PATCH_SIZE,PATCH_SIZE,offs,dims);n=idx(PATCH_SIZE,PATCH_SIZE+1,offs,dims)
       call circulate(id,1,[3*id+1,3*n+1,3*id+3],5)
       map(d)%rot(:,2,id+1)=map(d)%rot(:,1,id+1)
       map(d)%rot_written(2,id+1)=.true.
    case default
       call die('pentagon corner')
    end select
  end subroutine
  subroutine initialize
    integer::d,n,l,j,i
    if(allocated(plan)) deallocate(plan,map,curl)
    allocate(plan(level_start:level_end),map(size(grid)),curl(size(grid)))
    do d=1,size(grid)
       n=grid(d)%node%length
       allocate(map(d)%source(7,3,n),map(d)%rot(7,2,n),map(d)%written(3,n),map(d)%rot_written(2,n),curl(d)%key(7,3*n))
       map(d)%source=0;map(d)%rot=0;curl(d)%key=0
       do i=1,n
          do j=1,3
             map(d)%source(:,j,i)=[4,0,0,0,1,0,0]
          end do
       end do
    end do
    ! Only slots outside all curl writes/receives inherit initialized zero.
    coverage=.true.
    do l=level_end,level_start,-1
       do d=1,size(grid)
          compiling_domain=d
          do j=1,grid(d)%lev(l)%length
             call apply_onescale_to_patch(compile_curl,grid(d),grid(d)%lev(l)%elts(j),z_null,0,0)
          end do
       end do
       call exchange_edge_definitions(curl,l,.true.)
    end do
    coverage=.false.
    do d=1,size(grid)
       do i=1,size(curl(d)%key,2)
          if(curl(d)%key(1,i)==0) then
             curl(d)%key(:,i)=[4,0,0,0,1,0,0]
          else
             curl(d)%key(:,i)=0
          end if
       end do
    end do
    generation=inverse_operand_generation()
  end subroutine
  subroutine compile_plan(k,l)
    integer,intent(in)::k,l
    integer::d,n,a,s,id,e,t,ip,refs(5,57),p(5,8),r(5,3),one(5,1),empty(5,0),key(5),nkey(7)
    real(dp)::g(9),near
    logical,parameter::flux_needed(3,9)=reshape([ &
         .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
         .true.,.true.,.true.,.true.,.true.,.false.,.false.,.true.,.false., &
         .false.,.true.,.true.,.false.,.false.,.true.,.true.,.false.,.true.],[3,9])
    logical,parameter::pv_needed(3,9)=reshape([ &
         .true.,.true.,.true.,.true.,.true.,.true.,.true.,.true.,.true., &
         .true.,.true.,.true.,.true.,.true.,.false.,.false.,.true.,.false., &
         .false.,.true.,.true.,.false.,.false.,.true.,.true.,.false.,.true.],[3,9])
    compiling_level=l;compiling_layer=k
    one=program_add(plan(l)%code,block_catalog(1)%owner,OP_ZERO,empty,[real(dp)::],1);zero=one(:,1)
    do d=1,size(map)
       map(d)%written=.false.;map(d)%rot=0;map(d)%rot_written=.false.
    end do
    if(Laplace_rotu==2) then
       do d=1,size(grid)
          compiling_domain=d
          do s=1,grid(d)%lev(l)%length
             call apply_onescale_to_patch(compile_curl,grid(d),grid(d)%lev(l)%elts(s),z_null,0,0)
          end do
       end do
       call exchange_edge_definitions(curl,l)
       do d=1,size(grid)
          compiling_domain=d
          do s=1,grid(d)%lev(l)%length
             call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(s),compile_circulation)
          end do
          call apply_to_penta_d(compile_pentagon,grid(d),l,z_null)
       end do
    end if
    do d=1,size(native_velocity)
       compiling_domain=d
       associate(tape=>native_velocity(d)%level(l))
       do a=1,size(tape%action)
          s=tape%operand(a)
          select case(tape%action(a))
          case(VELOCITY_DIRECT)
             associate(stencil=>tape%direct(s))
             id=stencil%node(1)
             do e=1,3
                r(:,e)=zero
             end do
             if(stencil%active) then
                do t=1,57
                   refs(:,t)=zero
                end do
                do n=1,9
                   do e=1,3
                      if(flux_needed(e,n)) then
                         key=owner_mass_source(d,3*(stencil%node(n)-1)+e,l)
                         refs(:,3*(n-1)+e)=parent_reference(key)
                      end if
                      if(pv_needed(e,n)) then
                         nkey=primitive_result_source(d,stencil%node(n)-1,e,.false.,k)
                         refs(:,27+3*(n-1)+e)=reference(nkey)
                      end if
                   end do
                end do
                do t=1,8
                   p(:,t)=zero
                end do
                if(Laplace_divu/=0) then
                   do t=1,4
                      key=owner_diffusion_source(d,stencil%node(t),merge(4,1,Laplace_divu==2))
                      p(:,t)=parent_reference(key)
                   end do
                end if
                if(Laplace_rotu/=0) then
                   p(:,5)=vorticity(d,id,1,k)
                   p(:,6)=vorticity(d,id,2,k)
                   p(:,7)=vorticity(d,stencil%node(7),2,k)
                   p(:,8)=vorticity(d,stencil%node(5),1,k)
                end if
                near=0.0_dp
                if(Laplace_divu==1) then
                   do ip=1,12
                      if(geodesic(grid(d)%node%elts(id),penta_node(ip))<1.5*dx_avg(min_level)) near=1.0_dp
                   end do
                end if
                g=[stencil%length,grid(d)%pedlen%elts(3*(id-1)+1:3*id), &
                     real(Laplace_divu,dp),real(Laplace_rotu,dp),near]
                r=program_add(plan(l)%code,owner(d,id-1),OP_VELOCITY_PHYSICS,p,g,3)
                refs(:,55:57)=r
                r=program_add(plan(l)%code,owner(d,id-1),OP_VELOCITY_SOURCE,refs, &
                     [reshape(stencil%weights,[30]),stencil%length],3)
                if(any(r(3,:)==0)) then
                   write(*,*) 'missing direct source domain/node/level ',d,id,l
                   call die('undefined active source input')
                end if
             end if
             do e=1,3
                map(d)%source(:,e,id)=result_key(r(:,e));map(d)%written(e,id)=.true.
             end do
             end associate
          case(VELOCITY_RESTRICT)
             associate(stencil=>tape%restriction(s))
             do e=1,3
                if(.not.stencil%edge(e)) cycle
                p(:,1)=reference(map(d)%source(:,e,stencil%child))
                p(:,2)=reference(map(d)%source(:,e,stencil%neighbor(e)))
                one=program_add(plan(l)%code,owner(d,stencil%target-1),OP_ADD,p(:,1:2),[real(dp)::],1)
                if(one(3,1)==0) call die('undefined source restriction input')
                map(d)%source(:,e,stencil%target)=result_key(one(:,1));map(d)%written(e,stencil%target)=.true.
             end do
             end associate
          case default
             call die('source action')
          end select
       end do
       end associate
    end do
    call compile_deliveries(l)
  end subroutine
  function vorticity(d,id,e,k) result(ref)
    integer,intent(in)::d,id,e,k
    integer::ref(5),key(7)
    if(Laplace_rotu==2.and.map(d)%rot_written(e,id)) then
       key=map(d)%rot(:,e,id)
    else
       key=primitive_result_source(d,id-1,e,.true.,k,.false.)
    end if
    ref=reference(key)
  end function
  subroutine compile_deliveries(l)
    integer,intent(in)::l
    integer::d,id,e,n,pass,nr,j,nall,key(7)
    integer,allocatable::refs(:,:),rotrefs(:,:)
    do pass=1,2
       n=0;nr=0;nall=0
       do d=1,size(map)
          do id=1,grid(d)%node%length
             do e=1,3
                if(.not.map(d)%written(e,id)) cycle
                n=n+1
                if(pass==2) then
                   refs(:,n)=reference(map(d)%source(:,e,id));plan(l)%target(:,n)=[d,id,e]
                end if
             end do
             do e=1,2
                if(.not.map(d)%rot_written(e,id)) cycle
                nall=nall+1;key=map(d)%rot(:,e,id);j=0
                if(key(1)/=0) then
                   nr=nr+1;j=nr
                   if(pass==2) rotrefs(:,nr)=reference(key)
                end if
                if(pass==2) then
                   plan(l)%rot_target(:,nall)=[d,id,e,j];plan(l)%rot_key(:,nall)=key
                end if
             end do
          end do
       end do
       if(pass==1) allocate(refs(5,n),rotrefs(5,nr),plan(l)%target(3,n), &
            plan(l)%rot_target(4,nall),plan(l)%rot_key(7,nall))
    end do
    call program_compile(plan(l)%code)
    call program_read_compile(plan(l)%code,plan(l)%returned,refs)
    call program_read_compile(plan(l)%code,plan(l)%rotation,rotrefs)
  end subroutine
  subroutine execute_owner_velocity(k,l)
    integer,intent(in)::k,l
    integer::s,d,id,e,j
    real(dp)::coeff(3),saved,expected,actual
    if(.not.velocity_owner_enabled()) call die('execute without owner path')
    if(generation/=inverse_operand_generation()) then
       if(l/=level_end) call die('layout changed during source execution')
       call initialize
    end if
    if(.not.plan(l)%code%compiled) call compile_plan(k,l)
    coeff(1)=nu_scale(S_DIVU,k);coeff(2)=nu_scale(S_ROTU,k);coeff(3)=coeff(1)
    if(Laplace_divu==1) then
       saved=C_visc(S_DIVU,k);C_visc(S_DIVU,k)=0.05_dp
       coeff(3)=nu_scale(S_DIVU,k);C_visc(S_DIVU,k)=saved
    end if
    call program_execute(plan(l)%code,k,coeff,external)
    if(velocity_owner_check()) call program_read(plan(l)%code,plan(l)%returned)
    call program_read(plan(l)%code,plan(l)%rotation)
    do s=1,size(plan(l)%rot_target,2)
       d=plan(l)%rot_target(1,s);id=plan(l)%rot_target(2,s);e=plan(l)%rot_target(3,s);j=plan(l)%rot_target(4,s)
       if(velocity_owner_check().and.j>0) then
          actual=plan(l)%rotation%value(j);expected=grid(d)%vort%elts(2*(id-1)+e)
          if(transfer(actual,0_int64)/=transfer(expected,0_int64)) call die('independent rotational diffusion mismatch')
       end if
       call register_primitive_result(d,id-1,e,.true.,plan(l)%rot_key(:,s),k)
    end do
  end subroutine
  subroutine install_owner_velocity(k,l)
    integer,intent(in)::k,l
    integer::s,d,id,e
    real(dp)::actual,expected
    if(plan(l)%code%layer/=k) call die('source install layer')
    if(.not.velocity_owner_check()) return
    do s=1,size(plan(l)%target,2)
       d=plan(l)%target(1,s);id=plan(l)%target(2,s);e=plan(l)%target(3,s)
       actual=plan(l)%returned%value(s)
       if(velocity_owner_check()) then
          expected=native_velocity(d)%source(e,id)
          if(transfer(actual,0_int64)/=transfer(expected,0_int64)) then
             write(*,*) 'source mismatch d/id/edge/k/l/actual/expected ',d,id,e,k,l,actual,expected
             call die('independent source mismatch')
          end if
          checked=checked+1
       end if
       native_velocity(d)%source(e,id)=actual
    end do
    if(k==1.and.l==level_start) write(*,'(a,i0,a,i0)') 'Block velocity owner rank=',rank,' checked=',checked
  end subroutine
  function owner_velocity_source_key(d,id,e) result(key)
    integer,intent(in)::d,id,e
    integer::key(5)
    if(.not.velocity_owner_enabled()) call die('source address without owner execution')
    if(map(d)%source(1,e,id)/=3) call die('active gradient lacks source producer')
    key=map(d)%source(2:6,e,id)
  end function
  function owner_velocity_source_value(key,k) result(value)
    integer,intent(in)::key(5),k
    real(dp)::value
    value=program_value(plan(key(5))%code,key,k)
  end function
end module
