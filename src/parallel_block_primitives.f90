module parallel_block_primitives_mod
  use parallel_block_primitive_registry_mod
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run
  use shared_mod, only: level_start,level_end,S_MASS,S_TEMP,scalars,compressible,N_BDRY
  use domain_mod, only: Domain,grid,idx
  use domain_ops_mod, only: apply_onescale_to_patch
  use parallel_block_primitive_walk_mod, only: walk_primitive_patch
  use parallel_block_primitive_kernels_mod, only: primitive_regular,primitive_southwest
  use parallel_block_columns_mod, only: owner_tendency_without_reference
  use parallel_block_columns_mod, only: columns_active,columns_oracle,column_consumer_address,owner_column_result
  use parallel_block_inverse_mod, only: inverse_operand_generation
  use parallel_block_resident_operands_mod
  implicit none
  private
  public :: primitive_owner_value
  public :: begin_primitive_phase,end_primitive_phase,consume_regular,consume_southwest
  public :: primitives_active,primitives_oracle,require_primitive_bernoulli,restrict_primitive_bernoulli
  public :: primitive_vorticity_source,primitive_vorticity_value
  public :: primitive_thermodynamic_source,primitive_thermodynamic_value
  integer, parameter :: N_META=47,N_GEO=129,N_OUT=39
  logical :: primitives_active=.false.,primitives_oracle=.false.
  type :: Primitive_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),ids(:,:),local_domain(:),order(:),meta(:,:),destination(:)
     real(dp), allocatable :: geometry(:,:),computed(:,:),returned(:,:)
     logical, allocatable :: southwest(:),has_pressure(:)
     type(Resident_Operand_Plan) :: input
     integer :: cursor=0,last=0,layer=0
     logical :: sw_consumed=.true.
  end type
  type :: Integer_Map
     integer, allocatable :: visit(:),source(:,:)
  end type
  type(Integer_Map), allocatable :: bernoulli_ready(:)
  type(Primitive_Plan), allocatable :: plan(:)
  integer(int64) :: generation=-1_int64,executed(2)=0_int64,checked(2)=0_int64
  integer :: active_layer=0,active_level=0
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': block primitives: ',message
    call abort_run
    error stop 'block primitives abort returned'
  end subroutine
  subroutine offsets(count,disp)
    integer, intent(in) :: count(:)
    integer, intent(out) :: disp(:)
    integer :: r
    disp=0
    do r=2,size(count)
       disp(r)=disp(r-1)+count(r-1)
    end do
  end subroutine
  subroutine compile_plans
    type(Integer_Map), allocatable :: visits(:)
    integer :: l,d,j,pass,n,r,i,pos,dest,column,ierr,meta(N_META),key(5)
    integer, allocatable :: local_meta(:,:),send_meta(:,:),count(:),input_keys(:,:)
    real(dp), allocatable :: local_geo(:,:),send_geo(:,:)
    if(allocated(plan)) deallocate(plan,bernoulli_ready)
    allocate(plan(level_start:level_end),visits(size(grid)),bernoulli_ready(size(grid)))
    do d=1,size(grid)
       allocate(visits(d)%visit(grid(d)%node%length),bernoulli_ready(d)%visit(grid(d)%node%length), &
            bernoulli_ready(d)%source(4,grid(d)%node%length))
       bernoulli_ready(d)%visit=0
    end do
    do l=level_start,level_end
       do d=1,size(grid)
          visits(d)%visit=0
       end do
       do pass=1,2
          n=0
          do d=1,size(grid)
             do j=1,grid(d)%lev(l)%length
                if(pass==2) call apply_onescale_to_patch(record_pressure,grid(d),grid(d)%lev(l)%elts(j),0,0,1)
                call walk_primitive_patch(grid(d),grid(d)%lev(l)%elts(j),record_operation)
             end do
          end do
          if(pass==1) then
             allocate(local_meta(N_META,n),local_geo(N_GEO,n))
             allocate(plan(l)%ids(0:7,n),plan(l)%local_domain(n),plan(l)%order(n), &
                  plan(l)%destination(n),plan(l)%southwest(n),plan(l)%has_pressure(n))
             plan(l)%southwest=.false.
          end if
       end do
       associate(q=>plan(l))
       allocate(q%sc(n_process),q%sd(n_process),q%rc(n_process),q%rd(n_process),count(n_process))
       q%sc=0
       do i=1,n
          call column_consumer_address(q%local_domain(i),q%ids(0,i),dest,column)
          q%sc(dest+1)=q%sc(dest+1)+1
       end do
       call MPI_Alltoall(q%sc,1,MPI_INTEGER,q%rc,1,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('primitive plan counts')
       call offsets(q%sc,q%sd);call offsets(q%rc,q%rd)
       allocate(send_meta(N_META,n),send_geo(N_GEO,n),q%meta(N_META,sum(q%rc)),q%geometry(N_GEO,sum(q%rc)))
       count=0
       do i=1,n
          call column_consumer_address(q%local_domain(i),q%ids(0,i),dest,column)
          r=dest+1;pos=q%sd(r)+count(r)+1;count(r)=count(r)+1
          q%order(i)=pos;q%destination(i)=dest
          send_meta(:,pos)=local_meta(:,i);send_meta(6,pos)=merge(1,0,q%southwest(i))
          send_geo(:,pos)=local_geo(:,i)
       end do
       call MPI_Alltoallv(send_meta,N_META*q%sc,N_META*q%sd,MPI_INTEGER, &
            q%meta,N_META*q%rc,N_META*q%rd,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('primitive operation addresses')
       call MPI_Alltoallv(send_geo,N_GEO*q%sc,N_GEO*q%sd,MPI_DOUBLE_PRECISION, &
            q%geometry,N_GEO*q%rc,N_GEO*q%rd,MPI_DOUBLE_PRECISION,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('primitive geometry')
       allocate(input_keys(5,8*sum(q%rc)))
       do i=1,size(q%meta,2)
          input_keys(:,8*(i-1)+1:8*i)=reshape(q%meta(8:N_META,i),[5,8])
       end do
       call compile_resident_operands(q%input,input_keys)
       allocate(q%computed(N_OUT,sum(q%rc)),q%returned(N_OUT,n))
       deallocate(input_keys,send_meta,send_geo,local_meta,local_geo,count)
       end associate
    end do
    generation=inverse_operand_generation()
  contains
    subroutine record_pressure(dom,x,y,k,offs,dims)
      type(Domain), intent(inout) :: dom
      integer, intent(in) :: x,y,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
      integer :: id
      id=idx(x,y,offs,dims)+1
      visits(dom%id+1)%visit(id)=visits(dom%id+1)%visit(id)+1
    end subroutine
    subroutine record_operation(ids,southwest)
      integer, intent(in) :: ids(0:7)
      logical, intent(in) :: southwest
      real(dp) :: length(3,0:7),pedlen(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
      integer :: t,id
      if(southwest) then
         if(pass==2) then
            if(n==0) call die('southwest without regular producer')
            if(any(plan(l)%ids(:,n)/=ids)) call die('southwest producer identity')
            plan(l)%southwest(n)=.true.
         end if
         return
      end if
      n=n+1
      if(pass==1) return
      plan(l)%ids(:,n)=ids;plan(l)%local_domain(n)=d
      call column_consumer_address(d,ids(0),dest,column)
      if(visits(d)%visit(ids(0)+1)==0) column=0
      plan(l)%has_pressure(n)=column>0
      meta(1:7)=[rank,d,ids(0),column,visits(d)%visit(ids(0)+1),0,l]
      do t=0,7
         key=resident_operand_key(d,ids(t));meta(8+5*t:12+5*t)=key
         id=ids(t)
         length(:,t)=grid(d)%len%elts(3*id+1:3*id+3)
         pedlen(:,t)=grid(d)%pedlen%elts(3*id+1:3*id+3)
         area(:,t)=grid(d)%areas%elts(id+1)%part
         triarea(:,t)=grid(d)%triarea%elts(2*id+1:2*id+2)
         coriolis(:,t)=grid(d)%coriolis%elts(2*id+1:2*id+2)
      end do
      local_meta(:,n)=meta
      local_geo(:,n)=[reshape(length,[24]),reshape(pedlen,[24]),reshape(area,[48]), &
           reshape(triarea,[16]),reshape(coriolis,[16]),grid(d)%areas%elts(ids(0)+1)%hex_inv]
    end subroutine
  end subroutine

  subroutine initialize_registry(k)
    integer,intent(in)::k
    integer::l,i,d,e
    call reset_primitive_registry(k,generation)
    if(.not.registry_needs_declarations) return
    do l=lbound(plan,1),ubound(plan,1)
       associate(q=>plan(l))
       do i=1,size(q%order)
          d=q%local_domain(i)
          do e=1,3
             call declare_primitive_result(d,q%ids(0,i),e,.false.)
          end do
          do e=1,2
             call declare_primitive_result(d,q%ids(0,i),e,.true.)
          end do
          call declare_primitive_result(d,q%ids(4,i),1,.true.)
          call declare_primitive_result(d,q%ids(3,i),2,.true.)
          if(q%southwest(i)) then
             call declare_primitive_result(d,q%ids(4,i),1,.false.)
             call declare_primitive_result(d,q%ids(7,i),2,.false.)
             call declare_primitive_result(d,q%ids(3,i),3,.false.)
             do e=1,2
                call declare_primitive_result(d,q%ids(7,i),e,.true.)
             end do
          end if
       end do
       end associate
    end do
    call finish_primitive_declarations
  end subroutine

  subroutine begin_primitive_phase(k,l)
    integer, intent(in) :: k,l
    integer :: i,j,a,ierr,d,key(7),e
    real(dp) :: rho(2,0:7),u(3,0:7),length(3,0:7),pedlen(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
    real(dp) :: qe(3,0:7),vort(2,0:7),ke,b,context(2),dual(3),pressure(5)
    call end_primitive_phase
#ifdef BLOCK_PRIMITIVES
    primitives_active=columns_active
#endif
    if(.not.primitives_active) return
    if(.not.compressible) call die('incompressible producer unsupported')
    primitives_oracle=columns_oracle
#ifdef BLOCK_PRIMITIVES_ORACLE
    primitives_oracle=.true.
#endif
    if(generation/=inverse_operand_generation()) call compile_plans
    if(l==level_end) then
       call initialize_registry(k)
       do i=1,size(bernoulli_ready)
          bernoulli_ready(i)%visit=0
       end do
    end if
    active_layer=k;active_level=l
    associate(q=>plan(l))
    q%cursor=0;q%last=0;q%layer=k;q%sw_consumed=.true.
    call fetch_resident_operands(q%input,k)
    do i=1,size(q%meta,2)
       length=reshape(q%geometry(1:24,i),[3,8]);pedlen=reshape(q%geometry(25:48,i),[3,8])
       area=reshape(q%geometry(49:96,i),[6,8]);triarea=reshape(q%geometry(97:112,i),[2,8])
       coriolis=reshape(q%geometry(113:128,i),[2,8])
       do j=0,7
          a=q%input%address(8*(i-1)+j+1)
          rho(:,j)=q%input%recv(1:2,a);u(:,j)=q%input%recv(3:5,a)
       end do
       pressure=0.0_dp
       if(q%meta(4,i)>0) pressure=owner_column_result(q%meta(1,i),q%meta(4,i),q%meta(5,i),k,l)
       qe=0.0_dp;vort=0.0_dp
       call primitive_regular(rho(1,:),u,length,pedlen,area,q%geometry(129,i),triarea,coriolis, &
            pressure(2),pressure(5),qe,vort,ke,b,context,dual)
       q%computed(1:17,i)=[qe(:,0),vort(:,0),vort(1,4),vort(2,3),ke,b,context,dual, &
            u(1,4)*pedlen(1,4),u(2,7)*pedlen(2,7),u(3,3)*pedlen(3,3)]
       q%computed(18:22,i)=0.0_dp
       if(q%meta(6,i)==1) then
          call primitive_southwest(rho(1,:),u,length,area,triarea,coriolis,context,qe,vort)
          q%computed(18:22,i)=[qe(1,4),qe(2,7),qe(3,3),vort(:,7)]
          executed(2)=executed(2)+1_int64
       end if
       q%computed(23:38,i)=reshape(rho,[16])
       q%computed(39,i)=pressure(4)
       executed(1)=executed(1)+1_int64
    end do
    call MPI_Alltoallv(q%computed,N_OUT*q%rc,N_OUT*q%rd,MPI_DOUBLE_PRECISION, &
         q%returned,N_OUT*q%sc,N_OUT*q%sd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('primitive result delivery')
    do i=1,size(q%order)
       d=q%local_domain(i)
       key=[1,q%destination(i),rank,q%order(i)-q%sd(q%destination(i)+1),0,l,0]
       do e=1,3
          key(5)=e;call register_primitive_result(d,q%ids(0,i),e,.false.,key,k)
       end do
       do e=1,2
          key(5)=3+e;call register_primitive_result(d,q%ids(0,i),e,.true.,key,k)
       end do
       key(5)=6;call register_primitive_result(d,q%ids(4,i),1,.true.,key,k)
       key(5)=7;call register_primitive_result(d,q%ids(3,i),2,.true.,key,k)
       if(q%southwest(i)) then
          key(5)=18;call register_primitive_result(d,q%ids(4,i),1,.false.,key,k)
          key(5)=19;call register_primitive_result(d,q%ids(7,i),2,.false.,key,k)
          key(5)=20;call register_primitive_result(d,q%ids(3,i),3,.false.,key,k)
          do e=1,2
             key(5)=20+e;call register_primitive_result(d,q%ids(7,i),e,.true.,key,k)
          end do
       end if
    end do
    end associate
  end subroutine

  function primitive_owner_value(key,k) result(value)
    integer,intent(in)::key(5),k
    real(dp)::value
    integer::l,s
    ! owner, origin, relative packet, output row, producing level
    l=key(5)
    if(key(1)/=rank.or.l<lbound(plan,1).or.l>ubound(plan,1)) call die('owner result address')
    if(plan(l)%layer/=k) call die('owner result layer')
    if(key(2)<0.or.key(2)>=n_process) call die('owner result origin')
    if(key(3)<1.or.key(3)>plan(l)%rc(key(2)+1)) call die('owner result packet')
    if(key(4)<1.or.key(4)>N_OUT) call die('owner result row')
    s=plan(l)%rd(key(2)+1)+key(3);value=plan(l)%computed(key(4),s)
  end function

  function primitive_vorticity_source(d,id,component,l) result(source)
    integer, intent(in) :: d,id,component,l
    integer :: source(2),i
    source=0
    if(component<1.or.component>2) call die('vorticity component')
    associate(q=>plan(l))
    do i=size(q%order),1,-1
       if(q%local_domain(i)/=d) cycle
       if(q%southwest(i).and.q%ids(7,i)==id) then
          source=[i,20+component];return
       end if
       if(component==1.and.q%ids(4,i)==id) then
          source=[i,6];return
       end if
       if(component==2.and.q%ids(3,i)==id) then
          source=[i,7];return
       end if
       if(q%ids(0,i)==id) then
          source=[i,3+component];return
       end if
    end do
    end associate
    call die('pentagon vorticity has no current primitive producer')
  end function

  real(dp) function primitive_vorticity_value(source,l) result(value)
    integer, intent(in) :: source(2),l
    if(.not.primitives_active.or.l/=active_level) call die('vorticity producer epoch')
    if(source(1)<1.or.source(1)>size(plan(l)%order)) call die('vorticity producer address')
    if(source(2)<4.or.source(2)>22) call die('vorticity result component')
    value=plan(l)%returned(source(2),plan(l)%order(source(1)))
  end function

  subroutine exact(actual,expected,label)
    real(dp), intent(in) :: actual(:),expected(:)
    character(*), intent(in) :: label
    if(size(actual)/=size(expected)) call die('oracle extent')
    if(any(transfer(actual,[0_int64],size(actual))/=transfer(expected,[0_int64],size(expected)))) then
       write(*,*) 'actual:',actual,' expected:',expected
       call die('independent mismatch: '//label)
    end if
  end subroutine
  subroutine consume_regular(dom,ids,k,data)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: ids(0:7),k
    real(dp), intent(out) :: data(38)
    integer :: n,pos,id
    if(.not.primitives_active.or.k/=active_layer) call die('regular consumer epoch')
    associate(q=>plan(active_level))
    if(.not.q%sw_consumed) call die('southwest consumer was skipped')
    q%cursor=q%cursor+1;n=q%cursor
    if(n>size(q%order)) call die('extra regular consumer')
    if(q%local_domain(n)/=dom%id+1.or.any(q%ids(:,n)/=ids)) call die('regular consumer order')
    q%last=n;q%sw_consumed=.not.q%southwest(n)
    pos=q%order(n);data=q%returned(1:38,pos);id=ids(0)
    if(primitives_oracle) then
       call exact([dom%qe%elts(3*id+1:3*id+3),dom%vort%elts(2*id+1:2*id+2), &
            dom%vort%elts(2*ids(4)+1),dom%vort%elts(2*ids(3)+2),dom%ke%elts(id+1)],data(1:8),'regular primitive')
       if(q%has_pressure(n)) call exact([dom%bernoulli%elts(id+1)],data(9:9),'Bernoulli')
       checked(1)=checked(1)+1_int64
    end if
    if(.not.owner_tendency_without_reference()) then
    dom%qe%elts(3*id+1:3*id+3)=data(1:3)
    dom%vort%elts(2*id+1:2*id+2)=data(4:5)
    dom%vort%elts(2*ids(4)+1)=data(6);dom%vort%elts(2*ids(3)+2)=data(7)
    dom%ke%elts(id+1)=data(8)
    end if
    bernoulli_ready(dom%id+1)%visit(id+1)=merge(1,0,q%has_pressure(n))
    if(q%has_pressure(n)) then
       if(.not.owner_tendency_without_reference()) dom%bernoulli%elts(id+1)=data(9)
       bernoulli_ready(dom%id+1)%source(:,id+1)= &
            [q%destination(n),rank,pos-q%sd(q%destination(n)+1),active_level]
    end if
    end associate
  end subroutine
  subroutine require_primitive_bernoulli(d,ids)
    integer, intent(in) :: d,ids(:)
    if(.not.primitives_active) return
    if(any(bernoulli_ready(d)%visit(ids)==0)) call die('live Bernoulli consumer has no current producer')
  end subroutine
  subroutine restrict_primitive_bernoulli(d,parent,child)
    integer, intent(in) :: d,parent,child
    if(.not.primitives_active) return
    call require_primitive_bernoulli(d,[child])
    bernoulli_ready(d)%visit(parent)=bernoulli_ready(d)%visit(child)
    bernoulli_ready(d)%source(:,parent)=bernoulli_ready(d)%source(:,child)
  end subroutine
  function primitive_thermodynamic_source(d,id) result(source)
    integer, intent(in) :: d,id ! one based node
    integer :: source(4)
    call require_primitive_bernoulli(d,[id])
    if(.not.primitives_active) call die('thermodynamic source outside transaction')
    source=bernoulli_ready(d)%source(:,id)
  end function
  function primitive_thermodynamic_value(source,k) result(value)
    integer, intent(in) :: source(3),k ! origin rank, relative operation, level
    real(dp) :: value(2)
    integer :: r,p,l
    r=source(1)+1;l=source(3)
    if(.not.primitives_active.or.k/=active_layer) call die('thermodynamic operand epoch')
    if(l<level_start.or.l>level_end.or.r<1.or.r>n_process) call die('thermodynamic source address')
    associate(q=>plan(l))
    if(q%layer/=k) call die('stale thermodynamic layer')
    if(source(2)<1.or.source(2)>q%rc(r)) call die('thermodynamic operation address')
    p=q%rd(r)+source(2)
    if(q%meta(4,p)==0) call die('thermodynamic operand has no pressure producer')
    value=[q%computed(9,p),q%computed(39,p)]
    end associate
  end function
  subroutine consume_southwest(dom,ids,k,data)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: ids(0:7),k
    real(dp), intent(out) :: data(38)
    integer :: n,pos
    if(.not.primitives_active.or.k/=active_layer) call die('southwest consumer epoch')
    associate(q=>plan(active_level))
    n=q%last
    if(n==0.or.q%sw_consumed) call die('southwest consumer has no regular producer')
    if(q%local_domain(n)/=dom%id+1.or.any(q%ids([0,3,4,7],n)/=ids([0,3,4,7]))) &
         call die('southwest consumer identity')
    q%sw_consumed=.true.;pos=q%order(n);data=q%returned(1:38,pos)
    if(primitives_oracle) then
       call exact([dom%qe%elts(3*ids(4)+1),dom%qe%elts(3*ids(7)+2),dom%qe%elts(3*ids(3)+3), &
            dom%vort%elts(2*ids(7)+1:2*ids(7)+2)],data(18:22),'southwest primitive')
       checked(2)=checked(2)+1_int64
    end if
    if(.not.owner_tendency_without_reference()) then
    dom%qe%elts(3*ids(4)+1)=data(18);dom%qe%elts(3*ids(7)+2)=data(19);dom%qe%elts(3*ids(3)+3)=data(20)
    dom%vort%elts(2*ids(7)+1:2*ids(7)+2)=data(21:22)
    end if
    end associate
  end subroutine
  subroutine end_primitive_phase
    if(.not.primitives_active) return
    associate(q=>plan(active_level))
    if(q%cursor/=size(q%order).or..not.q%sw_consumed) call die('incomplete primitive consumers')
    end associate
    if(active_level==level_start) write(*,'(a,i0,a,4(1x,i0))') &
         'Block primitives rank=',rank,' regular/SW executed/checked=',executed,checked
    primitives_active=.false.
  end subroutine
end module
