module parallel_block_pentagons_mod
  use parallel_block_primitive_registry_mod
  ! Pentagon dependencies refer to primitive result packets, never Domain
  ! vorticity inputs. The independent reference remains in ops_mod.
  use parallel_block_columns_mod, only: owner_tendency_without_reference
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run
  use shared_mod, only: N_BDRY,level_start,level_end,IJMINUS,IJPLUS,IPLUSJMINUS,IMINUSJPLUS
  use domain_mod, only: Domain,grid,idx
  use domain_ops_mod, only: apply_to_penta_d
  use patch_mod, only: PATCH_SIZE
  use parallel_block_primitives_mod, only: primitives_active,primitives_oracle, &
       primitive_vorticity_source,primitive_vorticity_value
  use parallel_block_columns_mod, only: column_consumer_address
  use parallel_block_primitive_kernels_mod, only: primitive_pentagon
  use parallel_block_resident_operands_mod
  use parallel_block_inverse_mod, only: inverse_operand_generation
  implicit none
  private
  public :: pentagon_owner_value
  public :: begin_pentagon_phase,end_pentagon_phase,consume_pentagon
  integer, parameter :: NG=105
  type :: Pentagon_Plan
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),order(:),identity(:,:),ids(:,:),corner(:),source(:,:,:),destination(:)
     real(dp), allocatable :: geometry(:,:),incoming(:,:),outgoing(:,:),result(:,:)
     type(Resident_Operand_Plan) :: input
     integer :: cursor=0,layer=0
  end type
  type(Pentagon_Plan), allocatable :: plan(:)
  logical :: active=.false.
  integer :: active_level=0
  integer(int64) :: generation=-1_int64,executed=0_int64,checked=0_int64
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': block pentagons: ',message
    call abort_run
    error stop 'block pentagons abort returned'
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
  subroutine slots(c,qnode,qedge,vnode,vedge,prior_node,prior_edge)
    integer, intent(in) :: c
    integer, intent(out) :: qnode(2),qedge(2),vnode(2),vedge(2),prior_node(2),prior_edge(2)
    select case(c)
    case(IJMINUS)
       qnode=[4,3];qedge=[1,3];vnode=[7,7];vedge=[1,2];prior_node=[4,3];prior_edge=[1,2]
    case(IPLUSJMINUS)
       qnode=[0,7];qedge=[1,2];vnode=[7,3];vedge=[1,2];prior_node=[0,7];prior_edge=[1,2]
    case(IMINUSJPLUS)
       qnode=[0,7];qedge=[3,2];vnode=[7,4];vedge=[2,1];prior_node=[0,7];prior_edge=[2,1]
    case(IJPLUS)
       qnode=[0,0];qedge=[1,3];vnode=[0,0];vedge=[1,2];prior_node=[4,3];prior_edge=[1,2]
    case default
       call die('unknown corner')
    end select
  end subroutine
  subroutine compile_plans
    integer :: l,d,pass,n,dest,column,r,i,pos,ierr,j,qnode(2),qedge(2),vnode(2),vedge(2),pn(2),pe(2)
    integer, allocatable :: local_keys(:,:),send_keys(:,:),recv_keys(:,:),count(:),send_corner(:),destinations(:)
    real(dp), allocatable :: local_geo(:,:),send_geo(:,:)
    if(allocated(plan)) deallocate(plan)
    allocate(plan(level_start:level_end))
    do l=level_start,level_end
       do pass=1,2
          n=0
          do d=1,size(grid)
             call apply_to_penta_d(record,grid(d),l,0)
          end do
          if(pass==1) then
             allocate(plan(l)%identity(3,n),plan(l)%ids(0:7,n),plan(l)%source(2,2,n), &
                  local_keys(40,n),local_geo(NG,n),destinations(n))
          end if
       end do
       ! A later correction must not silently consume an earlier correction
       ! as if it were still a regular result. Reject such a dependency until
       ! represented as an explicit additional phase.
       do i=1,n
          call slots(plan(l)%identity(3,i),qnode,qedge,vnode,vedge,pn,pe)
          do j=1,i-1
             if(plan(l)%identity(1,j)/=plan(l)%identity(1,i)) cycle
             call reject_prior_write(i,j,pn,pe)
          end do
       end do
       associate(q=>plan(l))
       allocate(q%sc(n_process),q%sd(n_process),q%rc(n_process),q%rd(n_process),count(n_process),q%order(n))
       q%destination=destinations
       q%sc=0
       do i=1,n
          q%sc(destinations(i)+1)=q%sc(destinations(i)+1)+1
       end do
       call MPI_Alltoall(q%sc,1,MPI_INTEGER,q%rc,1,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('counts')
       call offsets(q%sc,q%sd);call offsets(q%rc,q%rd)
       allocate(send_keys(40,n),recv_keys(40,sum(q%rc)),send_corner(n),q%corner(sum(q%rc)), &
            send_geo(NG,n),q%geometry(NG,sum(q%rc)))
       count=0
       do i=1,n
          r=destinations(i)+1;pos=q%sd(r)+count(r)+1;count(r)=count(r)+1;q%order(i)=pos
          send_keys(:,pos)=local_keys(:,i);send_geo(:,pos)=local_geo(:,i);send_corner(pos)=q%identity(3,i)
       end do
       call MPI_Alltoallv(send_keys,40*q%sc,40*q%sd,MPI_INTEGER,recv_keys,40*q%rc,40*q%rd,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('operand keys')
       call MPI_Alltoallv(send_corner,q%sc,q%sd,MPI_INTEGER,q%corner,q%rc,q%rd,MPI_INTEGER,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('corners')
       call MPI_Alltoallv(send_geo,NG*q%sc,NG*q%sd,MPI_DOUBLE_PRECISION, &
            q%geometry,NG*q%rc,NG*q%rd,MPI_DOUBLE_PRECISION,comm,ierr)
       if(ierr/=MPI_SUCCESS) call die('geometry')
       call compile_resident_operands(q%input,reshape(recv_keys,[5,8*sum(q%rc)]))
       allocate(q%incoming(2,sum(q%rc)),q%outgoing(4,sum(q%rc)),q%result(4,n))
       deallocate(local_keys,local_geo,destinations,send_keys,recv_keys,send_corner,send_geo,count)
       end associate
    end do
    generation=inverse_operand_generation()
  contains
    subroutine record(dom,p,c,offs,dims,zlev)
      type(Domain), intent(inout) :: dom
      integer, intent(in) :: p,c,zlev,offs(N_BDRY+1),dims(2,N_BDRY+1)
      integer :: ids(0:7),t,id,key(5),x,y
      real(dp) :: length(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
      n=n+1
      if(pass==1) return
      x=0;y=0
      if(c==IPLUSJMINUS.or.c==IJPLUS) x=PATCH_SIZE
      if(c==IMINUSJPLUS.or.c==IJPLUS) y=PATCH_SIZE
      ids=[idx(x,y,offs,dims),idx(x,y+1,offs,dims),idx(x+1,y,offs,dims), &
           idx(x,y-1,offs,dims),idx(x-1,y,offs,dims),idx(x+1,y+1,offs,dims),idx(x,y,offs,dims), &
           idx(x-1,y-1,offs,dims)]
      ! Match the exact operand slots of each original pentagon expression.
      select case(c)
      case(IJMINUS)
         ids(5)=ids(0)
      case(IPLUSJMINUS)
         ids(1)=ids(0)
      case(IMINUSJPLUS)
         ids(2)=ids(0)
      case(IJPLUS)
         ids(5:7)=ids(0)
      end select
      plan(l)%ids(:,n)=ids;plan(l)%identity(:,n)=[d,p,c]
      call column_consumer_address(d,ids(0),dest,column);destinations(n)=dest
      call slots(c,qnode,qedge,vnode,vedge,pn,pe)
      do t=1,2
         plan(l)%source(:,t,n)=primitive_vorticity_source(d,ids(pn(t)),pe(t),l)
      end do
      do t=0,7
         key=resident_operand_key(d,ids(t));local_keys(5*t+1:5*t+5,n)=key
         id=ids(t)
         length(:,t)=dom%len%elts(3*id+1:3*id+3);area(:,t)=dom%areas%elts(id+1)%part
         triarea(:,t)=dom%triarea%elts(2*id+1:2*id+2);coriolis(:,t)=dom%coriolis%elts(2*id+1:2*id+2)
      end do
      local_geo(:,n)=[reshape(length,[24]),reshape(area,[48]),reshape(triarea,[16]),reshape(coriolis,[16]),0.0_dp]
    end subroutine
    subroutine reject_prior_write(current,previous,pn,pe)
      integer, intent(in) :: current,previous,pn(2),pe(2)
      integer :: qn(2),qe(2),vn(2),ve(2),dummy1(2),dummy2(2),a,b
      call slots(plan(l)%identity(3,previous),qn,qe,vn,ve,dummy1,dummy2)
      do a=1,2
         do b=1,2
            if(plan(l)%ids(pn(a),current)==plan(l)%ids(vn(b),previous).and.pe(a)==ve(b)) &
                 call die('ordered pentagon dependency requires another producer phase')
         end do
      end do
    end subroutine
  end subroutine
  subroutine begin_pentagon_phase(k,l)
    integer, intent(in) :: k,l
    integer :: i,j,a,pos,ierr,qn(2),qeidx(2),vn(2),ve(2),pn(2),pe(2),key(7),d
    real(dp), allocatable :: prior(:,:)
    real(dp) :: rho(0:7),u(3,0:7),length(3,0:7),area(6,0:7),triarea(2,0:7),coriolis(2,0:7)
    real(dp) :: qe(3,0:7),vort(2,0:7)
    call end_pentagon_phase
    if(.not.primitives_active) return
    if(generation/=inverse_operand_generation()) call compile_plans
    active=.true.;active_level=l
    associate(q=>plan(l))
    q%cursor=0;q%layer=k
    allocate(prior(2,size(q%order)))
    do i=1,size(q%order)
       pos=q%order(i)
       do j=1,2
          prior(j,pos)=primitive_vorticity_value(q%source(:,j,i),l)
       end do
    end do
    call MPI_Alltoallv(prior,2*q%sc,2*q%sd,MPI_DOUBLE_PRECISION,q%incoming,2*q%rc,2*q%rd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('prior vorticity delivery')
    call fetch_resident_operands(q%input,k)
    do i=1,size(q%corner)
       do j=0,7
          a=q%input%address(8*(i-1)+j+1);rho(j)=q%input%recv(1,a);u(:,j)=q%input%recv(3:5,a)
       end do
       length=reshape(q%geometry(1:24,i),[3,8]);area=reshape(q%geometry(25:72,i),[6,8])
       triarea=reshape(q%geometry(73:88,i),[2,8]);coriolis=reshape(q%geometry(89:104,i),[2,8])
       call slots(q%corner(i),qn,qeidx,vn,ve,pn,pe)
       qe=0.0_dp;vort=0.0_dp
       do j=1,2
          vort(pe(j),pn(j))=q%incoming(j,i)
       end do
       call primitive_pentagon(q%corner(i),rho,u,length,area,triarea,coriolis,qe,vort)
       q%outgoing(:,i)=[qe(qeidx(1),qn(1)),qe(qeidx(2),qn(2)),vort(ve(1),vn(1)),vort(ve(2),vn(2))]
       executed=executed+1_int64
    end do
    call MPI_Alltoallv(q%outgoing,4*q%rc,4*q%rd,MPI_DOUBLE_PRECISION, &
         q%result,4*q%sc,4*q%sd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('result delivery')
    do i=1,size(q%order)
       d=q%identity(1,i)
       key=[2,q%destination(i),rank,q%order(i)-q%sd(q%destination(i)+1),0,l,0]
       call slots(q%identity(3,i),qn,qeidx,vn,ve,pn,pe)
       do j=1,2
          key(5)=j;call register_primitive_result(d,q%ids(qn(j),i),qeidx(j),.false.,key,k)
          key(5)=j+2;call register_primitive_result(d,q%ids(vn(j),i),ve(j),.true.,key,k)
       end do
    end do
    end associate
  end subroutine
  function pentagon_owner_value(key,k) result(value)
    integer,intent(in)::key(5),k
    real(dp)::value
    integer::l,s
    l=key(5)
    if(key(1)/=rank.or.l<lbound(plan,1).or.l>ubound(plan,1)) call die('owner result address')
    if(plan(l)%layer/=k) call die('owner result layer')
    if(key(2)<0.or.key(2)>=n_process) call die('owner result origin')
    if(key(3)<1.or.key(3)>plan(l)%rc(key(2)+1)) call die('owner result packet')
    if(key(4)<1.or.key(4)>4) call die('owner result row')
    s=plan(l)%rd(key(2)+1)+key(3);value=plan(l)%outgoing(key(4),s)
  end function
  subroutine consume_pentagon(dom,p,c)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: p,c
    integer :: n,pos,j,qn(2),qe(2),vn(2),ve(2),pn(2),pe(2),ids(0:7)
    real(dp) :: actual(4),expected(4)
    if(.not.active) call die('consumer epoch')
    associate(q=>plan(active_level))
    q%cursor=q%cursor+1;n=q%cursor
    if(n>size(q%order)) call die('extra consumer')
    if(any(q%identity(:,n)/=[dom%id+1,p,c])) call die('consumer order')
    ids=q%ids(:,n);pos=q%order(n);expected=q%result(:,pos)
    call slots(c,qn,qe,vn,ve,pn,pe)
    if(primitives_oracle) then
       actual=[dom%qe%elts(3*ids(qn(1))+qe(1)),dom%qe%elts(3*ids(qn(2))+qe(2)), &
            dom%vort%elts(2*ids(vn(1))+ve(1)),dom%vort%elts(2*ids(vn(2))+ve(2))]
       if(any(transfer(actual,[0_int64],4)/=transfer(expected,[0_int64],4))) then
          write(*,*) 'actual:',actual,' expected:',expected
          call die('independent pentagon mismatch')
       end if
       checked=checked+1_int64
    end if
    if(.not.owner_tendency_without_reference()) then
    do j=1,2
       dom%qe%elts(3*ids(qn(j))+qe(j))=expected(j)
       dom%vort%elts(2*ids(vn(j))+ve(j))=expected(j+2)
    end do
    end if
    end associate
  end subroutine
  subroutine end_pentagon_phase
    if(.not.active) return
    if(plan(active_level)%cursor/=size(plan(active_level)%order)) call die('unconsumed result')
    if(active_level==level_start) write(*,'(a,i0,a,2(1x,i0))') &
         'Block pentagons rank=',rank,' executed/checked=',executed,checked
    active=.false.
  end subroutine
end module
