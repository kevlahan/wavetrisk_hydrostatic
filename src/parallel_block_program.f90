module parallel_block_program_mod
  ! Ordered numerical dependencies. Geometry compiles addresses; final block
  ! owners execute and retain values. A reference names one immutable write,
  ! never a mutable Domain field. Component sign represents oriented aliases.
  use kind_mod, only: dp
  use iso_fortran_env, only: int64
  use mpi_f08
  use arch_mod, only: rank,n_process,comm,abort_run
  use ieee_arithmetic, only: ieee_is_finite
  use parallel_block_resident_operands_mod
  use parallel_block_velocity_mod, only: velocity_qperp,velocity_source
  use parallel_block_primitive_kernels_mod, only: velocity_diffusion_source
  use parallel_block_mass_mod, only: Mass_Restriction,mass_restriction_value
  implicit none
  private
  public :: program_value
  abstract interface
     function External_Operand(code,key,k) result(value)
       import dp
       integer,intent(in)::code,key(5),k
       real(dp)::value
     end function
  end interface
  public :: Owner_Program,Owner_Read,program_add,program_compile,program_execute,program_read_compile,program_read
  public :: OP_RESIDENT,OP_GRADIENT,OP_GRADIENT_SW,OP_VELOCITY_FLUX,OP_DIVERGENCE,OP_RESTRICTION
  public :: OP_ADD,OP_CURL,OP_CIRCULATION,OP_VELOCITY_PHYSICS,OP_VELOCITY_SOURCE
  integer,parameter :: OP_ADD=10,OP_CURL=11,OP_CIRCULATION=12,OP_VELOCITY_PHYSICS=13,OP_VELOCITY_SOURCE=14
  public :: OP_MASS_DIVERGENCE
  integer, parameter :: OP_MASS_DIVERGENCE=9
  public :: OP_ZERO
  integer, parameter :: OP_ZERO=8
  public :: OP_SCALAR_FLUX
  integer, parameter :: OP_SCALAR_FLUX=7
  integer, parameter :: OP_RESIDENT=1,OP_GRADIENT=2,OP_GRADIENT_SW=3,OP_VELOCITY_FLUX=4, &
       OP_DIVERGENCE=5,OP_RESTRICTION=6
  type :: Operation
     integer :: owner=0,origin=0,id=0,code=0,depth=0,nout=0,key(5)=0,input_start=0,source=0
     integer, allocatable :: input(:,:)
     real(dp), allocatable :: geometry(:)
     real(dp) :: value(7)=0.0_dp
     logical :: ready=.false.
  end type
  type :: Owner_Read
     integer, allocatable :: sc(:),sd(:),rc(:),rd(:),source(:,:),address(:),orientation(:)
     real(dp), allocatable :: send(:),received(:),value(:)
  end type
  type :: Owner_Program
     type(Operation), allocatable :: local(:),owned(:)
     integer :: n=0,maxdepth=0,layer=0
     integer, allocatable :: table(:)
     type(Owner_Read), allocatable :: read(:)
     type(Resident_Operand_Plan) :: resident
     logical :: compiled=.false.
  end type
contains
  subroutine die(message)
    character(*), intent(in) :: message
    write(*,'(a,i0,2a)') 'Rank ',rank,': owner program: ',message
    call abort_run
    error stop 'owner program abort returned'
  end subroutine
  subroutine offsets(c,d)
    integer,intent(in)::c(:)
    integer,intent(out)::d(:)
    integer::i
    d=0
    do i=2,size(c)
       d(i)=d(i-1)+c(i-1)
    end do
  end subroutine
  function program_add(p,owner,code,input,geometry,nout,key) result(ref)
    type(Owner_Program),intent(inout)::p
    integer,intent(in)::owner,code,input(:,:),nout
    real(dp),intent(in)::geometry(:)
    integer,intent(in),optional::key(5)
    integer::ref(5,nout),depth,i
    type(Operation),allocatable::tmp(:)
    if(p%compiled) call die('append to compiled program')
    if(owner<0.or.owner>=n_process.or.nout<1.or.nout>7.or.size(input,1)/=5) call die('operation shape')
    ! Preserve undefined geometry-only writes as undefined definitions. A live
    ! consumer is rejected by program_read_compile; no numerical value is made.
    if(any(input(3,:)<=0).or.any(input(4,:)==0)) then
       ref=0;return
    end if
    if(.not.allocated(p%local)) allocate(p%local(1024))
    if(p%n==size(p%local)) then
       allocate(tmp(2*size(p%local)))
       tmp(1:p%n)=p%local
       call move_alloc(tmp,p%local)
    end if
    depth=0
    if(size(input,2)>0) depth=maxval(input(5,:))+1
    p%n=p%n+1
    associate(o=>p%local(p%n))
    o%owner=owner;o%origin=rank;o%id=p%n;o%code=code;o%depth=depth;o%nout=nout
    o%input=input;o%geometry=geometry
    if(present(key)) o%key=key
    do i=1,nout
       ref(:,i)=[owner,rank,p%n,i,depth]
    end do
    end associate
  end function
  integer function hash(origin,id,n) result(h)
    integer,intent(in)::origin,id,n
    integer(int64)::bits
    ! Producer IDs form long contiguous ranges on each geometry rank. Merely
    ! translating those ranges by origin creates long linear-probing clusters
    ! after redistribution. Mix the bits before reducing to the table size.
    ! Every multiply stays within signed int64, including maximum default IDs.
    bits=modulo(2654435761_int64*int(origin,int64)+int(id,int64),4294967296_int64)
    bits=modulo(ieor(bits,ishft(bits,-16))*73244475_int64,4294967296_int64)
    bits=modulo(ieor(bits,ishft(bits,-16))*73244475_int64,4294967296_int64)
    bits=ieor(bits,ishft(bits,-16))
    h=int(modulo(bits,int(n,int64)))+1
  end function
  integer function locate(p,origin,id) result(s)
    type(Owner_Program),intent(in)::p
    integer,intent(in)::origin,id
    integer::h,t
    h=hash(origin,id,size(p%table))
    do t=1,size(p%table)
       s=p%table(h)
       if(s==0) call die('unknown producer identity')
       if(p%owned(s)%origin==origin.and.p%owned(s)%id==id) return
       h=modulo(h,size(p%table))+1
    end do
    call die('producer lookup overflow')
  end function
  subroutine program_compile(p)
    type(Owner_Program),intent(inout)::p
    integer,allocatable::sc(:),sd(:),rc(:),rd(:),gc(:),gd(:),hc(:),hd(:),si(:),ri(:),ic(:),jc(:),keys(:,:),refs(:,:)
    real(dp),allocatable::sg(:),rg(:)
    integer::i,r,a,b,n,h,t,depth,ierr,ni,ng,ns,cursor
    if(p%compiled) call die('duplicate compile')
    allocate(sc(n_process),sd(n_process),rc(n_process),rd(n_process),gc(n_process),gd(n_process), &
         hc(n_process),hd(n_process),ic(n_process),jc(n_process))
    sc=0;gc=0
    do i=1,p%n
       r=p%local(i)%owner+1
       sc(r)=sc(r)+12+size(p%local(i)%input)
       gc(r)=gc(r)+size(p%local(i)%geometry)
    end do
    call MPI_Alltoall(sc,1,MPI_INTEGER,rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('operation counts')
    call MPI_Alltoall(gc,1,MPI_INTEGER,hc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('geometry counts')
    call offsets(sc,sd);call offsets(rc,rd);call offsets(gc,gd);call offsets(hc,hd)
    allocate(si(sum(sc)),ri(sum(rc)),sg(sum(gc)),rg(sum(hc)))
    ic=sd;jc=gd
    do i=1,p%n
       associate(o=>p%local(i))
       r=o%owner+1;a=ic(r);b=jc(r);ni=size(o%input,2);ng=size(o%geometry)
       si(a+1:a+12)=[o%origin,o%id,o%code,o%depth,o%nout,ni,ng,o%key]
       si(a+13:a+12+5*ni)=reshape(o%input,[5*ni]);sg(b+1:b+ng)=o%geometry
       ic(r)=a+12+5*ni;jc(r)=b+ng
       end associate
    end do
    call MPI_Alltoallv(si,sc,sd,MPI_INTEGER,ri,rc,rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('operation delivery')
    call MPI_Alltoallv(sg,gc,gd,MPI_DOUBLE_PRECISION,rg,hc,hd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('geometry delivery')
    a=0;n=0
    do while(a<size(ri))
       n=n+1;a=a+12+5*ri(a+6)
    end do
    if(a/=size(ri)) call die('operation packet extent')
    allocate(p%owned(n),p%table(max(1,2*n+1)))
    p%table=0;a=0;b=0;depth=0;ns=0
    do i=1,n
       associate(o=>p%owned(i))
       o%owner=rank;o%origin=ri(a+1);o%id=ri(a+2);o%code=ri(a+3);o%depth=ri(a+4)
       o%nout=ri(a+5);ni=ri(a+6);ng=ri(a+7);o%key=ri(a+8:a+12)
       o%input=reshape(ri(a+13:a+12+5*ni),[5,ni]);o%geometry=rg(b+1:b+ng)
       a=a+12+5*ni;b=b+ng;depth=max(depth,o%depth)
       if(o%code==OP_RESIDENT) then
          ns=ns+1;o%source=ns
       end if
       h=hash(o%origin,o%id,size(p%table))
       do while(p%table(h)/=0)
          t=p%table(h)
          if(p%owned(t)%origin==o%origin.and.p%owned(t)%id==o%id) call die('duplicate producer')
          h=modulo(h,size(p%table))+1
       end do
       p%table(h)=i
       end associate
    end do
    if(allocated(p%local)) deallocate(p%local)
    allocate(keys(5,ns));ns=0
    do i=1,n
       if(p%owned(i)%code/=OP_RESIDENT) cycle
       ns=ns+1;keys(:,ns)=p%owned(i)%key
    end do
    call compile_resident_operands(p%resident,keys)
    call MPI_Allreduce(depth,p%maxdepth,1,MPI_INTEGER,MPI_MAX,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('program depth')
    allocate(p%read(0:p%maxdepth))
    do depth=0,p%maxdepth
       ni=0
       do i=1,n
          if(p%owned(i)%depth==depth) ni=ni+size(p%owned(i)%input,2)
       end do
       allocate(refs(5,ni));cursor=0
       do i=1,n
          associate(o=>p%owned(i))
          if(o%depth/=depth) cycle
          ni=size(o%input,2);o%input_start=cursor
          refs(:,cursor+1:cursor+ni)=o%input;cursor=cursor+ni
          end associate
       end do
       call program_read_compile(p,p%read(depth),refs,depth)
       deallocate(refs)
    end do
    p%compiled=.true.
    write(*,'(a,i0,a,3(i0,1x))') 'Block owner program rank=',rank,' local/owned/depth=',p%n,n,p%maxdepth
  end subroutine
  subroutine program_read_compile(p,q,refs,before)
    type(Owner_Program),intent(in)::p
    type(Owner_Read),intent(out)::q
    integer,intent(in)::refs(:,:)
    integer,intent(in),optional::before
    integer,allocatable::packed(:,:),incoming(:,:),count(:),table(:),unique(:,:),mapping(:),order(:)
    integer::i,j,r,a,h,s,n,u,ierr,identity(5)
    integer(int64)::hv
    n=size(refs,2)
    allocate(q%sc(n_process),q%sd(n_process),q%rc(n_process),q%rd(n_process),count(n_process), &
         q%address(n),q%orientation(n),table(max(1,2*n+1)),unique(5,n),mapping(n),order(n))
    table=0;u=0
    do i=1,n
       identity=refs(:,i);q%orientation(i)=sign(1,identity(4));identity(4)=abs(identity(4))
       if(identity(3)<=0.or.identity(4)<=0) call die('undefined consumer input')
       if(present(before)) then
          if(identity(5)>=before) call die('noncausal consumer')
       end if
       hv=0
       do j=1,5
          hv=modulo(65599_int64*hv+identity(j),2147483647_int64)
       end do
       h=int(modulo(hv,int(size(table),int64)))+1
       do
          s=table(h)
          if(s==0) then
             u=u+1;unique(:,u)=identity;table(h)=u;s=u;exit
          end if
          if(all(unique(:,s)==identity)) exit
          h=modulo(h,size(table))+1
       end do
       mapping(i)=s
    end do
    q%sc=0
    do i=1,u
       r=unique(1,i)+1
       if(r<1.or.r>n_process) call die('consumer owner')
       q%sc(r)=q%sc(r)+1
    end do
    call MPI_Alltoall(q%sc,1,MPI_INTEGER,q%rc,1,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('read counts')
    call offsets(q%sc,q%sd);call offsets(q%rc,q%rd)
    allocate(packed(4,u),incoming(4,sum(q%rc)),q%source(2,sum(q%rc)))
    count=0
    do i=1,u
       r=unique(1,i)+1;a=q%sd(r)+count(r)+1;count(r)=count(r)+1
       packed(:,a)=unique(2:5,i);order(i)=a
    end do
    do i=1,n
       q%address(i)=order(mapping(i))
    end do
    call MPI_Alltoallv(packed,4*q%sc,4*q%sd,MPI_INTEGER,incoming,4*q%rc,4*q%rd,MPI_INTEGER,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('read identities')
    do i=1,size(incoming,2)
       s=locate(p,incoming(1,i),incoming(2,i));j=incoming(3,i)
       if(j<1.or.j>p%owned(s)%nout.or.incoming(4,i)/=p%owned(s)%depth) call die('producer shape/version')
       q%source(:,i)=[s,j]
    end do
    allocate(q%send(sum(q%rc)),q%received(u),q%value(n))
  end subroutine
  subroutine program_read(p,q)
    type(Owner_Program),intent(in)::p
    type(Owner_Read),intent(inout)::q
    integer::i,s,j,ierr
    if(p%layer<=0) call die('read before execution')
    do i=1,size(q%send)
       s=q%source(1,i);j=q%source(2,i)
       if(.not.p%owned(s)%ready) call die('producer not ready in current layer')
       q%send(i)=p%owned(s)%value(j)
    end do
    call MPI_Alltoallv(q%send,q%rc,q%rd,MPI_DOUBLE_PRECISION, &
         q%received,q%sc,q%sd,MPI_DOUBLE_PRECISION,comm,ierr)
    if(ierr/=MPI_SUCCESS) call die('result delivery')
    do i=1,size(q%value)
       q%value(i)=q%orientation(i)*q%received(q%address(i))
    end do
  end subroutine
  subroutine program_execute(p,k,coefficients,external)
    type(Owner_Program),intent(inout)::p
    integer,intent(in)::k
    real(dp),intent(in),optional::coefficients(:)
    procedure(External_Operand),optional::external
    real(dp),allocatable::coeff(:)
    integer::depth,i,a,n
    if(.not.p%compiled) call die('execute before compile')
    coeff=[0.0_dp,0.0_dp]
    if(present(coefficients)) coeff=coefficients
    p%layer=k;p%owned%ready=.false.
    call fetch_resident_operands(p%resident,k)
    do depth=0,p%maxdepth
       call program_read(p,p%read(depth))
       do i=1,size(p%owned)
          associate(o=>p%owned(i))
          if(o%depth/=depth) cycle
          a=o%input_start;n=size(o%input,2)
          if(o%code==OP_RESIDENT) then
             o%value=p%resident%recv(:,p%resident%address(o%source))
          else if(o%code>=100) then
             if(.not.present(external)) call die('missing external producer callback')
             if(o%nout/=1.or.size(o%input,2)/=0) call die('external producer shape')
             o%value(1)=external(o%code,o%key,k)
          else
             call evaluate(o%code,p%read(depth)%value(a+1:a+n),o%geometry,coeff,o%value(1:o%nout))
          end if
          if(.not.all(ieee_is_finite(o%value(1:o%nout)))) call die('nonfinite result')
          o%ready=.true.
          end associate
       end do
    end do
  end subroutine
  function program_value(p,key,k) result(value)
    type(Owner_Program),intent(in)::p
    integer,intent(in)::key(5),k
    real(dp)::value
    integer::s
    if(.not.p%compiled.or.p%layer/=k.or.key(1)/=rank) call die('external result epoch/owner')
    s=locate(p,key(2),key(3))
    if(.not.p%owned(s)%ready) call die('external result not ready')
    if(abs(key(4))<1.or.abs(key(4))>p%owned(s)%nout) call die('external result component')
    value=sign(1,key(4))*p%owned(s)%value(abs(key(4)))
  end function
  subroutine evaluate(code,x,g,coeff,y)
    integer,intent(in)::code
    real(dp),intent(in)::x(:),g(:),coeff(:)
    real(dp),intent(out)::y(:)
    type(Mass_Restriction)::s
    integer::i
    real(dp)::grad,u(3),qperp(3),nu
    select case(code)
    case(OP_ADD)
       y(1)=x(1)+x(2)
    case(OP_CURL)
       y(1)=-(x(1)-x(2))/g(1)
    case(OP_CIRCULATION)
       u=x*g(1:3)
       select case(nint(g(5)))
       case(1)
          y(1)=(u(1)+u(2)+u(3))/g(4)
       case(2)
          y(1)=-(u(1)+u(2)+u(3))/g(4)
       case(3)
          y(1)=(u(1)-u(2)+u(3))/g(4)
       case(4)
          y(1)=(-u(1)+u(2)+u(3))/g(4)
       case(5)
          y(1)=(u(1)-u(2)-u(3))/g(4)
       case default
          call die('circulation expression')
       end select
    case(OP_VELOCITY_PHYSICS)
       nu=coeff(1)
       if(nint(g(9))==1) nu=coeff(3)
       y=velocity_diffusion_source(x(1:4),x(5:8),g(1:3),g(4:6),nint(g(7)),nint(g(8)),nu,coeff(2))
    case(OP_VELOCITY_SOURCE)
       qperp=velocity_qperp(reshape(x(1:27),[3,9]),reshape(x(28:54),[3,9]),reshape(g(1:30),[5,2,3]))
       y=velocity_source(qperp,x(55:57),g(31:33),.true.)
    case(OP_ZERO)
       y=0.0_dp
    case(OP_GRADIENT,OP_GRADIENT_SW)
       y=[x(2)-x(1),x(1)-x(3),x(4)-x(1)]
       if(code==OP_GRADIENT_SW) y=-y
       y=y/g(1:3)*g(4:6)
    case(OP_VELOCITY_FLUX)
       y=x*g
    case(OP_MASS_DIVERGENCE)
       y(1)=-(x(1)-x(2)+x(3)-x(4)+x(5)-x(6))*g(1)
    case(OP_DIVERGENCE)
       y(1)=(x(1)-x(2)+x(3)-x(4)+x(5)-x(6))*g(1)
    case(OP_SCALAR_FLUX)
       i=nint(g(9));y(1)=0.0_dp
       if(nint(g(7))/=0) then
          grad=(x(5)-x(4))/g(i)
          if(i==2) grad=(x(4)-x(5))/g(i)
          y(1)=(-1)**nint(g(7))*coeff(nint(g(8)))*grad*g(i+3)
       end if
       y(2)=(x(3)*g(i+3))*(0.5*(x(1)+x(2)))+y(1)
    case(OP_RESTRICTION)
       s%small=reshape([(i,i=1,36)],[3,3,4]);s%partial_flux=reshape([(i,i=37,54)],[6,3])
       s%partial_node=reshape([(i,i=1,12)],[4,3]);s%coarse_node=reshape([(i,i=13,30)],[6,3])
       s%weight=reshape(g(1:24),[3,2,4]);s%area=reshape(g(25:30),[2,3])
       s%overlap=reshape(g(31:42),[4,3]);s%coarse=reshape(g(43:57),[5,3]);s%edge=g(58:60)>0
       y=mass_restriction_value(s,x(1:54),x(55:84))
    case default
       call die('unknown numerical operation')
    end select
  end subroutine
end module
