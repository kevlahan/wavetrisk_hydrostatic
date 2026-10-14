module parallel_block_rhs_compile_mod
  ! Geometry-only front end. Original physical storage addresses are retained
  ! here for exact remapping; these are not runtime numerical Domain views.
  use kind_mod, only : dp
  use shared_mod, only : N_BDRY,EDGE,TRIAG,NORTH,EAST,SOUTH,WEST,NORTHEAST,SOUTHWEST, &
       level_start,level_end,z_null,IJPLUS,IJMINUS,IPLUSJMINUS,IMINUSJPLUS,eps,radius, &
       N_CHDRN,RESTRCT,penta_node,dx_avg,min_level
  use patch_mod, only : PATCH_SIZE
  use domain_mod, only : Domain,idx
  use domain_ops_mod, only : apply_onescale_to_patch,apply_to_penta_d,apply_interscale_to_patch3
  use parallel_block_primitive_walk_mod, only : Primitive_Walk,walk_primitive_patch
  use parallel_block_mass_mod, only : Mass_Level,Mass_Divergence
  use parallel_block_velocity_mod, only : Native_Velocity_Workspace,Velocity_Level_Program,Velocity_Gradient_Stencil
  use geom_mod, only : geodesic
  use parallel_block_rhs_mod, only : Rhs_Primitive_Stencil,Rhs_Pentagon_Stencil
  use parallel_block_rhs_diffusion_mod, only : Rhs_Diffusion_Level,Rhs_Curl_Stencil
  implicit none
  private
  public :: Rhs_Domain_Operators,compile_rhs_domain_operators,merge_rhs_domain_operators,compile_rhs_dynamics
  public :: Rhs_Dynamics_Level,Rhs_Velocity_Physics
  type :: Rhs_Patch_Phase
     integer :: first_producer=0,last_producer=0,first_pressure=0,last_pressure=0
  end type
  type :: Rhs_Pressure_Level
     type(Rhs_Patch_Phase), allocatable :: patch(:)
     integer, allocatable :: node(:),producer_ordinal(:)
  end type
  type :: Rhs_Rotation_Level
     type(Rhs_Curl_Stencil), allocatable :: curl(:)
     type(Rhs_Pentagon_Stencil), allocatable :: pentagon(:)
  end type
  type :: Rhs_Velocity_Physics
     integer :: node(4)=0,vorticity(4)=0
     real(dp) :: length(3)=0.0_dp,dual_length(3)=0.0_dp
     logical :: pentagon=.false.
  end type
  type :: Rhs_Dynamics_Level
     type(Mass_Divergence), allocatable :: divergence(:)
     type(Velocity_Level_Program) :: source
     type(Rhs_Velocity_Physics), allocatable :: physics(:)
     integer, allocatable :: injection(:,:)
  end type
  type :: Rhs_Domain_Operators
     integer :: domain_id=-1,nnode=0
     type(Rhs_Dynamics_Level), allocatable :: dynamics(:)
     type(Velocity_Gradient_Stencil), allocatable :: gradient(:)
     integer, allocatable :: surface_visits(:)
     type(Rhs_Pressure_Level), allocatable :: pressure(:)
     type(Rhs_Diffusion_Level), allocatable :: diffusion(:)
     type(Rhs_Rotation_Level), allocatable :: rotation(:)
  end type
  type(Rhs_Domain_Operators), pointer :: compiling=>null()
  integer :: current_level=0,pass=0,used=0
contains
  subroutine compile_rhs_domain_operators(dom,mass,output)
    type(Domain), intent(inout) :: dom
    type(Mass_Level), intent(in) :: mass(level_start:)
    type(Rhs_Domain_Operators), intent(out), target :: output
    type(Primitive_Walk) :: walk
    integer :: l,p,j,t,n,s
    if (associated(compiling)) error stop 'RHS geometry compiler is not reentrant'
    output%domain_id=dom%id;output%nnode=dom%node%length
    allocate(output%diffusion(level_start:level_end),output%rotation(level_start:level_end))
    allocate(output%surface_visits(output%nnode),output%pressure(level_start:level_end))
    output%surface_visits=0
    compiling=>output
    do p=3,dom%patch%length
       call apply_onescale_to_patch(compile_surface_visits,dom,p-1,z_null,0,1)
    end do
    do l=level_start,level_end
       current_level=l;output%diffusion(l)%level=l
       n=0
       do j=1,dom%lev(l)%length
          p=dom%lev(l)%elts(j)
          call walk_primitive_patch(dom,p,operations=walk)
          n=n+count(.not.walk%southwest(:walk%count))
       end do
       allocate(output%diffusion(l)%producer(n),output%pressure(l)%producer_ordinal(n));s=0
       output%pressure(l)%producer_ordinal=[(j,j=1,n)]
       allocate(output%pressure(l)%patch(dom%lev(l)%length))
       do j=1,dom%lev(l)%length
          p=dom%lev(l)%elts(j)
          call walk_primitive_patch(dom,p,operations=walk)
          output%pressure(l)%patch(j)%first_producer=s+1
          do t=1,walk%count
             if (.not.walk%southwest(t)) s=s+1
             call primitive_geometry(dom,walk%ids(:,t),walk%southwest(t),output%diffusion(l)%producer(s))
          end do
       end do
       output%diffusion(l)%restriction=mass(l)%restriction
       do j=1,dom%lev(l)%length-1
          output%pressure(l)%patch(j)%last_producer=output%pressure(l)%patch(j+1)%first_producer-1
       end do
       if(dom%lev(l)%length>0) output%pressure(l)%patch(dom%lev(l)%length)%last_producer=s
       do pass=1,2
          used=0
          do j=1,dom%lev(l)%length
             if(pass==2) output%pressure(l)%patch(j)%first_pressure=used+1
             call apply_onescale_to_patch(compile_divergence,dom,dom%lev(l)%elts(j),z_null,0,1)
             if(pass==2) output%pressure(l)%patch(j)%last_pressure=used
          end do
          if (pass==1) allocate(output%diffusion(l)%divergence(used))
       end do
       output%pressure(l)%node=output%diffusion(l)%divergence%target
       do pass=1,2
          used=0
          do j=1,dom%lev(l)%length
             call apply_onescale_to_patch(compile_curl,dom,dom%lev(l)%elts(j),z_null,0,0)
          end do
          if (pass==1) allocate(output%rotation(l)%curl(used))
       end do
       do pass=1,2
          used=0
          call apply_to_penta_d(compile_pentagon,dom,l,z_null)
          if (pass==1) allocate(output%rotation(l)%pentagon(used))
       end do
    end do
    nullify(compiling)
  end subroutine

  subroutine compile_rhs_dynamics(dom,mass,velocity,output)
    ! Reuse the accepted topology compiler's exact mask predicates and mixed
    ! direct/restriction order. No numerical values are copied from its workspace.
    type(Domain), intent(inout) :: dom
    type(Mass_Level), intent(in) :: mass(level_start:)
    type(Native_Velocity_Workspace), intent(in) :: velocity
    type(Rhs_Domain_Operators), intent(out), target :: output
    integer :: l,s,n,p,j,c,ip,q
    call compile_rhs_domain_operators(dom,mass,output)
    allocate(output%dynamics(level_start:level_end))
    output%gradient=velocity%gradient
    compiling=>output
    do l=level_start,level_end
       current_level=l
       associate(a=>output%dynamics(l))
       a%divergence=mass(l)%divergence
       a%source=velocity%level(l)
       n=count(a%source%direct%active)
       allocate(a%physics(n));q=0
       do s=1,size(a%source%direct)
          if(.not.a%source%direct(s)%active) cycle
          q=q+1
          associate(v=>a%source%direct(s),f=>a%physics(q))
          n=v%node(1)
          f%node=v%node(:4)
          f%vorticity=[2*n-1,2*n,2*v%node(5)-1,2*v%node(7)]
          f%length=v%length;f%dual_length=dom%pedlen%elts(3*n-2:3*n)
          do ip=1,12
             if(geodesic(dom%node%elts(n),penta_node(ip))>=1.5*dx_avg(min_level)) cycle
             f%pentagon=.true.;exit
          end do
          end associate
       end do
       do pass=1,2
          used=0
          if(l<level_end) then
             do j=1,dom%lev(l)%length
                p=dom%lev(l)%elts(j)
                do c=1,N_CHDRN
                   if(dom%patch%elts(p+1)%children(c)<=0) cycle
                   call apply_interscale_to_patch3(compile_injection,dom,p,c,z_null,0,1)
                end do
             end do
          end if
          if(pass==1) allocate(a%injection(2,used))
       end do
       end associate
    end do
    nullify(compiling)
  end subroutine

  subroutine compile_injection(dom,p,ip,jp,ic,jc,k,op,dp_,oc,dc)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: p,ip,jp,ic,jc,k,op(N_BDRY+1),dp_(2,N_BDRY+1),oc(N_BDRY+1),dc(2,N_BDRY+1)
    integer :: n
    n=idx(ip,jp,op,dp_)+1
    if(dom%mask_n%elts(n)<RESTRCT) return
    used=used+1
    if(pass==2) compiling%dynamics(current_level)%injection(:,used)=[n,idx(ic,jc,oc,dc)+1]
  end subroutine

  subroutine compile_surface_visits(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: n
    n=idx(i,j,offs,dims)+1
    compiling%surface_visits(n)=compiling%surface_visits(n)+1
  end subroutine

  subroutine primitive_geometry(dom,ids,southwest,p)
    type(Domain), intent(in) :: dom
    integer, intent(in) :: ids(0:7)
    logical, intent(in) :: southwest
    type(Rhs_Primitive_Stencil), intent(inout) :: p
    integer, parameter :: vnode(10)=[0,1,4,0,3,7,4,0,2,3],component(10)=[1,1,1,2,2,2,2,3,3,3]
    integer, parameter :: dual(6)=[1,4,8,3,6,10]
    integer :: j,a
    if (southwest) then
       if (p%target/=ids(0)+1) error stop 'RHS SW producer did not follow its direct producer'
       p%southwest=.true.;p%sw_mass=ids(7)+1
       p%sw_velocity=[3*ids(7)+1,3*ids(7)+3]
       p%sw_edge=reshape([1,ids(4)+1,2,ids(7)+1,3,ids(3)+1],[2,3])
       p%sw_vorticity=[2*ids(7)+1,2*ids(7)+2]
       p%sw_geometry(1:2)=dom%len%elts(p%sw_velocity)
       p%sw_geometry(3:8)=[dom%areas%elts(ids(7)+1)%part(1),dom%areas%elts(ids(3)+1)%part(3), &
            dom%areas%elts(ids(0)+1)%part(5),dom%areas%elts(ids(7)+1)%part(2), &
            dom%areas%elts(ids(0)+1)%part(4),dom%areas%elts(ids(4)+1)%part(6)]
       p%sw_geometry(9:10)=dom%coriolis%elts(p%sw_vorticity)
       p%sw_geometry(11:12)=dom%triarea%elts(p%sw_vorticity)
       return
    end if
    p%target=ids(0)+1;p%mass=ids(:5)+1
    p%vorticity=[2*ids(0)+1,2*ids(0)+2,2*ids(4)+1,2*ids(3)+2]
    do j=1,10
       a=3*ids(vnode(j))+component(j);p%velocity(j)=a;p%geometry(j)=dom%len%elts(a)
    end do
    p%geometry(11:16)=dom%pedlen%elts(p%velocity(dual))
    p%geometry(17:28)=[dom%areas%elts(ids(0)+1)%part(1),dom%areas%elts(ids(2)+1)%part(3), &
         dom%areas%elts(ids(5)+1)%part(5),dom%areas%elts(ids(0)+1)%part(2), &
         dom%areas%elts(ids(5)+1)%part(4),dom%areas%elts(ids(1)+1)%part(6), &
         dom%areas%elts(ids(4)+1)%part(1),dom%areas%elts(ids(0)+1)%part(3), &
         dom%areas%elts(ids(1)+1)%part(5),dom%areas%elts(ids(3)+1)%part(2), &
         dom%areas%elts(ids(2)+1)%part(4),dom%areas%elts(ids(0)+1)%part(6)]
    p%geometry(29:32)=dom%coriolis%elts(p%vorticity)
    p%geometry(33)=dom%areas%elts(ids(0)+1)%hex_inv
    p%geometry(34:37)=dom%triarea%elts(p%vorticity)
  end subroutine

  subroutine compile_divergence(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: id,w,s,sw
    used=used+1
    if (pass==1) return
    id=idx(i,j,offs,dims);w=idx(i-1,j,offs,dims);s=idx(i,j-1,offs,dims);sw=idx(i-1,j-1,offs,dims)
    associate(p=>compiling%diffusion(current_level)%divergence(used))
    p%target=id+1;p%edge=[3*id+1,3*w+1,3*sw+2,3*id+2,3*id+3,3*s+3]
    p%inverse_area=dom%areas%elts(id+1)%hex_inv
    end associate
  end subroutine

  subroutine compile_curl(dom,i,j,k,offs,dims)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: i,j,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer :: id,w,s
    used=used+1
    if (pass==1) return
    id=idx(i,j,offs,dims);w=idx(i-1,j,offs,dims);s=idx(i,j-1,offs,dims)
    associate(p=>compiling%rotation(current_level)%curl(used))
    p%target=id+1;p%vorticity=[2*id+1,2*id+2,2*w+1,2*s+2]
    p%dual_length=dom%pedlen%elts(3*id+1:3*id+3)
    p%diagonal=p%dual_length(2)>eps(radius)
    end associate
  end subroutine

  subroutine compile_pentagon(dom,p,c,offs,dims,k)
    type(Domain), intent(inout) :: dom
    integer, intent(in) :: p,c,k,offs(N_BDRY+1),dims(2,N_BDRY+1)
    integer, parameter :: dx(7)=[0,0,1,0,-1,1,-1],dy(7)=[0,1,0,-1,0,1,-1]
    integer :: i,j,t,node(7),n
    used=used+1
    if (pass==1) return
    i=0;j=0
    if (c==IJPLUS.or.c==IPLUSJMINUS) i=PATCH_SIZE
    if (c==IJPLUS.or.c==IMINUSJPLUS) j=PATCH_SIZE
    do t=1,7
       if ((c==IJPLUS.and.t>=6).or.(c==IJMINUS.and.t==6).or. &
            (c==IPLUSJMINUS.and.t==2).or.(c==IMINUSJPLUS.and.t==3)) then
          node(t)=idx(i,j,offs,dims)+1
       else
          node(t)=idx(i+dx(t),j+dy(t),offs,dims)+1
       end if
    end do
    associate(q=>compiling%rotation(current_level)%pentagon(used))
    q%corner=c;q%mass=node;q%pv_node=node
    do t=1,7
       n=node(t)-1;q%velocity(:,t)=[3*n+1,3*n+2,3*n+3];q%vorticity(:,t)=[2*n+1,2*n+2]
       q%length(:,t)=dom%len%elts(q%velocity(:,t));q%triangle_area(:,t)=dom%triarea%elts(q%vorticity(:,t))
       q%coriolis(:,t)=dom%coriolis%elts(q%vorticity(:,t));q%parts(:,t)=dom%areas%elts(node(t))%part
    end do
    end associate
  end subroutine
  subroutine merge_rhs_domain_operators(input,base,output,with_dynamics)
    ! Validation/assembly bridge: distinct physical storage namespaces remain
    ! distinct when several source Domains contribute to a rank's workspace.
    type(Rhs_Domain_Operators), intent(in) :: input(:)
    integer, allocatable, intent(out) :: base(:)
    type(Rhs_Domain_Operators), intent(out) :: output
    type(Rhs_Domain_Operators) :: part
    logical, optional, intent(in) :: with_dynamics
    logical :: complete
    integer :: d,l,s,e,b,np,nn,ng,gradient_used
    integer,allocatable :: count(:,:),used(:,:)
    allocate(base(size(input)+1))
    base(1)=0
    do d=1,size(input)
       base(d+1)=base(d)+input(d)%nnode
    end do
    output%nnode=base(size(base))
    complete=.false.
    if(size(input)>0) complete=allocated(input(1)%dynamics)
    if(present(with_dynamics)) complete=with_dynamics
    ! Count every destination once. Data copies remain linear in the assembled
    ! program size even when one rank receives many physical Domains.
    allocate(count(15,level_start:level_end),used(15,level_start:level_end));count=0;used=0;ng=0
    do d=1,size(input)
       if(allocated(input(d)%dynamics).neqv.complete) error stop 'RHS inconsistent dynamics assembly'
       if(allocated(input(d)%gradient)) ng=ng+size(input(d)%gradient)
       do l=level_start,level_end
          count(1,l)=count(1,l)+size(input(d)%pressure(l)%patch)
          count(2,l)=count(2,l)+size(input(d)%pressure(l)%node)
          count(3,l)=count(3,l)+size(input(d)%pressure(l)%producer_ordinal)
          count(4,l)=count(4,l)+size(input(d)%diffusion(l)%producer)
          count(5,l)=count(5,l)+size(input(d)%diffusion(l)%restriction)
          count(6,l)=count(6,l)+size(input(d)%diffusion(l)%divergence)
          count(7,l)=count(7,l)+size(input(d)%rotation(l)%curl)
          count(8,l)=count(8,l)+size(input(d)%rotation(l)%pentagon)
          if(complete) then
             count(9,l)=count(9,l)+size(input(d)%dynamics(l)%divergence)
             count(10,l)=count(10,l)+size(input(d)%dynamics(l)%physics)
             count(11,l)=count(11,l)+size(input(d)%dynamics(l)%injection,2)
             count(12,l)=count(12,l)+size(input(d)%dynamics(l)%source%direct)
             count(13,l)=count(13,l)+size(input(d)%dynamics(l)%source%restriction)
             count(14,l)=count(14,l)+size(input(d)%dynamics(l)%source%action)
             count(15,l)=count(15,l)+size(input(d)%dynamics(l)%source%operand)
          end if
       end do
    end do
    allocate(output%diffusion(level_start:level_end),output%rotation(level_start:level_end))
    allocate(output%pressure(level_start:level_end),output%surface_visits(output%nnode))
    if(complete) allocate(output%dynamics(level_start:level_end),output%gradient(ng))
    gradient_used=0
    do l=level_start,level_end
       output%diffusion(l)%level=l
       allocate(output%pressure(l)%patch(count(1,l)))
       allocate(output%pressure(l)%node(count(2,l)))
       allocate(output%pressure(l)%producer_ordinal(count(3,l)))
       allocate(output%diffusion(l)%producer(count(4,l)))
       allocate(output%diffusion(l)%restriction(count(5,l)))
       allocate(output%diffusion(l)%divergence(count(6,l)))
       allocate(output%rotation(l)%curl(count(7,l)))
       allocate(output%rotation(l)%pentagon(count(8,l)))
       if(complete) then
          allocate(output%dynamics(l)%divergence(count(9,l)))
          allocate(output%dynamics(l)%physics(count(10,l)))
          allocate(output%dynamics(l)%injection(2,count(11,l)))
          allocate(output%dynamics(l)%source%direct(count(12,l)))
          allocate(output%dynamics(l)%source%restriction(count(13,l)))
          allocate(output%dynamics(l)%source%action(count(14,l)))
          allocate(output%dynamics(l)%source%operand(count(15,l)))
       end if
    end do
    do d=1,size(input)
       part=input(d);b=base(d)
       if(allocated(part%dynamics).neqv.allocated(output%dynamics)) error stop 'RHS inconsistent dynamics assembly'
       if(allocated(part%gradient)) then
          do s=1,size(part%gradient)
             part%gradient(s)%node(1)=part%gradient(s)%node(1)+b
             if(part%gradient(s)%active) part%gradient(s)%node(2:)=part%gradient(s)%node(2:)+b
          end do
          output%gradient(gradient_used+1:gradient_used+size(part%gradient))=part%gradient
          gradient_used=gradient_used+size(part%gradient)
       end if
       output%surface_visits(base(d)+1:base(d+1))=part%surface_visits
       do l=level_start,level_end
          np=used(4,l);nn=used(2,l)
          part%pressure(l)%node=part%pressure(l)%node+b
          part%pressure(l)%patch%first_producer=part%pressure(l)%patch%first_producer+np
          part%pressure(l)%patch%last_producer=part%pressure(l)%patch%last_producer+np
          part%pressure(l)%patch%first_pressure=part%pressure(l)%patch%first_pressure+nn
          part%pressure(l)%patch%last_pressure=part%pressure(l)%patch%last_pressure+nn
          output%pressure(l)%patch(used(1,l)+1:used(1,l)+size(part%pressure(l)%patch))=part%pressure(l)%patch
          used(1,l)=used(1,l)+size(part%pressure(l)%patch)
          output%pressure(l)%node(used(2,l)+1:used(2,l)+size(part%pressure(l)%node))=part%pressure(l)%node
          used(2,l)=used(2,l)+size(part%pressure(l)%node)
          output%pressure(l)%producer_ordinal(used(3,l)+1:used(3,l)+size(part%pressure(l)%producer_ordinal))= &
               part%pressure(l)%producer_ordinal
          used(3,l)=used(3,l)+size(part%pressure(l)%producer_ordinal)
          do s=1,size(part%diffusion(l)%producer)
             associate(p=>part%diffusion(l)%producer(s))
             p%target=p%target+b;p%mass=p%mass+b;p%velocity=p%velocity+3*b;p%vorticity=p%vorticity+2*b
             if (p%southwest) then
                p%sw_mass=p%sw_mass+b;p%sw_velocity=p%sw_velocity+3*b
                p%sw_edge(2,:)=p%sw_edge(2,:)+b;p%sw_vorticity=p%sw_vorticity+2*b
             end if
             end associate
          end do
          do s=1,size(part%diffusion(l)%restriction)
             associate(p=>part%diffusion(l)%restriction(s))
             p%target=p%target+b;p%small=p%small+3*b
             do e=1,3
                if (.not.p%edge(e)) cycle
                p%partial_flux(:,e)=p%partial_flux(:,e)+3*b
                p%partial_node(:,e)=p%partial_node(:,e)+b;p%coarse_node(:,e)=p%coarse_node(:,e)+b
             end do
             end associate
          end do
          do s=1,size(part%diffusion(l)%divergence)
             associate(p=>part%diffusion(l)%divergence(s))
             p%target=p%target+b;p%edge=p%edge+3*b
             end associate
          end do
          do s=1,size(part%rotation(l)%curl)
             associate(p=>part%rotation(l)%curl(s))
             p%target=p%target+b;p%vorticity=p%vorticity+2*b
             end associate
          end do
          do s=1,size(part%rotation(l)%pentagon)
             associate(p=>part%rotation(l)%pentagon(s))
             p%mass=p%mass+b;p%velocity=p%velocity+3*b;p%vorticity=p%vorticity+2*b;p%pv_node=p%pv_node+b
             end associate
          end do
          if(allocated(part%dynamics)) call merge_dynamics(part%dynamics(l),b,output%dynamics(l),used(9:15,l))
          output%diffusion(l)%producer(used(4,l)+1:used(4,l)+size(part%diffusion(l)%producer))=part%diffusion(l)%producer
          used(4,l)=used(4,l)+size(part%diffusion(l)%producer)
          output%diffusion(l)%restriction(used(5,l)+1:used(5,l)+size(part%diffusion(l)%restriction))=part%diffusion(l)%restriction
          used(5,l)=used(5,l)+size(part%diffusion(l)%restriction)
          output%diffusion(l)%divergence(used(6,l)+1:used(6,l)+size(part%diffusion(l)%divergence))=part%diffusion(l)%divergence
          used(6,l)=used(6,l)+size(part%diffusion(l)%divergence)
          output%rotation(l)%curl(used(7,l)+1:used(7,l)+size(part%rotation(l)%curl))=part%rotation(l)%curl
          used(7,l)=used(7,l)+size(part%rotation(l)%curl)
          output%rotation(l)%pentagon(used(8,l)+1:used(8,l)+size(part%rotation(l)%pentagon))=part%rotation(l)%pentagon
          used(8,l)=used(8,l)+size(part%rotation(l)%pentagon)
       end do
    end do
    if(any(used/=count)) error stop 'RHS assembly count differs'
    if(gradient_used/=ng) error stop 'RHS gradient assembly count differs'
  end subroutine
  subroutine merge_dynamics(part,b,out,used)
    use parallel_block_velocity_mod, only : VELOCITY_DIRECT,VELOCITY_RESTRICT
    type(Rhs_Dynamics_Level), intent(inout) :: part,out
    integer, intent(in) :: b
    integer,intent(inout) :: used(7)
    integer :: s,n
    do s=1,size(part%divergence)
       part%divergence(s)%target=part%divergence(s)%target+b
       if(part%divergence(s)%active) part%divergence(s)%edge=part%divergence(s)%edge+3*b
    end do
    do s=1,size(part%physics)
       part%physics(s)%node=part%physics(s)%node+b
       part%physics(s)%vorticity=part%physics(s)%vorticity+2*b
    end do
    part%injection=part%injection+b
    do s=1,size(part%source%direct)
       part%source%direct(s)%node(1)=part%source%direct(s)%node(1)+b
       if(part%source%direct(s)%active) part%source%direct(s)%node(2:)=part%source%direct(s)%node(2:)+b
    end do
    do s=1,size(part%source%restriction)
       associate(p=>part%source%restriction(s))
       p%target=p%target+b;p%child=p%child+b;p%neighbor=p%neighbor+b
       end associate
    end do
    do s=1,size(part%source%action)
       select case(part%source%action(s))
       case(VELOCITY_DIRECT)
          n=used(4)
       case(VELOCITY_RESTRICT)
          n=used(5)
       case default
          error stop 'RHS invalid velocity action'
       end select
       part%source%operand(s)=part%source%operand(s)+n
    end do
    out%source%direct(used(4)+1:used(4)+size(part%source%direct))=part%source%direct
    used(4)=used(4)+size(part%source%direct)
    out%source%restriction(used(5)+1:used(5)+size(part%source%restriction))=part%source%restriction
    used(5)=used(5)+size(part%source%restriction)
    out%source%action(used(6)+1:used(6)+size(part%source%action))=part%source%action
    used(6)=used(6)+size(part%source%action)
    out%source%operand(used(7)+1:used(7)+size(part%source%operand))=part%source%operand
    used(7)=used(7)+size(part%source%operand)
    out%divergence(used(1)+1:used(1)+size(part%divergence))=part%divergence
    used(1)=used(1)+size(part%divergence)
    out%physics(used(2)+1:used(2)+size(part%physics))=part%physics
    used(2)=used(2)+size(part%physics)
    out%injection(:,used(3)+1:used(3)+size(part%injection,2))=part%injection
    used(3)=used(3)+size(part%injection,2)
  end subroutine
end module
